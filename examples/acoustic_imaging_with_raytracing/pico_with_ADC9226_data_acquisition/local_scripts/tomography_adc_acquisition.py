"""
Master Tomography Acquisition System — AD9226 ADC + Pico 2W
============================================================

Replaces tomography_master_acquisition.py from the old Hantek 6022BE pipeline.

Continuously iterates across all multiplexed transducer pairs, captures
waveforms via the Pico 2W + AD9226 ADC, computes precision Time-of-Flight,
and exports results to timestamped CSV files.

Algorithms (carried over from the Hantek version):
  1. Hilbert Transform Envelope Detection (eliminates cycle skipping)
  2. Sub-sample Parabolic Interpolation (sub-100ns ToF precision at 10 MSPS)
  3. Exponential Moving Average (EMA) baseline tracking with ±25% acceptance window
  4. Adaptive Inter-Pulse Cooldown (+10% safety factor)

Key improvements over the Hantek version:
  - Hardware-deterministic trigger-to-sample synchronization (no random timing)
  - 12-bit ADC resolution (vs Hantek's 8-bit)
  - 10 MSPS sample rate (100ns sample interval, vs Hantek's ~1 MSPS effective)
  - No USB scope driver issues, no exclusive access conflicts

Usage:
    python tomography_adc_acquisition.py
    python tomography_adc_acquisition.py --tof 70.0
    python tomography_adc_acquisition.py --tof 70.0 --tolerance 0.30
"""

import sys
import os
import time
import csv
import argparse
from datetime import datetime
import numpy as np
from scipy.signal import hilbert

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface


class TomographyADCAcquisitionMaster:
    """
    Full automated tomography acquisition using the Pico 2W + AD9226 ADC.
    
    Equivalent to TomographyAcquisitionMaster but uses the deterministic
    PIO + DMA ADC capture instead of the Hantek 6022BE USB scope.
    """
    
    def __init__(self, vref=1.0, safety_factor=0.10, min_cooldown_ms=5.0,
                 serial_port=None, runs_dir=None,
                 initial_target_tof_us=70.0, tolerance_pct=0.25,
                 pre_trigger_offset_us=5.0, gain_code=None):
        """
        Args:
            vref: AD9226 reference voltage (±V).
            safety_factor: Percentage added to observed ToF for inter-pulse cooldown.
            min_cooldown_ms: Minimum acoustic ring-down dissipation time in ms.
            serial_port: COM port for Pico (None for auto-detect).
            runs_dir: Directory for saving CSV run data.
            initial_target_tof_us: Expected ToF for grounding the acceptance filter.
            tolerance_pct: ±percentage around baseline ToF for acceptance window.
            pre_trigger_offset_us: Estimated time from capture start to TX pulse.
            gain_code: Optional MCP4161 volatile wiper code (0-256). If None,
                       leave the firmware's current wiper setting unchanged.
        """
        self.adc = PicoADCInterface(port=serial_port, vref=vref)
        self.vref = vref
        self.safety_factor = safety_factor
        self.min_cooldown_ms = min_cooldown_ms
        self.pre_trigger_offset_us = pre_trigger_offset_us
        self.gain_code = gain_code
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.runs_dir = runs_dir if runs_dir else os.path.join(script_dir, "runs")
        
        # User-grounded adaptive dynamic window filtering
        self.initial_target_tof_us = initial_target_tof_us
        self.tolerance_pct = tolerance_pct
        self.active_baseline_tof = {}  # Per-pair EMA running baseline
        
        # Transducer pairs (matches pico firmware MUX channels)
        self.transducer_pairs = [
            (1, 2),
            #(2, 1),
            #(3, 4), (4, 3),
            #(5, 6), (6, 5),
            #(7, 8), (8, 7),
        ]
        
        self.tof_results = {}
        self.last_observed_tof_us = {}
    
    def connect(self):
        """Connects to the Pico 2W ADC system."""
        print("=== Initializing AD9226 ADC Tomography Acquisition System ===")
        
        if not self.adc.connect():
            print("[X] Pico 2W connection failed.")
            return False
        
        status = self.adc.get_status()
        if status:
            print(f"[+] Pico status: {status}")

        if self.gain_code is not None:
            if not self.adc.set_gain(self.gain_code):
                print("[X] Could not configure the MCP4161 wiper.")
                self.adc.close()
                return False
        
        print(f"[+] AD9226 ADC system connected and ready.")
        print(f"    MCP4161 volatile wiper code: {self.adc.last_wiper_code}")
        print(f"    Resolution: 12-bit | V_ref: ±{self.vref}V")
        return True
    
    def calculate_precision_tof_us(self, volts, dt_us, tx=1, rx=2):
        """
        User-Grounded Adaptive Dynamic Window ToF Extraction.
        
        Pipeline:
          1. Remove DC offset from the waveform
          2. Compute Hilbert analytic signal envelope
          3. Estimate noise floor from pre-trigger baseline
          4. Define search window around active baseline ToF (±tolerance)
          5. Find envelope peak within the grounded search window
          6. Sub-sample parabolic interpolation for precision
          7. Validate against acceptance window, update EMA baseline
        
        Returns: ToF in µs, or None if rejected.
        """
        if volts is None or len(volts) == 0:
            return None
        
        n = len(volts)
        pair_key = (tx, rx)
        
        # Current running baseline (or initial target)
        current_baseline = self.active_baseline_tof.get(
            pair_key, self.initial_target_tof_us
        )
        
        # Dynamic ±tolerance acceptance window
        min_valid_tof = current_baseline * (1.0 - self.tolerance_pct)
        max_valid_tof = current_baseline * (1.0 + self.tolerance_pct)
        
        # 1. Remove DC
        v_ac = volts - np.mean(volts)
        
        # 2. Hilbert envelope
        envelope = np.abs(hilbert(v_ac))
        
        # 3. Noise floor from pre-trigger baseline (~first 50 samples)
        baseline_end = min(50, n // 4)
        noise_floor = np.std(envelope[:baseline_end])
        
        # 4. Search window grounded around baseline ToF
        #    The echo arrives at pre_trigger_offset + ToF into the capture
        search_start_us = self.pre_trigger_offset_us + min_valid_tof
        search_end_us = self.pre_trigger_offset_us + max_valid_tof
        
        search_start_idx = int(search_start_us / dt_us)
        search_end_idx = min(int(search_end_us / dt_us), n)
        
        if search_start_idx >= search_end_idx or search_start_idx >= n:
            return None
        
        # 5. Find envelope peak in search window
        search_region = envelope[search_start_idx:search_end_idx]
        if len(search_region) == 0:
            return None
        
        local_peak_idx = np.argmax(search_region)
        global_peak_idx = search_start_idx + local_peak_idx
        peak_envelope = envelope[global_peak_idx]
        
        # SNR check
        snr = peak_envelope / noise_floor if noise_floor > 1e-12 else float('inf')
        if snr < 3.0:
            return None
        
        # 6. Sub-sample parabolic interpolation
        delta = 0.0
        if 0 < local_peak_idx < len(search_region) - 1:
            y0 = search_region[local_peak_idx - 1]
            y1 = search_region[local_peak_idx]
            y2 = search_region[local_peak_idx + 1]
            denom = y0 - 2.0 * y1 + y2
            if abs(denom) > 1e-15:
                delta = 0.5 * (y0 - y2) / denom
        
        # ToF = peak_time_in_capture - pre_trigger_offset
        peak_time_us = (global_peak_idx + delta) * dt_us
        measured_tof_us = peak_time_us - self.pre_trigger_offset_us
        
        # 7. Validate against acceptance window
        if min_valid_tof <= measured_tof_us <= max_valid_tof:
            # Update running baseline via EMA (85/15 blend)
            updated_baseline = 0.85 * current_baseline + 0.15 * measured_tof_us
            self.active_baseline_tof[pair_key] = updated_baseline
            return max(0.0, measured_tof_us)
        else:
            # Outlier rejected
            return None
    
    def calculate_adaptive_cooldown_sec(self, tx, rx, measured_tof_us):
        """
        Adaptive inter-pulse cooldown:
          Cooldown = (Observed_ToF × (1 + safety_factor)) + ring_down_time
        """
        if measured_tof_us is not None and measured_tof_us > 0:
            safe_tof_us = measured_tof_us * (1.0 + self.safety_factor)
            cooldown_sec = (safe_tof_us / 1e6) + (self.min_cooldown_ms / 1000.0)
        else:
            cooldown_sec = 0.010  # Fallback 10ms
        return cooldown_sec
    
    def execute_tomography_scan(self, run_folder_name=None, continuous=True):
        """
        Continuously iterates across all transducer pairs, captures ADC waveforms,
        computes precision ToF, and exports data to CSV files.
        
        Press Ctrl+C to stop acquisition cleanly.
        """
        now = datetime.now()
        timestamp = now.strftime("%Y_%m_%d__%H_%M_%S")
        
        if run_folder_name is None:
            run_folder_name = f"adc_tomography_{timestamp}"
        
        target_run_dir = os.path.join(self.runs_dir, run_folder_name)
        os.makedirs(target_run_dir, exist_ok=True)
        
        tof_csv_path = os.path.join(target_run_dir, f"tof_data_{timestamp}.csv")
        scope_csv_path = os.path.join(target_run_dir, f"waveform_data_{timestamp}.csv")
        
        print(f"\n---> Starting AD9226 ADC Tomography Acquisition...")
        print(f"     [Algorithm] Hilbert Envelope + Sub-sample Parabolic Interpolation")
        print(f"     [Safety] +{int(self.safety_factor*100)}% ToF Safety Factor | "
              f"Min Cooldown: {self.min_cooldown_ms} ms")
        print(f"     [Baseline] {self.initial_target_tof_us:.1f} µs ± "
              f"{int(self.tolerance_pct*100)}% acceptance window")
        print(f"     [Run Dir] {target_run_dir}")
        print(f"     [Mode] Continuous (Ctrl+C to stop)")
        
        scan_index = 1
        
        try:
            with (open(tof_csv_path, "w", newline="") as f_tof,
                  open(scope_csv_path, "w", newline="") as f_scope):
                
                tof_writer = csv.writer(f_tof)
                scope_writer = csv.writer(f_scope)
                
                # Headers
                tof_writer.writerow([
                    "Scan_Index", "Timestamp_s", "RayPath_ID",
                    "Tx_Channel", "Rx_Channel", "TimeOfFlight_us",
                    "AdaptiveCooldown_ms", "Method", "Status"
                ])
                scope_writer.writerow([
                    "Scan_Index", "Timestamp_s", "RayPath_ID",
                    "Tx_Channel", "Rx_Channel",
                    "SampleIndex", "Time_us", "ADC_Volts"
                ])
                
                f_tof.flush()
                os.fsync(f_tof.fileno())
                f_scope.flush()
                os.fsync(f_scope.fileno())
                
                while True:
                    print(f"\n--- [Scan #{scan_index}] Processing "
                          f"{len(self.transducer_pairs)} ray paths ---")
                    
                    for ray_id, (tx, rx) in enumerate(self.transducer_pairs, start=1):
                        ray_timestamp = time.time()
                        print(f"  [Scan #{scan_index} | Ray #{ray_id}] "
                              f"Tx {tx} -> Rx {rx}...", end="")
                        
                        # 1. Capture ADC waveform
                        volts, sample_rate, dt_us = self.adc.trigger_and_capture(
                            tx=tx, rx=rx
                        )
                        
                        # 2. Save raw waveform to CSV
                        if volts is not None:
                            n_pts = len(volts)
                            for i in range(n_pts):
                                t_us = i * dt_us
                                scope_writer.writerow([
                                    scan_index, f"{ray_timestamp:.3f}",
                                    ray_id, tx, rx, i,
                                    f"{t_us:.4f}", f"{volts[i]:.6f}"
                                ])
                            f_scope.flush()
                            os.fsync(f_scope.fileno())
                        
                        # 3. Compute precision ToF
                        tof = self.calculate_precision_tof_us(
                            volts, dt_us, tx=tx, rx=rx
                        ) if volts is not None else None
                        
                        if tof is not None:
                            status = "OK"
                            self.last_observed_tof_us[(tx, rx)] = tof
                            cooldown_s = self.calculate_adaptive_cooldown_sec(
                                tx, rx, tof
                            )
                            print(f" ToF: {tof:.3f} µs | "
                                  f"Cooldown: {cooldown_s*1000:.1f} ms")
                        else:
                            tof = -1.0
                            status = "NO_SIGNAL"
                            cooldown_s = 0.010
                            print(f" [NO_SIGNAL]")
                        
                        self.tof_results[(tx, rx)] = tof
                        
                        tof_writer.writerow([
                            scan_index, f"{ray_timestamp:.3f}",
                            ray_id, tx, rx,
                            f"{tof:.4f}" if tof >= 0 else "N/A",
                            f"{cooldown_s*1000:.2f}",
                            "HilbertEnvelope", status
                        ])
                        f_tof.flush()
                        os.fsync(f_tof.fileno())
                        
                        # 4. Adaptive cooldown
                        time.sleep(cooldown_s)
                    
                    scan_index += 1
                    if not continuous:
                        break
        
        except KeyboardInterrupt:
            print(f"\n[SUCCESS] Acquisition stopped cleanly by user.")
        
        finally:
            self.adc.close()
            print(f"  - ToF data saved to: '{tof_csv_path}'")
            print(f"  - Waveforms saved to: '{scope_csv_path}'")
        
        return self.tof_results


def parse_gain_code(value):
    """Argparse type for the MCP4161's 257-position wiper."""
    code = int(value)
    if not 0 <= code <= 256:
        raise argparse.ArgumentTypeError("gain code must be from 0 to 256")
    return code


def main():
    parser = argparse.ArgumentParser(
        description='AD9226 ADC Automated Tomography Acquisition',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--tof', type=float, default=70.0,
                        help='Expected baseline ToF in µs (default: 70.0)')
    parser.add_argument('--tolerance', type=float, default=0.25,
                        help='ToF acceptance tolerance as fraction (default: 0.25 = ±25%%)')
    parser.add_argument('--safety', type=float, default=0.10,
                        help='Safety factor for cooldown (default: 0.10 = +10%%)')
    parser.add_argument('--cooldown', type=float, default=5.0,
                        help='Minimum cooldown in ms (default: 5.0)')
    parser.add_argument('--vref', type=float, default=1.0,
                        help='AD9226 reference voltage in V (default: 1.0)')
    parser.add_argument('--gain-code', type=parse_gain_code, default=None,
                        metavar='0..256',
                        help='set MCP4161 volatile wiper (unchanged by default)')
    parser.add_argument('--pre-trigger', type=float, default=5.0,
                        help='Pre-trigger offset in µs (default: 5.0)')
    parser.add_argument('--port', type=str, default=None,
                        help='COM port for Pico (default: auto-detect)')
    parser.add_argument('--single', action='store_true',
                        help='Run a single scan instead of continuous')
    args = parser.parse_args()
    
    print("\n=======================================================")
    print("      AD9226 ADC AUTOMATED TOMOGRAPHY ACQUISITION      ")
    print("=======================================================\n")
    
    tof = args.tof
    print(f"[+] Baseline ToF: {tof:.1f} µs | "
          f"Acceptance: {tof*(1-args.tolerance):.1f} – {tof*(1+args.tolerance):.1f} µs\n")
    
    master = TomographyADCAcquisitionMaster(
        vref=args.vref,
        safety_factor=args.safety,
        min_cooldown_ms=args.cooldown,
        serial_port=args.port,
        initial_target_tof_us=tof,
        tolerance_pct=args.tolerance,
        pre_trigger_offset_us=args.pre_trigger,
        gain_code=args.gain_code,
    )
    
    if master.connect():
        master.execute_tomography_scan(continuous=not args.single)


if __name__ == "__main__":
    main()
