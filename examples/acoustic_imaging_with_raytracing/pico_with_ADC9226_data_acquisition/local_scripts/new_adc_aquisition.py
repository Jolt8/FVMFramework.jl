import sys
import os
import time
import csv
import argparse
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface


class VisualTomographyADCAcquisitionMaster:
    """
    Modified automated tomography acquisition that plots every single capture
    dynamically, allowing the user to visually inspect the Time-of-Flight 
    envelope detection in real-time.
    """
    def __init__(self, vref=2.5, safety_factor=0.10, min_cooldown_ms=5.0,
                 serial_port=None, runs_dir=None,
                 initial_target_tof_us=70.0, tolerance_pct=0.25,
                 pre_trigger_offset_us=130.0, gain_code=None,
                 sample_rate_hz=5_000_000, capture_time_us=200,
                 plot_pause_sec=2.0):
        
        self.adc = PicoADCInterface(port=serial_port, vref=vref)
        self.vref = vref
        self.safety_factor = safety_factor
        self.min_cooldown_ms = min_cooldown_ms
        self.pre_trigger_offset_us = pre_trigger_offset_us
        self.gain_code = gain_code
        
        self.sample_rate_hz = sample_rate_hz
        self.capture_time_us = capture_time_us
        self.n_samples = int((self.capture_time_us / 1_000_000.0) * self.sample_rate_hz)
        self.plot_pause_sec = plot_pause_sec
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.runs_dir = runs_dir if runs_dir else os.path.join(script_dir, "runs")
        
        self.initial_target_tof_us = initial_target_tof_us
        self.tolerance_pct = tolerance_pct
        self.active_baseline_tof = {}
        
        self.transducer_pairs = [
            (1, 2),
            # Add more pairs here as needed
            # (2, 1),
            # (3, 4), (4, 3),
        ]
        
        self.tof_results = {}
        self.last_observed_tof_us = {}
        
        # Enable interactive plotting
        plt.ion()
        self.fig, self.ax = plt.subplots(figsize=(12, 6))

    def connect(self):
        print("=== Initializing Visual AD9226 Tomography Acquisition ===")
        if not self.adc.connect():
            print("[X] Pico 2W connection failed.")
            return False
        
        status = self.adc.get_status()
        if status:
            print(f"[+] Pico status: {status}")

        if self.gain_code is not None:
            self.adc.set_gain(self.gain_code)
            
        print(f"\nConfiguring ADC for {self.capture_time_us} µs ({self.n_samples} samples at {self.sample_rate_hz/1e6:.1f} MSPS requested)...")
        self.adc.configure(sample_rate_hz=self.sample_rate_hz, n_samples=self.n_samples)
        
        # Pull actual achieved hardware configuration
        self.sample_rate_hz = self.adc.last_sample_rate
        self.n_samples = self.adc.last_n_samples
        print(f"[+] Actual ADC Rate: {self.sample_rate_hz/1e6:.6f} MSPS")
        
        print(f"[+] AD9226 visual acquisition connected and configured.")
        return True
    
    def calculate_precision_tof_us(self, volts, dt_us, tx=1, rx=2):
        if volts is None or len(volts) == 0:
            return None
        
        n = len(volts)
        pair_key = (tx, rx)
        
        t0_sample_index = self.adc.last_trigger_sample_index
        if t0_sample_index is None:
            print("[X] Capture rejected: Firmware did not provide t0_sample_index.")
            return None
            
        # AD9226 has a fixed 7-cycle pipeline delay. The data appearing at index i
        # corresponds to the analog voltage sampled 7 clock cycles earlier.
        PIPELINE_DELAY_SAMPLES = 7
        
        # Hardware-synchronized time array
        time_us = (np.arange(n) - t0_sample_index - PIPELINE_DELAY_SAMPLES) * dt_us
        
        v_ac = volts - np.mean(volts)
        envelope = np.abs(hilbert(v_ac))
        
        # Noise floor calculation
        if t0_sample_index > 20:
            noise_floor = np.std(envelope[:t0_sample_index-5])
        else:
            noise_floor = np.std(envelope[-100:])
            
        if noise_floor < 1e-6:
            noise_floor = 1e-6
        
        # Search window for the ultrasonic wave
        # Avoid the immediate ringdown of the TX pulse by starting a bit later.
        search_start_us = 30.0
        search_end_us = min(500.0, time_us[-1])
        
        valid_indices = np.where((time_us >= search_start_us) & (time_us <= search_end_us))[0]
        
        global_peak_idx = None
        peak_time_us = None
        measured_tof_us = None
        snr = 0
        status_text = "No Peak Found / Window Invalid"
        
        if len(valid_indices) > 0:
            search_region = envelope[valid_indices]
            local_peak_idx = np.argmax(search_region)
            global_peak_idx = valid_indices[local_peak_idx]
            peak_envelope = envelope[global_peak_idx]
            
            snr = peak_envelope / noise_floor if noise_floor > 1e-12 else float('inf')
            if snr >= 2.0:
                delta = 0.0
                if 0 < local_peak_idx < len(search_region) - 1:
                    y0 = search_region[local_peak_idx - 1]
                    y1 = search_region[local_peak_idx]
                    y2 = search_region[local_peak_idx + 1]
                    denom = y0 - 2.0 * y1 + y2
                    if abs(denom) > 1e-15:
                        delta = 0.5 * (y0 - y2) / denom
                
                peak_time_us = time_us[global_peak_idx] + (delta * dt_us)
                measured_tof_us = peak_time_us
                
                status_text = f"Valid ToF: {measured_tof_us:.3f} µs (SNR {snr:.1f})"
            else:
                status_text = f"Rejected Peak (Low SNR: {snr:.1f})"

        # Dynamic plotting
        self.ax.clear()
        self.ax.plot(time_us, volts, label='Raw Signal (Two\'s Complement)', color='blue', alpha=0.5)
        #self.ax.plot(time_us, envelope, label='Hilbert Envelope', color='red', linewidth=1.5)
        
        # Highlight t=0 (TX pulse)
        self.ax.axvline(0, color='orange', linestyle='-', linewidth=2, label='TX Pulse (Hardware t=0)')
        
        self.ax.axvspan(search_start_us, search_end_us, color='green', alpha=0.2, label='Search Window')
        
        if peak_time_us is not None:
            self.ax.axvline(peak_time_us, color='purple', linestyle='--', label=f'Peak @ {peak_time_us:.2f} µs')
            self.ax.plot(peak_time_us, envelope[global_peak_idx], 'x', color='black', markersize=10)
            
        self.ax.set_title(f"Scan TX {tx} -> RX {rx} | {status_text}")
        self.ax.set_xlabel("Time Relative to TX Pulse (µs)")
        self.ax.set_ylabel("Voltage (V)")
        self.ax.set_xlim(time_us[0], time_us[-1])
        self.ax.set_ylim(-self.vref, self.vref)
        self.ax.legend(loc="upper right")
        self.ax.grid(True)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
        
        # Pause to let the user see the plot and slow down sampling
        plt.pause(self.plot_pause_sec)
        
        return max(0.0, measured_tof_us) if measured_tof_us is not None else None

    def execute_tomography_scan(self, run_folder_name=None, continuous=True):
        now = datetime.now()
        timestamp = now.strftime("%Y_%m_%d__%H_%M_%S")
        
        if run_folder_name is None:
            run_folder_name = f"adc_visual_tomography_{timestamp}"
        
        target_run_dir = os.path.join(self.runs_dir, run_folder_name)
        os.makedirs(target_run_dir, exist_ok=True)
        
        tof_csv_path = os.path.join(target_run_dir, f"tof_data_{timestamp}.csv")
        
        print(f"\n---> Starting Visual AD9226 Tomography Acquisition...")
        print(f"     [Mode] Visual Continuous (Ctrl+C to stop)")
        
        scan_index = 1
        
        try:
            with open(tof_csv_path, "w", newline="") as f_tof:
                tof_writer = csv.writer(f_tof)
                tof_writer.writerow([
                    "Scan_Index", "Timestamp_s", "RayPath_ID",
                    "Tx_Channel", "Rx_Channel", "TimeOfFlight_us",
                    "Method", "Status"
                ])
                f_tof.flush()
                
                while True:
                    print(f"\n--- [Scan #{scan_index}] ---")
                    
                    for ray_id, (tx, rx) in enumerate(self.transducer_pairs, start=1):
                        ray_timestamp = time.time()
                        print(f"  [Tx {tx} -> Rx {rx}] Capturing...", end="", flush=True)
                        
                        pretrigger_samples = int(self.pre_trigger_offset_us * self.sample_rate_hz / 1e6)
                        volts, actual_rate, dt_us = self.adc.trigger_and_capture(
                            tx=tx, rx=rx, pre_trigger_samples=pretrigger_samples
                        )
                        
                        tof = self.calculate_precision_tof_us(volts, dt_us, tx=tx, rx=rx) if volts is not None else None
                        
                        if tof is not None:
                            status = "OK"
                            self.last_observed_tof_us[(tx, rx)] = tof
                            print(f" ToF: {tof:.3f} µs")
                        else:
                            tof = -1.0
                            status = "NO_SIGNAL/REJECTED"
                            print(f" [{status}]")
                        
                        self.tof_results[(tx, rx)] = tof
                        
                        tof_writer.writerow([
                            scan_index, f"{ray_timestamp:.3f}",
                            ray_id, tx, rx,
                            f"{tof:.4f}" if tof >= 0 else "N/A",
                            "HilbertEnvelope", status
                        ])
                        f_tof.flush()
                        os.fsync(f_tof.fileno())
                    
                    scan_index += 1
                    if not continuous:
                        break
        
        except KeyboardInterrupt:
            print(f"\n[SUCCESS] Acquisition stopped cleanly by user.")
        
        finally:
            self.adc.close()
            print(f"  - ToF data saved to: '{tof_csv_path}'")
        
        return self.tof_results


def main():
    parser = argparse.ArgumentParser(description='Visual AD9226 Automated Tomography Acquisition')
    parser.add_argument('--tof', type=float, default=70.0, help='Expected baseline ToF in µs (default: 70.0)')
    parser.add_argument('--tolerance', type=float, default=0.20, help='ToF acceptance tolerance as fraction (default: 0.25)')
    parser.add_argument('--vref', type=float, default=2.5, help='AD9226 reference voltage in V (default: 2.5)')
    parser.add_argument('--pre-trigger', type=float, default=0.0, help='Pre-trigger offset in µs (default: 130.0)')
    parser.add_argument('--port', type=str, default=None, help='COM port for Pico')
    parser.add_argument('--single', action='store_true', help='Run a single scan')
    parser.add_argument('--pause', type=float, default=2.0, help='Seconds to pause between captures to view the plot')
    parser.add_argument('--capture-time', type=float, default=100.0, help='Capture time in microseconds')
    parser.add_argument('--sample-rate', type=int, default=2_500_000, help='Sample rate in Hz')
    args = parser.parse_args()
    
    master = VisualTomographyADCAcquisitionMaster(
        vref=args.vref,
        serial_port=args.port,
        initial_target_tof_us=args.tof,
        tolerance_pct=args.tolerance,
        pre_trigger_offset_us=args.pre_trigger,
        plot_pause_sec=args.pause,
        capture_time_us=args.capture_time,
        sample_rate_hz=args.sample_rate
    )
    
    if master.connect():
        master.execute_tomography_scan(continuous=not args.single)
        # Keep plot open at the end if it was a single scan
        if args.single:
            plt.ioff()
            plt.show()

if __name__ == "__main__":
    main()
