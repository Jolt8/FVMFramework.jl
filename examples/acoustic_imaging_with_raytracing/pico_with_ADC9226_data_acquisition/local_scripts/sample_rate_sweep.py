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

class SampleRateSweep:
    def __init__(self, vref=2.5, serial_port=None, runs_dir=None,
                 pre_trigger_offset_us=0.0, capture_time_us=100.0,
                 num_captures=1):
        self.adc = PicoADCInterface(port=serial_port, vref=vref)
        self.vref = vref
        self.pre_trigger_offset_us = pre_trigger_offset_us
        self.capture_time_us = capture_time_us
        self.num_captures = num_captures
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.runs_dir = runs_dir if runs_dir else os.path.join(script_dir, "runs")
        
        self.tx = 1
        self.rx = 2
        
        self.rates_to_test = np.linspace(
            2_000_000, 5_000_000, 20, dtype=int
        )
        
    def connect(self):
        print("=== Initializing Sample Rate Sweep ===")
        if not self.adc.connect():
            print("[X] Pico 2W connection failed.")
            return False
        return True
        
    def calculate_precision_tof_us(self, volts, dt_us):
        if volts is None or len(volts) == 0:
            return None
            
        n = len(volts)
        t0_sample_index = self.adc.last_trigger_sample_index
        if t0_sample_index is None:
            return None
            
        # AD9226 has a fixed 7-cycle pipeline delay. The data appearing at index i
        # corresponds to the analog voltage sampled 7 clock cycles earlier.
        PIPELINE_DELAY_SAMPLES = 7
        
        # Hardware-synchronized time array
        time_us = (np.arange(n) - t0_sample_index - PIPELINE_DELAY_SAMPLES) * dt_us
        
        v_ac = volts - np.mean(volts)
        envelope = np.abs(hilbert(v_ac))
        
        if t0_sample_index > 20:
            noise_floor = np.std(envelope[:t0_sample_index-5])
        else:
            noise_floor = np.std(envelope[-100:])
            
        if noise_floor < 1e-6:
            noise_floor = 1e-6
            
        search_start_us = 30.0
        search_end_us = min(500.0, time_us[-1])
        
        valid_indices = np.where((time_us >= search_start_us) & (time_us <= search_end_us))[0]
        
        measured_tof_us = None
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
                
                measured_tof_us = time_us[global_peak_idx] + (delta * dt_us)
                
        return max(0.0, measured_tof_us) if measured_tof_us is not None else None

    def execute_sweep(self, run_folder_name=None):
        now = datetime.now()
        timestamp = now.strftime("%Y_%m_%d__%H_%M_%S")
        
        if run_folder_name is None:
            run_folder_name = f"sample_rate_sweep_{timestamp}"
            
        target_run_dir = os.path.join(self.runs_dir, run_folder_name)
        os.makedirs(target_run_dir, exist_ok=True)
        
        csv_path = os.path.join(target_run_dir, f"sweep_data_{timestamp}.csv")
        
        results = []
        
        try:
            with open(csv_path, "w", newline="") as f_csv:
                writer = csv.writer(f_csv)
                writer.writerow(["SampleRate_Hz", "Mean_ToF_us", "StdDev_ToF_us", "Successful_Captures", "Raw_ToFs_us"])
                
                print("\n---> Starting Sample Rate Sweep...")
                
                for rate in self.rates_to_test:
                    n_samples = int((self.capture_time_us / 1_000_000.0) * rate)
                    
                    if n_samples < 50:
                        print(f"\n[Rate: {rate/1e6:.1f} MSPS] Skipping (too few samples for capture window: {n_samples})")
                        continue
                        
                    print(f"\n[Rate: {rate/1e6:.1f} MSPS] Configuring for {n_samples} samples...")
                    self.adc.configure(sample_rate_hz=rate, n_samples=n_samples)
                    
                    # Pull actual achieved hardware configuration
                    actual_rate = self.adc.last_sample_rate
                    print(f"[+] Actual ADC Rate: {actual_rate/1e6:.6f} MSPS")
                    
                    # Pre-trigger samples calculation using ACTUAL rate
                    pretrigger_samples = int(self.pre_trigger_offset_us * actual_rate / 1e6)
                    
                    tofs = []
                    for i in range(self.num_captures):
                        volts, actual_rate, dt_us = self.adc.trigger_and_capture(
                            tx=self.tx, rx=self.rx, pre_trigger_samples=pretrigger_samples
                        )
                        tof = self.calculate_precision_tof_us(volts, dt_us) if volts is not None else None
                        
                        if tof is not None:
                            tofs.append(tof)
                        time.sleep(0.05) # Small pause to avoid overwhelming serial/mcp
                        
                    if len(tofs) > 0:
                        mean_tof = np.mean(tofs)
                        std_tof = np.std(tofs)
                        success_rate = len(tofs)
                    else:
                        mean_tof = 0.0
                        std_tof = 0.0
                        success_rate = 0
                        
                    print(f"  -> Captured {success_rate}/{self.num_captures} valid signals.")
                    if success_rate > 0:
                        print(f"  -> Mean ToF: {mean_tof:.3f} µs (StdDev: {std_tof:.3f} µs)")
                    
                    writer.writerow([actual_rate, mean_tof, std_tof, success_rate, str(tofs)])
                    f_csv.flush()
                    
                    if success_rate > 0:
                        results.append({
                            'rate_mhz': actual_rate / 1e6,
                            'mean': mean_tof,
                            'std': std_tof
                        })
                        
        except KeyboardInterrupt:
            print("\n[!] Sweep stopped by user.")
        finally:
            self.adc.close()
            print(f"\n[+] Sweep data saved to: '{csv_path}'")
            
        if results:
            self.plot_results(results, target_run_dir, timestamp)
            
    def plot_results(self, results, target_run_dir, timestamp):
        rates = [r['rate_mhz'] for r in results]
        means = [r['mean'] for r in results]
        stds = [r['std'] for r in results]
        
        plt.figure(figsize=(10, 6))
        
        # Plot with error bars
        plt.errorbar(rates, means, yerr=stds, fmt='-o', capsize=5, capthick=2, 
                    ecolor='red', markerfacecolor='blue', markersize=6, linewidth=2)
                    
        plt.title('Time of Flight vs. ADC Sample Rate')
        plt.xlabel('ADC Sample Rate (MSPS)')
        plt.ylabel('Average ToF (µs)')
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Annotate std dev above points if there's enough space, or offset
        for i, std in enumerate(stds):
            # Only annotate if std dev is > 0.001 to avoid clutter on perfect runs
            if std > 0.001:
                plt.annotate(f'±{std:.3f}', 
                             (rates[i], means[i]), 
                             textcoords="offset points", 
                             xytext=(0,10), 
                             ha='center', 
                             fontsize=9)
            
        plot_path = os.path.join(target_run_dir, f"sweep_plot_{timestamp}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"[+] Plot saved to: '{plot_path}'")
        
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='ADC Sample Rate ToF Sweep')
    parser.add_argument('--vref', type=float, default=2.5, help='AD9226 reference voltage in V (default: 2.5)')
    parser.add_argument('--pre-trigger', type=float, default=0.0, help='Pre-trigger offset in µs (default: 0.0)')
    parser.add_argument('--port', type=str, default=None, help='COM port for Pico')
    parser.add_argument('--capture-time', type=float, default=100.0, help='Capture time in microseconds')
    parser.add_argument('--captures', type=int, default=10, help='Number of captures per sample rate')
    args = parser.parse_args()
    
    sweep = SampleRateSweep(
        vref=args.vref,
        serial_port=args.port,
        pre_trigger_offset_us=args.pre_trigger,
        capture_time_us=args.capture_time,
        num_captures=args.captures
    )
    
    if sweep.connect():
        sweep.execute_sweep()

if __name__ == "__main__":
    main()
