import sys
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# Add current dir to path to import pico_adc_interface
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface

def debug_twos_complement(tx=1, rx=2, vref=1.0, pre_trigger_offset_us=5.0, initial_target_tof_us=70.0, tolerance_pct=0.25, capture_time_us=400, sample_rate_hz=5_000_000):
    print("=== Two's Complement ADC Debugging ===")
    
    adc = PicoADCInterface(vref=vref)
    if not adc.connect():
        print("[X] Could not connect to Pico ADC.")
        return
        
    n_samples = int((capture_time_us / 1_000_000.0) * sample_rate_hz)
    print(f"\nConfiguring ADC for {capture_time_us} µs ({n_samples} samples at {sample_rate_hz/1e6:.1f} MSPS)...")
    adc.configure(sample_rate_hz=sample_rate_hz, n_samples=n_samples)

    print(f"\nTriggering capture (TX {tx} -> RX {rx})...")
    volts_original, sample_rate, dt_us = adc.trigger_and_capture(tx=tx, rx=rx)
    
    if volts_original is None:
        print("[X] Capture failed.")
        adc.close()
        return
        
    print(f"[+] Capture successful! {len(volts_original)} samples at {sample_rate/1e6:.1f} MSPS.")
    
    # -------------------------------------------------------------------------
    # Fix Two's Complement Interpretation
    # -------------------------------------------------------------------------
    # Back out the raw 12-bit unsigned code from PicoADCInterface's math
    # PicoADCInterface does: volts = (adc_codes - 2048.0) / 2048.0 * vref
    adc_codes = np.round((volts_original / vref) * 2048.0 + 2048.0).astype(int)
    
    # Interpret as 12-bit Two's Complement
    signed_codes = adc_codes.copy()
    
    # If the MSB (bit 11, value 2048) is set, it's a negative number in two's complement.
    signed_codes[signed_codes >= 2048] -= 4096
    
    # Convert back to actual voltages
    # 12-bit signed ranges from -2048 to +2047
    volts = (signed_codes / 2048.0) * vref
    # -------------------------------------------------------------------------
    
    # Process signal
    n = len(volts)
    time_us = np.arange(n) * dt_us
    
    # 1. Remove DC
    v_ac = volts - np.mean(volts)
    
    # 2. Hilbert envelope
    envelope = np.abs(hilbert(v_ac))
    
    # 3. Noise floor from pre-trigger baseline (~first 50 samples)
    baseline_end = min(50, n // 4)
    noise_floor = np.std(envelope[:baseline_end])
    
    # 4. Search window grounded around baseline ToF
    min_valid_tof = initial_target_tof_us * (1.0 - tolerance_pct)
    max_valid_tof = initial_target_tof_us * (1.0 + tolerance_pct)
    
    search_start_us = pre_trigger_offset_us + min_valid_tof
    search_end_us = pre_trigger_offset_us + max_valid_tof
    
    search_start_idx = int(search_start_us / dt_us)
    search_end_idx = min(int(search_end_us / dt_us), n)
    
    global_peak_idx = None
    peak_time_us = None
    measured_tof_us = None
    
    if search_start_idx < search_end_idx and search_start_idx < n:
        search_region = envelope[search_start_idx:search_end_idx]
        if len(search_region) > 0:
            local_peak_idx = np.argmax(search_region)
            global_peak_idx = search_start_idx + local_peak_idx
            peak_envelope = envelope[global_peak_idx]
            
            snr = peak_envelope / noise_floor if noise_floor > 1e-12 else float('inf')
            if snr >= 3.0:
                delta = 0.0
                if 0 < local_peak_idx < len(search_region) - 1:
                    y0 = search_region[local_peak_idx - 1]
                    y1 = search_region[local_peak_idx]
                    y2 = search_region[local_peak_idx + 1]
                    denom = y0 - 2.0 * y1 + y2
                    if abs(denom) > 1e-15:
                        delta = 0.5 * (y0 - y2) / denom
                
                peak_time_us = (global_peak_idx + delta) * dt_us
                measured_tof_us = peak_time_us - pre_trigger_offset_us
                
                print(f"[+] ToF detected: {measured_tof_us:.3f} µs (SNR: {snr:.1f})")
            else:
                print(f"[-] Peak found but rejected (SNR {snr:.1f} < 3.0)")
    else:
        print("[-] Invalid search window indices.")

    adc.close()

    # Plot the results
    plt.figure(figsize=(12, 6))
    
    # Plot the original misinterpreted signal (light gray, so we can see how bad it was)
    plt.plot(time_us, volts_original, label='Original (Offset Binary interpreted)', color='gray', alpha=0.3)
    
    # Plot the fixed signal
    plt.plot(time_us, volts, label='Fixed (Two\'s Complement)', color='blue', alpha=0.7)
    
    # Plot the envelope of the fixed signal
    plt.plot(time_us, envelope, label='Hilbert Envelope (Fixed)', color='red', linewidth=1.5)
    
    # Show the search window
    plt.axvspan(search_start_us, search_end_us, color='green', alpha=0.2, label='Search Window')
    
    if peak_time_us is not None:
        plt.axvline(peak_time_us, color='purple', linestyle='--', label=f'Peak Time: {peak_time_us:.2f} µs')
        plt.plot(peak_time_us, envelope[global_peak_idx], 'x', color='black', markersize=10, label=f'ToF: {measured_tof_us:.2f} µs')
        
    plt.title(f'ADC Two\'s Complement Debugging (TX {tx} -> RX {rx})')
    plt.xlabel('Time (µs)')
    plt.ylabel('Voltage (V)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    debug_twos_complement()
