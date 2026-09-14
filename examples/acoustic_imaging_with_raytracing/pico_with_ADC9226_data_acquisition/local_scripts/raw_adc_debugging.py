import sys
import os
import time
import numpy as np

# Add current dir to path to import pico_adc_interface
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pico_adc_interface import PicoADCInterface

def print_binary_diagnostics(tx=1, rx=2, vref=2.5, n_samples=100, sample_rate_hz=100000):
    print("=== Raw ADC Binary Debugging ===")
    
    adc = PicoADCInterface(vref=vref)
    if not adc.connect():
        print("[X] Could not connect to Pico ADC.")
        return

    print(f"\nConfiguring ADC for {n_samples} samples at {sample_rate_hz/1e6:.1f} MSPS...")
    adc.configure(sample_rate_hz=sample_rate_hz, n_samples=n_samples)

    print(f"\nTriggering capture (TX {tx} -> RX {rx})...")
    volts, actual_rate, dt_us = adc.trigger_and_capture(tx=tx, rx=rx)
    
    if volts is None:
        print("[X] Capture failed.")
        adc.close()
        return

    # Back-calculate the raw 16-bit integer codes that the interface received
    # Formula in interface: volts = (adc_codes - 2048.0) / 2048.0 * self.vref
    # So: adc_codes = (volts / vref) * 2048.0 + 2048.0
    
    adc_codes = np.round((volts / vref) * 2048.0 + 2048.0).astype(int)
    
    print("\n--- ADC Data Printout ---")
    print(f"Total samples captured: {len(adc_codes)}")
    print(f"Sample Rate: {actual_rate} Hz (dt = {dt_us:.3f} us)")
    print(f"{'Idx':>5} | {'Volts':>8} | {'Raw Int':>7} | {'Hex':>6} | {'Binary (16-bit)':>19}")
    print("-" * 55)
    
    # Print the first N samples
    limit = min(n_samples, len(adc_codes))
    for i in range(limit):
        code = adc_codes[i]
        v = volts[i]
        
        # Format as 16-bit binary string (e.g. 0000 1111 0000 1111)
        # Note: We constrain to 16 bits for display even if it overflows 
        # (though it shouldn't if frombuffer used <u2)
        code_16 = code & 0xFFFF
        bin_str = f"{code_16:016b}"
        bin_fmt = f"{bin_str[:4]} {bin_str[4:8]} {bin_str[8:12]} {bin_str[12:16]}"
        
        hex_str = f"0x{code_16:04X}"
        
        print(f"{i:5d} | {v:8.3f} | {code:7d} | {hex_str:6s} | {bin_fmt}")
        
    print("-" * 55)
    
    # Statistical Summary
    print("\n--- Summary Statistics ---")
    min_code, max_code = np.min(adc_codes), np.max(adc_codes)
    print(f"Min Code: {min_code:5d} (0x{min_code & 0xFFFF:04X}) -> {np.min(volts):.3f} V")
    print(f"Max Code: {max_code:5d} (0x{max_code & 0xFFFF:04X}) -> {np.max(volts):.3f} V")
    
    # Check if bits are stuck
    bitwise_or = 0
    bitwise_and = 0xFFFF
    for code in adc_codes:
        bitwise_or |= (code & 0xFFFF)
        bitwise_and &= (code & 0xFFFF)
        
    print(f"\nBitwise OR of all samples:  {bitwise_or:016b} (Any 0 here means bit is ALWAYS 0)")
    print(f"Bitwise AND of all samples: {bitwise_and:016b} (Any 1 here means bit is ALWAYS 1)")
    
    if bitwise_or == 0:
        print("WARNING: All data is completely zero. AD9226 may not be powered or clocked.")
    elif bitwise_and == 0xFFFF:
        print("WARNING: All data is completely ones (0xFFFF). ADC pins might be floating high.")
    else:
        # Check for typical floating pin patterns
        toggling_bits = bitwise_or ^ bitwise_and
        print(f"Toggling bits:              {toggling_bits:016b} (1 = bit changes state, 0 = stuck)")
        
    adc.close()

if __name__ == '__main__':
    # You can change n_samples here to see more or fewer samples.
    # By default, looking at the first 100 samples is enough to see binary patterns.
    print_binary_diagnostics(n_samples=100)
