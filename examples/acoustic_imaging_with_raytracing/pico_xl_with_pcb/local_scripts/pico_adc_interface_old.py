"""
Pico 2W + AD9226 ADC Serial Interface
======================================

PC-side interface for communicating with the Pico 2W running pico_adc_tomography.py.
Handles auto-detection of the Pico COM port, serial protocol, binary data reception,
and ADC code-to-voltage conversion.

Replaces hantek_interface_testing.py from the old Hantek 6022BE pipeline.

Usage:
    from pico_adc_interface import PicoADCInterface

    adc = PicoADCInterface()
    if adc.connect():
        volts, sample_rate, dt_us = adc.trigger_and_capture(tx=1, rx=2)
        print(f"Captured {len(volts)} samples at {sample_rate/1e6:.1f} MSPS")
        adc.close()
"""

import time
import struct
import numpy as np

try:
    import serial
    import serial.tools.list_ports
    SERIAL_AVAILABLE = True
except ImportError:
    SERIAL_AVAILABLE = False
    print("[!] PySerial not installed. Run: pip install pyserial")


class PicoADCInterface:
    """
    Serial interface to Pico 2W + AD9226 ADC firmware.
    
    Handles:
    - Auto-detection of Pico 2W COM port (VID 0x2E8A)
    - PING/PONG connection verification
    - TRIG command with binary ADC data reception
    - CONF command for runtime reconfiguration
    - GAIN command for MCP4161 volatile wiper control
    - 12-bit unsigned ADC code -> voltage conversion
    
    The AD9226 outputs unsigned 12-bit codes:
        0    = -V_ref (most negative)
        2048 = 0V     (midscale)
        4095 = +V_ref (most positive)
    
    Default V_ref = 1.0V (AD9226 internal reference), but this depends on
    the analog front-end configuration. Adjust via the vref parameter.
    """
    
    PICO_VID = "2E8A"  # Raspberry Pi USB Vendor ID
    
    def __init__(self, port=None, baudrate=115200, timeout=2.0, vref=1.0):
        """
        Args:
            port: COM port string (e.g., "COM3"). None for auto-detect.
            baudrate: Serial baud rate (irrelevant for USB CDC, but required by PySerial).
            timeout: Serial read timeout in seconds.
            vref: AD9226 reference voltage. ADC output spans -vref to +vref.
        """
        self.ser = None
        self.port = port
        self.baudrate = baudrate
        self.timeout = timeout
        self.vref = vref
        
        # Cached from last capture
        self.last_sample_rate = None
        self.last_n_samples = None
        self.last_wiper_code = None
    
    def connect(self):
        """
        Connects to the Pico 2W over USB serial.
        Auto-detects the COM port if not specified.
        Verifies connection with PING/PONG handshake.
        
        Returns: True if connected and verified, False otherwise.
        """
        if not SERIAL_AVAILABLE:
            print("[X] PySerial not available.")
            return False
        
        # Try specified port first
        if self.port:
            if self._try_connect(self.port):
                return True
        
        # Auto-detect Pico serial port
        ports = serial.tools.list_ports.comports()
        for p in ports:
            hwid_upper = p.hwid.upper() if p.hwid else ""
            desc_upper = p.description.upper() if p.description else ""
            
            if (self.PICO_VID in hwid_upper or 
                "PICO" in desc_upper or 
                "MICROPYTHON" in desc_upper or
                "BOARD IN FS MODE" in desc_upper):
                if self._try_connect(p.device):
                    return True
        
        print("[X] Pico 2W USB serial port not detected.")
        print("    Checklist:")
        print("    1. Is the Pico 2W powered and running pico_adc_tomography.py?")
        print("    2. Is Thonny, PuTTY, or another serial monitor closed?")
        print("    3. Is the USB cable a data cable (not charge-only)?")
        return False
    
    def _try_connect(self, port):
        """Attempts to connect to a specific port and verify with PING."""
        try:
            self.ser = serial.Serial(port, self.baudrate, timeout=self.timeout)
            self.port = port
            time.sleep(0.1)  # Let USB CDC settle
            
            # Flush any boot messages from the Pico
            self.ser.reset_input_buffer()
            
            # Verify with PING/PONG
            if self.ping():
                print(f"[+] Connected to Pico 2W on {port}")
                return True
            else:
                print(f"[*] Found device on {port} but PING failed")
                self.ser.close()
                self.ser = None
                return False
                
        except Exception as e:
            if self.ser:
                try:
                    self.ser.close()
                except Exception:
                    pass
                self.ser = None
            return False
    
    def ping(self):
        """
        Sends PING and expects PONG response.
        Returns True if the Pico responds correctly.
        """
        if not self.ser or not self.ser.is_open:
            return False
        
        try:
            self.ser.reset_input_buffer()
            self.ser.write(b"PING\n")
            self.ser.flush()
            
            # Read response lines (may have boot messages before PONG)
            deadline = time.time() + 2.0
            while time.time() < deadline:
                line = self.ser.readline().decode('utf-8', errors='ignore').strip()
                if line == "PONG":
                    return True
            return False
        except Exception:
            return False
    
    def trigger_and_capture(self, tx=1, rx=2):
        """
        Commands the Pico to fire a TX pulse and capture the ADC response.
        
        Args:
            tx: TX transducer channel (MUX select)
            rx: RX transducer channel (MUX select)
        
        Returns:
            (volts_array, sample_rate_hz, dt_us) on success
            (None, None, None) on failure
            
            volts_array: numpy float64 array of voltage values
            sample_rate_hz: actual sample rate used
            dt_us: time between samples in microseconds
        """
        if not self.ser or not self.ser.is_open:
            print("[X] Not connected")
            return None, None, None
        
        try:
            # Flush input buffer
            self.ser.reset_input_buffer()
            
            # Send trigger command
            cmd = f"TRIG {tx} {rx}\n"
            self.ser.write(cmd.encode('ascii'))
            self.ser.flush()
            
            # The Pico acknowledges TRIG before acquiring and packing the data.
            # Packing 20,000 reversed samples in MicroPython can take several
            # seconds, so keep reading until DATA or an explicit firmware error.
            header = None
            last_line = None
            deadline = time.time() + 15.0
            while time.time() < deadline:
                line = self._read_line_timeout(
                    timeout=max(0.05, deadline - time.time())
                )
                if line is None:
                    break
                last_line = line
                if line.startswith("DATA"):
                    header = line
                    break
                if line.startswith("ERR"):
                    print(f"[X] Pico capture error: {line}")
                    return None, None, None

            if header is None:
                print(f"[X] Expected DATA header, last response: {last_line!r}")
                return None, None, None
            
            parts = header.split()
            if len(parts) < 3:
                print(f"[X] Malformed DATA header: {header}")
                return None, None, None
            
            n_samples = int(parts[1])
            sample_rate = int(parts[2])
            
            # Read binary payload: n_samples * 2 bytes (uint16 LE)
            n_bytes = n_samples * 2
            raw_bytes = self._read_exact(n_bytes, timeout=5.0)
            if raw_bytes is None or len(raw_bytes) != n_bytes:
                got = len(raw_bytes) if raw_bytes else 0
                print(f"[X] Expected {n_bytes} binary bytes, got {got}")
                return None, None, None
            
            # Read END marker
            end_line = self._read_line_timeout(timeout=2.0)
            if end_line is None or "END" not in end_line:
                print(f"[*] Warning: expected END marker, got: {repr(end_line)}")
                # Continue anyway — we have the data
            
            # Unpack uint16 LE array
            adc_codes = np.frombuffer(raw_bytes, dtype='<u2').astype(np.float64)
            
            # Convert 12-bit unsigned to voltage
            # AD9226: 0 = -Vref, 2048 = 0V, 4095 = +Vref
            volts = (adc_codes - 2048.0) / 2048.0 * self.vref
            
            self.last_sample_rate = sample_rate
            self.last_n_samples = n_samples
            
            dt_us = 1.0 / (sample_rate / 1e6)
            
            return volts, sample_rate, dt_us
            
        except Exception as e:
            print(f"[X] Capture error: {e}")
            return None, None, None
    
    def configure(self, sample_rate_hz=10_000_000, n_samples=20000):
        """
        Reconfigures the Pico's ADC capture parameters.
        
        Args:
            sample_rate_hz: ADC sample rate in Hz (100kHz - 65MHz)
            n_samples: number of samples per capture (10 - 50000)
        
        Returns: True if acknowledged, False otherwise.
        """
        if not self.ser or not self.ser.is_open:
            print("[X] Not connected")
            return False
        
        try:
            self.ser.reset_input_buffer()
            cmd = f"CONF {sample_rate_hz} {n_samples}\n"
            self.ser.write(cmd.encode('ascii'))
            self.ser.flush()
            
            response = self._read_line_timeout(timeout=2.0)
            if response and "CONF_ACK" in response:
                parts = response.split()
                if len(parts) >= 3:
                    self.last_sample_rate = int(parts[1])
                    self.last_n_samples = int(parts[2])
                    print(f"[+] Reconfigured: {self.last_sample_rate/1e6:.1f} MSPS, "
                          f"{self.last_n_samples} samples")
                return True
            else:
                print(f"[X] CONF failed, response: {repr(response)}")
                return False
        except Exception as e:
            print(f"[X] Configure error: {e}")
            return False

    def set_gain(self, wiper_code):
        """
        Sets the MCP4161 volatile wiper code.

        The accepted range is 0-256 (257 positions). Whether a larger code
        produces more or less amplifier gain depends on how P0A/P0W/P0B are
        connected on the analog board. The command does not write MCP4161 NVM.

        Returns: True if acknowledged, False otherwise.
        """
        if not self.ser or not self.ser.is_open:
            print("[X] Not connected")
            return False

        if (self.last_sample_rate == sample_rate_hz
                and self.last_n_samples == n_samples):
            print(
                f"[+] ADC already configured: {sample_rate_hz/1e6:.1f} MSPS, "
                f"{n_samples} samples; CONF skipped"
            )
            return True

        if isinstance(wiper_code, bool) or not isinstance(wiper_code, int):
            raise ValueError("wiper_code must be an integer from 0 to 256")
        if not 0 <= wiper_code <= 256:
            raise ValueError("wiper_code must be from 0 to 256")

        if self.last_wiper_code == wiper_code:
            print(f"[+] MCP4161 wiper already set to {wiper_code}; write skipped")
            return True

        try:
            self.ser.reset_input_buffer()
            self.ser.write(f"GAIN {wiper_code}\n".encode("ascii"))
            self.ser.flush()

            response = self._read_line_timeout(timeout=2.0)
            if response == f"GAIN_ACK {wiper_code}":
                self.last_wiper_code = wiper_code
                print(f"[+] MCP4161 wiper set to {wiper_code}")
                return True

            print(f"[X] GAIN failed, response: {response!r}")
            return False
        except Exception as e:
            print(f"[X] Gain configuration error: {e}")
            return False

    # The firmware accepts WIPER as an alias; expose the same terminology here.
    set_wiper = set_gain
    
    def get_status(self):
        """Queries the Pico for current status."""
        if not self.ser or not self.ser.is_open:
            return None
        
        try:
            self.ser.reset_input_buffer()
            self.ser.write(b"STATUS\n")
            self.ser.flush()
            response = self._read_line_timeout(timeout=2.0)
            if response and response.startswith("STATUS"):
                for field in response.split()[1:]:
                    if "=" not in field:
                        continue
                    key, value = field.split("=", 1)
                    try:
                        if key == "rate":
                            self.last_sample_rate = int(value)
                        elif key == "n":
                            self.last_n_samples = int(value)
                        elif key == "wiper":
                            self.last_wiper_code = int(value)
                    except ValueError:
                        pass
            return response
        except Exception:
            return None
    
    def _read_line_timeout(self, timeout=2.0):
        """Reads a single line from serial with a timeout."""
        self.ser.timeout = timeout
        try:
            line = self.ser.readline()
            if line:
                return line.decode('utf-8', errors='ignore').strip()
            return None
        except Exception:
            return None
    
    def _read_exact(self, n_bytes, timeout=5.0):
        """Reads exactly n_bytes from the serial port."""
        self.ser.timeout = timeout
        data = bytearray()
        deadline = time.time() + timeout
        
        while len(data) < n_bytes and time.time() < deadline:
            remaining = n_bytes - len(data)
            chunk = self.ser.read(remaining)
            if chunk:
                data.extend(chunk)
            else:
                break
        
        return bytes(data) if len(data) == n_bytes else None
    
    def close(self):
        """Closes the serial connection."""
        if self.ser and self.ser.is_open:
            try:
                self.ser.close()
                print("[+] Serial connection closed.")
            except Exception:
                pass
            self.ser = None


def main():
    """Quick connectivity test."""
    print("=== Pico 2W + AD9226 ADC Interface Test ===")
    
    adc = PicoADCInterface()
    if not adc.connect():
        return
    
    status = adc.get_status()
    if status:
        print(f"[+] Pico status: {status}")
    
    print("\nCapturing single waveform...")
    volts, rate, dt_us = adc.trigger_and_capture(tx=1, rx=2)
    
    if volts is not None:
        print(f"[+] Captured {len(volts)} samples at {rate/1e6:.1f} MSPS")
        print(f"    Time step: {dt_us:.3f} µs")
        print(f"    Voltage range: {np.min(volts):.4f} V to {np.max(volts):.4f} V")
        print(f"    Mean: {np.mean(volts):.4f} V, Std: {np.std(volts):.4f} V")
    else:
        print("[X] Capture failed")
    
    adc.close()


if __name__ == "__main__":
    main()
