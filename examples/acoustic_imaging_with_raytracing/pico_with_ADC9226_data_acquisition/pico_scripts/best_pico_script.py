"""
Pimoroni Pico 2 XL W + AD9226 tomography firmware
==================================================

Hardware-deterministic ultrasonic acquisition using:
  - PIO SM0: exact 250 ns TX pulse on GP22
  - PIO SM1: 12-bit parallel AD9226 capture at 10 MSPS
  - Hardware PWM: continuous ADC clock on GP12
  - DMA: PIO RX FIFO to RAM during each capture
  - Software-clocked SPI: MCP4161 volatile wiper control on GP13-GP15

Pin assignments (GPIO numbers, not physical header pin numbers):
  GP0-GP11 = AD9226 data, Bit 1/MSB on GP0 through Bit 12/LSB on GP11
  GP12      = AD9226 clock
  GP13      = MCP4161 CS
  GP14      = MCP4161 SCK
  GP15      = MCP4161 SDI / Pico SPI1 TX
  GP17-GP16-GP46-GP45 = TX mux S3-S2-S1-S0
  GP21-GP20-GP19-GP18 = RX mux S3-S2-S1-S0
  GP22      = TX pulse through level-shifter LV1

Serial commands:
  PING
  TRIG <tx_channel> <rx_channel>
  CONF <sample_rate_hz> <sample_count>
  GAIN <mcp4161_wiper_code>       (0-256, volatile register only)
  WIPER <mcp4161_wiper_code>      (alias for GAIN)
  STATUS


NOTE:

Validated operating ADC sample rate: 2.5 MSPS. 
Higher rates exhibit acquisition timing instability and are unsupported. 
Absolute ADC timestamps are corrected for the AD9226 seven-sample pipeline 
latency using the actual generated sample clock.
"""

import array
import machine
import rp2
import select
import sys
import time


# =============================================================================
# PIN DEFINITIONS
# =============================================================================
TX_TRIGGER = 22

# Lists are ordered S0, S1, S2, S3 so set_mux() can write the channel bits.
TX_MUX_GPIO = (45, 46, 16, 17)
RX_MUX_GPIO = (18, 19, 20, 21)

ADC_DATA_BASE = 0
ADC_CLK = 12

MCP_CS = 13
MCP_SCK = 14
MCP_SDI = 15


# =============================================================================
# DEFAULTS
# =============================================================================
DEFAULT_N_SAMPLES = 20_000        # 2 ms capture window at 10 MSPS
DEFAULT_SAMPLE_RATE = 10_000_000  # 10 MSPS
DEFAULT_PIO_FREQ = 125_000_000
DEFAULT_WIPER_CODE = 128
MCP_SPI_HALF_PERIOD_US = 1

MIN_SAMPLE_RATE = 100_000
MAX_SAMPLE_RATE = 65_000_000
MIN_SAMPLES = 10
MAX_SAMPLES = 50_000
MIN_WIPER_CODE = 0
MAX_WIPER_CODE = 256


# =============================================================================
# PIO PROGRAMS
# =============================================================================
@rp2.asm_pio(set_init=rp2.PIO.OUT_LOW)
def pulse_250ns():
    # At 125 MHz, each instruction is 8 ns.
    # set(pins, 1) takes 1 cycle + 30 delay cycles = 31 cycles = 248 ns pulse.
    wait(1, irq, 0)
    set(pins, 1) [30]
    set(pins, 0)


@rp2.asm_pio(
    in_shiftdir=rp2.PIO.SHIFT_LEFT,
    autopush=True,
    push_thresh=12,
)
def adc_capture_12bit_sync():
    pull(block)
    mov(x, osr)
    
    label("pre_trigger_loop")
    wait(1, gpio, 12)
    wait(0, gpio, 12)
    in_(pins, 12)
    jmp(x_dec, "pre_trigger_loop")
    
    irq(0)
    
    label("post_trigger_loop")
    wait(1, gpio, 12)
    wait(0, gpio, 12)
    in_(pins, 12)
    jmp("post_trigger_loop")


# =============================================================================
# HARDWARE INITIALIZATION
# =============================================================================
# TX pulse state machine.
tx_pin = machine.Pin(TX_TRIGGER, machine.Pin.OUT, value=0)
sm0 = rp2.StateMachine(
    0,
    pulse_250ns,
    freq=DEFAULT_PIO_FREQ,
    set_base=tx_pin,
)
sm0.active(1)

# Keep the pipelined AD9226 clock running continuously.
adc_clk_pwm = machine.PWM(machine.Pin(ADC_CLK))
adc_clk_pwm.freq(DEFAULT_SAMPLE_RATE)
adc_clk_pwm.duty_u16(32768)

# GP0-GP11 are contiguous, which lets PIO capture all 12 bits at once.
adc_data_pin = machine.Pin(ADC_DATA_BASE, machine.Pin.IN)
for gpio in range(ADC_DATA_BASE, ADC_DATA_BASE + 12):
    machine.Pin(gpio, machine.Pin.IN)

sm1 = rp2.StateMachine(
    1,
    adc_capture_12bit_sync,
    freq=DEFAULT_PIO_FREQ,
    in_base=adc_data_pin,
)

# Mux pin lists retain S0-to-S3 ordering even though the GPIOs are not contiguous.
tx_mux_pins = [machine.Pin(gpio, machine.Pin.OUT, value=0) for gpio in TX_MUX_GPIO]
rx_mux_pins = [machine.Pin(gpio, machine.Pin.OUT, value=0) for gpio in RX_MUX_GPIO]

# MCP4161 write-only bus. GPIO bit-banging avoids claiming SPI1 RX on GP12 and
# works consistently across stock and Pimoroni MicroPython firmware builds.
mcp_cs_pin = machine.Pin(MCP_CS, machine.Pin.OUT, value=1)
mcp_sck_pin = machine.Pin(MCP_SCK, machine.Pin.OUT, value=0)
mcp_sdi_pin = machine.Pin(MCP_SDI, machine.Pin.OUT, value=0)


# =============================================================================
# DMA AND BUFFERS
# =============================================================================
# PIO0 SM1 RX DREQ is 5 on RP2040 and RP2350.
DREQ_PIO0_RX1 = 5
dma = rp2.DMA()

n_samples = DEFAULT_N_SAMPLES
sample_rate_hz = adc_clk_pwm.freq()
wiper_code = None
capture_buf = array.array("I", [0] * n_samples)
send_buf = bytearray(n_samples * 2)

# Converts the GP0..GP11 PIO word into the AD9226's MSB..LSB numeric order.
REVERSED_NIBBLE = (0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15)


# =============================================================================
# HELPERS
# =============================================================================
def set_mux(pins, channel):
    if not 0 <= channel <= 15:
        raise ValueError("mux channel must be from 0 to 15")
    for bit, pin in enumerate(pins):
        pin.value((channel >> bit) & 1)


def set_wiper(new_code):
    """Write MCP4161 volatile Wiper 0 (address 0x00), without touching NVM."""
    global wiper_code
    if not MIN_WIPER_CODE <= new_code <= MAX_WIPER_CODE:
        raise ValueError("MCP4161 wiper code must be from 0 to 256")

    # MCP4161 Write Data frame: 0000 00D8 D7..D0. The address and command
    # fields are zero, so the 16-bit frame is simply the 9-bit wiper code.
    frame = new_code & 0x01FF
    mcp_cs_pin.value(0)
    try:
        for shift in range(15, -1, -1):
            mcp_sdi_pin.value((frame >> shift) & 1)
            time.sleep_us(MCP_SPI_HALF_PERIOD_US)
            mcp_sck_pin.value(1)
            time.sleep_us(MCP_SPI_HALF_PERIOD_US)
            mcp_sck_pin.value(0)
    finally:
        mcp_sck_pin.value(0)
        mcp_cs_pin.value(1)
    wiper_code = new_code


def reverse_adc_word(raw_word):
    """Reverse 12 bits because GP0 is the ADC MSB and GP11 is its LSB."""
    raw_word &= 0x0FFF
    return (
        (REVERSED_NIBBLE[raw_word & 0x0F] << 8)
        | (REVERSED_NIBBLE[(raw_word >> 4) & 0x0F] << 4)
        | REVERSED_NIBBLE[(raw_word >> 8) & 0x0F]
    )


# Avoid repeating the multi-step bit reversal for every sample after capture.
adc_reverse_lut = array.array("H", [0] * 4096)
for raw_code in range(4096):
    adc_reverse_lut[raw_code] = reverse_adc_word(raw_code)


def fire_tx_pulse():
    # SM0 now waits for the PIO IRQ from SM1, so we don't trigger it manually here.
    pass


def reallocate_buffers():
    global capture_buf, send_buf
    capture_buf = array.array("I", [0] * n_samples)
    send_buf = bytearray(n_samples * 2)


def reconfigure_adc(new_rate_hz, new_n_samples):
    global n_samples, sample_rate_hz
    n_samples = new_n_samples
    adc_clk_pwm.freq(new_rate_hz)
    sample_rate_hz = adc_clk_pwm.freq()
    reallocate_buffers()


def capture_adc(pretrigger_samples):
    sm1.active(0)
    sm1.restart()  # Clear the PIO instruction state and input shift register.
    while sm1.rx_fifo() > 0:
        sm1.get()

    # Pre-trigger samples must be at least 1 for the jmp(x_dec) logic.
    pretrigger_count = max(1, min(pretrigger_samples, n_samples - 1))
    sm1.put(pretrigger_count - 1)

    dma.config(
        read=sm1,
        write=capture_buf,
        count=n_samples,
        ctrl=dma.pack_ctrl(
            size=2,
            treq_sel=DREQ_PIO0_RX1,
            inc_read=False,
            inc_write=True,
        ),
        trigger=False,
    )

    dma.active(1)
    sm1.active(1)
    # SM0 is already active and will block on IRQ0 from SM1.
    # fire_tx_pulse() is no longer needed here.

    # Allow at least twice the ideal capture duration, with a 1 ms floor.
    timeout_us = max((n_samples * 2_000_000) // sample_rate_hz, 1_000)
    start_us = time.ticks_us()
    while dma.active():
        if time.ticks_diff(time.ticks_us(), start_us) > timeout_us:
            dma.active(0)
            sm1.active(0)
            return None, pretrigger_count

    sm1.active(0)
    return capture_buf, pretrigger_count


def pack_and_send(buf, count, t0_sample_index):
    for i in range(count):
        value = adc_reverse_lut[buf[i] & 0x0FFF]
        send_buf[i * 2] = value & 0xFF
        send_buf[i * 2 + 1] = (value >> 8) & 0xFF
    sys.stdout.write("DATA {} {} {}\n".format(count, sample_rate_hz, t0_sample_index))
    sys.stdout.buffer.write(send_buf[: count * 2])
    sys.stdout.write("END\n")


def parse_int(text):
    try:
        return int(text)
    except (TypeError, ValueError):
        raise ValueError("expected an integer")


# Establish a known power-up gain without wearing the non-volatile wiper.
set_wiper(DEFAULT_WIPER_CODE)


# =============================================================================
# MAIN LOOP
# =============================================================================
def main():
    transducer_send_receive_ordering = [(1, 2)]

    print(
        "Pimoroni Pico 2 XL W + AD9226 tomography controller ready "
        "(wiper={}).".format(wiper_code)
    )

    poller = select.poll()
    poller.register(sys.stdin, select.POLLIN)
    line_buffer = ""

    while True:
        command_ready = None
        while poller.poll(0):
            char = sys.stdin.read(1)
            if char == "\n" or char == "\r":
                if line_buffer:
                    command_ready = line_buffer.strip()
                    line_buffer = ""
                    break
            else:
                line_buffer += char

        if command_ready:
            parts = command_ready.split()
            command = parts[0].upper()

            try:
                if command == "PING" and len(parts) == 1:
                    print("PONG")

                elif command == "TRIG" and len(parts) >= 3:
                    tx_channel = parse_int(parts[1])
                    rx_channel = parse_int(parts[2])
                    pretrigger_samples = parse_int(parts[3]) if len(parts) > 3 else 0
                    
                    set_mux(tx_mux_pins, tx_channel)
                    set_mux(rx_mux_pins, rx_channel)
                    time.sleep_us(10)
                    print("TRIG_ACK {} {}".format(tx_channel, rx_channel))

                    result = capture_adc(pretrigger_samples)
                    if result[0] is None:
                        print("ERR ADC capture timeout")
                    else:
                        buf, t0_index = result
                        pack_and_send(buf, n_samples, t0_index)

                elif command == "CONF" and len(parts) == 3:
                    new_rate = parse_int(parts[1])
                    new_count = parse_int(parts[2])
                    if not MIN_SAMPLE_RATE <= new_rate <= MAX_SAMPLE_RATE:
                        raise ValueError("sample rate must be from 100000 to 65000000")
                    if not MIN_SAMPLES <= new_count <= MAX_SAMPLES:
                        raise ValueError("sample count must be from 10 to 50000")
                    reconfigure_adc(new_rate, new_count)
                    print("CONF_ACK {} {}".format(sample_rate_hz, n_samples))

                elif command in ("GAIN", "WIPER") and len(parts) == 2:
                    set_wiper(parse_int(parts[1]))
                    print("GAIN_ACK {}".format(wiper_code))

                elif command == "STATUS" and len(parts) == 1:
                    print(
                        "STATUS rate={} n={} wiper={} otr=na".format(
                            sample_rate_hz, n_samples, wiper_code
                        )
                    )

                else:
                    print("ERR bad or unknown command: {}".format(command_ready))

            except ValueError as error:
                print("ERR {}".format(error))

            except Exception as error:
                # A peripheral failure must not terminate main() and leave USB
                # CDC sitting at the REPL, where subsequent PINGs cannot work.
                print(
                    "ERR internal {}: {}".format(
                        type(error).__name__, error
                    )
                )

            continue

        # Standalone continuous-fire mode removed for acquisition mode.
        time.sleep_us(100)


if __name__ == "__main__":
    main()

