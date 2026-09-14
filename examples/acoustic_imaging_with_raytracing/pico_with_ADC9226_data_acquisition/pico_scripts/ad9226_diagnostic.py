"""AD9226 standalone diagnostic for the Pimoroni Pico 2 XL W wiring."""

import machine
import rp2
import time


ADC_DATA_BASE = 0
ADC_CLK = 12
REVERSED_NIBBLE = (0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15)


def reverse_adc_word(raw_word):
    raw_word &= 0x0FFF
    return (
        (REVERSED_NIBBLE[raw_word & 0x0F] << 8)
        | (REVERSED_NIBBLE[(raw_word >> 4) & 0x0F] << 4)
        | REVERSED_NIBBLE[(raw_word >> 8) & 0x0F]
    )


def read_gpio_word(data_pins):
    raw_word = 0
    for gpio_offset, pin in enumerate(data_pins):
        raw_word |= pin.value() << gpio_offset
    return reverse_adc_word(raw_word)


print("=" * 60)
print("  AD9226 STANDALONE TEST - PICO 2 XL W PINOUT")
print("=" * 60)

print("\n[1] Starting 10 MHz clock on GP12...")
clk = machine.PWM(machine.Pin(ADC_CLK))
clk.freq(10_000_000)
clk.duty_u16(32768)
print("    OK. Reported freq={} Hz".format(clk.freq()))

print("[2] Waiting 100 ms for the ADC pipeline to fill...")
time.sleep_ms(100)
print("    OK.")

print("[3] Direct GPIO reads (GP0=MSB through GP11=LSB):")
data_pins = [
    machine.Pin(ADC_DATA_BASE + offset, machine.Pin.IN) for offset in range(12)
]
for trial in range(3):
    code = read_gpio_word(data_pins)
    print("    Trial {}: code={:4d} (0x{:03X})".format(trial, code, code))
    time.sleep_ms(10)


@rp2.asm_pio(in_shiftdir=rp2.PIO.SHIFT_LEFT, autopush=True, push_thresh=12)
def adc_read():
    wait(1, gpio, 12)
    wait(0, gpio, 12)
    in_(pins, 12)


print("[4] Setting up the clock-synchronized PIO capture...")
sm = rp2.StateMachine(
    0,
    adc_read,
    freq=125_000_000,
    in_base=machine.Pin(ADC_DATA_BASE, machine.Pin.IN),
)
sm.active(1)

deadline = time.ticks_add(time.ticks_ms(), 2_000)
while sm.rx_fifo() == 0 and time.ticks_diff(deadline, time.ticks_ms()) > 0:
    time.sleep_ms(10)

if sm.rx_fifo() == 0:
    print("    ERROR: FIFO remained empty for 2 seconds.")
    print("    Check that the GP12 clock reaches the AD9226 and the Pico input.")
else:
    print("    PIO is capturing. Reading 20 samples:")
    for sample_index in range(20):
        sample_deadline = time.ticks_add(time.ticks_us(), 1_000)
        while (
            sm.rx_fifo() == 0
            and time.ticks_diff(sample_deadline, time.ticks_us()) > 0
        ):
            pass

        if sm.rx_fifo() == 0:
            print("      Sample {:2d}: FIFO timeout".format(sample_index))
            continue

        code = reverse_adc_word(sm.get())
        volts = (code - 2048.0) / 2048.0
        print(
            "      Sample {:2d}: code={:4d} (0x{:03X}) | volts={:+.4f} V".format(
                sample_index, code, code, volts
            )
        )

sm.active(0)

print("\n[5] Final direct GPIO code:")
code = read_gpio_word(data_pins)
print("    code={:4d} (0x{:03X})".format(code, code))

clk.deinit()
print("\n" + "=" * 60)
print("  TEST COMPLETE")
print("=" * 60)
