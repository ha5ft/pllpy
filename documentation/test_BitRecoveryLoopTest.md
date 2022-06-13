![pllpy logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# BitRecoveryLoopTest class

```python
# This class implements a testbed for the BitRecoveryLoop class.
# It uses the following classes from loops.py and algorithms.py:
#   - BitRecoveryLoop class
#       - SquareWaveNCO
#       - GatedIntegrator
#       - EdgePhaseDetector
#       - Type2LoopFilter
#   - SquareWaveNCO
#   - RectPulsePAMGen
#   - ComplexNoiseGen
#   - ComplexChannel
# During the tests we could
#   At instance construction time:
#       - Set the bit pattern to be used
#       - Set the base bit frequency
#       - Set the sampling rate
#       - Set the threshold for the edge phase detector
#   At runtime
#       - Add phase jump to the PAM generation
#       - Add jump to the bit rate
#       - Change the loop bandwidth
#       - Change the signal and noise gain to change the SNR
# The algorithm of the computation of the Kp and Kiparameters:
#   BL    : SSB noise bandwidth of the loop
#   phiPM : desired phase margine in the loop
#   Ts    : sampling time at the input of the loop
#   N     : decimation factor in the loop
#   For the Type2 loop filter:
#       rho = tan((phiPM/90.0)*pi / 2)
#       K = 4 * BL * rho / (1 + rho)
#       w0 = K / rho
#       Kp = K
#       Ki = w0 * Ts * N
#
# Handling in the frequency jump
#   At the constructor
#       dwacc = 0.0
#   At the beginning of run()
#       dwacc += df * 2.0 * math.pi
#   dwacc is connected to the frequency control of the bit NCO
#
# Handling the phase jump
#   In the run() function in the first call of the next() function of the bitnco
#   the dteta parameter is passed. For the following call 0.0 is passed as the
#   dteta parameter.
#
# The block diagram of the testbed:
#
#          fs  w0bit Abit=1.0          N=1
#           |   |    |                  |
#         |------------|bnext[n]   |------------|       |------------|
# dwacc   |  Square    |>--------->|    Rect    |       |  Complex   |
# >>----->|   Wave     |           |    Pulse   |>->|   |  Noise     |>->|
#         |    NCO     |>-----     |    PAMGen  |   |   |  Gen       |   |
#         |------------|bclk[n]    |------------|   |   |------------|   |
# dteta     |   |                       |           |                    |
# >>------->|   |                      bmode        |                    |
#            teta0=0                                | pam[n]             |
#                                                   |                    |
#                                                   |                    |
#                                  |------------|   |                    |
#                 x[n]             |  Complex   |<-<|                    |
#    |<---------------------------<|  Channel   |                        |
#    |                             |            |<-----------------------|
#    |                             |------------|     noise[n]
#    |                                |      |
#    |                               Gs      Gn
#    |
#    |
#    |
#    |
#    |    |---------------------------------------------------|
#    |>-->|                                                   |
#         |                                                   |
#         |                                                   |>--------> Iy[n]
#   Kp    |               BitRecoveryLoop                     |
#   ----->|                                                   |>-------> clk[n]
#   Ki    |                                                   |
#   ----->|                                                   |>--------> v[n]
#         |                                                   |
#         |                                                   |
#         |---------------------------------------------------|
#           |   |       |              |
#          fs  f0   threshold     openloop
#
#
#
# BL      |-----------------------|
# >>----->|                       |     
# phiPM   |                       |>------ Kp
# >>----->| Kp and Ki computation |
#         |                       |>------ Ki
#         |                       |
#         |-----------------------|
#
# The constructor:
#   BitRecoveryLoopTest(fs, f0, length, bmode=0, threshold=1e-2)
#   Parameters:
#       fs      : sampling frequency in Hz
#       f0      : bit frequency in Hz
#       length  : maximum length of the simulation in samples
#       bmode   : bit generating mode
#                 mode == 0 : random bit sequence
#                 mode == 1 : {-1.0,1.0,-1.0,1.0,...} series
#                 mode == 2 : {1.0,1.0,1.0,1.0,...} series
#                 mode == 3 : {-1.0,-1.0,-1.0,-1.0,...} series
#                 mode == 4 : {-1.0,-1.0,1.0,1.0,...} series
#                 default : 0
#       threshold:The threshold for the EdgePhaseDetector in the loop
#                 default : 1e-2
#
# The run() function
#   Inputs:
#       samples: run length in samples
#       BL     : SSB noise bandwidth of the loop in Hz
#       phiPM  : desired phase margin in the loop in degree
#       dteta  : phase jump in the carrier at the beginning of the run in radian
#                default : 0.0
#       df     : frequency jump in the carrier at the beginning of the run in Hz
#                default : 0.0
#       Gs     : signal gain
#                default : 1.0
#       Gn     : noise gain
#                default : 0
#       openloop:open loop mode enable flag
#                0 : closed loop runs the test with closed loop
#                1 : open loop runs the test with open loop
#                default : 0
#                if this flag change from the value in the previous run it rests
#                some parts of the loop.
#   Output:
#       none   : The run function each time it finished processing 10 * N samples
#                print a '.' character on the console.
#       
# The show functions:
#   show_inout(begin, end)
#       shows graphically the internal signals of the test and theloop
#           
#   Input:    show() function
#       begin : beginning of the signal sample range
#       end   : end of the signal sample range
#
#   Output:   show() function
#       none  : The function graphically displays the internal signals of the
#               loops or prints some signal statistict to the console.     
#
```
# Running scenarios

## Test with different bit sequences

### -1,1,-1,1,... bit sequence and phase step of 1 rad

```python
> python3
>>> test=BitRecoveryLoopTest(200000, 160, 10000000, bmode=1, threshold=1e-2)
>>> test.run(400000,  4,  65.6, dteta=1.0, df=0.0, Gs=1.0, Gn=0.0, openloop=0)
>>> test.show_error(0,400000)
>>> exit()
```

See the [results](scenario_bit_recovery_-1,1_bit_sequence_phase_step.md)

### random bit sequence and phase step of 1 rad

```python
> python3
>>> test=BitRecoveryLoopTest(200000, 160, 10000000, bmode=0, threshold=1e-2)
>>> test.run(400000,  4,  65.6, dteta=1.0, df=0.0, Gs=1.0, Gn=0.0, openloop=0)
>>> test.show_error(0,400000)
>>> exit()
```

See the [results](scenario_bit_recovery_random_bit_sequence_phase_step.md)

### -1,1,-1,1,... bit sequence and frequency step of 2 Hz

```python
> python3
>>> test=BitRecoveryLoopTest(200000, 160, 10000000, bmode=1, threshold=1e-2)
>>> test.run(400000,  4,  65.6, dteta=0.0, df=2.0, Gs=1.0, Gn=0.0, openloop=0)
>>> test.show_error(0,400000)
>>> exit()
```

See the [results](scenario_bit_recovery_-1,1_bit_sequence_frequency_step.md)

### random bit sequence and frequency step of 2 Hz

```python
> python3
>>> test=BitRecoveryLoopTest(200000, 160, 10000000, bmode=0, threshold=1e-2)
>>> test.run(400000,  4,  65.6, dteta=0.0, df=2.0, Gs=1.0, Gn=0.0, openloop=0)
>>> test.show_error(0,400000)
>>> exit()
```

See the [results](scenario_bit_recovery_random_bit_sequence_frequency_step.md)

Go back to the [tests page](tests.md)\
Go back to the [start page](../README.md)
