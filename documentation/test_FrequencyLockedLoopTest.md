![pllpy logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# FrequencyLockedLoop class
```python
# This class implements a testbed for the FrequencyLockedLoop class.
# It uses the following classes from loops.py and algorithms.py:
#   - FrequencyLoop class
#       - ComplexNCO
#       - ComplexMixer
#       - FFTFrequencyDetector
#       - Type2LoopFilter
#   - ComplexNCO
#   - RectPulsePAMGen
#   - ComplexChannel
# During the tests we could
#   At instance construction time:
#       - Set the base carrier frequency
#       - Set the sampling rate
#       - Set the size of the FFT
#       - Set the carrier power threshold
#   At runtime
#       - Add a frequency jump to the carrier frequency
#       - Change the loop bandwidth
#       - change the signal and noise gain to change the SNR
#
# The algorithm of the computation of the Kp and Ki parameters:
# We use Type2 loop filter
#   BL    : SSB noise bandwidth of the loop
#   phiPM : desired phase margine in the loop
#   Ts    : sampling time at the input of the loop
#   Nfft  : FFT size
#   For the Type2 loop filter:
#       rho = tan((phiPM/90.0)*pi / 2)
#       K = 4 * BL * rho / (1 + rho)
#       w0 = K / rho
#       Kp = K
#       Ki = w0 * Ts * Nfft
#
# Handling in the frequency jump
#   At the constructor
#       dwacc = 0.0
#   At the beginning of run()
#       dwacc += df * 2.0 * math.pi
#   dwacc is connected to the frequency control of the carrier NCO
#
# The block diagram of the testbed:
#
#          fs  w0   Abit=1.0  
#           |   |    |                    
#         |------------|                                |------------|
# dwacc   |  Complex   |carrier[n]                      |  Complex   |
# >>----->|            |>-------------------------->|   |  Noise     |>->|
#         |    NCO     |                            |   |  Gen       |   |
#         |------------|                            |   |------------|   |
# dteta=0   |   |    |                              |                    |
# --------->|   |    |                              |                    |
#            teta0=0 |                              | pam[n]             |
#                 inverQ=0                          |                    |
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
#    |>-->|                                                   |>-------> ps[n]
#         |                                                   |>-------> pn[n]
#         |                                                   |
#   Kp    |               FrequencyLockedLoop                 |>-------> Iz[n]
#   ----->|                                                   |>-------> Qz[n]
#   Ki    |                                                   |
#   ----->|                                                   |>--------> e[n]
#         |                                                   |>--------> c[n]
#         |                                                   |>---> binidx[n]
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
#
# The constructor:
#   FrequencyLockedLoopTest(fs, f0, Nfft, length, threshold=1e-2)
#   Parameters:
#       fs      : sampling frequency in Hz
#       f0      : bit frequency in Hz
#       Nfft    : the length of the FFT
#       length  : maximum length of the simulation in samples
#       threshold:The threshold for the FFTFrequencyDetector in the loop
#                 default : 1e-2
#
# The run() function:
#   run(samples, BL, phiPM, df=0.0, Gs=1.0, Gn=0.0, openloop=0)
#   Inputs:
#       samples: run length in samples
#       BL     : SSB noise bandwidth of the loop in Hz
#       phiPM  : desired phase margin in the loop in degree
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
#       shows graphically the internal signals of the test
#   show_loop(begin, end)
#       shows graphically the internal signals of the loop
#           
#   Input:    show() functions
#       begin : beginning of the signal sample range
#       end   : end of the signal sample range
#
#   Output:   show() function
#       none  : The function graphically displays the internal signals of the loop
#               or prints some signal statistict to the console.     
#
```

# Running scenarios

## First experience with frequency locked loop

```python
>python3
>>> from test import *
>>> test=FrequencyLockedLoopTest(2e5, 1e4, 4096,    2e7, threshold=1e-2)
>>> test.run(800000, 4.0,  65.6, df=1e4, Gs=1.0, Gn=0.0, openloop=0)
>>> test.show_loop(0,800000)
>>> exit()
```
See the [results](scenario_first_experience_with_fll.md)

Go back to the [tests page](tests.md)\
Go back to the [start page](../README.md)
