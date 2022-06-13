![pllpy logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# PLLorCostasLoopTest class
```python
# This class implements a test for the PLLorCostasLoop class in the loops.py file
# It uses the following classes from loops.py and algorithms.py:
#   - PLLorCostasLoop from loops.py which uses from algorithms.py the following classes
#       - ComplexMixer
#       - ComplexMovingAverageDecimator
#       - ComplexVariableMovingAverageFilter
#       - PhaseDetector
#       - Type2LoopFilter
#       - Type3LoopFilter
#       - ComplexNCO
#       - LockDetector
#   - SweeperGen from algoritms.py
#   - ComplexNCO from algoritms.py
#   - ComplexNoiseGen from algoritms.py
#   - ComplexChannel from algoritms.py
#   - SquareWaveNCO from algoritms.py
#   - RectPulsePAMGen from algoritms.py
#   - BPSKModulator from algoritms.py
# In the test the input signal to the PLL or Castas loop could be:
#   - a unmodulated carrier
#   - an unmodulated carrier with given level of white noise
#   - BPSK modulated carrier
#   - BPSK modulated carrier with given level of white noise
# For the carrier we could
#   - apply a given phase jump
#   - apply a frequency jump
#   - apply a triangle shape frequency sweep
#   - a sine function shape frequency sweep
# Instances of this class could be constructed for PLL test or Costas Loop test. The pdmode
# parameter passed to the constructor determines for wich test the instance will be generated.
# The PAM generator (SquareWaveNCO and RectPulsePAMGen) will operate only in the CostasLoop mode.
# Wheter the instance of the PLLorCostasLoop class operates in PLL or CostasLoop mode depend
# on the pdmode parameter passed to its constructor. This behavior can not be changed at
# runtime.
# At runtime we could:
#   - change the the loops noise bandwidth. From the BL and phiPM (phase margine) parameters the
#     Kp and Ki parameters will be computed in the run function.
#   - change the power of the added white noise in the complex channel
#   - apply phase jump on the PLLs input signal
#   - apply frequency jump on the loops input signal
#   - start a triangle or sine shape frequency change on the loops input signal. The sweep
#     will be one sweep period long. The sweep period and sweep amplitude could be controlled
#     at instance construction time.
# The run() function could be invoked multiple times.
# After the run() function completes we could investigate the result using the show() functions.
#
# The algorithm of the computation of the Kp and Ki parameters:
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
#   For the Type3 loop filter:
#       phi = (phiPM + 90) / 2
#       rho = tan((phiPM/90.0)*pi / 2)
#       K = 4 * BL * (2 * rho - 1) / (2 * rho + 3)
#       w0 = K / rho
#       Kp = K
#       Ki = w0 * Ts * N
#
# Handling in the frequency jump
#   At the constructor
#       dwacc = 0.0
#   At the beginning of run()
#       dwacc += df * 2.0 * math.pi
#   dwacc is added to the sweep gwnwrator output to form the frequency control of the carrier NCO
#
# Handling the phase jump
#   In the run() function in the first call of the next() function of the bitnco
#   the dteta parameter is passed. For the following call 0.0 is passed as the
#   dteta parameter.
#
# The block diagramm of the test is the following:
#
# >>------------------------>|
# dwacc                      |
#               w0sw         |          A=1.0                                                    Gn
#           fs   |  Asw      |       fs   |    w0                             Gs=1.0   |<--------<<
#           |    |   |       |       |    |    | carrier[n]         modcarr[n]   |     |
#   swend |------------|     |     |------------|       |------------|       |------------|
#   -----<|            |   |---|   |  Complex   |>----->|   BPSK     |>----->|  Complex   | x[n]
# swstart | SweeperGen |>->| + |>->|    NCO     |       |            |       |            |>----->|
# >>----->|            |   |---|   |  (Carrier) |   |>->| Modulator  |   |>->|  Channel   |       |
#         |------------|           |------------|   |   |------------|   |   |------------|       |
# swmode         |          sweep[n] | |      |     |                    |                        |
# >>------------>|                   | |  invertQ=0 |                    |                        |
#                                    | |            |                    |                        |
# >>-------------------------------->| teta0=0      |pam[n]              |                        |
# dteta                                             |                    |                        |
#                                                   |                    |noise[n]                |
#          fs  w0bit Abit=1.0          N=1          |                    |                        |
#           |   |    |                  |           |                    |                        |
#         |------------|bnext[n]   |------------|   |   |------------|   |                        |
# dwbit=0 |  Complex   |>--------->|    Rect    |   |   |  Complex   |   |                        |
# >>----->|            |           |    Pulse   |>->|   |  Noise     |>->|                        |
#         |    NCO     |>-----     |    PAMGen  |       |  Gen       |                            |
#         |------------|bclk[n]    |------------|       |------------|                            |
# dteta     |   |    |                  |                                                         |
# >>------->|   |  invertQ=0           bmode                                                      |
#            teta0=0                                                                              |
#                                                                                                 |
#    |<------------------------------------------------------------------------------------------<|
#    |
#    |    |---------------------------------------------------|
#    |>-->|                                                   |>--------> z[n]
#   Kp    |                                                   |
#   ----->|                                                   |>--------> Iu[k]
#   Ki    |               PLLorCostasLoop                     |>--------> Qu[k]
#   ----->|                                                   |>--------> e[k]
#         |                                                   |>--------> c[k]
#         |                                                   |
#         |                                                   |
# >>----->|                                                   |>--------> vk[n]
# pdmode  |---------------------------------------------------|
#           |   |   |   |   |                  |       |   |
#          fs  w0   N  Nfir bw               lfsel  pdtype mu
#
#
#
# BL      |-----------------------|
# >>----->|                       |     
# phiPM   |                       |>------ Kp
# >>----->| Kp and Ki computation |
#   lfsel |                       |>------ Ki
#   ----->|                       |
#         |-----------------------|
#
# The constructor:
#   PLLorCostasLoopTest(fs, f0, N, Nfir, bw, f0sw, Asw, f0bit, length, lfsel=1, pdtype=0, bmode=0, mu=1.0)
#   Parameters passed to the constructor:
#       fs     : sampling frequency in Hz
#       f0     : nominal frequency of the carrier and the ComplexNCO in the loop in Hz
#       N      : decimation factor in the loop
#       Nfir   : Length of the FIR decimator
#       bw     : Cutoff frequency of the FIR decimator in Hz
#       f0sw   : Frequency of the sweep waveform
#       Asw    : Amplitude of the sweep waveform
#       f0bit  : bit frequency in Hz
#       length : maximum lenth of the simulation in samples
#       lfsel  : loop filter selector
#                0 : type2
#                1 : type3
#                default is type3
#       pdtype : type of the phase detector
#                0 : Carrier
#                1 : BPSK
#                2 : QPSK
#                default is Carrier
#       bmode  : bit generating mode
#                mode == 0 : random bit sequence
#                mode == 1 : {-1.0,1.0,-1.0,1.0,...} series
#                mode == 2 : {1.0,1.0,1.0,1.0,...} series
#                mode == 3 : {-1.0,-1.0,-1.0,-1.0,...} series
#                mode == 4 : {-1.0,-1.0,1.0,1.0,...} series
#                default : 0
#       mu     : safety factor for phase detector
#                default value: 1.0
# The run fuction:
#   run(samples, BL, phiPM, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
#   Input:     run() function
#       samples: run length in samples
#       BL     : SSB noise bandwidth of the loop in Hz
#       phiPM  : desired phase margin in the loop in degree
#       M      : legth of the variable moving average filrter
#                default : 1
#       dteta  : phase jump in the carrier at the beginning of the run in radian
#                default : 0.0
#       df     : frequency jump in the carrier at the beginning of the run in Hz
#                default : 0.0
#       Gs     : signal gain
#                defaukt : 1.0
#       Gn     : noise gain
#                default : 0
#       startsw: sweep start flag, strats a sweep at the beginning of the run if
#                sweep not already in progress
#                default : 0
#       swmode : sweep mode
#                0 : triangle wave
#                1 : sine wave
#                default : 0
#       pdmode : phase detector mode selector
#                0 : linear mode
#                1 : pricipal value mode
#                default : 0
#       openloop:open loop mode enable flag
#                0 : closed loop runs the test with closed loop
#                1 : open loop runs the test with open loop
#                if this flag change from the value in the previous run it rests
#                default : 0
#                some parts of the loop.
#
#   Output:    run() function
#       none   : The run function each time it finished processing 10 * N samples
#                print a character on the console.
#                o : open loop, no sweep
#                O : open loop, sweep
#                c : closed loop, no sweep
#                C : closed loop, sweep
#
# The show functions:
#   show_inout(begin, end)
#       shows the input and output signals of the loop
#   show_loop(begin, end)
#       shows internal signals of the loop
#   show_stat(begin, end)
#       print signal statistics on the console
#
#   Input:    show() function
#       begin : beginning of the signal sample range
#       end   : end of the signal sample range
#
#   Output:   show() function
#       none  : The function graphically displays the internal signals of the
#               loop or prints some signal statistict to the console.     
```
# Running scenarios

## The effect of delay in the loop

In this scenario we test the effect of a delay in the loop. If you increase the delay at some point the loop start to become instable and finally it wont lock.
In these tests we increase the delay by increasing the legth of the FIR decimator in the loop. We do the test with type2 loop filter.

BL=4Hz, N=1250, Nfir=5000 : it means 2 decimated sample delay becouse the group delay is half the FIR length.
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 20, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 400000,  4,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,400000)
>>> exit()
```

BL=4Hz, N=1250, Nfir=20000 : 8 sample delay. Here the oscillation starts.
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 20000, 20, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 400000,  4,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,400000)
>>> exit()
```

BL=4, N=1250, Nfir=35000 : 14 sample delay. Long oscillation, but stil locks. Long simulation run
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 35000, 300, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 2000000,  4,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,2000000)
>>> exit()
```

BL=4Hz, N=1250, Nfir=40000 : 16 sample delay. Do not locks. Long simulation run.
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 40000, 300, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 2000000,  4,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,2000000)
>>> exit()
```

See the [results](scenario_effects_of_delay.md)

## Loop bandwidth and delay

The allowable delay depends on the loop bandwidth BL. It is inversely proportional with BL.
In these scenario we keep the delay constant and increase BL until the loop wont lock. We do the tests with type2 loop filter.

BL=4Hz, N=1250, Nfir=5000
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 20, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 400000,  4,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,400000)
>>> exit()
```

BL=16Hz, N=1250, Nfir=5000 : oscillation starts
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 300, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 400000,  16,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,400000)
>>> exit()
```
BL=24Hz, N=1250, Nfir=5000 : Long oscillation, but stil locks. Long simulation run
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 300, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 2000000,  24,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,2000000)
>>> exit()
```

BL=28Hz, N=1250, Nfir=5000 : Do not locks. Long simulation run.
```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 300, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 2000000,  28,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,2000000)
>>> exit()
```

See the [results](scenario_loop_bandwidth_and_delay.md)

## Tracking behaviors of the PLL

We do this tests with type3 loop filter.

### Phase step
We start the simulation with 0 phase and frequency error. Run the loop for 1 second. After that we apply a phase step and run the loop for 2 seconds.
BL=4Hz, dteta=1rad

```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 20, 0.025, 1000,  40.0,    2e7, lfsel=1, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 200000,  4,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.run( 600000,  4,  65.6, M=1, dteta=1.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,800000)
>>> exit())
```

See the [results](scenario_tracking_phase_step.md)

### Frequency step

We start the simulation with 0 phase and frequency error. Run the loop for 1 second. After that we apply a frequency step and run the loop for 2 seconds.
BL=4Hz, df=8Hz

```python
> python3
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 20, 0.025, 1000,  40.0,    2e7, lfsel=1, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 200000,  4,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.run( 600000,  4,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,800000)
>>> exit())
```

See the [results](scenario_tracking_frequency_step.md)

### Linear frequency change

We start the simulation with 0 phase and frequency error. Run the loop for 1 second. After that we start a triangular frequency sweep. Very long simulation.
BL=4Hz, fsw=0.025Hz, Asw=1000Hz : this means 100 Hz/sec linear frequency change

```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 20, 0.025, 1000,  40.0,    2e7, lfsel=1, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 200000,  4,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.run( 8200000,  4,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=1, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,8400000)
>>> exit()
```

See the [results](scenario_lin_frequency_change.md)

### Sine wave frequency change

We start the simulation with 0 phase and frequency error. Run the loop for 1 second. After that we start a sine wave frequency sweep. Very long simulation.  Due to the non zero frequency acceleration the phase error is not zero. It is a sine wave too. The maximum phase error is inversely proportiona with the third power of BL.
BL=4Hz, fsw=0.025Hz, Asw=1000Hz.

```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 20, 0.025, 1000,  40.0,    2e7, lfsel=1, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 200000,  4,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=1, pdmode=0, openloop=0)
>>> test.run( 8200000,  4,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=1, swmode=1, pdmode=0, openloop=0)
>>> test.show_loop(0,8400000)
>>> exit()
```

BL=8Hz, fsw=0.025Hz, Asw=1000Hz. Here the maximumphase error is 1/8 of the phase error in the previous simulation.

```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 20, 0.025, 1000,  40.0,    2e7, lfsel=1, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 200000,  8,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=1, pdmode=0, openloop=0)
>>> test.run( 8200000,  4,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=1, swmode=1, pdmode=0, openloop=0)
>>> test.show_loop(0,8400000)
>>> exit()
```

See the [results](scenario_sin_frequency_change.md)

### Loop bandwidth change

```python
> python3
>>> from tests import *
>>> test=PLLorCostasLoopTest(2e5, 1e4, 1250, 5000, 300, 0.025, 1000,  40.0,    2e7, lfsel=0, pdtype=0, bmode=1, mu=1.0)
>>> test.run( 600000,  4.0,  65.6, M=1, dteta=0.0, df=8.0, Gs=1.0, Gn=31.6, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.run( 400000,  0.125,  65.6, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=31.6, startsw=0, swmode=0, pdmode=0, openloop=0)
>>> test.show_loop(0,1000000)
>>> exit()
```

See the [results](scenario_loop_bandwidth_change.md)


Go back to the [tests page](tests.md)\
Go back to the [start page](../README.md)
