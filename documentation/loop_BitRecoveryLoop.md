![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# BitRecoveryLoop class

```python
# This class implements a PLL and a Costas Loop.
# This class uses classes from algorithms.py
# This loop behaves as a phase locked loop or a costas loop depending on how
# its phase detector is working. The behavior of the phase detector is controlled
# by the pdtype parameter passed to the constructor of the class.
# The loop could synchronize to
#   - unmodulated carrier
#   - BPSK modulated carrier
#   - QPSK modulated carrier
# The loop has a Type2 or Type3 loop filter selectable at contruction time. It has
# a decimator. The dacimation factor could be selected at construction time. The
# bandwidth of the lowpass filter and the loop filter as well as the phase margine
# of the loop filter could be changed at runtime.
# The class implements a lock detector described in the algorithms.py file.
# In the current implementation the loop outputs the output signals of its NCO and
# the lock signal.
# The diagramm of the loop:
#
#                                       N                         
#                                       |                         
#                                       |                                 
#                                 |----------|vk[n]   
#      |----------|               |          |>--->   
# x[n] |  Complex |     y[n]      | Complex  |       
# >>-->|          |>------------->| FIR      |>----------------------------->|
#      |   Mixer  |               | Decimator|      u[k]                     |
#      |----------|               |----------|                               |
#            | z1[n]                                                         |
# z[n]       |                                                          u[k] |
# <<--------<|z[n]                                                           |
#            |                * Depending on lfsel                           |
#         A  |  teta0           Type2 or Type3                  pdmode       |
#         |  |  |               Filter is used                     |         |
# izc[n]|---------|               |----------|                |----------|   |
# <<---<| Complex |               |Type2 or 3|                |          |   |
#       |         |<-------------<|   Loop   |<--------------<|  Phase   |<-<|
# <<---<|   NCO   |<-|     c[k]   |  Filter* |      e[k]      | Detector |
# qzc[n]|---------|  |            |----------|                |----------|
#         |  |  |    |               |    |                        |
#        fs  |  w0  dteta=0.0        Kp   Ki                      pdtype
#            |
#         invertQ=1
#
# On the diagramm the signals x,y,u,v,z,s and mean[l] are complex signals represented
# by their inphase (real part) and quadrature phase (imaginary part) components.
# The sample rate for the various indexes:
#   n : fs
#   k : fs/N
#   l : fs/(N*Nl)
# The PLL could use Type2 or Type3 loop filter. The selected type of loopfilter
# will be instantiated in the constructor.
# Parameters passed to the constructor:
#   fs     : the input sampling frequency
#   w0     : nominal angular frequency of the carrier and the ComplexNCO
#   N      : decimation factor in the PLL loop
#   Nfir   : length of the fIR decimator
#   bw     : cutoff frequency of the FIR decinmator in Hz
#   lfsel  : loop filter selector
#            0 : type2
#            1 : type3
#            default is type3
#   pdtype : type of the phase detector
#            0 : Carrier
#            1 : BPSK
#            2 : QPSK
#            default is Carrier
#   mu     : safety factor for phase detector
#            default value: 1.0
#
# Inputs: (next())
#   Ix : current sample of the inphase input
#   Qx : current sample of the quadrature phase input
#   Kp : current value of the Kp parameter of the loop filter
#   Ki : current value of the Ki parameter of the loop filter
#   pdmode : phase detector mode selector
#            0 : linear mode
#            1 : pricipal value mode
#            default is the linear mode
#   openloop : open loop enable flag
#              0 : closed loop
#              1 : open loop
#
# Outputs: (next())
#   Iz      : current sample of the NCO inphase output
#   Qz      : current sample of the NCO quadrature phase output
#   Iu     : current sample of the inphase signal at the phase detector input
#   Qu     : current sample of the quadrature phase signal at the phase detector input
#   e      : current sample of the phase error
#   c      : current sample of the NCO control signal
#   vk     : valid flag for the Tu,Qu,e,c samples
```

back to the [loops page](loops.md)\
back to the [start page](../README.md)
