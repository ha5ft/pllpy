![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# PhaseDetector class

```python
# This class implements a linear phase detector
# It compute the angle of the input complex number.
# The input is represented by a inphase and quadrature phase real number pair.
# The standard mathematical functions are able to cumpute only the principal
# value of the angle. If we suppose, that the input is a rotating complex
# vector we could unwrap the wrapping occurs in the output of the standard
# mathematical functions at pi or -pi. So we could determine the tru linear
# phase difference.
# Algorithm
#   The input of the detector is a complex vector in IQ reprezentation
#       Ix[n] = cos(phi[n])
#       Qx[n] = sin(phi[n])
#   We create a complex variable
#           s = Ix + jQx
#   Depending on the the detector type (pdtype) we do the following
#       In case of carrier mode:
#           s = s
#       In case of BPSK mode
#           s = s * s
#           So we have 2 * phi angle, this maps the BPSK modulation angels
#           to 0
#       In case of QPSK mode
#           s = s * s * s * s
#           so we have 4 * phi angle, which maps all the QPSK modulation
#           angel to pi
#   We compute the principal value of the phase angle
#       phi = phase(s)
#   The detector algorithm unwrap the wraping of the angle of the complex vector
#   by observing the angle difference of two consecutive samples
#
#          pi |   /|                  pi |    |\
#             |  + |x[n-1]               |    | +x[n]
#             | /  |                     |    |  \
#             |/   |     /              \|    |   \
#       ------/----|----/-----     ------\----|----\-----
#            /|    |   /                 |\   |     \
#             |    |  /                  | \  |
#             |    | + x[n]              |  + |x[n-1]
#         -pi |    |/                -pi |   \|
#
#
#           increasing angle          decreasing angle
#           wrapping at pi            wrapping at -pi
#           phi[n-1]>pi/2             phi[n-1]<-pi/2
#           phi[n]<-pi/2              phi[n]>pi/2
#           phi[n]-phi[n-1]<-mu*pi    phi[n]-phi[n-1]>mu*pi
#           unwrap+=2*pi              unwrap-=2*pi
#   Above mu>=1 is used to lower the false unwrapping du to noise. It depends
#   on the ration of the sampling rate and the signal frequency. In case of
#   low signal frequency the change of the principal value ot wrapping is
#   close to 2*pi, in case of high signal frequency the change of the principal
#   value is close to pi.
#   The return value:
#       in principal angle mode:
#           phi, in case of carrier type detector
#           phi / 2, in case of BPSK type detector
#           phi / 4 - pi / 4, in case of QPSK type detector
#               here we have to substract pi/4 becouse the modulation angles
#               have been maped to pi.
#       in linear mode: phi + unwrap
#   The algorithm works only if the sampling rate is higher than the critical
#   sampling rate.
# Constructor
#   detector=PhaseDetector(pdtype=0, mu=1.0)
#   Parameter passed to the constructor:
#       pdtype:  type of the phase detector
#                0 : Carrier
#                1 : BPSK
#                2 : QPSK
#                default is Carrier
#       mu : safety factor
#            default value: 1.0
# next method
#   detector.next(Ix, Qx, pdmode=0)
#   Parameters:
#       Ix   :  current sample of the inphase input signal
#       Qx   :  current sample of the quadrature phase input signal
#       pdmode :mode of angle computation
#               0 - principal value mode
#               1 - linear mode
#               The default value is 0.
#               It is a named parameter.
#   Return:
#       e  : phase of the input complex vextor
# In case of very low SNR the principal value mode is preferred.
```

back to the [algorithms page](algorithms.md)\
back to the [loops page](loops.md)\
back to the [start page](../README.md)
