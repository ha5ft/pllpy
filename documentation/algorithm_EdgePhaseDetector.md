![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# EdgePhaseDetector class

```python
# This class implements a phase detector which detect the edge position of a
# square wave signal relative to the edge of the reference square wave
# signal.
# Algorithm:
#   The square wave signal
#            |  A, if      nT+dT <= t < n*T+T/2+dT
#     x(t) = |
#           | -A, if n*T+T/2+dT <= t < (n+1)*T+dT
#     T is the period of the signal
#     dT is the time shift relative to the reference signal
#    The reference signal is a square wave signal
#            |  B, if         mTr <= t < m*Tr + Tr/2
#     y(t) = |
#            | -B, if m*Tr + Tr/2 <= t < (m+1)*Tr
#     Tr is the period of the reference signal and Tr~T
#
#   The algoritm's'inputs:
#              (m+1)*Tr
#      Ix[m] : integral(x(t))
#               m*Tr
#
#              (m+1)Tr-Tr/2
#      Qx[m] : integral(x(t))
#              m*Tr-Tr/2
#
#    Computing the error:
#             | -sign(I[m-1])*Q if sign(I[m-1]) != sign(I[m])
#      e[n] = |
#             | 0 otherwise
#     IMPORTANT: If there is no tranzition edge between the two consecutive
#                interwall the error will bw zero.
#                To lessen the effect of noise in the implementation we set
#                certain signal to 0 if they are close to zero. We use this
#                in computing sign of the signal samples. The sign will be
#                0 if the sample is close to zero. This will result a
#                zero error instead of a wrong error.
# This algorithm works for other kind of signals, for example dirac delta
# trains shaped by rised cosine or root rised cosine filter.
#
# Constructor
#   detector=EdgePhaseDetector(threshold)
#   Parameters passed to the constructor:
#       threshold : threshold valu around zero
# next method
#   detector.next(I, Q)
#   Parameters:
#       I : current sample of the inphase integral
#       Q : current sample of the quadrature phase integral
# Output:
#   e  : phase error
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
