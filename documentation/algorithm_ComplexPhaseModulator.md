![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# ComplexPhaseModulator class

```python
# This class phase modulates a complex sinusoid carrier.
# The modulating signal is a real signal and it represents the modulating phase.
# The complex carries is represenred by the corresponding inphase and quadrature
# phase real signals. The complex output is represented too with the corresponding
# inphase and quadrature phase real signals.
# Algorithm:
#   Iy[n] = Re(exp(j*phi[n])) * Ix[n] - Im(exp(j*phi[n])) * Qx[n]
#   Qy[n] = Re(exp(j*phi[n])) * Qx[n] + Im(exp(j*phi[n])) * Ix[n]
# Constructor:
#   modulator=ComplexPhaseModulator()
#   Parameters passed to the constructor:
#       None
# next method:
#   modulator.next(Ix, Qx, phi)
#   Parameters:
#       Ix : inphase sample of the carrier.
#            It should be sample of a sine wave.
#       Qx : quadrature phase sample of the carrier.
#            It should be sample of a sine wave.
#       phi: sample of the modulating phase.
#   Return:
#       Iy : inphase sample of the modulated signal.
#       Qy : quadrature phase sample of the modulated signal.
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
