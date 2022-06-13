![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# BPSKModulator class

```python
# This class generate a BPSK modulated signal from
#   - a complex sinusoid carrier and
#   - an PAM encoded binary signal
# The two level of the PAM pulse corresponding to the two logic bit level is
# -1.0 or 1.0. To generate such a signal we could use dirac pulses of
# -1.0 or 1.0 value for the bits and we use a pulse shaping filter. The simples
# pulse shaping filter has a rectangular impulse response. The mathematical
# background of the used algoritm is that changing the sign of a sine wave is
# equivalent to changing its phase by pi.
# Algorithm:
#   bpsk[n] = bit[n] * carrier[n]
#   bpsk[n] and carrier[n] are complex signals i IQ representation
# Constructor:
#   modulator=BPSKModulator()
#   Parameters passed to the constructor:
#       None
# next method:
#   modulator.next(Icarrier, Qcarrier, pam)
#   Parameters:
#       Icarrier : inphase sample of the carrier.
#                  It should be sample of a sine wave.
#       Qcarrier : quadrature phase sample of the carrier.
#                  It should be sample of a sine wave.
#       pam      : sample of the PAM signal.
#                  It should have -1.0 or 1.0 value.
#   Return:
#       bpsk : sample of the BPSK modulated signal.
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
