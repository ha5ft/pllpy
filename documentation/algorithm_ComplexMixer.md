![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# ComplexMixer class

```python
# This class implements a complex mixer
# The input and output signals are in the inphase and quarature phase form.
# Algorithm:
#   Iy = Is * Ilo - Qs * Qlo
#   Qy = Is * Qlo + Qs * Ilo
# Constructor:
#   mixer=ComplexMixer()
#   Parameters passed to the constructor:
#       None
# next method
#   mixer.(Is, Qs, Ilo, Qlo)
#   Parameters:
#       Is  : Inphase sample of the input signal
#       Qs  : Quadrature phase sample of the input signal
#       Ilo : Inphase sample of the local oscillator signal
#       Qlo : Quadrature phase sample of the local oscillator signal
#   Return:
#       Iy  : Inphase output sample
#       Qy  : Quadrature phase output sample
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
