![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# ComplexChannel class

```python
# This class a noisy channel for  complex signal.
# The noise is Gaussian noise.
# Using the channel one could control the signal to noise ratio.
# The input and output signals are complex signal in the inphase and quadrature
# representation.
# Algorithm:
#   Iy[n] = Gs * Isl[n] + Gn * In[n]
#   Qy[n] = Gs * Qsl[n] + Gn * Qn[n]
# Constructor:
#   channel=ComplexChannel8
#   Parameters passed to the constructor:
#       None
# next method
#   channel.next(Ix, Qx, In, Qn, Gs, Gn)
#   Parameters:
#       Ix : Input inphase signal sample
#       Qx : Input quadrature phase signal sample
#       In : Input inphase noise sample
#       Qn : Input quadrature phase noise sample
#       Gs : Signal gain
#       Gn : Noise gain
#            If the variance of the input signal and the input noise is 1 then
#            SNR = 10*log(Gs^2/Gn^2)
#   Return:
#       Iy : Output inphase sample
#       Qy : Output quadrature sample
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
