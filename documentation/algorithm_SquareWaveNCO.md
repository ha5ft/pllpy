![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# SquareWaveNCO class

```python
# This class generate square wave tone.
# In addition it generates an impulse at the rising edge
# and another impulse at the falling edge
# The algorithm:
#   phi[0] = teta0
#   alpha  = phi[n-1] + (w0+dw[n])/fs + dteta[n]
#
#                   / 1 if teta0 >= 0.0
#   y[-1] = y[0] = |
#                   \ 0 if teta0 < 0.0
#
#            / 1 if phi[n-1] < 2*pi  and alpha >= 2*pi
#           /  1 if phi[n-1] > -2*pi and alpha <= -2*pi
#   y[n] = |   y[n-1]
#           \  0 if phi[n-1] < pi  and alpha >= pi
#            \ 0 if phi[n-1] > -pi and alpha < -pi
#
#
#   izc[0] = 0
#
#
#             / 1, if y[n-1] < y and y[n] >= 0
#   izc[n] = |  
#             \ 0, otherwise
#
#
#   qzc[0] = 0
#
#
#             / 1, if y[n-1] < *0 and y[n] >= 0
#   qzc[n] = |  
#             \ 0, otherwise
#
#   if phi crosses 2*pi or -2*pi then it is decresed or increased by 2*pi
#
#             / alpha - 2*pi if alpha >= 2*pi
#   phi[n] = |  alpha        otherwise
#             \ alpha + 2*pi if alpha <= -2*pi
#
# Constructor:
#   nco=SquareWaveNCO(fs, w0, A, teta0=0.0, invertQ=0)
#   Parameters:
#       A           : absolute value of the output
#       fs          : sampling frequency
#       w0          : angular base frequency of the NCO
#       teta0       : starting phase of the NCO
#       invertQ     : If not zero, the sign of Qy will be changed to its oposite
#                     This will causes the frequency of the z signal to change
#                     to the -(w0 + dw) value. This is a permanent change.
# next method
#   nco.next(dw=0.0, dteta=0.0) :                     
#   Parameters:
#       dw       : angular frequency offset
#       dteta    : phase offset
#                  This input usually should be 0.0.
#                  A single nonzero sample will causes a phase jump in the
#                  output signal.
#                  A step function in this input will produce an angular
#                  frequency jump in the output. The jump value is dteta*fs.
#   Return:
#       Iy       : output inphase signal sample
#       Qy       : output quadrature phase signal sample
#       izc      : output zero crossing flag sample
#       qzc      : output pi crossing flag sample
```

back to the [algorithms page](algorithms.md)
back to the [start page](../README.md)
