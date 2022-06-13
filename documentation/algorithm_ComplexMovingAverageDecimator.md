![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# ComplexMovingAverageDecimator class

```python
# This class implement a decimator with a moving average filter.
# The input and output signals are comlex signal in inphase and quadrature phase
# representation.
# Algorithm:
#   It decimates by N
#   The mathematical model of the decimation:
#                       N-1
#       Iy[k] = (1/N) * Sum(I[k*N+i])
#                       i=0
#
#                       N-1
#       Qy[k] = (1/N) * Sum(Qx[k*N+i])
#                       i=0
#
#   The cutoff frequency of the moving average filter:
#       fcutoff ~ fs*0,443/N
#       fs : sampling frequency
#
#   In the implementation:
#       For N-1 invocation we update the running summs and output the last valid
#       output signal sample for the signal output and 0 for the valid flag output.
#       For the Nth invocation we finalize the running sum, divide them by N,
#       output the results on the signal output, output 1 on the valid flag output,
#       and clear the running sums and start a new cycle.
#       Initialization in the constructor:
#           i  = 0
#           Is = 0
#           Qs = 0
#           Iy = 0
#           Qy = 0
#           v  = 0
#       In the next() function
#
#                | Ix         if i == 0
#           Is = |
#                | Is + Ix    otherwise
#
#                | Qx         if i == 0
#           Qs = |
#                | Qs + Qx    otherwise
#
#                | Is/N       if i == N-1
#           Iy = |
#                | Iy         otherwise
#
#                | Qs/N       if i == N-1
#           Qy = |
#                | Qy         otherwise
#
#                | 1          if i == N-1
#           v  = |
#                | 0          otherwise
#
#                | 0          if i == N-1
#           i  = |
#                | i + 1      otherwise
#
# Constructor
#   decimator=ComplexMovingAverageDecimator(N)
#   Parameter passed to the constructor:
#       N : decimation factor
# next method
#   decimator.next(Ix, Qx)
#   Parameters:
#       Ix : current sample of the inphase input signal
#       Qx : current sample of the quadrature phase input signal
#   Return:
#       Iy : current sample of the inphase input signal
#       Qy : current sample of the quadrature phase input signal
#       v  : valid flag
#            v == 1 means new valid output sample is provided
#                   in the Iy, Qy output
#            v == 0 means the last valid output sample is
#                   provided on Iy, Qy outputs
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
