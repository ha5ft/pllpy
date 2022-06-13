![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# ComplexVariableMovingAverageFilter class

```python
# This class implement a variable length moving average filter.
# The length of the filter is limited by the maximum length parameter set in
# the constructor.
# For every invocation the filter provides new output samples, so no decimation.
# Algorithm:
#   The mathematical model of the filter:
#       M <= Mmax
#       Ix[k] = 0.0 for -Mmax < k <= -1
#       Qx[k] = 0.0 for -Mmax < k <= -1
#
#                       M-1
#       Iy[n] = (1/M) * Sum(Ix(n-i))
#                       i=0
#
#                       M-1
#       Qy[n] = (1/M) * Sum(Qx(n-i))
#                       i=0
#
#   The cutoff frequency of the moving average filter:
#       fcutoff ~ fs*0,443/M
#       fs : sampling frequency
#
#   In the implementation:
#       Mmax is the limit of the filter length
#       M is the actual filter length
#
#       Initialization in the constructor:
#           Is[k] = 0.0 for 0 <= k < Mmax
#           Qs[k] = 0.0 for 0 <= k < Mmax
#       In the nex() function:
#           Is[k] = Is[k-1] for k=Mmax-1 ....k=1
#           Qs[k] = Qs[k-1] for k=Mmax-1 ....k=1
#
#           Is[0] = Ix
#           Qs[0] = Qx
#
#                        M-1
#           Iy = (1/M) * Sum Is[i]
#                        i=0
#
#                        M-1
#           Qy = (1/M) * Sum Qs[i]
#                        i=0
# Constructor:ComplexVariableMovingAverageFilter
#   filter=(Mmax)
#   Parameter passed to the constructor:
#       Mmax : Maximum length of the filter
# next method
#   filter.next(Ix, Qx, M)
#   Parameters:
#       M    : current filter length
#              0 < M <= Mmax
#       Ix   : Current sample of the input inphase signal
#       Qx   : Current sample of the input quadrature phase signal
#   Return:
#       Iy   : Current sample of the output inphase signal
#       Qy   : Current sample of the output quadrature phase signal
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
