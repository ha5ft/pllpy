![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# ComplexMovingAverageFilter class

```python
# This class implement a fixed length complex moving average filter.
# Algorithm:
#   The mathematical model of the filter:
#       N is the length of the filter
#       Ix[k] = 0.0 for -N < k <= -1
#       Qx[k] = 0.0 for -N < k <= -1
#
#                       N-1
#       Iy[n] = (1/N) * Sum(Ix(n-i))
#                       i=0
#
#                       N-1
#       Qy[n] = (1/N) * Sum(Qx(n-i))
#                       i=0
#
#   The cutoff frequency of the moving average filter:
#       fcutoff ~ fs*0,443/M
#       fs : sampling frequency
#
#   This filter could be used to compute the sliding mean of a signal. The mean is
#   using N signal sample.
#
# Implementation:
#   The diagramm of the implementation:
#
#    Ix[n]            + |---|      +|---|       Iu[n]     |---|   Iy[n]
# >>------------------->| + |>----->| + |---------------->| * |>-------->>
#          |            |---|       |---|            |    |---|
#          |              |-          |+             |      |
#          |              |           |              |      |
#          |    |----|    |           |    |----|    |     1/N
#          |>---|z^-N|--->|           |<---|z^-1|<---|
#               |----| Ix[n-N]     Iu[n-1] |----|
#
#
#    Qx[n]            + |---|      +|---|       Qu[n]     |---|   Qy[n]
# >>------------------->| + |>----->| + |---------------->| * |>-------->>
#          |            |---|       |---|            |    |---|
#          |              |-          |+             |      |
#          |              |           |              |      |
#          |    |----|    |           |    |----|    |     1/N
#          |>---|z^-N|--->|           |<---|z^-1|<---|
#               |----| Qx[n-N]     Qu[n-1] |----|
#
#   The delay is implemented by using a circular buffer.
#
#   In the constructor:
#       i = 0
#       Is[k] = 0.0  for 0 <= k < N
#       Qs[k] = 0.0  for 0 <= k < N
#   In the run function:
#       Iu += Ix - Is[i]
#       Qu += Qx - Qs[i]
#       Iy = Is / N
#       Qy = Qs / N
#       Is[i] = Ix
#       Qs[i] = Qx
#           | 0    if i==N-1
#       i = |
#           | i+1  otherwise
# Constructor
#   filter=ComplexMovingAverageFilter(N)
#   Parameter passed to the constructor:
#       N : The length of the filter
# next method
#   Parameters:
#       Ix : current sample of the input inphase signal
#       Qx : current sample of the input quadrature phase signal
#   Return:
#       Iy : Currenr sample of the output inphase signal
#       Qy : Currenr sample of the output quadrature phase signal
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
