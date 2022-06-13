![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# MovingAverageFilter class

```python
# This class implement a fixed length real moving average filter.
# Algorithm:
#   The mathematical model of the filter:
#       N is the length of the filter
#       x[k] = 0.0 for -N < k <= -1
#
#                      N-1
#       y[n] = (1/N) * Sum(x(n-i))
#                      i=0
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
#     x[n]            + |---|      +|---|        u[n]     |---|    y[n]
# >>------------------->| + |>----->| + |---------------->| * |>-------->>
#          |            |---|       |---|            |    |---|
#          |              |-          |+             |      |
#          |              |           |              |      |
#          |    |----|    |           |    |----|    |     1/N
#          |>---|z^-N|--->|           |<---|z^-1|<---|
#               |----|  x[n-N]      u[n-1] |----|
#
#   The delay is implemented by using a circular buffer.
#
#   In the constructor:
#       i = 0
#       s[k] = 0.0  for 0 <= k < N
#   In the run function:
#       u += x - s[i]
#       y = s / N
#       s[i] = x
#           | 0    if i==N-1
#       i = |
#           | i+1  otherwise
# Constructor
#   filter=MovingAverageFilter(N)
#   Parameter passed to the constructor:
#       N : The length of the filter
# next method
#   filter.next(x)
#   Parameter:
#       x : current sample of the input signal
#   Return:
#       y : Currenr sample of the output signal
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
