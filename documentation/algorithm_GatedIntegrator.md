![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# GatedIntegrator class

```python
# This class implements a special integrator.
# It is realy a moving average filter which has not fixed averaging length. The
# number of samples to be averaged is determined by an external clock.
# The result is the sum of the input samples divided by the number of samples added
# The clock pulse is 1 sample long.
# A clock pulse:
#   - ends the current integration cycle
#   - starts the next integration cycle
# The integrator puts the result to the output when a cycle has been finished and
# outputs the same valu until the next result wiil be available.
# The integrator signals a new result by a one sample long pulse on the valid output.
# The integrator also outputs of the length of the integration cycle.
# If each of the cyle is N sample long then the integrator does a decimation by N.
#
# Algorithm:
#   The clock signal: c
#   Initialization in the constructor:
#       k = 1
#       u = 0.0
#       y = 0.0
#       N = 0
#       v = 0
#
#   In the next() function:       
#       The output:
#               | u/k    if c==1
#           y = |
#               | y      otherwise
#
#               | k      if c==1
#           N = |
#               | N      otherwise
#
#           v = c
#
#       The integration:
#               | x      if c==1
#           u = |
#               | u+x    otherwise
#
#       The integration length count signal:
#               | 1,     if c==1
#           k = |
#               | k+1,   otherwise
#
#   The diagramm of the algorithm:
#
#  x[n]     |---|            u[n]               |---| v[n]|----|     |-----|
# >>------->| + |>----------------------------->| / |>--->|z^-1|>--->|- \  |
#           |---|                          |    |---|     |----|     |1  \ |
#             |                            |      |                  |    -|>---------------------->>  y[n]
#             |      |-----|     |----|    |      |                  |0    |                |
#             |      |  / -|<---<|z^-1|<--<|      |              |>->|-    |                |
#             |      | /  0|     |----|           |              |   |-----|                |
#             |<----<|-    |                      |              |      |        |----|     |
#                    |    1|                      |              |<-----)-------<|z^-1|<---<|
#                    |    -|<--- 0.0              |                     |        |----|
#                    |-----|                      |                    c[n]
#                       |                         |
#                       |                         |
#                      c[n]                       |
#           |---|               k[n]              |       |----|     |-----|
#     1 >-->| + |>------------------------------->|>----->|z^-1|>--->|- \  |
#           |---|                          |              |----|     |1  \ |
#             |                            |                         |    -|>---------------------->> N[n]
#             |      |-----|     |----|    |                         |0    |                |
#             |      |  / -|<--->|z^-1|<--<|                     |>->|-    |                |
#             |      | /  0|     |----|                          |   |-----|                |
#             |<----<|-    |                                     |      |        |----|     |
#                    |    1|                                     |<-----)-------<|z^-1|<---<|
#                    |    -|<--- 0                                      |        |----|
#                    |-----|                                           c[n]
#                       |                          
#                       |                          
#    c[n]              c[n]                        
# >>---->                                                                                  c[n] >-->> v[n]
#
# Constructor
#   integrator=GatedIntegrator()
# Parameter passed to the constructor:
#   None
# next method
#   integrator.next(c, x)
#   Parameters:
#       x : Current sample of the inputs signal
#       c : Current sample of the clock signal
#   Return:
#       y : the last valid output sample
#       N : the length of the last valid integration in samples
#       v : valid flag
#           v == 1 means new valid output sample is provided
#                  on the y output
#           v == 0 means the last valid output sample is
#                  provided on y outputs
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
