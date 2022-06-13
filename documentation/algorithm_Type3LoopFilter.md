![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
# GatedIntegrator class

```python
# This class implements an type3 loop filter with parameter changing
# The code managing the parameter change ensures a smooth output at paramater
# changing.
# Mathematical modell:
#   Parameters:
#       BL    : SSB noise bandwidth
#       phiPM : Phase margin (in degree)
#       M     : Averiging length for computing the average value of the output
#   Derived parameters:
#       Kp    : the gain of the proportional part
#       Ki    : the gain of the integrator part
#   Parameter calculation:
#       ro = tg((phiPM+90°)/2)
#       Kp = 4*BL*((2*ro-1)/(2*ro+3))
#       w0 = K/ro
#       Ki = w0*Ts
#   The transfer function of the filter using Backward Euler integration
#
#       L(z) = (1+Ki/(1-1/z))*(1+Ki/(1-1/z))*Kp
#
#   Loop parameter change:
#       Ki could be change without any transient effect
#       If we change Kp we would have a jump at the output eve the input is unchanged.
#       So instead of using K as a
#           y = K * x
#       transfer function and always use a transfer curve going through the origo
#       we locally change the slope of the transfer curve as shown in the following
#       diagramm
#
#         y  ^
#            |                 /          .
#            |             K0 /         . K1
#            |               /        .
#            |              /       .
#            |             /      .
#            |            /     .
#            |__________ /___ .
#         y1 |          /   .|
#            |         /  .  |
#            |        / .    |
#         y0 |_______/.______|_____________
#            |     ./|       |
#            |   . / |       |
#            | .  /  |       |
#            |   /   |       |
#            |  /    |       |
#            | /     |       |
#            |/      |       |
#            --------|-------|-----------------> x
#                    x0      x1
#
#            y[n] = K0 * x[n]  n <= n0
#            y[n0+i] = K1 * (x[n0+i] - x[n0]) + y[n0]  i>0
#       The reason this works is that we change the parameters after the PLL achived lock
#       in order to reduce the PLL bandwidth and reduce the noise.
#       In the implementation we use average(y) in place of y[n] becouse when x is noisy
#       we still conuld have some jump at the output if we use y[n].
#
#   Implementation
#
#       ei[-1] = 0.0
#       Kp[-1] = 0.0
#       ccorr[-1] = 0.0
#       c[i] = 0.0   if -M<=i<0
#       cav[-1] = 0.0
#       ci[-1] = 0.0
#
#       ei[n] = ei[n-1] + Ki[n]*e[n]
#       cr[n] = ei[n] + e[n]
#
#
#              | 1    if Kp[n] != Kp[n-1]
#       s[n] = |
#              | 0    otherwise
#
#                  | cav[n-1] - Kp[n] * cr[n-1]       if s[n] = 1
#       ccorr[n] = |
#                  | ccorr[n-1]                       otherwise
#
#       c[n] = Kp * cr[n] + ccorr[n]
#
#       ci[n] = ci[n-1] + c[n] - c[n-M]
#       cav[n] = ci[n] / M
#       
#    The block diagramm of the filter with Backward Euler integration:
#
#                     Ki
#                     |
#   e[n]            |---|     |---|              ei1[n]      |---|
#   >>------------->| X |>--->| + |>------------------------>| + |>-->|
#     |             |---|     |---|                |         |---|    |
#     |                         |                  |           |      |
#     |                         |       |---|      |           |      |
#     |                         |<-----<| -1|<----<|           |      |
#     |                          ei[n-1]|z  |                  |      |
#     |                                 |---|                  |      |
#     |                                                        |      |
#     |>------------------------------------------------------>|      |
#                                    e[n]                             |
#                                                                     |
#                                                                     |
#                                                                     |
#     |<-------------------------------------------------------------<|
#     |
#     |                             
#     |               Ki
#     |               |
#     |             |---|     |---|              ei2[n]      |---|
#     |>----------->| X |>--->| + |>------------------------>| + |>-->|
#     |             |---|     |---|                |         |---|    |
#     |                         |                  |           |      |
#     |                         |       |---|      |           |      |
#     |                         |<-----<| -1|<----<|           |      |
#     |                         ei2[n-1]|z  |                  |      |
#     |                                 |---|                  |      |
#     |                                                        |      |
#     |>------------------------------------------------------>|      |
#                                                                     |
#                                                                     |
#                                                                     |
#                                    cr[n]                            |
#     |<-------------------------------------------------------------<|
#     |
#     |                             
#     |      c[n]     + |---|      +|---|       ci[n]     |---| cav[n]
#     |    |>---------->| + |>----->| + |---------------->| * |>------>
#     |    |            |---|       |---|            |    |---|
#     |    |              |-          |+             |      |
#     |    |              |           |              |      |
#     |    |    |----|    |           |    |----|    |     1/N
#     |    |>---|z^-M|--->|           |<---|z^-1|<---|
#     |    |    |----| cn[n-N]     ci[n-1] |----|
#     |    |
#     |    |<----------------------------------------------------<|
#     |                                                           |
#     |                                                           |
#     |    cr[n]     |---|                             |---|      |      c[n]
#     |>------------>| X |>--------------------------->| + |>---->|>----------------->>
#     |              |---|                             |---|      
#     |                |                                 |
#     |                |    |----| K[n-1]|---|           |
#     |                |>---|z^-1|>----->|   |           |
#     |   Kp[n]        |    |----|       |   | s[n]      |
#  >>-)--------------->|                 |!= |>-->|      |
#     |                |      Kp[n]      |   |    |      |
#     |                |>--------------->|   |    |      |
#     |                |                 |---|    |      |
#     |                |                          |      | ccorr[n]
#     | cr[n]          |                          |      |
#     |                |                          |      |  
#     |    |----|    |---|    |---|            |----|    |
#     |>-->|z^-1|>-->| X |>-->| + |>---------->|1\  |    |
#          |----|    |---|    |---|            |  \ |    |
#                               |              |   \|>-->|
#                               |              |    |    |
#   cav[n] |----|    cav[n-1]   |        |>--->|0   |    |
#   >----->|z^-1|>--------------|        |     |----|    |
#          |----|                        |               |
#                                        |     |----|    |
#                                        |<---<|z^-1|<--<|
#                                              |----|
#
# Constructor
#   filter=Type3LoopFilter(Ts, M)
#   Parameters passed to the constructor
#       Ts : sampling time
#       M  : length of the averaging
# next method
#   filter.next(e, Kp, Ki)
#   Parameters:
#       e  : current error sample
#       Kp : current Kp value
#       Ki : current Ki value
#   Return:
#       c  : current sample of the filter output                                         
#
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
