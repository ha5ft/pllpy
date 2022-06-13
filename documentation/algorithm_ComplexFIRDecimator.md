![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# ComplexFIRDecimator class

```python
# This class implements A FIR decimator
# Algorithm:
#   Mathematical model:
#                    M-1
#           Iy[k] =  Sum(Ix[k*D-M+j]*b[j])
#                    j=0
#
#                    M-1
#           Qy[k] =  Sum(Qx[k*D-M+j]*b[j])
#                    j=0
#       An equivalent impulse response to an 3 stages D=5 CIC decimator filter
#
#           b = [1.0, 3.0, 6.0, 10.0, 15.0, 18.0, 19.0, 18.0, 15.0, 10.0, 6.0,\
#                3.0, 1.0] / 125.0
#
#   Implementation:
#
#   The operation is cyclic, the cycle length is D.
#   In the first D-1 invocation of next()
#       - receive an input I,Q sample pair and stores it in the shift register
#       - outputs the last valid output sample
#       - outputs a zero value valid flag
#   In the Dth invocation of next()
#       - receive an input I,Q sample pair and stores it in the shift register
#       - does the convolution
#       - shits the shift registers by D
#       - outputs the result of the convolution
#       - outputs a 1 value valid flag
#   D : decimation factor
#   M : Impulese response length
#   b : impulse response
#   Is, Qs : shift registers of length M
#   In the constructor:
#       Is[j] = 0.0 for 0 <= j < M
#       Is[j] = 0.0 for 0 <= j < M
#       l = M - D
#       c = 0
#       Iy = 0.0
#       Qy = 0.0
#       v = 0
#   In the next() function
#       Is[l] = Ix
#       Qs[l] = Qx
#       in case of c < D-1
#           c += 1
#           l += 1
#           v = 0
#           output Iy,Qy,v
#       in case of c == D-1
#                 M-1
#           Iy =  Sum(Is[i]*b[i])
#                 i=0
#
#                 M-1
#           Qy =  Sum(Is[i]*b[i])
#                 i=0
#
#           in case of M > D
#               Is[j] = Is[j+D]  for 0 <= j < M-D
#               Is[j] = Is[j+D]  for 0 <= j < M-D
#
#           v = 1
#           l = M-D
#           c = 0
#           output Iy,Qy,v
#
#
#   The diagramm of the algorithm:
#
#   Ix[kN+j] 0 <= j < D
# >>-------------------------------------------------
#                                                   |
#                                                   |M-D+j
#                    0 |                  M-D|      |       | M-1
#                     |--------------------------------------|
#             Is[]    |     shift register of length M       |
#                     |--------------------------------------|
#                      |                                    | Iu[N-1]
#                      |                                    |
#                      |                                    |
#                     |--------------------------------------|
#                     |           M-1                        |
#                     |           Sum(b[i]*Is[i])            |--------------->>
#                     |           i=0                        |         Iy[k]
#                     |--------------------------------------|          
#                      |                                    |
#                b[0]  |                                    | b[M-1]
#                      |                                    |
#                     |--------------------------------------|
#                     |              b[i] 0 <= i < M         |
#                     |--------------------------------------|
#                      |                                   |
#                b[0]  |                                   | b[M-1]
#                      |                                   |
#                     |--------------------------------------|
#                     |           M-1                        |
#                     |           Sum(b[i]*Qs[i])            |--------------->>
#                     |           i=0                        |         Iy[k]
#                     |--------------------------------------|          
#                      |                                    |
#                Qs[0] |                                    | Qs[M-1]
#                      |                                    |
#                     |--------------------------------------|
#                     |    shift register of length M        |
#                     |--------------------------------------|
#                   0 |                   M-D|      |       | M-1
#                                                   |M-D+j
#                                                   |
# >>-------------------------------------------------
#   Qx[kN+j] 0 <= j < D
#
# Constructor
#   decimator=ComplexFIRDecimator(b, D)
#   Parameters passed to the constructor:
#       D   : Decimation factor
#       b[] : Impulse response
#             It should be time reversed
#             It should be an numpy array with numpy float type of elements
#             The legth of the array defines M
# next method
#   decimator.next(Ix, Qx)
#   Parameters:
#       Ix : current sample of the inphase input
#       Qx : current sample of the quadrature phase input
#   Return:
#       Iy : current sample of the inphase output
#       Qy : current sample of the quadrature phase output
#       v  : valid flag
#           v == 1 means new valid output sample is provided
#                  in the Iy, Qy output
#           v == 0 means the last valid output sample is
#                  provided on Iy, Qy outputs
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
