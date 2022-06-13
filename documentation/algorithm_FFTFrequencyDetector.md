![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# GatedIntegrator class

```python
# This class implements a frequency detector using FFT
# The detector suppose that there is only a single signal in the frequency range
# visible for the detector. The detector determine the signal's frequency and
# power and the average noise power in a bin.
# The frequency range is [-fs/2,fs/2], where fs is the sampling frequency.
# The next() function is invoked for every new input sample.
# Until the function has reseived the necessary number of samples for the FFT
# it only stores the received sample and outputs the last results with a 0 value
# valid flag.
# After the necessary number of samples have been received the function processes
# the samples and outputs the new result with a 1 value valid flag.
# Algorithm:
#   Computing the FFT of the input samples
#   Computing the power spectrum
#   Finding the bin with maximum power
#   Correcting the bin power using the power of the two neighbor bins
#   Computing the total power by summing a bin powers
#   If the maximum bin power is above the threshold we say that there is carrier in
#       the frequency range
#   Correcting the total power using the powers of the signal.
#
# Constructor
#   detector=FFTFrequencyDetector(fs, N, threshold)
#   Parameters:
#       fs : samplink frequency
#       N   : FFT width
#       threshold  : factor for detecting carrier presece
#                    The carrier is a valid carrier if
#                       ps > threshold
#                    where ps is the signal power
# next method
#   detector.next(Ix, Qx)
#   Parameters:
#       Ix : current sample of the inphase input signal
#       Qx : current sample of the quadrature phase input signal
#   Return:
#       v    : valid flag
#            0 : last valid result
#            1 : new result
#       f    : frequency of the carrier
#       nbin : the index of the bin the carrier is in
#              positive index for positive frequency
#              negative index for negative frequency
#       pn   : average noise power of a bin
#       ps   : signal power (max bin power)
#       cvalid   : carrier presence flag
```

back to the [algorithms page](algorithms.md)\
back to the [start page](../README.md)
