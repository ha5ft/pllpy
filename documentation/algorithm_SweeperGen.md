![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# SweepGen class

```python
# The instances of this class following each start trigger generates a single
# periode of a triangle or a sine wave real signal. The DC value of the signals
# is 0.0.
# The amplitude and the angular frequency of the signal is controlled by
# parameters passed to the constructor.
# We use the output of this sweep generator as the control signal of an NCO to
# test the frequency tracking capabilities of a PLL loop.
# Algoritm:
#   Until a start trigger it outputs 0.0
#   At the start trigger it smples the mode input and uses the sampled value
#   for the current wavwform.
#   After receiving a start trigger it starts outputing the waveform determined
#   by the mode.
#   During the waveform output it ignores the start trigger
#   After one perode id finished it continues to output 0.0 and looks for other
#   start trigger.
#   During the waveform output:
#       phi[0] = 0
#       alpha  = phi[n-1] + w0/fs
#                 / 0       if alpha >= 2*pi
#       phi[n] = |  
#                 \ phi[n]  otherwise
#
#       In case of triangle wave mode:
#
#                  | A * (phi[n]/(pi/2))        if 0 <= phi[n] < pi/2
#           y[n] = | A * ((pi-phi[n])/(pi/2))   if pi/2 <= phi[n] < 3*pi/2
#                  | A * ((phi[n]-2*pi)/(pi/2)) if 3*pi/2 <= phi[n] < 2*pi
#
#       In case of sine wave mode:
#
#           y[n] = A * sin(phi[n])
#
# Constructor:
#   sweeper=SweeperGen(fs, w0, A)
#   Parameters passed to the constructor:
#       fs : sampling frequency
#       w0 : angular frequency of the waveform
#       A  : amplitude of the waveform
# next method
#   sweeper.(start=0, mode=0)
#   Parameters:
#    start : start trigger
#            a single sample of value 1 starts a new period
#            if a sweep is in progress no new sweep wil be started
#    mode  : waveform type
#            0 : triangle wave
#            1 : sine wave
#  Return:
#   y   : output waveform sample
#   end : a single sample of value 1 outputed at the end of the period
```

back to the [algorithms page](algorithms.md)
back to the [start page](../README.md)
