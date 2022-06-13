![pllpy logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# Loops in the pllpy software

The pllpy workbench implemented several loops. The loops are composit algorithms which use instances of the algorithm classes.The loops have been implemented as python classes. All the loops are in the ***loops.py*** file. The workbench has the following loops:
* [PLLorCostasLoop](loop_PLLorCostasLoop.md)
* [FrequencyLockedLoop](loop_FrequencyLockedLoop.md)
* [BitRecoveryLoop](loop_BitRecoveryLoop.md)

## PLLorCostaLoop
This loop demonstraits the working of an ordinary PLL or a Costastas loop.
Depending on the mode parameter it is
* a PLL for unmodulated carrier
* Costas loop for BPSK modulated carrier
* Costas loop for QPSK modulated carrier
The mode parameter determine the working mode of the phase detector algorithm. All the other parts of the loop are not depend on this parameter. It is important that this loop is different from the usual formalizm of the Costas loops, but the working of the loop are the same. See more in the documentation of the [PhaseDetector](algorithm_PhaseDetector.md) class. You could learn the working and the behaviors of the PLL experimenting with this class.

## FrequencyLockedLoop
This loop demonstaits the ideas behind a frequency locked loop. It could be used to synchronize to the frequency of a signal which have distict spectral line at its frequency which is above the other spectral parts of the signal. The loop could be used together with a PLL to provide the initial frequency lock in a wide frequency range.

## BitRecoveryLoop
This loop implements one of the possible many bit recovery algorithms. Its input is a rectanguler PAM baseband signal. The algorithm synchronizes to the frequency and phase of the input signal.

# Unified loop interface

All the loops has the same external interface.

## Contstructor

It is used to create an instance of the loop class. usually several parameters should be passed to the constructor. Those are the parameters that can not be changed at runtime.

## next() method

This is the work method of the loop. All the loops are working on ***sample by sample bases***. They receives the current time samples of their input signals and output the current time samples of their output signals. In some cases the output signals are sampled with different sampling rate. In this case the output samples are associated with ***valid new output flags***. These flags are set when a new valid output sample has been outputed. So, there is no block processing in the algorithm. We choose this to reduce the processing delays.
The ***next() method*** could receive not only the imput samples but several parameters too. Those are the runtime changable parameters.

## reset() method

This method reinitialises all the instant variables to the values set by the constructor.

Go back to the [start page](../README.md)
