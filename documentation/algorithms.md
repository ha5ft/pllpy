![sdrflow logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# Algoritms in the pllpy software

The pllpy workbench has many of the algorithms you need to implement various PLLs. The algorithms have been implemented as python classes. All the algorithms are in the ***algorithms.py*** file. The algorithms belong to different categories.
* Signal generators
* Decimators
* Low pass filters and integrators
* Phase and frequency detectors
* Loop filters
* Lock and signal power detectors
* Various other algorithms

Each of the above categories contain several algorithms.

## Signal generators
These algorithms are signal sources. They are used in the loop and in the test classes. We have the followin generators:
* [ComplexNCO](algorithm_ComplexNCO.md)
* [ComplexNoiseGenerator](algorithm_ComplexNoiseGenerator.md)
* [RectPulsePAMGen](algorithm_RectPulsePAMGen.md)
* [SquareWaveNCO](algorithm_SquareWaveNCO.md)
* [SweeperGen](algorithm_SweeperGen.md)

## Decimators
These algorithms reduce the sample rate
* [ComplexFIRDecimator](algorithm_ComplexFIRDecimator.md)
* [ComplexMovingAverageDecimator](algorithm_ComplexMovingAverageDecimator.md)

## Low pass filters and integrators
These algorithms low pass filter they input signals. The gated integrator is a moving average filter which length is not fixed rather determined by an external clock signal.
* [ComplexMovingAverageFilter](algorithm_ComplexMovingAverageFilter.md)
* [ComplexVariableMovingAverageFilter](algorithm_ComplexVariableMovingAverageFilter.md)
* [MovingAverageFilter](algorithm_MovingAverageFilter.md)
* [GatedIntegrator](algorithm_GatedIntegrator.md)

## Phase and frequency detectors
These algorithms determine the phase or frequency of they input signal.
* [EdgePhaseDetector](algorithm_EdgePhaseDetector.md)
* [FFTFrequencyDetector](algorithm_FFTFrequencyDetector.md)
* [PhaseDetector](algorithm_PhaseDetector.md)

## Loop filters
Special filters. They determines the dynamic behavoirs and the noise bandwidth of the loops.
* [Type2LoopFilter](algorithm_Type2LoopFilter.md)
* [Type3LoopFilter](algorithm_Type3LoopFilter.md)

## Various other algorithms
* [BPSKModulator](algorithm_BPSKModulator.md)
* [ComplexChannel](algorithm_ComplexChannel.md)
* [ComplexMixer](algorithm_ComplexMixer.md)
* [ComplexPhaseModulator](algorithm_ComplexPhaseModulator.md)

# Unified algorithm interface

All the algorithms has the same external interface.

## Contstructor

It is used to create an instance of the algorithm class. usually several parameters should be passed to the constructor. Those are the parameters that can not be changed at runtime.

## next() method

This is the work method of the algorithm. All the algorithms are working on ***sample by sample bases***. They receives the current time samples of their input signals and output the current time samples of their output signals. In some cases, for example in case of decimators the algorithm, the algorithm needs samples for several consecutive sampling times to generate an output sample. In this case until it gets the necessary number of samples it outputs the last computed output samples and signals the output with a ***not valid new sample flag***. When it gets the last necessary samples it computes the new output and signals the new output with a ***valid new sample flag***. So there is no block processing in the algorithm. We choose this to reduce the processing delays.
The ***next() method*** could receive not only the imput samples but several parameters too. Those are the runtime changable parameters.

## reset() method

This method reinitialises all the instant variables to the values set by the constructor.

Go back to the [start page](../README.md)
