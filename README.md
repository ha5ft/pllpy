![pllpy logo](documentation/images/pllpy_logo.svg  "logo")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
## Important documents you should read
- [Getting started with pllpy](documentation/get_started.md)
- [Try the included PLL test](documentation/try_PLLorCostasLoopTest.md)
- [Try the included FLL test](documentation/try_FrequencyLockedLoopTest.md)
- [Try the included Bit Recovery Loop test](documentation/try_BitRecoveryLoopTest.md)
- [Algorithms](documentation/algorithms.md)
- [Loops](documentation/loops.md)
- [Tests](documentation/tests.md)

## What is pllpy?
Pllpy is a software workbench, which could be used to experiment with software defined PLLs. It has the following features:

* It has a number of algorithm classes for implementing various loops
  * ***generators*** classes for signal sources
  * ***mixers*** class for frequency translation
  * ***decimators*** classes for changing the sample rate
  * ***low pass filters*** classes to reduce signal bandwidth
  * ***phase and frequency detectors*** classes for detecting phase error
  * ***loop filters*** classes to define the dynamic behavior of the loops
  * ***lock and power detectors*** classes to define the dynamic behavior of the loops
* 3 loop classes have been implemented
  * ***PLLorCostasLoop*** a loop for recovering not modulated and BPSK or QPSK modulated carriers
  * ***FrequencyLockedLoop*** for locking for the frequency of a  not modulated carrier. It could be used together with a PLL loop for the initial frequency acquisition.
  * ***BitRecoveryLoop*** for bit synchronization for a rectangular PAM signal.
* Test classes for the above loops.
* It uses the following python infrastructure
  * python 3
  * numpy
  * scipy
  * matplotlib.pyplot

It is Open Source software. Licensed under GNU GPL 3 or any later version.
The current version has been tested on ***64 bit Ubuntu 20.04 sytem***

## The steps you should take
First of all you should download the test bench. You should follow the instructions in the [Getting started with pllpy](documentation/get_started.md) document.
Next you could try the PLLorCostasLoopTest class following the instructions in the [Try the PLLorCostasLoopTest class](documentation/try_PLLorCostasLoopTest.md) document. You could try the FrequencyLockedLoopTest class following the instructions in the [Try the FrequencyLockedLoopTest class](documentation/try_FrequencyLockedLoopTest.md) document. Finally you could try the BitRecoveryLoopTest class following the instructions in the [Try the BitRecoveryLoopTest class](documentation/try_BitRecoveryLoopTest.md) document.
If the test with the suggested parameters are working you should start experimenting with your own parameters.
All of the classes have detailed documentation in the source code.

Finally you could start to roll out your own test and loop classes following the patterns of the existing source code.

If you have any questions, you could contact me at <ha5ft.jani@freemail.hu>.
