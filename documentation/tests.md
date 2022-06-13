![pllpy logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# Tests in the pllpy software

The pllpy workbench implemented several tests. There is test for each loop classes. The tests use an istance of the coresponding loop class an instances of the algorithm classes to provide input signals for the test. All the tests are in the ***test.py*** file. The workbench has the following tests:
* [PLLorCostasLoopTest](test_PLLorCostasLoopTest.md)
* [FrequencyLockedLoopTest](test_FrequencyLockedLoopTest.md)
* [BitRecoveryLoopTest](test_BitRecoveryLoopTest.md)

## PLLorCostaLoopTest
This class provide all the input signals you need to test a PLL or a costas loop and to inverstigate they dynamic behaviors and tracking capabilities.

## FrequencyLockedLoopTest
This class is a test bed for the instance of the FrequencyLockedLoop class. It provide a carrier with added white noise as an input signal. You can learn about the working of the frequency locked loop using this class.

## BitRecoveryLoopTest
This class runs the an instance of the BitRecovery class. It provide a baseband rectangular pulse PAM signal with optionally added noise as the input signal. You could investigat the non trivial working of the bit recovery loop using this test class.

# Unified test interface

All the tests has the same external interface.

## Contstructor

It is used to create an instance of the test class. usually several parameters should be passed to the constructor. Those are the parameters that can not be changed at runtime.

## run(samples,...) method

This is the work method of the test. It runs the test for a given number of input samples. It could provide run time changable parameters for the test.

## show(begin, end) methods

This methods graphically show you the various signals in the test class or in the used loop class. Somtimes the show method provide you statistics (mean and variance) of the signals. As an input you have to provide the sample range you are interested in. The range is for the loop input sample rate. If there are decimated signals, they are shown at the decimated rate. In the diagrams the x axis shows the sample index not the sample time.

Go back to the [start page](../README.md)
