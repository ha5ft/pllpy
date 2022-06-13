![pllpy logo](images/pllpy_logo.svg  "pllpy")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
# Trying the FrequencyLockedLoopTest class

This class is a test bench for the FrequencyLockedLoop class.
In addition of the test instance of the FrequencyLockedLoop class it implements
* a carrier generator,
* a noise generator,
* channel simulator which adds noise to the carrier,
* signal visualization functions for checking signals in the test as well as in the loop class.

## Run your first simulation

You could run an instance of the class from the python3 consol. You shoud run it from a ***command terminal in the pllpy directory***.

```bash
# start the sdrflow runtime_application
> python3
Python 3.8.10 (default, Mar 15 2022, 12:22:08)
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
```
Now you are in python console and you could use it for the testing
```python
# You are in the python3 console
# First you have to import the code of the class
>>> from test import *
# You have to create an instance of the class
>>> test=FrequencyLockedLoopTest(2e5, 1e4, 4096,    2e7, threshold=1e-2)
>>> test.run(800000, 4.0,  65.6, df=1e4, Gs=1.0, Gn=0.0, openloop=0)
BL = 4.0
phiPM = 65.6
Kp = 11.007002311039455
Ki = 0.10225659266991194
.................................................................................................
n = 800000
l = 195
t = 4.0
exec time = 33.99460768699646
>>> test.show_loop(0,800000)
n rate sample range : 0 : 800000
l rate sample range : 1 : 195
```
At this point you should se the following 4 diagramms

![frequency error](results/fll_df1e4_N4096_BL4_phiPM65.6_phase_error.png "frequency error")
![bin index](results/fll_df1e4_N4096_BL4_phiPM65.6_bin_index.png "bin index")
![NCO control](results/fll_df1e4_N4096_BL4_phiPM65.6_NCO_control.png "NCO control")
![bin powers](results/fll_df1e4_N4096_BL4_phiPM65.6_bin_powers.png "bin powers")

After you close the diagramm windows

```python
>>> exit()
# You left the python3 consol and you are back in the linux shell
>
```

## Printouts on the python console

You have simulated an FLL locking process.
You simulated 4 second and the simulatur run 13.9 sec. So it is an order of magnitude slower than real time.
The run() function printed on the consol
* BL, the SSB noise bandwidth
* phiPM, the phase margin
* Kp, the computed proportional gain
* Ki, the computed integrator gain divided by the proportional gain
* a series of c characters showing that closed loop simulation is running without sweep
* n, the number of input samples used in the simulation
* l, the number of decimated sample at the frequency detector output
* t, the time length of the simulation
* exec time, the execution time of the simulation

## Parameters of the simulation

The important parameters of the simulation:
|Parameter|Value|Description|
|---------|-----|-----------|
|fs|2e5Hz|sampling rate|
|f0|1e4Hz|carrier frequency|
|Nfft|4096|decimation factor|
|BL|4Hz|SSB noise bandwidth|
|phiPM|65.6|phase margin|
|df|1e4Hz|frequency step|

Go back to the [get started page](get_started.md)\
Go back to the [start page](../README.md)
