![sdrflow logo](images/pllpy_logo.svg  "sdrflow")

**The author:** ***Dr.Janos Selmeczi, HA5FT***. You could reach me at <ha5ft.jani@freemail.hu>
***
## Getting started with pllpy
### Platform requirements
You should have an ***64 bit Linux system*** to try the pllpy workbench. I used the ***ubuntu 20.04.4*** system for the development. You need ***python3, numpy, scipy and matplotlib.pyplot***. I used ***python3 v3.8.10***, ***numpy v1.17.4***, ***scipy v1.3.3*** and ***matplotlib v3.1.2*** for the development. All this modules came from the ubuntu repository.
### Install prerequisites
You should have some packeges installed on your system. On ubuntu or debian systems you could installed them using the folowing command:

```bash
> sudo  apt-get install git python3 python3-numpy python3-scipy python3-matplotlib
```

### Get the workbench source

You could download the software using git.

```bash
# create a work directory if you does not have one
> mkdir work
# go to this directory
> cd work
# clone the pllpy repository from github
> git clone https://github.com/ha5ft/pllpy
```
you will get the workbench in the ***pllpy*** directory in your working directory.

### Start experimenting with the workbench

You do not have to compile the source code. You just have to use the python3 console to test parts of the code. The code provides three test classes. You should read the instruction for running those test
* [Try the included PLL test](try_PLLorCostasLoopTest.md)
* [Try the included FLL test](try_FrequencyLockedLoopTest.md)
* [Try the included Bit Recovery Loop test](try_BitRecoveryLoopTest.md)

or you could go back to the [start page](../README.md)
