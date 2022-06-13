# ******************************************************************************
# *                             Tests
# ******************************************************************************
# * File name:      tests.py
# * Platform:       ubuntu 20.04 64 bit
# *                 python3 v3.8.10
# *                 numpy v1.17.4
# *                 scipy v1.3.3
# *                 matplotlib v3.1.2
# *                 algorithms.py
# * Author:			Copyright (C) Dr. Selmeczi János, ha5ft: original version
# ******************************************************************************
# *                             Licencing
# ******************************************************************************
# *  This file is part of the pllpy project.
# *
# *  Pllpy is free software: you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  Pllpy is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with Pllpy.  If not, see <http://www.gnu.org/licenses/>.
# ******************************************************************************
# *							Revision history
# ******************************************************************************
# *	Author				Date		Comment
# ******************************************************************************
# *	Selmeczi János		12.06.2022	Original version
# *
# ******************************************************************************

# ==============================================================================
#                               Description
# ==============================================================================
# This file contains classes implementing various tests for the cslasses in the
# algorithms.py and loops.py files.
# All the tests have a common interface:
#   Constructor: It creates an instance of the used classes.
#                It receives the parameters for the test
#                It allocate the necessary resources for the test
#   run()      : It runs the test for a given number of signal sample with a set
#                of runtime changable parameters
#                It stores the signal samples (even for the signals internal to
#                the class being tested) for later presentation.
#                While running display the progress of the test.
#   show()     : Present the test results in numerical and/or graphical form for
#                a given sinal sample range.
# ==============================================================================

from algorithms import *
from loops import *
import time

# -----------------------------------------------------------------------------
#                              PLLorCostasLoopTest
# -----------------------------------------------------------------------------
# This class implements a test for the PLLorCostasLoop class in the loops.py file
# It uses the following classes from loops.py and algorithms.py:
#   - PLLorCostasLoop from loops.py which uses from algorithms.py the following classes
#       - ComplexMixer
#       - ComplexMovingAverageDecimator
#       - ComplexVariableMovingAverageFilter
#       - PhaseDetector
#       - Type2LoopFilter
#       - Type3LoopFilter
#       - ComplexNCO
#       - LockDetector
#   - SweeperGen from algoritms.py
#   - ComplexNCO from algoritms.py
#   - ComplexNoiseGen from algoritms.py
#   - ComplexChannel from algoritms.py
#   - SquareWaveNCO from algoritms.py
#   - RectPulsePAMGen from algoritms.py
#   - BPSKModulator from algoritms.py
# In the test the input signal to the PLL or Castas loop could be:
#   - a unmodulated carrier
#   - an unmodulated carrier with given level of white noise
#   - BPSK modulated carrier
#   - BPSK modulated carrier with given level of white noise
# For the carrier we could
#   - apply a given phase jump
#   - apply a frequency jump
#   - apply a triangle shape frequency sweep
#   - a sine function shape frequency sweep
# Instances of this class could be constructed for PLL test or Costas Loop test. The pdmode
# parameter passed to the constructor determines for wich test the instance will be generated.
# The PAM generator (SquareWaveNCO and RectPulsePAMGen) will operate only in the CostasLoop mode.
# Wheter the instance of the PLLorCostasLoop class operates in PLL or CostasLoop mode depend
# on the pdmode parameter passed to its constructor. This behavior can not be changed at
# runtime.
# At runtime we could:
#   - change the the loops noise bandwidth. From the BL and phiPM (phase margine) parameters the
#     Kp and Ki parameters will be computed in the run function.
#   - change the power of the added white noise in the complex channel
#   - apply phase jump on the PLLs input signal
#   - apply frequency jump on the loops input signal
#   - start a triangle or sine shape frequency change on the loops input signal. The sweep
#     will be one sweep period long. The sweep period and sweep amplitude could be controlled
#     at instance construction time.
# The run() function could be invoked multiple times.
# After the run() function completes we could investigate the result using the show() functions.
#
# The algorithm of the computation of the Kp and Ki parameters:
#   BL    : SSB noise bandwidth of the loop
#   phiPM : desired phase margine in the loop
#   Ts    : sampling time at the input of the loop
#   N     : decimation factor in the loop
#   For the Type2 loop filter:
#       rho = tan((phiPM/90.0)*pi / 2)
#       K = 4 * BL * rho / (1 + rho)
#       w0 = K / rho
#       Kp = K
#       Ki = w0 * Ts * N
#   For the Type3 loop filter:
#       phi = (phiPM + 90) / 2
#       rho = tan((phiPM/90.0)*pi / 2)
#       K = 4 * BL * (2 * rho - 1) / (2 * rho + 3)
#       w0 = K / rho
#       Kp = K
#       Ki = w0 * Ts * N
#
# Handling in the frequency jump
#   At the constructor
#       dwacc = 0.0
#   At the beginning of run()
#       dwacc += df * 2.0 * math.pi
#   dwacc is added to the sweep gwnwrator output to form the frequency control of the carrier NCO
#
# Handling the phase jump
#   In the run() function in the first call of the next() function of the bitnco
#   the dteta parameter is passed. For the following call 0.0 is passed as the
#   dteta parameter.
#
# The block diagramm of the test is the following:
#
# >>------------------------>|
# dwacc                      |
#               w0sw         |          A=1.0                                                    Gn
#           fs   |  Asw      |       fs   |    w0                             Gs=1.0   |<--------<<
#           |    |   |       |       |    |    | carrier[n]         modcarr[n]   |     |
#   swend |------------|     |     |------------|       |------------|       |------------|
#   -----<|            |   |---|   |  Complex   |>----->|   BPSK     |>----->|  Complex   | x[n]
# swstart | SweeperGen |>->| + |>->|    NCO     |       |            |       |            |>----->|
# >>----->|            |   |---|   |  (Carrier) |   |>->| Modulator  |   |>->|  Channel   |       |
#         |------------|           |------------|   |   |------------|   |   |------------|       |
# swmode         |          sweep[n] | |      |     |                    |                        |
# >>------------>|                   | |  invertQ=0 |                    |                        |
#                                    | |            |                    |                        |
# >>-------------------------------->| teta0=0      |pam[n]              |                        |
# dteta                                             |                    |                        |
#                                                   |                    |noise[n]                |
#          fs  w0bit Abit=1.0          N=1          |                    |                        |
#           |   |    |                  |           |                    |                        |
#         |------------|bnext[n]   |------------|   |   |------------|   |                        |
# dwbit=0 |  Complex   |>--------->|    Rect    |   |   |  Complex   |   |                        |
# >>----->|            |           |    Pulse   |>->|   |  Noise     |>->|                        |
#         |    NCO     |>-----     |    PAMGen  |       |  Gen       |                            |
#         |------------|bclk[n]    |------------|       |------------|                            |
# dteta     |   |    |                  |                                                         |
# >>------->|   |  invertQ=0           bmode                                                      |
#            teta0=0                                                                              |
#                                                                                                 |
#    |<------------------------------------------------------------------------------------------<|
#    |
#    |    |---------------------------------------------------|
#    |>-->|                                                   |>--------> z[n]
#   Kp    |                                                   |
#   ----->|                                                   |>--------> Iu[k]
#   Ki    |               PLLorCostasLoop                     |>--------> Qu[k]
#   ----->|                                                   |>--------> e[k]
#         |                                                   |>--------> c[k]
#         |                                                   |
#         |                                                   |
# >>----->|                                                   |>--------> vk[n]
# pdmode  |---------------------------------------------------|
#           |   |   |   |   |                  |       |   |
#          fs  w0   N  Nfir bw               lfsel  pdtype mu
#
#
#
# BL      |-----------------------|
# >>----->|                       |     
# phiPM   |                       |>------ Kp
# >>----->| Kp and Ki computation |
#   lfsel |                       |>------ Ki
#   ----->|                       |
#         |-----------------------|
#
# The constructor:
#   PLLorCostasLoopTest(fs, f0, N, Nfir, bw, f0sw, Asw, f0bit, length, lfsel=1, pdtype=0, bmode=0, mu=1.0)
#   Parameters passed to the constructor:
#       fs     : sampling frequency in Hz
#       f0     : nominal frequency of the carrier and the ComplexNCO in the loop in Hz
#       N      : decimation factor in the loop
#       Nfir   : Length of the FIR decimator
#       bw     : Cutoff frequency of the FIR decimator in Hz
#       f0sw   : Frequency of the sweep waveform
#       Asw    : Amplitude of the sweep waveform
#       f0bit  : bit frequency in Hz
#       length : maximum lenth of the simulation in samples
#       lfsel  : loop filter selector
#                0 : type2
#                1 : type3
#                default is type3
#       pdtype : type of the phase detector
#                0 : Carrier
#                1 : BPSK
#                2 : QPSK
#                default is Carrier
#       bmode  : bit generating mode
#                mode == 0 : random bit sequence 
#                mode == 1 : {-1.0,1.0,-1.0,1.0,...} series
#                mode == 2 : {1.0,1.0,1.0,1.0,...} series
#                mode == 3 : {-1.0,-1.0,-1.0,-1.0,...} series
#                mode == 4 : {-1.0,-1.0,1.0,1.0,...} series
#                default : 0
#       mu     : safety factor for phase detector
#                default value: 1.0
# The run fuction:
#   run(samples, BL, phiPM, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0)
#   Input:     run() function
#       samples: run length in samples
#       BL     : SSB noise bandwidth of the loop in Hz
#       phiPM  : desired phase margin in the loop in degree
#       M      : legth of the variable moving average filrter
#                default : 1
#       dteta  : phase jump in the carrier at the beginning of the run in radian
#                default : 0.0
#       df     : frequency jump in the carrier at the beginning of the run in Hz
#                default : 0.0
#       Gs     : signal gain
#                defaukt : 1.0
#       Gn     : noise gain
#                default : 0
#       startsw: sweep start flag, strats a sweep at the beginning of the run if 
#                sweep not already in progress
#                default : 0
#       swmode : sweep mode
#                0 : triangle wave
#                1 : sine wave
#                default : 0
#       pdmode : phase detector mode selector
#                0 : linear mode
#                1 : pricipal value mode
#                default : 0
#       openloop:open loop mode enable flag
#                0 : closed loop runs the test with closed loop
#                1 : open loop runs the test with open loop
#                if this flag change from the value in the previous run it rests 
#                default : 0
#                some parts of the loop.
#
#   Output:    run() function
#       none   : The run function each time it finished processing 10 * N samples
#                print a character on the console.
#                o : open loop, no sweep
#                O : open loop, sweep
#                c : closed loop, no sweep
#                C : closed loop, sweep
#
# The show functions:
#   show_inout(begin, end)
#       shows the input and output signals of the loop
#   show_loop(begin, end)
#       shows internal signals of the loop
#   show_stat(begin, end)
#       print signal statistics on the console
#
#   Input:    show() function
#       begin : beginning of the signal sample range
#       end   : end of the signal sample range
#
#   Output:   show() function
#       none  : The function graphically displays the internal signals of the 
#               loop or prints some signal statistict to the console.     
# -----------------------------------------------------------------------------
#

class PLLorCostasLoopTest :
    def __init__(self, fs, f0, N, Nfir, bw, f0sw, Asw, f0bit, length, lfsel=1, pdtype=0, bmode=0, mu=1.0) :
        if ((f0 >= fs / 2) or (f0 <= -fs / 2)) :
            sys.exit("PLLorCostasLoopTest: f0 is out of allowed range")
        if (N <= 0) :
            sys.exit("PLLorCostasLoopTest: N should be greater than zero")
        if (Nfir <= 0) :
            sys.exit("PLLorCostasLoopTest: Nfir should be greater than zero")
        if not ((lfsel == 0) or (lfsel == 1)) :
            sys.exit("PLLorCostasLoopTest: invalid lfsel value")
        if (mu < 1.0) :
            sys.exit("PLLorCostasLoopTest: mu shoud be greater than 1.0")
        if (length <= 0) :
            sys.exit("PLLorCostasLoopTest: length should be greater than zero")
        if not ((bmode >= 0) and (bmode <= 5)) :
            sys.exit("PLLorCostasLoopTest: bmode is out of range")
            
        self.twopi = 2.0 * math.pi
        self.fs = fs
        self.w0 = self.twopi * f0
        self.N = N
        self.Nfir = Nfir
        self.bw = bw
        self.w0sw = self.twopi * f0sw
        self.Asw = self.twopi * Asw
        self.w0bit = self.twopi * f0bit
        self.length = int(length)
        self.lfsel = lfsel
        self.pdtype = pdtype
        self.bmode = bmode
        self.mu = mu   
        self.Mmax = 64    
        
        self.A = 1.0
        self.teta0 = 0.0
        self.invertQ = 0
        self.Abit = 1.0
        self.teta0bit = 0.0
        self.invertQbit = 0
        self.Nb = 1
        self.dwbit= 0.0
        self.BL = 0
        self.phiPM = 0
        self.Ts0 = 1.0 / self.fs
        
        self.n = 0
        self.k = 0
        self.nprint = 0
        self.printlimit = 10 * self.N
        self.swrun = 0
        self.dwacc = 0.0
        self.time = 0.0
        
        self.klength = math.floor(self.length / self.N) + 1
        self.sweep = np.zeros(self.length, np.float)
        self.Icarrier = np.zeros(self.length, np.float)
        self.pam = np.zeros(self.length, np.float)
        self.Ix = np.zeros(self.length, np.float)
        self.Qx = np.zeros(self.length, np.float)
        self.Iz = np.zeros(self.length, np.float)
        self.Iu = np.zeros(self.klength, np.float)
        self.Qu = np.zeros(self.klength, np.float)
        self.e = np.zeros(self.klength, np.float)
        self.c = np.zeros(self.klength, np.float)
        self.cref = np.zeros(self.klength, np.float)
        
        self.carriernco = ComplexNCO(self.fs, self.w0, self.A, self.teta0, self.invertQ)
        self.bitnco = ComplexNCO(self.fs, self.w0bit, self.Abit, self.teta0bit, self.invertQbit)
        self.sweeper = SweeperGen(self.fs, self.w0sw, self.Asw)
        self.pamgen = RectPulsePAMGen(self.bmode, self.Nb)
        self.modulator = BPSKModulator()
        self.noisegen = ComplexNoiseGen()
        self.channel = ComplexChannel()
        self.loop = PLLorCostasLoop(self.fs, self.w0, self.N, self.Nfir, self.bw, self.lfsel, self.pdtype, self.mu)
        
    def run(self, samples, BL, phiPM, M=1, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, startsw=0, swmode=0, pdmode=0, openloop=0) :
        if (M < 1) :
            M = 1
        elif (M > self.Mmax) :
            M = self.Mmax
        if (self.n + int(samples)) > self.length :
            print("Too many samples are requested")
        else :
            if ((BL != self.BL) or (phiPM != self.phiPM)) :
                if (self.lfsel == 0) :
                    rho = math.tan((phiPM / 90.0) * math.pi / 2.0)
                    K = 4.0 * BL * rho / (1.0 + rho)
                    omega0 = K / rho
                    self.Kp = K
                    self.Ki = omega0 * self.Ts0 * self.N
                else :
                    phi = (phiPM + 90.0) / 2.0
                    rho = math.tan((phi / 90.0) * math.pi / 2.0)
                    K = 4.0 * BL * (2.0 * rho - 1.0) / (2.0 * rho + 3.0)
                    omega0 = K / rho
                    self.Kp = K
                    self.Ki = omega0 * self.Ts0 * self.N
                self.BL = BL
                self.phiPM = phiPM
                print("BL =", BL)
                print("phiPM =", self.phiPM)
                print("Kp =", self.Kp)
                print("Ki =", self.Ki)
            self.dwacc += 2.0 * math.pi * df
            if (startsw == 1) :
                if (self.swrun == 1) :
                    print("Sweeper is already running")
                    startsw = 0
                else :
                    self.swrun = 1
            start_time = time.time()
            n = 0
            i = 0
            while (i < int(samples)) :
                sweep,swend = self.sweeper.next(startsw, swmode)
                startsw = 0
                sweep += self.dwacc
                Icarrier,Qcarrier,__,__ = self.carriernco.next(sweep, dteta)
                dteta = 0.0
                if (self.pdtype != 0) :
                    __,__,bnext,__ = self.bitnco.next(self.dwbit,self.teta0bit)
                    pam = self.pamgen.next(bnext)
                else :
                    pam = 1.0
                Imodcarrier,Qmodcarrier = self.modulator.next(Icarrier, Qcarrier, pam)
                Inoise,Qnoise = self.noisegen.next()
                Ix,Qx = self.channel.next(Imodcarrier, Qmodcarrier, Inoise, Qnoise, Gs, Gn)
                
                Iz,Qz,Iu,Qu,e,c,vk= self.loop.next(Ix, Qx, self.Kp, self.Ki, M, pdmode, openloop)
                
                self.Icarrier[self.n] = Icarrier
                self.pam[self.n] = pam
                self.Ix[self.n] = Ix
                self.Qx[self.n] = Qx
                self.Iz[self.n] = Iz
                if (vk == 1) :
                    self.Iu[self.k] = Iu
                    self.Qu[self.k] = Qu
                    self.e[self.k] = e
                    self.c[self.k] = c
                    self.cref[self.k] = sweep # resampling, becouse c is sampled at rate k
                    self.k += 1
                i += 1
                self.n += 1
                if (swend == 1) :
                    self.swrun = 0
                if (open == 1) :
                    if (self.swrun == 0) :
                        statuschar = 'o'
                    else :
                        statuschar = 'O'
                else :
                    if (self.swrun == 0) :
                        statuschar = 'c'
                    else :
                        statuschar = 'C'
                self.nprint += 1
                if (self.nprint >= self.printlimit) :
                    sys.stdout.write(statuschar)
                    sys.stdout.flush()
                    self.nprint = 0
            self.time = self.n * self.Ts0
            print("")
            print("n =", self.n)
            print("k =", self.k)
            print("t =", self.time)
            print("exec time =", (time.time() - start_time))

    def show_inout(self, begin, end) :
        if (begin < 0) or (begin >= end) or (end > self.length) :
            print("show_inout() : invalid index range")
        else :
            kbegin = math.floor(begin / self.N) + 1
            kend = math.floor(end / self.N)
            print("n rate sample range :",begin,":",end)
            print("k rate sample range :",kbegin,":",kend)
            
            n = np.linspace(begin, end, end - begin, endpoint=False)
            k = np.linspace(kbegin, kend, kend - kbegin, endpoint=False)
            
            plt.plot(n,self.pam[begin:end],"C1",label="pam[n]")
            plt.legend(loc='lower right')
            plt.title("PAM baseband signal")
            plt.grid(True)
            plt.figure()
            
            plt.plot(k,self.cref[kbegin:kend],"C2",label="cref[k]")
            plt.legend(loc='lower right')
            plt.title("Sweep frequency at k-rate")
            plt.grid(True)
            plt.figure()
            
            plt.plot(n,self.Ix[begin:end] + 2.0,"C1",label="Ix[n]")
            plt.plot(n,self.Qx[begin:end],"C2",label="Qx[n]")
            plt.legend(loc='lower right')
            plt.title("PLL input")
            plt.grid(True)
            plt.figure()

            plt.plot(n,self.Icarrier[begin:end] + 2.0,"C1",label="Icarrier[n]")
            plt.plot(n,self.Iz[begin:end],"C2",label="Iz[n]")
            plt.legend(loc='lower right')
            plt.title("Carrier and recovered carrier")
            plt.grid(True)
            
            plt.show()
    
    def show_loop(self, begin, end) :
        if (begin < 0) or (begin >= end) or (end > self.length) :
            print("show_pll() : invalid index range")
        else :
            kbegin = math.floor(begin / self.N) + 1
            kend = math.floor(end / self.N)
            print("n rate sample range =",begin,":",end)
            print("k rate sample range =",kbegin,":",kend)
            
            k = np.linspace(kbegin, kend, kend - kbegin, endpoint=False)
            
            plt.plot(k,self.Iu[kbegin:kend],"C1",label="Iu[k]")
            plt.plot(k,self.Qu[kbegin:kend],"C2",label="Qu[k]")
            plt.legend(loc='lower right')
            plt.title("Phase detector input")
            plt.grid(True)
            plt.figure()
            
            plt.plot(k,self.e[kbegin:kend],"C1",label="e")
            plt.legend(loc='lower right')
            plt.title("Phase error")
            plt.grid(True)
            plt.figure()
            
            plt.plot(k,(0.5 / math.pi) * self.c[kbegin:kend],"C1",label="c")
            plt.plot(k,(0.5 / math.pi) * self.cref[kbegin:kend],"C2",label="cref")
            plt.legend(loc='lower right')
            plt.title("Frequency control")
            plt.grid(True)
                        
            plt.show()
            

    def show_stat(self, begin, end) :
        if (begin < 0) or (begin >= end) or (end > self.length) :
            print("show_stat() : Invalid index range")
        else :
            kbegin = math.floor(begin / self.N) + 1
            kend = math.floor(end / self.N)
            print("n rate sample range =",begin,":",end)
            print("k rate sample range =",kbegin,":",kend)
            
            N = end - begin
            kN = kend - kbegin
            
            Px = 0.0
            for i in range(begin,end) :
                Ic = self.Ix[i]
                Qc = self.Qx[i]
                Px += Ic * Ic + Qc * Qc
            Px /= N
            
            Pu = 0.0
            for i in range(kbegin,kend) :
                Ic = self.Iu[i]
                Qc = self.Qu[i]
                Pu += Ic * Ic + Qc * Qc
            Pu /= N
            
            MeanIu = 0.0
            VarIu = 0.0
            MeanQu = 0.0
            VarQu = 0.0
            Meane = 0.0
            Vare = 0.0
            Meanc = 0.0
            Varc = 0.0
            for i in range(kbegin,kend) :
                Iu = self.Iu[i]
                Qu = self.Qu[i]
                e = self.e[i]
                c = self.c[i]
                MeanIu += Iu
                MeanQu += Qu
                VarIu += Iu * Iu
                VarQu += Qu * Qu
                Meane += e
                Vare += e * e
                Meanc += c
                Varc += c * c
            MeanIu /= kN
            VarIu /= kN
            VarIu -= MeanIu * MeanIu
            MeanQu /= kN
            VarQu /= kN
            VarQu -= MeanQu * MeanQu
            Meane /= kN
            Vare /= kN
            Vare -= Meane * Meane
            Meanc /= kN
            Varc /= kN
            Varc -= Meanc * Meanc
            
            print("PLL input signal power :",Px)
            print("")
            print("Power after decimator :",Pu)
            print("")
            print("mean(Iu) :",MeanIu)
            print("var(Iu) :",VarIu)
            print("")
            print("mean(Qu) :",MeanQu)
            print("var(Qu) :",VarQu)
            print("")
            print("mean(e) :",Meane)
            print("var(e) :",Vare)
            print("")
            print("mean(c) :",Meanc)
            print("var(c) :",Varc)
            print("")


# -----------------------------------------------------------------------------
#                              BitRecoveryLoopTest
# -----------------------------------------------------------------------------
# This class implements a testbed for the BitRecoveryLoop class.
# It uses the following classes from loops.py and algorithms.py:
#   - BitRecoveryLoop class
#       - SquareWaveNCO
#       - GatedIntegrator
#       - EdgePhaseDetector
#       - Type2LoopFilter
#   - SquareWaveNCO
#   - RectPulsePAMGen
#   - ComplexNoiseGen
#   - ComplexChannel
# During the tests we could
#   At instance construction time:
#       - Set the bit pattern to be used
#       - Set the base bit frequency
#       - Set the sampling rate
#       - Set the threshold for the edge phase detector
#   At runtime
#       - Add phase jump to the PAM generation
#       - Add jump to the bit rate
#       - Change the loop bandwidth
#       - Change the signal and noise gain to change the SNR
# The algorithm of the computation of the Kp and Kiparameters:
#   BL    : SSB noise bandwidth of the loop
#   phiPM : desired phase margine in the loop
#   Ts    : sampling time at the input of the loop
#   N     : decimation factor in the loop
#   For the Type2 loop filter:
#       rho = tan((phiPM/90.0)*pi / 2)
#       K = 4 * BL * rho / (1 + rho)
#       w0 = K / rho
#       Kp = K
#       Ki = w0 * Ts * N
#
# Handling in the frequency jump
#   At the constructor
#       dwacc = 0.0
#   At the beginning of run()
#       dwacc += df * 2.0 * math.pi
#   dwacc is connected to the frequency control of the bit NCO
#
# Handling the phase jump
#   In the run() function in the first call of the next() function of the bitnco
#   the dteta parameter is passed. For the following call 0.0 is passed as the
#   dteta parameter.
#
# The block diagram of the testbed:
#
#          fs  w0bit Abit=1.0          N=1 
#           |   |    |                  | 
#         |------------|bnext[n]   |------------|       |------------|
# dwacc   |  Square    |>--------->|    Rect    |       |  Complex   |
# >>----->|   Wave     |           |    Pulse   |>->|   |  Noise     |>->|
#         |    NCO     |>-----     |    PAMGen  |   |   |  Gen       |   |
#         |------------|bclk[n]    |------------|   |   |------------|   |
# dteta     |   |                       |           |                    |
# >>------->|   |                      bmode        |                    |
#            teta0=0                                | pam[n]             |
#                                                   |                    |
#                                                   |                    |
#                                  |------------|   |                    |
#                 x[n]             |  Complex   |<-<|                    |
#    |<---------------------------<|  Channel   |                        |
#    |                             |            |<-----------------------|
#    |                             |------------|     noise[n]
#    |                                |      |
#    |                               Gs      Gn
#    |
#    |
#    |
#    |
#    |    |---------------------------------------------------|
#    |>-->|                                                   |
#         |                                                   |
#         |                                                   |>--------> Iy[n]
#   Kp    |               BitRecoveryLoop                     |
#   ----->|                                                   |>-------> clk[n]
#   Ki    |                                                   |
#   ----->|                                                   |>--------> v[n]
#         |                                                   |
#         |                                                   |
#         |---------------------------------------------------|
#           |   |       |              | 
#          fs  f0   threshold     openloop
#
#
#
# BL      |-----------------------|
# >>----->|                       |     
# phiPM   |                       |>------ Kp
# >>----->| Kp and Ki computation |
#         |                       |>------ Ki
#         |                       |
#         |-----------------------|
#
# The constructor:
#   BitRecoveryLoopTest(fs, f0, length, bmode=0, threshold=1e-2)
#   Parameters:
#       fs      : sampling frequency in Hz
#       f0      : bit frequency in Hz
#       length  : maximum length of the simulation in samples
#       bmode   : bit generating mode
#                 mode == 0 : random bit sequence 
#                 mode == 1 : {-1.0,1.0,-1.0,1.0,...} series
#                 mode == 2 : {1.0,1.0,1.0,1.0,...} series
#                 mode == 3 : {-1.0,-1.0,-1.0,-1.0,...} series
#                 mode == 4 : {-1.0,-1.0,1.0,1.0,...} series
#                 default : 0
#       threshold:The threshold for the EdgePhaseDetector in the loop
#                 default : 1e-2
#
# The run() function
#   Inputs:
#       samples: run length in samples
#       BL     : SSB noise bandwidth of the loop in Hz
#       phiPM  : desired phase margin in the loop in degree
#       dteta  : phase jump in the carrier at the beginning of the run in radian
#                default : 0.0
#       df     : frequency jump in the carrier at the beginning of the run in Hz
#                default : 0.0
#       Gs     : signal gain
#                default : 1.0
#       Gn     : noise gain
#                default : 0
#       openloop:open loop mode enable flag
#                0 : closed loop runs the test with closed loop
#                1 : open loop runs the test with open loop
#                default : 0
#                if this flag change from the value in the previous run it rests
#                some parts of the loop.
#   Output:
#       none   : The run function each time it finished processing 10 * N samples
#                print a '.' character on the console.
#       
# The show functions:
#   show_inout(begin, end) 
#       shows graphically the internal signals of the test and theloop
#           
#   Input:    show() function
#       begin : beginning of the signal sample range
#       end   : end of the signal sample range
#
#   Output:   show() function
#       none  : The function graphically displays the internal signals of the
#               loops or prints some signal statistict to the console.     
#
# -----------------------------------------------------------------------------

class BitRecoveryLoopTest :
    def __init__(self, fs, f0, length, bmode=0, threshold=1e-2) :
        self.fs = fs
        self.w0 = 2.0 * math.pi * f0
        self.f0 = f0
        self.length = int(length)
        self.bmode = bmode
        if self.bmode == 0 :
            self.gain = 2.0
        else :
            self.gain = 1.0
        self.threshold = threshold
        self.N = math.floor(fs/f0) + 1 # number of samples during bit time
        
        self.bclk = np.zeros(self.length, np.float)
        self.clk = np.zeros(self.length, np.float)
        self.izc = np.zeros(self.length, np.float)
        self.qzc = np.zeros(self.length, np.float)
        self.pam = np.zeros(self.length, np.float)
        self.x = np.zeros(self.length, np.float)
        self.Iy = np.zeros(self.length, np.float)
        self.Qy = np.zeros(self.length, np.float)
        self.e = np.zeros(self.length, np.float)
        self.c = np.zeros(self.length, np.float)
        self.cref = np.zeros(self.length, np.float)
        
        self.bitnco = SquareWaveNCO(self.fs, self.w0, 1.0, 0.0)
        self.pamgen = RectPulsePAMGen(self.bmode, 1)
        self.noisegen = ComplexNoiseGen()
        self.channel = ComplexChannel()
        self.loop = BitRecoveryLoop(self.fs, self.f0, self.threshold)
        
        self.n = 0
        self.k = 0
        self.nprint = 0
        self.printlimit = 10 * self.N
        self.dwacc = 0.0
        self.time = 0.0
        
        self.BL = 0
        self.phiPM = 0
        self.Ts0 = 1.0 / self.fs


    def run(self, samples, BL, phiPM, dteta=0.0, df=0.0, Gs=1.0, Gn=0.0, openloop=0) :
        if (self.n + int(samples)) > self.length :
            print("Too many samples are requested")
        else :
            if ((BL != self.BL) or (phiPM != self.phiPM)) :
                rho = math.tan((phiPM / 90.0) * math.pi / 2.0)
                K = 4.0 * BL * rho / (1.0 + rho)
                omega0 = K / rho
                self.Kp = K
                self.Ki = omega0 * self.Ts0 * self.N
                self.BL = BL
                self.phiPM = phiPM
                print("BL =", BL)
                print("phiPM =", self.phiPM)
                print("Kp =", self.Kp)
                print("Ki =", self.Ki)
            self.dwacc += 2.0 * math.pi * df
            start_time = time.time()
            i = 0
            while (i < int(samples)) :
                bclk,bnext,qnext = self.bitnco.next(self.dwacc, dteta)
                dteta = 0.0
                pam = self.pamgen.next(bnext)
                noise,qn = self.noisegen.next()
                x,qx = self.channel.next(pam, 0.0, noise, 0.0, Gs, Gn)
                clk, izc, qzc, Iy, Qy, e, c, v = self.loop.next(x, self.Kp, self.Ki, openloop, self.gain)
                self.bclk[self.n] = bclk
                self.clk[self.n] = clk
                self.izc[self.n] = izc
                self.qzc[self.n] = qzc
                self.pam[self.n] = pam
                self.x[self.n] = x
                self.Iy[self.n] = Iy
                self.Qy[self.n] = Qy
                self.e[self.n] = e
                self.c[self.n] = c
                self.cref[self.n] = self.dwacc
                self.n +=1
                i += 1
                if v == 1 :
                    self.k += 1
                self.nprint += 1
                if (self.nprint >= self.printlimit) :
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    self.nprint = 0
            self.time = self.n * self.Ts0
            print("")
            print("n =", self.n)
            print("k =", self.k)
            print("t =", self.time)
            print("exec time =", (time.time() - start_time))

    def show(self, begin, end) :
        if (begin < 0) or (begin >= end) or (end > self.length) :
            print("show() : invalid index range")
        else :
            print("n rate sample range :",begin,":",end)
            
            n = np.linspace(begin, end, end - begin, endpoint=False)
            
            plt.plot(n,self.bclk[begin:end] + 2.0,"C1",label="bclk[n]")
            plt.plot(n,self.clk[begin:end],"C2",label="clk[n]")
            plt.legend(loc='upper left')
            plt.title("Clocks")
            plt.grid(True)
            plt.figure()
            
            plt.plot(n,self.x[begin:end] + 2.0,"C1",label="x[n]")
            plt.plot(n,self.izc[begin:end],"C2",label="izc[n]")
            plt.plot(n,self.qzc[begin:end],"C3",label="qzc[n]")
            plt.legend(loc='upper left')
            plt.title("Integrators inputs")
            plt.grid(True)
            plt.figure()

            plt.plot(n,self.pam[begin:end] + 3.0,"C1",label="pam[n]")
            plt.plot(n,self.Iy[begin:end],"C2",label="Iy[n]")
            plt.legend(loc='upper left')
            plt.title("Baseband")
            plt.grid(True)
            plt.figure()
            
            plt.plot(n,self.Iy[begin:end] + 3.0,"C1",label="Iy[n]")
            plt.plot(n,self.Qy[begin:end],"C2",label="Qy[n]")
            plt.legend(loc='upper left')
            plt.title("Integrators output")
            plt.grid(True)
            plt.figure()
            
            plt.plot(n,self.e[begin:end],"C1",label="e[n]")
            plt.legend(loc='upper right')
            plt.title("Error")
            plt.grid(True)
            plt.figure()
            
            plt.plot(n,self.c[begin:end],"C1",label="c[n]")
            plt.plot(n,self.cref[begin:end],"C2",label="cref[n]")
            plt.legend(loc='upper right')
            plt.title("NCO control")
            plt.grid(True)
            
            plt.show()

    def show_error(self, begin, end) :
        if (begin < 0) or (begin >= end) or (end > self.length) :
            print("show() : invalid index range")
        else :
            print("n rate sample range :",begin,":",end)
            
            n = np.linspace(begin, end, end - begin, endpoint=False)
            
            plt.plot(n,self.e[begin:end],"C1",label="e[n]")
            plt.legend(loc='upper right')
            plt.title("Error")
            plt.grid(True)
            plt.figure()
            
            plt.plot(n,self.c[begin:end],"C1",label="c[n]")
            plt.plot(n,self.cref[begin:end],"C2",label="cref[n]")
            plt.legend(loc='upper right')
            plt.title("NCO control")
            plt.grid(True)
            
            plt.show()

            

# -----------------------------------------------------------------------------
#                              FrequencyLockedLoopTest
# -----------------------------------------------------------------------------
# This class implements a testbed for the FrequencyLockedLoop class.
# It uses the following classes from loops.py and algorithms.py:
#   - FrequencyLoop class
#       - ComplexNCO
#       - ComplexMixer
#       - FFTFrequencyDetector
#       - Type2LoopFilter
#   - ComplexNCO
#   - RectPulsePAMGen
#   - ComplexChannel
# During the tests we could
#   At instance construction time:
#       - Set the base carrier frequency
#       - Set the sampling rate
#       - Set the size of the FFT
#       - Set the carrier power threshold
#   At runtime
#       - Add a frequency jump to the carrier frequency
#       - Change the loop bandwidth
#       - change the signal and noise gain to change the SNR
#
# The algorithm of the computation of the Kp and Ki parameters:
# We use Type2 loop filter
#   BL    : SSB noise bandwidth of the loop
#   phiPM : desired phase margine in the loop
#   Ts    : sampling time at the input of the loop
#   Nfft  : FFT size
#   For the Type2 loop filter:
#       rho = tan((phiPM/90.0)*pi / 2)
#       K = 4 * BL * rho / (1 + rho)
#       w0 = K / rho
#       Kp = K
#       Ki = w0 * Ts * Nfft
#
# Handling in the frequency jump
#   At the constructor
#       dwacc = 0.0
#   At the beginning of run()
#       dwacc += df * 2.0 * math.pi
#   dwacc is connected to the frequency control of the carrier NCO
#
# The block diagram of the testbed:
#
#          fs  w0   Abit=1.0  
#           |   |    |                    
#         |------------|                                |------------|
# dwacc   |  Complex   |carrier[n]                      |  Complex   |
# >>----->|            |>-------------------------->|   |  Noise     |>->|
#         |    NCO     |                            |   |  Gen       |   |
#         |------------|                            |   |------------|   |
# dteta=0   |   |    |                              |                    |
# --------->|   |    |                              |                    |
#            teta0=0 |                              | pam[n]             |
#                 inverQ=0                          |                    |
#                                                   |                    |
#                                  |------------|   |                    |
#                 x[n]             |  Complex   |<-<|                    |
#    |<---------------------------<|  Channel   |                        |
#    |                             |            |<-----------------------|
#    |                             |------------|     noise[n]
#    |                                |      |
#    |                               Gs      Gn
#    |
#    |
#    |
#    |
#    |    |---------------------------------------------------|
#    |>-->|                                                   |>-------> ps[n]
#         |                                                   |>-------> pn[n]
#         |                                                   |
#   Kp    |               FrequencyLockedLoop                 |>-------> Iz[n]
#   ----->|                                                   |>-------> Qz[n]
#   Ki    |                                                   |
#   ----->|                                                   |>--------> e[n]
#         |                                                   |>--------> c[n]
#         |                                                   |>---> binidx[n]
#         |---------------------------------------------------|
#           |   |       |              | 
#          fs  f0   threshold     openloop
#
#
#
# BL      |-----------------------|
# >>----->|                       |     
# phiPM   |                       |>------ Kp
# >>----->| Kp and Ki computation |
#         |                       |>------ Ki
#         |                       |
#         |-----------------------|
#
#
# The constructor:
#   FrequencyLockedLoopTest(fs, f0, Nfft, length, threshold=1e-2)
#   Parameters:
#       fs      : sampling frequency in Hz
#       f0      : bit frequency in Hz
#       Nfft    : the length of the FFT
#       length  : maximum length of the simulation in samples
#       threshold:The threshold for the FFTFrequencyDetector in the loop
#                 default : 1e-2
#
# The run() function:
#   run(samples, BL, phiPM, df=0.0, Gs=1.0, Gn=0.0, openloop=0)
#   Inputs:
#       samples: run length in samples
#       BL     : SSB noise bandwidth of the loop in Hz
#       phiPM  : desired phase margin in the loop in degree
#       df     : frequency jump in the carrier at the beginning of the run in Hz
#                default : 0.0
#       Gs     : signal gain
#                default : 1.0
#       Gn     : noise gain
#                default : 0
#       openloop:open loop mode enable flag
#                0 : closed loop runs the test with closed loop
#                1 : open loop runs the test with open loop
#                default : 0
#                if this flag change from the value in the previous run it rests 
#                some parts of the loop.
#   Output:
#       none   : The run function each time it finished processing 10 * N samples
#                print a '.' character on the console.
#       
# The show functions:
#   show_inout(begin, end) 
#       shows graphically the internal signals of the test
#   show_loop(begin, end) 
#       shows graphically the internal signals of the loop
#           
#   Input:    show() functions
#       begin : beginning of the signal sample range
#       end   : end of the signal sample range
#
#   Output:   show() function
#       none  : The function graphically displays the internal signals of the loop
#               or prints some signal statistict to the console.     
#
# -----------------------------------------------------------------------------

class FrequencyLockedLoopTest :
    def __init__(self, fs, f0, Nfft, length, threshold=1e-2) :
        self.fs = fs
        self.f0 = f0
        self.w0 = 2.0 * math.pi * f0
        self.Nfft = Nfft
        self.length = int(length)
        self.llength = math.floor(self.length / self.Nfft) + 1
        self.threshold = threshold
        self.Ts = 1.0 / self.fs
        self.Tfft = self.Ts * self.Nfft

        self.Icarrier = np.zeros(self.length, np.float)
        self.Qcarrier = np.zeros(self.length, np.float)
        self.Iz = np.zeros(self.length, np.float)
        self.Qz = np.zeros(self.length, np.float)
        self.e = np.zeros(self.llength, np.float)
        self.c = np.zeros(self.llength, np.float)
        self.ps = np.zeros(self.llength, np.float)
        self.pn = np.zeros(self.llength, np.float)
        self.binidx = np.zeros(self.llength, np.float)
        self.cref = np.zeros(self.llength, np.float)

        self.BL = 0.0
        self.phiPM = 0.0
        self.Kp = 0.0
        self.Ki = 0.0
        self.A = 1.0
        self.teta0 = 0
        self.invertQ = 0
        self.dwacc = 0.0
        
        self.n = 0
        self.l = 0
        self.nprint = 0
        self.printlimit = 2 * self.Nfft
        self.time = 0.0

        self.carriernco = ComplexNCO(self.fs, self.w0, self.A, self.teta0, self.invertQ)
        self.noisegen = ComplexNoiseGen()
        self.channel = ComplexChannel()
        self.loop = FrequencyLockedLoop(self.fs, self.f0, self.Nfft, self.threshold)

    def run(self, samples, BL, phiPM, df=0.0, Gs=1.0, Gn=0.0, openloop=0) :
        if (self.n + int(samples)) > self.length :
            print("Too many samples are requested")
        else :
            if ((BL != self.BL) or (phiPM != self.phiPM)) :
                rho = math.tan((phiPM / 90.0) * math.pi / 2.0)
                K = 4.0 * BL * rho / (1.0 + rho)
                omega0 = K / rho
                self.Kp = K
                self.Ki = omega0 * self.Tfft
                self.BL = BL
                self.phiPM = phiPM
                print("BL =", BL)
                print("phiPM =", self.phiPM)
                print("Kp =", self.Kp)
                print("Ki =", self.Ki)
            self.dwacc += 2.0 * math.pi * df
            start_time = time.time()
            i = 0
            while (i < int(samples)) :
                Icarrier,Qcarrier,__,__ = self.carriernco.next(self.dwacc, 0.0)
                Inoise,Qnoise = self.noisegen.next()
                Ix,Qx = self.channel.next(Icarrier, Qcarrier, Inoise, Qnoise, Gs, Gn)
                fftvalid,\
                self.Iz[self.n],self.Qz[self.n],\
                e,c,\
                ps,pn,\
                binidx = self.loop.next(Ix, Qx, self.Kp, self.Ki, openloop)
                self.Icarrier[self.n] = Icarrier
                self.Icarrier[self.n] = Icarrier
                if (fftvalid == 1) :
                    self.e[self.l] = e
                    self.c[self.l] = c
                    self.ps[self.l] = ps
                    self.pn[self.l] = pn
                    self.binidx[self.l] = binidx
                    self.cref[self.l] = self.dwacc
                    self.l += 1
                i += 1
                self.n += 1
                self.nprint += 1
                if (self.nprint >= self.printlimit) :
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    self.nprint = 0
            self.time = self.n * self.Ts
            print("")
            print("n =", self.n)
            print("l =", self.l)
            print("t =", self.time)
            print("exec time =", (time.time() - start_time))


    def show_inout(self, begin, end) :
        if (begin < 0) or (begin >= end) or (end > self.length) :
            print("show() : invalid index range")
        else :
            print("n rate sample range :",begin,":",end)
            
            n = np.linspace(begin, end, end - begin, endpoint=False)
            
            plt.plot(n,self.Icarrier[begin:end] + 2.0,"C1",label="Icarrier[n]")
            plt.plot(n,self.Iz[begin:end],"C2",label="Iz[n]")
            plt.legend(loc='upper left')
            plt.title("Input and recoveerd carrier")
            plt.grid(True)
            
            plt.show()


    def show_loop(self, begin, end) :
        if (begin < 0) or (begin >= end) or (end > self.length) :
            print("show() : invalid index range")
        else :
            lbegin = math.floor(begin / self.Nfft) + 1
            lend = math.floor(end / self.Nfft)
            
            print("n rate sample range :",begin,":",end)
            print("l rate sample range :",lbegin,":",lend)
            
            l = np.linspace(lbegin, lend, lend - lbegin, endpoint=False)
            
            plt.plot(l,self.e[lbegin:lend],"C1",label="e[l]")
            plt.legend(loc='lower right')
            plt.title("Frequency error")
            plt.grid(True)
            plt.figure()
            
            plt.plot(l,self.binidx[lbegin:lend],"C1",label="binidx[l]")
            plt.legend(loc='lower right')
            plt.title("Carrier bin index")
            plt.grid(True)
            plt.figure()
            
            plt.plot(l,self.c[lbegin:lend],"C1",label="c[l]")
            plt.plot(l,self.cref[lbegin:lend],"C2",label="cref[l]")
            plt.legend(loc='lower right')
            plt.title("Frequency control")
            plt.grid(True)
            plt.figure()
            
            plt.plot(l,self.ps[lbegin:lend],"C1",label="ps[l]")
            plt.plot(l,self.pn[lbegin:lend],"C2",label="pn[l]")
            plt.legend(loc='lower right')
            plt.title("Bin powers")
            plt.grid(True)
            
            plt.show()



