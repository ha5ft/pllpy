# ******************************************************************************
# *                             Loops
# ******************************************************************************
# * File name:      loops.py
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
# This file contains python classes implementing PLL a related loops. You could
# find the following loops here:
#   PLLorCostasLoop (a loop which behaves as a phase locked loop or a costas loop
#       depending on how its phase detector is working.
#   FLL (frequency locked loop)
#   BitSyncLoop (Loop for bit synchronization for PAM signals)
# All the classes has a common interface:
#   Constructor: It creates an instance of the classes.
#                It receives the parameters for the instance
#                It allocate the necessary resources for the instance
#   next()     : It runs the algorithm for one signal sample
#                It receives the samples of the input signals
#                It outputs the samples of the output signals
#                Optionally it could receive run-time changable parameters
#   debug()    : Calling after next it returns the values of several internal
#                variables.
#   reset()    : It reinitialize the instance state
# ==============================================================================

from algorithms import *

# -----------------------------------------------------------------------------
#                              PLLorCostasLoop
# -----------------------------------------------------------------------------
# This class implements a PLL and a Costas Loop.
# This class uses classes from algorithms.py
# This loop behaves as a phase locked loop or a costas loop depending on how 
# its phase detector is working. The behavior of the phase detector is controlled
# by the pdtype parameter passed to the constructor of the class.
# The loop could synchronize to
#   - unmodulated carrier
#   - BPSK modulated carrier
#   - QPSK modulated carrier
# The loop has a Type2 or Type3 loop filter selectable at contruction time. It has
# a decimator. The dacimation factor could be selected at construction time. The
# bandwidth of the lowpass filter and the loop filter as well as the phase margine
# of the loop filter could be changed at runtime.
# The class implements a lock detector described in the algorithms.py file.
# In the current implementation the loop outputs the output signals of its NCO and
# the lock signal.
# The diagramm of the loop:
#
#                                       N                     M    
#                                       |                     |   
#                                       |                     |           
#                                 |----------|vk[n]      |----------|
#      |----------|               |          |>--->      | Variable |
# x[n] |  Complex |     y[n]      | Complex  |           | Moving   |
# >>-->|          |>------------->| FIR      |>----------| Average  |------->|
#      |   Mixer  |               | Decimator|    v[k]   | Filter   |        |
#      |----------|               |----------|           |----------|        |
#            | z1[n]                                          |              |
# z[n]       |                                                |         u[k] |
# <<--------<|z[n]                                          Mmax             |
#            |                * Depending on lfsel                           |
#         A  |  teta0           Type2 or Type3                  pdmode       |
#         |  |  |               Filter is used                     |         |
# izc[n]|---------|               |----------|                |----------|   |
# <<---<| Complex |               |Type2 or 3|                |          |   |
#       |         |<-------------<|   Loop   |<--------------<|  Phase   |<-<|
# <<---<|   NCO   |<-|     c[k]   |  Filter* |      e[k]      | Detector |
# qzc[n]|---------|  |            |----------|                |----------|
#         |  |  |    |               |    |                        |
#        fs  |  w0  dteta=0.0        Kp   Ki                      pdtype
#            |
#         invertQ=1
#
# On the diagramm the signals x,y,u,v,z,s and mean[l] are complex signals represented
# by their inphase (real part) and quadrature phase (imaginary part) components.
# The sample rate for the various indexes:
#   n : fs
#   k : fs/N
#   l : fs/(N*Nl)
# The PLL could use Type2 or Type3 loop filter. The selected type of loopfilter
# will be instantiated in the constructor.
# Parameters passed to the constructor:
#   fs     : the input sampling frequency
#   w0     : nominal angular frequency of the carrier and the ComplexNCO
#   N      : decimation factor in the PLL loop
#   Nfir   : length of the fIR decimator
#   bw     : cutoff frequency of the FIR decinmator in Hz
#   lfsel  : loop filter selector
#            0 : type2
#            1 : type3
#            default is type3
#   pdtype : type of the phase detector
#            0 : Carrier
#            1 : BPSK
#            2 : QPSK
#            default is Carrier
#   mu     : safety factor for phase detector
#            default value: 1.0
#
# Inputs: (next())
#   Ix : current sample of the inphase input
#   Qx : current sample of the quadrature phase input
#   Kp : current value of the Kp parameter of the loop filter
#   Ki : current value of the Ki parameter of the loop filter
#   M  : actual length of the filter
#        default : 1
#   pdmode : phase detector mode selector
#            0 : linear mode
#            1 : pricipal value mode
#            default is the linear mode
#   openloop : open loop enable flag
#              0 : closed loop
#              1 : open loop
# 
# Outputs: (next())
#   Iz      : current sample of the NCO inphase output
#   Qz      : current sample of the NCO quadrature phase output
#   Iu     : current sample of the inphase signal at the phase detector input
#   Qu     : current sample of the quadrature phase signal at the phase detector input
#   e      : current sample of the phase error
#   c      : current sample of the NCO control signal
#   vk     : valid flag for the Tu,Qu,e,c samples
# -----------------------------------------------------------------------------

class PLLorCostasLoop :
    def __init__(self, fs, w0, N, Nfir, bw, lfsel=1, pdtype=0, mu=1.0) :
        if w0 >= math.pi * fs :
            sys.exit("PLL: w0 is too high")
        if N <= 0:
            sys.exit("PLL: zero or negative N value")
        if Nfir <= 0:
            sys.exit("PLL: zero or negative Nfir value")
#        if Ml <= 0:
#            sys.exit("PLL: zero or negative Ml value")
        if not ((lfsel == 0) or (lfsel == 1)) :
            sys.exit("PLL: wrong lfsel value")
        if not ((pdtype == 0) or (pdtype == 1) or (pdtype == 2)) :
            sys.exit("PLL: wrong pmode value")
        # store parameters for later use in instance variables
        self.fs = fs
        self.w0 = w0
        self.N = N
        self.Nfir = Nfir
        self.bw = bw
        self.lfsel = lfsel
        self.pdtype = pdtype
        self.mu = mu
        # creating and initializing internal signals
        self.Iy = 0.0
        self.Qy = 0.0
        self.Iu = 0.0
        self.Qu = 0.0
        self.vk = 0     # valid flag for u, triggers sampling rate k functions
        self.Iv = 0.0
        self.Qv = 0.0
        self.e = 0.0
        self.c = 0.0
        self.Iz = 0.0
        self.Qz = 0.0
        self.izc = 0
        self.qzc = 0
        # create instance variables
        self.n = 0    # sample count at n-rate
        self.k = 0    # sample count at k-rate
        self.Mmax = 256 # maximum length of the complex variable moving average filter
        self.Ts = 1.0 / self.fs
        self.A = 1.0
        self.teta0 = 0.0
        self.invertQ = 1
        self.openloop = 0
        # create algorithm's instances
        self.mixer = ComplexMixer()
        self.b = ss.firwin(self.Nfir, bw, window='blackmanharris', fs=self.fs)
        self.decimator = ComplexFIRDecimator(self.b, self.N)
        self.varfilter = ComplexVariableMovingAverageFilter(self.Mmax)
        self.detector = PhaseDetector(self.pdtype, self.mu)
        if self.lfsel == 0 :
            self.loopfilter = Type2LoopFilter(self.Ts, self.Mmax)
            print("PLLorCostasLoop: Type2 loop filter")
        else :
            self.loopfilter = Type3LoopFilter(self.Ts, self.Mmax)
            print("PLLorCostasLoop: Type3 loop filter")
        self.nco = ComplexNCO(self.fs, self.w0, self.A, self.teta0, self.invertQ)
        
    def next(self, Ix, Qx, Kp, Ki, M=1, pdmode=0, openloop=0) :
        self.vk = 0
        self.vl = 0
        if (M < 1) :
            M = 1
        elif (M > self.Mmax) :
            M = self.Mmax
        if (openloop == 1) :
            c = 0
            if (self.openloop != 1) :
                self.detector.reset()
                self.loopfilter.reset()
        else :
            c = self.c
            if (self.openloop != 0) :
                self.detector.reset()
                self.loopfilter.reset()
        self.openloop = openloop
        self.Iz,self.Qz,self.izc,self.qzc = self.nco.next(c, 0.0)
        self.Iy,self.Qy = self.mixer.next(Ix, Qx, self.Iz, self.Qz)
        Iv,Qv,self.vk = self.decimator.next(self.Iy, self.Qy)
        if (self.vk == 1) :
            self.Iv = Iv
            self.Qv = Qv
            self.Iu,self.Qu = self.varfilter.next(self.Iv, self.Qv, M)
            self.e = self.detector.next(self.Iu, self.Qu, pdmode)
            self.c = self.loopfilter.next(self.e, Kp, Ki)
            self.k += 1
        self.n += 1
        return self.Iz, self.Qz, self.Iu, self.Qu, self.e, self.c, self.vk
        

    def reset(self) :
        # reset signals
        self.Iy = 0.0
        self.Qy = 0.0
        self.Iu = 0.0
        self.Qu = 0.0
        self.vk = 0     # valid flag for u, triggers samplin rate k functions
        self.e = 0.0
        self.c = 0.0
        self.Iz = 0.0
        self.Qz = 0.0
        self.izc = 0
        self.qzc = 0
        # reset instance variables
        self.n = 0    # sample count at n-rate
        self.k = 0    # sample count at k-rate
        # reset algorithms instances
        self.mixer.reset()
        self.decimator.reset()
        self.detector.reset()
        self.loopfilter.reset()
        self.nco.reset()




# -----------------------------------------------------------------------------
#                              FLL
# -----------------------------------------------------------------------------
# This class implements an FLL
# It is simmilar to the PLL with the following differences:
#   - No complex variable moving averiging filter
#   - Uses FFT frequency detector
#   - Only Type2 loop filter available
#   - There is an additional integrator between the loop filter and the NCO
#     we need this becouse the output of the error detector is frequency not
#     phase as in a PLL. This integrator is implemented as inline code.
# This class uses classes from algorithms.py
# The diagramm of the loop:
#
#                                                               Nfft  fs
#                                                                |    |
#                                                             |----------|
#                                                    threshold|          |>-->> binidx[l]
#      |----------|                                  -------->|   FFT    |>-->> pn[n]
# x[n] |  Complex |     y[n]                                  | Frequency|>-->> ps[n]
# >>-->|          |>----------------------------------------->| Detector |----> fftvalid
#      |   Mixer  |                                           |          |-->|
#      |----------|                                           |----------|   |
#           |                                                      |         |
# z[n]      |                                                   threshold    |->> e[n]
# <<-------<|      invertQ=1                                                  |
#           |         |                     |----|                           |
#           |      fs | w0             |>-->|z^-1|>-->| c[l-1]           e[l]|
#           |       | |  |             |    |----|    |                      |
#           |    |----------|          |              |       |----------|   |
#           |<--<| Complex  |          |            |---|     |  Type2   |   |
#                |          |<--------<|<----------<| + |<---<|   Loop   |<-<|
#                |   NCO    |          |c[l]        |---| r[l]|  Filter  |
#                |----------|          |                      |----------|
#                   | |  |             |                         |    |
#                 A=1 |teta0=0         |                         Kp   Ki
#                     |                |
#                   dteta=0            |                        
#                                      |
# <<----------------------------------<|
# c[n]
#
# On the diagramm the signals x,y,u are complex signals represented
# by their inphase and quadrature phase components.
# The sample rate for the various indexes:
#   n : fs
#   l : fs/(Nfft)
# Parameters passed to the constructor:
#   fs          : the input sampling frequency
#   f0          : base frequency of the ComplexNCO
#   Nfft           : Width of the FFT
#   threshold   : For detecting signal presence
# Inputs: (next())
#   Ix : current sample of the inphase input
#   Qx : current sample of the quadrature phase input
#   Kp : current value of the Kp parameter of the loop filter
#   Ki : current value of the Ki parameter of the loop filter
#   openloop :
# Outputs: (next())
#   Iz     : current sample of the NCO inphase output
#   Qz     : current sample of the NCO quadrature phase output
#   e      : current sample of the frequency error
#   c      : current sample of the NCO control signal
#   ps     : current sample of the carrier power
#   pn     : current sample of the average noise power in bins
#   binidx : current sample of the index of the bin containing the carrier
# -----------------------------------------------------------------------------
#

class FrequencyLockedLoop :
    def __init__(self,fs, f0, Nfft, threshold) :
        self.fs = fs
        self.Ts = 1.0 / fs
        self.f0 = f0
        self.w0 = 2.0 * math.pi * f0
        self.Nfft = Nfft
        self.Tfft = self.Ts * self.Nfft
        self.threshold = threshold
        self.Nloop = 256
        self.A = 1.0
        self.teta0 = 0
        self.invertQ = 1
        self.openloop = 0
        
        self.Iz = 0.0
        self.Qz = 0.0
        self.Iy = 0.0
        self.Qy = 0.0
        self.fftvalid = 0
        self.e = 0
        self.binidx = 0
        self.ps = 0.0
        self.pn = 0.0
        self.cvalid = 0
        self.r = 0.0
        self.c = 0.0
        
        self.l = 0
        
        self.mixer = ComplexMixer()
        self.detector = FFTFrequencyDetector(self.fs, self.Nfft, self.threshold)
        self.loopfilter = Type2LoopFilter(self.Ts, self.Nloop)
        self.nco = ComplexNCO(self.fs, self.w0, self.A, self.teta0, self.invertQ)
        
    def next(self, Ix, Qx, Kp, Ki, openloop=0) :
        self.fftvalid = 0
        if (openloop == 1) :
            c = 0
            if (self.openloop != 1) :
                self.detector.reset()
                self.loopfilter.reset()
        else :
            c = self.c
            if (self.openloop != 0) :
                self.detector.reset()
                self.loopfilter.reset()
        self.openloop = openloop
        self.Iz,self.Qz,self.izc,self.qzc = self.nco.next(c, 0.0)
        self.Iy,self.Qy = self.mixer.next(Ix, Qx, self.Iz, self.Qz)
        self.fftvalid,\
        e,binidx,\
        pn,ps,\
        cvalid = self.detector.next(self.Iy, self.Qy)
        if (self.fftvalid == 1) :
            self.e = 2.0 * math.pi * e # frequency detector output is in Hz
            self.binidx = binidx
            self.ps = ps
            self.pn = pn
            self.cvalid = cvalid
            self.r = self.loopfilter.next(self.e, Kp, Ki)
            self.c += self.Tfft * self.r # backward Euler integration
            self.l += 1
        return self.fftvalid, self.Iz, self.Qz, self.e, self.c, self.ps, self.pn, self.binidx

                

# -----------------------------------------------------------------------------
#                              BitRecoveryLoop
# -----------------------------------------------------------------------------
# This class implements the bit time recovery for PAM signals.
# Mathematical modell:
# For simplicity I show the model for input signal generated by a pulse shaping
# filter having rectangular impulse response, but the model could be easily 
# extended for PAM signal generated by rised cosine or root rised cosine pulse
# shaping filter or by a pulse shaping filter having one period rised cosine 
# impulse response. This later is used to generate PAM signal for PSK31.
#
#   Ts : sampling time, Ts = 1/fs, fs is the sampling frequency
#   Tb : bit time, Tb = 1/fb
#
#   x[n] = b[k]*A if k*Tb-Tb/2+dTb <= n*Ts < k*Tb+Tb/2+dTb
#
#   A > 0, the amplitude
#   Tb is the bit duration
#   dTb starting delay (starting phase)
#   b[k] is the bit value, i.e. -1 or 1
#        / -1 if bit[n]==1
#   b[k]=|
#        \ 1 if bit[n]==0
#   bit[k] is the logic bit value, i.e. 0 or 1
#
# We use a square wave NCO. Its output is clk[n] and the zero crossibg signals: izc, qzc.
# The kth period of the output is T[k]=1/fnco[k].
# Its is implemented by an instance of the SquareWaveNCO class.
#
# The izc and qzc signals control two gated integrator, which integrate the x[n] signal.
#   t[0] = 0
#          k-1
#   t[k] = Sum(T[k])
#          i=0
#   Iy[k]=Sum{n:t[k]<=nTs<t[k+1]}(x[n])
#
#   Qy[k]=Sum{n:t[k]-T[k-1]/2<=nTs<t[k]+T[k]/2}(x[n])
#
#   Ts=1/fs
#   The gated integrators are implemented by instances of the GatedIntegrator class.
#
# The outputs of the gated integrators are used to compute the phase error
#          / -sign(Iy[k-1])*Qy[k] if sign(Iy[k-1]) != sign(Iy[k])
#   e[k] = |
#          \ 0 otherwise
#   The error detector is implemented by an instance of the EdgePhaseDetector class
# The error is multiplied by the gain parameter to compensate for the fact, that 
#   in case of random bit sequence we do finde an edge in all of the integration periods.
# The error sinal is filtered by the loop filter, which is a Type2 filter.
# It is controlled by
#   Kp : the proportional gain
#   Ki : the integrator gain
# The loop filter is implemented by an instance of the Type2LoopFilter class
#
#  The output c[k] of the loop filter controls the NCO frequency
#   fnco[k] = f0+c[k]
#
# We use the following algorithm blocks:
#   SquareWaveNCO
#   GatedIntegrator
#   EdgePhaseDetector
#   Type2LoopFilter
#
# The block diagram of the loop:
#
#                                                                      Iy[n]
#                                                   |>--------------------->>
#                                                   |
#                  |-----------------------|        |  |---------|   >----->>
#                  |                       |  Iy[k] |  |         |     clk[n]
#         |>------>|    GatedIntegrator    |>-------|--|         |
#         |        |                       |           |         |   >----->>
#         |        |-----------------------|           |         |     izc[n]
#         |               |                            |         |
#         |               | izc[n]                     |         |
#         |               |                            |         |
#         |          fs   | A=1.0       Kp   Ki        |         |
#         |          |    |   |          |   |         |         |
#         |        |------------|      |-------|       |  Edge   |
#   x[n]  |  clk[n]|            | c[k] | Type2 |  e[k] |         | threshold
# >>----->|  <----<| SquareWave |<----<| Loop  |<-----<|  Phase  |<---------
#         |        |     NCO    |      | Filter|       |         |
#         |        |------------|      |-------|       | Detector|
#         |          |    |   |                        |         |
#         |         w0    | teta0=0                    |         |
#         |               |                            |         |
#         |               | qzc[n]                     |         |
#         |               |                            |         |
#         |        |-----------------------|           |         |
#         |        |                       |   Qy[k]   |         |
#         |>------>|    GatedIntegrator    |>--------->|         |
#                  |                       |           |         |
#                  |-----------------------|           |---------|
#
# -----------------------------------------------------------------------------

class BitRecoveryLoop :
    def __init__(self, fs, f0, threshold=1e-2) :
        self.fs = fs
        self.w0 = 2.0 * math.pi * f0
        self.A = 1.0
        self.teta0 = 0.0
        self.threshold = threshold
        self.Iy = 0.0
        self.Qy = 0.0
        self.c = 0.0
        self.e = 0.0
        self.Mlf = 256
        
        
        self.nco = SquareWaveNCO(self.fs, self.w0, self.A, self.teta0)
        self.iIntegrator = GatedIntegrator()
        self.qIntegrator = GatedIntegrator()
        self.detector = EdgePhaseDetector(self.threshold)
        self.loopfilter = Type2LoopFilter(1.0 / fs, self.Mlf)
        
    def next(self, x, Kp, Ki, openloop=0, gain=1.0) :
        if openloop == 0 :
            c = self.c
        else :
            c = 0.0
        clk, izc, qzc = self.nco.next(c, 0.0)
        Iy, iN, iv = self.iIntegrator.next(izc, x)
        Qy, qN, qv = self.qIntegrator.next(qzc, x)
        if (qv == 1) :
            self.Qy = Qy
        if (iv == 1) :
            self.Iy = Iy
            self.e = self.detector.next(self.Iy, self.Qy)
            self.e *= gain
            self.c = self.loopfilter.next(self.e, Kp, Ki)
        return clk, izc, qzc, Iy, Qy, self.e, self.c, iv
    
    
    def reset(self) :
        self.Iy = 0.0
        self.Qy = 0.0
        self.c = 0.0
        self.e = 0.0
        self.nco.reset()
        self.iIntegrator.reset()
        self.qIntegrator.reset()
        self.detector.reset()
        self.loopfilter.reset()



