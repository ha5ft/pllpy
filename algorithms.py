# ******************************************************************************
# *                             Algorithms
# ******************************************************************************
# * File name:      algorithms.py
# * Platform:       ubuntu 20.04 64 bit
# *                 python3 v3.8.10
# *                 numpy v1.17.4
# *                 scipy v1.3.3
# *                 matplotlib v3.1.2
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
# *	Selmeczi János		09.05.2022	Original version
# *
# ******************************************************************************

# ==============================================================================
#                               Description
# ==============================================================================
# This file contains python classes implementing the elementary algorithms used
# in the pllpy peojects. You could find the following algorithms here:
#   ComplexNCO
#   SquareWaveNCO
#   TriangleAndSinePeriodGen
#   GatedRectPAMPulseGen
#   PSK31BaseBandGen
#   BPSKModulator
#   ComplexPhaseModulator
#   ComplexNoiseGen
#   ComplexChannel
#   ComplexMixer
#   ComplexMovingAverageDecimator
#   FIRDecimator
#   ComplexVariableMovingAverageFilter
#   MovingAverageFilter
#   ComplexMovinAverageFilter
#   GatedIntegrator
#   LockDetector
#   PhaseDetector
#   CostasPhaseDetector
#   EdgePhaseDetector
#   FFTFrequencyDetector
#   Type2LoopFilter
#   TypeI3LoopFilter
# All the classes has a common interface:
#   Constructor: It creates an instance of the classes.
#                It receives the parameters for the instance
#                It allocate the necessary resources for the instance
#   next()     : It runs the algorithm for one signal sample
#                It receives the samples of the input signals
#                It outputs the samples of the output signals
#                Optionally it could receive run-time changable parameters
#   reset()    : It reinitialize the instance state
# ==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as ss
import sys
import math
import cmath
from numpy.random import default_rng
import random
import time

# ------------------------------------------------------------------------------
#                        ComplexNCO
# ------------------------------------------------------------------------------
# This class generate complex exponential tone represented by the inphase and
# quadrature phase real signal pair. In addition it generates an impulse when
# the phase crosses zero and another impulse when the phase crosses pi.
#
# The algorithm:
#   Iy[n] = A*cos(phi[n])
#   Qy[n] = A*sin(phi[n])
#   phi[0] = teta0
#   alpha  = phi[n-1] + (w0+dw[n])/fs + dteta[n]
#
#   izc[0] = 0
#            |  1, if Iy[n-1] < 0 and Iy[n] >= 0
#   izc[n] = |  
#            |  0, otherwise
#
#   qzc[0] = 0
#            |  1, if Iy[n-1] < 0 and Iy[n] >= 0
#   qzc[n] = |  
#            |  0, otherwise
#
#   if phi crosses 2*pi or -2*pi then it is decreased or increased by 2*pi
#            |  alpha - 2*pi if alpha >= 2*pi
#   phi[n] = |  alpha        otherwise
#            |  alpha + 2*pi if alpha <= -2*pi
#
# Constructor:
#   nco=ComplexNCO(fs, w0, A, teta0=0.0, invertQ=0)
#   Parameters:
#       A           : absolute value of the output
#       fs          : sampling frequency
#       w0          : angular base frequency of the NCO
#       teta0       : starting phase of the NCO
#       invertQ     : If not zero, the sign of Qy will be changed to its oposite
#                     This will causes the frequency of the z signal to change
#                     to the -(w0 + dw) value. This is a permanent change.
# next method
#   nco.next(dw=0.0, dteta=0.0) :                     
#   Parameters:
#       dw       : angular frequency offset
#       dteta    : phase offset
#                  This input usually should be 0.0.
#                  A single nonzero sample will causes a phase jump in the
#                  output signal.
#                  A step function in this input will produce an angular
#                  frequency jump in the output. The jump value is dteta*fs.
#   Return:
#       Iy       : output inphase signal sample
#       Qy       : output quadrature phase signal sample
#       izc      : output zero crossing flag sample
#       qzc      : output pi crossing flag sample
# -----------------------------------------------------------------------------

class ComplexNCO :
    def __init__(self, fs, w0, A, teta0=0.0, invertQ=0) :
        self.A = A
        self.fs = fs
        self.w0 = w0
        self.teta0=teta0
        self.invertQ = invertQ
        self.izc = 0
        self.qzc = 0
        if invertQ != 0 :
            self.signQ = -1.0
        else :
            self.signQ = 1.0
        self.twopi = 2.0 * math.pi
        # move teta0 to (-2*pi,2*pi) intervall
        if self.teta0 >= self.twopi :
            while (self.teta0 >= self.twopi) :
                self.teta0 -= self.twopi
        if self.teta0 <= -self.twopi :
            while (self.teta0 <= -self.twopi) :
                self.teta0 += self.twopi
        self.phi = self.teta0
        # compute Iy[0] and Qy[0]
        self.Iy = self.A * math.cos(self.phi)
        self.Qy = self.A * math.sin(self.phi)

# This function outputs Iy[n], Qy[n], izc[n], qzc[n] computed in the previous call
# and computes Iy[n+1], Qy[n+1], izc[n+1], qzc[n+1] values
# For n=0 the values has been computed in the constructor or in the reset() function

    def next(self, dw=0.0, dteta=0.0) :
        # save the previously computed samples for the return
        izc = self.izc
        qzc = self.qzc
        Iy = self.Iy
        Qy = self.Qy
        phi = self.phi
        # clear the zero crossings flags becouse we look for if the flags should
        # be set to 1 only
        self.izc = 0
        self.qzc = 0
        # we compute the next phase value
        self.phi += (self.w0 + dw) / self.fs + dteta
        # move the phase to the (-2*pi,2*pi) intervall and
        # compute the next sample for the zero crossing signals
        if self.phi >= self.twopi :
            self.phi -= self.twopi
        elif self.phi <= -self.twopi :
            self.phi += self.twopi
        # we compute the next signal samples
        self.Iy = self.A * math.cos(self.phi)
        self.Qy = self.signQ * self.A * math.sin(self.phi)
        # we compute the next zero crossing samples
        if (Iy < 0) and (self.Iy >= 0) :
            self.izc = 1
        if (Iy > 0) and (self.Iy <= 0) :
            self.qzc = 1
        
        return Iy, Qy, izc, qzc
    
    def reset(self) :
        # reset the phase
        self.phi = self.teta0
        # compute Iy[0] and Qy[0]
        self.Iy = self.A * math.cos(self.phi)
        self.Qy = self.A * math.sin(self.phi)
        # compute the zero crossing signal starting samples
        if self.teta0 == 0.0 :
            self.izc = 1
        else :
            self.izc = 0
        if (self.teta0 == math.pi) or (self.teta0 == -math.pi) :
            self.qzc = 1
        else :
            self.qzc = 0


# ------------------------------------------------------------------------------
#                        SquareWaveNCO
# ------------------------------------------------------------------------------
# This class generate square wave tone.
# In addition it generates an impulse at the rising edge
# and another impulse at the falling edge
# The algorithm:
#   phi[0] = teta0
#   alpha  = phi[n-1] + (w0+dw[n])/fs + dteta[n]
#
#                   / 1 if teta0 >= 0.0
#   y[-1] = y[0] = |
#                   \ 0 if teta0 < 0.0
#
#            / 1 if phi[n-1] < 2*pi  and alpha >= 2*pi
#           /  1 if phi[n-1] > -2*pi and alpha <= -2*pi
#   y[n] = |   y[n-1]
#           \  0 if phi[n-1] < pi  and alpha >= pi
#            \ 0 if phi[n-1] > -pi and alpha < -pi
#
#
#   izc[0] = 0
#
#
#             / 1, if y[n-1] < y and y[n] >= 0
#   izc[n] = |  
#             \ 0, otherwise
#
#
#   qzc[0] = 0
#
#
#             / 1, if y[n-1] < *0 and y[n] >= 0
#   qzc[n] = |  
#             \ 0, otherwise
#
#   if phi crosses 2*pi or -2*pi then it is decresed or increased by 2*pi
#
#             / alpha - 2*pi if alpha >= 2*pi
#   phi[n] = |  alpha        otherwise
#             \ alpha + 2*pi if alpha <= -2*pi
#
# Constructor:
#   nco=SquareWaveNCO(fs, w0, A, teta0=0.0, invertQ=0)
#   Parameters:
#       A           : absolute value of the output
#       fs          : sampling frequency
#       w0          : angular base frequency of the NCO
#       teta0       : starting phase of the NCO
#       invertQ     : If not zero, the sign of Qy will be changed to its oposite
#                     This will causes the frequency of the z signal to change
#                     to the -(w0 + dw) value. This is a permanent change.
# next method
#   nco.next(dw=0.0, dteta=0.0) :                     
#   Parameters:
#       dw       : angular frequency offset
#       dteta    : phase offset
#                  This input usually should be 0.0.
#                  A single nonzero sample will causes a phase jump in the
#                  output signal.
#                  A step function in this input will produce an angular
#                  frequency jump in the output. The jump value is dteta*fs.
#   Return:
#       Iy       : output inphase signal sample
#       Qy       : output quadrature phase signal sample
#       izc      : output zero crossing flag sample
#       qzc      : output pi crossing flag sample
# -----------------------------------------------------------------------------

class SquareWaveNCO :
    def __init__(self, fs, w0, A, teta0=0.0) :
        self.A = A
        self.fs = fs
        self.w0 = w0
        self.teta0 = teta0
        self.twopi = 2.0 * math.pi
        self.izc = 0
        self.qzc = 0
        self.e = 0
        # move teta0 to (-2*pi,2*pi) intervall
        if self.teta0 >= self.twopi :
            while (self.teta0 >= self.twopi) :
                self.teta0 -= self.twopi
        if self.teta0 <= -self.twopi :
            while (self.teta0 <= -self.twopi) :
                self.teta0 += self.twopi
        self.phi = self.teta0
        # we compute y[0]
        if (self.phi >= 0.0) and (self.phi) < math.pi :
            self.y = self.A
        elif (self.phi <= -math.pi) and (self.phi > -self.twopi) :
            self.y = self.A
        else :
            self.y = 0.0

    def next(self, dw=0.0, dteta=0.0) :
        # save the previously computed samples for the return
        izc = self.izc
        qzc = self.qzc
        y = self.y
        phi = self.phi
        # zero the zero crossings flags becouse we look for only if the flags should
        # be set to 1
        self.izc = 0
        self.qzc = 0
        # we compute the next phase value
        self.phi += (self.w0 + dw) / self.fs + dteta
        if self.phi >= self.twopi :
            self.phi -= self.twopi
        elif self.phi <= -self.twopi :
            self.phi += self.twopi
        # we compute the next signal sample
        if (self.phi >= 0.0) and (self.phi < math.pi) :
            self.y = self.A
        elif (self.phi <= -math.pi) and (self.phi > -self.twopi) :
            self.y = self.A
        else :
            self.y = 0.0
        # we compute the next zero crossing samples
        if (y == 0) and (self.y != 0) :
            self.izc = 1
        if (y != 0) and (self.y == 0) :
            self.qzc = 1
        return y, izc, qzc
    
    def reset(self) :
        # reset the phase
        self.phi = self.teta0
        # compute the zero crossing signal starting samples
        if self.teta0 == 0.0 :
            self.izc = 1
        else :
            self.izc = 0
        if (self.teta0 == math.pi) or (self.teta0 == -math.pi) :
            self.qzc = 1
        else :
            self.qzc = 0
        # we compute y[0]
        if (self.phi >= 0.0) and (self.phi) < math.pi :
            self.y = self.A
        elif (self.phi <= -math.pi) and (self.phi > -self.twopi) :
            self.y = self.A
        else :
            self.y = 0.0


# ------------------------------------------------------------------------------
#                           SweeperGen
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

class SweeperGen:
    def __init__(self, fs, w0, A):
        self.fs = fs
        if w0 > 0 :
            self.w0 = w0
        else :
            sys.exit("SweepGen: Negative or zero w0 is not allowed")
        self.A = A
        self.phi = 0.0
        self.started = 0
        self.end = 0
        self.mode = 0
        self.halfpi = math.pi / 2.0
        self.twothirdpi = 3.0 * math.pi / 2.0
        self.twopi = 2.0 * math.pi

    def next(self, start=0, mode=0):
        if self.started == 0 :
            if start != 0 :
                self.mode = mode
                self.started = 1
        end = self.end
        y = 0.0
        if self.end != 0 :
            self.end = 0
            self.started = 0
        if self.started != 0 :
            if self.mode == 0 :
                if self.phi < self.halfpi :
                    y = self.A * (self.phi / self.halfpi)
                elif (self.halfpi <= self.phi) and (self.phi < self.twothirdpi) :
                    y = self.A * (math.pi - self.phi) / self.halfpi
                elif (self.twothirdpi <= self.phi) and (self.phi < self.twopi) :
                    y = self.A * (self.phi - self.twopi) / self.halfpi
                else :
                    y = 0.0
            else :
                y = self.A * math.sin(self.phi)
            alpha = self.phi + self.w0 / self.fs
            if alpha >= self.twopi :
                self.phi = 0.0
                self.end = 1
            else :
                self.phi = alpha
        else :
            y = 0.0
        return y, end

    def reset(self) :
        self.phi = 0.0
        self.started = 0
        self.end = 0
        self.mode = 0
    
# ------------------------------------------------------------------------------
#                              RectPulsePAMGen
# ------------------------------------------------------------------------------
# This class generate samples of an PAM signal with rectangular pulse shape.
#   - The two logic level are 1.0 and -1.0.
#   - 5 different bit series could be generated. One of them is a random series.
#   - The bit length is controlled by a bit clock trigger signal and the
#     N parameter. A new bit value is generated for every Nth clock pulse.
# Algorithm:
#   Bit length control:
#       The clock pulse is a single sample long.
#       For every Nth clock pulse :
#           a new bit value is generated and a sample corresponding
#           to this new bit value is emitted
#       otherwise:
#           a sample corresponding to the last generated bit value
#           is emitted
#   Bit value generation:
#       5 different bit sequence could be generated.
#       The sequence being used is selected by the mode parameter, which is passed
#       to the constructor.
#       The constructor generates the starting bit value.
#       A new bit value is generated for every Nth clock pulse.
#       The bit sequences:
#           mode == 0 : random bit sequence generated by using random.random()
#                       numpy
#           mode == 1 : {-1.0,1.0,-1.0,1.0,...} series
#           mode == 2 : {1.0,1.0,1.0,1.0,...} series
#           mode == 3 : {-1.0,-1.0,-1.0,-1.0,...} series
#           mode == 4 : {-1.0,-1.0,1.0,1.0,...} series
# Constructor:
#   pamgen=RectPulsePAMGen(self, mode, N)
#   Parameters passed to the constructor:
#       mode : selects the bit sequence
#       N    : determines the bit length in clock pulses
# next method
#   pamgen.next(clock)
#   Parameters:
#       clock : a single sample long clock pulse
#               for every Nth clock pulse a new bit value is generated
#   Returnt:
#       b : current sample of the NRZ signal
# ------------------------------------------------------------------------------

class RectPulsePAMGen:
    def __init__(self, mode, N):
        if not ((0 <= mode) and (mode <= 4)) :
            sys.exit("NRZGen: Wrong mode value")
        self.mode = mode
        if N < 0 :
            sys.exit("NRZGen: Wrong N value")
        self.N = N
        self.pattern_length = 4
        self.pattern = [[0.0,0.0,0.0,0.0],\
                        [-1.0,1.0,-1.0,1.0],\
                        [1.0,1.0,1.0,1.0],\
                        [-1.0,-1.0,-1.0,-1.0],\
                        [-1.0,-1.0,1.0,1.0]]
        self.pidx = 0
        self.nclk = 0
        if self.mode == 0 :
            self.b = -2.0 * round(random.random(), 0) + 1.0
        else :
            self.b = self.pattern[self.mode][self.pidx]

    def next(self, clock) :
        if clock == 1 :
            self.nclk += 1
            if self.nclk >= self.N :
                self.pidx += 1
                if self.pidx >= 4 :
                    self.pidx = 0
                if self.mode == 0 :
                    self.b = -2.0 * round(random.random(), 0) + 1.0
                else :
                    self.b = self.pattern[self.mode][self.pidx]
        return self.b

    def reset(self) :
        self.pidx = 0
        self.nclk = 0
        if self.mode == 0 :
            self.b = -2.0 * round(random.random(), 0) + 1.0
        else :
            self.b = self.pattern[self.mode][self.pidx]        

# ------------------------------------------------------------------------------
#                               BPSKModulator
# ------------------------------------------------------------------------------
# This class generate a BPSK modulated signal from
#   - a complex sinusoid carrier and
#   - an PAM encoded binary signal
# The two level of the PAM pulse corresponding to the two logic bit level is
# -1.0 or 1.0. To generate such a signal we could use dirac pulses of 
# -1.0 or 1.0 value for the bits and we use a pulse shaping filter. The simples
# pulse shaping filter has a rectangular impulse response. The mathematical
# background of the used algoritm is that changing the sign of a sine wave is
# equivalent to changing its phase by pi.
# Algorithm:
#   bpsk[n] = bit[n] * carrier[n]
#   bpsk[n] and carrier[n] are complex signals i IQ representation
# Constructor:
#   modulator=BPSKModulator()
#   Parameters passed to the constructor:
#       None
# next method:
#   modulator.next(Icarrier, Qcarrier, pam)
#   Parameters:
#       Icarrier : inphase sample of the carrier.
#                  It should be sample of a sine wave.
#       Qcarrier : quadrature phase sample of the carrier.
#                  It should be sample of a sine wave.
#       pam      : sample of the PAM signal.
#                  It should have -1.0 or 1.0 value.
#   Return:
#       bpsk : sample of the BPSK modulated signal.
# ------------------------------------------------------------------------------

class BPSKModulator:
    def __init__(self):
        pass
    
    def next(self, Icarrier, Qcarrier, pam) :
        Ibpsk = pam * Icarrier
        Qbpsk = pam * Qcarrier
        return Ibpsk, Qbpsk
        
    def reset(self) :
        pass


# ------------------------------------------------------------------------------
#                               ComplexPhaseModulator
# ------------------------------------------------------------------------------
# This class phase modulates a complex sinusoid carrier.
# The modulating signal is a real signal and it represents the modulating phase.
# The complex carries is represenred by the corresponding inphase and quadrature
# phase real signals. The complex output is represented too with the corresponding
# inphase and quadrature phase real signals.
# Algorithm:
#   Iy[n] = Re(exp(j*phi[n])) * Ix[n] - Im(exp(j*phi[n])) * Qx[n]
#   Qy[n] = Re(exp(j*phi[n])) * Qx[n] + Im(exp(j*phi[n])) * Ix[n]
# Constructor:
#   modulator=ComplexPhaseModulator()
#   Parameters passed to the constructor:
#       None
# next method:
#   modulator.next(Ix, Qx, phi)
#   Parameters:
#       Ix : inphase sample of the carrier.
#            It should be sample of a sine wave.
#       Qx : quadrature phase sample of the carrier.
#            It should be sample of a sine wave.
#       phi: sample of the modulating phase.
#   Return:
#       Iy : inphase sample of the modulated signal.
#       Qy : quadrature phase sample of the modulated signal.
# ------------------------------------------------------------------------------

class ComplexPhaseModulator:
    def __init__(self):
        pass
    
    def next(self, Ix, Qx, phi) :
        c = cmath.exp(phi * 1j)
        Iy = c.real * Ix - c.imag * Qx
        Qy = c.real * Qx + c.imag * Ix
        return Iy, Qy
        
    def reset(self) :
        pass

# ------------------------------------------------------------------------------
#                             ComplexNoiseGen
# ------------------------------------------------------------------------------
# The ComplexNoiseGen class provides complex noise signal represented by the
# inphase and quadrature phase real signal pair.
# It uses the random generator from the nympy library.
# The noise
#   - Gaussian noise
#   - The inphase and quadrature phase noise is independent.
#   - The mean value of the inphase and quadrature phase noise is 0.0
#   - The variance of the inphase and quadrature phase noise is 0.5
#   - The variance of the complex noise is 1.0
#   - The samplig frequency is : fs
# Algorithm:
#   We use the noise generator generated by default.rng() from numpy.random to
#   generate noise sample. We use one noise geberator for generating the inphase
#   samples and another one for generating the quadrature phase samples. This
#   way the inphase and quadrature phase noise will be independent.
# Constructor:
#   noise=ComplexNoiseGen()
#   Parameter passed to the constructor:
#       None
# next method
#   noise.next()
#   Parameters:
#       None
#   Return:
#       Iy[n]       : output inphase noise sample
#       Qy[n]       : output quadrature phase noise sample
# ------------------------------------------------------------------------------

class ComplexNoiseGen:
    def __init__(self):
        self.Igen = default_rng()
        self.Qgen = default_rng()
        self.sigma = 1.0 / math.sqrt(2.0) # standard deviation for variance=1/2
        
    def next(self):
        Iy = self.Igen.normal(0.0, self.sigma, 1)
        Qy = self.Igen.normal(0.0, self.sigma, 1)
        return Iy, Qy
    
    def reset(self):
        self.Igen = default_rng()
        self.Qgen = default_rng()
        
# ------------------------------------------------------------------------------
#                           ComplexChannel
# ------------------------------------------------------------------------------
# This class a noisy channel for  complex signal.
# The noise is Gaussian noise.
# Using the channel one could control the signal to noise ratio.
# The input and output signals are complex signal in the inphase and quadrature
# representation.
# Algorithm:
#   Iy[n] = Gs * Isl[n] + Gn * In[n]
#   Qy[n] = Gs * Qsl[n] + Gn * Qn[n]
# Constructor:
#   channel=ComplexChannel8
#   Parameters passed to the constructor:
#       None
# next method
#   channel.next(Ix, Qx, In, Qn, Gs, Gn)
#   Parameters:
#       Ix : Input inphase signal sample
#       Qx : Input quadrature phase signal sample
#       In : Input inphase noise sample
#       Qn : Input quadrature phase noise sample
#       Gs : Signal gain
#       Gn : Noise gain
#            If the variance of the input signal and the input noise is 1 then
#            SNR = 10*log(Gs^2/Gn^2)
#   Return:
#       Iy : Output inphase sample
#       Qy : Output quadrature sample
# ------------------------------------------------------------------------------

class ComplexChannel:
    def __init__(self) :
        pass
    
    def next(self, Ix, Qx, In, Qn, Gs, Gn) :
        Iy = Gs * Ix + Gn * In
        Qy = Gs * Qx + Gn * Qn
        return Iy, Qy

# ------------------------------------------------------------------------------
#                           ComplexMixer
# ------------------------------------------------------------------------------
# This class implements a complex mixer
# The input and output signals are in the inphase and quarature phase form.
# Algorithm:
#   Iy = Is * Ilo - Qs * Qlo
#   Qy = Is * Qlo + Qs * Ilo
# Constructor:
#   mixer=ComplexMixer()
#   Parameters passed to the constructor:
#       None
# next method
#   mixer.(Is, Qs, Ilo, Qlo)
#   Parameters:
#       Is  : Inphase sample of the input signal
#       Qs  : Quadrature phase sample of the input signal
#       Ilo : Inphase sample of the local oscillator signal
#       Qlo : Quadrature phase sample of the local oscillator signal
#   Return:
#       Iy  : Inphase output sample
#       Qy  : Quadrature phase output sample
# ------------------------------------------------------------------------------

class ComplexMixer:
    def __init__(self) :
        pass
    
    def next(self, Is, Qs, Ilo, Qlo) :
        Iy = Is * Ilo - Qs * Qlo
        Qy = Is * Qlo + Qs * Ilo
        return Iy, Qy

    def reset(self) :
        pass

# -----------------------------------------------------------------------------
#                       ComplexMovingAverageDecimator
# -----------------------------------------------------------------------------
# This class implement a decimator with a moving average filter.
# The input and output signals are comlex signal in inphase and quadrature phase
# representation.
# Algorithm:
#   It decimates by N
#   The mathematical model of the decimation:
#                       N-1
#       Iy[k] = (1/N) * Sum(I[k*N+i]) 
#                       i=0
#
#                       N-1
#       Qy[k] = (1/N) * Sum(Qx[k*N+i]) 
#                       i=0
#
#   The cutoff frequency of the moving average filter:
#       fcutoff ~ fs*0,443/N
#       fs : sampling frequency
#
#   In the implementation:
#       For N-1 invocation we update the running summs and output the last valid
#       output signal sample for the signal output and 0 for the valid flag output.
#       For the Nth invocation we finalize the running sum, divide them by N,
#       output the results on the signal output, output 1 on the valid flag output,
#       and clear the running sums and start a new cycle.
#       Initialization in the constructor:
#           i  = 0 
#           Is = 0
#           Qs = 0
#           Iy = 0
#           Qy = 0
#           v  = 0
#       In the next() function
#
#                | Ix         if i == 0
#           Is = |
#                | Is + Ix    otherwise
#
#                | Qx         if i == 0
#           Qs = |
#                | Qs + Qx    otherwise
#
#                | Is/N       if i == N-1
#           Iy = |
#                | Iy         otherwise
#
#                | Qs/N       if i == N-1
#           Qy = |
#                | Qy         otherwise
#
#                | 1          if i == N-1
#           v  = |
#                | 0          otherwise
#
#                | 0          if i == N-1
#           i  = |
#                | i + 1      otherwise
#
# Constructor
#   decimator=ComplexMovingAverageDecimator(N)
#   Parameter passed to the constructor:
#       N : decimation factor
# next method
#   decimator.next(Ix, Qx)
#   Parameters:
#       Ix : current sample of the inphase input signal
#       Qx : current sample of the quadrature phase input signal
#   Return:
#       Iy : current sample of the inphase input signal
#       Qy : current sample of the quadrature phase input signal
#       v  : valid flag
#            v == 1 means new valid output sample is provided
#                   in the Iy, Qy output
#            v == 0 means the last valid output sample is 
#                   provided on Iy, Qy outputs 
# -----------------------------------------------------------------------------

class ComplexMovingAverageDecimator:
    def __init__(self, N) :
        if (N < 0) :
            sys.exit("ComplexMovingAverageDecimator: Wrong N value")
        self.N = N
        self.Nm1 = N - 1
        self.i = 0
        self.Is = 0.0
        self.Qs = 0.0
        self.Iy = 0.0
        self.Qy = 0.0
        self.v = 0
        
    def next(self, Ix, Qx) :
        if (self.i == 0) :
            self.Is = Ix
            self.Qs = Qx
        else :
            self.Is += Ix
            self.Qs += Qx
        if (self.i == self.Nm1) :
            self.Iy = self.Is / self.N
            self.Qy = self.Qs / self.N
            self.v = 1
            self.i = 0
        else :
            self.v = 0
            self.i += 1
        return self.Iy, self.Qy, self.v
        
    def reset(self) :
        self.i = 0
        self.Is = 0.0
        self.Qs = 0.0
        self.Iy = 0.0
        self.Qy = 0.0
        self.v = 0

# -----------------------------------------------------------------------------
#                         Complex FIR Decimator
# -----------------------------------------------------------------------------
# This class implements A FIR decimator
# Algorithm:
#   Mathematical model:
#                    M-1
#           Iy[k] =  Sum(Ix[k*D-M+j]*b[j])
#                    j=0
#
#                    M-1
#           Qy[k] =  Sum(Qx[k*D-M+j]*b[j])
#                    j=0
#       An equivalent impulse response to an 3 stages D=5 CIC decimator filter
#
#           b = [1.0, 3.0, 6.0, 10.0, 15.0, 18.0, 19.0, 18.0, 15.0, 10.0, 6.0,\
#                3.0, 1.0] / 125.0
#
#   Implementation:
#
#   The operation is cyclic, the cycle length is D.
#   In the first D-1 invocation of next()
#       - receive an input I,Q sample pair and stores it in the shift register
#       - outputs the last valid output sample
#       - outputs a zero value valid flag
#   In the Dth invocation of next()
#       - receive an input I,Q sample pair and stores it in the shift register
#       - does the convolution
#       - shits the shift registers by D
#       - outputs the result of the convolution
#       - outputs a 1 value valid flag
#   D : decimation factor
#   M : Impulese response length
#   b : impulse response
#   Is, Qs : shift registers of length M
#   In the constructor:
#       Is[j] = 0.0 for 0 <= j < M
#       Is[j] = 0.0 for 0 <= j < M
#       l = M - D
#       c = 0
#       Iy = 0.0
#       Qy = 0.0
#       v = 0
#   In the next() function
#       Is[l] = Ix
#       Qs[l] = Qx
#       in case of c < D-1
#           c += 1
#           l += 1
#           v = 0
#           output Iy,Qy,v
#       in case of c == D-1
#                 M-1
#           Iy =  Sum(Is[i]*b[i])
#                 i=0
#
#                 M-1
#           Qy =  Sum(Is[i]*b[i])
#                 i=0
#
#           in case of M > D
#               Is[j] = Is[j+D]  for 0 <= j < M-D
#               Is[j] = Is[j+D]  for 0 <= j < M-D
#
#           v = 1
#           l = M-D
#           c = 0
#           output Iy,Qy,v
#
#
#   The diagramm of the algorithm:
#
#   Ix[kN+j] 0 <= j < D
# >>-------------------------------------------------
#                                                   |
#                                                   |M-D+j
#                    0 |                  M-D|      |       | M-1
#                     |--------------------------------------|
#             Is[]    |     shift register of length M       |
#                     |--------------------------------------|
#                      |                                    | Iu[N-1]
#                      |                                    |
#                      |                                    |
#                     |--------------------------------------| 
#                     |           M-1                        |
#                     |           Sum(b[i]*Is[i])            |--------------->>
#                     |           i=0                        |         Iy[k]
#                     |--------------------------------------|          
#                      |                                    |
#                b[0]  |                                    | b[M-1]
#                      |                                    |
#                     |--------------------------------------|
#                     |              b[i] 0 <= i < M         |
#                     |--------------------------------------|
#                      |                                   |
#                b[0]  |                                   | b[M-1]
#                      |                                   |
#                     |--------------------------------------| 
#                     |           M-1                        |
#                     |           Sum(b[i]*Qs[i])            |--------------->>
#                     |           i=0                        |         Iy[k]
#                     |--------------------------------------|          
#                      |                                    |
#                Qs[0] |                                    | Qs[M-1]
#                      |                                    |
#                     |--------------------------------------|
#                     |    shift register of length M        |
#                     |--------------------------------------|
#                   0 |                   M-D|      |       | M-1
#                                                   |M-D+j
#                                                   |
# >>-------------------------------------------------
#   Qx[kN+j] 0 <= j < D
#
# Constructor
#   decimator=ComplexFIRDecimator(b, D)
#   Parameters passed to the constructor:
#       D   : Decimation factor
#       b[] : Impulse response
#             It should be time reversed
#             It should be an numpy array with numpy float type of elements
#             The legth of the array defines M
# next method
#   decimator.next(Ix, Qx)
#   Parameters:
#       Ix : current sample of the inphase input
#       Qx : current sample of the quadrature phase input
#   Return:
#       Iy : current sample of the inphase output
#       Qy : current sample of the quadrature phase output
#       v  : valid flag
#           v == 1 means new valid output sample is provided 
#                  in the Iy, Qy output
#           v == 0 means the last valid output sample is
#                  provided on Iy, Qy outputs 
# -----------------------------------------------------------------------------


class ComplexFIRDecimator:
    def __init__(self, b, D) :
        if D < 0 :
            sys.exit("FIRDecimator: D must be positive")
        self.b = b
        self.D = D
        self.Dm1 = self.D -1
        self.M = len(b)
        if self.D > self.M :
            sys.exit("ComplexFIRDecimator: D is too large")
        self.MmD = self.M - self.D
        self.Iy = 0.0
        self.Qy = 0.0
        self.v = 0
        self.l = self.MmD
        self.c = 0
        self.Is = np.zeros(self.M, np.float)
        self.Qs = np.zeros(self.M, np.float)
    
    def next(self, Ix, Qx) :
        self.Is[self.l] = Ix
        self.Qs[self.l] = Qx
        if self.c < self.Dm1 :
            self.c += 1
            self.l += 1
            self.v = 0
        else :
            self.Iy = 0.0
            self.Qy = 0.0
            for j in range (0, self.M) :
                self.Iy += self.b[j] * self.Is[j]
                self.Qy += self.b[j] * self.Qs[j]
            if self.M > self.D :
                for j in range(0, self.MmD) :
                    self.Is[j] = self.Is[j + self.D]
                    self.Qs[j] = self.Qs[j + self.D]
            self.c = 0
            self.l = self.MmD
            self.v = 1
        return self.Iy, self.Qy, self.v
    
    def reset(self) :
        self.Iy = 0.0
        self.Qy = 0.0
        self.v = 0
        self.l = self.MmD
        self.c = 0
        self.Iu = np.zeros(self.M, np.float)
        self.Qu = np.zeros(self.M, np.float)
        
        

# -----------------------------------------------------------------------------
#                       ComplexVariableMovingAverageFilter
# -----------------------------------------------------------------------------
# This class implement a variable length moving average filter.
# The length of the filter is limited by the maximum length parameter set in
# the constructor.
# For every invocation the filter provides new output samples, so no decimation.
# Algorithm:
#   The mathematical model of the filter:
#       M <= Mmax
#       Ix[k] = 0.0 for -Mmax < k <= -1
#       Qx[k] = 0.0 for -Mmax < k <= -1
#
#                       M-1
#       Iy[n] = (1/M) * Sum(Ix(n-i)) 
#                       i=0
#
#                       M-1
#       Qy[n] = (1/M) * Sum(Qx(n-i)) 
#                       i=0
#
#   The cutoff frequency of the moving average filter:
#       fcutoff ~ fs*0,443/M
#       fs : sampling frequency
#
#   In the implementation:
#       Mmax is the limit of the filter length
#       M is the actual filter length
#
#       Initialization in the constructor:
#           Is[k] = 0.0 for 0 <= k < Mmax
#           Qs[k] = 0.0 for 0 <= k < Mmax
#       In the nex() function:
#           Is[k] = Is[k-1] for k=Mmax-1 ....k=1
#           Qs[k] = Qs[k-1] for k=Mmax-1 ....k=1
#
#           Is[0] = Ix
#           Qs[0] = Qx
#
#                        M-1
#           Iy = (1/M) * Sum Is[i]
#                        i=0
#
#                        M-1
#           Qy = (1/M) * Sum Qs[i]
#                        i=0
# Constructor:ComplexVariableMovingAverageFilter
#   filter=(Mmax)
#   Parameter passed to the constructor:
#       Mmax : Maximum length of the filter
# next method
#   filter.next(Ix, Qx, M)
#   Parameters:
#       M    : current filter length
#              0 < M <= Mmax
#       Ix   : Current sample of the input inphase signal
#       Qx   : Current sample of the input quadrature phase signal
#   Return:
#       Iy   : Current sample of the output inphase signal
#       Qy   : Current sample of the output quadrature phase signal
# -----------------------------------------------------------------------------

class ComplexVariableMovingAverageFilter:
    def __init__(self, Mmax) :
        if (Mmax < 0) :
            sys.exit("ComplexVariableMovingAverageFilter: Wrong Mmax value")
        self.Mmax = Mmax
        self.Is = np.zeros(self.Mmax, np.float)
        self.Qs = np.zeros(self.Mmax, np.float)
        
    def next(self, Ix, Qx, M) :
        if (M < 0) or (M > self.Mmax) :
            sys.exit("ComplexVariableMovingAverageFilter: Wrong M value")
        for i in reversed(range(1, self.Mmax)) :
            self.Is[i] = self.Is[i - 1]
            self.Qs[i] = self.Qs[i - 1]
        self.Is[0] = Ix
        self.Qs[0] = Qx
        Iy = 0.0
        Qy = 0.0
        for i in range(0,M) :
            Iy += self.Is[i]
            Qy += self.Qs[i]
        Iy /= M
        Qy /= M
        return Iy, Qy
        
    def reset(self) :
        for i in range(0,self.Mmax) :
            self.Is[i] = 0.0
            self.Qs[i] = 0.0


# -----------------------------------------------------------------------------
#                            MovingAverageFilter
# -----------------------------------------------------------------------------
# This class implement a fixed length real moving average filter.
# Algorithm:
#   The mathematical model of the filter:
#       N is the length of the filter
#       x[k] = 0.0 for -N < k <= -1
#
#                      N-1
#       y[n] = (1/N) * Sum(x(n-i)) 
#                      i=0
#
#   The cutoff frequency of the moving average filter:
#       fcutoff ~ fs*0,443/M
#       fs : sampling frequency
#
#   This filter could be used to compute the sliding mean of a signal. The mean is
#   using N signal sample.
#
# Implementation:
#   The diagramm of the implementation:
#
#     x[n]            + |---|      +|---|        u[n]     |---|    y[n]
# >>------------------->| + |>----->| + |---------------->| * |>-------->>
#          |            |---|       |---|            |    |---|
#          |              |-          |+             |      |
#          |              |           |              |      |
#          |    |----|    |           |    |----|    |     1/N
#          |>---|z^-N|--->|           |<---|z^-1|<---|
#               |----|  x[n-N]      u[n-1] |----|
#
#   The delay is implemented by using a circular buffer.
#
#   In the constructor:
#       i = 0
#       s[k] = 0.0  for 0 <= k < N
#   In the run function:
#       u += x - s[i]
#       y = s / N
#       s[i] = x
#           | 0    if i==N-1
#       i = |
#           | i+1  otherwise
# Constructor
#   filter=MovingAverageFilter(N)
#   Parameter passed to the constructor:
#       N : The length of the filter
# next method
#   filter.next(x)
#   Parameter:
#       x : current sample of the input signal
#   Return:
#       y : Currenr sample of the output signal
# -----------------------------------------------------------------------------

class MovingAverageFilter :
    def __init__(self, N) :
        if (N < 0) :
            sys.exit("MovingAverageFilter: Wrong N value")
        self.N = N
        self.s = np.zeros(self.N, np.float)
        self.u = 0.0
        self.i = 0
        
    def next(self, x) :
        u += x - s[self.i]
        y = u / N
        s[self.i] = x
        if self.i == N -1 :
            self.i = 0
        else :
            self.i += 1
        return y
        
    def reset(self) :
        self.u = 0.0
        self.i = 0
        for j in range(0,self.N) :
            self.s[j] = 0.0

    
# -----------------------------------------------------------------------------
#                           ComplexMovingAverageFilter
# -----------------------------------------------------------------------------
# This class implement a fixed length complex moving average filter.
# Algorithm:
#   The mathematical model of the filter:
#       N is the length of the filter
#       Ix[k] = 0.0 for -N < k <= -1
#       Qx[k] = 0.0 for -N < k <= -1
#
#                       N-1
#       Iy[n] = (1/N) * Sum(Ix(n-i)) 
#                       i=0
#
#                       N-1
#       Qy[n] = (1/N) * Sum(Qx(n-i)) 
#                       i=0
#
#   The cutoff frequency of the moving average filter:
#       fcutoff ~ fs*0,443/M
#       fs : sampling frequency
#
#   This filter could be used to compute the sliding mean of a signal. The mean is
#   using N signal sample.
#
# Implementation:
#   The diagramm of the implementation:
#
#    Ix[n]            + |---|      +|---|       Iu[n]     |---|   Iy[n]
# >>------------------->| + |>----->| + |---------------->| * |>-------->>
#          |            |---|       |---|            |    |---|
#          |              |-          |+             |      |
#          |              |           |              |      |
#          |    |----|    |           |    |----|    |     1/N
#          |>---|z^-N|--->|           |<---|z^-1|<---|
#               |----| Ix[n-N]     Iu[n-1] |----|
#
#
#    Qx[n]            + |---|      +|---|       Qu[n]     |---|   Qy[n]
# >>------------------->| + |>----->| + |---------------->| * |>-------->>
#          |            |---|       |---|            |    |---|
#          |              |-          |+             |      |
#          |              |           |              |      |
#          |    |----|    |           |    |----|    |     1/N
#          |>---|z^-N|--->|           |<---|z^-1|<---|
#               |----| Qx[n-N]     Qu[n-1] |----|
#
#   The delay is implemented by using a circular buffer.
#
#   In the constructor:
#       i = 0
#       Is[k] = 0.0  for 0 <= k < N
#       Qs[k] = 0.0  for 0 <= k < N
#   In the run function:
#       Iu += Ix - Is[i]
#       Qu += Qx - Qs[i]
#       Iy = Is / N
#       Qy = Qs / N
#       Is[i] = Ix
#       Qs[i] = Qx
#           | 0    if i==N-1
#       i = |
#           | i+1  otherwise
# Constructor
#   filter=ComplexMovingAverageFilter(N)
#   Parameter passed to the constructor:
#       N : The length of the filter
# next method
#   Parameters:
#       Ix : current sample of the input inphase signal
#       Qx : current sample of the input quadrature phase signal
#   Return:
#       Iy : Currenr sample of the output inphase signal
#       Qy : Currenr sample of the output quadrature phase signal
# -----------------------------------------------------------------------------

class ComplexMovingAverageFilter :
    def __init__(self, N) :
        if (N < 0) :
            sys.exit("MovingAverageFilter: Wrong N value")
        self.N = N
        self.Is = np.zeros(self.N, np.float)
        self.Qs = np.zeros(self.N, np.float)
        self.Iu = 0.0
        self.Qu = 0.0
        self.i = 0
        self.valid = 0
        
    def next(self, Ix, Qx) :
        self.Iu += Ix - self.Is[self.i]
        self.Qu += Qx - self.Qs[self.i]
        Iy = self.Iu / self.N
        Qy = self.Qu / self.N
        self.Is[self.i] = Ix
        self.Qs[self.i] = Qx
        if self.i == self.N -1 :
            self.i = 0
            self.valid = 1
        else :
            self.i += 1
        return Iy, Qy, self.valid
        
    def reset(self) :
        self.i = 0
        self.valid = 0
        self.Iu = 0.0
        self.Qu = 0.0
        for j in range(0,self.N) :
            self.Is[j] = 0.0
            self.Qs[j] = 0.0

    
# -----------------------------------------------------------------------------
#                            GatedIntegrator
# -----------------------------------------------------------------------------
# This class implements a special integrator.
# It is realy a moving average filter which has not fixed averaging length. The
# number of samples to be averaged is determined by an external clock.
# The result is the sum of the input samples divided by the number of samples added
# The clock pulse is 1 sample long.
# A clock pulse:
#   - ends the current integration cycle 
#   - starts the next integration cycle
# The integrator puts the result to the output when a cycle has been finished and
# outputs the same valu until the next result wiil be available.
# The integrator signals a new result by a one sample long pulse on the valid output.
# The integrator also outputs of the length of the integration cycle.
# If each of the cyle is N sample long then the integrator does a decimation by N.
#
# Algorithm:
#   The clock signal: c
#   Initialization in the constructor:
#       k = 1
#       u = 0.0
#       y = 0.0
#       N = 0
#       v = 0
#
#   In the next() function:       
#       The output:
#               | u/k    if c==1
#           y = |
#               | y      otherwise
#
#               | k      if c==1
#           N = |
#               | N      otherwise
#
#           v = c
#
#       The integration:
#               | x      if c==1
#           u = |
#               | u+x    otherwise
#
#       The integration length count signal:
#               | 1,     if c==1
#           k = |
#               | k+1,   otherwise
#
#   The diagramm of the algorithm:
#
#  x[n]     |---|            u[n]               |---| v[n]|----|     |-----|
# >>------->| + |>----------------------------->| / |>--->|z^-1|>--->|- \  |
#           |---|                          |    |---|     |----|     |1  \ |
#             |                            |      |                  |    -|>---------------------->>  y[n]
#             |      |-----|     |----|    |      |                  |0    |                |
#             |      |  / -|<---<|z^-1|<--<|      |              |>->|-    |                |
#             |      | /  0|     |----|           |              |   |-----|                |
#             |<----<|-    |                      |              |      |        |----|     |
#                    |    1|                      |              |<-----)-------<|z^-1|<---<|
#                    |    -|<--- 0.0              |                     |        |----|
#                    |-----|                      |                    c[n]
#                       |                         |
#                       |                         |
#                      c[n]                       |
#           |---|               k[n]              |       |----|     |-----|
#     1 >-->| + |>------------------------------->|>----->|z^-1|>--->|- \  |
#           |---|                          |              |----|     |1  \ |
#             |                            |                         |    -|>---------------------->> N[n]
#             |      |-----|     |----|    |                         |0    |                |
#             |      |  / -|<--->|z^-1|<--<|                     |>->|-    |                |
#             |      | /  0|     |----|                          |   |-----|                |
#             |<----<|-    |                                     |      |        |----|     |
#                    |    1|                                     |<-----)-------<|z^-1|<---<|
#                    |    -|<--- 0                                      |        |----|
#                    |-----|                                           c[n]
#                       |                          
#                       |                          
#    c[n]              c[n]                        
# >>---->                                                                                  c[n] >-->> v[n]
#
# Constructor
#   integrator=GatedIntegrator()
# Parameter passed to the constructor:
#   None
# next method
#   integrator.next(c, x)
#   Parameters:
#       x : Current sample of the inputs signal
#       c : Current sample of the clock signal
#   Return:
#       y : the last valid output sample
#       N : the length of the last valid integration in samples
#       v : valid flag
#           v == 1 means new valid output sample is provided 
#                  on the y output
#           v == 0 means the last valid output sample is
#                  provided on y outputs
# -----------------------------------------------------------------------------

class GatedIntegrator:
    def __init__(self) :
        self.k = 1
        self.u = 0.0
        self.y = 0.0
        self.N = 0
        self.v = 0
        
    def next(self, c, x) :
        if (c == 1) :
            self.y = self.u / self.k
            self.N = self.k
            self.u = x
            self.k = 1
        else :
            self.u += x
            self.k += 1
        self.v = c
        return self.y, self.N, self.v
        
    def reset(self) :
        self.k = 1
        self.u = 0.0
        self.y = 0.0
        self.N = 0
        self.v = 0
    
# -----------------------------------------------------------------------------
#                         Phase detector
# -----------------------------------------------------------------------------
# This class implements a linear phase detector
# It compute the angle of the input complex number.
# The input is represented by a inphase and quadrature phase real number pair.
# The standard mathematical functions are able to cumpute only the principal
# value of the angle. If we suppose, that the input is a rotating complex
# vector we could unwrap the wrapping occurs in the output of the standard
# mathematical functions at pi or -pi. So we could determine the tru linear
# phase difference.
# Algorithm
#   The input of the detector is a complex vector in IQ reprezentation
#       Ix[n] = cos(phi[n])
#       Qx[n] = sin(phi[n])
#   We create a complex variable
#           s = Ix + jQx
#   Depending on the the detector type (pdtype) we do the following
#       In case of carrier mode:
#           s = s
#       In case of BPSK mode
#           s = s * s
#           So we have 2 * phi angle, this maps the BPSK modulation angels
#           to 0
#       In case of QPSK mode
#           s = s * s * s * s
#           so we have 4 * phi angle, which maps all the QPSK modulation
#           angel to pi
#   We compute the principal value of the phase angle
#       phi = phase(s)
#   The detector algorithm unwrap the wraping of the angle of the complex vector
#   by observing the angle difference of two consecutive samples
#
#          pi |   /|                  pi |    |\
#             |  + |x[n-1]               |    | +x[n]
#             | /  |                     |    |  \
#             |/   |     /              \|    |   \
#       ------/----|----/-----     ------\----|----\-----
#            /|    |   /                 |\   |     \
#             |    |  /                  | \  |
#             |    | + x[n]              |  + |x[n-1]
#         -pi |    |/                -pi |   \|
#
#
#           increasing angle          decreasing angle
#           wrapping at pi            wrapping at -pi
#           phi[n-1]>pi/2             phi[n-1]<-pi/2
#           phi[n]<-pi/2              phi[n]>pi/2
#           phi[n]-phi[n-1]<-mu*pi    phi[n]-phi[n-1]>mu*pi
#           unwrap+=2*pi              unwrap-=2*pi
#   Above mu>=1 is used to lower the false unwrapping du to noise. It depends
#   on the ration of the sampling rate and the signal frequency. In case of
#   low signal frequency the change of the principal value ot wrapping is
#   close to 2*pi, in case of high signal frequency the change of the principal
#   value is close to pi.
#   The return value:
#       in principal angle mode:
#           phi, in case of carrier type detector
#           phi / 2, in case of BPSK type detector
#           phi / 4 - pi / 4, in case of QPSK type detector
#               here we have to substract pi/4 becouse the modulation angles
#               have been maped to pi. 
#       in linear mode: phi + unwrap
#   The algorithm works only if the sampling rate is higher than the critical
#   sampling rate.
# Constructor
#   detector=PhaseDetector(pdtype=0, mu=1.0)
#   Parameter passed to the constructor:
#       pdtype:  type of the phase detector
#                0 : Carrier
#                1 : BPSK
#                2 : QPSK
#                default is Carrier
#       mu : safety factor
#            default value: 1.0
# next method
#   detector.next(Ix, Qx, pdmode=0)
#   Parameters:
#       Ix   :  current sample of the inphase input signal
#       Qx   :  current sample of the quadrature phase input signal
#       pdmode :mode of angle computation
#               0 - principal value mode
#               1 - linear mode
#               The default value is 0.
#               It is a named parameter.
#   Return:
#       e  : phase of the input complex vextor
# In case of very low SNR the principal value mode is preferred.
# -----------------------------------------------------------------------------

class PhaseDetector:
    def __init__(self, pdtype=0, mu=1.0):
        if not ((pdtype == 0) or (pdtype == 1) or (pdtype == 2)) :
            sys.exit("PhaseDetector: invalid pdtype value")
        self.previous = 0.0
        self.unwrap = 0.0
        self.pdtype = pdtype
        self.mu = mu
        self.pi = math.pi
        self.halfpi = self.pi / 2.0
        self.quaterpi = math.pi / 4.0
        if (pdtype == 0) :
            self.d = 1.0
        elif (pdtype == 1) :
            self.d = 2.0
        else :
            self.d = 4.0
    
    def next(self, Ix, Qx, pdmode=0):
        e = 0.0
        s = complex(Ix, Qx)
        if (self.pdtype == 1) or (self.pdtype == 2) :
            s = s * s
        if (self.pdtype == 2) :
            s = s * s
        phi = cmath.phase(s)
        if (self.previous > self.halfpi) and (phi < -self.halfpi) and (phi - self.previous) < self.mu * self.pi :
            self.unwrap += 2.0 * math.pi
        elif (self.previous > -self.halfpi) and (phi > self.halfpi) and (phi - self.previous) > self.mu * self.pi :
            self.unwrap -= 2.0 * math.pi
#        if (phi - self.previous) > self.mu * math.pi :
#            self.unwrap -= 2.0 * math.pi
#        elif (phi - self.previous) < -self.mu * math.pi :
#            self.unwrap += 2.0 * math.pi
        self.previous = phi
        if pdmode == 0 :
            e = phi + self.unwrap
        else :
            e = phi
        e /= self.d
        if self.pdtype == 2 :
            e -= self.quarterpi
        return e
    
    def reset(self):
        self.previous = 0.0
        self.unwrap = 0.0


    
# -----------------------------------------------------------------------------
#                         Edge Phase detector
# -----------------------------------------------------------------------------
# This class implements a phase detector which detect the edge position of a
# square wave signal relative to the edge of the reference square wave
# signal.
# Algorithm:
#   The square wave signal
#            |  A, if      nT+dT <= t < n*T+T/2+dT
#     x(t) = |
#           | -A, if n*T+T/2+dT <= t < (n+1)*T+dT
#     T is the period of the signal
#     dT is the time shift relative to the reference signal
#    The reference signal is a square wave signal
#            |  B, if         mTr <= t < m*Tr + Tr/2
#     y(t) = |
#            | -B, if m*Tr + Tr/2 <= t < (m+1)*Tr
#     Tr is the period of the reference signal and Tr~T
#
#   The algoritm's'inputs:
#              (m+1)*Tr
#      Ix[m] : integral(x(t))
#               m*Tr
#
#              (m+1)Tr-Tr/2
#      Qx[m] : integral(x(t))
#              m*Tr-Tr/2
#
#    Computing the error:
#             | -sign(I[m-1])*Q if sign(I[m-1]) != sign(I[m])
#      e[n] = | 
#             | 0 otherwise
#     IMPORTANT: If there is no tranzition edge between the two consecutive
#                interwall the error will bw zero.
#                To lessen the effect of noise in the implementation we set
#                certain signal to 0 if they are close to zero. We use this
#                in computing sign of the signal samples. The sign will be
#                0 if the sample is close to zero. This will result a
#                zero error instead of a wrong error.
# This algorithm works for other kind of signals, for example dirac delta
# trains shaped by rised cosine or root rised cosine filter.
# 
# Constructor
#   detector=EdgePhaseDetector(threshold)
#   Parameters passed to the constructor:
#       threshold : threshold valu around zero
# next method
#   detector.next(I, Q)
#   Parameters:
#       I : current sample of the inphase integral
#       Q : current sample of the quadrature phase integral
# Output:
#   e  : phase error
# -----------------------------------------------------------------------------

class EdgePhaseDetector :
    def __init__(self, threshold) :
        self.threshold = threshold
        self.I = 0.0
    
    def next(self, I, Q) :
        if (abs(Q) <= self.threshold) :
            sQ = 0
        elif Q > 0 :
            sQ = 1
        else :
            sQ = -1
        
        if (abs(self.I) <= self.threshold) :
            sIb = sQ
        elif self.I > 0 :
            sIb = 1
        else :
            sIb = -1
        
        if (abs(I) <= self.threshold) :
            sIf = -sQ
        elif I > 0 :
            sIf = 1
        else :
            sIf = -1
        
        if (sIf != sIb) :
            e = -sIb * Q
        else :
            e = 0
        self.I = I
        
        e *= math.pi # to correct the detector gain
        
        return e
    
    def reset(self) :
        self.I = 0.0
    
    
# -----------------------------------------------------------------------------
#                      FFTFrequencyDetector
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------

class FFTFrequencyDetector :
    def __init__(self, fs, N, threshold) :
        if N <= 0:
            sys.exit("FFTFrequencyDetector: zero or egative N value")
        self.fs = fs
        self.N = N
        self.Nm1 = N - 1
        self.threshold = threshold
        self.bw = self.fs / self.N
        self.Is = np.zeros(self.N, np.float)
        self.Qs = np.zeros(self.N, np.float)
        self.l = 0
        self.v = 0
        self.f = 0
        self.pn = 0.0
        self.ps = 0.0
        self.cvalid = 0
        self.k = 0
    
    def next(self, Ix, Qx) :
        self.Is[self.l] = Ix
        self.Qs[self.l] = Qx
        k = 0
        if self.l < self.Nm1 :
            self.l +=1
            self.v = 0
        else :
            self.l = 0
            self.v = 1
            absX = (1.0 / self.N) * np.abs(np.fft.fftshift(np.fft.fft(self.Is + 1j * self.Qs)))
            PX = absX * absX
            self.pn = 0.0
            self.ps = 0.0
            self.cvalid = 0
            for i in range(0, self.N) :
                if PX[i] < 1e-20 :
                    PX[i] = 1e-20
                self.pn += PX[i]
                if PX[i] > self.ps :
                    self.ps = PX[i]
                    k = i
            # correct the signal power with the power leaked into the neighboring bins
            if (k == 0) :
                self.ps += PX[1]
            elif (k == self.Nm1) :
                self.ps += PX[self.Nm1 - 1]
            else :
                self.ps += PX[k-1] + PX[k + 1]
            self.pn -= self.ps  # substract signal power from the total noise power
            self.pn /= self.N   # average power represent the noise power in a bin
            if self.ps > self.threshold :
                self.cvalid = 1
            k -= self.N / 2
            self.f = k * self.bw
            self.k = k
            
        return self.v, self.f, self.k, self.pn, self.ps, self.cvalid
    
    def reset(self) :
        self.l = 0
        self.v = 0
        self.f = 0
        self.pn = 0.0
        self.ps = 0.0
        self.cp = 0
        for j in range(0,self.N) :
            self.Is[j] = 0.0
            self.Qs[j] = 0.0



# -----------------------------------------------------------------------------
#                         Type2Loopfilter
# -----------------------------------------------------------------------------
# This class implements an type2 loop filter with parameter changing 
# The code managing the parameter change ensures a smooth output at paramater
# changing.
# Mathematical modell:
#   Parameters:
#       BL    : SSB noise bandwidth
#       phiPM : Phase margin (in degree)
#       M     : Averiging length for computing the average value of the output
#   Derived parameters:
#       Kp    : the gain of the proportional part
#       Ki    : the gain of the integrator part
#   Parameter calculation:
#       ro = tg(phiPM)
#       Kp = 4*BL*(ro/(1+ro))
#       w0 = K/ro
#       Ki = w0*Ts
#   The transfer function of the filter using Backward Euler integration
#
#       L(z) = (1+Ki/(1-1/z))*Kp
#
#   Loop parameter change:
#       Ki could be change without any transient effect
#       If we change Kp we would have a jump at the output eve the input is unchanged.
#       So instead of using K as a
#           y = K * x
#       transfer function and always use a transfer curve going through the origo
#       we locally change the slope of the transfer curve as shown in the following
#       diagramm
#
#         y  ^
#            |                 /          .
#            |             K0 /         . K1
#            |               /        .
#            |              /       .
#            |             /      .
#            |            /     .
#            |__________ /___ .
#         y1 |          /   .|
#            |         /  .  |
#            |        / .    |
#         y0 |_______/.______|_____________
#            |     ./|       |
#            |   . / |       |
#            | .  /  |       |
#            |   /   |       |
#            |  /    |       |
#            | /     |       |
#            |/      |       |
#            --------|-------|-----------------> x
#                    x0      x1
#
#            y[n] = K0 * x[n]  n <= n0
#            y[n0+i] = K1 * (x[n0+i] - x[n0]) + y[n0]  i>0 
#       The reason this works is that we change the parameters after the PLL achived lock
#       in order to reduce the PLL bandwidth and reduce the noise.
#       In the implementation we use average(y) in place of y[n] becouse when x is noisy
#       we still conuld have some jump at the output if we use y[n].
#
#   Implementation
#
#       ei[-1] = 0.0
#       Kp[-1] = 0.0
#       ccorr[-1] = 0.0
#       c[i] = 0.0   if -M<=i<0
#       cav[-1] = 0.0
#       ci[-1] = 0.0
#
#       ei[n] = ei[n-1] + Ki[n]*e[n]
#       cr[n] = ei[n] + e[n]
#
#
#              | 1    if Kp[n] != Kp[n-1]
#       s[n] = |
#              | 0    otherwise
#
#                  | cav[n-1] - Kp[n] * cr[n-1]       if s[n] = 1
#       ccorr[n] = |
#                  | ccorr[n-1]                       otherwise
#
#       c[n] = Kp * cr[n] + ccorr[n]
#
#       ci[n] = ci[n-1] + c[n] - c[n-M]
#       cav[n] = ci[n] / M
#       
#    The block diagramm of the filter with Backward Euler integration:
#
#                     Ki
#                     |
#   e[n]            |---|     |---|              ei[n]       |---| 
#   >>------------->| X |>--->| + |>------------------------>| + |>-->|
#     |             |---|     |---|                |         |---|    |
#     |                         |                  |           |      |
#     |                         |       |---|      |           |      |
#     |                         |<-----<| -1|<----<|           |      |
#     |                          ei[n-1]|z  |                  |      |
#     |                                 |---|                  |      |
#     |                                                        |      |
#     |>------------------------------------------------------>|      |
#                                    e[n]                             |
#                                                                     |
#                                                                     |
#                                    cr[n]                            |
#     |<-------------------------------------------------------------<|
#     |
#     |                             
#     |      c[n]     + |---|      +|---|       ci[n]     |---| cav[n]
#     |    |>---------->| + |>----->| + |---------------->| * |>------>
#     |    |            |---|       |---|            |    |---|
#     |    |              |-          |+             |      |
#     |    |              |           |              |      |
#     |    |    |----|    |           |    |----|    |     1/N
#     |    |>---|z^-M|--->|           |<---|z^-1|<---|
#     |    |    |----| cn[n-N]     ci[n-1] |----|
#     |    |
#     |    |<----------------------------------------------------<|
#     |                                                           |
#     |                                                           |
#     |    cr[n]     |---|                             |---|      |      c[n]
#     |>------------>| X |>--------------------------->| + |>---->|>----------------->>
#     |              |---|                             |---|      
#     |                |                                 |
#     |                |    |----| K[n-1]|---|           |
#     |                |>---|z^-1|>----->|   |           |
#     |   Kp[n]        |    |----|       |   | s[n]      |
#  >>-)--------------->|                 |!= |>-->|      |
#     |                |      Kp[n]      |   |    |      |
#     |                |>--------------->|   |    |      |
#     |                |                 |---|    |      |
#     |                |                          |      | ccorr[n]
#     | cr[n]          |                          |      |
#     |                |                          |      |  
#     |    |----|    |---|    |---|            |----|    |
#     |>-->|z^-1|>-->| X |>-->| + |>---------->|1\  |    |
#          |----|    |---|    |---|            |  \ |    |
#                               |              |   \|>-->|
#                               |              |    |    |
#   cav[n] |----|    cav[n-1]   |        |>--->|0   |    |
#   >----->|z^-1|>--------------|        |     |----|    |
#          |----|                        |               |
#                                        |     |----|    |
#                                        |<---<|z^-1|<--<|
#                                              |----|
#
# Constructor
#   filter=Type2LoopFilter(Ts, M)
#   Parameters passed to the constructor
#       Ts : sampling time
#       M  : length of the averaging
# next method
#   filter.next(e, Kp, Ki)
#   Parameters:
#       e  : current error sample
#       Kp : current Kp value
#       Ki : current Ki value
#   Return:
#       c  : current sample of the filter output                                         
#
# -----------------------------------------------------------------------------


class Type2LoopFilter:
    def __init__(self, Ts, M):
        if M <= 0:
            sys.exit("Type2LoopFilter: zero or egative M value")
        self.M = M
        self.Ts = Ts
        self.ci = 0.0
        self.ei = 0.0
        self.Kp = 0.0
        self.cav = 0.0
        self.ccorr = 0.0   
        self.cr = 0.0   # the raw output of the filter     
        self.c = np.zeros(self.M, np.float)
        self.i = 0
    
    def next(self, e, Kp, Ki):
        if (Kp != self.Kp) and (self.Kp != 0.0) :
            self.ccorr = self.cav - Kp * self.cr
        self.Kp = Kp
        self.ei += Ki * e
        self.cr = self.ei + e
        c = self.Kp * self.cr + self.ccorr
        self.ci += c - self.c[self.i]
        self.cav = self.ci / self.M
        self.c[self.i] = c
        if self.i == self.M - 1 :
            self.i = 0
        else :
            self.i += 1
        return c
    
    def reset(self):
        self.ci = 0.0
        self.ei = 0.0
        self.Kp = 0
        self.cav = 0.0
        self.ccorr = 0.0
        self.cr = 0.0
        self.i = 0
        for j in range(0,self.M) :
            self.c[j] = 0.0
    

# -----------------------------------------------------------------------------
#                         Type3Loopfilter
# -----------------------------------------------------------------------------
# This class implements an type3 loop filter with parameter changing 
# The code managing the parameter change ensures a smooth output at paramater
# changing.
# Mathematical modell:
#   Parameters:
#       BL    : SSB noise bandwidth
#       phiPM : Phase margin (in degree)
#       M     : Averiging length for computing the average value of the output
#   Derived parameters:
#       Kp    : the gain of the proportional part
#       Ki    : the gain of the integrator part
#   Parameter calculation:
#       ro = tg((phiPM+90°)/2)
#       Kp = 4*BL*((2*ro-1)/(2*ro+3))
#       w0 = K/ro
#       Ki = w0*Ts
#   The transfer function of the filter using Backward Euler integration
#
#       L(z) = (1+Ki/(1-1/z))*(1+Ki/(1-1/z))*Kp
#
#   Loop parameter change:
#       Ki could be change without any transient effect
#       If we change Kp we would have a jump at the output eve the input is unchanged.
#       So instead of using K as a
#           y = K * x
#       transfer function and always use a transfer curve going through the origo
#       we locally change the slope of the transfer curve as shown in the following
#       diagramm
#
#         y  ^
#            |                 /          .
#            |             K0 /         . K1
#            |               /        .
#            |              /       .
#            |             /      .
#            |            /     .
#            |__________ /___ .
#         y1 |          /   .|
#            |         /  .  |
#            |        / .    |
#         y0 |_______/.______|_____________
#            |     ./|       |
#            |   . / |       |
#            | .  /  |       |
#            |   /   |       |
#            |  /    |       |
#            | /     |       |
#            |/      |       |
#            --------|-------|-----------------> x
#                    x0      x1
#
#            y[n] = K0 * x[n]  n <= n0
#            y[n0+i] = K1 * (x[n0+i] - x[n0]) + y[n0]  i>0 
#       The reason this works is that we change the parameters after the PLL achived lock
#       in order to reduce the PLL bandwidth and reduce the noise.
#       In the implementation we use average(y) in place of y[n] becouse when x is noisy
#       we still conuld have some jump at the output if we use y[n].
#
#   Implementation
#
#       ei[-1] = 0.0
#       Kp[-1] = 0.0
#       ccorr[-1] = 0.0
#       c[i] = 0.0   if -M<=i<0
#       cav[-1] = 0.0
#       ci[-1] = 0.0
#
#       ei[n] = ei[n-1] + Ki[n]*e[n]
#       cr[n] = ei[n] + e[n]
#
#
#              | 1    if Kp[n] != Kp[n-1]
#       s[n] = |
#              | 0    otherwise
#
#                  | cav[n-1] - Kp[n] * cr[n-1]       if s[n] = 1
#       ccorr[n] = |
#                  | ccorr[n-1]                       otherwise
#
#       c[n] = Kp * cr[n] + ccorr[n]
#
#       ci[n] = ci[n-1] + c[n] - c[n-M]
#       cav[n] = ci[n] / M
#       
#    The block diagramm of the filter with Backward Euler integration:
#
#                     Ki
#                     |
#   e[n]            |---|     |---|              ei1[n]      |---| 
#   >>------------->| X |>--->| + |>------------------------>| + |>-->|
#     |             |---|     |---|                |         |---|    |
#     |                         |                  |           |      |
#     |                         |       |---|      |           |      |
#     |                         |<-----<| -1|<----<|           |      |
#     |                          ei[n-1]|z  |                  |      |
#     |                                 |---|                  |      |
#     |                                                        |      |
#     |>------------------------------------------------------>|      |
#                                    e[n]                             |
#                                                                     |
#                                                                     |
#                                                                     |
#     |<-------------------------------------------------------------<|
#     |
#     |                             
#     |               Ki
#     |               |
#     |             |---|     |---|              ei2[n]      |---| 
#     |>----------->| X |>--->| + |>------------------------>| + |>-->|
#     |             |---|     |---|                |         |---|    |
#     |                         |                  |           |      |
#     |                         |       |---|      |           |      |
#     |                         |<-----<| -1|<----<|           |      |
#     |                         ei2[n-1]|z  |                  |      |
#     |                                 |---|                  |      |
#     |                                                        |      |
#     |>------------------------------------------------------>|      |
#                                                                     |
#                                                                     |
#                                                                     |
#                                    cr[n]                            |
#     |<-------------------------------------------------------------<|
#     |
#     |                             
#     |      c[n]     + |---|      +|---|       ci[n]     |---| cav[n]
#     |    |>---------->| + |>----->| + |---------------->| * |>------>
#     |    |            |---|       |---|            |    |---|
#     |    |              |-          |+             |      |
#     |    |              |           |              |      |
#     |    |    |----|    |           |    |----|    |     1/N
#     |    |>---|z^-M|--->|           |<---|z^-1|<---|
#     |    |    |----| cn[n-N]     ci[n-1] |----|
#     |    |
#     |    |<----------------------------------------------------<|
#     |                                                           |
#     |                                                           |
#     |    cr[n]     |---|                             |---|      |      c[n]
#     |>------------>| X |>--------------------------->| + |>---->|>----------------->>
#     |              |---|                             |---|      
#     |                |                                 |
#     |                |    |----| K[n-1]|---|           |
#     |                |>---|z^-1|>----->|   |           |
#     |   Kp[n]        |    |----|       |   | s[n]      |
#  >>-)--------------->|                 |!= |>-->|      |
#     |                |      Kp[n]      |   |    |      |
#     |                |>--------------->|   |    |      |
#     |                |                 |---|    |      |
#     |                |                          |      | ccorr[n]
#     | cr[n]          |                          |      |
#     |                |                          |      |  
#     |    |----|    |---|    |---|            |----|    |
#     |>-->|z^-1|>-->| X |>-->| + |>---------->|1\  |    |
#          |----|    |---|    |---|            |  \ |    |
#                               |              |   \|>-->|
#                               |              |    |    |
#   cav[n] |----|    cav[n-1]   |        |>--->|0   |    |
#   >----->|z^-1|>--------------|        |     |----|    |
#          |----|                        |               |
#                                        |     |----|    |
#                                        |<---<|z^-1|<--<|
#                                              |----|
#
# Constructor
#   filter=Type3LoopFilter(Ts, M)
#   Parameters passed to the constructor
#       Ts : sampling time
#       M  : length of the averaging
# next method
#   filter.next(e, Kp, Ki)
#   Parameters:
#       e  : current error sample
#       Kp : current Kp value
#       Ki : current Ki value
#   Return:
#       c  : current sample of the filter output                                         
#
# -----------------------------------------------------------------------------


class Type3LoopFilter:
    def __init__(self, Ts, M):
        if M <= 0:
            sys.exit("Type2LoopFilter: zero or negative M value")
        self.M = M
        self.Ts = Ts
        self.ci = 0.0
        self.ei1 = 0.0
        self.ei2 = 0.0
        self.Kp = 0
        self.cav = 0.0
        self.ccorr = 0.0
        self.cr = 0.0
        self.c = np.zeros(self.M, np.float)
        self.i = 0
    
    def next(self, e, Kp, Ki):
        if Kp != self.Kp :
            # at the first invocation or after reset() ccorr will be zero
            self.ccorr = self.cav - Kp * self.cr
        self.Kp = Kp
        self.ei1 += Ki * e
        self.ei2 += Ki * (self.ei1 + e)
        self.cr = self.ei2 + (self.ei1 + e)
        c = self.Kp * self.cr + self.ccorr
        self.ci += c - self.c[self.i]
        self.cav = self.ci / self.M
        self.c[self.i] = c
        if self.i == self.M - 1 :
            self.i = 0
        else :
            self.i += 1
        return c
    
    def reset(self):
        self.ci = 0.0
        self.ei1 = 0.0
        self.ei2 = 0.0
        self.Kp = 0
        self.cav = 0.0
        self.ccorr = 0.0
        self.cr = 0.0
        self.i = 0
        for j in range(0,self.M) :
            self.c[j] = 0.0
    


