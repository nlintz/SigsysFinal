import sys
import os
import pylab
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from Analysis import *


def test_noHands():
    S = SongAnalyzer('./samples/nohands.wav')
    S.plotSpectrogram()


def test_a440():
    S = SongAnalyzer('./samples/a440.wav')
    S.plotSpectrogram()


def test_a440_fft():
    S = SongAnalyzer('./samples/a440.wav')
    S.plotFFT()


def test_downsampling(n):
    S = SongAnalyzer('./samples/nohands.wav')
    S.plotFFTDownSample(n)


def test_partial_spectrum(start=30, stop=40):
    S = SongAnalyzer('./samples/nohands.wav')
    S.plotPartialSpectrogram(start, stop)


def test_partial_fft(start=30, stop=31):
    S = SongAnalyzer('./samples/nohands.wav')
    S.plotPartialFFT(start, stop)


def test_wavwrite():
    S = SongModifier('./samples/nohands.wav')
    S.writeWAV(S.data, './processed/nohands.wav')


def test_writeInverseFFT(startTime, stopTime):
    S = SongModifier('./samples/nohands.wav')
    AMP = S.partialRealFFT(startTime, stopTime)[1]
    IFFT = S.inverseRealFFT(AMP)
    S.writeWAV(IFFT, './processed/nohandsTEST.wav')


def testLPF():
    S = SongModifier('./samples/nohands.wav')
    FFT = S.partialRealFFT(10, 20)
    lpf = S.lowPassFFT(FFT, 22000)
    IFFT = S.inverseRealFFT(lpf)
    S.writeWAV(IFFT, './processed/nofiltertest.wav')


def testPlotLPF(start=30, stop=31, cutoff=100):
    S = SongModifier('./samples/nohands.wav')
    FFT = S.partialRealFFT(start, stop)
    lpf = S.lowPassFFT(FFT, cutoff)
    IFFT = S.inverseRealFFT(lpf)
    inputSignal = IFFT
    S.plotSpectrogram(inputSignal)


def testPlotHPF(start=30, stop=31, cutoff=5000):
    S = SongModifier('./samples/nohands.wav')
    FFT = S.partialRealFFT(start, stop)
    lpf = S.highPassFFT(FFT, cutoff)
    IFFT = S.inverseRealFFT(lpf)
    inputSignal = IFFT
    S.plotFFT(inputSignal)


def main():
    # testPlotLPF()
    test_partial_spectrum()

if __name__ == "__main__":
    main()
