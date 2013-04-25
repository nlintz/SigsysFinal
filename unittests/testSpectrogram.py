import sys
import os
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


def test_partial_spectrum(start=30, stop=35):
    S = SongAnalyzer('./samples/nohands.wav')
    S.plotPartialSpectrogram(start, stop)


def test_partial_fft():
    start = 30
    stop = 30.05
    S = SongAnalyzer('./samples/nohands.wav')
    S.plotPartialFFT(start, stop)


def test_wavwrite():
    S = SongModifier('./samples/nohands.wav')
    S.writeWAV(S.data, './processed/nohands.wav')


def test_writeInverseFFT():
    S = SongModifier('./samples/nohands.wav')
    FFT = S.partialRealFFT(0, 10)[1]
    IFFT = S.inverseRealFFT(FFT)
    S.writeWAV(IFFT, './processed/nohandsTEST.wav')


def main():
    test_writeInverseFFT()

if __name__ == "__main__":
    main()
