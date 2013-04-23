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
    S = SongModifier('./samples/a440.wav')
    S.writeWAV(S.data, './processed/a440.wav')


def main():
    test_wavwrite()

if __name__ == "__main__":
    main()
