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


def main():
    test_a440_fft()

if __name__ == "__main__":
    main()
