from Analysis import *
import numpy


def main(start=60, stop=90):
    highCutoff = 600
    lowCutoff = 5000
    vocalBoostLowCutoff = 24000
    vocalBoostHighCutoff = 5000
    boostFactor = 5

    yesterday = SongModifier('./samples/yesterday.wav')
    FFT = yesterday.partialRealFFT(start, stop)
    bandPassedSignal = yesterday.bandStopFFT(FFT, lowCutoff, highCutoff)
    vocalBoostSignal = yesterday.bandBoostFFT(bandPassedSignal, vocalBoostLowCutoff, vocalBoostHighCutoff, boostFactor)
    IFFT = yesterday.inverseRealFFT(vocalBoostSignal)
    yesterday.writeWAV(FFT[1], './processed/naiveFilter.wav')
    yesterday.plotPartialSpectrogram(0, 5, IFFT, tempo=120)

if __name__ == "__main__":
    main(10, 15)
