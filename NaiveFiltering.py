from Analysis import *
import numpy


def main(start=60, stop=70):
    highCutoff = 2048
    lowCutoff = 512
    vocalBoostLowCutoff = 20000
    vocalBoostHighCutoff = 1000
    boostFactor = 2

    yesterday = SongModifier('./samples/yesterday.wav')
    FFT = yesterday.partialRealFFT(start, stop)
    bandPassedSignal = yesterday.bandStopFFT(FFT, lowCutoff, highCutoff)
    #vocalBoostSignal = yesterday.bandBoostFFT(bandPassedSignal, vocalBoostLowCutoff, vocalBoostHighCutoff, boostFactor)

    IFFT = yesterday.inverseRealFFT(bandPassedSignal)
    yesterday.writeWAV(IFFT, './processed/naiveFilter2.wav')
    yesterday.plotPartialSpectrogram(10, 25, yesterday.data, tempo=120)
    yesterday.plotSpectrogram(IFFT, tempo=120)

if __name__ == "__main__":
    main(10, 20)