import spect
import scikits.audiolab as audiolab


def analyzeWAV(inputFile):
    """
    inputFile = .wav audiofile
    returns array of audiodata and the sampling rate
    """
    data, fs, nbits = audiolab.wavread(inputFile)
    samplingRate = fs
    return [data, samplingRate]


def main():
    d = analyzeWAV('./samples/a440.wav')
    spect.plotSpectrum(d[0], d[1])


if __name__ == "__main__":
    main()
