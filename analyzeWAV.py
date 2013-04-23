import scikits.audiolab as audiolab


def analyzeWAV(inputFile):
    """
    inputFile = .wav audiofile
    returns array of audiodata and the sampling rate
    """
    data, fs, nbits = audiolab.wavread(inputFile)
    samplingRate = fs
    return [data, samplingRate]
