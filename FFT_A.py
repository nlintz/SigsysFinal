import scipy
import wave
import struct
import numpy
import pylab

from scipy.io import wavfile


def fftTest(inputFile, outputFile):
    rate, data = wavfile.read(inputFile)
    print 'wav read'
    filtereddata = numpy.fft.rfft(data, axis=0)
    print 'fft done'
    print data
    filteredwrite = numpy.fft.irfft(filtereddata, axis=0)
    print 'filteredwrite'
    print filteredwrite
    wavfile.write(outputFile, rate, filteredwrite)


def main():
    fftTest('./samples/nohands.wav', './processed/nohandsTest.wav')


if __name__ == "__main__":
    main()
