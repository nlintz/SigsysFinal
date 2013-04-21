import scipy
import wave
import struct
import numpy
import pylab

from scipy.io import wavfile

rate, data = wavfile.read('./a440.wav')
filtereddata = numpy.fft.rfft(data, axis=0)

print data

filteredwrite = numpy.fft.irfft(filtereddata, axis=0)

print filteredwrite

wavfile.write('TestFiltered.wav', rate, filteredwrite)
