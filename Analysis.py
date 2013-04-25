from pylab import plot, show, xlabel, ylabel
from scipy import fft, ifft
import scipy.fftpack
import pylab
import analyzeWAV
from scikits.audiolab import Format, Sndfile
import numpy


class SongAnalyzer(object):
    """
    this object is used to perform ffts and spectrogram
    analysis on songs
    """
    def __init__(self, songPath):
        self.song = songPath
        self.analysis = analyzeWAV.analyzeWAV(self.song)
        self.data = self.analysis[0]
        self.samplingRate = self.analysis[1]

    def plotFFT(self, inputSignal=None):
        """
        plots the fft of an array
        inputSignal = numpy array
        samplingRate = sampling frequency for FFT
        """
        if inputSignal is None:
            inputSignal = self.data
        samplingRate = self.samplingRate
        n = len(inputSignal)  # length of the signal
        FFT = abs(fft(inputSignal))/n  # fft computing and normalization
        FFT = FFT[range(n/2)]  # we only want the positive signals, fft returns the positive and negative frequencies
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        freqs = freqs[range(n/2)]  # only plot the positive frequencies

        plot(freqs, abs(FFT), 'r')  # plotting the spectrum
        xlabel('Freq (Hz)')
        ylabel('|Y(freq)|')
        show()

    def FFT(self):
        inputSignal = self.data
        samplingRate = self.samplingRate
        FFT = numpy.fft.fft(inputSignal)  # fft computing and normalization
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        return [freqs, FFT]

    def partialFFT(self, timeStart, timeEnd):
        startIndex = timeStart * self.samplingRate
        stopIndex = timeEnd * self.samplingRate
        samplingRate = self.samplingRate
        inputSignal = self.data[startIndex:stopIndex]
        FFT = numpy.fft.fft(inputSignal, axis=0)  # fft computing and normalization
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        return [freqs, FFT]

    def realFFT(self):
        inputSignal = self.data
        samplingRate = self.samplingRate
        FFT = numpy.fft.rfft(inputSignal, axis=0)
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        return [freqs, FFT]

    def partialRealFFT(self, timeStart, timeEnd):
        startIndex = timeStart * self.samplingRate
        stopIndex = timeEnd * self.samplingRate
        samplingRate = self.samplingRate
        inputSignal = self.data[startIndex:stopIndex]
        amplitudes = numpy.fft.rfft(inputSignal, axis=0)  # fft computing and normalization
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        return [freqs, amplitudes]

    def plotPartialFFT(self, timeStart, timeEnd, inputSignal=None):
        """
        plots the fft of an array
        inputSignal = numpy array
        samplingRate = sampling frequency for FFT
        """
        startIndex = timeStart * self.samplingRate
        stopIndex = timeEnd * self.samplingRate
        samplingRate = self.samplingRate
        if inputSignal is None:
            inputSignal = self.data[startIndex:stopIndex]
        n = len(inputSignal)  # length of the signal
        FFT = abs(fft(inputSignal))/n  # fft computing and normalization
        FFT = FFT[range(n/2)]  # we only want the positive signals, fft returns the positive and negative frequencies
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        freqs = freqs[range(n/2)]  # only plot the positive frequencies

        plot(freqs, abs(FFT), 'r')  # plotting the spectrum
        xlabel('Freq (Hz)')
        ylabel('|Y(freq)|')
        show()

    def plotFFTDownSample(self, downsamplingPercent, inputSignal=None):
        """
        downsamplingPercent = percent of samples you want ie 10 = every 10th sample
        """
        if inputSignal is None:
            inputSignal = self.data
        samplingRate = self.samplingRate
        n = len(inputSignal)  # length of the signal
        FFT = abs(fft(inputSignal))/n  # fft computing and normalization
        FFT = FFT[range(n/2)]  # we only want the positive signals, fft returns the positive and negative frequencies
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        freqs = freqs[range(n/2)]  # only plot the positive frequencies

        plot(freqs[::downsamplingPercent], abs(FFT)[::downsamplingPercent], 'r')  # plotting the spectrum
        xlabel('Freq (Hz)')
        ylabel('|Y(freq)|')
        show()

    def inverseFFT(self, amplitudes):
        return numpy.fft.ifft(amplitudes, axis=0)

    def inverseRealFFT(self, amplitudes):
        return numpy.fft.irfft(amplitudes, axis=0)

    def plotSpectrogram(self, inputSignal=None, tempo=120):
        if inputSignal is None:
            inputSignal = self.data
        NFFT = int(8.0 * self.samplingRate / float(tempo))
        noverlap = NFFT >> 2
        if (len(inputSignal.shape) == 2):
            Pxx, freqs, bins, im = pylab.specgram(inputSignal[:, 0], NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            pylab.show()
        else:
            Pxx, freqs, bins, im = pylab.specgram(inputSignal, NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            pylab.show()

    def spectrogram(self, inputSignal=None, tempo=120):
        if inputSignal is None:
            inputSignal = self.data
        NFFT = int(8.0 * self.samplingRate / float(tempo))
        noverlap = NFFT >> 2
        if (len(inputSignal.shape) == 2):
            Pxx, freqs, bins, im = pylab.specgram(inputSignal[:, 0], NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            return [Pxx, freqs, bins]
        else:
            Pxx, freqs, bins, im = pylab.specgram(inputSignal, NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            return [Pxx, freqs, bins]

    def partialSpectrogram(self, timeStart, timeEnd, inputSignal=None, tempo=120):
        if inputSignal is None:
            inputSignal = self.data
        startIndex = timeStart * self.samplingRate
        stopIndex = timeEnd * self.samplingRate
        spect = self.spectrogram(inputSignal)
        Pxx = spect[0][startIndex:stopIndex]
        freqs = spect[0][startIndex:stopIndex]
        bins = spect[0][startIndex:stopIndex]
        return [Pxx, freqs, bins]

    def plotPartialSpectrogram(self, timeStart, timeEnd, inputSignal=None, tempo=120):
        if inputSignal is None:
            inputSignal = self.data
        startIndex = timeStart * self.samplingRate
        stopIndex = timeEnd * self.samplingRate
        NFFT = int(8.0 * self.samplingRate / float(tempo))
        noverlap = NFFT >> 2
        if (len(inputSignal.shape) == 2):
            Pxx, freqs, bins, im = pylab.specgram(inputSignal[startIndex:stopIndex, 0], NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            pylab.show()
        else:
            Pxx, freqs, bins, im = pylab.specgram(inputSignal[startIndex:stopIndex], NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            pylab.show()

    def changeSong(self, songPath):
        self.__init__(songPath)


class SongModifier(SongAnalyzer):
    def __init__(self, songPath):
        SongAnalyzer.__init__(self, songPath)

    def lowPassFFT(self, inputSignal, cutoff):  # input signal needs to be an array like [freqs, amplitudes]
        freqs = inputSignal[0]
        amplitudes = inputSignal[1]
        for i in range(len(amplitudes)):
            if freqs[i] > cutoff:
                amplitudes[i] = 0
        return amplitudes

    def highPassFFT(self, inputSignal, cutoff):  # input signal needs to be an array like [freqs, amplitudes]
        freqs = inputSignal[0]
        amplitudes = inputSignal[1]
        for i in range(len(amplitudes)):
            if freqs[i] < cutoff:
                amplitudes[i] = 0
        return amplitudes

    def writeWAV(self, inputSignal, filename):
        format = Format('wav')
        if (len(inputSignal.shape) == 2):
            f = Sndfile(filename, 'w', format, 2, self.samplingRate)
            f.write_frames(inputSignal)
            f.close()
        else:
            f = Sndfile(filename, 'w', format, 1, self.samplingRate)
            f.write_frames(inputSignal)
            f.close()

