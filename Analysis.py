from pylab import plot, show, xlabel, ylabel
from scipy import fft
import scipy.fftpack
import pylab
import analyzeWAV


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

    def plotFFT(self):
        """
        plots the fft of an array
        inputSignal = numpy array
        samplingRate = sampling frequency for FFT
        """
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
        n = len(inputSignal)  # length of the signal
        FFT = abs(fft(inputSignal))/n  # fft computing and normalization
        FFT = FFT[range(n/2)]  # we only want the positive signals, fft returns the positive and negative frequencies
        freqs = scipy.fftpack.fftfreq(inputSignal.size, 1.0/samplingRate)
        freqs = freqs[range(n/2)]  # only plot the positive frequencies
        return [freqs, abs(FFT)]

    def plotFFTDownSample(self, downsamplingPercent):
        """
        downsamplingPercent = percent of samples you want ie 10 = every 10th sample
        """
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

    def plotSpectrogram(self, tempo=120):
        NFFT = int(8.0 * self.samplingRate / float(tempo))
        noverlap = NFFT >> 2
        if (len(self.data.shape) == 2):
            Pxx, freqs, bins, im = pylab.specgram(self.data[:, 0], NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            pylab.show()
        else:
            Pxx, freqs, bins, im = pylab.specgram(self.data, NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            pylab.show()

    def spectrogram(self, tempo=120):
        NFFT = int(8.0 * self.samplingRate / float(tempo))
        noverlap = NFFT >> 2
        if (len(self.data.shape) == 2):
            Pxx, freqs, bins, im = pylab.specgram(self.data[:, 0], NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            return [Pxx, freqs, bins]
        else:
            Pxx, freqs, bins, im = pylab.specgram(self.data, NFFT=NFFT, Fs=self.samplingRate, noverlap=noverlap)
            return [Pxx, freqs, bins]

    def changeSong(self, songPath):
        self.__init__(songPath)
