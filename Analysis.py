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

    def plotSpectrogram(self):
        Pxx, freqs, bins, im = pylab.specgram(self.data, Fs=self.samplingRate)
        pylab.show()

    def spectrogram(self):
        Pxx, freqs, bins, im = pylab.specgram(self.data, NFFT = 256, Fs=self.samplingRate, noverlap=128,)
        return [Pxx, freqs]

    def changeSong(self, songPath):
        self.__init__(songPath)


def main():
    S = SongAnalyzer('./samples/nohands.wav')
    S.plotSpectrogram()

if __name__ == "__main__":
    main()
