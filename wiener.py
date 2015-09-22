import numpy as np


def edgeSmoothing(amplitude, smoothPoints=None, smoothFraction=0.05):
    N = len(amplitude)
    if smoothPoints is None:
        n = int(N * smoothFraction)
    else:
        n = smoothPoints
    amp = amplitude.copy()
    amp[0] = 0
    amp[1:n-1] *= (np.cos(np.arange(n-2)*np.pi/2/(n-2))**2)[::-1]
    amp[-n+1:-1] *= (np.cos(np.arange(n-2)*np.pi/2/(n-2))**2)
    amp[-1] = 0
    return amp


def noise(signal, SNR):
    return np.real(np.ifft(SNR / (SNR + 1) * np.fft(signal[:])))


def deconvolution(signal, SNR, response=None, response_fft=None):
    "Return a wiener filter deconvolution"
    if response_fft is not None:
        R = response_fft
    elif response is not None:
        R = np.fft.fft(response)
    else:
        R = 1

    return np.real(np.fft.ifft(SNR * np.conj(R) / (SNR * np.abs(R)**2 + 1.) *
                   np.fft.fft(signal[:])))


class Deconcolver:
    def __init__(self, snr, response):
        self._snr = snr
        self._response = response
        R = np.fft.fft(response)
        snr_R2 = snr * np.abs(R)**2

        self._G = (1. / R) * snr_R2 / (snr_R2 + 1)

    def deconvolve(self, signal):
        return np.real(np.fft.ifft(self._G * np.fft.fft(signal)))
