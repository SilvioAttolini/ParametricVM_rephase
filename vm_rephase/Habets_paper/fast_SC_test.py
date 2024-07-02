from scipy.signal import stft
from helpers.auto_and_cross_spectrum import auto_and_cross_power_spectrum
from helpers.create_ground_sinc import create_ground_sinc
from helpers.plots import *
from scipy.io import wavfile


def fast_SC_test():
    """ """

    params = {'Fs': 16000,
              'winLength': 4096,
              'sigLenSecs': 5,
              'd': 0.2,  # 0.078520,
              }
    save_spatial_coherence = ""

    file_name = 'my_implementation/habets_result.wav'
    Fs, stereo_data = wavfile.read(file_name)
    data1 = stereo_data[:, 0]  # Left channel
    data2 = stereo_data[:, 1]  # Right channel

    fAx_stft, tAx_stft, stft_of_data1 = stft(data1, fs=Fs, window='hann', nperseg=256, noverlap=128)
    fAx_stft, tAx_stft, stft_of_data2 = stft(data2, fs=Fs, window='hann', nperseg=256, noverlap=128)

    # plot_signal_time(params, imported_rir_time, vm, save_signal_time, False)
    # plot_signal_stft(tAx_stft, fAx_stft, complete_vms_stft[vm, :, :], vm, save_signal_stft, False)

    print("Computing Spatial Coherence...")

    # create current signals
    A = stft_of_data1
    B = stft_of_data2

    # compute the complex spatial coherence
    phi_na_nb = auto_and_cross_power_spectrum(A, B)
    phi_n_n = auto_and_cross_power_spectrum(A, A)  # ASSUMPTION: IT IS EQUAL TO (B,B)
    Gamma_nab = phi_na_nb / phi_n_n  # phi_na_nb / np.sqrt(phi_n_n)

    # keep just the real part to plot
    real_Gamma_nab = np.real(Gamma_nab)

    sinc = create_ground_sinc(fAx_stft, params)
    plot_spatial_coherence(fAx_stft, real_Gamma_nab, sinc, 0, 1, save_spatial_coherence, True)

    return


# def main():
#     fast_SC_test()
#
#
# if __name__ == '__main__':
#     main()