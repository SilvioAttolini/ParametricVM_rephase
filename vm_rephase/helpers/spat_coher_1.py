from scipy.io.wavfile import write
from helpers.auto_and_cross_spectrum import auto_and_cross_power_spectrum
from helpers.plots import *


def spat_coher_1(complete_vms_time, complete_vms_stft, n_Vms, fAx_stft, sinc,
                 save_spatial_coherence, save_stereo, params):
    micA = 0
    micB = micA + 1
    i = 0
    real_spat_coh_of_all_pairs = np.zeros((n_Vms // 2, len(fAx_stft)))
    while micB < n_Vms:
        print(f"pair {micA + 1}-{micB + 1}")

        # create current signals
        A = complete_vms_stft[micA, :, :]
        B = complete_vms_stft[micB, :, :]

        # compute the complex spatial coherence
        phi_na_nb = auto_and_cross_power_spectrum(A, B)
        phi_n_n = auto_and_cross_power_spectrum(A, A)  # ASSUMPTION: IT IS EQUAL TO (B,B)
        Gamma_nab = phi_na_nb / phi_n_n  # phi_na_nb / np.sqrt(phi_n_n)

        # keep just the real part to plot
        real_Gamma_nab = np.real(Gamma_nab)

        plot_spatial_coherence(fAx_stft, real_Gamma_nab, sinc, micA, micB, save_spatial_coherence, False)

        real_spat_coh_of_all_pairs[i, :] = real_Gamma_nab

        # stereo audio of this pair
        signal_left = np.int16(complete_vms_time[micA, :] / np.max(np.abs(complete_vms_time[micA, :])) * 32767)
        signal_right = np.int16(complete_vms_time[micB, :] / np.max(np.abs(complete_vms_time[micB, :])) * 32767)
        stereo_signal = np.column_stack((signal_left, signal_right))
        write(f"{save_stereo}/stereo_signal_vms{micA + 1}_{micB + 1}", params['Fs'], stereo_signal)

        # loop
        micA = micB + 1
        micB = micA + 1
        i = i + 1

    # mean of the spatial coherence over all the pairs
    spatial_coherence_avg = np.mean(real_spat_coh_of_all_pairs, axis=0)

    # plot average spatial coherence
    plot_average_spatial_coherence(fAx_stft, spatial_coherence_avg, sinc, n_Vms, save_spatial_coherence, False)

    return
