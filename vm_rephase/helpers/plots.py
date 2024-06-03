import matplotlib.pyplot as plt
import numpy as np


def plot_complex(x_ax, real_part, imag_part, method, t, vm, save_directory, show):
    name_and_title = f"Phase @t={t} of vm {vm + 1}, method {method}"

    plt.figure(figsize=(10, 5))
    plt.subplot(2, 1, 1)  # 2 rows, 1 column, 1st subplot
    plt.plot(x_ax, real_part, label='Real Part')
    plt.title(f"{name_and_title}, Real part")
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Real value')
    plt.grid(True)
    plt.subplot(2, 1, 2)  # 2 rows, 1 column, 2nd subplot
    plt.plot(x_ax, imag_part, label='Imaginary Part', color='red')
    plt.title('Imaginary Part')
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Imaginary value')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{save_directory}/{name_and_title}.png', dpi=300)
    if show:
        plt.show()
    plt.close()


def plot_signal_time(params, imported_rir_time, vm, save_directory, show):
    # Create the time axis
    Ts = 1 / params['Fs']
    time_axis_sec = np.arange(0, len(imported_rir_time)) * Ts

    # index = np.where(time_axis_sec == 0.5)[0][0]
    # print(index)
    name_and_title = f"Complete signal of Vm RIR {vm + 1}"

    plt.plot(time_axis_sec, imported_rir_time)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    # plt.title(f'Diffuse contribution of Vm RIR {vm+1}')
    plt.title(f'Complete signal of Vm RIR {vm + 1}')
    plt.axvline(x=0.4, color='red', linestyle='--', label='x = 0.4')
    plt.axvline(x=0.5, color='red', linestyle='--', label='x = 0.5')
    plt.grid(True)
    plt.savefig(f'{save_directory}/{name_and_title}.png', dpi=300)
    if show:
        plt.show()
    plt.close()


def plot_signal_stft(t_ax, f_ax, imported_diff_contr_vm_rir_stft, vm, save_directory, show):
    name_and_title = f'STFT of Complete signal of Vm RIR {vm + 1}'

    plt.figure(figsize=(10, 6))
    plt.pcolormesh(t_ax, f_ax, np.log10(np.abs(imported_diff_contr_vm_rir_stft)), shading='gouraud')
    plt.colorbar(label='Intensity (dB)')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')
    plt.title(f'STFT of Complete signal of Vm RIR {vm + 1}')
    plt.tight_layout()
    plt.savefig(f'{save_directory}/{name_and_title}.png', dpi=300)
    if show:
        plt.show()
    plt.close()


def plot_sinc(f_ax, sinc, show):
    plt.plot(f_ax, sinc)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title(f'Sinc')
    plt.grid(True)
    if show:
        plt.show()
    plt.close()


def plot_spatial_coherence(freqs_gamma, real_Gamma_nab, sinc, micA, micB, save_directory, show):
    name_and_title = f'SpatialCoherence_{micA + 1}_{micB + 1}'

    plt.plot(freqs_gamma, real_Gamma_nab, label='SPA', color='blue')
    plt.plot(freqs_gamma, sinc, label='sinc', color='red')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Spatial Coherence')
    plt.title(f'Spatial Coherence of diffuse between mic {micA + 1} and mic {micB + 1}')
    plt.grid(True)
    plt.savefig(f'{save_directory}/{name_and_title}.png', dpi=300)
    if show:
        plt.show()
    plt.close()


def plot_average_spatial_coherence(freqs_gamma, spatial_coherence_avg, sinc, n_Vms, save_directory, show):
    name_and_title = f'SpatialCoherence_average'

    plt.plot(freqs_gamma, spatial_coherence_avg, label='SPA', color='blue')
    plt.plot(freqs_gamma, sinc, label='sinc', color='red')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Spatial Coherence')
    plt.title(f'Spatial Coherence averaged over the {n_Vms // 2} VM pairs, vs sinc')
    plt.grid(True)
    plt.savefig(f'{save_directory}/{name_and_title}.png', dpi=300)
    if show:
        plt.show()
    plt.close()
