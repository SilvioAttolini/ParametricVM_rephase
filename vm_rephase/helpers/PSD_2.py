import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, csd, coherence

# Define the signals and sampling frequency
# Example: Two sample sine waves with noise
fs = 2000  # Sampling frequency in Hz
t = np.linspace(0, 1, fs, endpoint=False)  # Time array for 1 second

# Define two signals with a common component and some noise
f_signal = 50  # Common signal frequency in Hz
x = np.sin(2 * np.pi * f_signal * t) # + 0.5 * np.random.normal(size=t.shape)  # Signal 1
y = np.sin(2 * np.pi * f_signal * t + np.pi/4) # + 0.5 * np.random.normal(size=t.shape)  # Signal 2 with a phase shift

# Parameters for the cross-PSD calculation
window = 'hamming'  # Window function (e.g., Hamming window)
nperseg = 1024  # Length of each segment
noverlap = 512  # Number of overlapping samples
nfft = 1024  # Number of FFT points

# Calculate the PSD using the Welch method
f, pxx = welch(x, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)
_, pyy = welch(y, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)

# Plot the PSD
plt.figure()
plt.semilogy(f, pxx)  # Use semilogarithmic scale for better visualization
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power/Frequency (dB/Hz)')
plt.title('Power Spectral Density')
plt.grid(True)
plt.show()

# Calculate the cross-PSD using the Welch method
f, Pxy = csd(x, y, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)

# Plot the cross-PSD
plt.figure()
plt.semilogy(f, np.abs(Pxy))  # Use semilogarithmic scale for better visualization
plt.xlabel('Frequency (Hz)')
plt.ylabel('Cross Power/Frequency (dB/Hz)')
plt.title('Cross Power Spectral Density')
plt.grid(True)
plt.show()

# Calculate the coherence
Cxy = coherence(x, y, fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft)[1]

# Plot the coherence
plt.figure()
plt.plot(f, Cxy)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Coherence')
plt.title('Coherence between x and y')
plt.grid(True)
plt.show()
