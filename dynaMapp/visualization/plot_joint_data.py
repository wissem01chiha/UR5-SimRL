# VIEW RECORDED DATA
#
# Author: Wissem Chiha
# Last Revision: 26-04-2024
"""
this file contain
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.signal import butter, filtfilt

# Clear all (Not needed in Python)

# Create directory if it doesn't exist
output_dir = "figures/view_recorded_data"
os.makedirs(output_dir, exist_ok=True)

data_file_name = "data/RunTrajectory1.csv"
data = pd.read_csv(data_file_name)
time = data.iloc[:, 0].to_numpy()

# Plot recorded position data
plt.figure(1)
for i in range(1, 8):
    plt.plot(time, data.iloc[:, i].to_numpy())
plt.legend(["Joint_1", "Joint_2", "Joint_3", "Joint_4",
            "Joint_5", "Joint_6", "Joint_7"])
plt.title("Joints Positions")
plt.xlabel("Time (seconds)")
plt.ylabel("Angular Displacement (rad)")
plt.savefig(os.path.join(output_dir, "Joints_Positions.png"))
plt.clf()

# Plot recorded velocity data
plt.figure(2)
for i in range(8, 15):
    plt.plot(time, data.iloc[:, i].to_numpy())
plt.legend(["Joint_1", "Joint_2", "Joint_3", "Joint_4",
            "Joint_5", "Joint_6", "Joint_7"])
plt.title("Joints Velocities")
plt.xlabel("Time (seconds)")
plt.ylabel("Angular Velocity (rad/sec)")
plt.savefig(os.path.join(output_dir, "Joints_Velocities.png"))
plt.clf()

# Plot recorded Torques from the sensor
plt.figure(3)
for i in range(15, 22):
    plt.plot(time, data.iloc[:, i].to_numpy())
plt.legend(["Joint_1", "Joint_2", "Joint_3", "Joint_4",
            "Joint_5", "Joint_6", "Joint_7"])
plt.title("Joints Torques")
plt.xlabel("Time (seconds)")
plt.ylabel("Torque Value (Nm)")
plt.savefig(os.path.join(output_dir, "Joints_Torques.png"))
plt.clf()

# Plot the Torque simulation variables
plt.figure(4)
plt.subplot(3, 1, 1)
for i in range(22, 29):
    plt.plot(time, data.iloc[:, i].to_numpy())
plt.legend(["Joint_1", "Joint_2", "Joint_3", "Joint_4",
            "Joint_5", "Joint_6", "Joint_7"])
plt.title("Joints Torques_Sim")
plt.xlabel("Time (seconds)")
plt.ylabel("Torque_Sim Value (Nm)")

# Plot the FFT of the joint sim torques
plt.subplot(3, 1, 2)
Fs = 200  # Sampling frequency
for i in range(22, 29):
    torque = data.iloc[:, i].to_numpy()
    L = len(torque)  # Length of signal
    Y = fft(torque)
    P2 = np.abs(Y / L)
    P1 = P2[:L // 2 + 1]
    P1[1:-1] = 2 * P1[1:-1]
    f = Fs * np.arange((L // 2) + 1) / L
    plt.plot(f, P1)
plt.title("Joints FFT simTorques")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.show()
# Apply a low pass Butterworth filter to joint sim_torques
fc = 0.3  # Cutoff frequency
fs = 100  # Sampling rate
plt.subplot(3, 1, 3)
for i in range(22, 29):
    b, a = butter(6, fc / (fs / 2))
    filtered_torques = filtfilt(b, a, data.iloc[:, i].to_numpy())
    plt.plot(time, filtered_torques)
plt.title("Joints filtered Torque_Sim")
plt.xlabel("Time (seconds)")
plt.ylabel("Filtered joints Torque_Sim")
plt.savefig(os.path.join(output_dir, "Joints_Torque_Sim.png"))
plt.show()
plt.clf()

# Plot joint torque_cur:
# couple calculated from the current sent to the motors
# @todo: how to get it?
plt.figure(6)
plt.subplot(3, 1, 1)
for i in range(29, 36):
    plt.plot(time, data.iloc[:, i].to_numpy())
plt.legend(["Joint_1", "Joint_2", "Joint_3", "Joint_4",
            "Joint_5", "Joint_6", "Joint_7"])
plt.title("Joints Computed Torque")
plt.xlabel("Time (seconds)")
plt.ylabel("Torque value")
plt.show()
Fs = 50  # Sampling frequency
b, a = butter(6, fc / (fs / 2))  # Butterworth filter of order 6

# Compute and plot FFT transformations of torque_Sim
plt.subplot(3, 1, 2)
for i in range(29, 36):
    filtered_torques = filtfilt(b, a, data.iloc[:, i].to_numpy())
    plt.plot(time, filtered_torques)
plt.title("Filtered Joints Computed Torque")
plt.xlabel("Time (seconds)")
plt.ylabel("Filtered Torque value")

plt.subplot(3, 1, 3)
for i in range(29, 36):
    filtered_torques_fft = filtfilt(b, a, data.iloc[:, i].to_numpy())
    L = len(filtered_torques_fft)
    Y = fft(filtered_torques_fft)
    P2 = np.abs(Y / L)
    P1 = P2[:L // 2 + 1]
    P1[1:-1] = 2 * P1[1:-1]
    f = Fs * np.arange((L // 2) + 1) / L
    plt.plot(f, P1)
plt.title("FFT of Joints Computed Torque")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.show()
plt.savefig(os.path.join(output_dir, "Joint_Torque_and_FFT.png"))
plt.clf()

# Plot the joint Temperature
plt.figure(7)
for i in range(36, 43):
    plt.plot(time, data.iloc[:, i].to_numpy())
plt.legend(["Joint_1", "Joint_2", "Joint_3", "Joint_4",
            "Joint_5", "Joint_6", "Joint_7"])
plt.title("Joints Temperature")
plt.xlabel("Time (seconds)")
plt.ylabel("Temperature")
plt.savefig(os.path.join(output_dir, "Joints_Temperature.png"))
plt.clf()

# Phase plan analysis: plot the joint velocity with respect to
# the joint positions
plt.figure(8)
for i in range(7):
    plt.plot(data.iloc[:, i + 1].to_numpy(), data.iloc[:, i + 8].to_numpy(), '*')
plt.legend(["Joint_1", "Joint_2", "Joint_3", "Joint_4",
            "Joint_5", "Joint_6", "Joint_7"])
plt.title("Joints Phase Plan trajectory")
plt.xlabel("Joint Position")
plt.ylabel("Joint velocity")
plt.savefig(os.path.join(output_dir, "Joints_PhasePlane.png"))
plt.clf()

# Bounded region within which the system's state remains confined over time.
# Closed-loop trajectories often suggest stability in the system.

# Estimate data parameters (to be implemented)
