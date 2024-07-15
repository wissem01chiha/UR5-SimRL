import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from scipy.stats import pearsonr

def joint_velocity(data_file_path, ndof, fc=0.2, fs=100, filter_order=2):
    """
    jointVelocity

    Inputs:
        data_file_path - Robot state variables measurements data file
        ndof           - Robot degree of freedom
        fc             - cutoff frequency
        fs             - sampling frequency
        filter_order   - low pass Butterworth filter order

    Returns: None
    """
    # Read the data file
    data = pd.read_csv(data_file_path)
    recorded_velocity = data.iloc[:, ndof+1:2*ndof+1].to_numpy()
    desired_velocity = data.iloc[:, 8*ndof+1:9*ndof+1].to_numpy()

    # Design the Butterworth filter
    b, a = butter(filter_order, fc / (fs / 2))

    # Apply the filter
    filter_recorded_velocity = filtfilt(b, a, recorded_velocity, axis=0)

    # Calculate the relative error
    relative_error = np.abs(filter_recorded_velocity - desired_velocity) / np.abs(desired_velocity) * 100
    relative_error[np.isinf(relative_error)] = 0

    # Calculate the correlation matrix
    velocity_correlation_matrix = np.corrcoef(filter_recorded_velocity, rowvar=False)

    # Create output directory if it doesn't exist
    output_dir = "figures/jointVelocity"
    os.makedirs(output_dir, exist_ok=True)

    # Plot the velocities
    plt.figure(1, figsize=(10, 6))
    for i in range(ndof):
        plt.subplot(3, 3, i+1)
        plt.plot(recorded_velocity[:, i], 'b-', label='recorded')
        plt.plot(desired_velocity[:, i], 'r--', label='desired')
        plt.plot(filter_recorded_velocity[:, i], 'g-.', label='filtered')
        plt.legend()
        plt.xlabel("Time (seconds)")
        plt.ylabel("Angular Velocity")
        plt.title(f'Joint {i+1}')
    plt.suptitle(f'Joints velocity, sampling frequency = {fs}, cutoff frequency = {fc}', fontsize=11)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_dir, 'joints_velocities.png'))

    # Plot the relative error
    plt.figure(2, figsize=(10, 6))
    for i in range(ndof):
        plt.subplot(3, 3, i+1)
        plt.plot(relative_error[:, i])
        plt.xlabel("Time (seconds)")
        plt.ylabel("Error (%)")
        plt.title(f'Joint {i+1}')
    plt.suptitle(f'Joints velocity relative error, filtered measurements, sampling = {fs}, cutoff frequency = {fc}', fontsize=11)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_dir, 'joints_velocity_relative_error.png'))

    # Print the correlation matrix
    print('Joints velocity correlation matrix =')
    print(velocity_correlation_matrix)

# Example usage
joint_velocity('data/RunTrajectory1.csv', 6)
