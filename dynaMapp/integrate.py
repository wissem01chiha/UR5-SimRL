import numpy as np

def integrate_time_vector(vector, sampling_rate):
    """
    Integrates a sampled vector with a given sampling rate to compute the time integral.

    Args:
        - vector         -  Input vector.
        - sampling_rate  -  Sampling frequency with which the vector was recorded.

    Returns:
        - result_vector (array): Integrated vector.
    """
    assert isinstance(sampling_rate, (int, float)) and sampling_rate > 0, \
        "Input validation failed: 'sampling_rate' must be a positive scalar."
    if sampling_rate > np.finfo(float).max:
        print("Warning: Sampling rate is too high; results may be inaccurate!")
        
    result_vector = np.zeros_like(vector)
    time_step = 1 / sampling_rate
    for i in range(len(vector)):
        if i == 0:
            result_vector[i] = vector[i]
        else:
            result_vector[i] = result_vector[i - 1] + 2 * vector[i]
    result_vector = (time_step / 2) * result_vector
    return result_vector
