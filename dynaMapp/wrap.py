import numpy as np

def wrap_to_interval(angle, min_value, max_value):
    """
    Wrap an angle to the interval [min_value, max_value].
    """
    range_size = max_value - min_value
    wrapped_angle = (angle - min_value) % range_size + min_value
    return wrapped_angle

def wrap_to_value(angle,value):
    """
    Wrap angle to the interval [-value, value].
    """
    wrapped_angle = angle - 2 * abs(value) * np.floor((angle + abs(value)) / (2 * abs(value)))
    return wrapped_angle

def wrap_to_pi(angle):
    """
    Wrap angle to the interval [-pi, pi]
    """
    wrapped_angle = angle - 2 * np.pi * np.floor((angle + np.pi) / (2 * np.pi))
    return wrapped_angle

def wrap_to_deg(angle):
    """
    Wrap angle to the interval [-180, 180]
    """
    wrapped_angle = angle - 360 * np.floor((angle + 180) / 360)
    return wrapped_angle
