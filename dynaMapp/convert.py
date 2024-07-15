import numpy as np

def convert_to_rad(deg_angle):
    """ Convert an angle from radians to degrees """
    return deg_angle * np.pi/180

def convert_to_deg(rad_angle):
    """ Convert an angle from radians to degrees """
    return rad_angle * 180 /np.pi

def convert_to_rad_s(rpm_velocity):
    """ Convert velocity from RPM to rad/s."""
    return rpm_velocity * (2 * np.pi / 60)

def convert_to_rpm(rad_s_velocity):
    """ Convert velocity from rad/s to RPM. """
    return rad_s_velocity * (60 / (2 * np.pi))

def convert_to_m_s(mm_s_value):
    """ Convert velcity from mm/s to m/s"""
    return mm_s_value