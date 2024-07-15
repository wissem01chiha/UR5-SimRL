"""
This module contain electrical models implmentaion of actuator systems
"""
import numpy as np
from scipy.signal import square
from scipy.signal import lti, lsim

class BLDC(object):
    
    def __init__(self) -> None:
        """
        Brushless Direct Current Motor Model Class.

        Args:
            - Jm (float)        - Robot inertia factor.
            - L (float)         - Armature inductance.
            - R (float)         - Armature resistance.
            - Ke (float)        - Constant speed.
            - D (float)         - Motor damping coefficient.
            - Kt (float)        - Motor current coefficient.
            - Tck (float)       - Motor cogging torque coefficients.
            - Ta (float)        - Motor mechanical disturbance coefficient A.
            - Tb (float)        - Motor mechanical disturbance coefficient B.
            - sampling (float)  - Sampling frequency of the PWM signal.
            - PWM_input (float) - PWM ESC controller input signal.

        Returns:
            - t (array)     - Simulation time vector.
            - Q (array)     - Rotor angular position vector [-pi, pi].
            - dQ_dt (array) - Rotor angular velocity vector.
            - Ia (array)    - Armature current vector.
            - Va (array)    - Motor phase voltage vector.
            - T (array)     - Motor total torque vector.
            - Tf (array)    - Motor viscous friction torque vector.
            - Tff (array)   - Motor friction torque with respect to Lugre model.
            - Tcogg (array) - Motor cogging torque vector.

        Ref:
            - Practical Modeling and Comprehensive System Identification of a BLDC 
            Motor - C.iang, X.Wang, Y.Ma, B.Xu - 2015.
        """
        # memebre variables 
        self.Jm = None
        self.L = None
        self.R = None
        self.kv =None
        self.ke= None
        
        self.PWM_input = None
        self.voltage = None 
        self.sampling = 500
        
        self.controller_gain = None 
        self.controller_time = None 
        
        self.time = None 
        self.tstep = None
        self.tspan = None
        
    @classmethod
    def _init__state(self, sample_number)-> None:
        self.tspan = (sample_number- 1) / self.sampling
        self.tstep = 1 / self.sampling
        self.time = np.arange(0, self.tspan + self.tstep, self.tstep)
    
    @classmethod
    def set_random_input(self)-> None:
        self.PWM_input = 1
         
    @classmethod
    def _init__controller(self, gain, time_factor) -> None:
        self.controller_gain = gain
        self.controller_time = time_factor
        
    @classmethod
    def lsim_controller(self, gain = 0.578, time_factor = 0.7)-> None:
        
        self._init__controller(gain, time_factor)
        G = lti([self.controller_gain], [self.controller_time, 1])
        if self.PWM_input == None:
            print("error controller input signal not intilized")
        else:
            self._init__state(len(self.PWM_input))
            self.time, self.voltage, _ = lsim(G, self.PWM_input, self.time)

    @classmethod
    def lsim(self):
        """
        
        """
        return
    
    @classmethod
    def nlsim(self):
        """
        Simulate the non linear model 
        """
        return 






 




def PWM(frequency, sampling_rate, duty_cycle, magnitude, duration, mode='multi_channel'):
    """
    Generates a Pulse Width Modulated (PWM) Signal.

    Parameters:
    frequency (list)     - Signal frequency in Hz.
    sampling_rate (list) - Sampling frequency in Hz.
    duty_cycle (list)    - Duty cycle value of the signal (0 to 1).
    magnitude (list)     - Signal max value.
    duration (list)      - Total signal time duration in seconds.
    mode (str)           - 'single_channel' or 'multi_channel' (default).

    Returns:
    signal (ndarray)     - PWM sampled signal for each channel.
    """
    channel_number = len(frequency)
    if not (channel_number == len(sampling_rate) == len(duty_cycle) == len(magnitude) == len(duration)):
        raise ValueError("Invalid number of channels parameters inputs!")
    
    for dc in duty_cycle:
        if dc > 1 or dc < 0:
            raise ValueError("Duty cycle should be within [0,1]!")
    for mag in magnitude:
        if mag <= 0:
            raise ValueError("Signal magnitude should be a positive value")

    if mode == 'multi_channel':
        t_min = None
        for i in range(channel_number):
            t_i = np.arange(0, duration[i], 1 / (sampling_rate[i] * frequency[i]))
            if t_min is None or len(t_i) < len(t_min):
                t_min = t_i
        
        signal = np.zeros((channel_number, len(t_min)))
        for i in range(channel_number):
            t_i = np.arange(0, duration[i], 1 / (sampling_rate[i] * frequency[i]))
            step_fn = magnitude[i] * square(2 * np.pi * frequency[i] * t_min, duty_cycle[i] * 100)
            signal[i, :len(t_i)] = np.maximum(step_fn[:len(t_i)], np.zeros_like(step_fn[:len(t_i)]))
    
    elif mode == 'single_channel':
        max_sampling_rate = max(sampling_rate)
        signal = np.array([])
        for i in range(channel_number):
            t_i = np.arange(0, duration[i], 1 / (max_sampling_rate * frequency[i]))
            signal_i = magnitude[i] * square(2 * np.pi * frequency[i] * t_i, duty_cycle[i] * 100)
            signal = np.concatenate((signal, signal_i))
        signal = np.maximum(signal, np.zeros_like(signal))
    
    else:
        raise ValueError("Unsupported PWM generation option!")
    
    return signal







