"""

© 2025 by Bode Plotter

MIT License

Copyright (c) [2025] [Bode Plotter]
Original code written by Simone Albano on 10/15/2023
Bode plot implementation for the handheld oscilloscope OWON HDS320S
Highly rewritten and modified between 01/20/2025 to 05/09/2025
Modifications were done by Bode Plotter with assistance from AI sourced from Google and Bing.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import usb.core
import usb.util
import numpy as np
from numpy.fft import fft, fftshift, fftfreq
import cmath
from scipy.interpolate import UnivariateSpline
import scipy.optimize as optimize
import dearpygui.dearpygui as dpg
import datetime
import time
import os
from pathlib import Path
import re  # Used for extracting the first number (integer or float) found in a string.
import multiprocessing
import queue
import threading
import json
import traceback
from collections import deque
from screeninfo import get_monitors
import platform
import warnings
import socket
import logging
import signal
import sys
import argparse

# 2025/01/30 
# /usr/lib/python3.12/contextlib.py:137: DeprecationWarning: anti_aliased keyword removed
#  return next(self.gen)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# 2025/03/19 Import the matplotlib class for bode-plotter
from bodeplots.managers import PlotManagerFFToscilloscope, PlotManagerMagnitudePhase, XYoscilloscope, PlotManagerSmithChart

from uuid import uuid4
import shutil  # Import for cleanup
import atexit  # Import for safe cleanup upon exit

# --- System variables --- ------------------------------------------------------------
available_channels = ['CH1', 'CH2']
probes_attenuation_ratio = ['1X', '10X']
channels_coupling_mode = ['AC', 'DC']  # AC or DC coupling (for Bode plots use AC coupling!)
sample_modes = ['SAMPle', 'PEAK']  # SAMPle mode is preferred
memory_depth_modes = ['4K', '8K']
AWG_output_impedance_modes = ['ON', 'OFF']
plot_win_settings = ['DearPyGui', 'MatPlotLib']
drag_line_values = ['Drag Lines ON', 'Drag Lines OFF']
points_spacing = ['Linear', 'Logarithmic']

# Global variables to control chart visibility and names:
# Index mapping: 0 -> XY, 1 -> MP, 2 -> FFT, 3 -> SC
chart_visibility = [True, True, True, False]
chart_names = ["XY Plot", "Bode Plots", "FFT/Oscilloscope", "Smith Chart"]

# Global system variables
global channel_in, channel_out, CH1_probe_attenuation_ratio, CH2_probe_attenuation_ratio
global CH1_coupling, CH2_coupling, Sample_command, DEPMEM, waveform_amplitude_V
global AWG_output_impedance, points_per_decade, start_decade, stop_decade, addpm
global point_resize_factor, vertical_scaling_factor, horizontal_scaling_factor
global nWaveSamples, FTTcorection, read_delay_ms, sample_delay_s, plot_win_disposition, LogFile
global JSONLogFile, plot_process, oscilloscope_OUT, oscilloscope_IN, oscilloscope
global is_playing, is_paused, play_speed, is_recording

global magnitude_phase_stop_event, fft_oscilloscope_stop_event, xy_stop_event, sc_stop_event
global magnitude_phase_ready_event, fft_oscilloscope_ready_event, xy_ready_event, sc_ready_event
global thread_magnitude_phase, thread_fft_oscilloscope, thread_xy_oscilloscope, thread_sc_oscilloscope
global CH1_amplitude_scales, CH2_amplitude_scales, MaxFreqAdjust
# Global thread objects
thread_magnitude_phase = None
thread_fft_oscilloscope = None
thread_xy_oscilloscope = None
thread_sc_oscilloscope = None

# Stop events for controlling thread termination
magnitude_phase_stop_event = threading.Event()
fft_oscilloscope_stop_event = threading.Event()
xy_stop_event = threading.Event()
sc_stop_event = threading.Event()

magnitude_phase_ready_event = threading.Event()
fft_oscilloscope_ready_event = threading.Event()
xy_ready_event = threading.Event()
sc_ready_event = threading.Event()

LogFile = "AUTO.csv"
JSONLogFile = "AUTO.json"
points_per_decade = 30
# Global variables to keep track of the signal state for play back of data
is_playing = False
is_paused = False
is_recording = False
play_speed = 1.0  # Default speed

time_bases_commands = [
    '2.0ns', '5.0ns', '10.0ns', '20.0ns', '50.0ns', '100ns', '200ns', '500ns', '1.0us', '2.0us', '5.0us', '10us',
    '20us', '50us', '100us', '200us', '500us', '1.0ms', '2.0ms', '5.0ms', '10ms', '20ms', '50ms', '100ms',
    '200ms', '500ms', '1.0s', '2.0s', '5.0s', '10s', '20s', '50s', '100s', '200s', '500s', '1000s'
]
time_bases_values = [
    0.000000002, 0.000000005, 0.00000001, 0.00000002, 0.00000005, 0.0000001, 0.0000002, 0.0000005, 0.000001,
    0.000002, 0.000005, 0.00001, 0.00002, 0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02,
    0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000
]
# Define a dictionary mapping each attenuation ratio to its corresponding
# amplitude scales (with both commands and values)
# 2025/04/24 Seeing 500uV with 10X and 1X this doesn't seem to be documented in HDS200_Series_SCPI_Protocol.pdf
#            Need to slow down SALe changes to have 200 ms between changes.
amplitude_scales = {
    '10X': {
        'commands': ['100V', '50.0V', '10.0V', '5.00V', '2.00V', '1.00V', '500mV', '200mV', '100mV'],
        'values': [100, 50, 10, 5, 2, 1, 0.5, 0.2, 0.1]
    },
    '1X': {
        'commands': ['10.0V', '5.00V', '2.00V', '1.00V', '500mV', '200mV', '100mV', '50.0mV', '20.0mV', '10.0mV'],
        'values': [10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    }
}
    
decades_list_string = [
    '100mHz', '1Hz', '10Hz', '100Hz', '1kHz', '10kHz', '100kHz', '1MHz', '10MHz'
]

decades_list = [
    1E-1, 1E0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6, 1E7
]

read_buffer_size = 100000  # buffer size in bytes (8192 or 8k)
horizontal_decades_n = 6  # 12 in total
vertical_decades_n = 4  # 8 in total

raw_frequencies_range = []
gain_linear = []
gain_Y = []
gain_X = []
phase_Y = []
phase_X = []
OSCx = []
OSCy = []
OSCxin = []
FFTx = []
FFTy = []
MxFFTx = []
MxFFTy = []
XYx = []
XYy = []
fft_magnitude = []
FFT_freqList = []
FFT_resultList = []

# --- end system variables --- --------------------------------------------------------


# --- Start Graphics Settings --- -----------------------------------------------------
win_vertical_border = 3
win_horizontal_border = 3
items_standard_width = 100
TitleBarHeight = 60

main_window_height = 1600
set_window_height = 1600

setY_pos = 100

# Adjust the viewport height minus the taskbar height (assuming taskbar height is 40 pixels)
# Get the display height
monitor = get_monitors()[0]
display_height = monitor.height

if set_window_height > display_height - 40:
    set_window_height = display_height - 40
    setY_pos=0

main_window_width = 1300
plot_window_height = main_window_height / 4 - win_horizontal_border
setting_window_height = main_window_height - plot_window_height
setting_window_width = 500
scrollbar_width = 20  # Adjust this value as needed for your scrollbar width
ScrollbarOffSet= -10
plot_window_width = main_window_width - setting_window_width
# Increase the viewport width to accommodate the scrollbar
if platform.system() == 'Windows':
    viewport_width = main_window_width + 2 * scrollbar_width
else:
    viewport_width = main_window_width + scrollbar_width
dpgWindow1 = main_window_width
dpgWindow2 = main_window_width 
dpgWindow3 = setting_window_width

# List of available themes
win_theme = ['Dark', 'Light']

# Global variables to store themes
dark_theme = None
light_theme = None

# --- End Graphics Settings --- -------------------------------------------------------
def signal_handler(sig, frame):
    logging.info("Gracefully shutting down...")
    is_playing = False
    is_paused = False
    is_recording = False

    # Signal all threads/processes to stop.
    shutdown_event.set()
    
    time.sleep(0.5)
    logging.info("Gracefully closing MatPlotLib windows now ...")
    
    # shutdown all running PlotManagerSmithChart instances.
    for osc in PlotManagerSmithChart.get_running_instances():
        try:
            osc.send_close_signal()
        except Exception as e:
            logging.error("Error stopping PlotManagerSmithChart instance: %s", e)

    # shutdown all running PlotManagerMagnitudePhase instances.
    for osc in PlotManagerMagnitudePhase.get_running_instances():
        try:
            osc.send_close_signal()
        except Exception as e:
            logging.error("Error stopping PlotManagerMagnitudePhase instance: %s", e)

    # shutdown all running plotting/oscilloscope instances.
    for osc in PlotManagerFFToscilloscope.get_running_instances():
        try:
            osc.send_close_signal()
        except Exception as e:
            logging.error("Error stopping PlotManagerFFToscilloscope instance: %s", e)
    
    # shutdown all running XYoscilloscope instances.
    for osc in XYoscilloscope.get_running_instances():
        try:
            osc.send_close_signal()
        except Exception as e:
            logging.error("Error stopping XYoscilloscope instance: %s", e)
            
    dpg.stop_dearpygui()  # Stop the Dear PyGui main loop
    # Shut down the logging system
    logging.shutdown()
    
# Global shutdown event for threads/processes
shutdown_event = threading.Event()

def print_to_terminal(data, color=None):
    """
    Add text data to the terminal with optional color and bold font, then auto-scroll.
    """
    dpg.add_text(data, parent='terminal', color=color)  # Default font

    dpg.set_y_scroll(item='terminal', value=-1)  # Auto-scroll to the end of the window
    dpg.focus_item(item='terminal')


def scroll_data_table():
    """
    Auto-scroll to the end of the data table and focus on it.
    """
    dpg.set_y_scroll(item='dataTableWindow', value=-1)  # Auto-scroll to the end of the window
    dpg.focus_item(item='dataTableWindow')

def oscilloscope_query(cmd):
    global sample_delay_s
    """
    Send a command to the oscilloscope and read the response.
    Retries until a valid response is received or an exception occurs.

    Args:
        cmd (str): The command to send to the oscilloscope.

    Returns:
        str: The response from the oscilloscope.
    """
    result = None
    while result is None:
        try:
            oscilloscope_OUT.write(cmd)
            result = oscilloscope_IN.read(read_buffer_size, read_delay_ms)
            if result is None:
                logging.info("Failed to get value from oscilloscope")
            else:
                return result.tobytes().decode('utf-8').strip()
        except Exception as e:
            logging.info("Oscilloscope command sent: " + str(cmd))
            logging.info(f"Exception in getting value from oscilloscope: {e}")
            oscilloscope.reset()
        time.sleep(sample_delay_s)

def oscilloscope_write(cmd):
    """
    Send a command to the oscilloscope without reading a response.

    Args:
        cmd (str): The command to send to the oscilloscope.
    """
    oscilloscope_OUT.write(cmd)

def oscilloscope_read(cmd):
    """
    Send a command to the oscilloscope and read the response.

    Args:
        cmd (str): The command to send to the oscilloscope.

    Returns:
        bytes: The raw response from the oscilloscope.
    """
    oscilloscope_OUT.write(cmd)
    result = oscilloscope_IN.read(read_buffer_size, read_delay_ms)
    return result


def sleep_with_abort(total_time, is_setup, check_interval=0.05):
    """
    Sleep for a total of total_time seconds, but with periodic checks.
    
    Args:
      total_time (float): Total time to sleep in seconds.
      is_setup (bool): Setup flag that affects the abort condition.
      check_interval (float): Time increments for each sleep check.
      
    Returns:
      bool: True if an abort condition occurred, otherwise False.
    """
    elapsed = 0.0
    while elapsed < total_time:
        # Abort if not recording and not in setup mode.
        # If is_setup is True, we are allowed to continue even if recording is off.
        if (not is_recording) and (not is_setup):
            return True  # Abort condition met.
        time.sleep(check_interval)
        elapsed += check_interval
    return False

def AWG_first_set_frequency(frequency, is_setup=False):
    """
    Set the Arbitrary Waveform Generator (AWG) frequency in Hz (0.1Hz to 25MHz for SINE waveform).

    Args:
      frequency (float): The desired frequency to set.
      is_setup (bool): Flag indicating if the function is in setup mode.

    Returns:
      tuple: (float measured_frequency, list PhaseList) after adjustment.
    """
    global sample_delay_s, horizontal_decades_n, is_recording, MaxFreqAdjust
    
    # Initialize variables
    Vpkpk = 10000
    ChkFREQ = 0.0
    FreqAdjust = 0
    PhaseNew = 0
    PhaseList = []
    
    # Adjust vertical scale
    # Retrieve the appropriate amplitude scale values for the selected attenuation
    if channel_out == 'CH2':
        amplitude_scales_values = amplitude_scales[CH2_probe_attenuation_ratio]['values']
        amplitude_scales_commands = amplitude_scales[CH2_probe_attenuation_ratio]['commands']
    if channel_out == 'CH1':
        amplitude_scales_values = amplitude_scales[CH1_probe_attenuation_ratio]['values']
        amplitude_scales_commands = amplitude_scales[CH1_probe_attenuation_ratio]['commands']
    # Measure frequency
    ChkFREQ = FREQ_out_to_float(oscilloscope_query(':MEASurement:{CH}:FREQuency?'.format(CH=channel_in)))
    CHin_amplitude_scales = get_amplitude_scale(CH1_probe_attenuation_ratio, channel_in, channel_in, amplitude_scales)
    CHout_amplitude_scales = get_amplitude_scale(CH2_probe_attenuation_ratio, channel_out, channel_in, amplitude_scales)
    if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
        return ChkFREQ
    set_amplitude_scale(channel_in, CHin_amplitude_scales)
    if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
        return ChkFREQ
    set_amplitude_scale(channel_out, CHout_amplitude_scales)
    if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
        return ChkFREQ

    # Loop to adjust frequency and minimize transient spikes
    while (int(ChkFREQ / 10 + 0.5) != int(frequency / 10 + 0.5)) and (MaxFreqAdjust > FreqAdjust):
        # Check if the process should be aborted
        if (not is_recording) and (not is_setup):
            break

        if (FreqAdjust) % 2 == 0:
            oscilloscope_OUT.write(':FUNCtion:FREQuency {}'.format(frequency))
        FreqAdjust += 1

        # Allow time for frequency change to settle using our abortable sleep
        if sleep_with_abort(FreqAdjust * sample_delay_s, is_setup):
            break

        Vpkpk = vertical_scale_to_float(get_pkpk_voltage(channel_out))
        closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (Vpkpk * vertical_scaling_factor))))
        set_amplitude_scale(channel_out, amplitude_scales_commands[closest_v_scale_index])
        if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
            break
        
        Vpkpk = vertical_scale_to_float(get_pkpk_voltage(channel_out))
        closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (Vpkpk * vertical_scaling_factor))))
        set_amplitude_scale(channel_out, amplitude_scales_commands[closest_v_scale_index])
        if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
            break        
    
        # Adjust horizontal scale
        closest_h_scale_index = time_bases_values.index(min(time_bases_values, key=lambda x: abs(x - ((1 / frequency) * horizontal_scaling_factor))))
        set_time_base(time_bases_commands[closest_h_scale_index])
        if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
            break
        
        # Measure frequency
        ChkFREQ = FREQ_out_to_float(oscilloscope_query(':MEASurement:{CH}:FREQuency?'.format(CH=channel_in)))
        #if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup):
        #    break
            
    # Return the final measured frequency and phase list
    return ChkFREQ


def AWG_set_frequency(frequency, is_setup=False):
    """
    Set the Arbitrary Waveform Generator (AWG) frequency in Hz (0.1Hz to 25MHz for SINE waveform).

    Args:
      frequency (float): The desired frequency to set.
      is_setup (bool): Flag indicating if the function is in setup mode.

    Returns:
      tuple: (float measured_frequency, list PhaseList) after adjustment.
    """
    global sample_delay_s, horizontal_decades_n, is_recording, MaxFreqAdjust
    
    # Initialize variables
    Vpkpk = 10000
    ChkFREQ = 0.0
    FreqAdjust = 0
    PhaseNew = 0
    PhaseList = []
    
    # Adjust vertical scale
    # Retrieve the appropriate amplitude scale values for the selected attenuation
    if channel_out == 'CH2':
        amplitude_scales_values = amplitude_scales[CH2_probe_attenuation_ratio]['values']
        amplitude_scales_commands = amplitude_scales[CH2_probe_attenuation_ratio]['commands']
    if channel_out == 'CH1':
        amplitude_scales_values = amplitude_scales[CH1_probe_attenuation_ratio]['values']
        amplitude_scales_commands = amplitude_scales[CH1_probe_attenuation_ratio]['commands']

    # Loop to adjust frequency and minimize transient spikes
    while (int(ChkFREQ / 10 + 0.5) != int(frequency / 10 + 0.5)) and (MaxFreqAdjust > FreqAdjust):
        # Check if the process should be aborted
        if (not is_recording) and (not is_setup):
            break

        if (FreqAdjust) % 2 == 0:
            oscilloscope_OUT.write(':FUNCtion:FREQuency {}'.format(frequency))
        FreqAdjust += 1

        # Allow time for frequency change to settle using our abortable sleep
        if sleep_with_abort(FreqAdjust * sample_delay_s, is_setup):
            break

        Vpkpk = vertical_scale_to_float(get_pkpk_voltage(channel_out))
        closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (Vpkpk * vertical_scaling_factor))))
        set_amplitude_scale(channel_out, amplitude_scales_commands[closest_v_scale_index])
        if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
            break     
        
        Vpkpk = vertical_scale_to_float(get_pkpk_voltage(channel_out))
        closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (Vpkpk * vertical_scaling_factor))))
        set_amplitude_scale(channel_out, amplitude_scales_commands[closest_v_scale_index])
        if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
            break        
    
        # Adjust horizontal scale
        closest_h_scale_index = time_bases_values.index(min(time_bases_values, key=lambda x: abs(x - ((1 / frequency) * horizontal_scaling_factor))))
        set_time_base(time_bases_commands[closest_h_scale_index])
        if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup): 
            break
        
        # Measure frequency
        ChkFREQ = FREQ_out_to_float(oscilloscope_query(':MEASurement:{CH}:FREQuency?'.format(CH=channel_in)))
        #if sleep_with_abort(sample_delay_s * FreqAdjust, is_setup):
        #    break

        # If conditions are met, perform the full sine analysis
        if int(ChkFREQ / 10 + 0.5) == int(frequency / 10 + 0.5):
            current_v_scaleout = vertical_scale_to_float(oscilloscope_query(':{}:SCALe?'.format(channel_out)))
            if sleep_with_abort(sample_delay_s, is_setup):
                break
            current_v_scalein = vertical_scale_to_float(oscilloscope_query(':{}:SCALe?'.format(channel_in)))
            if sleep_with_abort(sample_delay_s, is_setup):
                break
            raw_waveform_in = get_waveform(channel_in, current_v_scalein)
            if sleep_with_abort(sample_delay_s, is_setup):
                break
            raw_waveform_out = get_waveform(channel_out, current_v_scaleout)
            if sleep_with_abort(sample_delay_s, is_setup):
                break
            current_h_scale = horizontal_scale_to_float(oscilloscope_query(':HORIzontal:SCALe?'))
            # Compute full sine analysis (FSA)
            sampling_PERIOD = 1 / (current_h_scale / (300 / (horizontal_decades_n * 2)))
            FSA = full_sin_analysis(sampling_PERIOD, FTTcorection, raw_waveform_in, raw_waveform_out)
            PhaseNew = FSA['phase_difference_degrees']
            PhaseList.append(PhaseNew)
            
    # Return the final measured frequency and phase list
    return ChkFREQ, PhaseList

def set_time_base(period):
    """
    Set the oscilloscope time-base.
    
    Args:
        period (str): The desired time period to set on the oscilloscope.
    """
    oscilloscope_OUT.write(':HORizontal:SCALe {}'.format(period))

def set_amplitude_scale(channel, scale):
    channel = str(channel)
    scale = str(scale)
    """
    Set the oscilloscope amplitude scale.
    
    Args:
        channel (str): The channel to set the amplitude scale for.
        scale (str): The desired amplitude scale.
    """
    oscilloscope_OUT.write(':{CH}:SCALe {data}'.format(CH=channel, data=scale))

def get_pkpk_voltage(channel):
    """
    Get the peak-to-peak voltage for a given channel.
    
    Args:
        channel (str): The channel to measure the peak-to-peak voltage from.
        
    Returns:
        str: The peak-to-peak voltage.
    """
    return oscilloscope_query(':MEASurement:{CH}:PKPK?'.format(CH=channel)).lstrip('Vpp=')

def vertical_scale_to_float(voltage):
    """
    Convert the vertical scale voltage reading to a float (in volts).

    Args:
        voltage (str): The voltage reading as a string with units (e.g., "500mV", "500uV", "500µV", or "0.5V").

    Returns:
        float: The converted voltage in volts, scaled as appropriate.
    """
    # Extract the numeric part from the voltage string
    match = re.search(r"[-+]?\d*\.?\d+", voltage)
    if not match:
        return None
    value = float(match.group(0))
    
    # Check for units and convert accordingly
    if 'mV' in voltage:
        # millivolts: convert to volts (1mV = 1E-3 V)
        return value / 1E3
    elif 'uV' in voltage or 'µV' in voltage:
        # microvolts: convert to volts (1uV = 1E-6 V)
        return value / 1E6
    else:
        # If no unit is specified, assume the value is already in volts
        return value


def horizontal_scale_to_float(timescale):
    """
    Convert the horizontal scale time reading to a float.
    
    Args:
        timescale (str): The time scale reading as a string.
        
    Returns:
        float: The converted time scale as a float, scaled if necessary.
    """
    if 'ns' in timescale:
        timescale = timescale.rstrip('ns')
        return float(timescale) / 1E9
    elif 'us' in timescale:
        timescale = timescale.rstrip('us')
        return float(timescale) / 1E6
    elif 'ms' in timescale:
        timescale = timescale.rstrip('ms')
        return float(timescale) / 1E3
    else:
        timescale = timescale.rstrip('s')
        return float(timescale)

def FREQ_out_to_float(FREQ):
    """
    Convert the frequency reading to a float.
    
    Args:
        FREQ (str): The frequency reading as a string.
        
    Returns:
        float: The converted frequency as a float, scaled if necessary.
    """
    if 'OFF' in FREQ:
        return 0.0
    if 'Hz' in FREQ:
        FREQ = FREQ.lstrip('F=')
    else:
        return 0.0 
    if 'mHz' in FREQ:
        FREQ = FREQ.rstrip('mHz')
        return float(FREQ) * 1E-3
    elif 'kHz' in FREQ:
        FREQ = FREQ.rstrip('kHz')
        return float(FREQ) * 1E3
    elif 'MHz' in FREQ:
        FREQ = FREQ.rstrip('MHz')
        return float(FREQ) * 1E6
    elif 'Hz' in FREQ:
        FREQ = FREQ.rstrip('Hz')
        return float(FREQ)
    else:
        return 0.0

def get_waveform(channel, v_scale):
    """
    Retrieve waveform data from the oscilloscope.
    
    Args:
        channel (str): The channel to retrieve waveform data from.
        v_scale (float): The vertical scale for the waveform data.
        
    Returns:
        list: The waveform data points scaled appropriately.
    """

    # The first 4 bytes are discarded; in total, there are 600 points.
    rawdata = oscilloscope_read(':DATA:WAVE:SCREen:{}?'.format(channel))
    data_points = []
    
    # Process the raw data in 2-byte chunks
    for val in range(4, len(rawdata), 2):
        # Convert 2 bytes to signed integer using "little-endian"
        point = int.from_bytes([rawdata[val], rawdata[val + 1]], 'little', signed=True)
        # Scale the data point and add to the list
        data_points.append((point * v_scale) / point_resize_factor)
    
    return data_points

def play_callback(sender, app_data, user_data):
    global is_playing, is_paused
    if not is_playing:
        is_playing = True
        is_paused = False
        threading.Thread(target=PlayBack, daemon=True).start()
        dpg.configure_item('Play', label='Pause')
    elif is_paused:
        is_paused = False
        dpg.configure_item('Play', label='Pause')
    elif not is_paused:
        is_paused = True
        dpg.configure_item('Play', label='Resume')

def stop_callback(sender, app_data, user_data):
    global is_playing, is_paused, is_recording
    is_playing = False
    is_paused = False
    is_recording = False
    dpg.bind_item_theme('Stop', 'YellowButton')
    dpg.configure_item('Play', label='Play')

    
def speed_slider_callback(sender, app_data, user_data):
    global play_speed
    play_speed = app_data
    
def search_oscilloscope():
    """
    Search and connect to the OWON HDS320S Handheld Oscilloscope. 
    Set up the oscilloscope if found, otherwise attempt to reconnect.

    Note:
        This function assumes the oscilloscope's vendor and product IDs are known.
    """
    # Make global variables accessible
    global oscilloscope_OUT
    global oscilloscope_IN
    global oscilloscope

    if 'oscilloscope' in globals():
        print_to_terminal('HDS320S previously found and connected...')
        try:
            setup_oscilloscope()
        except Exception as e:
            print_to_terminal('HDS320S failed setup:')
            print_to_terminal({e})
            print_to_terminal('Clearing connection...')
            if 'oscilloscope' in globals():
                 # Delete the variable
                 del globals()['oscilloscope']
            if 'oscilloscope_IN' in globals():     
                 # Delete the variable
                 del globals()['oscilloscope_IN']
            if 'oscilloscope_OUT' in globals():     
                 # Delete the variable
                 del globals()['oscilloscope_OUT']
        return
    else:
        oscilloscope = usb.core.find(idVendor=0x5345, idProduct=0x1234)  
        if oscilloscope is None:
            print_to_terminal('HDS320S not found, please try again...')
            if 'oscilloscope' in globals():
                 # Delete the variable
                 del globals()['oscilloscope']
            if 'oscilloscope_IN' in globals():     
                 # Delete the variable
                 del globals()['oscilloscope_IN']
            if 'oscilloscope_OUT' in globals():     
                 # Delete the variable
                 del globals()['oscilloscope_OUT']

            return
        else:
            print_to_terminal('HDS320S found!')
            oscilloscope.set_configuration()
            
            # Get an endpoint instance
            config = oscilloscope.get_active_configuration()
            intf = config[(0, 0)]

            oscilloscope_OUT = usb.util.find_descriptor(
                intf,
                custom_match=lambda e: usb.util.endpoint_direction(e.bEndpointAddress) == usb.util.ENDPOINT_OUT
            )
            assert oscilloscope_OUT is not None

            oscilloscope_IN = usb.util.find_descriptor(
                intf,
                custom_match=lambda e: usb.util.endpoint_direction(e.bEndpointAddress) == usb.util.ENDPOINT_IN
            )
            assert oscilloscope_IN is not None

    # Print endpoint information for debugging
    print_to_terminal(oscilloscope_IN)
    print_to_terminal(oscilloscope_OUT)
            
    # Call the oscilloscope setup function
    try:
        setup_oscilloscope()
    except Exception as e:
        traceback.print_exc()
        print_to_terminal('HDS320S failed setup:')
        print_to_terminal({e})
        print_to_terminal('Clearing connection...')
        if 'oscilloscope' in globals():
             # Delete the variable
             del globals()['oscilloscope']
        if 'oscilloscope_IN' in globals():     
             # Delete the variable
             del globals()['oscilloscope_IN']
        if 'oscilloscope_OUT' in globals():     
             # Delete the variable
             del globals()['oscilloscope_OUT']

def get_amplitude_scale(probe_attenuation_ratio, channel, channel_in, amplitude_scales):
    # Define a mapping for each channel and ratio to the corresponding index
    if channel_in == 'CH1':
        mapping = {
            ('CH1', '1X'): 4,
            ('CH2', '1X'): 5,
            ('CH1', '10X'): 6,
            ('CH2', '10X'): 7,
        }
    if channel_in == 'CH2':
        mapping = {
            ('CH1', '1X'): 5,
            ('CH2', '1X'): 4,
            ('CH1', '10X'): 7,
            ('CH2', '10X'): 6,
        }
    try:
        index = mapping[(channel, probe_attenuation_ratio)]
    except KeyError:
        raise ValueError(f"Invalid combination: {channel} with attenuation {probe_attenuation_ratio}")

    return amplitude_scales[probe_attenuation_ratio]['commands'][index]

def setup_oscilloscope():
    """
    Configure the OWON HDS320S Handheld Oscilloscope and Arbitrary Waveform Generator (AWG).

    This function performs general device configuration, sets up channels, and initializes the AWG.
    """
    global LogFile, JSONLogFile, channel_in, channel_out
    global CH1_amplitude_scales, CH2_amplitude_scales

    # Programmatically select the "Terminal" tab
    select_tab('terminal_tab')

    # General device configuration
    print_to_terminal('\n --- Oscilloscope configurations --- \n')
    print_to_terminal(oscilloscope_query('*IDN?') + '\n')

    # Set up channels
    oscilloscope_write(':{CH}:DISPlay OFF'.format(CH=channel_out))
    print_to_terminal('Channel out status: ' + oscilloscope_query(':{CH}:DISPlay?'.format(CH=channel_out)).upper())
    oscilloscope_write(':{CH}:DISPlay OFF'.format(CH=channel_in))
    print_to_terminal('Channel in status: ' + oscilloscope_query(':{CH}:DISPlay?'.format(CH=channel_in)).upper())

    # Set AC coupling
    oscilloscope_write(':{CH}:COUPling {data}'.format(CH='CH1', data=CH1_coupling))
    print_to_terminal('Channel 1 coupling: ' + oscilloscope_query(':{CH}:COUPling?'.format(CH='CH1')).upper())
    oscilloscope_write(':{CH}:COUPling {data}'.format(CH='CH2', data=CH2_coupling))
    print_to_terminal('Channel 2 coupling: ' + oscilloscope_query(':{CH}:COUPling?'.format(CH='CH2')).upper())

    # Set attenuation mode
    oscilloscope_write(':{CH}:PROBe {data}'.format(CH='CH1', data=CH1_probe_attenuation_ratio))
    print_to_terminal('Channel 1 probe attenuation ratio: ' + oscilloscope_query(':{CH}:PROBe?'.format(CH='CH1')).upper())
    oscilloscope_write(':{CH}:PROBe {data}'.format(CH='CH2', data=CH2_probe_attenuation_ratio))
    print_to_terminal('Channel 2 probe attenuation ratio: ' + oscilloscope_query(':{CH}:PROBe?'.format(CH='CH2')).upper())
    
    # Turn on frequency and amplitude peak-to-peak measurements
    oscilloscope_write(':MEASurement:DISPlay ON')

    # Set acquisition mode
    oscilloscope_write(':ACQuire:MODE {}'.format(Sample_command))
    print_to_terminal('Acquisition mode: ' + oscilloscope_query(':ACQuire:MODE?').upper())

    # Set memory depth
    oscilloscope_write(':ACQuire:DEPMEM {}'.format(DEPMEM))
    print_to_terminal('Memory depth: ' + oscilloscope_query(':ACQuire:DEPMEM?').upper())


    oscilloscope_write(':TRIGger:SINGle:SWEEp AUTO')
    oscilloscope_write(':TRIGger:SINGLe:SOURce {CH}'.format(CH=channel_in))
    oscilloscope_write(':TRIGger:SINGle:COUPling AC')

    # Set the trigger to rising edge, VERY IMPORTANT!
    EdgeValue = oscilloscope_query(':TRIGger:SINGle:EDGe?')
    if EdgeValue == 'RISE':
        print_to_terminal('TRIGer EDGe: RISE', color=[0, 255, 0])  # Green
    elif EdgeValue == 'FALL':
        print_to_terminal('TRIGer EDGe: FALL', color=[255, 102, 102])  # Red
    else:
        print_to_terminal(f'TRIGer EDGe: {EdgeValue}', color=[255, 255, 255])  # Default (white)

    print_to_terminal('TRIGer STATus: ' + oscilloscope_query(':TRIGger:STATus?'))

    # Setup the AWG
    print_to_terminal('\n --- AWG configurations --- \n')

    # Turn off the AWG
    oscilloscope_write(':CHANnel OFF')
    print_to_terminal('Channel status: ' + oscilloscope_query(':CHANnel?').upper())

    # Set the output waveform: for Bode plots, the sine waveform is used
    oscilloscope_write(':FUNCtion SINE')
    print_to_terminal('Output waveform: ' + oscilloscope_query(':FUNCtion?').upper())

    # Set the waveform amplitude
    oscilloscope_write(':FUNCtion:AMPLitude {}'.format(waveform_amplitude_V))
    print_to_terminal('Waveform amplitude (pk-pk): ' + str(float(str(oscilloscope_query(':FUNCtion:AMPLitude?'))[0:8])) + 'V')

    # Set output impedance
    oscilloscope_write(':FUNCtion:LOAD {}'.format(AWG_output_impedance))
    print_to_terminal('High Output impedance?: ' + oscilloscope_query(':FUNCtion:LOAD?').upper() + ' Ω')

    # Set the waveform offset to zero
    oscilloscope_write(':FUNCtion:OFFSet 0')
    print_to_terminal('Waveform offset: ' + str(float(str(oscilloscope_query(':FUNCtion:OFFSet?'))[0:8])) + 'V\n')
    print_to_terminal('Now adjust both the vertical and horizontal scales...')

    # Set the vertical offset of both channels to zero
    oscilloscope_write(':CH1:OFFSet 0')
    oscilloscope_write(':CH2:OFFSet 0')
    
    # Turn on the device at the correct initial range
    oscilloscope_write(':{CH}:DISPlay ON'.format(CH='CH1'))
    oscilloscope_write(':{CH}:DISPlay ON'.format(CH='CH2'))
    oscilloscope_write(':CHANnel ON')
    
    CHin_amplitude_scales = get_amplitude_scale(CH1_probe_attenuation_ratio, channel_in, channel_in, amplitude_scales)
    CHout_amplitude_scales = get_amplitude_scale(CH2_probe_attenuation_ratio, channel_out, channel_in, amplitude_scales)
   
    time.sleep(sample_delay_s)
    set_amplitude_scale(channel_in, CHin_amplitude_scales)
    time.sleep(sample_delay_s)
    set_amplitude_scale(channel_out, CHout_amplitude_scales)
    time.sleep(sample_delay_s)
    channel_in_scale = oscilloscope_query(':{}:SCALe?'.format(channel_in))
    channel_out_scale = oscilloscope_query(':{}:SCALe?'.format(channel_out))

    print_to_terminal("channel_in: " + channel_in + " CHin_amplitude_scales: " + str(CHin_amplitude_scales))
    print_to_terminal("channel_out: " + channel_out + " CHout_amplitude_scales: " + str(CHout_amplitude_scales))
    
    AWG_first_set_frequency(decades_list[start_decade], True)

    MeasurementFREQ = oscilloscope_query(':MEASurement:{CH}:FREQuency?'.format(CH=channel_in))
    print_to_terminal('Channel out FREQuency: ' + MeasurementFREQ + " Calculate Value: " + str(round(FREQ_out_to_float(MeasurementFREQ), 0)) + '\n')



    
    if round(FREQ_out_to_float(MeasurementFREQ), 0) < 0.1:
        print_to_terminal('Issue with calculated channel out FREQuency.', color=[255, 102, 102])
        print_to_terminal('Recommend pressing Auto button on the HDS320S.', color=[255, 102, 102])
    elif EdgeValue == 'FALL':
        print_to_terminal('Issue Trigger is set to FALL needs to be set to RISE.', color=[255, 102, 102])
        print_to_terminal('Recommend pressing Trig -> F4 -> F1 to set to RISE.', color=[255, 102, 102])
    else:
        # Enable measurement button
        dpg.configure_item(item='START_MEASURE', enabled=True)
        dpg.bind_item_theme('START_MEASURE', 'GreenButton')
        dpg.bind_item_theme('SEARCH_OSCILLOSCOPE', 'DisabledButton')
        if LogFile == 'AUTO.csv' or os.path.isfile(LogFile):
            LogFile = create_filename_from_timestamp("HDS320S_Plotter", "csv")
            LogFileGUIlabel = "CSV Log File: " + LogFile
            dpg.configure_item('LOGFILE', label=LogFileGUIlabel)
        if JSONLogFile == 'AUTO.json' or os.path.isfile(JSONLogFile):
            JSONLogFile = create_filename_from_timestamp("HDS320S_Plotter", "json")
            JSONLogFileGUIlabel = "JSON Log File: " + JSONLogFile
            dpg.configure_item('JSONLOGFILEtag', label=JSONLogFileGUIlabel)
        # Disable all the control items
        control_items = [
            'Play', 'Stop'
        ]
        for item in control_items:
            dpg.configure_item(item=item, enabled=False)
        dpg.bind_item_theme('Play', 'DisabledButton')
        dpg.bind_item_theme('Stop', 'DisabledButton')

def start_mesurement():
    global LogFile, JSONLogFile, nWaveSamples, FTTcorection, is_recording
    global magnitude_phase_conn
    global xy_conn
    global sc_conn   
    global fft_oscilloscope_conn, addpm
    global magnitude_phase_stop_event, fft_oscilloscope_stop_event, xy_stop_event, sc_stop_event
    global magnitude_phase_ready_event, fft_oscilloscope_ready_event, xy_ready_event, sc_ready_event
    global thread_magnitude_phase, thread_fft_oscilloscope, thread_xy_oscilloscope, thread_sc_oscilloscope

    magnitude_phase_stop_event.clear()
    fft_oscilloscope_stop_event.clear()
    xy_stop_event.clear()
    sc_stop_event.clear()
    magnitude_phase_ready_event.clear()
    fft_oscilloscope_ready_event.clear()
    xy_ready_event.clear()
    sc_ready_event.clear()

    # Initialize plot managers and socket
    magnitude_phase_conn = None
    fft_oscilloscope_conn = None
    xy_conn = None
    sc_conn = None
    dpg.bind_item_theme('START_MEASURE', 'YellowButton')
    
    
    try:
        if plot_win_disposition == 'MatPlotLib':
            # Set environment variables for decades and points per decade
            os.environ["START_DECADE"] = str(decades_list[start_decade])
            os.environ["STOP_DECADE"] = str(decades_list[stop_decade])
            os.environ["POINTS_PER_DECADE"] = str(points_per_decade)
            
            # Start socket threads
            start_socket_threads()
 
 
            # Wait for a maximum of 10 seconds for all threads to start
            if not wait_for_threads(max_wait_time=40):
                logging.info("Aborting: Failed to start all threads.")
                return
            
            # Start plot managers
            
            # List to collect events of the plot processes that will be started.
            events = []
            
            # Start Magnitude/Phase (Bode Plots) if visible (chart_visibility[1])
            if chart_visibility[1]:
                plot_manager = PlotManagerMagnitudePhase()
                plot_manager.start_plot_process(decades_list[start_decade],
                                                decades_list[stop_decade],
                                                points_per_decade)
                events.append(magnitude_phase_ready_event)
            else:
                plot_manager = None
            
            # Start XY Oscilloscope if visible (chart_visibility[0])
            if chart_visibility[0]:
                xy_oscilloscope = XYoscilloscope()
                xy_oscilloscope.start_plot_process()
                events.append(xy_ready_event)
            else:
                xy_oscilloscope = None
            
            # Start FFT/Oscilloscope if visible (chart_visibility[2])
            if chart_visibility[2]:
                plot_manager_fft_oscilloscope = PlotManagerFFToscilloscope()
                plot_manager_fft_oscilloscope.start_plot_process(decades_list[start_decade],
                                                                decades_list[stop_decade],
                                                                points_per_decade)
                events.append(fft_oscilloscope_ready_event)
            else:
                plot_manager_fft_oscilloscope = None
            
            # Start Smith Chart if visible (chart_visibility[3])
            if chart_visibility[3]:
                plot_manager_sc = PlotManagerSmithChart()
                plot_manager_sc.start_plot_process()
                events.append(sc_ready_event)
            else:
                plot_manager_sc = None

            # Define a timeout duration (in seconds) for the servers to signal readiness.
            timeout_duration = 180  # seconds

            if events:
                logging.info("Waiting for all requested plot servers to be ready...")
                all_ready = wait_for_all_events(events, timeout_duration)
            
                # Check each event individually if the corresponding plot was started:
                if chart_visibility[1]:
                    if magnitude_phase_ready_event.is_set():
                        logging.info("Magnitude/Phase server is ready.")
                    else:
                        logging.info("Error: Magnitude/Phase server failed to start within the timeout.")
            
                if chart_visibility[0]:
                    if xy_ready_event.is_set():
                        logging.info("XY Oscilloscope server is ready.")
                    else:
                        logging.info("Error: XY Oscilloscope server failed to start within the timeout.")
            
                if chart_visibility[2]:
                    if fft_oscilloscope_ready_event.is_set():
                        logging.info("FFT Oscilloscope server is ready.")
                    else:
                        logging.info("Error: FFT Oscilloscope server failed to start within the timeout.")
            
                if chart_visibility[3]:
                    if sc_ready_event.is_set():
                        logging.info("SC (Smith Chart) server is ready.")
                    else:
                        logging.info("Error: SC (Smith Chart) server failed to start within the timeout.")
            else:
                logging.info("No plot managers were requested for MatPlotLib display.")


        else:
            # If not using Matplotlib, disable plot managers
            plot_manager = None
            plot_manager_fft_oscilloscope = None
            xy_oscilloscope = None
            plot_manager_sc = None
            magnitude_phase_conn = None
            fft_oscilloscope_conn = None
            xy_conn = None
            sc_oscilloscope_conn = None
            sc_conn = None
            
        if is_recording == False:
            is_recording = True
            # Enable all the control items
            control_items = [
                'Stop',
            ]
            for item in control_items:
                dpg.configure_item(item=item, enabled=True)
            dpg.bind_item_theme('Stop', 'RedButton')
            # Launch the background thread for processing
            processing_start_mesurement_threaded = threading.Thread(
                target=start_mesurement_threaded, 
                args=(plot_win_disposition, magnitude_phase_conn, fft_oscilloscope_conn, xy_conn, sc_conn),
                daemon=True
            )
            processing_start_mesurement_threaded.start()
    except Exception as e:
        logging.error(f"Error in start_mesurement: {e}")

def choose_phase(candidate_phase, current_phase, threshold):
    """Return candidate_phase if its difference from current_phase is below threshold, else None."""
    return candidate_phase if abs(candidate_phase - current_phase) < threshold else None

def is_between(val, bound1, bound2):
    """Checks whether val is between bound1 and bound2 (inclusive), regardless of their order."""
    if bound1 <= bound2:
        return bound1 <= val <= bound2
    else:
        return bound2 <= val <= bound1

def start_mesurement_threaded(plot_win_disposition,
                              magnitude_phase_conn, 
                              fft_oscilloscope_conn, xy_conn, sc_conn):
    """
    Start the measurement process for the oscilloscope and Arbitrary Waveform Generator (AWG).
    
    This function handles the UI setup, data clearing, and measurement initiation.
    """
    global LogFile, JSONLogFile, nWaveSamples, FTTcorection, is_recording, addpm, sample_delay_s
    
    OKrun = True
                                   
    # Programmatically select the "Data Table" tab
    select_tab('data_table_tab')
    
    # Clear all the plots
    try:
        # dpg.delete_item(item='MAG_SCATTER')
        dpg.delete_item(item='MAG_LINE')
        # dpg.delete_item(item='PHASE_SCATTER')
        dpg.delete_item(item='PHASE_LINE')
    except Exception as e:
        logging.info(f"Exception clearing plots: {e}")
    
    # Disable all the control items
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'START_MEASURE', 'JSONLOGFILEtag',
        'SAMPLE_WAVES', 'FFToffset'
    ]
    # if platform.system() != 'Windows':
    control_items.append('PLOT_WIN_SETTING')

    for item in control_items:
        dpg.configure_item(item=item, enabled=False)

    # Disable the checkboxes while running the plots/charts:
    for i in range(len(chart_names)):
        dpg.configure_item(f"checkbox_{i}", enabled=False)
        
    # dpg.bind_item_theme('START_MEASURE', 'DisabledButton')
    dpg.bind_item_theme('JSONLOGFILEtag', 'DisabledButton')
    
    # Clear all the data containers
    raw_frequencies_range.clear()
    gain_Y.clear()
    gain_X.clear()
    phase_Y.clear()
    phase_X.clear()
    MxFFTx.clear()
    MxFFTy.clear()
    
    # Set Magnitude/Phase/FFT plot Max range
    dpg.set_axis_limits('MAG_X', decades_list[start_decade], decades_list[stop_decade])
    dpg.set_axis_limits('PHASE_X', decades_list[start_decade], decades_list[stop_decade])
    dpg.set_axis_limits('FFT_X', decades_list[start_decade], decades_list[stop_decade])
    
    try:
        # Set the vertical offset of both channels to zero
        oscilloscope_write(':CH1:OFFSet 0')
        oscilloscope_write(':CH2:OFFSet 0')
        
        # Generate the complete test frequency range
        if dpg.get_value(item='POINTS_SPACING') == points_spacing[0]:
            for indx in range(start_decade, stop_decade, 1):
                current_frequency_range = np.linspace(decades_list[indx], decades_list[indx + 1], num=points_per_decade)
                for indy in current_frequency_range:
                    raw_frequencies_range.append(indy)
                # Remove the last, duplicated value
                if indx != (stop_decade - 1):
                    raw_frequencies_range.pop()
        else:
            for indx in range(start_decade, stop_decade, 1):
                current_frequency_range = np.geomspace(decades_list[indx], decades_list[indx + 1], num=points_per_decade, endpoint=True)
                for indy in current_frequency_range:
                    raw_frequencies_range.append(indy)
                # Remove the last, duplicated value
                if indx != (stop_decade - 1):
                    raw_frequencies_range.pop()
    
        # Cleanup old data from table in GUI
        dpg.delete_item(item='DataTable')
        with dpg.table(parent='dataTableWindow', header_row=True, resizable=True, policy=dpg.mvTable_SizingStretchProp, width=setting_window_width - 5, height=318, freeze_rows=1,
                       scrollY=True, scrollX=False, borders_outerH=True, borders_innerV=True, borders_innerH=True, borders_outerV=True, tag='DataTable'):
            dpg.add_table_column(label="Frequency")
            dpg.add_table_column(label="MeasFreq")
            dpg.add_table_column(label="Voltage")
            dpg.add_table_column(label="Phase")
            dpg.add_table_column(label="FFT Hz")
            dpg.add_table_column(label="FFT Meg")
        
        # Create a dictionary to hold input and output waves for each frequency
        data = {}
       
        # Open output CSV file for data writes.
        f = open(LogFile,'w')
        f.write('Frequency,MeasFreq,VpkpkMeter,Phase,"FFT Hz","FFT Meg"\n') 
       
        # Measurement loop
        # PhaseLast is defined and initialized before entering the loop to zero.
        PhaseLast = 0
        sign = 1
        deltas = deque(maxlen=3)  # Track the last 3 deltas
        dpg.bind_item_theme('START_MEASURE', 'DisabledButton')
        FreqAdjustOffset = 0
        AWG_first_set_frequency(raw_frequencies_range[0])
        
        for index, frequency in enumerate(raw_frequencies_range):
            # Stop-measurement check
            if not is_recording:
                break
            PhaseList = []
            # Set the target frequency once
            MeasurementFREQ, NewPhases = AWG_set_frequency(frequency)
            PhaseList.extend(NewPhases)
            if index == 0 and PhaseLast == 0:
                if len(PhaseList) != 0:
                    PhaseLast = sum(PhaseList) / len(PhaseList)
            
            # ---------------------------
            # PERFORM MEASUREMENT
            # ---------------------------
            # Query current scales
            current_v_scaleout = vertical_scale_to_float(oscilloscope_query(':{}:SCALe?'.format(channel_out)))
            current_v_scalein = vertical_scale_to_float(oscilloscope_query(':{}:SCALe?'.format(channel_in)))
            current_h_scale = horizontal_scale_to_float(oscilloscope_query(':HORIzontal:SCALe?'))
            # Compute the time range for the waveform acquisition
            raw_frequencies_x = np.linspace(-(horizontal_decades_n * current_h_scale), (horizontal_decades_n * current_h_scale), 300)
             
            # Ask for the complete datapoint array of the input channel
            list_of_in_arrays = []
            for _ in range(nWaveSamples):
                readPk = 0
                while readPk == 0:
                    raw_waveform_in = get_waveform(channel_in, current_v_scalein)
                    readPk = peak_to_peak(raw_waveform_in)
                list_of_in_arrays.append(raw_waveform_in)
            raw_waveform_in = average_sinusoidal_arrays(list_of_in_arrays)
            sampling_PERIOD = 1 / (current_h_scale / (300 / (horizontal_decades_n * 2)))            
            raw_frequencies_x_in = []
            for timestep in range(len(raw_waveform_in)):
                raw_frequencies_x_in.append(timestep * sampling_PERIOD)

            # Ask for the complete datapoint array of the output channel
            list_of_out_arrays = []
            for _ in range(nWaveSamples):
                readPk = 0
                attempts = 0
                max_attempts = 3
                while readPk == 0 and attempts < max_attempts:
                    raw_waveform_out = get_waveform(channel_out, current_v_scaleout)
                    readPk = peak_to_peak(raw_waveform_out)
                    attempts += 1
                list_of_out_arrays.append(raw_waveform_out)
            raw_waveform_out = average_sinusoidal_arrays(list_of_out_arrays)
            raw_frequencies_x_out = []
            for timestep in range(len(raw_waveform_out)):
                raw_frequencies_x_out.append(timestep * sampling_PERIOD)
            
            # Compute full sine analysis (FSA)
            FSA = full_sin_analysis(sampling_PERIOD, FTTcorection, raw_waveform_in, raw_waveform_out)
            FFT_freq = FSA["fft_frequencies"][:len(FSA["fft_frequencies"]) // 2]
            FFT_result = FSA["amplitude_db_signal2"][:len(FSA["fft_frequencies"]) // 2]
            FFT_max = FSA['max_amplitude_db_signal2']
            FFT_maxfreq = FSA['angular_frequency_signal2'] / (2 * np.pi)
            
            # Attempt candidate phases in a prioritized list.
            # Candidate 1: Direct FFT measurement.
            candidate_fft = FSA['phase_difference_degrees']
            new_phase = choose_phase(candidate_fft, PhaseLast, addpm)
            
            # Candidate 2: Zero-cross measurement (after normalizing).
            if new_phase is None:
                candidate_zero_cross = normalize_phase(FSA['average_delta_time'] * MeasurementFREQ * 360)
                new_phase = choose_phase(candidate_zero_cross, PhaseLast, addpm)
                
            # Candidate 3: Average phase from PhaseList if available.
            if new_phase is None and PhaseList:
                if len(PhaseList) >= 2:
                    candidate_avg = (PhaseList[-1] + PhaseList[-2]) / 2.0
                else:
                    candidate_avg = PhaseList[-1]
                new_phase = choose_phase(candidate_avg, PhaseLast, addpm)

            # Fallback logic: if no candidate was acceptable from current measurement.
            if new_phase is None:
                # Safely get the next frequency from raw_frequencies_range
                next_frequency = raw_frequencies_range[index + 1] if (index + 1) < len(raw_frequencies_range) else None
                if next_frequency is not None:
                    # Measure the next frequency (using the same AWG_set_frequency call)
                    MeasurementFREQ_next, NextNewPhases = AWG_set_frequency(next_frequency)
                    if NextNewPhases:
                        if len(PhaseList) >= 2:
                            next_candidate_avg = (NextNewPhases[-1] + NextNewPhases[-2]) / 2.0
                        else:
                            next_candidate_avg = NextNewPhases[-1]
                        # Check if the current fallback (candidate_fft) lies between PhaseLast and the next candidate average.
                        if is_between(candidate_fft, PhaseLast, next_candidate_avg):
                            new_phase = candidate_fft
                        else:
                            # Split the difference between PhaseLast and the next candidate average.
                            new_phase = (PhaseLast + next_candidate_avg) / 2.0
                    else:
                        new_phase = candidate_fft  # Fallback to candidate_fft if no NextNewPhases
                else:
                    new_phase = candidate_fft  # Use candidate_fft if there's no next frequency
            
            # Update PhaseLast and Phase_Diff_degree based on selected candidate.
            PhaseLast = new_phase
            Phase_Diff_degree = new_phase

            # --- Update after data read from oscilloscope ---
            OSCvRange = 1
            OSCvRange1 = 1
            OSCvRange2 = 1
            
            VpkpkMeterOut = vertical_scale_to_float(get_pkpk_voltage(channel_out))
            VpkpkMeterIn = vertical_scale_to_float(get_pkpk_voltage(channel_in))

            # Optionally handle the case where a command wasn't found
            if VpkpkMeterOut is None or VpkpkMeterIn is None:
                logging.warning(
                    "One of the DMrange values was not found in amplitude_scales_commands: VpkpkMeterIn=%s, VpkpkMeterOut=%s",
                    VpkpkMeterIn,
                    VpkpkMeterOut
                )
            else:
                OSCvRange = max(VpkpkMeterIn, VpkpkMeterOut)
             
            graph_processing(OSCvRange,
                             FFT_maxfreq,
                             FFT_max,
                             FFT_freq,
                             FFT_result,
                             raw_frequencies_x_in,
                             raw_waveform_in,
                             raw_frequencies_x_out,
                             raw_waveform_out)                   
            

            if VpkpkMeterIn > 0:          
                gY = (20 * np.log10(VpkpkMeterOut / VpkpkMeterIn)) if VpkpkMeterOut > 0 else float('-inf')
            else:
                gY = (20 * np.log10(VpkpkMeterOut / waveform_amplitude_V)) if VpkpkMeterOut > 0 else float('-inf')
            gX = MeasurementFREQ
            pY = Phase_Diff_degree
            pX = frequency
            
            graph_processing_Bode(gX, gY, pX, pY)

            # Convert the array to a list
            list_FFT_freq = check_and_convert_to_list(FFT_freq)
            list_FFT_result = check_and_convert_to_list(FFT_result)          
            list_raw_frequencies_x_in = check_and_convert_to_list(raw_frequencies_x_in)
            list_raw_waveform_in = check_and_convert_to_list(raw_waveform_in)
            list_raw_frequencies_x_out = check_and_convert_to_list(raw_frequencies_x_out)
            list_raw_waveform_out = check_and_convert_to_list(raw_waveform_out)     
            data[frequency] = {
                            "OSCvRange": OSCvRange,
                            "FFT_maxfreq": FFT_maxfreq,
                            "FFT_max": FFT_max,
                            "FFT_freq": list_FFT_freq,
                            "FFT_result": list_FFT_result,                   
                            "raw_frequencies_x_in": list_raw_frequencies_x_in,
                            "raw_waveform_in": list_raw_waveform_in,
                            "raw_frequencies_x_out": list_raw_frequencies_x_out,
                            "raw_waveform_out": list_raw_waveform_out,
                            "VpkpkMeter": VpkpkMeterOut,
                            "gX": gX,
                            "gY": gY,
                            "pX": pX,
                            "pY": pY
            }
            
            if plot_win_disposition == 'MatPlotLib':
                if magnitude_phase_conn:
                    # Send data to Magnitude/Phase PlotManager socket
                    if all(key in data[frequency] and data[frequency][key] is not None for key in ['gX', 'gY', 'pX', 'pY']):
                        try:
                            payload = json.dumps({
                                "gain_X": round(data[frequency]['gX']),  # Single value, rounded to whole number
                                "gain_Y": round(data[frequency]['gY'], 4),  # Single value, rounded to 4 decimal places
                                "phase_X": round(data[frequency]['pX']),  # Single value, rounded to whole number
                                "phase_Y": round(data[frequency]['pY'], 4)  # Single value, rounded to 4 decimal places
                            })

                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 2560):  # Ensure `range` is not shadowed
                                chunk = payload[start:start+2560]
                                magnitude_phase_conn.sendall(chunk.encode('utf-8'))

                            # Send end-of-message marker
                            magnitude_phase_conn.sendall(b"<END>")  # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to Magnitude/Phase plot: {e}")
                    else:
                        logging.error("Magnitude/Phase connection is None. Skipping data transmission.")
                    
                if fft_oscilloscope_conn:  
                    # Send data to FFT Oscilloscope PlotManager socket
                    if all(key in data[frequency] and data[frequency][key] is not None for key in ['FFT_maxfreq', 'FFT_max', 'FFT_freq', 'FFT_result', 'raw_frequencies_x_in', 'raw_waveform_in', 'raw_frequencies_x_out', 'raw_waveform_out']):
                        try:
                            payload = json.dumps({
                                "FFTxMax": round(data[frequency]['FFT_maxfreq']),  # No decimals
                                "FFTmaxVal": round(data[frequency]['FFT_max'], 4),  # Round to whole number
                                "FFTx": [round(x) for x in data[frequency]['FFT_freq']],  # Round each element to 4 decimals
                                "FFTy": [round(y, 4) for y in data[frequency]['FFT_result']],  # Round each element to whole numbers
                                "OSCxin": [round(x) for x in data[frequency]['raw_frequencies_x_in']],  # Round to 4 decimals
                                "OSCyin": [round(y, 4) for y in data[frequency]['raw_waveform_in']],  # Round each element to whole numbers
                                "OSCxout": [round(x) for x in data[frequency]['raw_frequencies_x_out']],  # Round to 4 decimals
                                "OSCyout": [round(y, 4) for y in data[frequency]['raw_waveform_out']]  # Round each element to whole numbers
                            })
                            
                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 32768):  # Ensure `range` is not shadowed
                                chunk = payload[start:start+32768]
                                fft_oscilloscope_conn.sendall(chunk.encode('utf-8'))              
                            # Send end-of-message marker
                            fft_oscilloscope_conn.sendall(b"<END>")  # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to FFT Oscilloscope plot: {e}")
                
                if xy_conn:
                    # Send data to XY Oscilloscope PlotManager socket
                    if all(key in data[frequency] and data[frequency][key] is not None for key in ['raw_waveform_in', 'raw_waveform_out']):
                        try:
                            payload = json.dumps({
                                "raw_waveform_in": [round(value, 4) for value in data[frequency]['raw_waveform_in']],  # Round to 4 decimal places
                                "raw_waveform_out": [round(value, 4) for value in data[frequency]['raw_waveform_out']]  # Round to 4 decimal places
                            })
                            
                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 20480):  # Ensure `range` is not shadowed
                                chunk = payload[start:start+20480]
                                xy_conn.sendall(chunk.encode('utf-8'))

                            # Send end-of-message marker
                            xy_conn.sendall(b"<END>")  # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to XY Oscilloscope plot: {e}")
                    else:
                        logging.error("XY Oscilloscope connection is None. Skipping data transmission.")

                if sc_conn:
                    # Send data to Smith Chart PlotManager socket
                    if (all(key in data[frequency] and data[frequency][key] is not None 
                            for key in ['gX', 'gY', 'pX', 'pY']) ):
                        try:
                            payload = json.dumps({
                                "frequency": [round(data[frequency]['pX'])],  # Wrap in a list
                                "magnitude_db": [round(data[frequency]['gY'], 4)],  # Wrap in a list
                                "phase_degrees": [round(data[frequency]['pY'], 4)]  # Wrap in a list
                            })    
                            
                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 2560):  # Ensure `range` is not shadowed
                                chunk = payload[start:start+2560]
                                sc_conn.sendall(chunk.encode('utf-8'))

                            # Send end-of-message marker
                            sc_conn.sendall(b"<END>")  # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to Smith Chart plot: {e}")
                    else:
                        logging.error("Smith Chart connection is None. Skipping data transmission.")
                        
            with dpg.table_row(parent='DataTable'):
                dpg.add_text(str(int(round(pX, 0))))
                dpg.add_text(str(int(round(gX, 0))))
                dpg.add_text(str(round(VpkpkMeterOut, 3)))
                dpg.add_text(str(round(pY, 3)))
                dpg.add_text(str(int(round(FFT_maxfreq, 0))))
                dpg.add_text(str(round(FFT_max, 2)))
            
            # Scroll to the bottom
            num_rows = len(dpg.get_item_children('DataTable', slot=1))
            scroll_pos = num_rows * 25  # Adjust the multiplier if row height is different
            dpg.set_y_scroll('DataTable', scroll_pos)
            
            f.write(str(int(round(pX, 0))))
            f.write(",")
            f.write(str(int(round(gX, 0))))
            f.write(",")
            f.write(str(round(VpkpkMeterOut, 3)))
            f.write(",")
            # f.write(str(round(OrgPhaseCal, 3)))
            f.write(str(round(pY, 3)))
            f.write(",")
            f.write(str(int(round(FFT_maxfreq, 0))))
            f.write(",")
            f.write(str(round(FFT_max, 2)))
            f.write("\n")

    except Exception as e:
        traceback.print_exc()
        print_to_terminal('HDS320S failed to complete the Measurements...')
        OKrun = False
       
    # Reset Frequency of Oscilloscope
    AWG_first_set_frequency(decades_list[start_decade], True)
    # oscilloscope_OUT.write(':FUNCtion:FREQuency {}'.format(decades_list[start_decade]))
    
    # Re-enable all the control elements
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'JSONLOGFILEtag',
        'SAMPLE_WAVES', 'FFToffset'
    ]
    # if platform.system() != 'Windows':
    control_items.append('PLOT_WIN_SETTING')
    
    for item in control_items:
        dpg.configure_item(item=item, enabled=True)
    dpg.bind_item_theme('SEARCH_OSCILLOSCOPE', 'GreenButton')
    dpg.bind_item_theme('JSONLOGFILEtag', 'GreenButton')
    is_recording = False

    # Disable all the control items
    control_items = [
        'Stop'
    ]
    for item in control_items:
        dpg.configure_item(item=item, enabled=False)

    # Enable the checkboxes while running the plots/charts:
    for i in range(len(chart_names)):
        dpg.configure_item(f"checkbox_{i}", enabled=True)
    
    dpg.bind_item_theme('Stop', 'DisabledButton')    
    # Programmatically select the "Data Table" tab
    select_tab('data_table_tab')    
    with open(JSONLogFile, 'w') as json_file:
        json.dump(data, json_file, indent=1)
    f.close()
    JSONLogFile = 'AUTO.json'
    LogFile = 'AUTO.csv'
    LogFileGUIlabel = "CSV Log File: " + LogFile
    JSONLogFileGUIlabel = "JSON Log File: " + JSONLogFile
    # Start post-processing if run was successful
    if OKrun:
        post_processing(gain_X, gain_Y, phase_X, phase_Y)
    if plot_win_disposition == 'MatPlotLib':
        # Close all connections
        close_connection(magnitude_phase_conn, "Magnitude/Phase", 2)
        close_connection(fft_oscilloscope_conn, "FFT Oscilloscope", 2)
        close_connection(xy_conn, "XY Oscilloscope", 2)
        close_connection(sc_conn, "SC Oscilloscope", 2)
        PlayBack_cleanup()

    
def processing_thread(data, plot_win_disposition,
                      magnitude_phase_conn, fft_oscilloscope_conn, xy_conn, sc_conn):
    """Synchronized thread for handling data processing and sending to separate sockets for Matplotlib updates."""
    try:
        dpg.bind_item_theme('Play', 'GreenButton')
        for index, frequency in enumerate(list(data.keys())):
            # **Step 1: Pause Loop**
            pause_loop()

            # Handle stopping playback
            if not is_playing and not is_paused:
                break  # Exit the loop if playback stops

            # **Step 2: Process Graph Data**
            graph_processing(data[frequency]['OSCvRange'],
                             data[frequency]['FFT_maxfreq'],
                             data[frequency]['FFT_max'],
                             data[frequency]['FFT_freq'],
                             data[frequency]['FFT_result'],
                             data[frequency]['raw_frequencies_x_in'],
                             data[frequency]['raw_waveform_in'],
                             data[frequency]['raw_frequencies_x_out'],
                             data[frequency]['raw_waveform_out'])

            graph_processing_Bode(data[frequency]['gX'],
                                  data[frequency]['gY'],
                                  data[frequency]['pX'],
                                  data[frequency]['pY'])

            # **Step 3: Send Data to Sockets**
            if plot_win_disposition == 'MatPlotLib':
                if magnitude_phase_conn:
                    # Send data to Magnitude/Phase PlotManager socket
                    if all(key in data[frequency] and data[frequency][key] is not None for key in ['gX', 'gY', 'pX', 'pY']):
                        try:
                            payload = json.dumps({
                                "gain_X": round(data[frequency]['gX']),  # Single value, rounded to whole number
                                "gain_Y": round(data[frequency]['gY'], 4),  # Single value, rounded to 4 decimal places
                                "phase_X": round(data[frequency]['pX']),  # Single value, rounded to whole number
                                "phase_Y": round(data[frequency]['pY'], 4)  # Single value, rounded to 4 decimal places
                            })
                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 2560):  # Ensure `range` is not shadowed
                                chunk = payload[start:start+2560]
                                magnitude_phase_conn.sendall(chunk.encode('utf-8'))

                            # Send end-of-message marker
                            magnitude_phase_conn.sendall(b"<END>")  # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to Magnitude/Phase plot: {e}")
                    else:
                        logging.error("Magnitude/Phase connection is None. Skipping data transmission.")
                    
                if fft_oscilloscope_conn:  
                    # Send data to FFT Oscilloscope PlotManager socket
                    if all(key in data[frequency] and data[frequency][key] is not None for key in ['FFT_maxfreq', 'FFT_max', 'FFT_freq', 'FFT_result', 'raw_frequencies_x_in', 'raw_waveform_in', 'raw_frequencies_x_out', 'raw_waveform_out']):
                        try:
                            payload = json.dumps({
                                "FFTxMax": round(data[frequency]['FFT_maxfreq']),  # No decimals
                                "FFTmaxVal": round(data[frequency]['FFT_max'], 4),  # Round to whole number
                                "FFTx": [round(x) for x in data[frequency]['FFT_freq']],  # Round each element to 4 decimals
                                "FFTy": [round(y, 4) for y in data[frequency]['FFT_result']],  # Round each element to whole numbers
                                "OSCxin": [round(x) for x in data[frequency]['raw_frequencies_x_in']],  # Round to 4 decimals
                                "OSCyin": [round(y, 4) for y in data[frequency]['raw_waveform_in']],  # Round each element to whole numbers
                                "OSCxout": [round(x) for x in data[frequency]['raw_frequencies_x_out']],  # Round to 4 decimals
                                "OSCyout": [round(y, 4) for y in data[frequency]['raw_waveform_out']]  # Round each element to whole numbers
                            })
                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 32768):  # Ensure `range` is not shadowed
                                chunk = payload[start:start+32768]
                                fft_oscilloscope_conn.sendall(chunk.encode('utf-8'))              
                            # Send end-of-message marker
                            fft_oscilloscope_conn.sendall(b"<END>")  # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to FFT Oscilloscope plot: {e}")
                
                if xy_conn:
                    # Send data to XY Oscilloscope PlotManager socket
                    if all(key in data[frequency] and data[frequency][key] is not None for key in ['raw_waveform_in', 'raw_waveform_out']):
                        try:
                            payload = json.dumps({
                                "raw_waveform_in": [round(value, 4) for value in data[frequency]['raw_waveform_in']],  # Round to 4 decimal places
                                "raw_waveform_out": [round(value, 4) for value in data[frequency]['raw_waveform_out']]  # Round to 4 decimal places
                            })
                            # xy_conn.sendall(payload.encode('utf-8'))
                            # xy_conn.sendall(b"<END>")
                            
                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 20480):  # Ensure `range` is not shadowed
                                chunk = payload[start:start+20480]
                                xy_conn.sendall(chunk.encode('utf-8'))

                            # Send end-of-message marker
                            xy_conn.sendall(b"<END>")  # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to XY Oscilloscope plot: {e}")
                    else:
                        logging.error("XY Oscilloscope connection is None. Skipping data transmission.")
                        
                if sc_conn:
                    # Send data to Smith Chart PlotManager socket
                    if (all(key in data[frequency] and data[frequency][key] is not None 
                            for key in ['gX', 'gY', 'pX', 'pY']) ):
                        try:
                            payload = json.dumps({
                                "frequency": [round(data[frequency]['pX'])],  # Wrap in a list
                                "magnitude_db": [round(data[frequency]['gY'], 4)],  # Wrap in a list
                                "phase_degrees": [round(data[frequency]['pY'], 4)]  # Wrap in a list
                            })                            
                            
                            # Split serialized data into chunks and send
                            for start in range(0, len(payload), 2560):
                                chunk = payload[start:start+2560]
                                sc_conn.sendall(chunk.encode('utf-8'))

                            # Send end-of-message marker
                            sc_conn.sendall(b"<END>") # Signal the end of the current payloadfor frequency {frequency}")
                        except Exception as e:
                            logging.error(f"Error sending data to Smith Chart plot: {e}")
                    else:
                        logging.error("Either required data fields are missing/inconsistent or the connection is None. Skipping data transmission.")
                        

            # **Step 4: Update DearPyGui Table**
            try:
                with dpg.table_row(parent='DataTable'):
                    dpg.add_text(str(int(round(data[frequency]['pX'], 0))))
                    dpg.add_text(str(int(round(data[frequency]['gX'], 0))))
                    dpg.add_text(str(round(data[frequency]['VpkpkMeter'], 3)))
                    dpg.add_text(str(round(data[frequency]['pY'], 3)))
                    dpg.add_text(str(int(round(data[frequency]['FFT_maxfreq'], 0))))
                    dpg.add_text(str(round(data[frequency]['FFT_max'], 2)))

                    # Scroll to the bottom
                    num_rows = len(dpg.get_item_children('DataTable', slot=1))
                    scroll_pos = num_rows * 25
                    dpg.set_y_scroll('DataTable', scroll_pos)
            except Exception as e:
                logging.error(f"Error updating DearPyGui table: {e}")

            # **Step 5: Control Playback Speed**
            time.sleep(play_speed)

    except Exception as e:
        logging.error(f"Error in processing thread: {e}")

    # Update Bode Plots in DearPyGui
    post_processing(gain_X, gain_Y, phase_X, phase_Y)
    
    # Close all connections
    close_connection(magnitude_phase_conn, "Magnitude/Phase", 2)
    close_connection(fft_oscilloscope_conn, "FFT Oscilloscope", 2)
    close_connection(xy_conn, "XY Oscilloscope", 2)
    close_connection(sc_conn, "SC Oscilloscope", 2)
    
# Close connections gracefully
def close_connection(conn, conn_name, timeout=5):
    try:
        if conn:
            # Send the "<CLOSE>" signal to the server
            conn.sendall(b"<CLOSE>")
            logging.info(f"Close signal sent to {conn_name} socket.")

            # Wait briefly for the server to process and respond
            conn.settimeout(timeout)

            # Retry mechanism for acknowledgment
            retries = 3
            for attempt in range(retries):
                try:
                    # Listen for acknowledgment from the server
                    response = conn.recv(128).decode('utf-8')  # Buffer size of 128 bytes
                    if response == "<CLOSED>":
                        logging.info(f"Acknowledgment received: Server successfully closed the {conn_name} connection.")
                        logging.info(f"Acknowledgment received: Server successfully closed the {conn_name} connection.")
                        break
                except socket.timeout:
                    logging.warning(f"Attempt {attempt + 1}: Timeout waiting for acknowledgment from {conn_name}.")
                    logging.info(f"Attempt {attempt + 1}: Timeout waiting for acknowledgment from {conn_name}.")
                if attempt == retries - 1:
                    logging.error(f"No acknowledgment after {retries} attempts. Proceeding with socket closure.")
                    logging.info(f"No acknowledgment after {retries} attempts. Proceeding with socket closure.")
            
            # Wait briefly to allow any remaining data to flush
            time.sleep(0.5)

            # Gracefully shut down the socket
            try:
                if conn.fileno() != -1:  # Check socket validity
                    conn.shutdown(socket.SHUT_RDWR)
                    logging.info(f"{conn_name} socket shutdown completed.")
                else:
                    logging.warning(f"{conn_name} socket already invalid or closed.")
            except Exception as e:
                logging.error(f"Error during socket shutdown for {conn_name}: {e}")
    except Exception as e:
        logging.error(f"Error closing {conn_name} socket: {e}")
    finally:
        # Ensure the socket is always closed
        try:
            if conn and conn.fileno() != -1:  # Check socket validity before closing
                conn.close()
                logging.info(f"{conn_name} socket fully closed.")
            else:
                logging.warning(f"{conn_name} socket already closed. Skipping close.")
        except Exception as e:
            logging.error(f"Error during socket close for {conn_name}: {e}")
        finally:
            # Explicitly set the connection variable to None
            if conn_name == "Magnitude/Phase":
                magnitude_phase_conn = None
            elif conn_name == "FFT Oscilloscope":
                fft_oscilloscope_conn = None
            elif conn_name == "XY Oscilloscope":
                xy_conn = None
            elif conn_name == "SC Oscilloscope":
                sc_conn = None
            logging.info(f"{conn_name} global variable reset to None.")

def clear_events(*events):
    """Utility function to clear multiple threading events."""
    for event in events:
        event.clear()

def handle_socket(stop_event, ready_event, port, connection_var_name):
    """Generic socket handling function for thread management."""
    global_vars = globals()
    socket_instance = None
    try:
        socket_instance = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        socket_instance.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        # socket_instance.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEPORT, 1)
        socket_instance.bind(('127.0.0.1', port))
        socket_instance.listen(1)
        socket_instance.settimeout(4)  # Set timeout for accept calls
        logging.info(f"Waiting for connection on socket (port {port})...")

        while not stop_event.is_set():
            if global_vars[connection_var_name] is None:
                try:
                    global_vars[connection_var_name], _ = socket_instance.accept()
                    logging.info(f"Connection established on port {port}.")
                except socket.timeout:
                    pass  # Timeout to recheck stop_event
                except Exception as e:
                    logging.error(f"Error handling socket on port {port}: {e}")
            else:
                if not ready_event.is_set():
                    ready_event.set()
                time.sleep(0.01)  # Small delay between loops
    finally:
        if global_vars.get(connection_var_name):
            global_vars[connection_var_name].close()
            global_vars[connection_var_name] = None
            logging.info(f"Connection on port {port} closed.")
        if socket_instance:
            try:
                socket_instance.shutdown(socket.SHUT_RDWR)
            except Exception as e:
                logging.warning(f"Error during shutdown: {e}")
            socket_instance.close()
            logging.info(f"Socket on port {port} fully closed.")

def start_socket_threads():
    """Initialize and start threads for handling different sockets."""
    global magnitude_phase_stop_event, fft_oscilloscope_stop_event, xy_stop_event, sc_stop_event
    global magnitude_phase_ready_event, fft_oscilloscope_ready_event, xy_ready_event, sc_ready_event
    global thread_magnitude_phase, thread_fft_oscilloscope, thread_xy_oscilloscope, thread_sc_oscilloscope

    # Clear stop and ready events for all socket types.
    clear_events(magnitude_phase_stop_event, fft_oscilloscope_stop_event, xy_stop_event, sc_stop_event)
    clear_events(magnitude_phase_ready_event, fft_oscilloscope_ready_event, xy_ready_event, sc_ready_event)

    # Retrieve base port from environment variable or use default.
    try:
        base_port = int(os.getenv("BASE_PORT", 5001))
    except ValueError:
        logging.error("Invalid BASE_PORT value in environment; falling back to 5001")
        base_port = 5001

    # Assign related ports dynamically.
    MP_PORT, FFT_PORT, XY_PORT, SC_PORT = base_port, base_port + 1, base_port + 2, base_port + 3

    logging.info(f"MP_PORT: {MP_PORT}, FFT_PORT: {FFT_PORT}, XY_PORT: {XY_PORT}, SC_PORT: {SC_PORT}")

    # Start threads for each socket.
    thread_magnitude_phase = threading.Thread(
        target=handle_socket,
        args=(magnitude_phase_stop_event, magnitude_phase_ready_event, MP_PORT, "magnitude_phase_conn"),
        daemon=True
    )
    thread_fft_oscilloscope = threading.Thread(
        target=handle_socket,
        args=(fft_oscilloscope_stop_event, fft_oscilloscope_ready_event, FFT_PORT, "fft_oscilloscope_conn"),
        daemon=True
    )
    thread_xy_oscilloscope = threading.Thread(
        target=handle_socket,
        args=(xy_stop_event, xy_ready_event, XY_PORT, "xy_conn"),
        daemon=True
    )
    thread_sc_oscilloscope = threading.Thread(
        target=handle_socket,
        args=(sc_stop_event, sc_ready_event, SC_PORT, "sc_conn"),
        daemon=True
    )

    # Start all threads.
    thread_magnitude_phase.start()
    thread_fft_oscilloscope.start()
    thread_xy_oscilloscope.start()
    thread_sc_oscilloscope.start()

    logging.info("All socket threads started.")

def close_socket(socket, name):
    """Utility to safely close a socket."""
    if socket:
        try:
            socket.close()
            logging.info(f"{name} socket closed.")
        except Exception as e:
            logging.error(f"Error closing {name} socket: {e}")

def wait_for_threads(max_wait_time=10):
    """Wait for all threads to start running within a given time limit."""
    start_time = time.time()
    while time.time() - start_time < max_wait_time:
        # Check if all threads are alive
        if (
            thread_magnitude_phase.is_alive() and
            thread_fft_oscilloscope.is_alive() and
            thread_xy_oscilloscope.is_alive() and
            thread_sc_oscilloscope.is_alive()
        ):
            logging.info("All threads are running. Proceeding to start plot managers...")
            return True
        else:
            logging.info("Waiting for threads to start...")
            time.sleep(0.5)  # Small delay before rechecking
    
    # If the loop exits without all threads alive
    logging.info("Error: One or more threads failed to start within the time limit.")
    return False

def wait_for_all_events(events, timeout):
    """Wait until all events in the list are set or until timeout.

    Args:
        events (list): A list of threading.Event or multiprocessing.Event objects.
        timeout (float): Maximum time in seconds to wait for all events.

    Returns:
        bool: True if all events are set before the timeout, otherwise False.
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        if all(event.is_set() for event in events):
            return True
        time.sleep(0.05)  # Polling interval. Adjust as needed.
    return False

def PlayBack():
    global JSONLogFile
    global is_playing, is_paused, play_speed
    global magnitude_phase_conn
    global xy_conn, sc_conn
    global fft_oscilloscope_conn
    global magnitude_phase_stop_event, fft_oscilloscope_stop_event, xy_stop_event, sc_stop_event
    global magnitude_phase_ready_event, fft_oscilloscope_ready_event, xy_ready_event, sc_ready_event
    global thread_magnitude_phase, thread_fft_oscilloscope, thread_xy_oscilloscope, thread_sc_oscilloscope

    magnitude_phase_stop_event.clear()
    fft_oscilloscope_stop_event.clear()
    xy_stop_event.clear()
    sc_stop_event.clear()
    magnitude_phase_ready_event.clear()
    fft_oscilloscope_ready_event.clear()
    xy_ready_event.clear()
    sc_ready_event.clear()

    logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')

    # Disable all the control items
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'START_MEASURE', 'JSONLOGFILEtag',
        'SAMPLE_WAVES', 'FFToffset'
    ]
    # if platform.system() != 'Windows':
    control_items.append('PLOT_WIN_SETTING')
    
    for item in control_items:
        dpg.configure_item(item=item, enabled=False)
    dpg.bind_item_theme('JSONLOGFILEtag', 'DisabledButton')
    dpg.bind_item_theme('SEARCH_OSCILLOSCOPE', 'DisabledButton')
    
    # Disable the checkboxes while running the plots/charts:
    for i in range(len(chart_names)):
        dpg.configure_item(f"checkbox_{i}", enabled=False)
    
    # Enable all the control items
    control_items = [
        'Stop'
    ]
    for item in control_items:
        dpg.configure_item(item=item, enabled=True)
    dpg.bind_item_theme('Play', 'YellowButton')
    dpg.bind_item_theme('Stop', 'RedButton') 
    # Clear all the plots
    try:
        # dpg.delete_item(item='MAG_SCATTER')
        dpg.delete_item(item='MAG_LINE')
        # dpg.delete_item(item='PHASE_SCATTER')
        dpg.delete_item(item='PHASE_LINE')
    except Exception as e:
        logging.info(f"Exception clearing plots: {e}")
    
    # Clear all the data containers
    raw_frequencies_range.clear()
    gain_Y.clear()
    gain_X.clear()
    phase_Y.clear()
    phase_X.clear()
    MxFFTx.clear()
    MxFFTy.clear()

    # Create a dictionary to hold input and output waves for each frequency
    data = {}
       
    try:
        with open(JSONLogFile, 'r') as json_file:
            data = json.load(json_file)
    except FileNotFoundError:
        logging.info(f"JSONLogFile not found: {JSONLogFile}")
        return
    except json.JSONDecodeError:
        logging.info(f"Error decoding JSON file: {JSONLogFile}")
        return

    my_range = list(data.keys())
    start_decade = decades_list.index(find_closest(decades_list, my_range[0]))
    dpg.set_value('START_DEC', decades_list_string[start_decade])
    stop_decade = decades_list.index(find_closest(decades_list, my_range[-1]))
    dpg.set_value('STOP_DEC', decades_list_string[stop_decade])
    # logging.info("Start: ", decades_list[start_decade], " End: ", decades_list[stop_decade])
    # logging.info("Start: ", decades_list_string[start_decade], " End: ", decades_list_string[stop_decade])
    
    # Set Magnitude/Phase/FFT plot Max range
    dpg.set_axis_limits('MAG_X', decades_list[start_decade], decades_list[stop_decade])
    dpg.set_axis_limits('PHASE_X', decades_list[start_decade], decades_list[stop_decade])
    dpg.set_axis_limits('FFT_X', decades_list[start_decade], decades_list[stop_decade])

    # Cleanup old data from table in GUI
    dpg.delete_item(item='DataTable')
    with dpg.table(parent='dataTableWindow', header_row=True, resizable=True, policy=dpg.mvTable_SizingStretchProp, width=setting_window_width - 5, height=318, freeze_rows=1,
                    scrollY=True, scrollX=False, borders_outerH=True, borders_innerV=True, borders_innerH=True, borders_outerV=True, tag='DataTable'):
        dpg.add_table_column(label="Frequency")
        dpg.add_table_column(label="MeasFreq")
        dpg.add_table_column(label="Voltage")
        dpg.add_table_column(label="Phase")
        dpg.add_table_column(label="FFT Hz")
        dpg.add_table_column(label="FFT Meg")

    # Programmatically select the "Data Table" tab
    select_tab('data_table_tab')
    
    # Initialize plot managers and socket
    magnitude_phase_conn = None
    fft_oscilloscope_conn = None
    xy_conn = None
    sc_conn = None
    
    try:
        if plot_win_disposition == 'MatPlotLib':
            # Set environment variables for decades and points per decade
            os.environ["START_DECADE"] = str(decades_list[start_decade])
            os.environ["STOP_DECADE"] = str(decades_list[stop_decade])
            os.environ["POINTS_PER_DECADE"] = str(points_per_decade)
            
            # Start socket threads
            start_socket_threads()
 
            # Conditionally start each plot process based on chart_visibility.
            # Index mapping: 0 -> XY Plot, 1 -> Bode Plots, 2 -> FFT/Oscilloscope, 3 -> Smith Chart
            events = []

            # Bode Plots (Magnitude/Phase): chart_visibility[1]
            if chart_visibility[1]:
                plot_manager = PlotManagerMagnitudePhase()
                plot_manager.start_plot_process(decades_list[start_decade], decades_list[stop_decade], points_per_decade)
                events.append(magnitude_phase_ready_event)
            else:
                plot_manager = None

            # XY Plot: chart_visibility[0]
            if chart_visibility[0]:
                xy_oscilloscope = XYoscilloscope()
                xy_oscilloscope.start_plot_process()
                events.append(xy_ready_event)
            else:
                xy_oscilloscope = None

            # FFT/Oscilloscope: chart_visibility[2]
            if chart_visibility[2]:
                plot_manager_fft_oscilloscope = PlotManagerFFToscilloscope()
                plot_manager_fft_oscilloscope.start_plot_process(decades_list[start_decade],
                                                                 decades_list[stop_decade],
                                                                 points_per_decade)
                events.append(fft_oscilloscope_ready_event)
            else:
                plot_manager_fft_oscilloscope = None

            # Smith Chart: chart_visibility[3] (not started in your configuration)
            if chart_visibility[3]:
                plot_manager_sc = PlotManagerSmithChart()
                plot_manager_sc.start_plot_process()
                events.append(sc_ready_event)
            else:
                plot_manager_sc = None

            timeout_duration = 180  # seconds
            logging.info("Waiting for all requested plot servers to be ready...")
            
            if events:
                all_ready = wait_for_all_events(events, timeout_duration)
            
                # Verify each individually if started
                if chart_visibility[1]:
                    if magnitude_phase_ready_event.is_set():
                        logging.info("Magnitude/Phase server is ready.")
                    else:
                        logging.info("Error: Magnitude/Phase server failed to start within the timeout.")

                if chart_visibility[0]:
                    if xy_ready_event.is_set():
                        logging.info("XY Oscilloscope server is ready.")
                    else:
                        logging.info("Error: XY Oscilloscope server failed to start within the timeout.")

                if chart_visibility[2]:
                    if fft_oscilloscope_ready_event.is_set():
                        logging.info("FFT Oscilloscope server is ready.")
                    else:
                        logging.info("Error: FFT Oscilloscope server failed to start within the timeout.")

                if chart_visibility[3]:
                    if sc_ready_event.is_set():
                        logging.info("SC (Smith Chart) server is ready.")
                    else:
                        logging.info("Error: SC (Smith Chart) server failed to start within the timeout.")
            else:
                logging.info("No plot processes were requested to start.")
                
        else:
            # If not using Matplotlib, disable plot managers
            plot_manager = None
            plot_manager_fft_oscilloscope = None
            xy_oscilloscope = None
            plot_manager_sc = None
            magnitude_phase_conn = None
            fft_oscilloscope_conn = None
            xy_conn = None
            sc_conn = None
        
                
        # Launch the background thread for processing
        processing_thread_instance = threading.Thread(
            target=processing_thread,
            args=(data, plot_win_disposition, magnitude_phase_conn, fft_oscilloscope_conn, xy_conn, sc_conn),
            daemon=True
        )
        processing_thread_instance.start()
        processing_thread_instance.join()

    except Exception as e:
        logging.error(f"Error in PlayBack: {e}")
    finally:
        logging.info("Final cleanup initiated.")
        stop_play_cleaup()
        is_playing = False
        dpg.configure_item('Play', label='Play')

        if plot_win_disposition == 'MatPlotLib':
            PlayBack_cleanup()
            logging.info("Cleanup complete. Ready for next playback run.")

def PlayBack_cleanup():
    # Set stop events and close connections
    magnitude_phase_stop_event.set()
    fft_oscilloscope_stop_event.set()
    xy_stop_event.set()
    sc_stop_event.set()
    
    close_socket(magnitude_phase_conn, "Magnitude/Phase")
    close_socket(fft_oscilloscope_conn, "FFT Oscilloscope")
    close_socket(xy_conn, "XY Oscilloscope")
    close_socket(sc_conn, "SC Oscilloscope")
    logging.info("Stop events set for all threads.")

    # Wait for threads to terminate
    if thread_magnitude_phase and thread_magnitude_phase.is_alive():
        thread_magnitude_phase.join()
    if thread_fft_oscilloscope and thread_fft_oscilloscope.is_alive():
        thread_fft_oscilloscope.join()
    if thread_xy_oscilloscope and thread_xy_oscilloscope.is_alive():
        thread_xy_oscilloscope.join()
    if thread_sc_oscilloscope and thread_sc_oscilloscope.is_alive():
        thread_sc_oscilloscope.join()
        
    # Clear events for next run
    magnitude_phase_stop_event.clear()
    fft_oscilloscope_stop_event.clear()
    xy_stop_event.clear()
    sc_stop_event.clear()
    magnitude_phase_ready_event.clear()
    fft_oscilloscope_ready_event.clear()
    xy_ready_event.clear()
    sc_ready_event.clear()
    
    # Reset data containers
    raw_frequencies_range.clear()
    gain_Y.clear()
    gain_X.clear()
    phase_Y.clear()
    phase_X.clear()
    MxFFTx.clear()
    MxFFTy.clear()


def stop_play_cleaup():
    # Re-enable all the control elements
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'JSONLOGFILEtag',
        'SAMPLE_WAVES', 'FFToffset'
        ]
    # if platform.system() != 'Windows':
    control_items.append('PLOT_WIN_SETTING')

    for item in control_items:
        dpg.configure_item(item=item, enabled=True)

    # Enable the checkboxes while running the plots/charts:
    for i in range(len(chart_names)):
        dpg.configure_item(f"checkbox_{i}", enabled=True)
        
    dpg.bind_item_theme('SEARCH_OSCILLOSCOPE', 'GreenButton')
    dpg.bind_item_theme('JSONLOGFILEtag', 'GreenButton')
    control_items = [
        'Stop'
        ]
    for item in control_items:
        dpg.configure_item(item=item, enabled=False)
    dpg.bind_item_theme('Stop', 'DisabledButton') 
    
def pause_loop():
    global is_playing, is_paused
    while is_playing:
        if not is_paused:
            control_items = [
                'Stop'
                ]
            for item in control_items:
                dpg.configure_item(item=item, enabled=True)
            dpg.bind_item_theme('Stop', 'RedButton') 
            return
        else:
            control_items = [
                'Stop'
                ]
            for item in control_items:
                dpg.configure_item(item=item, enabled=False)
            dpg.bind_item_theme('Stop', 'DisabledButton')     
            control_items = [
                'Play'
                ]
            for item in control_items:
                dpg.configure_item(item=item, enabled=True)
                dpg.bind_item_theme('Play', 'GreenButton')
            time.sleep(0.1)  # Small delay to prevent busy-waiting        

def check_and_convert_to_list(variable):
    if isinstance(variable, np.ndarray):
        # logging.info("The variable is a NumPy array.")
        list_variable = variable.tolist()
        # logging.info("NumPy array converted to list.")
        return list_variable
    elif isinstance(variable, list):
         # logging.info("The variable is a list.")
         return variable
    else:
        # logging.info("The variable is not a NumPy array or a list.")
        return None

def find_closest(list_of_numbers, target_number):
    """
    Finds the value in a list closest to a target number.

    Args:
        list_of_numbers: A list of numbers.
        target_number: The target number.

    Returns:
        The value in the list closest to the target number.
    """
    return min(list_of_numbers, key=lambda x: abs(x - float(target_number)))

# Function to calculate the peak-to-peak value of a signal
def peak_to_peak(signal):
    """
    Calculates the peak-to-peak value of a signal.

    Args:
        signal: A list of numerical values representing the signal.

    Returns:
        The peak-to-peak value of the signal, which is the difference between the maximum and minimum values.
    """
    max_val = max(signal)
    min_val = min(signal)
    return max_val - min_val

def average_sinusoidal_arrays(arrays):
    """
    Averages multiple arrays of sinusoidal data, throwing out data outside
    one standard deviation, and correcting for phase shift offset.

    Args:
        arrays: A list of arrays, where each array represents sinusoidal data.

    Returns:
        The final averaged array after outlier removal and phase shift correction.
    """

   # If only one waveform is provided, simply return it.
    if len(arrays) == 1:
        return np.array(arrays[0])
    
    # Find the array with the maximum length
    max_length = max(len(arr) for arr in arrays)

    # Interpolate all arrays to the maximum length
    interpolated_arrays = [
        np.interp(np.linspace(0, 1, max_length), np.linspace(0, 1, len(arr)), arr)
        for arr in arrays
    ]

    # Convert to numpy array for easier manipulation
    data = np.array(interpolated_arrays)

    # Calculate the mean and standard deviation for each point
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)

    # Create a mask to remove data outside one standard deviation
    mask = np.abs(data - mean) <= std

    # Correct for phase shift offset
    corrected_data = np.zeros_like(data)
    for i in range(data.shape[0]):
        correlation = np.correlate(data[i], mean, mode='full')
        shift = np.argmax(correlation) - (len(mean) - 1)
        corrected_data[i] = np.roll(data[i], -shift)

    # Apply the mask to the corrected data
    corrected_data[~mask] = np.nan

    # Calculate the final average, ignoring NaN values
    final_average = np.nanmean(corrected_data, axis=0)

    return final_average

"""
When performing an FFT on a sampled sinusoidal wave with multiple wavelengths, 
the magnitude in the FFT will appear larger than the peak-to-peak value of the 
wave because the FFT calculates the energy distributed across all frequency 
components within the signal, not just the peak amplitude, and this energy is 
spread out across multiple bins depending on how well the wave's frequency 
aligns with the FFT's frequency bins, leading to a larger apparent magnitude in 
the relevant frequency bin; to correct this, you need to normalize the FFT result 
by dividing by the number of samples in your signal to get a more accurate 
representation of the amplitude at a specific frequency. 

Key points to understand this phenomenon:

    Spectral Leakage:
    If the frequency of your sine wave doesn't perfectly align with a single FFT bin, 
    the energy will be distributed across multiple bins, making the peak magnitude 
    in the correct bin appear smaller than the actual amplitude of the wave. 

Normalization:
To get a more accurate amplitude representation, you need to divide the FFT magnitude 
by the number of samples in your signal, especially when comparing it to the peak-to-peak 
value of the original wave. 

How to correct the FFT magnitude:

    Calculate the scaling factor:
    Divide the number of samples in your signal by the number of frequency 
    bins in the FFT output.

    Apply the scaling factor:
    After performing the FFT, divide each magnitude value by the calculated scaling 
    factor to get a more accurate amplitude representation. 

2025/02/12
Analyze two sinusoidal arrays where sampling starts mid wave and are not exactly integer 
periodic using Python. Determine the normalize amplitude of FFT in dB, providing Hz 
values for FFT, and maximum value with frequency. Get angular frequency, phase, 
phase difference in degrees for FFT. Get gain between input and out put in db.  
Look at graphing the imaginary VS real VS the frequency domain use def with 
dictionary containing the analysis results
full_sin_analysis(sample_rate, full_sin_analysis, signal1, signal2):
"""
def normalize_phase(angle_deg):
    """Wraps an angle in degrees to the [-180, 180] range, treating 180° as a special case."""
    """ Maybe consider looking at last phase before 180 to determine if +/-. """
    if angle_deg % 360 == 180:
        # Return 180 if the original angle is exactly 180° modulo 360.
        return 180
    else:
        # Standard normalization to [-180, 180)
        return (angle_deg + 180) % 360 - 180

def find_zero_crossings(signal):
    """Find the indices of all rising zero-crossings in the signal."""
    zero_crossings = []
    for i in range(1, len(signal)):
        if signal[i - 1] <= 0 and signal[i] > 0:  # Rising edge zero-crossing
            zero_crossings.append(i)
    return zero_crossings
    deltas = deque(maxlen=3) 
def full_sin_analysis(sample_rate, analysis_factor, signal1, signal2):
    """
    Performs a full sinusoidal analysis on two signals, including time-domain zero-crossing analysis.
    The average_delta_time is computed after adjusting for a common trigger determined by the first
    rising zero crossing of the reference signal (signal1).
    """
    # Find zero-crossings for both signals.
    zc1 = find_zero_crossings(signal1)
    zc2 = find_zero_crossings(signal2)

    if not zc1:
        raise ValueError("Signal1 did not contain any zero-crossings; cannot compute trigger offset.")

    # Use first crossing in signal1 as the trigger.
    trigger_index = zc1[0]
    trigger_time = trigger_index / sample_rate

    # Compute zero-crossing times relative to the trigger for signal1.
    times1 = [(x / sample_rate) - trigger_time for x in zc1]
    # For signal2, consider only those crossings that occur at or after the trigger.
    times2 = [(x / sample_rate) - trigger_time for x in zc2 if x >= trigger_index]

    # Pair up the crossings (using the minimum available count):
    num_pairs = min(len(times1), len(times2))
    corrected_crossing_differences = []
    for i in range(num_pairs):
        diff = times1[i] - times2[i]
        corrected_crossing_differences.append(diff)

    # The average delta time is given by:
    average_delta_time = np.mean(corrected_crossing_differences)

    # --- FFT-based Analysis ---
    n = len(signal1)
    fft_signal1 = np.fft.fft(signal1)
    fft_signal2 = np.fft.fft(signal2)
    frequencies = np.fft.fftfreq(n, 1 / sample_rate)
    valid_indices = np.where(frequencies > 0)[0]

    amplitude_db_signal1 = analysis_factor + 20 * np.log10((np.abs(fft_signal1) / n) + 1e-10)
    amplitude_db_signal2 = analysis_factor + 20 * np.log10((np.abs(fft_signal2) / n) + 1e-10)

    max_index_signal1 = valid_indices[np.argmax(np.abs(fft_signal1[valid_indices]))]
    max_index_signal2 = valid_indices[np.argmax(np.abs(fft_signal2[valid_indices]))]

    # Compute phases and then normalize them.
    phase_signal1 = normalize_phase(np.angle(fft_signal1[max_index_signal1], deg=True))
    phase_signal2 = normalize_phase(np.angle(fft_signal2[max_index_signal2], deg=True))
    # The raw phase difference might also be outside the desired range; normalize it.
    phase_difference = normalize_phase(phase_signal2 - phase_signal1)

    # Max amplitude and frequency
    max_freq_signal1 = frequencies[max_index_signal1]
    max_freq_signal2 = frequencies[max_index_signal2]
    max_amplitude_db_signal1 = amplitude_db_signal1[max_index_signal1]
    max_amplitude_db_signal2 = amplitude_db_signal2[max_index_signal2]

    # Angular frequencies
    angular_frequency_signal1 = 2 * np.pi * max_freq_signal1
    angular_frequency_signal2 = 2 * np.pi * max_freq_signal2

    # Gain
    gain_db = max_amplitude_db_signal2 - max_amplitude_db_signal1

    # Results dictionary
    results = {
        "fft_frequencies": frequencies,
        "amplitude_db_signal1": amplitude_db_signal1,
        "amplitude_db_signal2": amplitude_db_signal2,
        "max_frequency_signal1": max_freq_signal1,
        "max_amplitude_db_signal1": max_amplitude_db_signal1,
        "max_frequency_signal2": max_freq_signal2,
        "max_amplitude_db_signal2": max_amplitude_db_signal2,
        "angular_frequency_signal1": angular_frequency_signal1,
        "angular_frequency_signal2": angular_frequency_signal2,
        "phase_signal1_degrees": phase_signal1,
        "phase_signal2_degrees": phase_signal2,
        "phase_difference_degrees": phase_difference,
        "average_delta_time": average_delta_time,  # Corrected average zero-crossing time difference
        "gain_db": gain_db,
        "fft_real_signal1": np.real(fft_signal1),
        "fft_imag_signal1": np.imag(fft_signal1),
        "fft_real_signal2": np.real(fft_signal2),
        "fft_imag_signal2": np.imag(fft_signal2),
    }
    return results

def filter_within_std(data):
    """
    Filters an array to include only values within one standard deviation of the mean.

    Args:
        data (list or numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: A new array containing only the values within one standard deviation of the mean.
    """
    data_np = np.array(data)
    mean = np.mean(data_np)
    std_dev = np.std(data_np)
    filtered_data = data_np[(data_np >= mean - std_dev) & (data_np <= mean + std_dev)]
    return filtered_data

 
def replace_below_threshold(data_list, threshold, replacement_value):
    """
    Converts a list to a NumPy array and replaces values below a threshold.

    Args:
        data_list (list or numpy.ndarray): The input list or array.
        threshold (float): The threshold value.
        replacement_value (float): The value to replace elements below the threshold.

    Returns:
        numpy.ndarray: A NumPy array with values below the threshold replaced.
    """
    data_array = np.array(data_list)
    data_array[data_array < threshold] = replacement_value
    return data_array

def replace_above_threshold(data_list, threshold, replacement_value):
    """
    Converts a list to a NumPy array and replaces values above a threshold.

    Args:
        data_list (list or numpy.ndarray): The input list or array.
        threshold (float): The threshold value.
        replacement_value (float): The value to replace elements above the threshold.

    Returns:
        numpy.ndarray: A NumPy array with values above the threshold replaced.
    """
    data_array = np.array(data_list)
    data_array[data_array > threshold] = replacement_value
    return data_array
    
def graph_processing(OSCvRange, FFTxMax, FFTmaxVal, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout):
    """
    Processes FFT and oscilloscope data, and updates GUI plots accordingly.

    Args:
        FFTxMax: Maximum X value of the FFT.
        FFTmaxVal: Maximum value of the FFT.
        FFTx: X values of the FFT.
        FFTy: Y values of the FFT.
        OSCxin: Input X values from the oscilloscope.
        OSCyin: Input Y values from the oscilloscope.
        OSCxout: Output X values from the oscilloscope.
        OSCyout: Output Y values from the oscilloscope.
    """

    MxFFTx.append(FFTxMax)
    MxFFTy.append(FFTmaxVal)
    
    if len(MxFFTx) == len(MxFFTy) and len(MxFFTx) > 0:
        nsamples = len(MxFFTy)
        dpg.set_value('MxFTT_SCATTEROUT', [list(MxFFTx[-nsamples:]), list(MxFFTy[-nsamples:])])
    else:
        logging.info("FTT Plot data issue seen MaxFFTG V at ", len(MxFFTy))

    threshold = -120
    new_value = -120
    FFTythreshold = replace_below_threshold(FFTy, threshold, new_value)
    
    if len(FFTx) == len(FFTythreshold) and len(FFTx) > 0:
        nsamples = len(FFTx)
        dpg.set_value('FTT_SCATTEROUT', [list(FFTx[-nsamples:]), list(FFTythreshold[-nsamples:])])
        dpg.configure_item(item='DRAG_LINE_FFT_X', default_value=np.median(FFTx))
        dpg.configure_item(item='DRAG_LINE_FFT_X', default_value=np.median(FFTx))
        if len(MxFFTy) > 0:
            dpg.configure_item(item='DRAG_LINE_FFT_Y', default_value=np.median(filter_within_std(MxFFTy)))
        dpg.set_axis_limits_auto(axis='FFT_Y')
        dpg.fit_axis_data(axis='FFT_Y')
    else:
        logging.info("FTT Plot data issue seen at ", len(FFTx))
    # Plotting the oscilloscope data
    if len(OSCyout) == len(OSCxout) and len(OSCyout) > 0:
        nsamples = len(OSCyout)
        dpg.set_value('OSC_SCATTEROUT', [list(OSCxout[-nsamples:]), list(OSCyout[-nsamples:])])
    
    if len(OSCyin) == len(OSCxin) and len(OSCyin) > 0:
        nsamples = len(OSCyin)
        dpg.set_value('OSC_SCATTERIN', [list(OSCxin[-nsamples:]), list(OSCyin[-nsamples:])])
        dpg.configure_item(item='DRAG_LINE_OSC_X', default_value=np.median(OSCxin))
        if len(OSCyin) > 0:
            dpg.configure_item(item='DRAG_LINE_OSC_Y', default_value=np.median(OSCyin))
        RangeOffset = 0.75
        dpg.set_axis_limits('OSC_Y', -RangeOffset * OSCvRange, RangeOffset * OSCvRange)
        dpg.set_axis_limits_auto(axis='OSC_X')
        dpg.fit_axis_data(axis='OSC_X')
    
    if len(OSCyout) == len(OSCyin) and len(OSCyin) > 0:
        dpg.set_value('XY_SCATTEROUT', [list(OSCyin[-nsamples:]), list(OSCyout[-nsamples:])])
        dpg.set_axis_limits('XY_Y', -RangeOffset * OSCvRange, RangeOffset * OSCvRange)
        dpg.set_axis_limits('XY_X', -RangeOffset * OSCvRange, RangeOffset * OSCvRange)
    else:
        logging.info("XY Plot data issue seen at " + str(len(FFTx)))

def graph_processing_Bode(gX, gY, pX, pY):
    """
    Processes Bode plot data and updates GUI plots accordingly.

    Args:
        gX: X values for gain.
        gY: Y values for gain.
        pX: X values for phase.
        pY: Y values for phase.
    """
    
    gain_Y.append(gY)
    gain_X.append(gX)

    # pY = replace_above_threshold(pYb, threshold, new_value)
    phase_Y.append(pY)
    phase_X.append(pX)
    
    # Plotting the magnitude
    if len(gain_X) == len(gain_Y) and len(gain_X) > 0:
        nsamples = len(gain_X)
        dpg.set_value('MAG_SCATTER', [list(gain_X[-nsamples:]), list(gain_Y[-nsamples:])])
        dpg.configure_item(item='DRAG_LINE_MAG_X', default_value=np.median(gain_X))
        dpg.configure_item(item='DRAG_LINE_MAG_Y', default_value=np.median(gain_Y))
        dpg.set_axis_limits_auto(axis='MAG_Y')
        dpg.fit_axis_data(axis='MAG_Y')
    
    # Plotting the phase
    if len(phase_X) == len(phase_Y) and len(phase_X) > 0:
        nsamples = len(phase_X)
        dpg.set_value('PHASE_SCATTER', [list(phase_X[-nsamples:]), list(phase_Y[-nsamples:])])
        dpg.configure_item(item='DRAG_LINE_PHASE_X', default_value=np.median(phase_X))
        dpg.configure_item(item='DRAG_LINE_PHASE_Y', default_value=np.median(phase_Y))

def post_processing(gain_X, gain_Y, phase_X, phase_Y):
    """
    Performs post-processing of the plots and updates GUI elements accordingly.
    """
    # Ensure gain_X is sorted and strictly increasing
    sorted_indices = np.argsort(gain_X)
    gain_X = np.array(gain_X)[sorted_indices]
    gain_Y = np.array(gain_Y)[sorted_indices]

    # Remove duplicates in gain_X
    unique_indices = np.diff(gain_X) > 0
    gain_X = gain_X[np.concatenate(([True], unique_indices))]
    gain_Y = gain_Y[np.concatenate(([True], unique_indices))]

    # Plotting the magnitude
    if len(gain_X) >= 4:  # Ensure enough points for spline
        magnitude_spline = UnivariateSpline(gain_X, gain_Y, k=3, s=0)
        try:
            dpg.delete_item('MAG_LINE')  # Ensure the tag is unique
        except Exception:
            pass
        dpg.add_line_series(x=gain_X, y=magnitude_spline(gain_X), parent='MAG_Y', tag='MAG_LINE', label='Interpolated waveform')
        dpg.configure_item(item='DRAG_LINE_MAG_X', default_value=np.median(gain_X))
        dpg.configure_item(item='DRAG_LINE_MAG_Y', default_value=np.median(gain_Y))
        dpg.set_axis_limits_auto(axis='MAG_Y')
        dpg.fit_axis_data(axis='MAG_Y')
        # dpg.set_axis_limits_auto(axis='MAG_X')
        # dpg.fit_axis_data(axis='MAG_X')
    else:
        logging.warning("Insufficient points for Magnitude Spline. Skipping fit.")
        
    # Ensure phase_X is sorted and strictly increasing
    sorted_indices = np.argsort(phase_X)
    phase_X = np.array(phase_X)[sorted_indices]
    phase_Y = np.array(phase_Y)[sorted_indices]

    # Remove duplicates in phase_X
    unique_indices = np.diff(phase_X) > 0
    phase_X = phase_X[np.concatenate(([True], unique_indices))]
    phase_Y = phase_Y[np.concatenate(([True], unique_indices))]

    # Outlier detection based on local deviation
    window_size = 12  # Number of neighbors to consider on each side
    threshold = 10   # Threshold in degrees

    filtered_indices = []
    for i in range(len(phase_Y)):
        # Determine the neighborhood indices
        start = max(0, i - window_size)
        end = min(len(phase_Y), i + window_size + 1)
        neighbors = np.concatenate([phase_Y[start:i], phase_Y[i + 1:end]])

        if len(neighbors) > 0:
            local_median = np.median(neighbors)
            if abs(phase_Y[i] - local_median) <= threshold:
                filtered_indices.append(i)

    # Filter phase_X and phase_Y
    filtered_phase_X = phase_X[filtered_indices]
    filtered_phase_Y = phase_Y[filtered_indices]

    # Plotting the filtered phase
    if len(filtered_phase_X) >= 4:  # Ensure enough points for spline
        phase_spline = UnivariateSpline(filtered_phase_X, filtered_phase_Y, k=3, s=0)
        try:
            dpg.delete_item('PHASE_LINE')  # Ensure the tag is unique
        except Exception:
            pass
        dpg.add_line_series(x=filtered_phase_X, y=phase_spline(filtered_phase_X), parent='PHASE_Y', tag='PHASE_LINE', label='Interpolated waveform')
        dpg.configure_item(item='DRAG_LINE_PHASE_X', default_value=np.median(filtered_phase_X))
        dpg.configure_item(item='DRAG_LINE_PHASE_Y', default_value=np.median(filtered_phase_Y))
        # dpg.set_axis_limits_auto(axis='PHASE_X')
        # dpg.fit_axis_data(axis='PHASE_X')
    else:
        logging.warning("Insufficient points for Phase Spline. Skipping fit.")

def set_drag_lines():
    """
    Sets the visibility of drag lines in the GUI based on the selected option.
    """
    
    if dpg.get_value(item='DRAG_LINES_COMBO') == drag_line_values[1]:
        dpg.configure_item(item='DRAG_LINE_MAG_X', show=False)
        dpg.configure_item(item='DRAG_LINE_MAG_Y', show=False)
        dpg.configure_item(item='DRAG_LINE_PHASE_X', show=False)
        dpg.configure_item(item='DRAG_LINE_PHASE_Y', show=False)
    else:
        dpg.configure_item(item='DRAG_LINE_MAG_X', show=True)
        dpg.configure_item(item='DRAG_LINE_MAG_Y', show=True)
        dpg.configure_item(item='DRAG_LINE_PHASE_X', show=True)
        dpg.configure_item(item='DRAG_LINE_PHASE_Y', show=True)

def callbackCSV(sender, app_data, user_data):
    global LogFile
    """
    Handles the callbackCSV for file selection.

    Args:
        sender: The sender of the event.
        app_data: The data associated with the event, containing file path information.
        user_data: Additional user data passed to the callbackCSV.
    """
    LogFile = app_data.get("file_path_name")
    LogFileGUIlabel = "CSV Log File: " + LogFile
    dpg.configure_item('LOGFILE', label=LogFileGUIlabel)

def callbackJSON(sender, app_data, user_data):
    global JSONLogFile
    """
    Handles the callbackJSON for file selection.

    Args:
        sender: The sender of the event.
        app_data: The data associated with the event, containing file path information.
        user_data: Additional user data passed to the callbackJSON.
    """
    JSONLogFile = app_data.get("file_path_name")
    LogFileGUIlabel = "JSON Log File: " + JSONLogFile
    dpg.configure_item('JSONLOGFILEtag', label=LogFileGUIlabel)
    dpg.configure_item(item='Play', enabled=True)
    dpg.bind_item_theme('Play', 'GreenButton')
    dpg.configure_item(item='PlaySpeed', enabled=True)
    
def create_filename_from_timestamp(prefix="file", extension="txt"):
    """
    Creates a filename using a timestamp.

    Args:
        prefix (str): The prefix for the filename.
        extension (str): The file extension.

    Returns:
        str: The generated filename with the timestamp.
    """
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{prefix}_{timestamp}.{extension}"
    return filename


def stop_exit():
    # Close all the processes and exit
    os._exit(0)
    
# Function to switch tabs programmatically
def select_tab(tab_name):
    dpg.set_value('tab_bar', tab_name)

# Function to create themes
def create_themes():
    global dark_theme, light_theme
    
    # -- Theme Settings for GUI ---
    with dpg.theme(tag='Dark') as dark_theme:
        with dpg.theme_component(0):
            # Set theme colors for various elements
            # sourced from https://raw.githubusercontent.com/hoffstadt/DearPyGui_Ext/refs/heads/master/dearpygui_ext/themes.py
            dpg.add_theme_color(dpg.mvThemeCol_Text                   , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TextDisabled           , (0.50 * 255, 0.50 * 255, 0.50 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg               , (0.06 * 255, 0.06 * 255, 0.06 * 255, 0.94 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg                , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PopupBg                , (0.08 * 255, 0.08 * 255, 0.08 * 255, 0.94 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Border                 , (0.43 * 255, 0.43 * 255, 0.50 * 255, 0.50 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_BorderShadow           , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg                , (0.16 * 255, 0.29 * 255, 0.48 * 255, 0.54 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered         , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.40 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive          , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.67 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBg                , (0.04 * 255, 0.04 * 255, 0.04 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive          , (0.16 * 255, 0.29 * 255, 0.48 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgCollapsed       , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.51 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg              , (0.14 * 255, 0.14 * 255, 0.14 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarBg            , (0.02 * 255, 0.02 * 255, 0.02 * 255, 0.53 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrab          , (0.31 * 255, 0.31 * 255, 0.31 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabHovered   , (0.41 * 255, 0.41 * 255, 0.41 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabActive    , (0.51 * 255, 0.51 * 255, 0.51 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_CheckMark              , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrab             , (0.24 * 255, 0.52 * 255, 0.88 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrabActive       , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Button                 , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.40 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered          , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive           , (0.06 * 255, 0.53 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Header                 , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.31 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_HeaderHovered          , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.80 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_HeaderActive           , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Separator              , (0.43 * 255, 0.43 * 255, 0.50 * 255, 0.50 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorHovered       , (0.10 * 255, 0.40 * 255, 0.75 * 255, 0.78 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorActive        , (0.10 * 255, 0.40 * 255, 0.75 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ResizeGrip             , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.20 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ResizeGripHovered      , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.67 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ResizeGripActive       , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.95 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Tab                    , (0.18 * 255, 0.35 * 255, 0.58 * 255, 0.86 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabHovered             , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.80 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabActive              , (0.20 * 255, 0.41 * 255, 0.68 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocused           , (0.07 * 255, 0.10 * 255, 0.15 * 255, 0.97 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocusedActive     , (0.14 * 255, 0.26 * 255, 0.42 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_DockingPreview         , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.70 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_DockingEmptyBg         , (0.20 * 255, 0.20 * 255, 0.20 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotLines              , (0.61 * 255, 0.61 * 255, 0.61 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotLinesHovered       , (1.00 * 255, 0.43 * 255, 0.35 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotHistogram          , (0.90 * 255, 0.70 * 255, 0.00 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotHistogramHovered   , (1.00 * 255, 0.60 * 255, 0.00 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableHeaderBg          , (0.19 * 255, 0.19 * 255, 0.20 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableBorderStrong      , (0.31 * 255, 0.31 * 255, 0.35 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableBorderLight       , (0.23 * 255, 0.23 * 255, 0.25 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableRowBg             , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableRowBgAlt          , (1.00 * 255, 1.00 * 255, 1.00 * 255, 0.06 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TextSelectedBg         , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.35 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_DragDropTarget         , (1.00 * 255, 1.00 * 255, 0.00 * 255, 0.90 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_NavHighlight           , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_NavWindowingHighlight  , (1.00 * 255, 1.00 * 255, 1.00 * 255, 0.70 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_NavWindowingDimBg      , (0.80 * 255, 0.80 * 255, 0.80 * 255, 0.20 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ModalWindowDimBg       , (0.80 * 255, 0.80 * 255, 0.80 * 255, 0.35 * 255))
            dpg.add_theme_color(dpg.mvPlotCol_FrameBg                 , (1.00 * 255, 1.00 * 255, 1.00 * 255, 0.07 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_PlotBg                  , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.50 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_PlotBorder              , (0.43 * 255, 0.43 * 255, 0.50 * 255, 0.50 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_LegendBg                , (0.08 * 255, 0.08 * 255, 0.08 * 255, 0.94 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_LegendBorder            , (0.43 * 255, 0.43 * 255, 0.50 * 255, 0.50 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_LegendText              , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_TitleText               , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_InlayText               , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisBg                  , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisBgActive            , (0.06 * 255, 0.53 * 255, 0.98 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisBgHovered           , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisGrid                , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisText                , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_Selection               , (1.00 * 255, 0.60 * 255, 0.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_Crosshairs              , (1.00 * 255, 1.00 * 255, 1.00 * 255, 0.50 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvNodeCol_NodeBackground          , (50, 50, 50, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_NodeBackgroundHovered   , (75, 75, 75, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_NodeBackgroundSelected  , (75, 75, 75, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_NodeOutline             , (100, 100, 100, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_TitleBar                , (41, 74, 122, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_TitleBarHovered         , (66, 150, 250, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_TitleBarSelected        , (66, 150, 250, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_Link                    , (61, 133, 224, 200), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_LinkHovered             , (66, 150, 250, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_LinkSelected            , (66, 150, 250, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_Pin                     , (53, 150, 250, 180), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_PinHovered              , (53, 150, 250, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_BoxSelector             , (61, 133, 224, 30), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_BoxSelectorOutline      , (61, 133, 224, 150), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_GridBackground          , (40, 40, 50, 200), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_GridLine                , (200, 200, 200, 40), category=dpg.mvThemeCat_Nodes)
 
            
        # Theme for line series
        with dpg.theme_component(dpg.mvLineSeries):
            dpg.add_theme_color(dpg.mvPlotCol_Line, (51, 51, 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_style(dpg.mvPlotStyleVar_LineWeight, 3, category=dpg.mvThemeCat_Plots)

        # Theme for scatter series
        with dpg.theme_component(dpg.mvScatterSeries):
            dpg.add_theme_color(dpg.mvPlotCol_Line, (255, 0, 0), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_style(dpg.mvPlotStyleVar_Marker, dpg.mvPlotMarker_Circle, category=dpg.mvThemeCat_Plots)
            dpg.add_theme_style(dpg.mvPlotStyleVar_MarkerSize, 3, category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_MarkerOutline, (255, 0, 0), category=dpg.mvThemeCat_Plots)

        # Create a theme for the plot (green and blue)
        with dpg.theme(tag="plot_theme_green_blue"):
            with dpg.theme_component(dpg.mvStemSeries):
                dpg.add_theme_color(dpg.mvPlotCol_Line, (150, 255, 0), category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_Marker, dpg.mvPlotMarker_Diamond, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_MarkerSize, 3, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_color(dpg.mvPlotCol_MarkerOutline, (150, 255, 0), category=dpg.mvThemeCat_Plots)

            with dpg.theme_component(dpg.mvScatterSeries):
                dpg.add_theme_color(dpg.mvPlotCol_Line, (60, 150, 200), category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_Marker, dpg.mvPlotMarker_Square, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_MarkerSize, 3, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_color(dpg.mvPlotCol_MarkerOutline, (60, 150, 200), category=dpg.mvThemeCat_Plots)

        # Create a theme for the plot (orange and yellow)
        with dpg.theme(tag="plot_theme_org_yellow"):
            with dpg.theme_component(dpg.mvStemSeries):
                dpg.add_theme_color(dpg.mvPlotCol_Line, (250, 128, 0), category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_Marker, dpg.mvPlotMarker_Diamond, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_MarkerSize, 3, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_color(dpg.mvPlotCol_MarkerOutline, (250, 128, 0), category=dpg.mvThemeCat_Plots)

            with dpg.theme_component(dpg.mvScatterSeries):
                dpg.add_theme_color(dpg.mvPlotCol_Line, (242, 255, 46), category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_Marker, dpg.mvPlotMarker_Square, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_style(dpg.mvPlotStyleVar_MarkerSize, 3, category=dpg.mvThemeCat_Plots)
                dpg.add_theme_color(dpg.mvPlotCol_MarkerOutline, (242, 255, 46), category=dpg.mvThemeCat_Plots)
                
        # Green button theme
        with dpg.theme(tag='GreenButton') as green_button_theme:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (0, 255, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (0, 200, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (0, 150, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (0, 0, 0, 255))

        # Red button theme
        with dpg.theme(tag='RedButton') as red_button_theme:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (255, 0, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (200, 0, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (150, 0, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (255, 255, 255, 255))

        # Yellow button theme
        with dpg.theme(tag='YellowButton') as yellow_button_theme:
            with dpg.theme_component(dpg.mvButton):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (255, 255, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (200, 200, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (150, 150, 0, 255))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (0, 0, 0, 255))

        # Disabled button theme
        with dpg.theme(tag='DisabledButton') as disabled_button_theme:
            with dpg.theme_component(dpg.mvButton, enabled_state=False):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (128, 128, 128, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (128, 128, 128, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (128, 128, 128, 255))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (255, 255, 255, 255))
        
        # Create a theme and add a specific style for the scrollbar size
        with dpg.theme(tag='WinScrollbar') as window_theme:
            with dpg.theme_component(dpg.mvAll):
                dpg.add_theme_style(dpg.mvStyleVar_ScrollbarSize, scrollbar_width)
        
        # Create a theme for the plot with custom padding values.
        with dpg.theme(tag="PlotPadding") as plot_theme:
            with dpg.theme_component(dpg.mvPlot):
                # For a vector (x, y) value, pass each component separately.
                dpg.add_theme_style(dpg.mvPlotStyleVar_PlotPadding, 20, 0, category=dpg.mvThemeCat_Plots)
       
    # Define a theme with the tag 'Light'
    with dpg.theme(tag='Light') as light_theme:
        # Apply the theme component with ID 0
        with dpg.theme_component(0):
            # Set theme colors for various elements
            # sourced from https://raw.githubusercontent.com/hoffstadt/DearPyGui_Ext/refs/heads/master/dearpygui_ext/themes.py
            dpg.add_theme_color(dpg.mvThemeCol_Text                   , (0.00 * 255, 0.00 * 255, 0.00 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TextDisabled           , (0.60 * 255, 0.60 * 255, 0.60 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg               , (0.94 * 255, 0.94 * 255, 0.94 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg                , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PopupBg                , (1.00 * 255, 1.00 * 255, 1.00 * 255, 0.98 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Border                 , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.30 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_BorderShadow           , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg                , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered         , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.40 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive          , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.67 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBg                , (0.96 * 255, 0.96 * 255, 0.96 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive          , (0.82 * 255, 0.82 * 255, 0.82 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgCollapsed       , (1.00 * 255, 1.00 * 255, 1.00 * 255, 0.51 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg              , (0.86 * 255, 0.86 * 255, 0.86 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarBg            , (0.98 * 255, 0.98 * 255, 0.98 * 255, 0.53 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrab          , (0.69 * 255, 0.69 * 255, 0.69 * 255, 0.80 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabHovered   , (0.49 * 255, 0.49 * 255, 0.49 * 255, 0.80 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabActive    , (0.49 * 255, 0.49 * 255, 0.49 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_CheckMark              , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrab             , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.78 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrabActive       , (0.46 * 255, 0.54 * 255, 0.80 * 255, 0.60 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Button                 , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.40 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered          , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive           , (0.06 * 255, 0.53 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Header                 , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.31 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_HeaderHovered          , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.80 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_HeaderActive           , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Separator              , (0.39 * 255, 0.39 * 255, 0.39 * 255, 0.62 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorHovered       , (0.14 * 255, 0.44 * 255, 0.80 * 255, 0.78 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorActive        , (0.14 * 255, 0.44 * 255, 0.80 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ResizeGrip             , (0.35 * 255, 0.35 * 255, 0.35 * 255, 0.17 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ResizeGripHovered      , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.67 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ResizeGripActive       , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.95 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_Tab                    , (0.76 * 255, 0.80 * 255, 0.84 * 255, 0.93 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabHovered             , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.80 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabActive              , (0.60 * 255, 0.73 * 255, 0.88 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocused           , (0.92 * 255, 0.93 * 255, 0.94 * 255, 0.99 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocusedActive     , (0.74 * 255, 0.82 * 255, 0.91 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_DockingPreview         , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.22 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_DockingEmptyBg         , (0.20 * 255, 0.20 * 255, 0.20 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotLines              , (0.39 * 255, 0.39 * 255, 0.39 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotLinesHovered       , (1.00 * 255, 0.43 * 255, 0.35 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotHistogram          , (0.90 * 255, 0.70 * 255, 0.00 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_PlotHistogramHovered   , (1.00 * 255, 0.45 * 255, 0.00 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableHeaderBg          , (0.78 * 255, 0.87 * 255, 0.98 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableBorderStrong      , (0.57 * 255, 0.57 * 255, 0.64 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableBorderLight       , (0.68 * 255, 0.68 * 255, 0.74 * 255, 1.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableRowBg             , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TableRowBgAlt          , (0.30 * 255, 0.30 * 255, 0.30 * 255, 0.09 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_TextSelectedBg         , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.35 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_DragDropTarget         , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.95 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_NavHighlight           , (0.26 * 255, 0.59 * 255, 0.98 * 255, 0.80 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_NavWindowingHighlight  , (0.70 * 255, 0.70 * 255, 0.70 * 255, 0.70 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_NavWindowingDimBg      , (0.20 * 255, 0.20 * 255, 0.20 * 255, 0.20 * 255))
            dpg.add_theme_color(dpg.mvThemeCol_ModalWindowDimBg       , (0.20 * 255, 0.20 * 255, 0.20 * 255, 0.35 * 255))
            dpg.add_theme_color(dpg.mvPlotCol_FrameBg                 , (1.00 * 255, 1.00 * 255, 1.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_PlotBg                  , (0.42 * 255, 0.57 * 255, 1.00 * 255, 0.13 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_PlotBorder              , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_LegendBg                , (1.00 * 255, 1.00 * 255, 1.00 * 255, 0.98 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_LegendBorder            , (0.82 * 255, 0.82 * 255, 0.82 * 255, 0.80 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_LegendText              , (0.00 * 255, 0.00 * 255, 0.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_TitleText               , (0.00 * 255, 0.00 * 255, 0.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_InlayText               , (0.00 * 255, 0.00 * 255, 0.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisBg                  , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisBgActive            , (0.06 * 255, 0.53 * 255, 0.98 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisBgHovered           , (0.26 * 255, 0.59 * 255, 0.98 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisGrid                , (0.00 * 255, 0.00 * 255, 0.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_AxisText                , (0.00 * 255, 0.00 * 255, 0.00 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_Selection               , (0.82 * 255, 0.64 * 255, 0.03 * 255, 1.00 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_Crosshairs              , (0.00 * 255, 0.00 * 255, 0.00 * 255, 0.50 * 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvNodeCol_NodeBackground          , (240, 240, 240, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_NodeBackgroundHovered   , (240, 240, 240, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_NodeBackgroundSelected  , (240, 240, 240, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_NodeOutline             , (100, 100, 100, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_TitleBar                , (248, 248, 248, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_TitleBarHovered         , (209, 209, 209, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_TitleBarSelected        , (209, 209, 209, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_Link                    , (66, 150, 250, 100), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_LinkHovered             , (66, 150, 250, 242), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_LinkSelected            , (66, 150, 250, 242), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_Pin                     , (66, 150, 250, 160), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_PinHovered              , (66, 150, 250, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_BoxSelector             , (90, 170, 250, 30), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_BoxSelectorOutline      , (90, 170, 250, 150), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_GridBackground          , (225, 225, 225, 255), category=dpg.mvThemeCat_Nodes)
            dpg.add_theme_color(dpg.mvNodeCol_GridLine                , (180, 180, 180, 100), category=dpg.mvThemeCat_Nodes)
        
        # Theme for line series
        with dpg.theme_component(dpg.mvLineSeries):
            dpg.add_theme_color(dpg.mvPlotCol_Line, (51, 51, 255), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_style(dpg.mvPlotStyleVar_LineWeight, 3, category=dpg.mvThemeCat_Plots)

        # Theme for scatter series
        with dpg.theme_component(dpg.mvScatterSeries):
            dpg.add_theme_color(dpg.mvPlotCol_Line, (255, 0, 0), category=dpg.mvThemeCat_Plots)
            dpg.add_theme_style(dpg.mvPlotStyleVar_Marker, dpg.mvPlotMarker_Circle, category=dpg.mvThemeCat_Plots)
            dpg.add_theme_style(dpg.mvPlotStyleVar_MarkerSize, 3, category=dpg.mvThemeCat_Plots)
            dpg.add_theme_color(dpg.mvPlotCol_MarkerOutline, (255, 0, 0), category=dpg.mvThemeCat_Plots)

        # Theme for a disabled button
        with dpg.theme_component(dpg.mvButton, enabled_state=False):
            dpg.add_theme_color(dpg.mvThemeCol_Text, [150, 150, 150])


# Function to switch theme based on user selection
def switch_theme(sender, app_data):
    selected_theme = dpg.get_value('WIN_THEME')
    if selected_theme == 'Dark':
        dpg.bind_theme(dark_theme)
    else:
        dpg.bind_theme(light_theme)

def on_viewport_close():
    logging.info("Dear PyGUI viewport closed.")
    stop_event.set()  # Signal to stop any running processes or threads

# Always initialize stop_event
stop_event = multiprocessing.Event()

# Callback function that updates CH_OUT based on the selected input channel
def update_ch_out(sender, app_data):
    # app_data contains the new selection from CH_IN
    if app_data == available_channels[0]:
        new_channel = available_channels[1]
    else:
        new_channel = available_channels[0]
    # Update the text widget using its tag
    dpg.set_value("CH_OUT", new_channel)

def checkbox_callback(sender, app_data, user_data):
    """Update chart visibility flags based on checkbox interactions."""
    index = user_data["index"]
    chart_visibility[index] = app_data
    print(f"{chart_names[index]} visibility set to {app_data}")

def radio_button_callback(sender, app_data):
    """
    Toggle visibility of the checkbox group and its placeholder based on the plot engine selection.
    When MatPlotLib is selected, show the checkboxes.
    Otherwise, hide them and show the spacing placeholder.
    """
    if app_data == "MatPlotLib":
        dpg.show_item("checkbox_group")
        dpg.hide_item("spacing_placeholder")
    else:
        dpg.hide_item("checkbox_group")
        dpg.show_item("spacing_placeholder")


# -- gui settings --- ---------------------------------------------
def main():
    global channel_in, channel_out, CH1_probe_attenuation_ratio, CH2_probe_attenuation_ratio
    global CH1_coupling, CH2_coupling, Sample_command, DEPMEM, waveform_amplitude_V
    global AWG_output_impedance, points_per_decade, start_decade, stop_decade, addpm
    global point_resize_factor, vertical_scaling_factor, horizontal_scaling_factor
    global nWaveSamples, FTTcorection, read_delay_ms, sample_delay_s, plot_win_disposition
    global is_playing, is_paused, play_speed, is_recording, working_dir
    global CH1_amplitude_scales, CH2_amplitude_scales, MaxFreqAdjust


    # --- Set working path --- ------------------------------------------------------------
    try:
        documents_path = os.path.join(os.path.expanduser("~"), "Documents")
        if os.path.exists(documents_path):
            os.chdir(documents_path)
        else:
            logging.error("Documents directory not found. Defaulting to home directory.")
            os.chdir(os.path.expanduser("~"))
    except Exception as e:
        logging.error(f"Failed to set working directory: {e}")
        # Fallback to a default safe directory (e.g., home directory)
        # os.chdir(os.path.expanduser("~"))
    # --- End Set working path --- --------------------------------------------------------
    
    logging.info("Main application running normally")
    # Normal application logic goes here
    
    # Register the signal handler in the main thread
    signal.signal(signal.SIGINT, signal_handler)
    
    dpg.create_context()

    # Create themes
    create_themes()
    
    # Create a file dialog for selecting files
    with dpg.file_dialog(directory_selector=False, show=False, callback=callbackCSV, id="file_dialog_id1", width=700, height=400):
        dpg.add_file_extension(".csv", color=(106, 127, 235, 255), custom_text="[csv]")
    with dpg.file_dialog(directory_selector=False, show=False, callback=callbackJSON, id="file_dialog_id2", width=700, height=400):
        dpg.add_file_extension(".json", color=(106, 127, 235, 255), custom_text="[json]")        

    with dpg.window(tag='main_window', pos=(0, 0), width=dpgWindow1, no_bring_to_front_on_focus=True, autosize=True, height=dpgWindow1, no_title_bar=True, no_move=True, no_resize=False):
        dpg.bind_item_theme('main_window', 'WinScrollbar') 
        with dpg.child_window(tag='scrolling_container', width=dpgWindow2, height=main_window_height, horizontal_scrollbar=False, border=False):
            with dpg.child_window(tag='controlwin', width=dpgWindow3, height=setting_window_height - win_vertical_border, pos=(0, win_vertical_border), border=False, menubar=False): 
                dpg.add_text("Oscilloscope settings:")
                # Add various combo boxes for different settings
                with dpg.group(horizontal=True):
                    dpg.add_combo(tag='CH_IN', items=available_channels, label='Input Channel', default_value=available_channels[0], width=items_standard_width, callback=update_ch_out)
                    # Add the output channel input field, locked for editing via readonly=True.
                    dpg.add_input_text(tag="CH_OUT", label="Output Channel", default_value=available_channels[1], width=items_standard_width, readonly=True)
                dpg.add_combo(tag='CH1_ATTENUATION_RATIO', items=probes_attenuation_ratio, label='Channel 1 Probe attenuation ratio', default_value=probes_attenuation_ratio[1], width=items_standard_width)
                dpg.add_combo(tag='CH2_ATTENUATION_RATIO', items=probes_attenuation_ratio, label='Channel 2 Probe attenuation ratio', default_value=probes_attenuation_ratio[1], width=items_standard_width)
                dpg.add_combo(tag='CH1_COUPLING_MODE', items=channels_coupling_mode, label='Channel 1 Coupling mode', default_value=channels_coupling_mode[0], width=items_standard_width)
                dpg.add_combo(tag='CH2_COUPLING_MODE', items=channels_coupling_mode, label='Channel 2 Coupling mode', default_value=channels_coupling_mode[0], width=items_standard_width)
                dpg.add_combo(tag='SAMPL_MODE', items=sample_modes, label='Oscilloscope Sample mode', default_value=sample_modes[0], width=items_standard_width)
                dpg.add_combo(tag='DEPMEM', items=memory_depth_modes, label='Oscilloscope Memory depth', default_value=memory_depth_modes[1], width=items_standard_width)
                dpg.add_input_float(tag='AWG_OUT_VOLTAGE', label='AWG pk-pk output voltage', min_value=0, max_value=5, min_clamped=True, max_clamped=True, default_value=1, width=items_standard_width)
                dpg.add_combo(tag='HIGH_Z', items=AWG_output_impedance_modes, label='AWG High output impedance', default_value=AWG_output_impedance_modes[0], width=items_standard_width)
                dpg.add_text('High output impedance = OFF -> Z_out = 50 Ohm\nThe AWG pk-pk output voltage will be doubled!')

                # Add inputs for points and spacing settings
                dpg.add_text('\nPoints/decade   Points spacing')
                with dpg.group(horizontal=True):
                    dpg.add_input_int(tag='POINTS_X_DEC', min_value=0, min_clamped=True, default_value=30, width=items_standard_width)
                    dpg.add_combo(tag='POINTS_SPACING', items=points_spacing, default_value=points_spacing[1], width=120)
                dpg.add_combo(tag='START_DEC', items=decades_list_string, label='Start frequency', default_value=decades_list_string[3], width=items_standard_width)
                dpg.add_combo(tag='STOP_DEC', items=decades_list_string, label='Stop frequency', default_value=decades_list_string[7], width=items_standard_width)

                # Add graphics settings
                dpg.add_text('\nGraphics setting:')
                with dpg.group(horizontal=True):
                    dpg.add_radio_button(tag='WIN_THEME', label='Window theme', items=win_theme, default_value=win_theme[0], callback=switch_theme)
                    # if platform.system() == 'Windows':
                    dpg.add_radio_button(tag='PLOT_WIN_SETTING', label='Plot Engine', items=plot_win_settings, default_value=plot_win_settings[0], callback=radio_button_callback)
                    dpg.add_radio_button(tag='DRAG_LINES_COMBO', items=drag_line_values, default_value=drag_line_values[1], callback=set_drag_lines)


                # Group for chart visibility checkboxes, initially visible since MatPlotLib is default.
                with dpg.group(horizontal=True, tag="checkbox_group"):
                    for i, name in enumerate(chart_names):
                        dpg.add_checkbox(
                            label=name,
                            default_value=chart_visibility[i],
                            callback=checkbox_callback,
                            user_data={"index": i},
                            tag=f"checkbox_{i}"
                        )
                dpg.hide_item("checkbox_group")  # Hide checkboxes initially
                # Add a spacing placeholder to maintain layout spacing
                # (Here we show a text widget with a newline.
                # In practice, dpg.add_dummy(width=200, height=50, tag="spacing_placeholder") # for more precise control over the empty space.
                dpg.add_text("\n", tag="spacing_placeholder")
                 
                # Add advanced settings
                dpg.add_text('\nAdvanced settings:')
                dpg.add_input_float(tag='POINT_SCALE_COEFF', label='Point scale coefficient', min_value=0, min_clamped=True, default_value=5850, width=items_standard_width)
                dpg.add_input_float(tag='V_SCALE_COEFF', label='Vertical scale calibration coeff.', min_value=0, min_clamped=True, default_value=0.33, width=items_standard_width)
                dpg.add_input_float(tag='H_SCALE_COEFF', label='Horizontal scale calibration coeff.', min_value=0, min_clamped=True, default_value=0.80, width=items_standard_width)
                dpg.add_input_int(tag='OSCILL_TIMEOUT', label='Oscilloscope reading timeout (ms)', min_value=1200, min_clamped=True, default_value=1200, width=items_standard_width)
                dpg.add_input_float(tag='CODE_EXEC_PAUSE', label='Commands execution delay (s)', min_value=0.2, min_clamped=True, default_value=0.20, width=items_standard_width)

                # Add buttons for various actions
                dpg.add_input_int(tag='SAMPLE_WAVES', label='Number of times waves are sampled', min_value=1, min_clamped=True, default_value=1, width=items_standard_width)
                dpg.add_input_int(tag='FFToffset', label='FFT correction offset (db)', default_value=2, width=items_standard_width)
                # Acceptable Degree Delta Per Measurement
                dpg.add_input_int(tag='ADDPM', label='Acceptable Degree Delta Per Measurement', min_value=0, max_value=180, min_clamped=True, max_clamped=True, default_value=5, width=items_standard_width)
                dpg.add_input_int(tag='MAXTFREQADJUST', label='Max attempts to Adjust Freq.', min_value=3, max_value=7, min_clamped=True, max_clamped=True, default_value=3, width=items_standard_width)
                dpg.add_text('\n')
                with dpg.group(horizontal=True):
                    dpg.add_button(tag='SEARCH_OSCILLOSCOPE', label='Search and Setup Oscilloscope', callback=search_oscilloscope)
                    dpg.add_button(tag='START_MEASURE', label='Start measurements', callback=start_mesurement, enabled=False)
                    # Create the window with playback buttons and signal plot
                    dpg.add_button(tag='Play', label='Play', callback=play_callback, enabled=False)
                    dpg.add_button(tag='Stop', label='Stop', callback=stop_callback, enabled=False)
                dpg.add_slider_float(tag='PlaySpeed', label="Playback Speed", default_value=1.0, min_value=0.2, max_value=3.0, callback=speed_slider_callback, format="%.1f", enabled=False)
                dpg.add_button(tag='LOGFILE', label="CSV Log File: " + LogFile, callback=lambda: dpg.show_item("file_dialog_id1"))
                dpg.add_button(tag='JSONLOGFILEtag', label="JSON Log File: " + JSONLogFile, callback=lambda: dpg.show_item("file_dialog_id2"))
                dpg.add_text('\n')
                
                # Create tabs
                with dpg.tab_bar(tag='tab_bar'):
                    with dpg.tab(label="Terminal", tag='terminal_tab'):
                        with dpg.child_window(tag='terminal', width=setting_window_width, height=int(main_window_height / 4) - 80, label='Terminal', border=False):
                            pass
                    with dpg.tab(label="Data Table", tag='data_table_tab'):
                        with dpg.child_window(tag='dataTableWindow', width=setting_window_width, height=318, border=False, menubar=False):
                            pass

                            # Add a table with various columns
                            with dpg.table(header_row=True, resizable=True, policy=dpg.mvTable_SizingStretchProp, borders_outerH=True, borders_innerV=True, borders_innerH=True, borders_outerV=True, tag='DataTable'):
                                dpg.add_table_column(label="Frequency")
                                dpg.add_table_column(label="MeasFreq")
                                dpg.add_table_column(label="Voltage")
                                dpg.add_table_column(label="Phase")
                                dpg.add_table_column(label="FFT Hz")
                                dpg.add_table_column(label="FFT Meg")
        
            # Create a window for the magnitude Bode plot
            # with dpg.window(tag='MAG_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, 0), no_close=True, no_collapse=True, no_move=True, no_title_bar=True, no_resize=True):
            with dpg.child_window(tag='MAG_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, 0), border=False, menubar=False):
                with dpg.plot(tag='MAG_PLOT_GRAPH', label='Magnitude Bode Plot', height=plot_window_height, width=plot_window_width, crosshairs=False, anti_aliased=True, no_mouse_pos=False):
                    # Add X-axis with log scale
                    dpg.add_plot_axis(dpg.mvXAxis, label='Frequency (Hz)', scale=dpg.mvPlotScale_Log10, tag='MAG_X')
                    dpg.set_axis_limits('MAG_X', 100, 1000000)
                    # Add Y-axis for measured voltage gain
                    dpg.add_plot_axis(dpg.mvYAxis, label='Measured Voltage Gain (dB)', tag='MAG_Y')
                    # Add draggable lines for frequency and gain
                    dpg.add_drag_line(tag='DRAG_LINE_MAG_X', parent='MAG_PLOT_GRAPH', label="Frequency", color=[255, 255, 0, 255], default_value=0, thickness=2, show=False)
                    dpg.add_drag_line(tag='DRAG_LINE_MAG_Y', parent='MAG_PLOT_GRAPH', label="Gain dB", color=[255, 255, 0, 255], vertical=False, default_value=1, thickness=2, show=False)
                    # Add scatter series for gain data
                    dpg.add_scatter_series(x=gain_X, y=gain_Y, parent='MAG_Y', tag='MAG_SCATTER', label='{} points per decade'.format(points_per_decade))
                    dpg.add_plot_legend(location=2)
                    dpg.bind_item_theme("MAG_PLOT_GRAPH", "PlotPadding")

            # Create a window for the phase Bode plot
            # with dpg.window(tag='PHASE_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, plot_window_height), no_close=True, no_collapse=True, no_move=True, no_title_bar=True, no_resize=True):
            with dpg.child_window(tag='PHASE_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, plot_window_height), border=False, menubar=False):
                with dpg.plot(tag='PHASE_PLOT_GRAPH', label='Phase Bode Plot', height=plot_window_height, width=plot_window_width, crosshairs=False, anti_aliased=True, no_mouse_pos=False):
                    # Add X-axis with log scale
                    dpg.add_plot_axis(dpg.mvXAxis, label='Frequency (Hz)', scale=dpg.mvPlotScale_Log10, tag='PHASE_X')
                    dpg.set_axis_limits('PHASE_X', 100, 1000000)
                    # Add Y-axis for phase shift
                    dpg.add_plot_axis(dpg.mvYAxis, label='Phase shift (deg°)', tag='PHASE_Y')
                    dpg.set_axis_limits('PHASE_Y', -200, 200)
                    # Add draggable lines for frequency and phase
                    dpg.add_drag_line(tag='DRAG_LINE_PHASE_X', parent='PHASE_PLOT_GRAPH', label="Frequency", color=[255, 255, 0, 255], default_value=0, thickness=2, show=False)
                    dpg.add_drag_line(tag='DRAG_LINE_PHASE_Y', parent='PHASE_PLOT_GRAPH', label="Degrees", color=[255, 255, 0, 255], vertical=False, default_value=0, thickness=2, show=False)
                    # Add scatter series for phase data
                    dpg.add_scatter_series(x=phase_X, y=phase_Y, parent='PHASE_Y', tag='PHASE_SCATTER', label='{} points per decade'.format(points_per_decade))
                    dpg.add_plot_legend(location=2)
                    dpg.bind_item_theme("PHASE_PLOT_GRAPH", "PlotPadding")
                    
            # Create a window for the FFT plot
            # with dpg.window(tag='FFT_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, 2*plot_window_height), no_close=True, no_collapse=True, no_move=True, no_title_bar=True, no_resize=True):
            with dpg.child_window(tag='FFT_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, 2 * plot_window_height), border=False, menubar=False):
                with dpg.plot(tag='FFT_PLOT_GRAPH', label='FTT Plot', height=plot_window_height, width=plot_window_width, crosshairs=False, anti_aliased=True, no_mouse_pos=False):
                    dpg.add_plot_axis(dpg.mvXAxis, label='Frequency (Hz)', scale=dpg.mvPlotScale_Log10, tag='FFT_X')
                    # dpg.set_axis_limits doesn't work with scale=dpg.mvPlotScale_Log10 as zero is outside log10 input range.
                    dpg.set_axis_limits('FFT_X', 100, 1000000)
                    dpg.add_plot_axis(dpg.mvYAxis, label='Magnitude (dB)', tag='FFT_Y')
                    dpg.add_drag_line(tag='DRAG_LINE_FFT_X', parent='FFT_PLOT_GRAPH', label="Frequency", color=[255, 255, 0, 255], default_value=0, thickness=2, show=False)
                    dpg.add_drag_line(tag='DRAG_LINE_FFT_Y', parent='FFT_PLOT_GRAPH', label="Magnitude", color=[255, 255, 0, 255], vertical=False, default_value=0, thickness=2, show=False)
                    dpg.add_scatter_series(x=FFTx, y=FFTy, parent='FFT_Y', tag='FTT_SCATTEROUT', label='FFT of Voltage Signal Out')
                    dpg.add_scatter_series(x=MxFFTx, y=MxFFTy, parent='FFT_Y', tag='MxFTT_SCATTEROUT', label='MaxFFT V')
                    dpg.add_plot_legend(location=2)
                    
                    """
                    ===============   =============
                    Location String   Location Code
                    ===============   =============
                    'best'            
                    'upper right'     9
                    'upper left'      5
                    'lower left'      6
                    'lower right'     10
                    'right'           
                    'center left'     4 or 7
                    'center right'    8
                    'center center'   
                    'lower center'    2
                    'upper center'    1
                    'center'          0
                    ===============   =============
                    """
                    # Apply theme to series
                    dpg.bind_item_theme("FTT_SCATTEROUT", "plot_theme_green_blue")
                    dpg.bind_item_theme("FFT_PLOT_GRAPH", "PlotPadding")

            # Create a window for the Oscilloscope plot
            # with dpg.window(tag='OSC_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, 3*plot_window_height), no_close=True, no_collapse=True, no_move=True, no_title_bar=True, no_resize=True):
            with dpg.child_window(tag='OSC_PLOT_WIN', height=plot_window_height, width=plot_window_width, pos=(setting_window_width, 3 * plot_window_height), border=False, menubar=False):
                with dpg.plot(tag='OSC_PLOT_GRAPH', label='Oscilloscope Plot', height=plot_window_height, width=plot_window_width, crosshairs=False, anti_aliased=True, no_mouse_pos=False):
                    dpg.add_plot_axis(dpg.mvXAxis, label='Time', tag='OSC_X')
                    dpg.add_plot_axis(dpg.mvYAxis, label='Amplitude', tag='OSC_Y')
                    dpg.set_axis_limits('OSC_Y', -2, 2)
                    dpg.add_drag_line(tag='DRAG_LINE_OSC_X', parent='OSC_PLOT_GRAPH', label="Time", color=[255, 255, 0, 255], default_value=0, thickness=2, show=False)
                    dpg.add_drag_line(tag='DRAG_LINE_OSC_Y', parent='OSC_PLOT_GRAPH', label="Amplitude", color=[255, 255, 0, 255], vertical=False, default_value=1, thickness=2, show=False)
                    dpg.add_plot_legend()
                    dpg.add_scatter_series(x=OSCx, y=OSCy, parent='OSC_Y', tag='OSC_SCATTEROUT', label='Data Points OUT')
                    dpg.add_scatter_series(x=OSCx, y=OSCy, parent='OSC_Y', tag='OSC_SCATTERIN', label='Data Points IN')
                    dpg.add_plot_legend(location=2)
                    # Apply theme to series
                    dpg.bind_item_theme("OSC_SCATTEROUT", "plot_theme_green_blue")
                    dpg.bind_item_theme("OSC_SCATTERIN", "plot_theme_org_yellow")
                    dpg.bind_item_theme("OSC_PLOT_GRAPH", "PlotPadding")

            # Create a window for the XY plot
            #with dpg.window(tag='XY_PLOT_WIN', height=plot_window_height, width=setting_window_width, pos=(0, 3*plot_window_height), no_close=True, no_collapse=True, no_move=True, no_title_bar=True, no_resize=True):
            with dpg.child_window(tag='XY_PLOT_WIN', height=plot_window_height, width=setting_window_width, pos=(0, 3 * plot_window_height), border=False, menubar=False):
                with dpg.plot(tag='XY_PLOT_GRAPH', label='XY Plot', height=plot_window_height, width=setting_window_width, crosshairs=False, anti_aliased=True, no_mouse_pos=False):
                    dpg.add_plot_axis(dpg.mvXAxis, label='X = IN', tag='XY_X')
                    dpg.set_axis_limits('XY_X', -2, 2)
                    dpg.add_plot_axis(dpg.mvYAxis, label='Y = OUT', tag='XY_Y')
                    dpg.set_axis_limits('XY_Y', -2, 2)
                    dpg.add_plot_legend()
                    dpg.add_scatter_series(x=XYx, y=XYy, parent='XY_Y', tag='XY_SCATTEROUT', label='XY for IN OUT')

    # Apply the default theme (Dark)
    dpg.bind_theme(dark_theme)
    
    # Step starting button colors
    dpg.bind_item_theme('SEARCH_OSCILLOSCOPE', 'GreenButton')    
    dpg.bind_item_theme('START_MEASURE', 'DisabledButton')
    dpg.bind_item_theme('Play', 'DisabledButton')
    dpg.bind_item_theme('Stop', 'DisabledButton')
    dpg.bind_item_theme('JSONLOGFILEtag', 'GreenButton')
    
    # Create a viewport with specified settings
    dpg.create_viewport(title='HDS320S Magnitude/Phase Bode Plotter',
                        x_pos=100,
                        y_pos=setY_pos,
                        width=viewport_width + scrollbar_width,
                        height=set_window_height + scrollbar_width,
                        resizable=True,
                        max_height=main_window_height + scrollbar_width,
                        min_height=int(main_window_height/4),
                        max_width=viewport_width + scrollbar_width,
                        min_width=viewport_width + scrollbar_width,
                        clear_color=[150, 150, 150, 255]
    )

    dpg.setup_dearpygui()
    dpg.show_viewport()
    # dpg.set_primary_window('main', True)

    while dpg.is_dearpygui_running():
        # Create and update all the system variables
        channel_in = str(dpg.get_value(item='CH_IN'))
        channel_out = str(dpg.get_value(item='CH_OUT'))
        CH1_probe_attenuation_ratio = str(dpg.get_value(item='CH1_ATTENUATION_RATIO'))
        CH2_probe_attenuation_ratio = str(dpg.get_value(item='CH2_ATTENUATION_RATIO'))
        CH1_coupling = str(dpg.get_value(item='CH1_COUPLING_MODE'))
        CH2_coupling = str(dpg.get_value(item='CH2_COUPLING_MODE'))
        Sample_command = str(dpg.get_value(item='SAMPL_MODE'))
        DEPMEM = str(dpg.get_value(item='DEPMEM'))
        waveform_amplitude_V = float(dpg.get_value(item='AWG_OUT_VOLTAGE'))
        AWG_output_impedance = str(dpg.get_value(item='HIGH_Z'))
        points_per_decade = int(dpg.get_value(item='POINTS_X_DEC'))
        start_decade = int(decades_list_string.index(dpg.get_value(item='START_DEC')))
        stop_decade = int(decades_list_string.index(dpg.get_value(item='STOP_DEC')))
        point_resize_factor = float(dpg.get_value(item='POINT_SCALE_COEFF'))
        vertical_scaling_factor = float(dpg.get_value(item='V_SCALE_COEFF'))  # used for optimal vertical scale calibration
        horizontal_scaling_factor = float(dpg.get_value(item='H_SCALE_COEFF'))  # used for optimal horizontal scale calibration
        nWaveSamples = int(dpg.get_value(item='SAMPLE_WAVES'))  # used to set the number of times wave data is sampled to average out noise.
        FTTcorection = float(dpg.get_value(item='FFToffset'))  # 2025/02/27 Comparative Measurement of HDS320S with Rigol DS1054Z
        addpm =  dpg.get_value('ADDPM') # Acceptable Degree Delta Per Measurement default 10 with range 0 to 180.
        MaxFreqAdjust = dpg.get_value('MAXTFREQADJUST') # 2025/04/27 Maximum Frequency adjustment attepts to get a good phase value.
        read_delay_ms = int(dpg.get_value(item='OSCILL_TIMEOUT'))
        sample_delay_s = float(dpg.get_value(item='CODE_EXEC_PAUSE'))
        # Plot parameters
        # if platform.system() == 'Windows':
        plot_win_disposition = str(dpg.get_value(item='PLOT_WIN_SETTING'))
        #else:
        # plot_win_disposition = default_value=plot_win_settings[0]
        # Keep updating plot sizes
        dpg.set_item_height(item='MAG_PLOT_GRAPH', height=dpg.get_item_height(item='MAG_PLOT_WIN'))
        dpg.set_item_width(item='MAG_PLOT_GRAPH', width=dpg.get_item_width(item='MAG_PLOT_WIN'))
        dpg.set_item_height(item='PHASE_PLOT_GRAPH', height=dpg.get_item_height(item='PHASE_PLOT_WIN'))
        dpg.set_item_width(item='PHASE_PLOT_GRAPH', width=dpg.get_item_width(item='PHASE_PLOT_WIN'))

        dpg.render_dearpygui_frame()
    # Perform cleanup before destroying context
    on_viewport_close()
    dpg.destroy_context()
# --- End GUI settings --- ----------------------------------------

