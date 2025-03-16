"""

© 2025 by Bode Plotter

MIT License

Copyright (c) [2025] [Bode Plotter]
Original code written by Simone Albano on 10/15/2023
Bode plot implementation for the handheld oscilloscope OWON HDS320S
Highly rewritton and modified between 01/20/2025 to 03/15/2025
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
OUT OF OR IN CO

"""

import usb.core
import usb.util
import numpy as np
from numpy.fft import fft, fftshift, fftfreq
import cmath
import matplotlib
matplotlib.use("Qt5Agg")  # Switch to Qt5Agg backend for better performance
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import UnivariateSpline
import scipy.optimize as optimize
import dearpygui.dearpygui as dpg
import datetime
import time
import os
from pathlib import Path
import re  # Used for extracting the first number (integer or float) found in a string.
from scipy.optimize import curve_fit
import multiprocessing
import queue
import threading
import json
import traceback
from collections import deque
from screeninfo import get_monitors
import platform
import warnings
# 2025/01/30 
# /usr/lib/python3.12/contextlib.py:137: DeprecationWarning: anti_aliased keyword removed
#  return next(self.gen)
warnings.filterwarnings("ignore", category=DeprecationWarning)

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

# Global system variables
global channel_in, channel_out, CH1_probe_attenuation_ratio, CH2_probe_attenuation_ratio
global CH1_coupling, CH2_coupling, Sample_command, DEPMEM, waveform_amplitude_V
global AWG_output_impedance, points_per_decade, start_decade, stop_decade
global point_resize_factor, vertical_scaling_factor, horizontal_scaling_factor
global nWaveSamples, FTTcorection, read_delay_ms, sample_delay_s, plot_win_disposition, LogFile
global JSONLogFile, plot_process, oscilloscope_OUT, oscilloscope_IN, oscilloscope
global is_playing, is_paused, play_speed, is_recording

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
amplitude_scales_commands = [
    '100V', '50V', '10.0V', '5.00V', '2.00V', '1.00V', '500mV', '200mV', '100mV', '50.0mV', '20.0mV', '10.0mV'
]
amplitude_scales_values = [
    100, 50, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01
]
decades_list_string = [
    '1mHz', '1Hz', '10Hz', '100Hz', '1kHz', '10kHz', '100kHz', '1MHz', '10MHz'
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

# --- Set working path --- ------------------------------------------------------------
# Get the path to the user's Documents directory
documents_path = os.path.join(os.path.expanduser("~"), "Documents")

# Set the working directory to the Documents directory
os.chdir(documents_path)
# --- End Set working path --- --------------------------------------------------------

"""
# --- class PlotManager Magnitude Phase using FuncAnimation:--- -----------------------
class PlotManagerMagnitudePhase:
    def __init__(self):
        self.data_queue = multiprocessing.Queue(maxsize=100)  # Limited queue size for control
        self.command_queue = multiprocessing.Queue()  # Queue for commands
        self.stop_event = multiprocessing.Event()  # Event for stopping the process
        self.plot_process = None  # Plotting process
        # self.gain_X =([])  # Store all gain X data, ([]) Initialize as a NumPy array
        # self.gain_Y = ([])  # Store all gain Y data, ([]) Initialize as a NumPy array
        # self.phase_X = ([])  # Store all phase X data, ([]) Initialize as a NumPy array
        # self.phase_Y = ([])  # Store all phase Y data, ([]) Initialize as a NumPy array
        # Initialize variables as lists
        self.gain_X = []  # Store all gain X data
        self.gain_Y = []  # Store all gain Y data
        self.phase_X = []  # Store all phase X data
        self.phase_Y = []  # Store all phase Y data
    def create_plot(self, data_queue, command_queue, stop_event, start_decade, stop_decade, points_per_decade):
        try:
            plt.ion()
            figure, (magnitude, phase) = plt.subplots(2)
            
            figure.canvas.manager.set_window_title("Magnitude / Phase Bode Plot")

            # Configure magnitude plot
            magnitude.set_xlim([start_decade, stop_decade])
            magnitude.set_xscale('log')
            magnitude.set_title('Magnitude Bode Plot')
            magnitude.set_xlabel('Frequency (Hz)')
            magnitude.set_ylabel('Gain (dB)')
            magnitude.grid(which='both')

            # Configure phase plot
            phase.set_xlim([start_decade, stop_decade])
            phase.set_xscale('log')
            phase.set_title('Phase Bode Plot')
            phase.set_xlabel('Frequency (Hz)')
            phase.set_ylabel('Phase (Degrees)')
            phase.set_ylim([-200, 200])  # Fixed Y-axis range for the phase plot
            phase.grid(which='both')

            # Initialize scatter plots
            magnitude_scatter = magnitude.scatter([0], [0], color='red', label='Raw Data', s=20)
            phase_scatter = phase.scatter([0], [0], color='orange', label='Raw Data', s=20)

            figure.tight_layout()
            figure.canvas.draw()
        except Exception as e:
            print(f"Exception in create_plot for PlotManagerMagnitudePhase: {e}")
                  
        def update_plot(i):
            try:
                if not data_queue.empty():  # Only update if data is available
                    # Retrieve updated data
                    gain_X, gain_Y, phase_X, phase_Y = data_queue.get_nowait()
                    
                    # Reset the scatter plots with the first received data
                    magnitude_scatter.set_offsets(np.c_[gain_X, gain_Y])
                    phase_scatter.set_offsets(np.c_[phase_X, phase_Y])

                    # Update Y-axis limits dynamically for magnitude
                    if gain_Y:
                        magnitude.set_ylim([min(gain_Y) - 5, max(gain_Y) + 5])
            except queue.Empty:
                pass  # If the queue is empty, do nothing

            # Handle commands (e.g., trigger line plot)
            try:
                while not command_queue.empty():
                    command = command_queue.get_nowait()
                    if command == "trigger_line_plot":
                        # Clear previous line plot to avoid duplicate labels for magnitude
                        for line in magnitude.get_lines():
                            line.remove()

                        # Add the line plot for the magnitude
                        if gain_X and len(gain_X) >= 3:
                            # Ensure gain_X is sorted and strictly increasing
                            sorted_indices = np.argsort(gain_X)
                            gain_X = np.array(gain_X)[sorted_indices]
                            gain_Y = np.array(gain_Y)[sorted_indices]

                            # Remove duplicates in gain_X
                            unique_indices = np.diff(gain_X) > 0
                            gain_X = gain_X[np.concatenate(([True], unique_indices))]
                            gain_Y = gain_Y[np.concatenate(([True], unique_indices))]

                            magnitude_spline = UnivariateSpline(gain_X, gain_Y, k=3, s=0)
                            magnitude.plot(gain_X, magnitude_spline(gain_X), color='blue', label='Magnitude Spline')

                        # Clear previous line plot for the phase plot
                        for line in phase.get_lines():
                            line.remove()

                        # Add the line plot for the phase
                        if  phase_X and len(phase_X) >= 3:
                            # Ensure phase_X is sorted and strictly increasing
                            sorted_indices = np.argsort(phase_X)
                            phase_X = np.array(phase_X)[sorted_indices]
                            phase_Y = np.array(phase_Y)[sorted_indices]

                            # Remove duplicates in phase_X
                            unique_indices = np.diff(phase_X) > 0
                            phase_X = phase_X[np.concatenate(([True], unique_indices))]
                            phase_Y = phase_Y[np.concatenate(([True], unique_indices))]
                            
                            # Outlier detection based on local deviation
                            window_size = 5  # Number of neighbors to consider on each side
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

                            # Fit spline only to inliers
                            if len(filtered_phase_X) >= 3:  # Ensure enough points for spline
                                phase_spline = UnivariateSpline(filtered_phase_X, filtered_phase_Y, k=3, s=0)
                                phase.plot(filtered_phase_X, phase_spline(filtered_phase_X), color='blue', label='Phase Spline')

                        # Update legend for the magnitude plot and phase plot
                        magnitude.legend()
                        phase.legend()
            except queue.Empty:
                time.sleep(0.1)  # Sleep for 0.1 seconds when the queue is empty

            return magnitude_scatter, phase_scatter

        ani = FuncAnimation(figure, update_plot, interval=100, blit=False, cache_frame_data=False)
        plt.show(block=True)

    def start_plot_process(self, start_decade, stop_decade, points_per_decade):
        self.plot_process = multiprocessing.Process(
            target=self.create_plot,
            args=(self.data_queue, self.command_queue, self.stop_event, start_decade, stop_decade, points_per_decade)
        )
        self.plot_process.start()

    def set_data_point(self, new_gain_X, new_gain_Y, new_phase_X, new_phase_Y):
        # Append new data points to the accumulated dataset
        self.gain_X.append(new_gain_X)
        self.gain_Y.append(new_gain_Y)
        self.phase_X.append(new_phase_X)
        self.phase_Y.append(new_phase_Y)

        # Send the accumulated dataset to the data queue
        self.data_queue.put((self.gain_X, self.gain_Y, self.phase_X, self.phase_Y))

    def set_full_data(self, gain_X, gain_Y, phase_X, phase_Y):
        # Flush the data queue by processing all single points
        print("Flushing the data queue...")
        while not self.data_queue.empty():
            try:
                self.data_queue.get_nowait()
            except queue.Empty:
                break
            time.sleep(0.01)  # Short delay to allow smooth processing without blocking

        # Overwrite the accumulated dataset with the new full dataset
        self.gain_X = np.array(gain_X)
        self.gain_Y = np.array(gain_Y)
        self.phase_X = np.array(phase_X)
        self.phase_Y = np.array(phase_Y)

        # Send the full dataset to the data queue
        self.data_queue.put((gain_X, gain_Y, phase_X, phase_Y))

        # Add a command to trigger the line plot
        self.command_queue.put("trigger_line_plot")

    def stop_plot_process(self):
        if self.plot_process and self.plot_process.is_alive():
            self.stop_event.set()
            self.plot_process.join()

# --- End class PlotManager Magnitude Phase using FuncAnimation:--- -------------------
"""

# --- class PlotManager Magnitude Phase:--- -------------------------------------------
class PlotManagerMagnitudePhase:
    def __init__(self):    
        self.data_queue = multiprocessing.Queue(maxsize=100)  # Limited queue size for control
        self.command_queue = multiprocessing.Queue()  # Queue for commands
        self.stop_event = multiprocessing.Event()  # Event for stopping the process
        self.plot_process = None  # Plotting process
        # self.gain_X =([])  # Store all gain X data, ([]) Initialize as a NumPy array
        # self.gain_Y = ([])  # Store all gain Y data, ([]) Initialize as a NumPy array
        # self.phase_X = ([])  # Store all phase X data, ([]) Initialize as a NumPy array
        # self.phase_Y = ([])  # Store all phase Y data, ([]) Initialize as a NumPy array
        # Initialize variables as lists
        self.gain_X = []  # Store all gain X data
        self.gain_Y = []  # Store all gain Y data
        self.phase_X = []  # Store all phase X data
        self.phase_Y = []  # Store all phase Y data
        
    def create_plot(self, data_queue, command_queue, stop_event, start_decade, stop_decade, points_per_decade):
        try:
            plt.ion()
            figure, (magnitude, phase) = plt.subplots(2)
            
            figure.canvas.manager.set_window_title("Magnitude / Phase Bode Plot")

            # Configure magnitude plot
            magnitude.set_xlim([start_decade, stop_decade])
            magnitude.set_xscale('log')
            magnitude.set_title('Magnitude Bode Plot')
            magnitude.set_xlabel('Frequency (Hz)')
            magnitude.set_ylabel('Gain (dB)')
            magnitude.grid(which='both')

            # Configure phase plot
            phase.set_xlim([start_decade, stop_decade])
            phase.set_xscale('log')
            phase.set_title('Phase Bode Plot')
            phase.set_xlabel('Frequency (Hz)')
            phase.set_ylabel('Phase (Degrees)')
            phase.set_ylim([-200, 200])  # Fixed Y-axis range for the phase plot
            phase.grid(which='both')

            # Initialize scatter plots
            magnitude_scatter = magnitude.scatter([0], [0], color='red', label='Raw Data', s=20)
            phase_scatter = phase.scatter([0], [0], color='orange', label='Raw Data', s=20)

            figure.tight_layout()
            
            # Define the close event handler
            def on_close(event):
                # print("Window is closing...")
                stop_event.set()

            # Bind the close event handler to the figure
            figure.canvas.mpl_connect('close_event', on_close)

        except Exception as e:
            print(f"Exception in create_plot for PlotManagerMagnitudePhase: {e}")
                  
        while not stop_event.is_set():
            try:
                if not data_queue.empty():  # Only update if data is available
                    # Retrieve updated data
                    gain_X, gain_Y, phase_X, phase_Y = data_queue.get_nowait()
                    
                    # Reset the scatter plots with the first received data
                    magnitude_scatter.set_offsets(np.c_[gain_X, gain_Y])
                    phase_scatter.set_offsets(np.c_[phase_X, phase_Y])

                    # Update Y-axis limits dynamically for magnitude
                    if gain_Y:
                        magnitude.set_ylim([min(gain_Y) - 5, max(gain_Y) + 5])
                    figure.canvas.draw()
                    figure.canvas.flush_events()
            except queue.Empty:
                pass  # If the queue is empty, do nothing
                # Add a slight delay to prevent excessive CPU usage
                time.sleep(0.1)

            # Handle commands (e.g., trigger line plot)
            try:
                while not command_queue.empty():
                    command = command_queue.get_nowait()
                    if command == "trigger_line_plot":
                        # Clear previous line plot to avoid duplicate labels for magnitude
                        for line in magnitude.get_lines():
                            line.remove()

                        # Add the line plot for the magnitude
                        if gain_X and len(gain_X) >= 3:
                            # Ensure gain_X is sorted and strictly increasing
                            sorted_indices = np.argsort(gain_X)
                            gain_X = np.array(gain_X)[sorted_indices]
                            gain_Y = np.array(gain_Y)[sorted_indices]

                            # Remove duplicates in gain_X
                            unique_indices = np.diff(gain_X) > 0
                            gain_X = gain_X[np.concatenate(([True], unique_indices))]
                            gain_Y = gain_Y[np.concatenate(([True], unique_indices))]

                            magnitude_spline = UnivariateSpline(gain_X, gain_Y, k=3, s=0)
                            magnitude.plot(gain_X, magnitude_spline(gain_X), color='blue', label='Magnitude Spline')

                        # Clear previous line plot for the phase plot
                        for line in phase.get_lines():
                            line.remove()

                        # Add the line plot for the phase
                        if  phase_X and len(phase_X) >= 3:
                            # Ensure phase_X is sorted and strictly increasing
                            sorted_indices = np.argsort(phase_X)
                            phase_X = np.array(phase_X)[sorted_indices]
                            phase_Y = np.array(phase_Y)[sorted_indices]

                            # Remove duplicates in phase_X
                            unique_indices = np.diff(phase_X) > 0
                            phase_X = phase_X[np.concatenate(([True], unique_indices))]
                            phase_Y = phase_Y[np.concatenate(([True], unique_indices))]
                            
                            # Outlier detection based on local deviation
                            window_size = 5  # Number of neighbors to consider on each side
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

                            # Fit spline only to inliers
                            if len(filtered_phase_X) >= 3:  # Ensure enough points for spline
                                phase_spline = UnivariateSpline(filtered_phase_X, filtered_phase_Y, k=3, s=0)
                                phase.plot(filtered_phase_X, phase_spline(filtered_phase_X), color='blue', label='Phase Spline')

                        # Update legend for the magnitude plot and phase plot
                        magnitude.legend()
                        phase.legend()
                        figure.canvas.draw()
                        figure.canvas.flush_events()                       
            except queue.Empty:
                time.sleep(0.1)  # Sleep for 0.1 seconds when the queue is empty
            # Allow GUI events to process
            figure.canvas.flush_events()
            time.sleep(0.1)


    def start_plot_process(self, start_decade, stop_decade, points_per_decade):
        self.plot_process = multiprocessing.Process(
            target=self.create_plot,
            args=(self.data_queue, self.command_queue, self.stop_event, start_decade, stop_decade, points_per_decade)
        )
        self.plot_process.start()

    def set_data_point(self, new_gain_X, new_gain_Y, new_phase_X, new_phase_Y):
        # Append new data points to the accumulated dataset
        self.gain_X.append(new_gain_X)
        self.gain_Y.append(new_gain_Y)
        self.phase_X.append(new_phase_X)
        self.phase_Y.append(new_phase_Y)

        # Send the accumulated dataset to the data queue
        self.data_queue.put((self.gain_X, self.gain_Y, self.phase_X, self.phase_Y))

    def set_full_data(self, gain_X, gain_Y, phase_X, phase_Y):
        # Flush the data queue by processing all single points
        print("Flushing the data queue...")
        while not self.data_queue.empty():
            try:
                self.data_queue.get_nowait()
            except queue.Empty:
                break
            time.sleep(0.01)  # Short delay to allow smooth processing without blocking

        # Overwrite the accumulated dataset with the new full dataset
        self.gain_X = np.array(gain_X)
        self.gain_Y = np.array(gain_Y)
        self.phase_X = np.array(phase_X)
        self.phase_Y = np.array(phase_Y)

        # Send the full dataset to the data queue
        self.data_queue.put((gain_X, gain_Y, phase_X, phase_Y))

        # Add a command to trigger the line plot
        self.command_queue.put("trigger_line_plot")

    def stop_plot_process(self):
        if self.plot_process and self.plot_process.is_alive():
            self.stop_event.set()
            self.plot_process.join()

# --- End class PlotManager Magnitude Phase:--- ---------------------------------------

"""
# --- class PlotManager FFT oscilloscope using FuncAnimation:--- ----------------------

class PlotManagerFFToscilloscope:
    def __init__(self):
        self.data_queue = multiprocessing.Queue(maxsize=100)  # Limited queue size for control
        self.stop_event = multiprocessing.Event()  # Event for stopping the process
        self.plot_process = None  # Plotting process
        self.FFTxMax = ([])  # Store all Maximum X values of the FFT, ([]) Initialize as a NumPy array
        self.FFTmaxVal = ([])  # Store all Maximum values of the FFT, ([]) Initialize as a NumPy array
        self.scatter_updated = False  # Flag to indicate scatter plot update

    def create_plot(self, data_queue, stop_event, start_decade, stop_decade, points_per_decade):
        try:
            plt.ion()
            figure, (FFT, oscilloscope) = plt.subplots(2)
            
            figure.canvas.manager.set_window_title("Max FFT, FFT, and Oscilloscope In/Out Plot")

            # Configure FFT plot
            FFT.set_xlim([start_decade, stop_decade])
            FFT.set_xscale('log')
            FFT.set_title('FFT Plot')
            FFT.set_xlabel('Frequency (Hz)')
            FFT.set_ylabel('Gain (dB)')
            FFT.grid(which='both')

            # Configure oscilloscope plot
            oscilloscope.set_xlim([0, 1])  # Default time range; can be updated dynamically
            oscilloscope.set_xscale('linear')
            oscilloscope.set_title('Oscilloscope Plot')
            oscilloscope.set_xlabel('Time (s)')
            oscilloscope.set_ylabel('Voltage (V)')
            oscilloscope.grid(which='both')

           # Initialize scatter plots
            FFT_scatter = FFT.scatter([0], [0], color='blue', label='FFT Data', s=20)
            FFT_max_scatter = FFT.scatter([0], [0], color='red', label='FFT Max Data', s=20)
            oscilloscope_in_scatter = oscilloscope.scatter([0], [0], color='yellow', label='Input Data', s=20)
            oscilloscope_out_scatter = oscilloscope.scatter([0], [0], color='blue', label='Output Data', s=20)
            
            # Adding legend in the very left corner, outside the plot area
            FFT.legend(loc='lower left', bbox_to_anchor=(0, 0), borderaxespad=0, frameon=True)
            oscilloscope.legend(loc='lower left', bbox_to_anchor=(0, 0), borderaxespad=0, frameon=True)
            
            figure.tight_layout()
            figure.canvas.draw()
        except Exception as e:
            print(f"Exception in create_plot for PlotManagerFFToscilloscope: {e}")
 
        # Initialize update_count as a global variable
        global update_count
        update_count = -1
        
        def update_plot(i):
            global update_count  # To modify the update_count variable
            try:
                if not data_queue.empty():  # Only update if data is available
                    update_count += 1
                    
                    # Adjust the skip interval based on the current size of the queue
                    queue_size = data_queue.qsize()
                    skip_interval = min(30, 10 + max(1, queue_size // 50))
                    data = data_queue.get_nowait()
                     # Retrieve updated data
                    FFTxMax, FFTmaxVal, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout = map(np.array, data)
                    if FFTxMax.any() and FFTmaxVal.any():
                        FFT_max_scatter.set_offsets(np.c_[FFTxMax, FFTmaxVal])
                        # self.scatter_updated = True
                    # Skip plotting unless update_count is a multiple of skip_interval
                    if update_count % skip_interval == 0:
                        # Clear previous line plot to avoid duplicate labels for FFT
                        for line in FFT.get_lines():
                            line.remove()

                        if FFTx.any() and FFTy.any():
                            FFT_scatter.set_offsets(np.c_[FFTx, FFTy])
                            self.scatter_updated = True

                        if OSCxin.any() and OSCyin.any():
                            oscilloscope_in_scatter.set_offsets(np.c_[OSCxin, OSCyin])
                            self.scatter_updated = True

                        if OSCxout.any() and OSCyout.any():
                            oscilloscope_out_scatter.set_offsets(np.c_[OSCxout, OSCyout])
                            self.scatter_updated = True

                        # Update Y-axis limits dynamically for FFT
                        if FFTy.any() and FFTmaxVal.any():
                            FFTy_all = np.concatenate((FFTy, FFTmaxVal))
                            if FFTy_all.size:
                                FFT.set_ylim([min(FFTy_all) - 5, max(FFTy_all) + 5])

                        # Update X-axis limits dynamically for Oscilloscope
                        if OSCxin.any() or OSCxout.any():
                            oscilloscope_x_all = np.concatenate((OSCxin, OSCxout))
                            if oscilloscope_x_all.size:
                                oscilloscope.set_xlim([min(oscilloscope_x_all), max(oscilloscope_x_all)])
                else:
                    # Reset update_count to allow to print the next time there is data.
                    update_count = -1
            except queue.Empty:
                pass  # If the queue is empty, do nothing
                 
            except Exception as e:
                print(f"Error in handling command: {e}")
                
            # Handle commands (e.g., trigger line plot)
            try:

                if self.scatter_updated:
                    # Clear previous line plot to avoid duplicate labels for FFT
                    for line in FFT.get_lines():
                        line.remove()

                    # Add the line plot for the FFT
                    if FFTx.any() and FFTy.any() and len(FFTx) >= 6:
                        try:
                            # Ensure FFTx is sorted and strictly increasing
                            sorted_indices = np.argsort(FFTx)
                            FFTx = np.array(FFTx)[sorted_indices]
                            FFTy = np.array(FFTy)[sorted_indices]
                            unique_indices = np.diff(FFTx) > 0
                            FFTx = FFTx[np.concatenate(([True], unique_indices))]
                            FFTy = FFTy[np.concatenate(([True], unique_indices))]

                            FFT_spline = UnivariateSpline(FFTx, FFTy, k=3, s=0)
                            FFT.plot(FFTx, FFT_spline(FFTx), color='blue', label='FFT Spline')
                        except Exception as e:
                            print(f"Error in fitting spline for FFT data: {e}")

                    # Add the line plot for the FFT max data
                    if FFTxMax.any() and FFTmaxVal.any() and len(FFTxMax) >= 6:
                        try:
                            # Ensure FFTxMax is sorted and strictly increasing
                            sorted_indices = np.argsort(FFTxMax)
                            FFTxMax = np.array(FFTxMax)[sorted_indices]
                            FFTmaxVal = np.array(FFTmaxVal)[sorted_indices]
                            unique_indices = np.diff(FFTxMax) > 0
                            FFTxMax = FFTxMax[np.concatenate(([True], unique_indices))]
                            FFTmaxVal = FFTmaxVal[np.concatenate(([True], unique_indices))]

                            FFT_max_spline = UnivariateSpline(FFTxMax, FFTmaxVal, k=3, s=0)
                            FFT.plot(FFTxMax, FFT_max_spline(FFTxMax), color='red', label='FFT Max Spline')
                        except Exception as e:
                            print(f"Error in fitting spline for FFT max data: {e}")

                    # Clear previous line plot for the oscilloscope plot
                    for line in oscilloscope.get_lines():
                        line.remove()

                    # Add the line plot for the oscilloscope input
                    if OSCxin.any() and OSCyin.any() and len(OSCxin) >= 6:
                        oscilloscope_in_spline = UnivariateSpline(OSCxin, OSCyin, k=3, s=0)
                        oscilloscope.plot(OSCxin, oscilloscope_in_spline(OSCxin), color='yellow', label='Input Spline')

                    # Add the line plot for the oscilloscope output
                    if OSCxout.any() and OSCyout.any() and len(OSCxout) >= 6:
                        oscilloscope_out_spline = UnivariateSpline(OSCxout, OSCyout, k=3, s=0)
                        oscilloscope.plot(OSCxout, oscilloscope_out_spline(OSCxout), color='blue', label='Output Spline')

                    # Reset the scatter update flag
                    self.scatter_updated = False
            except queue.Empty:
                time.sleep(0.1)  # Sleep for 0.1 seconds when the queue is empty
                 
            except Exception as e:
                print(f"Error in handling command: {e}")
                
            return FFT_scatter, FFT_max_scatter, oscilloscope_in_scatter, oscilloscope_out_scatter

        ani = FuncAnimation(figure, update_plot, interval=100, blit=False, cache_frame_data=False)
        plt.show(block=True)

    def start_plot_process(self, start_decade, stop_decade, points_per_decade):
        self.plot_process = multiprocessing.Process(
            target=self.create_plot,
            args=(self.data_queue, self.stop_event, start_decade, stop_decade, points_per_decade)
        )
        self.plot_process.start()

    def set_data_point(self, new_FFTxMax, new_FFTmaxVal):
        # Append new data points to the accumulated dataset
        self.FFTxMax.append(new_FFTxMax)
        self.FFTmaxVal.append(new_FFTmaxVal)

    def set_full_data(self, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout):
        # Ensure the FFTxMax and FFTmaxVal arrays are not cleared
        # Ensure data is a NumPy array
        FFTx = np.array(FFTx)
        FFTy = np.array(FFTy)
        OSCxin = np.array(OSCxin)
        OSCyin = np.array(OSCyin)
        OSCxout = np.array(OSCxout)
        OSCyout = np.array(OSCyout)
        # Send the full dataset to the data queue
        self.data_queue.put((self.FFTxMax, self.FFTmaxVal, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout))

    def stop_plot_process(self):
        if self.plot_process and self.plot_process.is_alive():
            self.stop_event.set()
            self.plot_process.join()

# --- End class PlotManager FFT oscilloscope using FuncAnimation:--- ------------------
"""

# --- class PlotManager FFT oscilloscope:--- ------------------------------------------

class PlotManagerFFToscilloscope:
    def __init__(self):
        self.data_queue = multiprocessing.Queue(maxsize=100)  # Limited queue size for control
        self.stop_event = multiprocessing.Event()  # Event for stopping the process
        self.plot_process = None  # Plotting process
        self.FFTxMax = ([])  # Store all Maximum X values of the FFT, ([]) Initialize as a NumPy array
        self.FFTmaxVal = ([])  # Store all Maximum values of the FFT, ([]) Initialize as a NumPy array
        self.scatter_updated = False  # Flag to indicate scatter plot update
        
    def create_plot(self, data_queue, stop_event, start_decade, stop_decade, points_per_decade):
        try:
            plt.ion()
            figure, (FFT, oscilloscope) = plt.subplots(2)
            
            figure.canvas.manager.set_window_title("Max FFT, FFT, and Oscilloscope In/Out Plot")

            # Configure FFT plot
            FFT.set_xlim([start_decade, stop_decade])
            FFT.set_xscale('log')
            FFT.set_title('FFT Plot')
            FFT.set_xlabel('Frequency (Hz)')
            FFT.set_ylabel('Gain (dB)')
            FFT.grid(which='both')

            # Configure oscilloscope plot
            oscilloscope.set_xlim([0, 1])  # Default time range; can be updated dynamically
            oscilloscope.set_xscale('linear')
            oscilloscope.set_title('Oscilloscope Plot')
            oscilloscope.set_xlabel('Time (s)')
            oscilloscope.set_ylabel('Voltage (V)')
            oscilloscope.grid(which='both')

           # Initialize scatter plots
            FFT_scatter = FFT.scatter([0], [0], color='blue', label='FFT Data', s=20)
            FFT_max_scatter = FFT.scatter([0], [0], color='red', label='FFT Max Data', s=20)
            oscilloscope_in_scatter = oscilloscope.scatter([0], [0], color='yellow', label='Input Data', s=20)
            oscilloscope_out_scatter = oscilloscope.scatter([0], [0], color='blue', label='Output Data', s=20)
            
            # Adding legend in the very left corner, outside the plot area
            FFT.legend(loc='lower left', bbox_to_anchor=(0, 0), borderaxespad=0, frameon=True)
            oscilloscope.legend(loc='lower left', bbox_to_anchor=(0, 0), borderaxespad=0, frameon=True)
            
            figure.tight_layout()

        except Exception as e:
            print(f"Exception in create_plot for PlotManagerFFToscilloscope: {e}")
 
        # Initialize update_count as a global variable
        global update_count
        update_count = -1
        
        # Define the close event handler
        def on_close(event):
            # print("Window is closing...")
            stop_event.set()

        # Bind the close event handler to the figure
        figure.canvas.mpl_connect('close_event', on_close)
        
        while not stop_event.is_set():
            try:
                if not data_queue.empty():  # Only update if data is available
                    update_count += 1
                    
                    # Adjust the skip interval based on the current size of the queue
                    queue_size = data_queue.qsize()
                    skip_interval = min(30, 10 + max(1, queue_size // 50))
                    data = data_queue.get_nowait()
                     # Retrieve updated data
                    FFTxMax, FFTmaxVal, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout = map(np.array, data)
                    if FFTxMax.any() and FFTmaxVal.any():
                        FFT_max_scatter.set_offsets(np.c_[FFTxMax, FFTmaxVal])
                        # self.scatter_updated = True
                    # Skip plotting unless update_count is a multiple of skip_interval
                    if update_count % skip_interval == 0:
                        # Clear previous line plot to avoid duplicate labels for FFT
                        for line in FFT.get_lines():
                            line.remove()

                        if FFTx.any() and FFTy.any():
                            FFT_scatter.set_offsets(np.c_[FFTx, FFTy])
                            self.scatter_updated = True

                        if OSCxin.any() and OSCyin.any():
                            oscilloscope_in_scatter.set_offsets(np.c_[OSCxin, OSCyin])
                            self.scatter_updated = True

                        if OSCxout.any() and OSCyout.any():
                            oscilloscope_out_scatter.set_offsets(np.c_[OSCxout, OSCyout])
                            self.scatter_updated = True

                        # Update Y-axis limits dynamically for FFT
                        if FFTy.any() and FFTmaxVal.any():
                            FFTy_all = np.concatenate((FFTy, FFTmaxVal))
                            if FFTy_all.size:
                                FFT.set_ylim([min(FFTy_all) - 5, max(FFTy_all) + 5])

                        # Update X-axis limits dynamically for Oscilloscope
                        if OSCxin.any() or OSCxout.any():
                            oscilloscope_x_all = np.concatenate((OSCxin, OSCxout))
                            if oscilloscope_x_all.size:
                                oscilloscope.set_xlim([min(oscilloscope_x_all), max(oscilloscope_x_all)])
                else:
                    # Reset update_count to allow to print the next time there is data.
                    update_count = -1
                    # Add a slight delay to prevent excessive CPU usage
                    time.sleep(0.1)
            except queue.Empty:
                pass  # If the queue is empty, do nothing
                time.sleep(0.1)  # Sleep if the queue is empty                 
            except Exception as e:
                print(f"Error in handling command: {e}")
                
            # Handle commands (e.g., trigger line plot)
            try:

                if self.scatter_updated:
                    # Clear previous line plot to avoid duplicate labels for FFT
                    for line in FFT.get_lines():
                        line.remove()

                    # Add the line plot for the FFT
                    if FFTx.any() and FFTy.any() and len(FFTx) >= 6:
                        try:
                            # Ensure FFTx is sorted and strictly increasing
                            sorted_indices = np.argsort(FFTx)
                            FFTx = np.array(FFTx)[sorted_indices]
                            FFTy = np.array(FFTy)[sorted_indices]
                            unique_indices = np.diff(FFTx) > 0
                            FFTx = FFTx[np.concatenate(([True], unique_indices))]
                            FFTy = FFTy[np.concatenate(([True], unique_indices))]

                            FFT_spline = UnivariateSpline(FFTx, FFTy, k=3, s=0)
                            FFT.plot(FFTx, FFT_spline(FFTx), color='blue', label='FFT Spline')
                        except Exception as e:
                            print(f"Error in fitting spline for FFT data: {e}")

                    # Add the line plot for the FFT max data
                    if FFTxMax.any() and FFTmaxVal.any() and len(FFTxMax) >= 6:
                        try:
                            # Ensure FFTxMax is sorted and strictly increasing
                            sorted_indices = np.argsort(FFTxMax)
                            FFTxMax = np.array(FFTxMax)[sorted_indices]
                            FFTmaxVal = np.array(FFTmaxVal)[sorted_indices]
                            unique_indices = np.diff(FFTxMax) > 0
                            FFTxMax = FFTxMax[np.concatenate(([True], unique_indices))]
                            FFTmaxVal = FFTmaxVal[np.concatenate(([True], unique_indices))]

                            FFT_max_spline = UnivariateSpline(FFTxMax, FFTmaxVal, k=3, s=0)
                            FFT.plot(FFTxMax, FFT_max_spline(FFTxMax), color='red', label='FFT Max Spline')
                        except Exception as e:
                            print(f"Error in fitting spline for FFT max data: {e}")

                    # Clear previous line plot for the oscilloscope plot
                    for line in oscilloscope.get_lines():
                        line.remove()

                    # Add the line plot for the oscilloscope input
                    if OSCxin.any() and OSCyin.any() and len(OSCxin) >= 6:
                        oscilloscope_in_spline = UnivariateSpline(OSCxin, OSCyin, k=3, s=0)
                        oscilloscope.plot(OSCxin, oscilloscope_in_spline(OSCxin), color='yellow', label='Input Spline')

                    # Add the line plot for the oscilloscope output
                    if OSCxout.any() and OSCyout.any() and len(OSCxout) >= 6:
                        oscilloscope_out_spline = UnivariateSpline(OSCxout, OSCyout, k=3, s=0)
                        oscilloscope.plot(OSCxout, oscilloscope_out_spline(OSCxout), color='blue', label='Output Spline')

                    # Reset the scatter update flag
                    self.scatter_updated = False
                    figure.canvas.draw()
                    figure.canvas.flush_events()
            except queue.Empty:
                time.sleep(0.1)  # Sleep for 0.1 seconds when the queue is empty
                 
            except Exception as e:
                print(f"Error in handling command: {e}")
            # Allow GUI events to process
            figure.canvas.flush_events()
            time.sleep(0.1)


    def start_plot_process(self, start_decade, stop_decade, points_per_decade):
        self.plot_process = multiprocessing.Process(
            target=self.create_plot,
            args=(self.data_queue, self.stop_event, start_decade, stop_decade, points_per_decade)
        )
        self.plot_process.start()

    def set_data_point(self, new_FFTxMax, new_FFTmaxVal):
        # Append new data points to the accumulated dataset
        self.FFTxMax.append(new_FFTxMax)
        self.FFTmaxVal.append(new_FFTmaxVal)

    def set_full_data(self, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout):
        # Ensure the FFTxMax and FFTmaxVal arrays are not cleared
        # Ensure data is a NumPy array
        FFTx = np.array(FFTx)
        FFTy = np.array(FFTy)
        OSCxin = np.array(OSCxin)
        OSCyin = np.array(OSCyin)
        OSCxout = np.array(OSCxout)
        OSCyout = np.array(OSCyout)
        # Send the full dataset to the data queue
        self.data_queue.put((self.FFTxMax, self.FFTmaxVal, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout))

    def stop_plot_process(self):
        if self.plot_process and self.plot_process.is_alive():
            self.stop_event.set()
            self.plot_process.join()

# --- End class PlotManager FFT oscilloscope:--- --------------------------------------

"""
# --- class PlotManager XYoscilloscope using FuncAnimation:--- ------------------------
class XYoscilloscope:
    def __init__(self):
        self.data_queue = multiprocessing.Queue(maxsize=100)  # Limited queue size for control
        self.command_queue = multiprocessing.Queue()  # Queue for commands
        self.stop_event = multiprocessing.Event()  # Event for stopping the process
        self.plot_process = None  # Plotting process

    def create_plot(self, data_queue, command_queue, stop_event):
        try:
            plt.ion()
            figure, (XY_plot) = plt.subplots(1)

            # Set the window title
            figure.canvas.manager.set_window_title("Oscilloscope X-Y Plot")
            
            # Configure XY plot
            XY_plot.set_title('XY Oscilloscope Plot')
            XY_plot.set_xlabel('raw_waveform_in')
            XY_plot.set_ylabel('raw_waveform_out')
            XY_plot.grid(which='both')

            # Initialize scatter plots
            XY_scatter = XY_plot.scatter([0], [0], color='red', label='XY Data', s=20)
            
            # Adding legend in the very left corner, outside the plot area
            XY_plot.legend(loc='lower left', bbox_to_anchor=(0, 0), borderaxespad=0, frameon=True)




            figure.tight_layout()
            figure.canvas.draw()
            
            # Initialize update_count as a global variable
            global update_count
            update_count = -1
            
            def update_plot(i):
                global update_count  # To modify the update_count variable
                try:
                    if not data_queue.empty():  # Only update if data is available
                        update_count += 1
                        # Adjust the skip interval based on the current size of the queue
                        queue_size = data_queue.qsize()
                        skip_interval = min(30, 10 + max(1, queue_size // 50))
                        # Retrieve updated data
                        raw_waveform_in, raw_waveform_out = data_queue.get_nowait()
                        # Skip plotting unless update_count is a multiple of 5
                        # Skip plotting unless update_count is a multiple of skip_interval
                        if update_count % skip_interval == 0:
                            XY_scatter.set_offsets(np.c_[raw_waveform_in, raw_waveform_out])
                            
                            # Calculate padding
                            x_padding = (max(raw_waveform_in) - min(raw_waveform_in)) * 0.1
                            y_padding = (max(raw_waveform_out) - min(raw_waveform_out)) * 0.1

                            # Set plot limits with padding
                            XY_plot.set_xlim([min(raw_waveform_in) - x_padding, max(raw_waveform_in) + x_padding])
                            XY_plot.set_ylim([min(raw_waveform_out) - y_padding, max(raw_waveform_out) + y_padding])
                    else:
                        # Reset update_count to allow to print the next time there is data.
                        update_count = -1
                except queue.Empty:
                    pass  # If the queue is empty, do nothing

                return XY_scatter,

            ani = FuncAnimation(figure, update_plot, interval=100, blit=False, cache_frame_data=False)
            plt.show(block=True)
            
        except Exception as e:
            print(f"Exception in create_plot for XYoscilloscope: {e}")
            
    def start_plot_process(self):
        self.plot_process = multiprocessing.Process(
            target=self.create_plot,
            args=(self.data_queue, self.command_queue, self.stop_event)
        )
        self.plot_process.start()

    def set_data(self, raw_waveform_in, raw_waveform_out):
        # Ensure data is a NumPy array
        raw_waveform_in = np.array(raw_waveform_in)
        raw_waveform_out = np.array(raw_waveform_out)

        # Send the dataset to the data queue
        if not self.data_queue.full():
            self.data_queue.put((raw_waveform_in, raw_waveform_out))

        # Add a command to trigger the line plot
        self.command_queue.put("trigger_line_plot")

    def stop_plot_process(self):
        if self.plot_process and self.plot_process.is_alive():
            self.stop_event.set()
            self.plot_process.join()

# --- End class PlotManager XYoscilloscope using FuncAnimation:--- --------------------
"""

# --- class PlotManager XYoscilloscope:--- --------------------------------------------
class XYoscilloscope:
    def __init__(self):

        self.data_queue = multiprocessing.Queue(maxsize=100)  # Limited queue size for control
        self.command_queue = multiprocessing.Queue()  # Queue for commands
        self.stop_event = multiprocessing.Event()  # Event for stopping the process
        self.plot_process = None  # Plotting process

    def create_plot(self, data_queue, command_queue, stop_event):
        try:
            plt.ion()
            figure, (XY_plot) = plt.subplots(1)

            # Set the window title
            figure.canvas.manager.set_window_title("Oscilloscope X-Y Plot")
            
            # Configure XY plot
            XY_plot.set_title('XY Oscilloscope Plot')
            XY_plot.set_xlabel('raw_waveform_in')
            XY_plot.set_ylabel('raw_waveform_out')
            XY_plot.grid(which='both')

            # Initialize scatter plots
            XY_scatter = XY_plot.scatter([0], [0], color='red', label='XY Data', s=20)
            
            # Adding legend in the very left corner, outside the plot area
            XY_plot.legend(loc='lower left', bbox_to_anchor=(0, 0), borderaxespad=0, frameon=True)

            figure.tight_layout()
            
            # Initialize update_count as a global variable
            global update_count
            update_count = -1
            
            # Define the close event handler
            def on_close(event):
                # print("Window is closing...")
                stop_event.set()

            # Bind the close event handler to the figure
            figure.canvas.mpl_connect('close_event', on_close)
            
            while not stop_event.is_set():
                try:
                    if not data_queue.empty():  # Only update if data is available
                        update_count += 1
                        # Adjust the skip interval based on the current size of the queue
                        queue_size = data_queue.qsize()
                        skip_interval = min(30, 10 + max(1, queue_size // 50))
                        # Retrieve updated data
                        raw_waveform_in, raw_waveform_out = data_queue.get_nowait()
                        # Skip plotting unless update_count is a multiple of 5
                        # Skip plotting unless update_count is a multiple of skip_interval
                        if update_count % skip_interval == 0:
                            XY_scatter.set_offsets(np.c_[raw_waveform_in, raw_waveform_out])
                            
                            # Calculate padding
                            x_padding = (max(raw_waveform_in) - min(raw_waveform_in)) * 0.1
                            y_padding = (max(raw_waveform_out) - min(raw_waveform_out)) * 0.1

                            # Set plot limits with padding
                            XY_plot.set_xlim([min(raw_waveform_in) - x_padding, max(raw_waveform_in) + x_padding])
                            XY_plot.set_ylim([min(raw_waveform_out) - y_padding, max(raw_waveform_out) + y_padding])
                            
                            figure.canvas.draw()
                            figure.canvas.flush_events()
                    else:
                        # Reset update_count to allow to print the next time there is data.
                        update_count = -1
                        # Add a slight delay to prevent excessive CPU usage
                        time.sleep(0.1)

                except queue.Empty:
                    pass  # If the queue is empty, do nothing
                    time.sleep(0.1)  # Sleep if the queue is empty
                # Allow GUI events to process
                figure.canvas.flush_events()
                time.sleep(0.1)

        except Exception as e:
            print(f"Exception in create_plot for XYoscilloscope: {e}")
            
    def start_plot_process(self):
        self.plot_process = multiprocessing.Process(
            target=self.create_plot,
            args=(self.data_queue, self.command_queue, self.stop_event)
        )
        self.plot_process.start()

    def set_data(self, raw_waveform_in, raw_waveform_out):
        # Ensure data is a NumPy array
        raw_waveform_in = np.array(raw_waveform_in)
        raw_waveform_out = np.array(raw_waveform_out)

        # Send the dataset to the data queue
        if not self.data_queue.full():
            self.data_queue.put((raw_waveform_in, raw_waveform_out))

        # Add a command to trigger the line plot
        self.command_queue.put("trigger_line_plot")

    def stop_plot_process(self):
        if self.plot_process and self.plot_process.is_alive():
            self.stop_event.set()
            self.plot_process.join()

# --- End class PlotManager XYoscilloscope:--- ----------------------------------------

# --- Start class PhaseDiffProcessor:--- ----------------------------------------------

class PhaseDiffProcessor:
    def __init__(self, initial_evaluation_count, threshold, weights):
        # Initialize the processor with given parameters
        self.initial_evaluation_count = initial_evaluation_count
        self.threshold = threshold
        self.weights = weights
        self.phase_diff_values = deque(maxlen=len(weights))
        self.sign_locked = False
        self.locked_sign = None
    
    def weighted_moving_average(self, data):
        # Calculate the weighted moving average of the data
        return np.dot(data, self.weights[:len(data)]) / sum(self.weights[:len(data)])
    
    def majority_sign(self):
        # Determine the majority sign of the phase difference values
        return 1 if sum(np.sign(self.phase_diff_values)) >= 0 else -1
    
    def process_phase_diff(self, value):
        # Append the incoming value to phase_diff_values
        self.phase_diff_values.append(value)
        
        # If initial evaluation period is not complete, wait
        if len(self.phase_diff_values) < self.initial_evaluation_count:
            # print(f"Value: {value}, Waiting for initial evaluation period to complete.")
            return None
        
        # After initial evaluation period, lock the sign based on majority sign
        if len(self.phase_diff_values) == self.initial_evaluation_count and not self.sign_locked:
            self.locked_sign = self.majority_sign()
            self.sign_locked = True
            # print(f"Initial evaluation completed. Locked Sign: {'Positive' if self.locked_sign > 0 else 'Negative'}")
        
        # If sign is locked, check if the incoming value is outside the threshold
        if self.sign_locked:
            if abs(value) > self.threshold:
                # print(f"Value: {value}, Outside threshold. Locked Sign remains: {'Positive' if self.locked_sign > 0 else 'Negative'}")
                return self.locked_sign
            else:
                moving_avg = self.weighted_moving_average(self.phase_diff_values)
                self.locked_sign = 1 if moving_avg >= 0 else -1
                # print(f"Value: {value}, Within threshold. Updated Locked Sign: {'Positive' if self.locked_sign > 0 else 'Negative'}")

        return self.locked_sign
    
    def cleanup(self):
        # Reset the processor attributes to their initial state
        self.phase_diff_values.clear()
        self.sign_locked = False
        self.locked_sign = None
        # print("Processor cleaned up and ready for new data.")

# --- End class PhaseDiffProcessor:--- -----------------------------------------------
            
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
plot_window_height = main_window_height / 4 - win_vertical_border
setting_window_height = main_window_height - plot_window_height
setting_window_width = 500
plot_window_width = main_window_width - setting_window_width - win_horizontal_border
scrollbar_width = 15  # Adjust this value as needed for your scrollbar width
# Increase the viewport width to accommodate the scrollbar
viewport_width = main_window_width + scrollbar_width
# List of available themes
win_theme = ['Dark', 'Light']

# Global variables to store themes
dark_theme = None
light_theme = None

# Initialize graphics settings
def __init__(self):
    return

# --- End Graphics Settings --- -------------------------------------------------------

def print_to_terminal(data):
    """
    Add text data to the terminal and auto-scroll to the end.
    """
    dpg.add_text(data, parent='terminal')
    dpg.set_y_scroll(item='terminal', value=-1)  # Auto-scroll to the end of the window
    dpg.focus_item(item='terminal')

def scroll_data_table():
    """
    Auto-scroll to the end of the data table and focus on it.
    """
    dpg.set_y_scroll(item='dataTableWindow', value=-1)  # Auto-scroll to the end of the window
    dpg.focus_item(item='dataTableWindow')

def oscilloscope_query(cmd):
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
                print("Failed to get value from oscilloscope")
            else:
                return result.tobytes().decode('utf-8').strip()
        except Exception as e:
            print(f"Exception in getting value from oscilloscope: {e}")
            oscilloscope.reset()

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

def AWG_set_frequency(frequency):
    """
    Set the Arbitrary Waveform Generator (AWG) frequency in Hz (0.1Hz to 25MHz for SINE waveform).
    
    Args:
        frequency (float): The desired frequency to set.
        
    Returns:
        float: The measured frequency after adjustment.
    """
    global nWaveSamples
    
    # Initialize variables
    Chkpkpk = 2000
    ChkFREQ = 0.0
    FreqAdjust = 0
    
    # Loop to adjust frequency and minimize transient spikes
    while (Chkpkpk > 1000 or int(ChkFREQ / 10 + 0.5) != int(frequency / 10 + 0.5)) and (int(nWaveSamples) > FreqAdjust):
        FreqAdjust += 1
        oscilloscope_OUT.write(':FUNCtion:FREQuency {}'.format(frequency))
        
        # Allow time for frequency change to settle
        time.sleep(2 * sample_delay_s)

        # Pre-measurement to adjust ranges
        current_v_scale = vertical_scale_to_float(oscilloscope_query(':{}:SCALe?'.format(channel_out)))
        current_h_scale = horizontal_scale_to_float(oscilloscope_query(':HORIzontal:SCALe?'))
        
        # Compute the time range for the waveform acquisition
        raw_frequencies_x = np.linspace(-(horizontal_decades_n * current_h_scale), 
                                        (horizontal_decades_n * current_h_scale), 
                                        300)  # Constant x size of 300 points
        
        # Get waveform data from the output channel
        list_of_out_arrays = []
        for _ in range(nWaveSamples):
            raw_waveform_out = get_waveform(channel_out, current_v_scale)
            list_of_out_arrays.append(raw_waveform_out)

        # Average the collected waveforms
        raw_waveform_out = average_sinusoidal_arrays(list_of_out_arrays)
        Vpkpk_from_curve = peak_to_peak(raw_waveform_out)
        Vpkpk_measured = vertical_scale_to_float(get_pkpk_voltage(channel_out))
        
        # Determine the peak-to-peak voltage
        if Vpkpk_measured == np.nan:
            Vpkpk = Vpkpk_from_curve if Vpkpk_from_curve != np.nan else 0
        elif Vpkpk_from_curve == np.nan:
            Vpkpk = Vpkpk_measured
        else:
            Vpkpk = (Vpkpk_from_curve + Vpkpk_measured) / 2

        # Adjust vertical scale
        if Vpkpk < 0.05:
            closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (3 * Vpkpk * vertical_scaling_factor))))
        elif Vpkpk < 0.1:
            closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (2.5 * Vpkpk * vertical_scaling_factor))))
        elif Vpkpk < 0.2:
            closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (2 * Vpkpk * vertical_scaling_factor))))
        elif Vpkpk < 0.3:
            closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (1.5 * Vpkpk * vertical_scaling_factor))))
        else:
            closest_v_scale_index = amplitude_scales_values.index(min(amplitude_scales_values, key=lambda x: abs(x - (Vpkpk * vertical_scaling_factor))))

        set_amplitude_scale(channel_out, amplitude_scales_commands[closest_v_scale_index])

        # Adjust horizontal scale
        closest_h_scale_index = time_bases_values.index(min(time_bases_values, key=lambda x: abs(x - ((1 / frequency) * horizontal_scaling_factor))))
        set_time_base(time_bases_commands[closest_h_scale_index])

        # Measure frequency
        ChkFREQ = FREQ_out_to_float(oscilloscope_query(':MEASurement:{CH}:FREQuency?'.format(CH=channel_in)))
        Chkpkpk = Vpkpk_from_curve

    # Return the final measured frequency
    return ChkFREQ
def set_time_base(period):
    """
    Set the oscilloscope time-base.
    
    Args:
        period (str): The desired time period to set on the oscilloscope.
    """
    oscilloscope_OUT.write(':HORizontal:SCALe {}'.format(period))

def set_amplitude_scale(channel, scale):
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
    Convert the vertical scale voltage reading to a float.
    
    Args:
        voltage (str): The voltage reading as a string.
        
    Returns:
        float: The converted voltage as a float, scaled if necessary.
    """
    if 'mV' in voltage:
        match = re.search(r"[-+]?\d*\.?\d+", voltage)
        if match:
            return float(match.group(0)) / 1E3
        return None
    else:
        match = re.search(r"[-+]?\d*\.?\d+", voltage)
        if match:
            return float(match.group(0))
        return None

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

def setup_oscilloscope():
    """
    Configure the OWON HDS320S Handheld Oscilloscope and Arbitrary Waveform Generator (AWG).

    This function performs general device configuration, sets up channels, and initializes the AWG.
    """
    global LogFile, JSONLogFile

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

    # Set the trigger to rising edge, VERY IMPORTANT!

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
    print_to_terminal('Now adjust both the vertical and horizontal scales before performing any measurement...')

    # Turn on the device at the correct initial range
    oscilloscope_write(':{CH}:DISPlay ON'.format(CH='CH1'))
    oscilloscope_write(':{CH}:DISPlay ON'.format(CH='CH2'))
    oscilloscope_write(':CHANnel ON')
    AWG_set_frequency(decades_list[start_decade])
    MeasurementFREQ = oscilloscope_query(':MEASurement:{CH}:FREQuency?'.format(CH=channel_in))
    print_to_terminal('Channel out FREQuency: ' + MeasurementFREQ + " Calculate Value: " + str(round(FREQ_out_to_float(MeasurementFREQ), 0)) + '\n')

    if round(FREQ_out_to_float(MeasurementFREQ), 0) < 0.1:
        print_to_terminal('Issue with calculated channel out FREQuency.')
        print_to_terminal('Recommend pressing Auto button on the HDS320S.')
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
    if is_recording == False:
        is_recording = True
        threading.Thread(target=start_mesurement_threaded, daemon=True).start()
        # Enable all the control items
        control_items = [
            'Stop'
        ]
        for item in control_items:
            dpg.configure_item(item=item, enabled=True)
        dpg.bind_item_theme('Stop', 'RedButton')
        
def start_mesurement_threaded():
    """
    Start the measurement process for the oscilloscope and Arbitrary Waveform Generator (AWG).
    
    This function handles the UI setup, data clearing, and measurement initiation.
    """
    global LogFile, JSONLogFile, nWaveSamples, FTTcorection, is_recording
    
    OKrun = True

    initial_evaluation_count = 5
    threshold = 60
    weights = np.array([0.1, 0.15, 0.2, 0.25, 0.3])

    processor = PhaseDiffProcessor(initial_evaluation_count, threshold, weights)
                                   
    # Programmatically select the "Data Table" tab
    select_tab('data_table_tab')
    
    # Clear all the plots
    try:
        # dpg.delete_item(item='MAG_SCATTER')
        dpg.delete_item(item='MAG_LINE')
        # dpg.delete_item(item='PHASE_SCATTER')
        dpg.delete_item(item='PHASE_LINE')
    except Exception as e:
        print(f"Exception clearing plots: {e}")
    
    # Disable all the control items
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'START_MEASURE', 'JSONLOGFILEtag'
    ]
    if platform.system() != 'Windows':
        control_items.append('PLOT_WIN_SETTING')
    for item in control_items:
        dpg.configure_item(item=item, enabled=False)
    dpg.bind_item_theme('START_MEASURE', 'DisabledButton')
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
        with dpg.table(parent='dataTableWindow', header_row=True, resizable=True, policy=dpg.mvTable_SizingStretchProp, width=setting_window_width - 10, height=350, freeze_rows=1,
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
 
        if plot_win_disposition == 'MatPlotLib':
           # 2025/02/20 Fix code that was unstable and causing crashes, by having Matplotlib plot opened in a separate process          
           # If the operating system is Mint Linux might be okay to run
           # Check the operating system is not Windows.
               # Indent Mint Linux  commands here      
           if platform.system() != 'Windows':
               plot_manager = PlotManagerMagnitudePhase()
               plot_manager.start_plot_process(decades_list[start_decade], decades_list[stop_decade], points_per_decade)   

               plot_manager_fft_oscilloscope = PlotManagerFFToscilloscope()
               plot_manager_fft_oscilloscope.start_plot_process(decades_list[start_decade], decades_list[stop_decade], points_per_decade)
                   
               # Creating XY Oscilloscope plot
               xy_oscilloscope = XYoscilloscope()
               xy_oscilloscope.start_plot_process()
           
        # Measurement loop
        for index, frequency in enumerate(raw_frequencies_range):
            if is_recording == False:
                break
            # Set the current test frequency
            MeasurementFREQ = AWG_set_frequency(frequency)
 
            # Ask for the current vertical and horizontal scale (only the first time)
            current_v_scale = vertical_scale_to_float(oscilloscope_query(':{}:SCALe?'.format(channel_out)))
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
                # print("READin: ", peak_to_peak(raw_waveform_in))
            raw_waveform_in = average_sinusoidal_arrays(list_of_in_arrays)
            # print("READin AVG: ", peak_to_peak(raw_waveform_in))
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
                    raw_waveform_out = get_waveform(channel_out, current_v_scale)
                    readPk = peak_to_peak(raw_waveform_out)
                    attempts += 1
                list_of_out_arrays.append(raw_waveform_out)
                # print("READout: ", peak_to_peak(raw_waveform_out))
            raw_waveform_out = average_sinusoidal_arrays(list_of_out_arrays)
            # print("READout AVG: ", peak_to_peak(raw_waveform_out))
            raw_frequencies_x_out = []
            for timestep in range(len(raw_waveform_out)):
                raw_frequencies_x_out.append(timestep * sampling_PERIOD)

            FSA = full_sin_analysis(sampling_PERIOD, FTTcorection, raw_waveform_in, raw_waveform_out)
            FFT_freq = FSA["fft_frequencies"][:len(FSA["fft_frequencies"]) // 2]
            FFT_result = FSA["amplitude_db_signal2"][:len(FSA["fft_frequencies"]) // 2]
            FFT_max = FSA['max_amplitude_db_signal2']
            FFT_maxfreq = FSA['angular_frequency_signal2'] / (2 * np.pi)
  
            Phase_Diff_degree = FSA['phase_difference_degrees']
            
            sign = processor.process_phase_diff(Phase_Diff_degree)
            if sign is not None:
                if sign > 0:
                    Phase_Diff_degree = abs(Phase_Diff_degree)
                else:
                    Phase_Diff_degree = -abs(Phase_Diff_degree)
            
            # --- Update after data read from oscilloscope ---
            OSCvRange = 0.5
            OSCvRange1 = 0.5
            OSCvRange2 = 0.5
            DMrangeV1 = oscilloscope_query(':CH1:SCALe?')
            DMrangeV2 = oscilloscope_query(':CH2:SCALe?')
            for Voffset in range(len(amplitude_scales_commands)):
                if DMrangeV1 == amplitude_scales_commands[Voffset]:
                    OSCvRange1 = amplitude_scales_values[Voffset]
                if DMrangeV2 == amplitude_scales_commands[Voffset]:
                    OSCvRange2 = amplitude_scales_values[Voffset]
            OSCvRange = max(OSCvRange1, OSCvRange2)
             
            graph_processing(OSCvRange,
                             FFT_maxfreq,
                             FFT_max,
                             FFT_freq,
                             FFT_result,
                             raw_frequencies_x_in,
                             raw_waveform_in,
                             raw_frequencies_x_out,
                             raw_waveform_out)                   
            
            VpkpkMeter = vertical_scale_to_float(get_pkpk_voltage(channel_out))
            VpkpkMeterIn = vertical_scale_to_float(get_pkpk_voltage(channel_in))
            if VpkpkMeterIn > 0:          
                gY = (20 * np.log10(VpkpkMeter / VpkpkMeterIn)) if VpkpkMeter > 0 else float('-inf')
            else:
                gY = (20 * np.log10(VpkpkMeter / waveform_amplitude_V)) if VpkpkMeter > 0 else float('-inf')
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
                            "VpkpkMeter": VpkpkMeter,
                            "gX": gX,
                            "gY": gY,
                            "pX": pX,
                            "pY": pY
            }
            
            if plot_win_disposition == 'MatPlotLib':
                # If the operating system is Mint Linux might be okay to run
                # Check the operating system is not Windows.
                if platform.system() != 'Windows':
                    # plot_manager.set_data(gain_X, gain_Y, phase_X, phase_Y)
                    plot_manager.set_data_point(data[frequency]['gX'],
                                            data[frequency]['gY'],
                                            data[frequency]['pX'],
                                            data[frequency]['pY'])      
                    

                    # Indent Mint Linux  commands here
                    plot_manager_fft_oscilloscope.set_data_point(data[frequency]['FFT_maxfreq'],
                                                                data[frequency]['FFT_max'])
                    
                    plot_manager_fft_oscilloscope.set_full_data(data[frequency]['FFT_freq'],
                                                                data[frequency]['FFT_result'],
                                                                data[frequency]['raw_frequencies_x_in'],
                                                                data[frequency]['raw_waveform_in'],
                                                                data[frequency]['raw_frequencies_x_out'],
                                                                data[frequency]['raw_waveform_out'])
                            
                    raw_waveform_in = data[frequency]['raw_waveform_in']
                    raw_waveform_out = data[frequency]['raw_waveform_out']
                    xy_oscilloscope.set_data(raw_waveform_in, raw_waveform_out)
        
            with dpg.table_row(parent='DataTable'):
                dpg.add_text(str(int(round(pX, 0))))
                dpg.add_text(str(int(round(gX, 0))))
                dpg.add_text(str(round(VpkpkMeter, 3)))
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
            f.write(str(round(VpkpkMeter, 3)))
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
       
    # Re-enable all the control elements
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'JSONLOGFILEtag'
    ]
    if platform.system() != 'Windows':
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
    dpg.bind_item_theme('Stop', 'DisabledButton')    
    # Programmatically select the "Data Table" tab
    select_tab('data_table_tab')    
    # Cleanup after processing
    processor.cleanup()
    with open(JSONLogFile, 'w') as json_file:
        json.dump(data, json_file, indent=1)
    f.close()
    JSONLogFile = 'AUTO.json'
    LogFile = 'AUTO.csv'
    LogFileGUIlabel = "CSV Log File: " + LogFile
    JSONLogFileGUIlabel = "JSON Log File: " + JSONLogFile
    # Start post-processing if run was successful
    if OKrun:
        if plot_win_disposition == 'MatPlotLib':
            plot_manager.set_full_data(gain_X, gain_Y, phase_X, phase_Y)
        post_processing(gain_X, gain_Y, phase_X, phase_Y)
   
def PlayBack():
    global JSONLogFile
    global is_playing, is_paused, play_speed
    # Disable all the control items
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'START_MEASURE', 'JSONLOGFILEtag'
    ]
    if platform.system() != 'Windows':
        control_items.append('PLOT_WIN_SETTING')
    for item in control_items:
        dpg.configure_item(item=item, enabled=False)
    dpg.bind_item_theme('JSONLOGFILEtag', 'DisabledButton')
    dpg.bind_item_theme('SEARCH_OSCILLOSCOPE', 'DisabledButton')  
    # Enable all the control items
    control_items = [
        'Stop'
    ]
    for item in control_items:
        dpg.configure_item(item=item, enabled=True)
    dpg.bind_item_theme('Stop', 'RedButton') 
    # Clear all the plots
    try:
        # dpg.delete_item(item='MAG_SCATTER')
        dpg.delete_item(item='MAG_LINE')
        # dpg.delete_item(item='PHASE_SCATTER')
        dpg.delete_item(item='PHASE_LINE')
    except Exception as e:
        print(f"Exception clearing plots: {e}")
    
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
       
    # Open output JSON file for data read.
    with open(JSONLogFile, 'r') as json_file:
        data = json.load(json_file)
            
    range = list(data.keys())
    start_decade = decades_list.index(find_closest(decades_list, range[0]))
    dpg.set_value('START_DEC', decades_list_string[start_decade])
    stop_decade = decades_list.index(find_closest(decades_list, range[-1]))
    dpg.set_value('STOP_DEC', decades_list_string[stop_decade])
    # print("Start: ", decades_list[start_decade], " End: ", decades_list[stop_decade])
    # print("Start: ", decades_list_string[start_decade], " End: ", decades_list_string[stop_decade])
    
    # Set Magnitude/Phase/FFT plot Max range
    dpg.set_axis_limits('MAG_X', decades_list[start_decade], decades_list[stop_decade])
    dpg.set_axis_limits('PHASE_X', decades_list[start_decade], decades_list[stop_decade])
    dpg.set_axis_limits('FFT_X', decades_list[start_decade], decades_list[stop_decade])

    # Cleanup old data from table in GUI
    dpg.delete_item(item='DataTable')
    with dpg.table(parent='dataTableWindow', header_row=True, resizable=True, policy=dpg.mvTable_SizingStretchProp, width=setting_window_width - 10, height=350, freeze_rows=1,
                    scrollY=True, scrollX=False, borders_outerH=True, borders_innerV=True, borders_innerH=True, borders_outerV=True, tag='DataTable'):
        dpg.add_table_column(label="Frequency")
        dpg.add_table_column(label="MeasFreq")
        dpg.add_table_column(label="Voltage")
        dpg.add_table_column(label="Phase")
        dpg.add_table_column(label="FFT Hz")
        dpg.add_table_column(label="FFT Meg")

    # Programmatically select the "Data Table" tab
    select_tab('data_table_tab')
    
    if plot_win_disposition == 'MatPlotLib':
        # 2025/02/20 Fix code that was unstable and causing crashes, by having Matplotlib plot opened in a separate process
        # If the operating system is Mint Linux might be okay to run
        # Check the operating system is not Windows.
        if platform.system() != 'Windows':
            # Indent Mint Linux  commands here    
            plot_manager = PlotManagerMagnitudePhase()
            plot_manager.start_plot_process(decades_list[start_decade], decades_list[stop_decade], points_per_decade)
          
            plot_manager_fft_oscilloscope = PlotManagerFFToscilloscope()
            plot_manager_fft_oscilloscope.start_plot_process(decades_list[start_decade], decades_list[stop_decade], points_per_decade)
        
            # Creating XY Oscilloscope plot
            xy_oscilloscope = XYoscilloscope()
            xy_oscilloscope.start_plot_process()
         
    for index, frequency in enumerate(list(data.keys())):
        # print(frequency, " : ", data[frequency]['FFT_maxfreq'])
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
                              
        if plot_win_disposition == 'MatPlotLib':
            # If the operating system is Mint Linux might be okay to run
            # Check the operating system is not Windows.
            if platform.system() != 'Windows':
                # plot_manager.set_data(gain_X, gain_Y, phase_X, phase_Y)
                plot_manager.set_data_point(data[frequency]['gX'],
                                            data[frequency]['gY'],
                                            data[frequency]['pX'],
                                            data[frequency]['pY'])

                # Indent Mint Linux  commands here        
                plot_manager_fft_oscilloscope.set_data_point(data[frequency]['FFT_maxfreq'],
                                                            data[frequency]['FFT_max'])
                
                plot_manager_fft_oscilloscope.set_full_data(data[frequency]['FFT_freq'],
                                                            data[frequency]['FFT_result'],
                                                            data[frequency]['raw_frequencies_x_in'],
                                                            data[frequency]['raw_waveform_in'],
                                                            data[frequency]['raw_frequencies_x_out'],
                                                            data[frequency]['raw_waveform_out'])
                
                raw_waveform_in = data[frequency]['raw_waveform_in']
                raw_waveform_out = data[frequency]['raw_waveform_out']
                xy_oscilloscope.set_data(raw_waveform_in, raw_waveform_out)
        
        with dpg.table_row(parent='DataTable'):
            dpg.add_text(str(int(round(data[frequency]['pX'], 0))))
            dpg.add_text(str(int(round(data[frequency]['gX'], 0))))
            dpg.add_text(str(round(data[frequency]['VpkpkMeter'], 3)))
            dpg.add_text(str(round(data[frequency]['pY'], 3)))
            dpg.add_text(str(int(round(data[frequency]['FFT_maxfreq'], 0))))
            dpg.add_text(str(round(data[frequency]['FFT_max'], 2)))
            
            # Scroll to the bottom
            num_rows = len(dpg.get_item_children('DataTable', slot=1))
            scroll_pos = num_rows * 25  # Adjust the multiplier if row height is different
            dpg.set_y_scroll('DataTable', scroll_pos)
        time.sleep(play_speed)
        pause_loop()
        if is_playing == False and is_paused == False:
            stop_play_cleaup()
            return
    if plot_win_disposition == 'MatPlotLib':
        plot_manager.set_full_data(gain_X, gain_Y, phase_X, phase_Y)
    post_processing(gain_X, gain_Y, phase_X, phase_Y)
    stop_play_cleaup()
    is_playing = False
    dpg.configure_item('Play', label='Play')
    
def stop_play_cleaup():
    # Re-enable all the control elements
    control_items = [
        'CH_IN', 'CH_OUT', 'CH1_ATTENUATION_RATIO', 'CH2_ATTENUATION_RATIO',
        'CH1_COUPLING_MODE', 'CH2_COUPLING_MODE', 'SAMPL_MODE', 'DEPMEM',
        'AWG_OUT_VOLTAGE', 'HIGH_Z', 'POINTS_X_DEC', 'POINTS_SPACING',
        'START_DEC', 'STOP_DEC', 'POINT_SCALE_COEFF',
        'V_SCALE_COEFF', 'H_SCALE_COEFF', 'OSCILL_TIMEOUT', 'CODE_EXEC_PAUSE',
        'WIN_THEME', 'SEARCH_OSCILLOSCOPE', 'JSONLOGFILEtag'
        ]
    if platform.system() != 'Windows':
        control_items.append('PLOT_WIN_SETTING')
    for item in control_items:
        dpg.configure_item(item=item, enabled=True)
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
        # print("The variable is a NumPy array.")
        list_variable = variable.tolist()
        # print("NumPy array converted to list.")
        return list_variable
    elif isinstance(variable, list):
         # print("The variable is a list.")
         return variable
    else:
        # print("The variable is not a NumPy array or a list.")
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
def full_sin_analysis(sample_rate, analysis_factor, signal1, signal2):
    """
    Performs a full sinusoidal analysis on two signals.

    Args:
        sample_rate (float): The sample rate of the signals.
        analysis_factor (float): The factor to be added to the normalized amplitude.
        signal1 (list): The first signal.
        signal2 (list): The second signal.

    Returns:
        dict: A dictionary containing the analysis results.
    """
    
    n = len(signal1)
    duration = n / sample_rate
    time = np.linspace(0, duration, n, endpoint=False)
    
    fft_signal1 = np.fft.fft(signal1)
    fft_signal2 = np.fft.fft(signal2)
    frequencies = np.fft.fftfreq(n, 1/sample_rate)
    
    # Normalized amplitude in dB
    amplitude_db_signal1 = analysis_factor + 20 * np.log10((np.abs(fft_signal1) / n) + 1e-10)  # Add small value 1e-10 (-200 dB) to avoid log(0)
    amplitude_db_signal2 = analysis_factor + 20 * np.log10((np.abs(fft_signal2) / n) + 1e-10)  # Add small value 1e-10 (-200 dB) to avoid log(0)

    # Max amplitude and frequency
    max_index_signal1 = np.argmax(np.abs(fft_signal1[0:n//2]))
    max_freq_signal1 = frequencies[max_index_signal1]
    max_amplitude_db_signal1 = amplitude_db_signal1[max_index_signal1]
    
    max_index_signal2 = np.argmax(np.abs(fft_signal2[0:n//2]))
    max_freq_signal2 = frequencies[max_index_signal2]
    max_amplitude_db_signal2 = amplitude_db_signal2[max_index_signal2]

    # Angular frequency
    angular_frequency_signal1 = 2 * np.pi * max_freq_signal1
    angular_frequency_signal2 = 2 * np.pi * max_freq_signal2
    
    # Phase and phase difference
    phase_signal1 = np.angle(fft_signal1[max_index_signal1], deg=True)
    phase_signal2 = np.angle(fft_signal2[max_index_signal2], deg=True)
    phase_difference = phase_signal2 - phase_signal1

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
        "gain_db": gain_db,
        "fft_real_signal1": np.real(fft_signal1),
        "fft_imag_signal1": np.imag(fft_signal1),
        "fft_real_signal2": np.real(fft_signal2),
        "fft_imag_signal2": np.imag(fft_signal2)
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
        print("FTT Plot data issue seen MaxFFTG V at ", len(MxFFTy))

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
        print("FTT Plot data issue seen at ", len(FFTx))
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
        RangeOffset = 2
        dpg.set_axis_limits('OSC_Y', -RangeOffset * OSCvRange, RangeOffset * OSCvRange)
        dpg.set_axis_limits_auto(axis='OSC_X')
        dpg.fit_axis_data(axis='OSC_X')
    
    if len(OSCyout) == len(OSCyin) and len(OSCyin) > 0:
        dpg.set_value('XY_SCATTEROUT', [list(OSCyin[-nsamples:]), list(OSCyout[-nsamples:])])
        dpg.set_axis_limits('XY_Y', -RangeOffset * OSCvRange, RangeOffset * OSCvRange)
        dpg.set_axis_limits('XY_X', -RangeOffset * OSCvRange, RangeOffset * OSCvRange)
    else:
        print("XY Plot data issue seen at ", len(FFTx))

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
    magnitude_spline = UnivariateSpline(gain_X, gain_Y, k=3, s=0)
    dpg.add_line_series(x=gain_X, y=magnitude_spline(gain_X), parent='MAG_Y', tag='MAG_LINE', label='Interpolated waveform')
    dpg.configure_item(item='DRAG_LINE_MAG_X', default_value=np.median(gain_X))
    dpg.configure_item(item='DRAG_LINE_MAG_Y', default_value=np.median(gain_Y))
    dpg.set_axis_limits_auto(axis='MAG_Y')
    dpg.fit_axis_data(axis='MAG_Y')
    dpg.set_axis_limits_auto(axis='MAG_X')
    dpg.fit_axis_data(axis='MAG_X')

    # Ensure phase_X is sorted and strictly increasing
    sorted_indices = np.argsort(phase_X)
    phase_X = np.array(phase_X)[sorted_indices]
    phase_Y = np.array(phase_Y)[sorted_indices]

    # Remove duplicates in phase_X
    unique_indices = np.diff(phase_X) > 0
    phase_X = phase_X[np.concatenate(([True], unique_indices))]
    phase_Y = phase_Y[np.concatenate(([True], unique_indices))]

    # Outlier detection based on local deviation
    window_size = 5  # Number of neighbors to consider on each side
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
    if len(filtered_phase_X) >= 3:  # Ensure enough points for spline
        phase_spline = UnivariateSpline(filtered_phase_X, filtered_phase_Y, k=3, s=0)
        dpg.add_line_series(x=filtered_phase_X, y=phase_spline(filtered_phase_X), parent='PHASE_Y', tag='PHASE_LINE', label='Interpolated waveform')
        dpg.configure_item(item='DRAG_LINE_PHASE_X', default_value=np.median(filtered_phase_X))
        dpg.configure_item(item='DRAG_LINE_PHASE_Y', default_value=np.median(filtered_phase_Y))
        dpg.set_axis_limits_auto(axis='PHASE_X')
        dpg.fit_axis_data(axis='PHASE_X')


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

        # Theme for a disabled button
        with dpg.theme_component(dpg.mvButton, enabled_state=False):
            dpg.add_theme_color(dpg.mvThemeCol_Text, [150, 150, 150])

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

        # Disabled button theme
        with dpg.theme(tag='DisabledButton') as disabled_button_theme:
            with dpg.theme_component(dpg.mvButton, enabled_state=False):
                dpg.add_theme_color(dpg.mvThemeCol_Button, (128, 128, 128, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (128, 128, 128, 255))
                dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (128, 128, 128, 255))
                dpg.add_theme_color(dpg.mvThemeCol_Text, (255, 255, 255, 255))
       
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

        
# -- gui settings --- ---------------------------------------------
def main():
    global channel_in, channel_out, CH1_probe_attenuation_ratio, CH2_probe_attenuation_ratio
    global CH1_coupling, CH2_coupling, Sample_command, DEPMEM, waveform_amplitude_V
    global AWG_output_impedance, points_per_decade, start_decade, stop_decade
    global point_resize_factor, vertical_scaling_factor, horizontal_scaling_factor
    global nWaveSamples, FTTcorection, read_delay_ms, sample_delay_s, plot_win_disposition
    global is_playing, is_paused, play_speed, is_recording

    dpg.create_context()

    # Create themes
    create_themes()

    # Create a file dialog for selecting files
    with dpg.file_dialog(directory_selector=False, show=False, callback=callbackCSV, id="file_dialog_id1", width=700, height=400):
        dpg.add_file_extension(".csv", color=(106, 127, 235, 255), custom_text="[csv]")
    with dpg.file_dialog(directory_selector=False, show=False, callback=callbackJSON, id="file_dialog_id2", width=700, height=400):
        dpg.add_file_extension(".json", color=(106, 127, 235, 255), custom_text="[json]")        

    with dpg.window(tag='main_window', pos=(0, 0), width=main_window_width - scrollbar_width, no_bring_to_front_on_focus=True, autosize=True, height=main_window_height, no_title_bar=True, no_move=True, no_resize=False):
        with dpg.child_window(tag='scrolling_container', width=main_window_width, height=main_window_height, horizontal_scrollbar=False, border=False):
            # Create a Control window with specified settings
            # with dpg.window(tag='controlwin', width=setting_window_width, height=setting_window_height - win_vertical_border, no_resize=True, pos=(0, win_vertical_border), no_move=True, no_close=True, no_collapse=True, no_title_bar=True):
            with dpg.child_window(tag='controlwin', width=setting_window_width, height=setting_window_height - win_vertical_border, pos=(0, win_vertical_border), border=False, menubar=False): 
                dpg.add_text("Oscilloscope settings:")
                # Add various combo boxes for different settings
                dpg.add_combo(tag='CH_IN', items=available_channels, label='Input Channel', default_value=available_channels[0], width=items_standard_width)
                dpg.add_combo(tag='CH_OUT', items=available_channels, label='Output Channel', default_value=available_channels[1], width=items_standard_width)
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
                    if platform.system() != 'Windows':
                        dpg.add_radio_button(tag='PLOT_WIN_SETTING', label='Plot Engine', items=plot_win_settings, default_value=plot_win_settings[0])
                    dpg.add_radio_button(tag='DRAG_LINES_COMBO', items=drag_line_values, default_value=drag_line_values[1], callback=set_drag_lines)

                # Add advanced settings
                dpg.add_text('\nAdvanced settings:')
                dpg.add_input_float(tag='POINT_SCALE_COEFF', label='Point scale coefficient', min_value=0, min_clamped=True, default_value=5850, width=items_standard_width)
                dpg.add_input_float(tag='V_SCALE_COEFF', label='Vertical scale calibration coeff.', min_value=0, min_clamped=True, default_value=0.33, width=items_standard_width)
                dpg.add_input_float(tag='H_SCALE_COEFF', label='Horizontal scale calibration coeff.', min_value=0, min_clamped=True, default_value=0.80, width=items_standard_width)
                dpg.add_input_float(tag='OSCILL_TIMEOUT', label='Oscilloscope reading timeout (ms)', min_value=0, min_clamped=True, default_value=1200, width=items_standard_width)
                dpg.add_input_float(tag='CODE_EXEC_PAUSE', label='Commands execution delay (s)', min_value=0, min_clamped=True, default_value=0.50, width=items_standard_width)

                # Add buttons for various actions
                dpg.add_input_int(tag='SAMPLE_WAVES', label='Number of times waves are sampled', min_value=1, min_clamped=True, default_value=2, width=items_standard_width)
                dpg.add_input_int(tag='FFToffset', label='FFT correction offset (db)', default_value=2, width=items_standard_width)
                dpg.add_text('\n')
                with dpg.group(horizontal=True):
                    dpg.add_button(tag='SEARCH_OSCILLOSCOPE', label='Search and Setup DSO', callback=search_oscilloscope)
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
                        with dpg.child_window(tag='terminal', width=setting_window_width, height=int(main_window_height / 4) - 40, label='Terminal', border=False):
                            pass
                    with dpg.tab(label="Data Table", tag='data_table_tab'):
                        with dpg.child_window(tag='dataTableWindow', width=setting_window_width - 10, height=350, border=False, menubar=False):
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
			height=set_window_height,
			resizable=True,
			max_height=main_window_height,
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
        read_delay_ms = int(dpg.get_value(item='OSCILL_TIMEOUT'))
        sample_delay_s = float(dpg.get_value(item='CODE_EXEC_PAUSE'))
        # Plot parameters
        if platform.system() != 'Windows':
            plot_win_disposition = str(dpg.get_value(item='PLOT_WIN_SETTING'))
        else:
            plot_win_disposition = default_value=plot_win_settings[0]
        # Keep updating plot sizes
        dpg.set_item_height(item='MAG_PLOT_GRAPH', height=dpg.get_item_height(item='MAG_PLOT_WIN'))
        dpg.set_item_width(item='MAG_PLOT_GRAPH', width=dpg.get_item_width(item='MAG_PLOT_WIN'))
        dpg.set_item_height(item='PHASE_PLOT_GRAPH', height=dpg.get_item_height(item='PHASE_PLOT_WIN'))
        dpg.set_item_width(item='PHASE_PLOT_GRAPH', width=dpg.get_item_width(item='PHASE_PLOT_WIN'))

        # Limit the amplitude scales based on the chosen probe attenuation value (see datasheet for allowed values)
        if CH1_probe_attenuation_ratio == '1X':
            CH1_amplitude_scales = amplitude_scales_commands[2:]
            CH1_amplitude_scales = amplitude_scales_values[2:]
        elif CH1_probe_attenuation_ratio == '10X':
            CH1_amplitude_scales = amplitude_scales_commands[:9]
            CH1_amplitude_scales = amplitude_scales_values[:9]
        if CH2_probe_attenuation_ratio == '1X':
            CH2_amplitude_scales = amplitude_scales_commands[2:]
            CH2_amplitude_scales = amplitude_scales_values[2:]
        elif CH2_probe_attenuation_ratio == '10X':
            CH2_amplitude_scales = amplitude_scales_commands[:9]
            CH2_amplitude_scales = amplitude_scales_values[:9]

        dpg.render_dearpygui_frame()
    dpg.destroy_context()
# --- End GUI settings --- ----------------------------------------
