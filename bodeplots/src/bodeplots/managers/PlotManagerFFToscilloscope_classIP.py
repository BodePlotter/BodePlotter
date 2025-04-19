"""

© 2025 by Bode Plotter

MIT License

Copyright (c) [2025] [Bode Plotter]
Bode plot implementation for the handheld oscilloscope OWON HDS320S
Highly rewritten and modified between 01/20/2025 to 04/20/2025
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

import socket
import json
import logging
# import matplotlib.pyplot as plt
import numpy as np
import multiprocessing
import threading
from threading import Event
import weakref
from multiprocessing import Process, Event
from scipy.interpolate import UnivariateSpline
import os
import time
import random
import platform
import subprocess
import sys
import atexit

# --- class PlotManager FFT oscilloscope:--- ------------------------------------------
class PlotManagerFFToscilloscope:
    # Class-level registry to hold instances with weak references
    _instances = weakref.WeakSet()
    _instances_lock = threading.Lock()  # For thread-safety on the registry

    def __init__(self):
        """Initialize the PlotManagerFFToscilloscope."""
        # Initialize a logger specific to this class
        self.logger = logging.getLogger(__name__)  # Inherit global logging configuration
        self.logger.info("PlotManagerFFToscilloscope initialized.")
        

        # Initialize attributes
        self.stop_event = Event()  # Event for stopping the process
        self.plot_process = None  # Plotting process
        self.scatter_updated = False  # Flag to indicate scatter plot update
        self.FFTxMax = np.array([])  # Maximum X values of the FFT as NumPy array
        self.FFTmaxVal = np.array([])  # Maximum values of the FFT as NumPy array
        self.host = '127.0.0.1'  # Socket host
        try:
            self.base_port = int(os.getenv("BASE_PORT", 5001))
        except ValueError:
            # Optionally log the error or set a fallback base port
            self.base_port = 5001
        self.port = self.base_port + 1
        self.stop_reading = False  # Flag to stop reading from the socket

        # Register the new instance in a thread-safe manner
        with PlotManagerFFToscilloscope._instances_lock:
            PlotManagerFFToscilloscope._instances.add(self)

        # Log the initialization of the class
        self.logger.info(f"PlotManagerFFToscilloscope initialized on {self.host}:{self.port}.")



    def create_plot(self, stop_event, start_decade, stop_decade, points_per_decade, cache_dir=None):
        current_os = platform.system()
        if current_os == "Windows":
            if cache_dir is None:
                timestamp = int(time.time() * 1000)  # Get current time in milliseconds
                cache_dir = f"C:\\Temp\\matplotlib_cache_FFT_{timestamp}"
                os.makedirs(cache_dir, exist_ok=True)  # Create the directory if it doesn't exist

            # Update environment variable properly
            os.environ["MPLCONFIGDIR"] = cache_dir
            time.sleep(0.5)  # Short wait to ensure directory is ready 

        # Import Matplotlib
        import matplotlib
        matplotlib.use("qtagg")  # qtagg preferred interactive backend
        import matplotlib.pyplot as plt

        """Matplotlib process to receive data via sockets and render FFT/Oscilloscope plots."""
        try:
            self.logger.info("Starting FFT/Oscilloscope Plot process...")

            # Constants for retry logic
            MAX_RETRIES = 5
            INITIAL_RETRY_DELAY = 1  # Start with 1 second delay
            MAX_RETRY_DELAY = 16     # Maximum wait time between retries

            for attempt in range(MAX_RETRIES):
                try:
                    # Attempt to connect to the server
                    self.logger.info(f"Attempt {attempt + 1}: Connecting to socket server...")
                    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    client_socket.settimeout(0.1)  # Non-blocking socket
                    client_socket.connect((self.host, self.port))
                    self.logger.info("Connected to socket server!")
                    break  # Exit the retry loop if successful
                except socket.error as e:
                    self.logger.error(f"Connection attempt {attempt + 1} failed: {e}")

                    # Calculate the retry delay (exponential backoff with jitter)
                    retry_delay = min(INITIAL_RETRY_DELAY * (2 ** attempt), MAX_RETRY_DELAY)
                    jitter = random.uniform(0, 1)  # Add randomness to prevent synchronized retries
                    retry_delay += jitter

                    self.logger.info(f"Retrying in {retry_delay:.2f} seconds...")
                    time.sleep(retry_delay)
            else:
                # All retries failed
                self.logger.error("Failed to connect to socket server after multiple attempts.")

        
            plt.ion()

            # Adjust figsize based on your desired overall window size.
            # For example, if you want the entire window to be ~600 pixels high at 100 DPI:
            figure, (FFT, oscilloscope) = plt.subplots(2, figsize=(12, 6))
            
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
            self.logger.error(f"Error initializing plots: {e}")
            return
                
        # (After creating your figure and before entering the main event loop.)
        manager = plt.get_current_fig_manager()
        if hasattr(manager, 'window'):
            try:
                window = manager.window

                # Desired initial position (these values can be any starting point)
                desired_x, desired_y = 1450, 717

                # Get the current frame geometry (which includes window decorations)
                geom = window.frameGeometry()
                window_width = geom.width()
                window_height = geom.height()

                # Use PyQt to retrieve the primary screen's available geometry
                # (availableGeometry avoids taskbars/docks)
                from PyQt5.QtWidgets import QApplication
                app = QApplication.instance()
                if app is None:
                    app = QApplication(sys.argv)
                screen = app.primaryScreen()
                screen_geom = screen.availableGeometry()
                screen_width = screen_geom.width()
                screen_height = screen_geom.height()

                # Adjust x and y so the window remains within screen boundaries.
                # We ensure the window's right/bottom edges don't exceed the screen.
                x = max(0, min(desired_x, screen_width - window_width))
                y = max(0, min(desired_y, screen_height - window_height))

                window.move(x, y)
            except Exception as e:
                self.logger.error(f"Error setting window geometry: {e}")

        # Define the close event handler
        def on_close(event):
            self.logger.info("Plot window closed by user.")
            stop_event.set()

        figure.canvas.mpl_connect('close_event', on_close)

        buffer = ""  # Buffer for incoming data
        
        last_update_time = 0  # Initialize the timestamp
        UpdateSkipCanvasCounter = 0  # Initialize the counter to 0
        # Flag to control whether to stop reading from the socket
        stop_reading = False
        
        try:
            while not stop_event.is_set():
                current_time = time.time()
                try:
                    # Receive data from the socket
                    if not stop_reading and client_socket and client_socket.fileno() != -1:
                        chunk = client_socket.recv(32768).decode('utf-8')
                        buffer += chunk
                        
                    # Handle "<CLOSE>" signal
                    if not self.stop_reading and "<CLOSE>" in buffer:
                        # Set flag to stop further attempts to read from the socket
                        stop_reading = True
                        self.logger.info("Received CLOSE signal from client. Processing remaining data...")

                        # Clear the buffer to avoid leftover data
                        buffer = ""
                        self.logger.info("Buffer cleared after processing all complete messages.")

                        # Send acknowledgment back to the server
                        try:
                            client_socket.sendall(b"<CLOSED>")
                            self.logger.info("Sent CLOSED acknowledgment to client.")
                        except Exception as e:
                            self.logger.error(f"Error sending acknowledgment: {e}")
                            break
                        # Close the socket gracefully
                        try:
                            if client_socket and client_socket.fileno() != -1:  # Verify socket is valid
                                self.logger.info("Shutting down socket for buffered data transmission...")
                                client_socket.shutdown(socket.SHUT_RDWR)  # Graceful shutdown for write operations
                                time.sleep(0.1)  # Allow pending operations to complete
                            else:
                                self.logger.warning("Socket is invalid or already closed. Skipping shutdown.")
                                break
                        except socket.error as e:
                            self.logger.error(f"Error during socket shutdown: {e}")
                            break
                        finally:
                            try:
                                if client_socket and client_socket.fileno() != -1:  # Verify socket is still valid
                                    client_socket.close()
                                    self.logger.info("Socket closed.")
                                else:
                                    self.logger.warning("Socket already closed. Skipping close operation.")
                                    break
                            except Exception as e:
                                self.logger.error(f"Error during socket close: {e}")
                                break
                            finally:
                                client_socket = None
                                self.logger.info("Client socket set to None to mark it as unavailable.")

                    if "<END>" in buffer:
                        raw_data, buffer = buffer.split("<END>", 1)

                        try:
                            # Parse the received JSON data
                            data = json.loads(raw_data)
                            # self.logger.info(data)
                            received_FFTxMax = np.array(data["FFTxMax"])
                            received_FFTmaxVal = np.array(data["FFTmaxVal"])

                            self.FFTxMax = np.append(self.FFTxMax, received_FFTxMax)
                            self.FFTmaxVal = np.append(self.FFTmaxVal, received_FFTmaxVal)
                            
                            FFTx = np.array(data["FFTx"])
                            FFTy = np.array(data["FFTy"])
                            OSCxin = np.array(data["OSCxin"])
                            OSCyin = np.array(data["OSCyin"])
                            OSCxout = np.array(data["OSCxout"])
                            OSCyout = np.array(data["OSCyout"])

                            # Check if enough time has elapsed since the last update
                            if current_time - last_update_time > 0.4:
                                last_update_time = current_time  # Update the timestamp
                    
                                # Clear previous line plot to avoid duplicate labels for FFT
                                for line in FFT.get_lines():
                                    line.remove()
                                
                                # Update scatter plots
                                if self.FFTxMax.any() and self.FFTmaxVal.any():
                                    FFT_max_scatter.set_offsets(np.c_[self.FFTxMax, self.FFTmaxVal])
                                    self.scatter_updated = True

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
                                if FFTy.any() and self.FFTmaxVal.any():
                                    FFTy_all = np.concatenate((FFTy, self.FFTmaxVal))
                                    if FFTy_all.size:
                                        FFT.set_ylim([min(FFTy_all) - 5, max(FFTy_all) + 5])

                                # Update Y-axis limits dynamically for Oscilloscope
                                if OSCyin.any() or OSCyout.any():
                                    oscilloscope_y_all = np.concatenate((OSCyin, OSCyout))
                                    if oscilloscope_y_all.size:
                                        padding = 0.5
                                        min_val = min(oscilloscope_y_all)
                                        max_val = max(oscilloscope_y_all)
                                        ValswingY = max(abs(int(min_val - padding)), abs(int(max_val + padding))) + 1
                                        oscilloscope.set_ylim([-ValswingY, ValswingY])

                                # Update X-axis limits dynamically for Oscilloscope
                                if OSCxin.any() or OSCxout.any():
                                    oscilloscope_x_all = np.concatenate((OSCxin, OSCxout))
                                    if oscilloscope_x_all.size:
                                        oscilloscope.set_xlim([min(oscilloscope_x_all), max(oscilloscope_x_all)])

                        except json.JSONDecodeError:
                            self.logger.error("Error decoding JSON data.")
                except socket.timeout:
                    pass

                # Handle plotting updates
                if self.scatter_updated:
                    # Clear previous lines and add spline fits
                    for line in FFT.get_lines():
                        line.remove()

                    for line in oscilloscope.get_lines():
                        line.remove()

                    # Add spline fits for FFT and Oscilloscope data
                    try:
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
                                self.logger.info(f"Error in fitting spline for FFT data: {e}")

                        if self.FFTxMax.any() and self.FFTmaxVal.any() and len(self.FFTxMax) >= 6:
                            try:
                                # Ensure FFTxMax is sorted and strictly increasing
                                sorted_indices = np.argsort(self.FFTxMax)
                                self.FFTxMax = np.array(self.FFTxMax)[sorted_indices]
                                self.FFTmaxVal = np.array(self.FFTmaxVal)[sorted_indices]
                                unique_indices = np.diff(self.FFTxMax) > 0
                                self.FFTxMax = self.FFTxMax[np.concatenate(([True], unique_indices))]
                                self.FFTmaxVal = self.FFTmaxVal[np.concatenate(([True], unique_indices))]

                                FFT_max_spline = UnivariateSpline(self.FFTxMax, self.FFTmaxVal, k=3, s=0)
                                FFT.plot(self.FFTxMax, FFT_max_spline(self.FFTxMax), color='red', label='FFT Max Spline')
                            except Exception as e:
                                self.logger.info(f"Error in fitting spline for FFT max data: {e}")

                        # Clear previous line plot for the oscilloscope plot
                        for line in oscilloscope.get_lines():
                            line.remove()
                        
                        if OSCxin.any() and OSCyin.any() and len(OSCxin) >= 6:
                            oscilloscope_in_spline = UnivariateSpline(OSCxin, OSCyin, k=3, s=0)
                            oscilloscope.plot(OSCxin, oscilloscope_in_spline(OSCxin), color='yellow', label='Input Spline')

                        if OSCxout.any() and OSCyout.any() and len(OSCxout) >= 6:
                            oscilloscope_out_spline = UnivariateSpline(OSCxout, OSCyout, k=3, s=0)
                            oscilloscope.plot(OSCxout, oscilloscope_out_spline(OSCxout), color='blue', label='Output Spline')

                    except Exception as e:
                        self.logger.error(f"Error fitting splines: {e}")
                
                # Increment the counter each loop iteration
                UpdateSkipCanvasCounter += 1
                if self.scatter_updated or UpdateSkipCanvasCounter >= 10:
                    # Allow GUI events to process
                    figure.canvas.flush_events()
                    time.sleep(0.01)
                    self.scatter_updated = False
                    # Reset the counter after the canvas update
                    UpdateSkipCanvasCounter = 0

        except Exception as e:
            self.logger.error(f"Exception during plotting loop: {e}")
        finally:
            # Ensure buffered data is sent before closing the socket
            if client_socket:  # Check if socket is initialized
                try:
                    if client_socket.fileno() != -1:  # Ensure socket is valid
                        self.logger.info("Shutting down socket for buffered data transmission...")
                        client_socket.shutdown(socket.SHUT_WR)  # Graceful shutdown for write operations
                        time.sleep(0.1)  # Allow any pending data to flush
                    else:
                        self.logger.warning("Socket is already closed or invalid. Skipping shutdown.")
                except socket.error as e:
                    self.logger.error(f"Error during socket shutdown: {e}")
                finally:
                    try:
                        if client_socket.fileno() != -1:  # Ensure socket is still valid before closing
                            client_socket.close()  # Fully close the socket
                            self.logger.info("Socket closed.")
                        else:
                            self.logger.warning("Socket already closed. Skipping close.")
                    except Exception as e:
                        self.logger.error(f"Error during socket close: {e}")
                        
        # Allow GUI events to process
        figure.canvas.flush_events()
        time.sleep(0.01)
        
    def read_output(self, process):
        """Reads and prints output line by line without decoding."""
        for line in iter(process.stdout.readline, ""):  # Use "" instead of b""
            print(line.strip())  # Simply strip without decode
            
    def start_plot_process(self, start_decade, stop_decade, points_per_decade):
        self.logger.info("Launching FFT plot process...")
        os.environ["PLOT_TYPE"] = "FFT"  # Set environment variable for FFT plotting
        current_os = platform.system()

        if current_os == "Linux":
            # Linux: Use multiprocessing for forking with FFT arguments
            self.logger.info("Using multiprocessing for FFT plotting on Linux.")
            self.plot_process = multiprocessing.Process(
                target=self.create_plot,
                args=(self.stop_event, start_decade, stop_decade, points_per_decade)
            )
            self.plot_process.start()

            if self.plot_process.is_alive():
                self.logger.info("FFT plot process started successfully!")
            else:
                self.logger.warning("Failed to start FFT plot process!")

        elif current_os == "Windows":
            # Windows: Use subprocess to launch the bundled executable.
            self.logger.info("Using subprocess for FFT plotting on Windows.")
            
            
            # Override env with environment variables
            env = os.environ.copy()
            env["PLOT_TYPE"] = "FFT"  # Pass environment variable
            
            # Use sys.executable to point to the bundled exe
            exe_path = sys.executable  # This should be the path to your Briefcase generated exe
            exe_dir = os.path.dirname(exe_path)
            self.logger.info(f"Launching bundled exe: {exe_path} with working directory: {exe_dir}")
            
            # Ensure null_handle is defined (for example, to suppress child process output)
            null_handle = open(os.devnull, 'w')
            atexit.register(null_handle.close)
            # Define startupinfo to adjust process startup options
            SW_SHOWNORMAL = 1
            SW_SHOWMINIMIZED = 2
            SW_SHOWMAXIMIZED = 3
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW | subprocess.STARTF_USESTDHANDLES
            startupinfo.wShowWindow = SW_SHOWNORMAL
            startupinfo.hStdInput = os.dup(0)             # inherited or explicitly fetched
            startupinfo.hStdOutput = null_handle.fileno() # Ensure null_handle is defined as: null_handle = open(os.devnull, 'w')
            startupinfo.hStdError = null_handle.fileno()
            
            # Consider using Popen instead of run if you want non-blocking behavior.
            processFFT = subprocess.Popen(
                [exe_path, "--plot-type", "FFT",
                "--start-decade", str(start_decade), 
                "--stop-decade", str(stop_decade), 
                "--points-per-decade", str(points_per_decade),
                "--base-port", str(self.base_port)],
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                env=env,
                cwd=exe_dir,  # Set working directory
                text=True,  # Ensure output is readable (string, not bytes)
                # startupinfo=startupinfo,
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS,  # Proper Windows process handling
                start_new_session=True  # Linux-friendly process detachment
            )
                
            ## Capture and print any outputs/errors for debugging
            ## Call the `read_output` method correctly
            thread = threading.Thread(target=self.read_output, args=(processFFT,))
            thread.start()
                
            def cleanupFFT():
                if processFFT.poll() is None:  # If process is still running
                    processFFT.terminate()  # Gracefully stop
                    processFFT.wait()  # Ensure full cleanup

            atexit.register(cleanupFFT)
            
        else:
            self.logger.error(f"Unsupported OS: {current_os}")
            raise NotImplementedError(f"OS '{current_os}' is not supported.")
        
    def set_data_point(self, conn, new_FFTxMax, new_FFTmaxVal):
        """Send new FFT max data points via socket."""
        try:
            payload = json.dumps({
                "FFTxMax": new_FFTxMax,
                "FFTmaxVal": new_FFTmaxVal
            })
            conn.sendall(payload.encode('utf-8'))
            conn.sendall(b"<END>")
        except Exception as e:
            self.logger.error(f"Error sending data point to FFT Oscilloscope: {e}")

    def set_full_data(self, conn, FFTx, FFTy, OSCxin, OSCyin, OSCxout, OSCyout):
        """Send the full dataset via socket."""
        try:
            payload = json.dumps({
                "FFTxMax": self.FFTxMax,  # Keep accumulating FFT max data
                "FFTmaxVal": self.FFTmaxVal,
                "FFTx": FFTx,
                "FFTy": FFTy,
                "OSCxin": OSCxin,
                "OSCyin": OSCyin,
                "OSCxout": OSCxout,
                "OSCyout": OSCyout
            })
            conn.sendall(payload.encode('utf-8'))
            conn.sendall(b"<END>")
        except Exception as e:
            self.logger.error(f"Error sending full dataset to FFT Oscilloscope: {e}")

    def is_running(self):
        return self.plot_process is not None and self.plot_process.is_alive()

    def stop_plot_process(self, timeout=5):
        """Stop the plotting process."""
        if self.plot_process and self.plot_process.is_alive():
            self.logger.info("Stopping XY plot process...")
            self.stop_event.set()
            self.plot_process.join(timeout=timeout)

            if self.plot_process.is_alive():
                self.logger.warning("Warning: XY plot process did not terminate within the timeout period.")
            else:
                self.logger.info("XY plot process stopped successfully!")
        else:
            self.logger.info("XY plot process is not running or already stopped.")

    def send_close_signal(self):
        # If there's no process, report and return.
        if self.plot_process is None:
            self.logger.info("Plot process is not running (None).")
            return

        if self.is_running():
            # Signal the child process to stop.
            self.stop_event.set()
            try:
                self.plot_process.join(timeout=0.5)
            except Exception as e:
                self.logger.error("Error joining plot process: %s", e)
                # As a fallback, try terminating the process.
                try:
                    self.plot_process.terminate()
                    self.plot_process.join(timeout=5)
                except Exception as e:
                    self.logger.error("Error terminating plot process: %s", e)
            self.logger.info("Plot process stopped.")
        else:
            self.logger.info("Plot process is not running.")
            
    @classmethod
    def get_running_instances(cls):
        """Returns a list of running PlotManagerFFToscilloscope instances."""
        with cls._instances_lock:
            return [instance for instance in cls._instances if instance.is_running()]
# --- End class PlotManager FFT oscilloscope:--- --------------------------------------

