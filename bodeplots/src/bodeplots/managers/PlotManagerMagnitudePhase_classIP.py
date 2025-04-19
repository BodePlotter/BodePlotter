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
from multiprocessing import Process, Event
import threading
import weakref
from scipy.interpolate import UnivariateSpline
import os
import time
import random
import platform
import subprocess
import sys
import atexit

# --- class PlotManager Magnitude Phase:--- -------------------------------------------
class PlotManagerMagnitudePhase:
    # Class-level registry to hold instances with weak references
    _instances = weakref.WeakSet()
    _instances_lock = threading.Lock()  # For thread-safety on the registry

    def __init__(self):
        """Initialize the PlotManagerMagnitudePhase."""
        # Set up the logger using the global logging configuration.
        self.logger = logging.getLogger(__name__)
        self.logger.info("PlotManagerMagnitudePhase initialized.")

        # Initialize class attributes
        self.stop_event = Event()  # Event to stop the plotting process
        self.scatter_updated = False  # Flag to indicate scatter plot update
        self.plot_process = None  # Process for handling Matplotlib
        self.host = '127.0.0.1'  # Socket host
        try:
            self.base_port = int(os.getenv("BASE_PORT", 5001))
        except ValueError:
            # Optionally log the error or set a fallback base port
            self.base_port = 5001
        self.port = self.base_port
        self.stop_reading = False  # Flag to control whether to stop reading from the socket

        # Initialize variables as numpy arrays
        self.gain_X = np.array([])  # Store all gain X data
        self.gain_Y = np.array([])  # Store all gain Y data
        self.phase_X = np.array([])  # Store all phase X data
        self.phase_Y = np.array([])  # Store all phase Y data

        # Register the new instance in a thread-safe manner
        with PlotManagerMagnitudePhase._instances_lock:
            PlotManagerMagnitudePhase._instances.add(self)
            
    def sort_gain_data(self):
        """Sort gain_X and gain_Y together while maintaining their relationship."""
        if self.gain_X.size > 0 and self.gain_Y.size > 0:  # Ensure non-empty arrays
            sorted_indices = np.argsort(self.gain_X)  # Sort based on gain_X
            self.gain_X = self.gain_X[sorted_indices]
            self.gain_Y = self.gain_Y[sorted_indices]

    def sort_phase_data(self):
        """Sort phase_X and phase_Y together while maintaining their relationship."""
        if self.phase_X.size > 0 and self.phase_Y.size > 0:  # Ensure non-empty arrays
            sorted_indices = np.argsort(self.phase_X)  # Sort based on phase_X
            self.phase_X = self.phase_X[sorted_indices]
            self.phase_Y = self.phase_Y[sorted_indices]

    def perform_outlier_detection(self, phase_X, phase_Y, window_size=12, threshold=10):
        """Detect and filter outliers based on local deviation."""
        filtered_indices = []
        for i in range(phase_Y.size):
            start = max(0, i - window_size)
            end = min(phase_Y.size, i + window_size + 1)
            neighbors = np.concatenate([phase_Y[start:i], phase_Y[i + 1:end]])

            if neighbors.size > 0:
                local_median = np.median(neighbors)
                if abs(phase_Y[i] - local_median) <= threshold:
                    filtered_indices.append(i)

        return phase_X[filtered_indices], phase_Y[filtered_indices]

    def create_plot(self, stop_event, start_decade, stop_decade, points_per_decade, cache_dir=None):
        current_os = platform.system()
        if current_os == "Windows":
            if cache_dir is None:
                timestamp = int(time.time() * 1000)  # Get current time in milliseconds
                cache_dir = f"C:\\Temp\\matplotlib_cache_MP_{timestamp}"
                os.makedirs(cache_dir, exist_ok=True)  # Create the directory if it doesn't exist

            # Update environment variable properly
            os.environ["MPLCONFIGDIR"] = cache_dir
            time.sleep(0.5)  # Short wait to ensure directory is ready 

        # Import Matplotlib
        import matplotlib
        matplotlib.use("qtagg")  # qtagg preferred interactive backend
        import matplotlib.pyplot as plt

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

        # Configure Matplotlib settings
        plt.ion()
         
        # Adjust figsize based on your desired overall window size.
        # For example, if you want the entire window to be ~600 pixels high at 100 DPI:
        figure, (magnitude, phase) = plt.subplots(2, figsize=(12, 6))

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
        phase.set_ylim([-200, 200])  # Fixed y-axis range for phase
        phase.grid(which='both')

        # Initialize scatter plots
        magnitude_scatter = magnitude.scatter([], [], color='red', label='Raw Data', s=20)
        phase_scatter = phase.scatter([], [], color='orange', label='Raw Data', s=20)

        figure.tight_layout()

        # (After creating your figure and before entering the main event loop.)
        manager = plt.get_current_fig_manager()
        if hasattr(manager, 'window'):
            try:
                window = manager.window

                # Desired initial position (these values can be any starting point)
                desired_x, desired_y = 1450, 10

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

        # Handle window close event
        def on_close(event):
            self.logger.info("Plot window closed by user.")
            stop_event.set()

        figure.canvas.mpl_connect('close_event', on_close)
        self.logger.info("Figure initialized and close event connected.")

        # Data buffer
        buffer = ""
        
        last_update_time = 0  # Initialize the timestamp
        UpdateSkipCanvasCounter = 0  # Initialize the counter to 0
        
        try:
            while not stop_event.is_set():
                current_time = time.time()
                try:
                    # Receive data via socket
                    if not self.stop_reading and client_socket and client_socket.fileno() != -1:
                        chunk = client_socket.recv(2560).decode('utf-8')  # Large buffer size for efficiency
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
                        
                    # Process complete messages
                    if "<END>" in buffer:
                        raw_data, buffer = buffer.split("<END>", 1)

                        try:
                            # Parse the received JSON data
                            data = json.loads(raw_data)

                            # Insert new gain data and sort
                            self.gain_X = np.append(self.gain_X, data["gain_X"])
                            self.gain_Y = np.append(self.gain_Y, data["gain_Y"])
                            self.sort_gain_data()

                            # Insert new phase data and sort
                            self.phase_X = np.append(self.phase_X, data["phase_X"])
                            self.phase_Y = np.append(self.phase_Y, data["phase_Y"])
                            self.sort_phase_data()
                            
                        except json.JSONDecodeError:
                            self.logger.error("Error decoding JSON data. Skipping this message.")

                except socket.timeout:
                    pass  # No data received, continue loop

                # Check if enough time has elapsed since the last update
                if current_time - last_update_time > 0.4:
                    last_update_time = current_time  # Update the timestamp
                    # Update scatter plots
                    try:
                        if self.gain_X.size == self.gain_Y.size:
                            magnitude_scatter.set_offsets(np.c_[self.gain_X, self.gain_Y])
                        if self.phase_X.size == self.phase_Y.size:
                            phase_scatter.set_offsets(np.c_[self.phase_X, self.phase_Y])
                        self.scatter_updated = True
                    except Exception as e:
                        self.logger.error(f"Error updating scatter plots: {e}")

                    # Update Y-axis limits dynamically
                    if self.gain_Y.size > 0:
                        magnitude.set_ylim([np.min(self.gain_Y) - 5, np.max(self.gain_Y) + 5])
                    if self.phase_Y.size > 0:
                        phase.set_ylim([np.min(self.phase_Y) - 5, np.max(self.phase_Y) + 5])
                    if self.scatter_updated:
                        # Clear previous lines and add spline fits
                        for line in magnitude.get_lines():
                            line.remove()
                        for line in phase.get_lines():
                            line.remove()

                        try:
                            # Spline fitting for Magnitude
                            if self.gain_X.size >= 4:
                                magnitude_spline = UnivariateSpline(self.gain_X, self.gain_Y, k=3, s=0)
                                magnitude.plot(self.gain_X, magnitude_spline(self.gain_X), color='blue', label='Magnitude Spline')
                            else:
                                self.logger.warning("Insufficient points for Magnitude Spline. Skipping fit.")

                            # Spline fitting for Phase
                            if self.phase_X.size >= 4:
                                filtered_phase_X, filtered_phase_Y = self.perform_outlier_detection(
                                    self.phase_X, self.phase_Y
                                )
                                if filtered_phase_X.size >= 4:
                                    phase_spline = UnivariateSpline(filtered_phase_X, filtered_phase_Y, k=3, s=0)
                                    phase.plot(filtered_phase_X, phase_spline(filtered_phase_X), color='blue', label='Phase Spline')
                            else:
                                self.logger.warning("Insufficient points for Phase Spline. Skipping fit.")
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
            self.logger.error(f"Exception in create_plot: {e}")
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
        self.logger.info("Launching Magnitude/Phase (MP) plot process...")
        os.environ["PLOT_TYPE"] = "MP"  # Set environment variable for MP plotting
        current_os = platform.system()
        
        if current_os == "Linux":
            # Linux: Use multiprocessing for forking with additional MP arguments
            self.logger.info("Using multiprocessing for MP plotting on Linux.")
            self.plot_process = multiprocessing.Process(
                target=self.create_plot,
                args=(self.stop_event, start_decade, stop_decade, points_per_decade)
            )
            self.plot_process.start()
            
            if self.plot_process.is_alive():
                self.logger.info("MP plot process started successfully!")
            else:
                self.logger.warning("Failed to start MP plot process!")

        elif current_os == "Windows":
            # Windows: Use subprocess to launch the bundled executable.
            self.logger.info("Using subprocess for MP plotting on Windows.")
            
            
            # Override env with environment variables
            env = os.environ.copy()
            env["PLOT_TYPE"] = "MP"  # Pass environment variable
            
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
            processMP = subprocess.Popen(
                [exe_path, "--plot-type", "MP", 
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
            thread = threading.Thread(target=self.read_output, args=(processMP,))
            thread.start()
                
            def cleanupMP():
                if processMP.poll() is None:  # If process is still running
                    processMP.terminate()  # Gracefully stop
                    processMP.wait()  # Ensure full cleanup

            atexit.register(cleanupMP)
            
        else:
            self.logger.error(f"Unsupported OS: {current_os}")
            raise NotImplementedError(f"OS '{current_os}' is not supported.")


    def is_running(self):
        return self.plot_process is not None and self.plot_process.is_alive()

    def stop_plot_process(self):
        """Stop the plotting process."""
        if self.plot_process and self.plot_process.is_alive():
            self.logger.info("Stopping Magnitude/Phase plot process...")
            self.stop_event.set()
            self.plot_process.join()
            if not self.plot_process.is_alive():
                self.logger.info("Magnitude/Phase plot process stopped successfully.")
            else:
                self.logger.warning("Failed to terminate Magnitude/Phase plot process.")

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
        """Returns a list of running PlotManagerMagnitudePhase instances."""
        with cls._instances_lock:
            return [instance for instance in cls._instances if instance.is_running()]
# --- End class PlotManager Magnitude Phase:--- ---------------------------------------

