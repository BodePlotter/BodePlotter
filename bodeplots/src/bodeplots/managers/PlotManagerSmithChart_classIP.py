"""

© 2025 by Bode Plotter

MIT License

Copyright (c) [2025] [Bode Plotter]
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

import socket
import json
import logging
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
import skrf as rf # pip install scikit-rf

# --- class PlotManagerSmithChart: --- -------------------------------------------
class PlotManagerSmithChart:
    """
    This class sets up a Matplotlib figure displaying a Smith chart, connects as a 
    client to a socket server on port 5004 (or the port defined by BASE_PORT), and 
    updates the chart with received JSON-formatted Bode plot data. The JSON data should 
    include the keys:
      - "frequency": list of frequency values in Hz,
      - "magnitude_db": list of magnitudes (in dB),
      - "phase_degrees": list of phases (in degrees).
    
    It includes robust error handling, exponential backoff for connection retries, window 
    management via PyQt, and a graceful shutdown of the socket connection.
    """

    # Class-level registry to hold instances with weak references
    _instances = weakref.WeakSet()
    _instances_lock = threading.Lock()  # For thread-safety on the registry

    def __init__(self):
        """Initialize the PlotManagerSmithChart."""
        # Set up the logger using the global logging configuration.
        self.logger = logging.getLogger(__name__)
        self.logger.info("PlotManagerSmithChart initialized.")

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
        self.port = self.base_port + 3
        self.stop_reading = False  # Flag to control whether to stop reading from the socket
        
        # These will be set later in create_plot()
        self.ax = None

        # Register this instance in a thread-safe manner.
        with PlotManagerSmithChart._instances_lock:
            PlotManagerSmithChart._instances.add(self)

    def initialize_smith_chart(self, draw_labels=True, chart_type='z'):
        """Clear the axis and draw the Smith chart background using the specified options.
        
        Additionally, the subplot layout is adjusted so the top margin is set at 85% of the figure height.
        This prevents overlap between the title (manually positioned with y=0.95) and other chart elements.
        """
        if self.ax is not None:
            self.ax.clear()
            # Draw the Smith chart background using scikit-rf's smith plotting function
            rf.plotting.smith(ax=self.ax, draw_labels=draw_labels, chart_type=chart_type)
            self.ax.set_title('Smith Chart from Bode Plot Data with Frequency Annotations')
            # Adjust the overall subplot layout: set the top margin at 80% of the figure's height.
            fig = self.ax.get_figure()
            fig.subplots_adjust(top=0.80)

    def segment_smith_data(self, real_data, imag_data, radius=1.0):
        """
        Process data points for a Smith chart line.
        
        Parameters:
            real_data: list or array of real parts.
            imag_data: list or array of imaginary parts.
            radius: the defining radius of the Smith chart (default is 1).
            
        Returns:
            new_real: list of real parts with NaNs inserted to break segments.
            new_imag: list of imaginary parts with NaNs inserted accordingly.
        """
        if not real_data or not imag_data:
            return [], []
            
        new_real = [real_data[0]]
        new_imag = [imag_data[0]]

        def is_outside(x, y, r=radius):
            return x**2 + y**2 > r**2

        for i in range(1, len(real_data)):
            prev_out = is_outside(real_data[i-1], imag_data[i-1])
            curr_out = is_outside(real_data[i], imag_data[i])
            if prev_out or curr_out:
                new_real.append(real_data[i])
                new_imag.append(imag_data[i])
            else:
                new_real.append(np.nan)
                new_imag.append(np.nan)
                new_real.append(real_data[i])
                new_imag.append(imag_data[i])
        return new_real, new_imag
    
    def create_plot(self, stop_event, cache_dir=None):
        """
        Set up the interactive Matplotlib figure, position the window,
        connect to the socket server with retry logic, and process incoming JSON messages
        to update the Smith chart.
        """
        # Platform-specific cache directory settings
        current_os = platform.system()
        if current_os == "Windows":
            if cache_dir is None:
                timestamp = int(time.time() * 1000)
                cache_dir = f"C:\\Temp\\matplotlib_cache_SC_{timestamp}"
                os.makedirs(cache_dir, exist_ok=True)
            os.environ["MPLCONFIGDIR"] = cache_dir
            time.sleep(0.5)

        # Configure Matplotlib and create the interactive figure.
        import matplotlib
        matplotlib.use("qtagg")
        import matplotlib.pyplot as plt
        plt.ion()
        figure, self.ax = plt.subplots(figsize=(6, 6))
        # Initialize the Smith chart only once.
        self.initialize_smith_chart()
        figure.canvas.manager.set_window_title("Smith Chart from Bode Plot Data with Frequency Annotations")
        figure.tight_layout()

        # --- Positioning with PyQt (same as before) ---
        manager = plt.get_current_fig_manager()
        if hasattr(manager, 'window'):
            try:
                window = manager.window
                desired_x, desired_y = 2060, 1400
                geom = window.frameGeometry()
                window_width = geom.width()
                window_height = geom.height()

                from PyQt5.QtWidgets import QApplication
                app = QApplication.instance()
                if app is None:
                    app = QApplication(sys.argv)
                screen = app.primaryScreen()
                screen_geom = screen.availableGeometry()
                screen_width = screen_geom.width()
                screen_height = screen_geom.height()

                x = max(0, min(desired_x, screen_width - window_width))
                y = max(0, min(desired_y, screen_height - window_height))
                window.move(x, y)
            except Exception as e:
                self.logger.error(f"Error setting window geometry: {e}")

        def on_close(event):
            self.logger.info("Plot window closed by user.")
            stop_event.set()
        figure.canvas.mpl_connect('close_event', on_close)
        self.logger.info("Figure initialized and close event connected.")

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

        # --- Initialize Accumulators and Plot Line ---
        accumulated_real = []
        accumulated_imag = []
        accumulated_frequencies = []  # optional, if you want to accumulate frequency labels

        # Create an initial empty plot line for the accumulated data.
        line, = self.ax.plot([], [], 'b-o', label='S11 Data')
        self.ax.legend()

        buffer = ""
        last_update_time = 0  # Initialize the timestamp
        UpdateSkipCanvasCounter = 0  # Initialize the counter to 0
        update_counter = 0  # Tracks updates over time
        annotation_interval = 50  # Annotate every 50 updates
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
                            data = json.loads(raw_data)
                            # self.logger.info("Data: %s", data)
                            # Validate essential keys.
                            if (data.get("magnitude_db") is None or 
                                data.get("frequency") is None or 
                                data.get("phase_degrees") is None):
                                self.logger.error("One or more required keys (magnitude_db, frequency, phase_degrees) are missing. Skipping this message.")
                                continue

                            frequency = np.atleast_1d(data.get("frequency"))
                            magnitude_db = np.atleast_1d(data.get("magnitude_db"))
                            phase_degrees = np.atleast_1d(data.get("phase_degrees"))
                            
                            """
                            # Log the input values to help with debugging/verification
                            self.logger.info("Input Frequency Values: %s", frequency)
                            self.logger.info("Input Magnitude (dB): %s", magnitude_db)
                            self.logger.info("Input Phase (degrees): %s", phase_degrees)
                            """
                            
                            magnitude_linear = 10 ** (magnitude_db / 20)
                            phase_radians = np.radians(phase_degrees)
                            s11 = magnitude_linear * np.exp(1j * phase_radians)
                            
                            """
                            # Log computed intermediate values
                            self.logger.info("Computed Linear Magnitude: %s", magnitude_linear)
                            self.logger.info("Computed Phase (radians): %s", phase_radians)
                            self.logger.info("Computed s11 (complex S-parameter): %s", s11)
                            """
                            
                            # --- Accumulate new data ---
                            accumulated_real.extend(s11.real.tolist())
                            accumulated_imag.extend(s11.imag.tolist())
                            accumulated_frequencies.extend(frequency.tolist())  # if needed for annotations

                            # Check if enough time has elapsed since the last update
                            if current_time - last_update_time > 0.4:
                                last_update_time = current_time  # Update the timestamp

                                # Process the accumulated data to insert breaks where
                                # both consecutive points are inside the Smith chart.
                                plot_x, plot_y = self.segment_smith_data(accumulated_real, accumulated_imag)

                                # Update the existing line with the processed data.
                                line.set_data(plot_x, plot_y)

                                # Compute the current data limits from the accumulated data.
                                data_min_real = np.min(accumulated_real)
                                data_max_real = np.max(accumulated_real)
                                data_min_imag = np.min(accumulated_imag)
                                data_max_imag = np.max(accumulated_imag)

                                # --- Define the fixed Smith Chart boundaries ---
                                SMITH_REAL_BOUND = (-1, 1)
                                SMITH_IMAG_BOUND = (-1, 1)

                                # Determine if the new points exceed the fixed boundaries.
                                update_xlim = data_min_real < SMITH_REAL_BOUND[0] or data_max_real > SMITH_REAL_BOUND[1]
                                update_ylim = data_min_imag < SMITH_IMAG_BOUND[0] or data_max_imag > SMITH_IMAG_BOUND[1]

                                if update_xlim or update_ylim:
                                    new_xlim = (
                                        data_min_real if data_min_real < SMITH_REAL_BOUND[0] else SMITH_REAL_BOUND[0],
                                        data_max_real if data_max_real > SMITH_REAL_BOUND[1] else SMITH_REAL_BOUND[1]
                                    )
                                    new_ylim = (
                                        data_min_imag if data_min_imag < SMITH_IMAG_BOUND[0] else SMITH_IMAG_BOUND[0],
                                        data_max_imag if data_max_imag > SMITH_IMAG_BOUND[1] else SMITH_IMAG_BOUND[1]
                                    )
                                else:
                                    new_xlim = SMITH_REAL_BOUND
                                    new_ylim = SMITH_IMAG_BOUND

                                # --- Add margin to leave extra space for the title and chart labels ---
                                margin_fraction = 0.1  # 10% margin; adjust as needed

                                # Compute the padding for x and y axes.
                                x_range = new_xlim[1] - new_xlim[0]
                                y_range = new_ylim[1] - new_ylim[0]
                                pad_x = x_range * margin_fraction
                                pad_y = y_range * margin_fraction

                                # Update the axes limits with the added padding.
                                self.ax.set_xlim(new_xlim[0] - pad_x, new_xlim[1] + pad_x)
                                self.ax.set_ylim(new_ylim[0] - pad_y, new_ylim[1] + pad_y)

                                
                                                                
                                """
                                # Log the accumulated results
                                self.logger.info("Accumulated Real Parts: %s", accumulated_real)
                                self.logger.info("Accumulated Imaginary Parts: %s", accumulated_imag)
                                self.logger.info("Accumulated Frequencies: %s", accumulated_frequencies)
                                """
                                
                                update_counter += 1  # Increment counter with each new data update
                                
                                if update_counter % annotation_interval == 0:  # Annotate every 50 updates
                                    step = max(30, len(accumulated_frequencies) // 60)  # Dynamically adjust step size
                                    for i, (real, imag, f) in enumerate(zip(accumulated_real, accumulated_imag, accumulated_frequencies)):
                                        if i % step == 0:  # Annotate selected points
                                            if f < 1e3:
                                                label = f"{f:.0f} Hz"
                                            elif f < 1e6:
                                                label = f"{f/1e3:.0f} kHz"
                                            else:
                                                label = f"{f/1e6:.0f} MHz"

                                            self.ax.annotate(label, xy=(real, imag),
                                                            textcoords='offset points', xytext=(15, 15),
                                                            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

                                # Reset counter occasionally to prevent overflow
                                if update_counter > 1e6:
                                    update_counter = 0
                                    
                        except json.JSONDecodeError:
                            self.logger.error("Error decoding JSON data. Skipping this message.")
                            # Allow GUI events to process
                        
                except socket.timeout:
                    pass  # No data received, continue loop
                
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
            if client_socket:
                try:
                    if client_socket.fileno() != -1:
                        self.logger.info("Shutting down socket for buffered data transmission...")
                        client_socket.shutdown(socket.SHUT_WR)
                        time.sleep(0.1)
                    else:
                        self.logger.warning("Socket is already closed or invalid. Skipping shutdown.")
                except socket.error as e:
                    self.logger.error(f"Error during socket shutdown: {e}")
                finally:
                    try:
                        if client_socket.fileno() != -1:
                            client_socket.close()
                            self.logger.info("Socket closed.")
                        else:
                            self.logger.warning("Socket already closed. Skipping close.")
                    except Exception as e:
                        self.logger.error(f"Error during socket close: {e}")

            figure.canvas.flush_events()
            time.sleep(0.01)

        
    def read_output(self, process):
        """Reads and prints output line by line (process.stdout is in text mode)."""
        for line in iter(process.stdout.readline, ""):
            print(line.strip())
    
    def start_plot_process(self):
        """
        Launch the Smith Chart plot process.
        For Linux, we use multiprocessing to call create_plot().
        For Windows, we assume a bundled executable is launched via subprocess.
        
        Note: The extra parameters (start_decade, stop_decade, points_per_decade)
        are used only when launching via subprocess; the local create_plot() only 
        requires the stop_event (and optionally a cache_dir).
        """
        self.logger.info("Launching Smith Chart plot process...")
        os.environ["PLOT_TYPE"] = "SC"  # Set environment variable for Smith Chart plotting
        current_os = platform.system()
        
        if current_os == "Linux":
            self.logger.info("Using multiprocessing for Smith Chart plotting on Linux.")
            # Pass only the stop event since create_plot() accepts (stop_event, cache_dir=None)
            self.plot_process = multiprocessing.Process(
                target=self.create_plot,
                args=(self.stop_event,)
            )
            self.plot_process.start()
            
            if self.plot_process.is_alive():
                self.logger.info("Smith Chart plot process started successfully!")
            else:
                self.logger.warning("Failed to start Smith Chart plot process!")
    
        elif current_os == "Windows":
            self.logger.info("Using subprocess for Smith Chart plotting on Windows.")
            
            # Override environment variables
            env = os.environ.copy()
            env["PLOT_TYPE"] = "SC"
            
            # Use sys.executable; in a bundled application this should be the packaged exe.
            exe_path = sys.executable
            exe_dir = os.path.dirname(exe_path)
            self.logger.info(f"Launching bundled exe: {exe_path} with working directory: {exe_dir}")
            
            # Define a null handle to suppress child process output.
            null_handle = open(os.devnull, 'w')
            atexit.register(null_handle.close)
            
            # Define startupinfo if needed.
            SW_SHOWNORMAL = 1
            SW_SHOWMINIMIZED = 2
            SW_SHOWMAXIMIZED = 3
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW | subprocess.STARTF_USESTDHANDLES
            startupinfo.wShowWindow = SW_SHOWNORMAL
            startupinfo.hStdInput = os.dup(0)
            startupinfo.hStdOutput = null_handle.fileno()
            startupinfo.hStdError = null_handle.fileno()
            
            # Build command-line arguments. (The extra parameters are used if needed by the exe.)
            args = [exe_path, "--plot-type", "SC", "--base-port", str(self.base_port)]
                
            processSC = subprocess.Popen(
                args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
                cwd=exe_dir,  # Set working directory
                text=True,    # Output in text mode
                # startupinfo=startupinfo,  # Uncomment if needed
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS,
                start_new_session=True
            )
                
            # Start a thread to capture and print output for debugging.
            thread = threading.Thread(target=self.read_output, args=(processSC,))
            thread.start()
                
            def cleanupSC():
                if processSC.poll() is None:  # If process is still running
                    processSC.terminate()  # Gracefully stop
                    processSC.wait()       # Ensure proper cleanup
    
            atexit.register(cleanupSC)
        else:
            self.logger.error(f"Unsupported OS: {current_os}")
            raise NotImplementedError(f"OS '{current_os}' is not supported.")
    
    
    def is_running(self):
        """Check if the plot process is running."""
        return self.plot_process is not None and self.plot_process.is_alive()
    
    
    def stop_plot_process(self):
        """Stop the Smith Chart plotting process."""
        if self.plot_process and self.plot_process.is_alive():
            self.logger.info("Stopping Smith Chart plot process...")
            self.stop_event.set()
            self.plot_process.join()
            if not self.plot_process.is_alive():
                self.logger.info("Smith Chart plot process stopped successfully.")
            else:
                self.logger.warning("Failed to terminate Smith Chart plot process.")
    
    
    def send_close_signal(self):
        """
        Signal the child plotting process to stop and wait for its termination.
        If the child process does not turn off in a timely manner, attempt to terminate it.
        """
        if self.plot_process is None:
            self.logger.info("Plot process is not running (None).")
            return
    
        if self.is_running():
            self.stop_event.set()
            try:
                self.plot_process.join(timeout=0.5)
            except Exception as e:
                self.logger.error("Error joining plot process: %s", e)
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
        """Returns a list of running PlotManagerSmithChart instances."""
        with cls._instances_lock:
            return [instance for instance in cls._instances if instance.is_running()]
# --- End class PlotManagerSmithChart:--- ---------------------------------------

