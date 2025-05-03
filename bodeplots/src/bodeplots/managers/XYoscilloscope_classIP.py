"""

© 2025 by Bode Plotter

MIT License

Copyright (c) [2025] [Bode Plotter]
Bode plot implementation for the handheld oscilloscope OWON HDS320S
Highly rewritten and modified between 01/20/2025 to 05/04/2025
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
# import matplotlib.pyplot as plt
import numpy as np
import multiprocessing
import threading
import weakref
import os
import logging
import time
import random
import platform
import subprocess
import sys
import atexit

# --- class PlotManager XYoscilloscope:--- --------------------------------------------
class XYoscilloscope:
    # Class-level registry to hold instances with weak references
    _instances = weakref.WeakSet()
    _instances_lock = threading.Lock()  # For thread-safety on the registry

    def __init__(self):
        """Initialize the XYoscilloscope class."""
        # Set up the logger using the global logging configuration.
        self.logger = logging.getLogger(__name__)
        self.logger.info("Initializing XYoscilloscope.")

        # Initialize class attributes
        self.scatter_updated = False  # Flag to indicate scatter plot update
        self.host = '127.0.0.1'       # Socket host
        try:
            self.base_port = int(os.getenv("BASE_PORT", 5001))
        except ValueError:
            # Optionally log the error or set a fallback base port
            self.base_port = 5001
        self.port = self.base_port + 2
        self.stop_reading = False   # Flag to control whether to stop reading
        self.stop_event = multiprocessing.Event()  # Event for stopping the process
        self.plot_process = None    # Plotting process

        # Register the new instance in a thread-safe manner
        with XYoscilloscope._instances_lock:
            XYoscilloscope._instances.add(self)

        # Log initialization details
        self.logger.info(f"XYoscilloscope initialized with host {self.host} and port {self.port}.")

            
    def create_plot(self, stop_event, cache_dir=None):
        current_os = platform.system()
        if current_os == "Windows":
            if cache_dir is None:
                timestamp = int(time.time() * 1000)  # Get current time in milliseconds
                cache_dir = f"C:\\Temp\\matplotlib_cache_XY_{timestamp}"
                os.makedirs(cache_dir, exist_ok=True)  # Create the directory if it doesn't exist
            # Update environment variable properly
            os.environ["MPLCONFIGDIR"] = cache_dir
            time.sleep(0.5)  # Short wait to ensure directory is ready 

        # Import Matplotlib
        import matplotlib
        matplotlib.use("qtagg")  # qtagg preferred interactive backend
        import matplotlib.pyplot as plt

        # Remove the inherited process object in the child
        """Matplotlib process to receive data via sockets and render plots."""

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
        plt.ion()  # Enable interactive mode
        
        # Adjust figsize based on your desired overall window size.
        # For example, if you want the entire window to be ~600 pixels high at 100 DPI:
        self.figure, XY_plot = plt.subplots(figsize=(6, 6))
        
        XY_scatter = XY_plot.scatter([0], [0], color='red', label='XY Data', s=20)
        XY_plot.set_title("Oscilloscope X-Y Plot")
        XY_plot.set_xlabel("raw_waveform_in")
        XY_plot.set_ylabel("raw_waveform_out")
        XY_plot.grid(True)

        XY_plot.legend(loc='lower left', bbox_to_anchor=(0, 0), borderaxespad=0, frameon=True)
        self.figure.canvas.manager.set_window_title("Oscilloscope X-Y Plot")

        # (After creating your figure and before entering the main event loop.)
        manager = plt.get_current_fig_manager()
        if hasattr(manager, 'window'):
            try:
                window = manager.window

                # Desired initial position (these values can be any starting point)
                desired_x, desired_y = 1450, 1400

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

        # Define close event handler
        def on_close(event):
            self.logger.info("Matplotlib window closed.")
            stop_event.set()

        self.figure.canvas.mpl_connect('close_event', on_close)
        self.logger.info("Figure initialized and close event connected.")

        buffer = ""  # Buffer for incoming data chunks
        UpdateSkipCanvasCounter = 0  # Initialize the counter to 0
        
        try:
            while not stop_event.is_set():
                try:
                    if not self.stop_reading and client_socket and client_socket.fileno() != -1:
                        chunk = client_socket.recv(20480).decode('utf-8')
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
                            data = json.loads(raw_data)  # Parse JSON data
                            raw_waveform_in = data["raw_waveform_in"]
                            raw_waveform_out = data["raw_waveform_out"]
                            
                            # Validate that the data lengths are consistent
                            if len(raw_waveform_in) != len(raw_waveform_out):
                                self.logger.error(
                                    f"Inconsistent data lengths: raw_waveform_in ({len(raw_waveform_in)}) "
                                    f"!= raw_waveform_out ({len(raw_waveform_out)}). Skipping update."
                                )
                                continue  # Skip this update cycle if lengths don't match

                            x_padding = (max(raw_waveform_in) - min(raw_waveform_in)) * 0.1
                            y_padding = (max(raw_waveform_out) - min(raw_waveform_out)) * 0.1

                            XY_scatter.set_offsets(np.c_[raw_waveform_in, raw_waveform_out])
                            XY_plot.set_xlim([min(raw_waveform_in) - x_padding, max(raw_waveform_in) + x_padding])
                            XY_plot.set_ylim([min(raw_waveform_out) - y_padding, max(raw_waveform_out) + y_padding])
                            plt.draw()
                            self.scatter_updated = True
                        except json.JSONDecodeError:
                            self.logger.error("Error decoding JSON. Skipping this chunk.")
                            
                except socket.timeout:
                    pass


                # Increment the counter each loop iteration
                UpdateSkipCanvasCounter += 1
                if self.scatter_updated or UpdateSkipCanvasCounter >= 10:
                    # Allow GUI events to process
                    self.figure.canvas.flush_events()
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
                        
    def read_output(self, process):
        """Reads and prints output line by line without decoding."""
        for line in iter(process.stdout.readline, ""):  # Use "" instead of b""
            print(line.strip())  # Simply strip without decode
            
    def start_plot_process(self):
        self.logger.info("Launching XY plot process...")
        os.environ["PLOT_TYPE"] = "XY"  # Set environment variable for XY plotting
        current_os = platform.system()
        
        if current_os == "Linux":
            # Linux: Use multiprocessing for forking
            self.logger.info("Using multiprocessing for XY plotting on Linux.")
            self.plot_process = multiprocessing.Process(
                target=self.create_plot,
                args=(self.stop_event,)
            )
            self.plot_process.start()
            
            if self.plot_process.is_alive():
                self.logger.info("XY plot process started successfully!")
            else:
                self.logger.warning("Failed to start XY plot process!")
        
        elif current_os == "Windows":
            # Windows: Use subprocess to launch the bundled executable.
            self.logger.info("Using subprocess for XY plotting on Windows.")
            
            
            # Override env with environment variables
            env = os.environ.copy()
            env["PLOT_TYPE"] = "XY"  # Pass environment variable
            
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
            processXY = subprocess.Popen(
                [exe_path, "--plot-type", "XY",
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
            thread = threading.Thread(target=self.read_output, args=(processXY,))
            thread.start()
    
            def cleanupXY():
                if processXY.poll() is None:  # If process is still running
                    processXY.terminate()  # Gracefully stop
                    processXY.wait()  # Ensure full cleanup

            atexit.register(cleanupXY)
            
        else:
            self.logger.error(f"Unsupported OS: {current_os}")
            raise NotImplementedError(f"OS '{current_os}' is not supported.")
    


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

    def is_running(self):
        return self.plot_process is not None and self.plot_process.is_alive()


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
        """Returns a list of running XYoscilloscope instances."""
        with cls._instances_lock:
            return [instance for instance in cls._instances if instance.is_running()]
# --- End class PlotManager XYoscilloscope:--- ----------------------------------------
