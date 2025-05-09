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
import os
import sys
import platform
import multiprocessing
import datetime
import logging
import argparse
from uuid import uuid4
import shutil  # Import for cleanup
import atexit  # Import for safe cleanup upon exit
# 2025/03/19 Import the matplotlib class for bode-plotter
from bodeplots.managers import PlotManagerFFToscilloscope, PlotManagerMagnitudePhase, XYoscilloscope, PlotManagerSmithChart
from bodeplots.app import main
import tempfile
import time

def robust_cleanup(directory):
    # Note: Avoid using print if there's any chance sys.stdout might be closed.
    for _ in range(10):
        try:
            shutil.rmtree(directory, ignore_errors=False)
            return
        except Exception as e:
            time.sleep(0.2)
    shutil.rmtree(directory, ignore_errors=True)

# Moved over from __main__.py
# --- Command-line Parsing for Plot Type ---
parser = argparse.ArgumentParser(description="BodePlots Application")
parser.add_argument(
    "--plot-type",
    choices=["XY", "MP", "FFT", "SC", "BODEPLOTTER"],
    default=None,
    help="Select the type of plot to display. Options: XY, MP, FFT, SC, or BODEPLOTTER."
)

parser.add_argument(
    "--start-decade",
    type=float,
    required=False,
    help="Specify the starting decade for data processing."
)

parser.add_argument(
    "--stop-decade",
    type=float,
    required=False,
    help="Specify the ending decade for data processing."
)
    
parser.add_argument(
    "--points-per-decade",
    type=int,
    required=False,
    help="Specify the number of points per decade for the plot."
)

parser.add_argument(
    "--base-port",
    type=int,
    required=False,
    default=5001,
    help="Specify first IP port for + 0, + 1, + 2 ports."
)

args = parser.parse_args()

# Set the environment variable correctly
os.environ["BASE_PORT"] = str(args.base_port)  # Environment variables must be strings

# --- Helper Function to Ensure Log Directory Exists ---
def ensure_log_directory(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)  # Create the directory if it doesn't exist.
        except Exception as e:
            pass
        
# Logic to set plot type based on environment variable and/or command-line argument
if args.plot_type:
    plot_type = args.plot_type
elif os.environ.get("PLOT_TYPE"):
    plot_type = os.environ["PLOT_TYPE"]
else:
    plot_type = "BODEPLOTTER"

cache_dir_arg=None
# --- Logging Configuration ---
if sys.platform == "win32" or platform.system() == 'Windows':
    # Generate a unique identifier for the instance
    instance_id = str(uuid4())
    # Ensure separate working directories
    # working_dir = f"instance_data_{instance_id}"
    working_dir = os.path.join(tempfile.gettempdir(), f"BodePlotter_{plot_type}_logs_{instance_id}")
    # working_dir = f"C:\\Temp\\instance_data_{instance_id}"
    os.makedirs(working_dir, exist_ok=True)
    
    # Clean up the temp directory upon exiting
    # atexit.register(lambda: shutil.rmtree(working_dir, ignore_errors=True))  # Safe cleanup
    # Register robust_cleanup of working_dir to be last.
    atexit.register(lambda: robust_cleanup(working_dir))
    
    cache_dir_arg=working_dir
    # ensure_log_directory(working_dir)
    # Get the current date and time in a formatted string (e.g., "20250404_1421")
    # represents the date_hours, minutes, and seconds_milliseconds
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")

    # Update log file name to include the timestamp
    log_file = os.path.join(working_dir, f"{plot_type}_{timestamp}_{os.getpid()}_app.log")
    try:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            filename=log_file,
            filemode="w"  # Overwrite each time on startup.
        )
    except Exception as e:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s"
        )
        log_file = None  # Logging file not available

    # Set this flag to True to enable redirection, or False to disable it.
    ENABLE_STDOUT_REDIRECTION = True

    if ENABLE_STDOUT_REDIRECTION:
        # Redirect stdout/stderr only for the parent process.
        try:
            if log_file:
                log_file = os.path.join(working_dir, f"{plot_type}_{timestamp}_{os.getpid()}_app_stdout.log")
                sys.stdout = open(log_file, "w", buffering=1)
                # Register sys.stdout.close before robust_cleanup.
                atexit.register(sys.stdout.close)
                sys.stderr = sys.stdout
        except Exception as e:
            logging.error(f"Failed to redirect stdout/stderr: {e}")
    else:
        # Option 2: Redirection is disabled.
        # This block intentionally does nothing, so stdout and stderr remain unchanged.
        pass
    
    # ...then register logging shutdown:
    atexit.register(logging.shutdown)
   
else:
    # Disable existing handlers
    logging.getLogger().handlers = []

    # Configure basic logging settings (but without a stream handler)
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[logging.NullHandler()]  # Prevent default handlers
    )

    # Create and configure the logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Add a NullHandler to discard logs
    null_handler = logging.NullHandler()
    logger.addHandler(null_handler)


# Suppress verbose logging from Matplotlib.
logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
logging.getLogger("matplotlib.ticker").setLevel(logging.WARNING)

start_decade_arg=args.start_decade
stop_decade_arg=args.stop_decade
points_per_decade_arg=args.points_per_decade
# End moved over from __main__.py

if plot_type == "XY":
    if cache_dir_arg:
        os.environ["MPLCONFIGDIR"] = cache_dir_arg
    logging.info("Launching XYoscilloscope plot...")
    xy_oscilloscope = XYoscilloscope()
    xy_oscilloscope.create_plot(xy_oscilloscope.stop_event, cache_dir_arg)
    if sys.platform == "win32" or platform.system() == 'Windows':
        # Register the shutdown callback so that when the program exits,
        # the plotting process is signaled to stop.
        atexit.register(xy_oscilloscope.send_close_signal)
elif plot_type == "FFT":
    if cache_dir_arg:
        os.environ["MPLCONFIGDIR"] = cache_dir_arg
    logging.info("Launching FFToscilloscope plot...")
    plot_manager_fft = PlotManagerFFToscilloscope()
    plot_manager_fft.create_plot(plot_manager_fft.stop_event, start_decade_arg, stop_decade_arg, points_per_decade_arg, cache_dir_arg)
    if sys.platform == "win32" or platform.system() == 'Windows':
        # Register the shutdown callback so that when the program exits,
        # the plotting process is signaled to stop.
        atexit.register(plot_manager_fft.send_close_signal)
elif plot_type == "MP":
    if cache_dir_arg:
        os.environ["MPLCONFIGDIR"] = cache_dir_arg
    logging.info("Launching MagnitudePhase plot...")
    plot_manager_mp = PlotManagerMagnitudePhase()
    plot_manager_mp.create_plot(plot_manager_mp.stop_event, start_decade_arg, stop_decade_arg, points_per_decade_arg, cache_dir_arg)
    if sys.platform == "win32" or platform.system() == 'Windows':
        # Register the shutdown callback so that when the program exits,
        # the plotting process is signaled to stop.
        atexit.register(plot_manager_mp.send_close_signal)
elif plot_type == "SC":
    if cache_dir_arg:
        os.environ["MPLCONFIGDIR"] = cache_dir_arg
    logging.info("Launching SCoscilloscope plot...")
    plot_manager_sc = PlotManagerSmithChart()
    plot_manager_sc.create_plot(plot_manager_sc.stop_event, cache_dir_arg)
    if sys.platform == "win32" or platform.system() == 'Windows':
        # Register the shutdown callback so that when the program exits,
        # the plotting process is signaled to stop.
        atexit.register(plot_manager_sc.send_close_signal)
else:
    if __name__ == "__main__":
        main()

