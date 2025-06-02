"""

© 2025 by Bode Plotter

MIT License

Copyright (c) [2025] [Bode Plotter]
Bode plot implementation for the handheld oscilloscope OWON HDS320S
Highly rewritten and modified between 01/20/2025 to 06/01/2025
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
import time
import platform
import tempfile
import datetime
import shutil
import atexit
import re
from uuid import uuid4
import argparse
import logging
import psutil  # Ensure psutil is installed (pip install psutil)

# --- Import your plotting managers/applications ---
from bodeplots.managers import PlotManagerFFToscilloscope, PlotManagerMagnitudePhase, XYoscilloscope, PlotManagerSmithChart
from bodeplots.app import main

# --- Check if a log file is in use based on its embedded PID ---
def is_log_file_in_use(file_path):
    """
    Determines if a log file appears to be in use by extracting the process ID from the file name.
    The file name is assumed to follow the format:
      <plot_type>_<timestamp>_<pid>_app.log
    Where the timestamp includes underscores (e.g., "20250404_142100_123456").
    """
    file_name = os.path.basename(file_path)
    pattern = r"^(?P<plot_type>.+?)_(\d{8}_\d{6}_\d+)_(?P<pid>\d+)_app\.log$"
    m = re.match(pattern, file_name)
    if m:
        pid = int(m.group("pid"))
        try:
            return psutil.pid_exists(pid)
        except Exception as e:
            print(f"Error checking PID {pid}: {e}")
            return False
    return False

# --- Cleanup orphan log directories created by previous instances ---
def cleanup_old_log_dirs(dir_prefix="BodePlotter", max_age_seconds=3600):
    """
    Iterates over directories in the temp folder that start with the designated prefix and
    contain '_logs_' in their names. For directories older than max_age_seconds (default 1 hour),
    the function inspects their log files. If none of the log files appear to be in use,
    the directory is removed.
    """
    temp_dir = tempfile.gettempdir()
    now = time.time()
    for entry in os.listdir(temp_dir):
        entry_path = os.path.join(temp_dir, entry)
        if os.path.isdir(entry_path) and entry.startswith(f"{dir_prefix}_") and "_logs_" in entry:
            mod_time = os.path.getmtime(entry_path)
            if now - mod_time > max_age_seconds:
                in_use = False
                # Walk through to locate any log file ending with "_app.log"
                for root, _, files in os.walk(entry_path):
                    for file in files:
                        if file.endswith("_app.log"):
                            log_file_path = os.path.join(root, file)
                            if is_log_file_in_use(log_file_path):
                                in_use = True
                                break
                    if in_use:
                        break
                if not in_use:
                    try:
                        shutil.rmtree(entry_path)
                        print(f"Removed old log directory: {entry_path}")
                    except Exception as e:
                        print(f"Could not remove {entry_path}: {e}")

# --- Robust cleanup function to remove working directories on exit ---
def robust_cleanup(directory):
    """
    Attempt to remove the directory upon exit. Retries a few times before giving up.
    """
    for _ in range(10):
        try:
            shutil.rmtree(directory, ignore_errors=False)
            return
        except Exception as e:
            time.sleep(0.2)
    shutil.rmtree(directory, ignore_errors=True)

# --- Main entry point of the application ---
def main_entry():
    # --- Command-line Argument Parsing ---
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
        help="Specify first IP port for + 0, + 1, + 2, + 3 ports."
    )
    # New flag that controls whether log directories are left intact and, on Linux, logs output to console.
    parser.add_argument(
        "--leave-logs",
        action="store_true",
        default=False,
        help="Leave log files and temp directories intact on exit and, on Linux, output logs to console."
    )
    args = parser.parse_args()

    # Set base port as an environment variable.
    os.environ["BASE_PORT"] = str(args.base_port)

    # Determine the plot type from command-line arguments or environment.
    if args.plot_type:
        plot_type = args.plot_type
    elif os.environ.get("PLOT_TYPE"):
        plot_type = os.environ["PLOT_TYPE"]
    else:
        plot_type = "BODEPLOTTER"

    cache_dir_arg = None

    # --- Windows-Specific Logging Configuration and Cleanup ---
    if sys.platform == "win32" or platform.system() == "Windows":
        # For BODEPLOTTER, clean up old log directories (if any exist and are safe to remove).
        if plot_type == "BODEPLOTTER":
            cleanup_old_log_dirs()

        # Generate a unique identifier for this instance.
        instance_id = str(uuid4())

        # Build the working directory path within the temp folder.
        working_dir = os.path.join(tempfile.gettempdir(), f"BodePlotter_{plot_type}_logs_{instance_id}")
        os.makedirs(working_dir, exist_ok=True)

        # Register exit cleanup if the user hasn't requested to leave logs intact.
        if not args.leave_logs:
            atexit.register(lambda: robust_cleanup(working_dir))

        cache_dir_arg = working_dir

        # Create a timestamp with microsecond precision.
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        log_file = os.path.join(working_dir, f"{plot_type}_{timestamp}_{os.getpid()}_app.log")
        print(f"Logging to: {log_file}")

        try:
            logging.basicConfig(
                level=logging.DEBUG,
                format="%(asctime)s %(levelname)s: %(message)s",
                filename=log_file,
                filemode="w"  # Overwrite the log file each time.
            )
        except Exception as e:
            logging.basicConfig(
                level=logging.DEBUG,
                format="%(asctime)s %(levelname)s: %(message)s"
            )
            log_file = None

        # Optionally redirect stdout/stderr to a dedicated log.
        ENABLE_STDOUT_REDIRECTION = True
        if ENABLE_STDOUT_REDIRECTION:
            try:
                if log_file:
                    stdout_log_file = os.path.join(working_dir, f"{plot_type}_{timestamp}_{os.getpid()}_app_stdout.log")
                    sys.stdout = open(stdout_log_file, "w", buffering=1)
                    atexit.register(sys.stdout.close)
                    sys.stderr = sys.stdout
            except Exception as e:
                logging.error(f"Failed to redirect stdout/stderr: {e}")
        atexit.register(logging.shutdown)
    else:
        # Non-Windows platforms (typically Linux or macOS):
        if args.leave_logs:
            # With --leave-logs, output logs to console.
            logging.basicConfig(
                level=logging.DEBUG,
                format="%(asctime)s %(levelname)s: %(message)s"
            )
        else:
            # Otherwise, use a minimal logging configuration.
            logging.getLogger().handlers = []
            logging.basicConfig(
                level=logging.DEBUG,
                format="%(asctime)s %(levelname)s: %(message)s",
                handlers=[logging.NullHandler()]
            )
            logger = logging.getLogger()
            logger.setLevel(logging.DEBUG)
            logger.addHandler(logging.NullHandler())

    # --- Suppress Verbose Matplotlib Logging ---
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.ticker").setLevel(logging.WARNING)

    # Pass additional arguments to plotting routines.
    start_decade_arg = args.start_decade
    stop_decade_arg = args.stop_decade
    points_per_decade_arg = args.points_per_decade

    # --- Dispatch to the Appropriate Plotting Manager or Main App ---
    if plot_type == "XY":
        if cache_dir_arg:
            os.environ["MPLCONFIGDIR"] = cache_dir_arg
        logging.info("Launching XYoscilloscope plot...")
        xy_oscilloscope = XYoscilloscope()
        xy_oscilloscope.create_plot(xy_oscilloscope.stop_event, cache_dir_arg)
        if sys.platform == "win32" or platform.system() == "Windows":
            atexit.register(xy_oscilloscope.send_close_signal)
    elif plot_type == "FFT":
        if cache_dir_arg:
            os.environ["MPLCONFIGDIR"] = cache_dir_arg
        logging.info("Launching FFToscilloscope plot...")
        plot_manager_fft = PlotManagerFFToscilloscope()
        plot_manager_fft.create_plot(
            plot_manager_fft.stop_event,
            start_decade_arg,
            stop_decade_arg,
            points_per_decade_arg,
            cache_dir_arg
        )
        if sys.platform == "win32" or platform.system() == "Windows":
            atexit.register(plot_manager_fft.send_close_signal)
    elif plot_type == "MP":
        if cache_dir_arg:
            os.environ["MPLCONFIGDIR"] = cache_dir_arg
        logging.info("Launching MagnitudePhase plot...")
        plot_manager_mp = PlotManagerMagnitudePhase()
        plot_manager_mp.create_plot(
            plot_manager_mp.stop_event,
            start_decade_arg,
            stop_decade_arg,
            points_per_decade_arg,
            cache_dir_arg
        )
        if sys.platform == "win32" or platform.system() == "Windows":
            atexit.register(plot_manager_mp.send_close_signal)
    elif plot_type == "SC":
        if cache_dir_arg:
            os.environ["MPLCONFIGDIR"] = cache_dir_arg
        logging.info("Launching SCoscilloscope plot...")
        plot_manager_sc = PlotManagerSmithChart()
        plot_manager_sc.create_plot(plot_manager_sc.stop_event, cache_dir_arg)
        if sys.platform == "win32" or platform.system() == "Windows":
            atexit.register(plot_manager_sc.send_close_signal)
    else:
        if __name__ == "__main__":
            main()

if __name__ == "__main__":
    main_entry()

