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
from collections import deque
import numpy as np
import multiprocessing
import queue
import threading



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
        
        # If initial evaluation period is not complete, wait (2025/03/18 only wait for initial evaluation period only where self.sign_locked is False)
        if len(self.phase_diff_values) < self.initial_evaluation_count and not self.sign_locked:
            # logging.info(f"Value: {value}, Waiting for initial evaluation period to complete.")
            return None
        
        # After initial evaluation period, lock the sign based on majority sign
        if len(self.phase_diff_values) == self.initial_evaluation_count and not self.sign_locked:
            self.locked_sign = self.majority_sign()
            self.sign_locked = True
            # logging.info(f"Initial evaluation completed. Locked Sign: {'Positive' if self.locked_sign > 0 else 'Negative'}")
        
        # If sign is locked, check if the incoming value is outside the threshold
        if self.sign_locked:
            if abs(value) > self.threshold:
                # logging.info(f"Value: {value}, Outside threshold. Locked Sign remains: {'Positive' if self.locked_sign > 0 else 'Negative'}")
                # 2025/03/17 clear phase_diff_values after lock so once value are within threshold new evaluation will be made.
                self.phase_diff_values.clear()
                return self.locked_sign
            else:
                moving_avg = self.weighted_moving_average(self.phase_diff_values)
                self.locked_sign = 1 if moving_avg >= 0 else -1
                # logging.info(f"Value: {value}, Within threshold. Updated Locked Sign: {'Positive' if self.locked_sign > 0 else 'Negative'}")

        return self.locked_sign
    
    def cleanup(self):
        # Reset the processor attributes to their initial state
        self.phase_diff_values.clear()
        self.sign_locked = False
        self.locked_sign = None
        # logging.info("Processor cleaned up and ready for new data.")

# --- End class PhaseDiffProcessor:--- -----------------------------------------------
            
