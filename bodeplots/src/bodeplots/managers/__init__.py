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
from .PlotManagerFFToscilloscope_classIP import PlotManagerFFToscilloscope
from .PlotManagerMagnitudePhase_classIP import PlotManagerMagnitudePhase
from .XYoscilloscope_classIP import XYoscilloscope
from.PlotManagerSmithChart_classIP import PlotManagerSmithChart

# Optionally, define __all__ to specify what gets imported with "from managers import *"
__all__ = [
    "PlotManagerFFToscilloscope",
    "PlotManagerMagnitudePhase",
    "XYoscilloscope",
    "PlotManagerSmithChart"
]
