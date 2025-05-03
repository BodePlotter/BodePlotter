
---

## Frequency Response (Bode Plot)

---

For installer downloads and source code, please visit our GitHub Pages site:  
[BodePlotter](https://github.com/BodePlotter/BodePlotter)  
You can also browse the latest releases on GitHub:  
[Releases](https://github.com/BodePlotter/BodePlotter/releases)

---

## Table of Contents
1. [Features](#features)
2. [Example RC Band Pass Filter Circuit](#example-rc-band-pass-filter-circuit)
3. [Screen shots of bodeplots](#screen-shots-of-bodeplots)
4. [Software Usage](#software-usage)
5. [Command Line Options](#command-line-options)
6. [Recommended Display](#recommended-display)
7. [Oscilloscope Setup](#oscilloscope-setup)
8. [Install on Windows 11](#install-on-windows-11)
9. [Install USB Driver for Windows 11](#install-usb-driver-for-windows-11)
10. [Install the .deb Package on Mint Linux (22.1)](#install-the-deb-package-on-mint-linux-221)
11. [Remove the .deb Package on Mint Linux (22.1)](#remove-the-deb-package-on-mint-linux-221)
12. [Build the .deb Package](#build-the-deb-package)
13. [Windows Build Process](#windows-build-process)
14. [Issues](#issues)

---

## Features
1.  **Python Application with Briefcase**:
   - Generate Bode-Plots.
   - Accumulate Peak FFT values.
   - Display oscilloscope.
   - Display X-Y mode or Lissajous patterns.
   - Current implementation uses the OWON HDS320S.
2.  **OWON HDS320S Oscilloscope + Multi-Meter + Waveform Generator**:
   - Bandwidth: 200 MHz, Sample Rate: 1 GS/s.
   - Sine Wave Generator: 0.1 Hz to 30 MHz.
   - Outputs 300 samples from the screen for both channels via SCPI.
3.  **Additional Information**:
   - [OWON HDS300 Series Digital Oscilloscope](https://www.owonna.com/products_owon_hds300_series_digital_oscilloscope)
   - [SCPI Interface for OWON HDS320S](http://files.owon.com.cn/software/Application/HDS200_Series_SCPI_Protocol.pdf)

---

## Example RC Band Pass Filter Circuit

- This example uses a passive RC Band Pass Filter.
- Measured values of the RC band pass filter circuit using the HDS320S multi-meter function:
  - **R1**: Resistance of the lower cutoff resistor = 1.46kΩ
  - **C1**: Capacitance of the lower cutoff capacitor = 105.89nF
  - **R2**: Resistance of the upper cutoff resistor = 2.93kΩ
  - **C2**: Capacitance of the upper cutoff capacitor = 1.02nF

---

### Lower Cutoff Frequency (f1)

```math
f1 = 1 / (2 * π * R1 * C1)

Given:
R1 = 1.46 [kΩ] = 1460 [Ω]
C1 = 105.89 [nF] = 105.89 × 10^-9 [F]

f1 = 1 / (2 * π * 1460 * 105.89 × 10^-9) ≈ 1029.5 [Hz]
```

---

### Upper Cutoff Frequency (f2)

```math
f2 = 1 / (2 * π * R2 * C2)

Given:
R2 = 2.93 [kΩ] = 2930 [Ω]
C2 = 1.02 [nF] = 1.02 × 10^-9 [F]

f2 = 1 / (2 * π * 2930 * 1.02 × 10^-9) ≈ 53.254 [kHz]
```

---

## TINA-TIV9 Simulation of the Circuit

The circuit was simulated using TINA-TIV9 software. Below is the visual representation of the circuit setup:

![TINA-TIV9 Circuit Simulation](./media/image1.png)

---

### TINA-TIV9 Simulation Results

The simulation provided the following results, showcasing the circuit's behavior:

![TINA-TIV9 Simulation Results](./media/image2.jpeg)

---

### Bode Plots Using OWON HDS320S

The Bode plots were generated using the OWON HDS320S oscilloscope in conjunction with the Bode-Plots Python application:

![Bode Plots with OWON HDS320S](./media/image3.png)

---

### Screen shots of bodeplots

### HDS320S Magnitude/Phase Bode Plotter window before Measurements or playback

![HDS320S Magnitude/Phase Bode Plotter window before Measurements or playback](./media/Screenshot_001.png)

---

### HDS320S Magnitude/Phase Bode Plotter window Playback Dark Theme

![HDS320S Magnitude/Phase Bode Plotter window Playback Dark Theme](./media/Screenshot_002.png)

---

### HDS320S Magnitude/Phase Bode Plotter window Playback Light Theme

![HDS320S Magnitude/Phase Bode Plotter window Playback Light Theme](./media/Screenshot_003.png)

---

### HDS320S Magnitude/Phase Bode Plotter window 1MHz to 10MHz getting into the noise floor

![HDS320S Magnitude/Phase Bode Plotter window 1MHz to 10MHz getting into the noise floor](./media/Screenshot_004.png)

---

### Watch the HDS320S Magnitude/Phase Bode Plotter window and MatPlotLib windows Playback Video

[Watch the HDS320S Magnitude/Phase Bode Plotter window and MatPlotLib windows Video](./media/HDS320S-bode-plots-0.0.3.mp4)

---

## Software Usage
1. Launch the software via the command line (`bodeplots`) or from the menu.
2. Verify oscilloscope settings and input desired parameters.
3. Click **Search and Setup Oscilloscope** to initialize.
4. Click **Start Measurements** (takes about 15 to 30 minutes).
5. After completion, load the JSON log file for playback and analysis.

---

## Command Line Options

- All inputs are used for Windows 11 process spawning except for `--base-port`.
- If the `--base-port` option is not specified, local ports 5001, 5002, and 5003 will be used to communicate with the Matplotlib windows.

```bash
usage: bodeplots [-h] [--plot-type {XY,MP,FFT,BODEPLOTTER}] 
                      [--start-decade START_DECADE] 
                      [--stop-decade STOP_DECADE] 
                      [--points-per-decade POINTS_PER_DECADE] 
                      [--base-port BASE_PORT]

BodePlots Application

options:
  -h, --help            show this help message and exit
  --plot-type {XY,MP,FFT,BODEPLOTTER}
                        Select the type of plot to display. Options: XY, MP, FFT, or BODEPLOTTER.
  --start-decade START_DECADE
                        Specify the starting decade for data processing.
  --stop-decade STOP_DECADE
                        Specify the ending decade for data processing.
  --points-per-decade POINTS_PER_DECADE
                        Specify the number of points per decade for the plot.
  --base-port BASE_PORT
                        Specify the first IP port for communication; subsequent ports (+1, +2) will be used automatically.

````

---

## Recommended Display
1. 4K Display (4096x2160) preferred; minimum 1600x1200 resolution.
2. VirtualBox with Mint Linux 22.1 works well.

---

## Oscilloscope Setup
1. Connect the oscilloscope to your PC using a USB cable (preferably with Ferrite Rings to reduce noise).
2. Set HDS320S to default:
   - Press **System** → **F3 (Default Setting)** → Confirm.
3. Connect probes:
   - CH1 (yellow lead) to circuit input.
   - CH2 (blue lead) to circuit output.
   - Gen Out to input ground.

---

## Install on Windows 11
1. Download [Bode-Plots-0.0.3.msi](https://github.com/BodePlotter/BodePlotter/releases/download/0.0.3/Bode-Plots-0.0.3.msi) and double-click it in File Explorer.
2. Follow the installation prompts: **Next** → **Next** → **Finish**.
3. Open "Bode" from the Windows Start menu.
4. Adjust vertical height as needed.
5. [Install USB Driver for Windows 11](#install-usb-driver-for-windows-11) using Zadig.

---

## Install USB Driver for Windows 11

1. **Download and Install Zadig**  
   Download [zadig-2.9.exe](https://github.com/pbatard/libwdi/releases) to install the drivers for the OWON HDS320S.

2. **Allow Device Changes**  
   When prompted, allow Zadig to change Windows 11 device settings by clicking the **Yes** button.  
   ![Allow Zadig to change Windows 11 device settings](./media/zadig-001.png)

3. **Allow Updates**  
   When prompted, allow Zadig to check for updates by clicking the **Yes** button.  
   ![Allow Zadig to check for updates](./media/zadig-002.png)

4. **Detect USB Device**  
   Ensure your HDS320S is plugged in and powered on. Zadig should detect the USB device otherwise select from list.  
   ![Zadig detects USB device](./media/zadig-003.png)

4. **Detect USB Device**  
   Correctly selected USB device will have USB ID 5345 1234.  
   ![Zadig detects USB device](./media/zadig-004.png)

5. **Select the Correct Driver**  
   In the driver details, rename the device to **OWON HDS320S**. Then, use the drop-down arrow to select **libusb-win32 (v1.4.0.0)**.  
   ![Select libusb-win32 driver](./media/zadig-005.png)

6. **Install the Driver**  
   Click the **Install Driver** button and wait for the installation to complete.  
   ![Install Driver](./media/zadig-006.png)

7. **Successful Installation**  
   Zadig should report that the driver was installed successfully.  
   ![Driver installed successfully](./media/zadig-007.png)

---

## Install the .deb Package on Mint Linux (22.1)
1. Download: [bodeplots_0.0.3-1~linuxmint-xia_amd64.deb](https://github.com/BodePlotter/BodePlotter/releases/download/0.0.3/bodeplots_0.0.3-1.linuxmint-xia_amd64.deb).
2. Open a terminal and run:
   ```bash
   sudo apt install ./bodeplots_0.0.3-1~linuxmint-xia_amd64.deb
   ```

---

## Remove the .deb Package on Mint Linux (22.1)
1. Open a terminal and run:
   ```bash
   sudo apt purge bodeplots
   ```
---

## Build the .deb Package
1. Use Live Mint (22.1) via VirtualBox or a local installation.
2. Download: [bodeplots-source-code.tgz](https://github.com/BodePlotter/BodePlotter/releases/download/0.0.3/bodeplots-source-code.tgz).
3. Follow these steps:
   ```bash
   mkdir BodePlotter-ScopeFFT
   mv ~/Downloads/bodeplots-source-code.tgz BodePlotter-ScopeFFT
   cd BodePlotter-ScopeFFT
   apt update
   sudo apt install python3-virtualenv git
   virtualenv bodeplotter -p python
   source bodeplotter/bin/activate
   pip install --upgrade pip
   pip install briefcase numpy dearpygui matplotlib pyusb scipy screeninfo PyQt5
   tar zpxvf bodeplots-source-code.tgz
   cd bodeplots
   rm -rf build/bodeplots
   briefcase create
   briefcase update -r
   briefcase build
   briefcase package
   ```
   Package file location: `dist/bodeplots_0.0.3-1~linuxmint-xia_amd64.deb`.

---

## Windows Build Process
### Environment Setup for Windows 11
1.  **Install Latest Microsoft Visual C++ Redistributable**:
   - [Documentation](https://learn.microsoft.com/en-gb/cpp/windows/latest-supported-vc-redist?view=msvc-170#latest-microsoft-visual-c-redistributable-version).
   - [Download vc_redist.x64.exe](https://aka.ms/vs/17/release/vc_redist.x64.exe).
2.  **Install Git**:
   - [Git Download Page](https://git-scm.com/downloads/win).
   - [Direct Download](https://github.com/git-for-windows/git/releases/download/v2.48.1.windows.1/Git-2.48.1-64-bit.exe).
3.  Install Python and dependencies via PowerShell:
    ```PowerShell
    mkdir C:\BodePlotter-ScopeFFT
    cd C:\BodePlotter-ScopeFFT
    python # or https://www.python.org/downloads/release/python-3132/
    # https://www.python.org/ftp/python/3.13.2/python-3.13.2-amd64.exe
    python -m venv bodeplotter
    .\bodeplotter\Scripts\activate
    python.exe -m pip install --upgrade pip
    pip install briefcase numpy dearpygui matplotlib pyusb scipy screeninfo PyQt5
    ```
4.  Extract source code and build: [bodeplots-source-code.tgz](https://github.com/BodePlotter/BodePlotter/releases/download/0.0.3/bodeplots-source-code.tgz)
    ```PowerShell
    cd C:\BodePlotter-ScopeFFT
    .\bodeplotter\Scripts\activate
    cd .\bodeplots
    briefcase create windows app
    briefcase update -r
    briefcase update
    briefcase build -v
    briefcase run -u
    briefcase package windows
    ```
5.  Packaged MSI: `dist\Bode-Plots-0.0.3.msi`.

---

## Issues

1. **DearPyGui Maximize Button**  
   The maximize button in DearPyGui does not function correctly due to fixed window size constraints.

2. **Matplotlib Performance on Windows 11**  
   On Windows 11, Matplotlib may struggle to keep pace at higher playback speeds.
   - If Matplotlib is selected, starting measurements or playing may take more than 15 seconds for additional windows load.
     You may want to click the Pause button until the windows are loaded after seeing the first data points.
   - *Matplotlib windows support was added in version 0.0.2.*
   - Matplotlib windows are not closed automatically before DearPyGui exits, leaving log and temporary files in the user’s TEMP directory.
   - **Command Prompt:**  
     ```cmd
     explorer %TEMP%
     ```
   - **PowerShell:**  
     ```powershell
     explorer $env:TEMP
     ```
   - [Bug #55 with briefcase-windows-VisualStudio-template](https://github.com/beeware/briefcase-windows-VisualStudio-template/issues/55) that causes Matplotlib to fail due to I/O handling issues. This happens as a result of the Briefcase logging wrapper does not support handling multiple instances due to spawning Matplotlib windows. Note, using Matplotlib within threads makes DearPyGui unstable and therefore starting a new process was used in place of threads. The new process in windows results in spawning a new instance.
   - **Workaround:**  
     Using the `console_app` flag in `pyproject.toml` to open a blank console appears to disable logging, mitigating these issues.
     The advantage of the console is you can do ctrl+c to close all Bode-Plots windows. Recommend Stop first for active measurements or Play back.
	 If you would like to hid the console recommend making a bode-plots.vbs.
	 If issues running with Invalid character use Notepad++, go to Encoding -> Convert to ANSI and then save the file.
	 ```cmd
	 Set WshSell = CreateObject("WScript.Shell")
     WshSell.Run "cmd /c bodeplots.exe", 0, True
     Set WshSell = Nothing
	 ```
	 
	 Without the `console_app` flag to get Matplotlib windows to open. (Allow spawned instances or more than on cocurrent instance.) 
     ```powershell
     cd C:\BodePlotter-ScopeFFT
     .\bodeplotter\Scripts\activate
     briefcase.exe -v
     cd "$HOME\AppData\Local\Programs\Bode Plotter\Bode-Plots"
     .\Bode-Plots.exe
     ```

---

