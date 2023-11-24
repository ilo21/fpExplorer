# <div align="center">fpExplorer<br><br>![icon64x64](https://user-images.githubusercontent.com/87764674/174671214-01d6a9e9-39bc-4bd4-8a02-e519a0bd834f.png)<br> [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Created my free logo at LogoMakr.com](https://img.shields.io/badge/Created%20my%20free%20logo%20at-LogoMakr.com-blue)](https://logomakr.com/)
  <br></div>
## What is it?
fpExplorer is an open-source application that allows users to preview and perform basic analysis of fiber photometry data recorded using Tucker Davis Technologies (TDT) systems. To make the application useful to the researchers who do not use TDT system, we added support for generic csv files. <br>This is our attempt to standardize and simplify analysis approaches.<br>
<br>![example_github](https://github.com/ilo21/fpExplorer/assets/87764674/db595b39-110a-413f-9d57-6032ef10e104)<br>
<br>
![batch-z-score-options](https://github.com/ilo21/fpExplorer/assets/87764674/6570e972-dc5a-4d3e-8117-e6c91f756199)
<br>
## Usage
We provide users with a Python source code as well as an independent windows application that does not require any programming background. (See the detailed [Application Guide](https://github.com/ilo21/fpExplorer/blob/main/fpExplorer_src/Documentation/docs.pdf) and [Installation](#installation) instructions.)<br>
<br><br>
![Flow](https://user-images.githubusercontent.com/87764674/174672419-8a7a6296-88f5-40da-a291-fd0218cd0c15.png)
<br>
## Installation
- Windows application: <br>
The installation file for the desktop application is available [here](https://github.com/ilo21/fpExplorer/releases). Download the fpExplorer.Installer.exe file and run it.
## Development
- Install [Anaconda](https://www.anaconda.com/products/individual) Python Distribution.
- Create a separate environment for the development.
  - From [yml](https://github.com/ilo21/fpExplorer/blob/main/environment_info/FPenv.yml) file:<br>Download the repository. Open cmd within environment_info folder where [yml](https://github.com/ilo21/fpExplorer/blob/main/environment_info/FPenv.yml) file is located and type:
  ```
    conda env create -f FPenv.yml
  ```
  - From [specification](https://github.com/ilo21/fpExplorer/blob/main/environment_info/win10FPspec-file.txt) file:<br>Download the repository. Open cmd within environment_info folder where [spec-file](https://github.com/ilo21/fpExplorer/blob/main/environment_info/win10FPspec-file.txt) is located and type:
  ```
    conda create --name FP2 --file win10FPspec-file.txt
  ```
  - From scratch:
      - Open Terminal and type:
      ```
      conda create --name FP2 python=3.7
      ```
      - Change change your working environment to the new environment:
      ```
        conda activate FP2
      ```
      - Install required packeges one by one:
      ```
      pip install tdt
      pip install pandas
      pip install scipy
      pip install sklearn
      pip install matplotlib
      pip install PyQt5
      pip install pyqtgraph
      ```
- Change your working environment to the new environment (if you haven't already):
```
  conda activate FP2
```
- Navigate to fpExplorer_src folder and run:
```
  python fpExplorer.py
```
You can find more information about managing conda environments [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).
## Optional Peri Event Video Extraction
In order to be able to preview recorded signal around our events synchronized together with the recorded videos, we developed a small additional application: [fpVideoExplorer](https://github.com/ilo21/fpExplorer/tree/main/fpVideoExplorer_src). 
## Limitations
This application was developed to meet the needs of researchers at Linköping University (LiU), Sweden. We have tried to make fpExplorer easily accessable (no programming knowledge required) and broadly applicable (as long as datafiles were created with TDT's [Synapse software](https://www.tdt.com/component/synapse-software/)). However, we do realize that each research group might have their own, custom fiber photometry setups and unique experimental designs, and our analysis tool will not be able to handle all of the cases.
- The application was developed and extensively tested on Windows 10 platform. It also works on Windows 11. We have only minimally tested running it on MacOS and Linux systems.
- Our application requires a very proprietary folder structure that is characteristic to [TDT](https://www.tdt.com/docs/synapse/managing-data-for-your-lab/) software (Subject-> Experiment or Experiment->Structure)
- Currently supported events are only the ones that don't start with Cam or Tick (e.g., "PrtA 253", "Note 1")
## Contributors
- [Ilona Szczot](https://liu.se/en/employee/ilosz01) (development and analytical expertise)
- [Joost Wiskerke](https://liu.se/en/employee/joowi80) (development and analytical expertise)
- [David Barker](https://www.thebarkerlab.com/) (we were inspired by his group's [pMAT](https://github.com/djamesbarker/pMAT) (Photometry Modular Analysis Tool))
- [Patrick Mulholland](https://education.musc.edu/MUSCApps/facultydirectory/Mulholland-Patrick) (provided a modified polynomial fitting method)






