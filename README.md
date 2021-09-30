# <div align="center">fpAT<br><br>![icon64x64](https://user-images.githubusercontent.com/87764674/126611941-34d2f3c2-0f24-4517-a82a-4c01ac82599f.png)<br>(Fiber Photometry Analysis Tool)<br> [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)<br></div>
## What is it?
fpAT is an open-source application that allows users to preview and perform basic analysis of fiber photometry data recorded using Tucker Davis Technologies (TDT) systems. This is our attempt to standardize and simplify analysis approaches.<br>
<br>![visualize](https://user-images.githubusercontent.com/87764674/134684805-b8c2564b-4e95-4837-a278-9ed875adcb3e.PNG)
<br>
![batch-z-score-options](https://user-images.githubusercontent.com/87764674/134880894-ba74ccb2-c06d-445f-8d69-3ae24cbae589.PNG)
<br>
## Usage
We provide users with a Python source code as well as an independent windows application that does not require any programming background. (See the detailed [Application Guide](https://github.com/ilo21/fpAT/blob/main/fpAT_src/Documentation/docs.pdf) and [Installation](#installation) instructions.)<br>
<br><br>
![Flow drawio](https://user-images.githubusercontent.com/87764674/134685031-f346e347-2bb2-498e-b3c9-a205b7d8b93e.png)
<br>
## Installation
- Windows application: <br>
The installation file for the desktop application is available [here](https://github.com/ilo21/fpAT/releases). Download the fpATinstaller.exe file and run it.
## Development
- Install [Anaconda](https://www.anaconda.com/products/individual) Python Distribution.
- Create a separate environment for the development.
  - From [yml](https://github.com/ilo21/fpAT/blob/main/environment_info/FPenv.yml) file.<br>Download the repository. Open cmd within environment_info folder where [yml](https://github.com/ilo21/fpAT/blob/main/environment_info/FPenv.yml) file is located and type:
  ```
    conda env create -f FPenv.yml
  ```
  - From [specification](https://github.com/ilo21/fpAT/blob/main/environment_info/win10FPspec-file.txt) file<br>Download the repository. Open cmd within environment_info folder where [spec-file](https://github.com/ilo21/fpAT/blob/main/environment_info/win10FPspec-file.txt) is located and type:
  ```
    conda create --name FP --file win10FPspec-file.txt
  ```
- Lastly, change your working environment to the new environment:
```
  conda activate FP
```
- Navigate to fpAT_src folder and run:
```
  python fpAT.py
```
You can find more information about managing conda environments [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).
<br><br>*[Recent adjustments for the development on MacOS](https://github.com/ilo21/fpAT/tree/main/fpAT_MacOS)*<br>
## Limitations
This application was developed to meet the needs of researchers at Linköping University (LiU), Sweden. We have tried to make fpAT easily accessable (no programming knowledge required) and broadly applicable (as long as datafiles were created with TDT's Synapse software). However, we do realize that each research group might have their own, custom fiber photometry setups and unique experimental designs, and our analysis tool will not be able to handle all of the cases.
- The application was developed and tested only on Windows 10 platform.
- Our application requires a very proprietary folder structure that is characteristic to [TDT](https://www.tdt.com/docs/synapse/managing-data-for-your-lab/) software (Subject-> Experiment or Experiment->Structure)
- Currently supported events are only the ones that start with Prt or Note (e.g., "PrtA 253", "Note 1")
## Contributors
- Ilona Szczot (development and analytical expertise)
- Joost Wiskerke (development and analytical expertise)
- [David Barker](https://www.thebarkerlab.com/) (we were inspired by his group's [pMAT](https://github.com/djamesbarker/pMAT) (Photometry Modular Analysis Tool))
- Patrick Mulholland (provided a modified polynomial fitting method)









