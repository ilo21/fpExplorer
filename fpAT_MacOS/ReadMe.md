## Development on MacOS
<br>Tested on Big Sur
- Install [Anaconda](https://www.anaconda.com/products/individual) Python Distribution.
- Create a separate environment for the development.
    - Open Terminal and type:
    ```
    conda create --name FP python=3.7
  ```
- Change change your working environment to the new environment:
```
  conda activate FP
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
- Place fpAT4MacOS.py in the fpAT_src folder and run the app:
```
  python fpAT4MacOS.py
```
