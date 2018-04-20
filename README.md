# Combining range and direction for improved localizatioin

This repository contains the code presented at the IEEE International Conference on Acoustic Speech and Signal Processing 2018 in Calgary CA. 

Authors: Gilles Baechler', Frederike Duembgen', Golnooshsadat Elhami', Miranda Krekovic', Robin Scheibler, Adam Scholefield and Martin Vetterli. 

' equal contributions. 

Last modified: 17 April 2018

## Overview

This repository contains the following components: 
 
- Plots.ipynb, the ipython notebook used to reproduce all plots used in the paper. 
- run_simulation.py, the code used to run your own simulations. 
- test_run_simulation.py, some very simple tests to make sure the code runs.

## Requirements

This repository relies on a python toolbox called _pylocus_, which contains implementations of various angle- and distance based localization algorithms.   

The source code and installation instructions can be found at [https:/github.com/LCAV/pylocus](https:/github.com/LCAV/pylocus). 

The results were created using the release 1.0.1 of the package. If you use pip, you can install it and other required packages by running

```bash
pip install -r requirements.txt
``` 

## Contribute

Feel free to create issues or pull requests, or send an e-mail to frederike.duembgen@epfl.ch if you run into any difficulties while using this code base. 
