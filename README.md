# Combining range and angles for improved localizatioin

This repository contains the code presented at the IEEE International Conference on Acoustic Speech and Signal Processing 2018 in Calgary CA. 

Last modified: 17 April 2018

## Authors 

This work is a collaboration of a research group of the Laboratory of Audiovisual Communications (LCAV), EPFL, Switzerland. The first 4 authors have contributed equally to this research. 

- Gilles Baechler
- Frederike Duembgen
- Miranda Krekovic
- Golnooshsadat Elhami
- Robin Scheibler
- Adam Scholefield
- Martin Vetterli

## Overview

This repository contains the following components: 
 
- Plots.ipynb, the ipython notebook used to reproduce all plots used in the paper. 
- run_simulation.py, the code used to run your own simulations. 
- src/ contains all the modules created for the described methods. 

## Requirements

This repository relies on a python toolbox called _pylocus_, which contains implementations of various angle- and distance based localization algorithms.   

The source code and installation instructions can be found at [https:/github.com/lcav/pylocus]. 

The results were created using the release 0.0.1 of the package. If you use pip, you can install it and other required packages by running, from this folder level, 

```bash
pip install -r requirements.txt
``` 

## Contribute

Feel free to create issues or pull requests, or send an e-mail to frederike.duembgen@epfl.ch if you run into any difficulties. 



