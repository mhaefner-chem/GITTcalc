# GITTcalc
![GitHub Release](https://img.shields.io/github/v/release/mhaefner-chem/GITT_Analysis?include_prereleases) ![GitHub License](https://img.shields.io/github/license/mhaefner-chem/GITT_Analysis) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17986238.svg)](https://doi.org/10.5281/zenodo.17986238)

GITTcalc analyzes raw GITT data to extract the diffusion coefficients of the conducting ion based on the mini-review 'Principle and Applications of Galvanostatic Intermittent Titration Technique for Lithium-ion Batteries' by Jaeyoung Kim, Sangbin Park, Sunhyun Hwang, and Won-Sub Yoon. (DOI: https://doi.org/10.33961/jecst.2021.00836) and equation 16, in particular. It can either be used as a standalone program or inside OriginLab.

## Table of contents

- [How to use GITTcalc?](#how-to-use-gittcalc)
- [Requirements & Installation](#requirements-and-installation)

## How to use GITTcalc?
GITT_Analysis processes raw GITT data to obtain diffusion coefficients. For this, a file with the time-voltage-pairs from the measurement are required, as well as the area-normed mass of the active material in g/cm², the molar mass of the active material in g/mol, the density of the active material in g/cm³, and the contact area with the electrode during the measurement in cm². Additionally, the diffusion coefficents at different ion contents can be calculated if either capacity or specific capacity is provided in the same file as the time and voltage. This also requires a theoretical capacity for a hypothetical ion content of 1 (e.g., Li<sub>1</sub>NiO<sub>2</sub> for Li<sub>x</sub>NiO<sub>2</sub> or Na<sub>1</sub>CoO<sub>2</sub> for Na<sub>x</sub>CoO<sub>2</sub>) and the starting ion content.

If this program is used inside OriginLab, the raw and processed data are automatically output into a Workbook for further processing. If this program is used on its own, the processed data can be saved as a CSV-file to a location of the user's choosing. In either case, the extra properties for the active material and measurement are saved as INFO-file in plain text at the same location as the raw data to store them in case the data needs to be processed again at a later point.

## Requirements and installation

So far, the program has been successfully tested with `python 3.10` and `OriginPro 2024b 10.1.5.132`.
The source code for the program is obtained with the command

```console
$ git clone https://github.com/mhaefner-chem/GITTcalc
```

Running the standalone program with python requires the python packages `numpy`, `scipy`, `matplotlib` and all their dependencies. They can be installed via pip

```console
$ pip install numpy scipy matplotlib
```
or via conda as

```console
$ conda install -c conda-forge numpy scipy matplotlib
```

If run inside the OriginLab command window with

```console
run -pyf GITTcalc.py
```
the packages `numpy`, `matplotlib`, and `scipy` need to be installed via OriginLab's native package manager. For further information on how to install python packages in OriginLab, please consult their [website](https://www.originlab.com/doc/python/Python-Packages).


