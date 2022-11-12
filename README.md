quantumAllostery
==============================
[//]: # (Badges)
[![Documentation Status](https://readthedocs.org/projects/quantumallostery/badge/?version=latest)](https://quantumallostery.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/quantumAllostery/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/quantumAllostery/branch/main)


# quantumAllostery Package
## Table of Contents
1. **Overview**
    * Introduction
    * Purpose
2. **Installation**
    * Installing quantumAllostery
    * Prerequisites
    * File structure
3. **What is included?**
    * Library
    * Functionality
4. **Documentation**
    * Read the Docs
    * Examples
5. **Developer Guide**
    * GitHub refresher


## Overview
The quantumAllostery package was designed to automate the identification of charge-transfer events in ab-initio molecular dynamics simulations. Routine tasks can be easily automated using the functionality contained in the library.


## Installation
Install the package by running the follow commands inside the repository. This will perform a developmental version install. It is good practice to do this inside of a virtual environment. A yaml environmental file has been created to automate the installation of dependencies.

### Creating python environment
```bash
conda create -n qa
conda activate qa
conda env create -f environment.yml
```

### Setup developing environment
Remember to update your GitHub [ssh keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
```bash
git clone git@github.com:davidkastner/quantumAllostery.git
cd quantumAllostery
pip install -e .
```

## What is included?
### File structure
```
.
├── docs
├── qa
│   ├── qa           # Top-level script that interacts with the rest of the package
│   ├── process      # Processes the raw AIMD data
│   ├── predict      # Machine learning analysis
│   └── plot         # Automated plotting and vizualization 
└── ...
```


## Documentation
### Update the docs
```bash
make clean
make html
```


## Developer guide
### GitHub refresher
```bash
git status
git add .
git commit -m "Change a specific functionality"
git push -u origin main
```


### Copyright
Copyright (c) 2022, David W. Kastner


#### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
