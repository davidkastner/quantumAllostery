![Graphical Summary of README](docs/_static/header.webp)
quantumAllostery
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/davidkastner/quantumAllostery/workflows/CI/badge.svg)](https://github.com/davidkastner/quantumAllostery/actions?query=workflow%3ACI)
[![Documentation Status](https://readthedocs.org/projects/quantumallostery/badge/?version=latest)](https://quantumallostery.readthedocs.io/en/latest/?badge=latest)


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


## 1. Overview
The quantumAllostery package was designed to automate the identification of charge-transfer events in ab-initio molecular dynamics (AIMD) simulations. Routine tasks can be easily automated using the functionality contained in the library.


## 2. Installation
Install the package by running the follow commands inside the repository. This will perform a developmental version install. It is good practice to do this inside of a virtual environment. A yaml environmental file has been created to automate the installation of dependencies.

### Creating python environment
All the dependencies can be loaded together using the prebuilt environment.yml or environment_dev.yml files.
We provide two YAML files. The dev version contains additional packages for code maintenance.
If you are only going to be using the package run:
```bash
conda env create -f environment.yml
```
If you are going to be developing the package run:
```bash
conda env create -f environment_dev.yml
```

### Setup developing environment
Remember to update your GitHub [ssh keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
```bash
git clone git@github.com:davidkastner/quantumAllostery.git
cd quantumAllostery
pip install -e .
```


## 3. What is included?
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


## 4. Documentation
### Run the following commands to update the ReadTheDocs
```bash
make clean
make html
```


## 5. Developer guide
### GitHub refresher for those who would like to contribute
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
