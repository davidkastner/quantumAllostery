Getting Started
===============

*Welcome to quantumAllostery!*

1 Overview
----------

The quantumAllostery package was designed to automate the identification of charge-transfer events in ab-initio molecular dynamics simulations. Routine tasks can be easily automated using the functionality contained in the library.

2 Installation
--------------

Install the package by running the follow commands inside the repository. This will perform a developmental version install. It is good practice to do this inside of a virtual environment. A yaml environmental file has been created to automate the installation of dependencies.

**Creating a Conda environment**

All the dependencies can be loaded together using the prebuilt environment.yml files.
We provide two YAML files. The dev version contains additional packages for code maintenance.

::

    conda env create -f environment.yml

If you are going to be developing the package run:

::

    conda env create -f environment_dev.yml

**Setup developing environment**

::

    git clone git@github.com:davidkastner/quantumAllostery.git
    cd quantumAllostery
    pip install -e .


3 What is included?
-------------------

::
    
    .
    ├── docs
    ├── qa
    │   ├── qa           # Top-level script that interacts with the rest of the package
    │   ├── process      # Processes the raw AIMD data
    │   ├── predict      # Machine learning analysis
    │   └── plot         # Automated plotting and vizualization 
    └── ...


4 Documentation
---------------

Update ReadtheDocs.

::

    cd docs
    make clean
    make html

5 Developer guide
-----------------

**GitHub Refresher**

::

    git status
    git add .
    git commit -m "Change a specific functionality"
    git push -u origin main

