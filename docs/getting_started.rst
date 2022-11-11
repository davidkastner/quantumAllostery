Getting Started
===============

Welcome to quantumAllostery!
A workflow for identify quantum allostery trends in molecular dynamics data.

Installation
------------
quantumAllostery is build as a library for automating data generation and analysis.

::

    $ conda create --name qa
    $ conda activate qa

Install dependencies
--------------------
::

    $ conda install -c anaconda requests
    $ conda install -c anaconda beautifulsoup4
    $ conda install -c anaconda pandas
    $ conda install -c plotly plotly_express
    $ pip install sphinx sphinx_rtd_theme
    $ conda install -c conda-forge sphinx-autoapi
    $ pip install https://github.com/revitron/revitron-sphinx-theme/archive/master.zip

Install package
---------------
::

    $ pip install -e .

Run a pipeline
--------------
::
    
    $ python pipelines.py