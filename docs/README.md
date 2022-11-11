# quantumAllostery's Complete Build Documentation (Devs)

## Creating python environment
```bash
conda create -n qa
conda activate qa
```

## Install useful packages
```bash
conda create --name package_name
conda activate package_name
conda install numpy
conda install -c conda-forge matplotlib
conda install -c conda-forge scipy
conda install cookiecutter
conda install -c conda-forge git
conda install -c conda-forge black
conda install -c conda-forge mypy
conda install -c conda-forge flake8
conda install -c conda-forge pytest
conda install -c conda-forge pydantic
pip install pytest-cov
pip install sphinx sphinx_rtd_theme
pip install https://github.com/revitron/revitron-sphinx-theme/archive/master.zip
```

## Build from MOLSSI scientific standards
```bash
cookiecutter gh:molssi/cookiecutter-cms
```

## Setup developing environment
```bash
cd quantumAllostery
pip install -e .
```

## Compile static HTML pages
```bash
make html
```
