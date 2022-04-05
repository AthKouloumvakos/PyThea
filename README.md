# PyThea: A software package to reconstruct the 3D structure of CMEs and shock waves

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Version](https://img.shields.io/github/v/release/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases)
[![Release Date](https://img.shields.io/github/release-date/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![flake8](https://github.com/AthKouloumvakos/PyThea/actions/workflows/flake8.yml/badge.svg)
![pytest](https://github.com/AthKouloumvakos/PyThea/actions/workflows/pytest.yml/badge.svg)

![Logo](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_logo.png)

_PyThea_ is an open-source software package that can be used to reconstruct the 3D structure of Coronal Mass Ejections (CMEs) and shock waves and determine their kinematics using remote-sensing observations. The tool implements the Graduated Cylindrical Shell (GCS) model that can be used to reconstruct CMEs and two geometrical models, namely a spheroid and ellipsoid model to reconstruct shock waves. It also implements remote-sensing observations from multiple viewpoints such as the Solar and Heliospheric Observatory (SoHO) and Solar Terrestrial Relations Observatory (STEREO). An online preview of this tool is available at [https://athkouloumvakos.github.io/PyThea](https://athkouloumvakos.github.io/PyThea).

## 💾 Installation

_PyThea_ is written in Python >=3.8 and has some package requirements to run. These are listed in the requirements.txt and environment.yml files located in the main repository. To run this application locally, we recommend to create its own virtual environment in Python and install the dependent packages. Because of a range of dependencies that the different packages have, the simplest way to work with _PyThea_ is installing the packages with _conda_ and creating its own environment. Alternatively, you can create a virtual environment in Python inside the project folder and install the packages using ```pip```. Both ways are explained below.

First, download the latest release of _PyThea_ and unzip the file. Navigate (e.g., ```cd```) to the root directory of _PyThea_ which is the unzipped folder and follow the instructions to create the virtual environment and install the dependent packages:

**Recommended (conda)**

We create the virtual environment (see [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)) and install the dependent packages by doing the following in the terminal,

```python
conda env create -f environment.yml
conda info --envs
```

Then every time before using _PyThea_, you have to activate the environment and when finishing your work deactivate it using the following commands,

```python
# Activate the enviroment
conda activate PyThea

# Here you run _PyThea_ (see Run locally section)

# When you are done you can deactivate a virtual environment
conda deactivate
```

**Alternative (pip)**

You can create a virtual environment in Python inside the _PyThea_ project (root) folder using ```pip``` and by doing the following in the terminal,

```python
# Create the virtual environment in PyThea's root folder
python3 -m venv env

# Activate the environment
source env/bin/activate

# install the required packages using pip3
pip3 install -r requirements.txt

# When you are done you can deactivate a virtual environment
deactivate
```

Now you can run any part of the _PyThea_ (see Run locally section).

You may also add your _PyThea_ directory to the environment variable ```PYTHONPATH```. This is useful if you need to run _PyThea_ tests or when you need to run some of the package modules out of streamlit.

In the terminal use the following and change the \<PyTheaRootDir\> with your path.

```
export PYTHONPATH="${PYTHONPATH}:<PyTheaRootDir>/PyThea"
```

For a permanent solution, if you're using bash (on a Mac or GNU/Linux distribution), add the above line to your ```~/.bashrc``` file (changing the \<PyTheaRootDir\> with your path first).

## 🐾 Run locally the _PyThea_ application
Install the required Python packages and activate the environment as shown above for each case.

Note that: For the ```pip``` method the ```source env/bin/activate``` is done in the root directory the ```conda activate PyThea``` can be done everywhere.

Then navigate into the _PyThea_ package directory (e.g. ```cd <PyTheaRootDir>/PyThea```) and run the application with streamlit.

```
streamlit run PyThea_app.py
```
The application should now open in the default browser!

When you finish the work deactivate the environment as shown above for each case.

## 📙 Usage

Complete documentation of the _PyThea_ can be found in (under construction)

## 📦 Useful Python packages

- [SunPy](https://sunpy.org/): The community-developed, free and open-source solar data analysis environment for Python.
- [AstroPy](https://www.astropy.org/): The Astropy Project is a community effort to develop a single core package for Astronomy in Python.

## 📜 Acknowledging or Citing PyThea [![https://doi.org/10.5281/zenodo.5713659](https://zenodo.org/badge/DOI/10.5281/zenodo.5713659.svg)](https://doi.org/10.5281/zenodo.5713659)

If you use _PyThea_ for scientific work or research presented in a publication, please cite this package by including in the methods or acknowledgement section the following: "This research has made use of PyThea v?.?.?, an open-source and free Python package to reconstruct the 3D structure of CMEs and shock waves (cite software).". For the _PyThea_ software citation use the following: "Athanasios Kouloumvakos. (2021). PyThea: A software package to reconstruct the 3D structure of CMEs and shock waves. Zenodo. [https://doi.org/10.5281/zenodo.5713659](https://doi.org/10.5281/zenodo.5713659)". _PyThea_ has a strong dependency to SunPy and AstroPy Python packages, consider citing these packages as well. To acknowledge _PyThea_ in posters or talks include the project [logo](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_logo.png) or [icon](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_icon.png).

## ⓘ The mythology of Thea:

In Greek mythology, Thea, also called Euryphaessa "wide-shining", is the Titaness of sight and the shining light of the clear blue sky. Her brother/consort is Hyperion, a Titan and god of the sun, and together they are the parents of Helios (the Sun), Selene (the Moon), and Eos (the Dawn).
