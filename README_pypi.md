# PyThea: A software package to reconstruct the 3D structure of CMEs and shock waves

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Version](https://img.shields.io/github/v/release/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases)
[![Release Date](https://img.shields.io/github/release-date/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![flake8](https://github.com/AthKouloumvakos/PyThea/actions/workflows/flake8.yml/badge.svg)
![pytest](https://github.com/AthKouloumvakos/PyThea/actions/workflows/pytest.yml/badge.svg)

_PyThea_ is an open-source software package that can be used to reconstruct the 3D structure of Coronal Mass Ejections (CMEs) and shock waves and determine their kinematics using remote-sensing observations. The tool implements the Graduated Cylindrical Shell (GCS) model that can be used to reconstruct CMEs and two geometrical models, namely a spheroid and ellipsoid model to reconstruct shock waves. It also implements remote-sensing observations from multiple viewpoints such as the Solar and Heliospheric Observatory (SoHO) and Solar Terrestrial Relations Observatory (STEREO).

## 💾 Installation

We recommend we recommend to create a virtual environment for _PyThea_ and installing this using ```pip```.

Use the following commands to create a python virtual environment, install _PyThea_ with ```pip```, and then run the application locally.

```python
# Create a virtual environment
python3 -m venv PyThea

# Activate the environment
source PyThea/bin/activate

# install the required packages using pip3
pip3 install PyThea

# Run locally _PyThea_ application running PyThea streamlit (see next step)
PyThea streamlit

# When you are done you can deactivate the virtual environment
deactivate
```

If you use conda (or anaconda) you can also install _PyThea_ in a similar way, as follows:

```python
# Create a virtual environment
conda create --name PyThea python=3.9

# Activate the environment
conda activate PyThea

# install the required packages using pip3
pip3 install PyThea

# Run locally _PyThea_ application running PyThea streamlit (see next step)
PyThea streamlit

# When you are done you can deactivate the virtual environment
conda deactivate
```

## 🐾 Run locally the _PyThea_ application

After installing _PyThea_ software package you can run the application using the terminal.

If the environment is not active then use ```source PyThea/bin/activate``` to activate this and then run _PyThea_ with

```
PyThea streamlit
```

The application should now open in the default browser!

Deactivate the environment when you finish.

## 📙 Usage

Complete documentation of the _PyThea_ is under construction.

## 📦 Useful Python packages

- [SunPy](https://sunpy.org/): The community-developed, free and open-source solar data analysis environment for Python.
- [AstroPy](https://www.astropy.org/): The Astropy Project is a community effort to develop a single core package for Astronomy in Python.

## 📜 Acknowledging or Citing PyThea [![https://doi.org/10.5281/zenodo.5713659](https://zenodo.org/badge/DOI/10.5281/zenodo.5713659.svg)](https://doi.org/10.5281/zenodo.5713659)

If you use _PyThea_ for scientific work or research presented in a publication, please cite this package by including in the methods or acknowledgement section the following: "This research has made use of PyThea v?.?.?, an open-source and free Python package to reconstruct the 3D structure of CMEs and shock waves (cite software).". For the _PyThea_ software citation use the following: "Athanasios Kouloumvakos. (2021). PyThea: A software package to reconstruct the 3D structure of CMEs and shock waves. Zenodo. [https://doi.org/10.5281/zenodo.5713659](https://doi.org/10.5281/zenodo.5713659)". _PyThea_ has a strong dependency to SunPy and AstroPy Python packages, consider citing these packages as well. To acknowledge _PyThea_ in posters or talks include the project [logo](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_logo.png) or [icon](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_icon.png).

## ⓘ The mythology of Thea:

In Greek mythology, Thea, also called Euryphaessa "wide-shining", is the Titaness of sight and the shining light of the clear blue sky. Her brother/consort is Hyperion, a Titan and god of the sun, and together they are the parents of Helios (the Sun), Selene (the Moon), and Eos (the Dawn).
