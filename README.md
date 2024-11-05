# PyThea: A software package to reconstruct the 3D structure of CMEs and shock waves

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Version](https://img.shields.io/github/v/release/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases)
[![Release Date](https://img.shields.io/github/release-date/AthKouloumvakos/PyThea)](https://github.com/AthKouloumvakos/PyThea/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![flake8](https://github.com/AthKouloumvakos/PyThea/actions/workflows/flake8.yml/badge.svg)
![pytest](https://github.com/AthKouloumvakos/PyThea/actions/workflows/pytest.yml/badge.svg)
[![pypi](https://img.shields.io/pypi/v/PyThea?style=flat&logo=pypi)](https://pypi.org/project/PyThea/)

![Logo](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_logo.png)

_PyThea_ is an open-source software package that can be used to reconstruct the 3D structure of Coronal Mass Ejections (CMEs) and shock waves and determine their kinematics using remote-sensing observations. The tool implements the Graduated Cylindrical Shell (GCS) model that can be used to reconstruct CMEs and two geometrical models, namely a spheroid and ellipsoid model to reconstruct shock waves. It also implements remote-sensing observations from multiple viewpoints such as the Solar and Heliospheric Observatory (SoHO), Solar Terrestrial Relations Observatory (STEREO), and Parker Solar Probe.

## 💾 Installation

We recommend, creating a virtual environment for _PyThea_ and installing the package from PyPI using ```pip```.

If you use Anaconda or Miniconda (which we also recommend) you can create a virtual environment using ```conda``` and then install _PyThea_ using from PyPI using ```pip```. In the terminal do the following:

```python
# Create a virtual environment. Use python>3.9
conda create --name PyThea python=3.10

# Activate the environment
conda activate PyThea

# install the required packages using pip3
pip3 install PyThea

# Run locally the application with streamlit (see also the section "Run locally the PyThea application" )
PyThea streamlit

# When you are done you can deactivate the virtual environment
conda deactivate
```

If ```conda``` is not your favorite way of creating a virtual environment in python, then you can manually create it and install _PyThea_ with ```pip``` as previously shown. For example, you can do the following:

```python
# Create a virtual environment.
python3 -m venv PyThea

# Activate the environment
source PyThea/bin/activate

# install the required packages using pip3
pip3 install PyThea

# Run locally the application with streamlit (see also the section "Run locally the PyThea application" )
PyThea streamlit

# When you are done you can deactivate the virtual environment
deactivate
```

At the directory where the terminal is open, this method will create a folder named ```/PyThea``` and install _PyThea_ and all the required packages inside.

### ⬆️ Update

To update the package to the latest version activate the environment first and then:

```python
# Update _PyThea_ using pip3
pip3 install PyThea -U
```

One way to see witch version is installed in your environment is to open a python session and do:
```python
import PyThea
PyThea.__version__
```
You can also see the verion used and the latest version on the main page of the application.

## 🐾 Run locally the _PyThea_ application

After installing _PyThea_ software package you can run the application using the terminal.

If the environment is not active then use ```conda activate PyThea``` or ```source PyThea/bin/activate``` to activate this and then run _PyThea_ with

```
PyThea streamlit
```

The application should now open in the default browser!

Deactivate the environment when you finish.

If there is an error when running ```PyThea streamlit``` then you can manually run PyThea following these steps:

1) Locate where PyThea_app.py have been instaled. This is usually inside the anaconda3 enviroment folder. For example here ```~/opt/anaconda3/envs/PyThea/lib/python3.9/site-packages/PyThea```

2) Then run the PyThea_app.py using: ```streamlit run <PyTheaRootDir>/PyThea_app.py``` replacing  \<PyTheaRootDir\> with your path first.

## 📙 Usage

Complete documentation of the _PyThea_ can be found in [https://www.pythea.org/](https://www.pythea.org/).

## 📦 Useful Python packages

- [SunPy](https://sunpy.org/): The community-developed, free and open-source solar data analysis environment for Python.
- [AstroPy](https://www.astropy.org/): The Astropy Project is a community effort to develop a single core package for Astronomy in Python.
- [gcs_python](https://github.com/johan12345/gcs_python/): An implementation of the Graduated Cylindrical Shell model in python.

## 📜 Acknowledging or Citing PyThea [![https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/](https://img.shields.io/static/v1?label=Paper&message=Frontiers&color=red)](https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/) [![https://doi.org/10.5281/zenodo.5713659](https://zenodo.org/badge/DOI/10.5281/zenodo.5713659.svg)](https://doi.org/10.5281/zenodo.5713659)

If you use _PyThea_ for scientific work or research presented in a publication, please mention it in the main text and cite _PyThea_ paper (see [Paper](https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/) or ADS). Additionally, add in the methods or acknowledgements section the following: "This research has made use of PyThea v?.?.?, an open-source and free Python package to reconstruct the 3D structure of CMEs and shock waves (Zenodo: [https://doi.org/10.5281/zenodo.5713659](https://doi.org/10.5281/zenodo.5713659))." and changing the v?.?.? to the version you have used. To acknowledge _PyThea_ in posters or talks include the project [logo](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_logo.png) or [icon](https://github.com/AthKouloumvakos/PyThea/blob/master/docs/logo/pythea_icon.png). _PyThea_ has a strong dependency on SunPy and AstroPy Python packages, consider citing these packages as well.

## ⓘ The mythology of Thea:

In Greek mythology, Thea, also called Euryphaessa "wide-shining", is the Titaness of sight and the shining light of the clear blue sky. Her brother/consort is Hyperion, a Titan and god of the sun, and together they are the parents of Helios (the Sun), Selene (the Moon), and Eos (the Dawn).

## Development Support:

The lead author of this software package Athanasios Kouloumvakos acknowledges financial support from NASA Grant 80NSSC24K0071   for the further development and improvement of PyThea during 2024. This grant was part of the NASA Headquarters Heliophysics Tools and Methods Program in response to NASA ROSES–2022 (NNH22ZDA001N).
