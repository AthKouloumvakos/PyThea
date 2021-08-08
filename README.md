# PyThea a CME and shock wave forward modeling software package

_PyThea_ is an open-source software package that can be used to perform coronal mass ejection (CME) and shock wave forward modeling using remote-sensing observations. The tool implements various geometrical models that can be used to reconstruct CME and shock waves. It also implements remote-sensing observations from multiple viewpoints such as the Solar and Heliospheric Observatory (SoHO) and Solar Terrestrial Relations Observatory (STEREO). An online preview of this tool is available at [https://athkouloumvakos.github.io/PyThea](https://athkouloumvakos.github.io/PyThea).

## Installation

_PyThea_ is written in Python >=3.8 and has some package requirements, which are listed in the requirements.txt and environment.yml files. 
To run localy this application we recomend to create its own virtual enviroment in python.

**Recomended (conda)**

Because of a range of dependencies that packages have, the simplest way to work with _PyThea_ 
is in conda and to create its own environment using the ```conda env create```. If you already have Anaconda installed, then with the anaconda prompt active you can, ```cd``` the root directory of _PyThea_ and in your terminal do:

```python
# Create a virtual environment in python using conda
#  see https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
conda env create -f environment.yml
conda info --envs

# Activate the enviroment
conda activate pythea

# When you are done you can deactivate a virtual environment
conda deactivate
```

**Alternative (pip)**

```python
# You can create a virtual environment in python inside the project folder.
#  see https://docs.python.org/3/library/venv.html
python3 -m venv env

# Activate the enviroment
source env/bin/activate

# install the required packages using pip3
pip3 install -r requirements.txt

# When you are done you can deactivate a virtual environment
deactivate
```

Now you can run any part of the _PyThea_ (see Usage section).

### How to run localy the _PyThea_ application
Install the required python packages, and then run the application with streamlit. 
```
# Activate the enviroment as shown above
# You can now run the application using streamlit
streamlit run HDPt_app.py
```
The application should now open in the default browser!

# Usage

Complete documentation of the _PyThea_ can be found in (under construction)

## Usefull python packages 📦:
        
- [SunPy](https://sunpy.org/): The community-developed, free and open-source solar data analysis environment for Python.
- [AstroPy](https://www.astropy.org/): The Astropy Project is a community effort to develop a single core package for Astronomy in Python.

## The mythology of Thea

In Greek mythology, Thea, also called Euryphaessa "wide-shining", is the Titaness of sight and the shining light of the clear blue sky. Her brother/consort is Hyperion, a Titan and god of the sun, and together they are the parents of Helios (the Sun), Selene (the Moon), and Eos (the Dawn).
