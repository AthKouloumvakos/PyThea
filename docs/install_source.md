## 💾 Installation from source files

_PyThea_ is written in Python >=3.8 and has some package requirements to run. These are listed in the requirements.txt and environment.yml files located in the main repository. To run this application locally, we recommend to create its own virtual environment in Python and install the dependent packages. Because of a range of dependencies that the different packages have, the simplest way to work with _PyThea_ is installing the packages with _conda_ and creating its own environment. Alternatively, you can create a virtual environment in Python inside the project folder and install the packages using ```pip```. Both ways are explained below.

First, download the latest release of _PyThea_ from here: https://github.com/AthKouloumvakos/PyThea/releases/latest

Unzip the file and navigate (e.g., ```cd```) to the root directory of _PyThea_ which is in the unzipped folder. Follow the instructions to create the virtual environment and install the dependent packages:

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
