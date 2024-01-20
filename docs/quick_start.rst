
Getting Started
===============

Requirements
------------

Python 3.8 or later is required to install and run PyThea.

If you do not have Python installed already, use these `instructions <https://www.python.org/downloads>`_ to download and install it.

We recommend to install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or `Anaconda <https://www.anaconda.com/download>`_, that come with a suite of packages useful for scientific data analysis and ease the installation process. These two provide `Conda <https://docs.conda.io/en/latest/>`_ which is a an open-source package and environment management system that runs on Windows, macOS, and Linux.

Conda as a package manager helps you find and install packages. With ``Conda`` you can quickly install, run, and update packages and their dependencies and you can also easily create, save, load, and switch between environments on your computer. Following the instuctions bellow, you can use ``Conda`` to make a separate environment to run ``PyThea``, while you can continue running your usual version of Python in your normal environment.

.. _installing-pythea:

Installing PyThea
-----------------

.. tip::

  We recommend to create a virtual environment before installing ``PyThea``.

If you use `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or `Anaconda <https://www.anaconda.com/download>`_
you can install ``PyThea`` creating a virtual environment first with `Conda <https://docs.conda.io/en/latest/>`_
and then installing the package using ``pip``.

**Creating, Activating, and Deactivating a Conda Environment:**

To create a new Conda environment called "PyThea", and then activate it follow these steps:

1. Open a terminal or command prompt.

2. Use the following command to create a new Conda environment named "PyThea" and install Python:

   .. code-block:: bash

      conda create --name PyThea python=3.8

   You can replace "3.8" with the desired Python version.

3. Use the following command to activate the "PyThea" environment:

   .. code-block:: bash

      conda activate PyThea

   After executing this command, your prompt should change to indicate that the "PyThea" environment is active.

You can now install ``PyThea`` and work within this isolated environment when using this software package.

To deactivate the "PyThea" Conda environment and return to the base environment use,

   .. code-block:: bash

      conda deactivate

 The environment will be deactivated, and you will return to the base Conda environment.

**Installing PyThea:**

.. code-block:: bash

  # Create the virtual environment
  conda create --name PyThea python=3.9

  # Activate the environment
  conda activate PyThea

  # Install the required packages using pip
  pip install PyThea

.. warning::

  Currently install with ``conda install PyThea`` is not suppoted.

.. _runnig-pythea:

Run the Application
-------------------

After installing ``PyThea`` software package, you can run the application using the terminal.

Activate the enviroment and then run ``PyThea`` using the following,

.. code-block:: bash

  PyThea streamlit

The application will open in the default browser.

Deactivate the environment when you finish your work.

.. code-block:: bash

  # When you are done you can deactivate the virtual environment
  conda deactivate

.. _acknowledging-pythea:

Acknowledging or Citing
-----------------------

If you use PyThea for scientific work or research presented in a publication, please cite it by acknowledging in the main text its use and include the following citation:

.. note:: "Athanasios Kouloumvakos et al. (2022). PyThea: An open-source software package to perform 3D reconstruction of coronal mass ejections and shock waves, Front. Astron. Space Sci. 9:974137. (DOI: `10.3389/fspas.2022.974137 <https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/>`)".

This article published in Frontiers in Astronomy and Space Sciences as part of the research topic "Snakes on a Spaceship: An Overview of Python in Space Physics" and can be found here: https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/

Also include in the methods or acknowledgement section the following:

.. note:: "This research has made use of PyThea v?.?.?, an open-source and free Python package to reconstruct the 3D structure of CMEs and shock waves (Zenodo: https://doi.org/10.5281/zenodo.5713659).". |zenodo-badge|

.. |zenodo-badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5713659.svg
   :target: https://doi.org/10.5281/zenodo.5713659
   :alt: zenodo

changing the v?.?.? to the version you have used.

To acknowledge PyThea in posters or talks include the project logo or icon.

PyThea has a strong dependency to SunPy and AstroPy Python packages, consider citing these packages as well.
