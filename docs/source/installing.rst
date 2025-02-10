.. _installing-pythea:

Installing PyThea
=================

Requirements
------------

Python 3.10 or later is required to install and run PyThea.

If you do not have Python installed already, use these `instructions <https://www.python.org/downloads>`_ to download and install it.

We recommend to install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or `Anaconda <https://www.anaconda.com/download>`_, that come with a suite of packages useful for scientific data analysis and ease the installation process. These two provide `Conda <https://docs.conda.io/en/latest/>`_ which is a an open-source package and environment management system that runs on Windows, macOS, and Linux.

Conda as a package manager helps you find and install packages. With ``Conda`` you can quickly install, run, and update packages and their dependencies and you can also easily create, save, load, and switch between environments on your computer. Following the instuctions bellow, you can use ``Conda`` to make a separate environment to run ``PyThea``, while you can continue running your usual version of Python in your normal environment.

Create Virtual Environment
--------------------------

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

      conda create --name PyThea python=3.10

   You can replace "3.10" with the desired Python version.

3. Use the following command to activate the "PyThea" environment:

   .. code-block:: bash

      conda activate PyThea

   After executing this command, your prompt should change to indicate that the "PyThea" environment is active.

You can now install ``PyThea`` and work within this isolated environment when using this software package.

To deactivate the "PyThea" Conda environment and return to the base environment use,

   .. code-block:: bash

      conda deactivate

 The environment will be deactivated, and you will return to the base Conda environment.

Install with pip
----------------

.. code-block:: bash

  # Create the virtual environment
  conda create --name PyThea python=3.10

  # Activate the environment
  conda activate PyThea

  # Install the required packages using pip
  pip install PyThea

.. warning::

  Currently install with ``conda install PyThea`` is not suppoted.
