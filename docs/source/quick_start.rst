Getting Started
===============

Requirements
------------

Python 3.8 or later is required to install and run PyThea.

We recommend using `Anaconda <https://www.anaconda.com>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, that come with a suite of packages useful for scientific data analysis and ease the installation process.

Installation
------------

.. tip::

  We recommend to create a virtual environment before installing ``PyThea``.

If you use `Anaconda <https://www.anaconda.com>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
you can install ``PyThea`` creating a virtual environment first with `Conda <https://docs.conda.io/projects/conda/en/latest/>`_
and then installing the package using ``pip``. In the terminal do the following:

In the terminal this is done as follows:

.. code-block:: bash

  # Create the virtual environment
  conda create --name PyThea python=3.9

  # Activate the environment
  conda activate PyThea

  # install the required packages using pip
  pip install PyThea

.. warning::

  Currently install with ``conda install PyThea`` is not suppoted.

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
