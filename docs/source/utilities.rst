.. _utilities:

Utilities
=========

PyThea includes several utilities that facilitate various tasks without requiring inputs from the Streamlit components. Therefore these functions can be used outside of the main application allowing the users to leverage some of PyThea's functionality without relying on the main application interface. Most of the utilities focus on data processing and visualization. Bellow there a few examples that show how to use the utility fucntions to download, process and plot imaging data from various instruments, load and process JSON fitting files, and analyze kinematic data. This last feature is particularly useful for users who need to perform detailed kinematic analysis without engaging the main application and produce their final figures for publication. Many of these utility functions are also demonstrated through examples available in the example section of the documentation. These examples provide step-by-step guidance on how to effectively use the utilities for various tasks.

Imaging data
------------

This section introduces the utilities related to the imaging download, process, and virtualization and provides details of how the imaging data are managed within PyThea.

PyThea utilizes remote sensing data from several spacecraft and various imagers on board heliospheric missions that continuously collect solar imaging data.The available imagers are specified in the ``selected_imagers.py`` configuration file. Each imager has several configuration options used in the data downloading and processing pipeline. These configurations are stored in the ``selected_imagers.imager_dict`` which is a dictionary.

To retrieve the list of available imagers and their corresponding keys, use the following script:

.. code-block:: python

    >>> from PyThea.config import selected_imagers

    >>> for key in selected_imagers.imager_dict.keys():
    >>>     imager = selected_imagers.imager_dict[key]
    >>>     detector_or_wavelenght = imager['detector'] if 'detector' in imager else imager['wavelength']
    >>>     print(f'{imager["source"]}/{imager["instrument"]}-{detector_or_wavelenght}:  {key}')

The ``imager_dict`` is a dictionary where each primary key is a short label for a selected imager. These short labels are used in the fitting process as identifiers for the selected imagers and the information stored in the JSON fitting file. For each imager, ``imager_dict[key]`` is another dictionary containing definitions for the download and processing pipeline. For example, ``imager_dict['LC2']['fido']`` is a list of parameters used by SunPy's Fido function to download data from LASCO instrument and C2 detector and ``imager_dict['LC2']['process']`` provides the imaging processing configuration details of which are shown bellow.

For the first iteration (where key='LC2'), the script will print:

.. code-block:: python

    SOHO/LASCO-C2:  LC2
