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


Download and Load fits files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To download imaging data for a selected imager you can use the ``download_fits`` utility. The ``download_fits`` utility uses SunPy's ``Fido`` which is a unified data search and retrieval interface (see details in `Fido <https://docs.sunpy.org/en/stable/tutorial/acquiring_data/>`_). The function requires two inputs from the user, first is the time range of the imaging query and second is the key of the selected imager. The time range is a search attribute from ``sunpy.net.attrs`` that defines the start and end time of the query, for example ``attrs.Time('2021-10-28T15:30:00', '2021-10-28T15:30:00')``. The code below will download LASCO-C2 fits files for the selected time range:

.. code-block:: python

    >>> from sunpy.net import attrs as a
    >>> from PyThea.utils import download_fits

    >>> imager = 'LC2'
    >>> timerange = a.Time('2021-10-28T15:30:00', '2021-10-28T17:30:00')
    >>> files = download_fits(timerange, imager)

In this example, ``Fido`` searches and downloads imaging data from LASCO-C2 from the Virtual Solar Observatory using the attributes that are defined in the ``imager_dict['LC2']['fido']``. For LASCO-C2 the ``imager_dict['LC2']['fido']`` list of atributes define only the instrument (``a.Instrument.lasco``) and the detector (``a.Detector.c2``). For other imagers the list of attributes maybe more extended. For example, the ``imager_dict['LC2']['fido']`` for the SDO/AIA imager contains also the attributes for wavelength (``a.Wavelength(19.3 * u.nm)``) and sample rate (``a.Sample(1*u.minute)``) to be used in the data query.

The fetched fits files from this downloading process are stored in PyThea's database directory which is pre-configured to be the ``os.path.join(Path.home(), 'PyThea')`` directory. Then the imaging data are stored in the ``\data`` directory with the following sub-directory structure for each imager: ``{imager_prop["source"]}/{imager_prop["instrument"]}/{sub_path}``, where the ``sub_path`` is either ``imager_prop['detector']`` or ``imager_prop['wavelength']`` depending on the imager.

To load the fits files in a sequence of SunPy maps use the ``load_fits`` utility as shown in the code bellow:

.. code-block:: python

    >>> from PyThea.utils import load_fits

    >>> maps = load_fits(flies)

The ``load_fits`` utility checks if all the fits files can be loaded in a SunPy ``Map`` and if errors are detected then it removes the file. These errors may occur because, on rare occasions, some fits files are downloaded corrupted and the loading process fails. The returned ``maps`` are ``sunpy.map.MapSequence`` which is a series of Maps in a single object.

In the 'Gallery of Examples', you can find an example of the above downloading and loading process.

Process fits files
~~~~~~~~~~~~~~~~~~

The next step is to process the loaded maps using the ``single_imager_maps_process`` utility. This function filters, prepares, and processes the maps. The ``single_imager_maps_process`` takes as input the laded maps and options for the maps processing and returs the processed maps as ``sunpy.map.MapSequence``. The default options for the maps processing can be found in ``imager_dict[imager]['process']``, for each imager.

.. code-block:: python

    >>> from PyThea.config import selected_imagers

    >>> imager = 'LC2'
    >>> print(selected_imagers.imager_dict[imager]['process'])

.. code-block:: python
    {'dimensions': (<Quantity 1024. pix>, <Quantity 1024. pix>), 'polar': 'Clear', 'superpixel': 2}

According to the code above, the default configuration for the image processing of LASCO-C2 maps filters out images with dimensions other than 1024x1024, images that are not total brightness images and resamples the images to half the original dimension. The user can provide a different configuration than the default one, however, it is not advisable to change the configuration for any imager when running the application.

The image processing in ``single_imager_maps_process`` utility consists of three different processing steps:

* Filters the Maps (``filter_maps``)

The map filtering uses the ``filter_maps`` function which is a part of the ``maputils`` utilites. With this function, the Maps can be filtered based on their exposure time, data dimension, and polarization (for white-light images).

* Prepare the Maps (``prepare_maps``)

With the ``prepare_maps`` function, which is also a part of the ``maputils`` utilities, the Maps are prepared from their standard level to a higher level. These preparations include various calibrations and corrections, such as pointing, observer location, exposure time corrections and the preparation of polarization images into total brightness images. Then the Maps are normalized to their exposure time and if the images are coronagraphic the occulter is masked. Finally, the Maps are resampled with the SunPy's ``superpixel`` method.

* Process the sequence (``maps_sequence_processing``)

The final prepared Maps can be processed into running/base difference or plain image sequence maps using the ``maps_sequence_processing`` function which is part of the ``maputils``utilities as well.

The code below gives a small example of how the ``maps_sequence_processing`` function is being used.

.. code-block:: python

    >>> from PyThea.utils import single_imager_maps_process

    >>> processed_images = single_imager_maps_process(maps,
                                                      **selected_imagers.imager_dict[imager]['process'],
                                                      image_mode='Running Diff.',
                                                      diff_num=1)
