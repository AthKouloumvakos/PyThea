..
   PyThea documentation master file, created by
   sphinx-quickstart on Mon May  9 22:22:34 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyThea's Documentation
======================

.. image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
   :target: https://www.python.org/
   :alt: python

.. image:: https://img.shields.io/pypi/v/PyThea?style=flat&logo=pypi
   :target: https://pypi.org/project/PyThea/
   :alt: pypi

.. image:: https://img.shields.io/github/release-date/AthKouloumvakos/PyThea
   :target: https://github.com/AthKouloumvakos/PyThea/releases
   :alt: Release Date

.. image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0
   :alt: License: GPL v3

.. image:: https://github.com/AthKouloumvakos/PyThea/actions/workflows/flake8.yml/badge.svg

.. image:: https://github.com/AthKouloumvakos/PyThea/actions/workflows/pytest.yml/badge.svg

PyThea is an open-source software package that can be used to reconstruct the 3D structure of Coronal Mass Ejections (CMEs) and shock waves and determine their kinematics using remote-sensing observations. The tool implements the Graduated Cylindrical Shell (GCS) model that can be used to reconstruct CMEs and two geometrical models, namely a spheroid and ellipsoid model to reconstruct shock waves. It also implements remote-sensing observations from multiple viewpoints such as the Solar and Heliospheric Observatory (SoHO) and Solar Terrestrial Relations Observatory (STEREO).

Content
========

.. toctree::
   :maxdepth: 1

   quick_start
   users_guide

Acknowledging or Citing PyThea
==============================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5713659.svg
   :target: https://doi.org/10.5281/zenodo.5713659
   :alt: zenodo

If you use PyThea for scientific work or research presented in a publication, please cite it by acknowledging in the main text its use and include the following citation:

.. note:: "Athanasios Kouloumvakos et al. (2022). PyThea: An open-source software package to perform 3D reconstruction of coronal mass ejections and shock waves, Front. Astron. Space Sci. 9:974137. (DOI: 10.3389/fspas.2022.974137)".

This article published in Frontiers in Astronomy and Space Sciences as part of the research topic "Snakes on a Spaceship: An Overview of Python in Space Physics" and can be found here: https://www.frontiersin.org/articles/10.3389/fspas.2022.974137/

Also include in the methods or acknowledgement section the following:

.. note:: "This research has made use of PyThea v?.?.?, an open-source and free Python package to reconstruct the 3D structure of CMEs and shock waves (Zenodo: https://doi.org/10.5281/zenodo.5713659).".

changing the v?.?.? to the version you have used.

To acknowledge PyThea in posters or talks include the project logo or icon.

PyThea has a strong dependency to SunPy and AstroPy Python packages, consider citing these packages as well.

The mythology of Thea
=====================

In Greek mythology, Thea, also called Euryphaessa "wide-shining", is the Titaness of sight and the shining light of the clear blue sky. Her brother/consort is Hyperion, a Titan and god of the sun, and together they are the parents of Helios (the Sun), Selene (the Moon), and Eos (the Dawn).
