"""
    PyThea: A software package to reconstruct the 3D structure of CMEs and
    shock waves using multi-viewpoint remote-sensing observations.
    Copyright (C) 2021  Athanasios Kouloumvakos

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import collections.abc
import json
import os
import re

from PyThea.config import database_dir_default
from PyThea.config.selected_imagers import imager_dict
from PyThea.utils import download_fits


def create_nested_dict(k, v):
    if len(k) == 0:
        return v
    else:
        return {k[0]: create_nested_dict(k[1:], v)}


def load_nested_dict(k, d):
    if len(k) == 0:
        return d
    else:
        if k[0] in d:
            return load_nested_dict(k[1:], d[k[0]])
        else:
            return None


def update_nested_dict(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_nested_dict(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def get_fits_filenames_from_database(event_id, timerange, imager):

    db_args = [imager_dict[imager]['source'],
               imager_dict[imager]['instrument'],
               imager_dict[imager].get('detector', imager_dict[imager].get('wavelength'))]

    # name of the JSON file with the database
    file_id = re.search(r'D(.*?)(?=T)', event_id).group(1)
    filename = os.path.join(database_dir_default, 'data_manager', f'{file_id}_fits_filenames.json')

    if not os.path.isfile(filename):
        # create file if it doesn't exist
        with open(filename, 'w') as f:
            json.dump({}, f)

    # load data from file
    with open(filename, 'r') as f:
        database_data = json.load(f)

    downloaded_files = None

    if str(timerange) in database_data:
        downloaded_files = load_nested_dict(db_args, database_data[str(timerange)])

    if downloaded_files is None:
        downloaded_files = download_fits(timerange, imager)
        if not downloaded_files:
            return None
        nested_dict_ = create_nested_dict(db_args, downloaded_files.data)
        nested_dict_ = {str(timerange): nested_dict_}

        # dictionary stored in JSON file
        database_data = update_nested_dict(database_data, nested_dict_)

        # write updated data to file
        with open(filename, 'w') as f:
            json.dump(database_data, f, indent=4)

    return downloaded_files
