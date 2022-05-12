import os

from setuptools import find_packages, setup

# from setuptools.config import read_configuration

# extras = read_configuration('pyproject.toml')

# create home directory
if not os.path.isdir(os.path.join(os.environ['HOME'], '.PyThea')):
    os.mkdir(os.path.join(os.environ['HOME'], '.PyThea'))

with open('README_pypi.md', 'r') as f:
    long_description = f.read()

# get requirements
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='PyThea',
    use_scm_version={'write_to': os.path.join('PyThea', '_version.py')},

    description='PyThea: A software package to reconstruct the 3D structure of CMEs and shock waves',
    url='https://github.com/AthKouloumvakos/PyThea',
    long_description=long_description,
    long_description_content_type='text/markdown',

    author='Athanasios Kouloumvakos',
    author_email='athkouloumvakos@gmail.com',
    license='GPL-3.0',
    license_file='LICENSE.md',
    keywords=['science', 'solar physics', 'solar', 'sun', 'shock waves'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
    ],

    python_requires='>=3.8, <3.10',
    install_requires=requirements,

    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'pythea = PyThea.pythea_cli:main'
        ]
    },
    # extras_require=extras,
)
