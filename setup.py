import os

from setuptools import find_packages, setup
from setuptools.config import read_configuration

extras = read_configuration('setup.cfg')

# create home directory
if not os.path.isdir(os.path.join(os.environ['HOME'], '.PyThea')):
    os.mkdir(os.path.join(os.environ['HOME'], '.PyThea'))

with open('README.md', 'r') as f:
    long_description=f.read()

# get requirements
with open('requirements.txt') as f:
    requirements=f.read().splitlines()

setup(
    name='PyThea',
    use_scm_version={'write_to': os.path.join('PyThea', '_version.py')},
    # version="0.4.0",
    description='PyThea: A software package to reconstruct the 3D structure of CMEs and shock waves',
    url='https://github.com/AthKouloumvakos/PyThea',
    # long_description=long_description,
    long_description_content_type='text/markdown',

    author='Athanasios Kouloumvakos',
    author_email='athkouloumvakos@gmail.com',
    license='GPL-3.0',
    #license_file="LICENSE.md",

    python_requires='>=3.8',
    install_requires=requirements,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'pythea = PyThea.pythea_cli:main'
        ]
    },
    extras_require=extras,
)
