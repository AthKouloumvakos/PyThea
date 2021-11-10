import os
import sys

# https://stackoverflow.com/questions/20971619/ensuring-py-test-includes-the-application-directory-in-sys-path
# Make sure that the application source directory (this directory's parent) is on sys.path.
# or else before runing pytest 'export PYTHONPATH='${PYTHONPATH}:<PyTheaRootDir>/PyThea''

dir_pythea = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'PyThea'))
sys.path.insert(0, dir_pythea)
