import os, inspect
from distutils.core import setup, Extension
file_path = os.path.dirname(inspect.stack()[0][1])

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='fastpace',
    version='0.1.0',
    author="Hazem M. Kotb",
    author_email="hazem.mamdouh.kotb@gmail.com",
    description="FaSTPACE: A Fast and Scalable Tool for Peptide Alignment and Consensus Extraction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hkotb/fastpace",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    ext_modules=[
        Extension('fastpace', [
            os.path.join(file_path, 'src', 'fastpace_module_main.c'), 
            os.path.join(file_path, 'src', 'fastpace_utils.c'),
            os.path.join(file_path, 'src', 'scipy_interface.cpp')])
    ])
