import os
import numpy
from setuptools import setup, find_packages

CURR_PATH = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(CURR_PATH, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
setup(
    name="LibPrep",
    version="1.0",
    description='Python package to analyse and prepare libraries of chemical compounds for molecular simulations.',
    long_description=long_description,
    url="https://github.com/carlesperez94/lib_prep",
    author='Carles Perez Lopez',
    author_email='carlesperez94@gmail.com',
    license='',
    packages=find_packages(exclude=['docs', 'tests']),
    include_dirs=[numpy.get_include()],
    install_requires=['rdkit', 'ProDy==1.10'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research"
    ],
    project_urls={
    'Documentation': 'https://carlesperez94.github.io/lib_prep/',
    'Source': 'https://carlesperez94.github.io/lib_prep/',
'Tracker': 'https://github.com/carlesperez94/lib_prep/issues',
},
)