.. LibPrep documentation master file, created by
   sphinx-quickstart on Tue Nov 12 16:22:00 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
==========================

.. toctree::
   :maxdepth: 2

Conda (Recommended)
-------------------
If you want, you can create a conda environment to work with lib_prep::

	conda create -n lib_prep_env python=3.7
	conda activate lib_prep_env

Otherwise, you can just continue here by installing the package into your current conda environment::

        conda install -c carlesperez94 -c rdkit -c nostrumbiodiscovery lib_prep


Now, it should work. We strongly recommend you to run the tests, but is not mandatory. 

Just remember that if you want to logout of this conda environment run the following when you finish::

	conda deactivate


Source
------
It is not recommended to install lib_prep from source, because you will need to install rdkit on your own. But if you feel brave, just ensure
that you have rdkit installed in your sistem and run the setup.py in the root lib_prep directory::

	git clone https://github.com/carlesperez94/lib_prep.git
	cd lib_prep
	python setup.py install

If you have permissions problems try with::

	sudo python setup.py install

Test
----
First of all you need to install pytest into your conda envirnoment::

	conda install pytest

If you have pytest already installed, go to the testing directory and run test_all.py::

	cd lib_prep/tests
	pytest test_all.py



