from setuptools import setup, find_packages

VERSION = "0.1.0"
DESCRIPTION = "Find Representative: Find a genome that is representing all the genomes the best"
setup(
	name="find_representative",
	version=VERSION,
	author="Jordi Weijers",
	author_email="s2228483@vuw.leidenuniv.nl",
	description=DESCRIPTION,
	packages=find_packages(),
	entry_points={'console_scripts': ["getphylo=bin.main:main"]}
