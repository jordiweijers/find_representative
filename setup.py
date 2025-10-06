from setuptools import setup, find_packages

VERSION = "0.1.0"
DESCRIPTION = "Find the most representative genome from a set of similar genomes, based on the presence/absence of genes."
setup(
	name="find_representative",
	version=VERSION,
	author="Jordi Weijers",
	author_email="jweijers@ifam.uni-kiel.de",
	description=DESCRIPTION,
	packages=find_packages(),
	entry_points={'console_scripts': ["find_representative=bin.main:main"]}
)
