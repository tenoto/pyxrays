#!/usr/bin/env python

from setuptools import setup, find_packages
from pyxrays import __author__, __version__, __license__

setp(
	name		 = 'pyxrays',
	version 	 = __version__,
	description  = 'python scripts for general X-ray analyses'
	license 	 = __license__,
	author  	 = __author__,
	author_email = 'teruaki.enoto@gmail.com',
	url     	 = 'https://github.com/tenoto/pyxrays.git',
	keywords 	 = 'xray analysis',
	packages 	 = find_packages(),
	install_requires = []
	)