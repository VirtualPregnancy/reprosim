""" Reprosim library - routines for a modelling the lung

The Reprosim library is an advanced modelling library for models of the reproductive system.
"""

classifiers = """\
Development Status :: 4 - Beta
Intended Audience :: Developers
Intended Audience :: Education
Intended Audience :: Science/Research
License :: OSI Approved :: Apache Software License
Programming Language :: Python
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3.5
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: MacOS :: MacOS X
Topic :: Scientific/Engineering :: Medical Science Apps.
Topic :: Software Development :: Libraries :: Python Modules
"""

from setuptools import setup

doclines = __doc__#.split("\n")

setup(
  name='reprosim',
  version='@Reprosim_VERSION@',
  author='Reproduction and Development Group, Auckland Bioengineering Institute.',
  author_email='alys.clark@auckland.ac.nz',
  packages=['reprosim'],
  package_data={'reprosim': [@SETUP_PY_PACKAGE_FILES_STR@]},
  platforms=['any'],
  url='http://www.abi.auckland.ac.nz/en/about/our-research/development-and-reproductive-health.html',
  license='http://www.apache.org/licenses/LICENSE-2.0',
  description='Reprosim library of routines for modelling the reproductive system.',
  classifiers = filter(None, classifiers.split("\n")),
  long_description=doclines,
)
