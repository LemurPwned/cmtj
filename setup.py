from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

setup(name='pymtj',
      version='1.0',
      description='A Python3 package for Magnetic Tunnel Junctions in spintronics and magnetics',
      maintainer='LemurPwned',
      install_requires=[
            'Cython>=0.29.13',
            'matplotlib>=3.1.1',
            'numpy>=1.16.4',
            'pandas>=0.24.2'
      ],
      ext_modules=cythonize("pymtj/junction.pyx",
                            language_level="3",
                            compiler_directives={'language_level': "3"}),
      include_dirs=[
          numpy.get_include(),
      ],
      classifiers=["Programming Language :: Python :: 3"],
      packages=find_packages())
