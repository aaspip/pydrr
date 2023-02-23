#!/usr/bin/env python
# -*- encoding: utf8 -*-
import glob
import inspect
import io
import os

from setuptools import find_packages
from setuptools import setup


long_description = """
Source code: https://github.com/aaspip/pydrr""".strip() 

def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")).read()


from distutils.core import Extension
# defining calculator_module as an extension class instance.
# 'calculator' is extension's name and sources is a file name list.
calculator_module = Extension('calculator', sources=['pydrr/src/calculator.c'])

setup(
    name="pydrr",
    version="0.0.2.1",
    license='GNU General Public License, Version 3 (GPLv3)',
    description="A python package for the damped rank reduction (DRR) method and its variants. The DRR method has a variety of applications in seismology, including but not limited to seismic denoising, seismic reconstruction, seismic diffraction separation, constrained LSRTM, constrained FWI, etc.",
    long_description=long_description,
    author="pydrr developing team",
    author_email="chenyk2016@gmail.com",
    url="https://github.com/aaspip/pydrr",
    ext_modules=[calculator_module],
    packages=['pydrr'],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
    ],
    keywords=[
        "seismology", "exploration seismology", "array seismology", "denoising", "science", "signal-to-noise ratio", "damped rank reduction method"
    ],
    install_requires=[
        "numpy", "scipy", "matplotlib"
    ],
    extras_require={
        "docs": ["sphinx", "ipython", "runipy"]
    }
)
