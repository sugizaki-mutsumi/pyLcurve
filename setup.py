#!/usr/bin/env python

from setuptools import setup, Extension  #, find_packages


lxdata_module = Extension(
    '_lxdata',
    sources = ['interpolate.cpp', 'lxdata.cpp', 'lxdata.i'],
    swig_opts=['-c++'],
    libraries=["cfitsio"], 
)


setup(
    name = 'lxdata',
    version = '0.1',
    author      = "Mutsumi Sugizaki",
    description = "Module to calculate light curve integration",
    ext_modules = [lxdata_module],
    py_modules = ["lxdata"],
)
