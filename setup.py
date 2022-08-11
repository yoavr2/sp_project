from setuptools import setup, find_packages, Extension
import os

setup(name='mykmeanssp',
      version='1.0',
      description='kmeans algorithm',
      ext_modules=[Extension('mykmeanssp', sources=['spkmeansmodule.c'])])
      

