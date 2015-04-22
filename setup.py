#!/usr/bin/env python

from setuptools import setup, find_packages
import os

# Allow to run setup.py from another directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

version = '1.0.0'

setup(
    name='rosemary',
    version=version,
    url='https://github.com/siegelzero/rosemary',
    license='MIT',
    author='Kenneth Brown',
    author_email='siegel.zero@gmail.com',
    description='Algorithms for combinatorics and number theory',
    long_description="""rosemary is a suite of algorithms for computations in number theory, discrete mathematics, and
                        combinatorial optimization""",
    include_package_data=True,
    packages=find_packages(
        exclude=(
            '.*',
            'EGG-INFO',
            '*.egg-info',
            '_trial*',
            '*.tests',
            '*.tests.*',
            'tests.*',
            'tests',
            'examples.*',
            'examples*',
        )
    ),
    tests_require=[
        'coverage',
        'mock',
        'nose',
    ],
    test_suite='nose.collector',
)
