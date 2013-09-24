from distutils.core import setup
import sys
import os.path
from subprocess import check_output

if os.path.isdir('.git'):
    __version__ = check_output(["git", "describe", "--tags"])
    with open('gemini/__init__.py', 'w') as o:
        o.write("__version__ = '{0}'".format(__version__.strip()))
else:
    with open('gemini/__init__.py', 'w') as o:
        o.write("__version__ = '{0}'".format('0.1'))

setup(
        name = 'fastqp',
        provides = 'fastqp',
        version = __version__,
        author = 'Matthew Shirley',
        author_email = 'mdshw5@gmail.com',
        url = 'http://mattshirley.com',
        description = 'simple NGS read quality assessment using Python',
        license = 'MIT',
        py_modules = ['fastqp'],
        scripts = ['scripts/fastqp'],
        classifiers = [
                "Development Status :: 3 - Alpha",
                "License :: OSI Approved :: MIT License",
                "Environment :: Console",
                "Intended Audience :: Science/Research",
                "Natural Language :: English",
                "Operating System :: Unix",
                "Programming Language :: Python :: 2.7",
                "Topic :: Scientific/Engineering :: Bio-Informatics"
        ]
)
