from distutils.core import setup
import ConfigParser
import sys
import os.path
from subprocess import check_output

config = ConfigParser.RawConfigParser()
config.read('setup.cfg')

if os.path.isdir('.git'):
    __version__ = check_output(["git", "describe", "--tags"])
    config.set('version', '__version__', __version__)
    with open('setup.cfg', 'wb') as configfile:
        config.write(configfile)
else:
    __version__ = config.get('version', '__version__')

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
