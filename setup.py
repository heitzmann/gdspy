######################################################################
#                                                                    #
#  Copyright 2009-2017 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import sys
from setuptools import setup, Extension

with open('README.md') as fin:
    long_description = fin.read()

with open('gdspy/__init__.py') as fin:
    for line in fin:
        if line.startswith('__version__ ='):
            version = eval(line[14:])
            break

setup_requires = ['pytest-runner'] if \
    {'pytest', 'test', 'ptr'}.intersection(sys.argv) else []

setup(
    name='gdspy',
    version=version,
    author='Lucas Heitzmann Gabrielli',
    author_email='heitzmann@gmail.com',
    license='Boost Software License v1.0',
    url='https://github.com/heitzmann/gdspy',
    description='Python module for creating/importing/merging GDSII files.',
    long_description=long_description,
    keywords='GDSII CAD layout',
    packages=['gdspy'],
    package_dir={'gdspy': 'gdspy'},
    package_data={'gdspy': ['data/*']},
    ext_modules=[
        Extension('gdspy.boolext', ['gdspy/boolext.c']),
        Extension('gdspy.clipper', ['gdspy/clipper.cpp'])
    ],
    provides=['gdspy'],
    install_requires=['numpy'] + (['future']
                                  if sys.version_info.major < 3 else []),
    setup_requires=setup_requires,
    tests_require=['pytest'],
    platforms='OS Independent',
    classifiers=[
        'Development Status :: 4 - Beta', 'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved',
        'Operating System :: OS Independent', 'Programming Language :: C',
        'Programming Language :: C++', 'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Electronic Design Automation (EDA)'
    ],
    zip_safe=False)
