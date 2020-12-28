######################################################################
#                                                                    #
#  Copyright 2009 Lucas Heitzmann Gabrielli.                         #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import sys
import platform
from setuptools import setup, Extension
from distutils.version import LooseVersion

with open("README.md") as fin:
    long_description = fin.read()

with open("gdspy/__init__.py") as fin:
    for line in fin:
        if line.startswith("__version__ ="):
            version = eval(line[14:])
            break

setup_requires = []
if {"pytest", "test", "ptr"}.intersection(sys.argv):
    setup_requires.append("pytest-runner")
if "build_sphinx" in sys.argv:
    setup_requires.extend(["sphinx", "sphinx_rtd_theme"])

# Mac OS X Mojave C++ compile + linking arguments
extra_compile_args = []
extra_link_args = []
if platform.system() == "Darwin" and LooseVersion(platform.release()) >= LooseVersion(
    "17.7"
):
    extra_compile_args = ["-std=c++11", "-mmacosx-version-min=10.9"]
    extra_link_args = ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

setup(
    name="gdspy",
    version=version,
    author="Lucas Heitzmann Gabrielli",
    author_email="heitzmann@gmail.com",
    license="Boost Software License v1.0",
    url="https://github.com/heitzmann/gdspy",
    description="Python module for creating/importing/merging GDSII files.",
    long_description=long_description,
    keywords="GDSII CAD layout",
    packages=["gdspy"],
    package_dir={"gdspy": "gdspy"},
    package_data={"gdspy": ["data/*"]},
    ext_modules=[
        Extension(
            "gdspy.clipper",
            ["gdspy/clipper.cpp"],
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
        )
    ],
    provides=["gdspy"],
    install_requires=["numpy"] + (["future"] if sys.version_info.major < 3 else []),
    setup_requires=setup_requires,
    tests_require=["pytest"],
    platforms="OS Independent",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Manufacturing",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Boost Software License 1.0 (BSL-1.0)",
        "Operating System :: OS Independent",
        "Programming Language :: C",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Electronic Design Automation (EDA)",
    ],
    zip_safe=False,
)
