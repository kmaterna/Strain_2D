#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert, Jeremy Maurer, and Piyush Agram
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import glob
import os
import re

import numpy as np
from setuptools import Extension, find_packages, setup

# Parameter defs
CWD = os.getcwd()


def get_version():
    with open('version.txt', 'r') as f:
        m = re.match("""version=['"](.*)['"]""", f.read())

    assert m, "Malformed 'version.txt' file!"
    return m.group(1)


setup(
    name='2D_Strain',
    version=get_version(),
    description='This is the 2D_Strain package',
    package_dir={
        'tools': 'tools',
        '': 'tools'
    },
    packages=['tools'] + find_packages('tools'),
    scripts=[
        'tools/bin/strain_driver.py',
        'tools/bin/compare_driver.py',
    ],
    zip_safe=False,
)
