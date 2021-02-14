#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Kathryn Materna, Jeremy Maurer
# Copyright 2021. ALL RIGHTS RESERVED.
# United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import re
from setuptools import find_packages, setup

# Parameter defs
CWD = os.getcwd()


def get_version():
    with open('version.txt', 'r') as f:
        m = re.match("""version=['"](.*)['"]""", f.read())

    assert m, "Malformed 'version.txt' file!"
    return m.group(1)


setup(
    name='Strain_2D',
    version=get_version(),
    description='This is the Strain_2d package',
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
