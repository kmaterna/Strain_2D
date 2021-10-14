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
    description='This is the Strain_2D package',
    package_dir={
        'Strain_Tools': 'Strain_Tools',
        '': 'Strain_Tools'
    },
    packages=['Strain_Tools'] + find_packages('Strain_Tools'),
    scripts=[
        'Strain_Tools/bin/strain_rate_compute.py',
        'Strain_Tools/bin/strain_rate_comparison.py',
        'Strain_Tools/bin/moment_from_strain.py'
    ],
    zip_safe=False,
)
