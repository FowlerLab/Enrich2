#  Copyright 2016-2019 Alan F Rubin
#
#  This file is part of Enrich2.
#
#  Enrich2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Enrich2 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Enrich2.  If not, see <http://www.gnu.org/licenses/>.

import sys
from setuptools import setup, find_packages

setup(
    name="Enrich2",
    version="1.2.1",

    packages=find_packages(),

    entry_points={
        'console_scripts': ['enrich_cmd = enrich2.main:main_cmd'],
        'gui_scripts': ['enrich_gui = enrich2.main:main_gui'],
    },

    author="Alan F Rubin",
    author_email="alan.rubin@wehi.edu.au",
    description="Analysis program for calculating variant scores from "
    "deep mutational scanning data.",
    url="https://github.com/FowlerLab/Enrich2/",

    install_requires=[
        'numpy >= 1.10.4',
        'scipy >= 0.16.0 ',
        'pandas >= 0.18.0',
        'statsmodels >= 0.6.1',
        'matplotlib >= 1.4.3',
        'tables >= 3.2.0',
    ]
)
