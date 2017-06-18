#  Copyright 2016-2017 Alan F Rubin
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

requirements = [
    'numpy >= 1.10.4',
    'scipy >= 0.16.0',
    'pandas >= 0.18.0',
    'pandas < 0.20',
    'statsmodels >= 0.6.1',
    'matplotlib >= 1.4.3',
]

# workaround to deal with Enthought Canopy referring to the tables module
# as pytables
if 'Canopy' in sys.executable:
    requirements.append('pytables >= 3.2.0')
else:
    requirements.append('tables >= 3.2.0')

setup(
    name="Enrich2",
    version="1.1.0a",

    packages=find_packages(),
    package_data={
        'test': ['test_files/*/*'],
    },

    entry_points={
        'console_scripts': ['enrich_cmd = enrich2.main:main_cmd'],
        'gui_scripts': ['enrich_gui = enrich2.main:main_gui'],
    },

    author="Alan F Rubin",
    author_email="alan.rubin@wehi.edu.au",
    description="Analysis program for calculating variant scores from "
    "deep mutational scanning data.",
    url="https://github.com/FowlerLab/Enrich2/",

    install_requires=requirements,
)
