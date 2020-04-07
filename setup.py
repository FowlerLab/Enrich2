import sys
from setuptools import setup, find_packages

setup(
    name="Enrich2",
    version="1.3.1",

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
