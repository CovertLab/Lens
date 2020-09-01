import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

with open("requirements.txt", 'r') as requirements:
    install_requires = list(requirements.read().splitlines())

setup(
    name='wholecell-vivarium',
    version='0.0.54',
    packages=[
        'vivarium',
        'vivarium.analysis',
        'vivarium.core',
        'vivarium.compartments',
        'vivarium.parameters',
        'vivarium.plots',
        'vivarium.processes',
        'vivarium.reference_data',
        'vivarium.states',
        'vivarium.data',
        'vivarium.data.chromosomes',
        'vivarium.data.flat',
        'vivarium.data.flat.media',
        'vivarium.data.json_files',
        'vivarium.library'
    ],
    author='Eran Agmon, Ryan Spangler',
    author_email='eagmon@stanford.edu, spanglry@stanford.edu',
    url='https://github.com/CovertLab/vivarium',
    license='MIT',
    entry_points={
        'console_scripts': [
            'vivarium.analysis.analyze=vivarium.analysis.analyze:run',
            'vivarium.environment.boot=vivarium.environment.boot:run',
            'vivarium.environment.control=vivarium.environment.control:run',
            'vivarium.compartments=vivarium.compartments:run'
        ],
    },
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={
        'vivarium.data.flat': ['*.tsv', '*.fa'],
        'vivarium.data.flat.media': ['*.tsv'],
        'vivarium.reference_data': ['*.csv'],
        'vivarium.data.json_files': ['*.json']},
    include_package_data=True,
    install_requires=install_requires)
