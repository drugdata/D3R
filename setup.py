#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    "argparse",
    "lockfile",
    "psutil",
    "biopython",
    "xlsxwriter",
    "ftpretty",
    "python-dateutil",
    "easywebdav"
]

test_requirements = [
    "argparse",
    "lockfile",
    "psutil",
    "biopython",
    "xlsxwriter",
    "ftpretty",
    "python-dateutil",
    "mock",
    "easywebdav"
]

setup(
    name='d3r',
    version='1.6.8',
    description='Drug Design Data Resource CELPP Runner is an application to run the filtering, docking, and '
                'evaluation of new sequences from wwpdb',
    long_description=readme + '\n\n' + history,
    author='Christopher Churas <churas@ncmir.ucsd.edu>,'
           'Shaui Liu <shuailiu25@gmail.com>,'
           'Jeff Wagner <j5wagner@ucsd.edu>,'
           'Rob Swift <rvswift@ucsd.edu>',
    author_email='drugdesigndata@gmail.com',
    url='https://github.com/drugdata/D3R',
    packages=[
        'd3r', 'd3r.blast', 'd3r.filter', 'd3r.utilities', 'd3r.celpp', 
        'd3r.celppade'
    ],
    package_dir={'d3r':
                 'd3r'},
    include_package_data=True,
    install_requires=requirements,
    license="Other",
    zip_safe=False,
    keywords='d3r',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        'Environment :: Console',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    scripts = ['d3r/celpprunner.py', 'd3r/blastnfilter.py',
               'd3r/postanalysis.py', 'd3r/proteinligprep.py',
               'd3r/glidedocking.py', 'd3r/evaluate.py',
               'd3r/celppreports.py', 'd3r/vinadocking.py',
               'd3r/genchallengedata.py','d3r/chimera_proteinligprep.py',
               'd3r/getchallengedata.py','d3r/packdockingresults.py'
               ],
    test_suite='tests',
    tests_require=test_requirements
)
