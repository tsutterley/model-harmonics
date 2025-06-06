import os
from setuptools import setup, find_packages

# package description and keywords
description = ('Python tools for obtaining and working with model '
    'synthetic spherical harmonic coefficients for comparing with '
    'data from the NASA/DLR GRACE and NASA/GFZ GRACE Follow-on missions')
keywords = 'gravity synthetics, physical geodesy, spherical harmonics'
# get long_description from README.rst
with open("README.rst", mode='r', encoding='utf8') as fh:
    long_description = fh.read()
long_description_content_type = "text/x-rst"

# get install requirements
with open('requirements.txt', encoding='utf8') as fh:
    install_requires = fh.read().splitlines()

# get version
with open('version.txt', encoding='utf8') as fh:
    version = fh.read()

# list of all scripts to be included with package
scripts=[]
for s in ['scripts','GIA','OBP','reanalysis','seismic','SMB','TWS']:
    scripts.extend([os.path.join(s,f) for f in os.listdir(s) if f.endswith('.py')])

setup(
    name='model-harmonics',
    version=version,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    url='https://github.com/tsutterley/model-harmonics',
    author='Tyler Sutterley',
    author_email='tsutterl@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords=keywords,
    packages=find_packages(),
    install_requires=install_requires,
    scripts=scripts,
)
