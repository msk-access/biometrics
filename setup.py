#!/usr/bin/env python

"""The setup script."""

import os

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


def req_file(filename):
    """
    We're using a requirements.txt file so that pyup.io can use this for security checks
    :param filename:
    :return str:
    """
    with open(filename) as f:
        content = f.readlines()
        content = filter(lambda x: not x.startswith("#"), content)
    return [x.strip() for x in content]

with open(os.path.join(os.path.dirname(__file__), "biometrics/VERSION"), "r") as fh:
    __version__ = fh.read().strip()


setup(
    author="Charlie Murphy",
    author_email='murphyc4@mskcc.org',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    description="Package to generate sample based biometrics.",
    entry_points={
        'console_scripts': [
            'biometrics=biometrics.cli:main',
        ],
    },
    install_requires=req_file("requirements.txt"),
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='biometrics',
    name='biometrics',
    packages=find_packages(include=['biometrics', 'biometrics.*']),
    package_data={
        "": ['requirements.txt', 'requirements_dev.txt'],
    },
    test_suite='tests',
    url='https://github.com/msk-access/biometrics',
    version=__version__,
    zip_safe=False,
)
