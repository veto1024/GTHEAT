#!/usr/bin/python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
deps = ['matplotlib',
        'gt3',
        'scipy',
        'PyQT5']

setuptools.setup(
    name="GTHEAT",
    version="0.0.2",
    author="Jonathan J. Roveto",
    install_requires=deps,
    include_package_data=True,
    author_email="veto1024@gmail.com",
    description="GTHEAT - The Georgia Tech Heat and Particle Transport Analysis Code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/veto1024/gtheat/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    python_requires='>=3.8',
)

if __name__ == '__main__':
    pass
