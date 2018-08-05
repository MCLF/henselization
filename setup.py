from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="henselization",
    version="0.0.0",
    author="Julian Rüth",
    author_email="julian.rueth@fsfe.org",
    description="Henselizations in Sage",
    url="https://github.com/mclf/henselization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(include=['henselization*']),
    install_requires=['patchy', 'recursive-monkey-patch==0.4.0'],
    classifiers=(
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Operating System :: OS Independent")
)
