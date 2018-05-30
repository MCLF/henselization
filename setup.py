from setuptools import setup, find_packages

setup(
    name="henselization",
    packages=find_packages(include=['henselization*']),
    install_requires=['recursive-monkey-patch==0.3.0'],
    dependency_links=['https://github.com/saraedum/recursive-monkey-patch/tarball/import#egg=recursive-monkey-patch-0.3.0']
)


