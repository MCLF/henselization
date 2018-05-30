from setuptools import setup, find_packages

setup(
    name="henselization",
    packages=find_packages(include=['henselization*']),
    install_requires=['recursive-monkey-patch']
)


