from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="hostdesigner",
    version="0.1",
    description="HostDesigner python wrapper",
    author="Kutay B. Sezginel",
    author_email="kbs37@pitt.edu",
    install_requires=requirements,
    include_package_data=True,
    packages=find_packages()
)
