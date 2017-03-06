from setuptools import setup, find_packages

setup(
    name="hostdesigner",
    version="0.1",
    description="HostDesigner python wrapper",
    author="Kutay B. Sezginel",
    author_email="kbs37@pitt.edu",
    install_requires=[
        'numpy',
        'pyyaml',
        'tabulate',
        'nglview'
    ],
    include_package_data=True,
    packages=find_packages()
)
