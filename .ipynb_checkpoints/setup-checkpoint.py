from setuptools import setup, find_packages
import sys
import os
import contextlib

with open(os.devnull, 'w') as nullfile:
    with contextlib.redirect_stderr(nullfile):
        setup(
            name="graphimate",
            version="0.1.0",
            packages=find_packages(),
            install_requires=[
                "scanpy==1.9.3",
                "pandas==1.5.3",
                "numpy==1.22.1",
                "scipy==1.7.3",
                "matplotlib==3.7.1",
                "seaborn==0.12.2",
                "requests==2.31.0",
                "psutil==5.9.5",
                "Pillow==9.5.0",
                "moviepy",
                "fa2",
                "igraph",
            ],
        )