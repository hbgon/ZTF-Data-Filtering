from setuptools import setup, find_packages

requirements = [
    "numpy",
    "pandas",
    "tqdm",
    "astropy",
    "matplotlib",
    "astroquery"
]

setup(
    name="ZTF-Data-Filter", 
    version="1.2.0",
    description="A project for ZTF data filtering and analysis",
    author="Hygor B. Gonçalves",
    author_email="hygor.benati@ufrgs.br",
    packages=find_packages(),  
    install_requires=requirements, 
    entry_points={
        "console_scripts": [
            "ztfdatafiltering=ztfdatafiltering.main:main",  # Entry point
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)