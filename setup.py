"""Setup configuration for imgtsero package."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="imgtsero",
    version="0.1.0",
    author="Tavi Nathanson",
    description="A Python library for downloading and working with IMGT/HLA data files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tavinathanson/imgtsero",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=[
        # No external dependencies - uses only Python standard library
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
        ],
    },
    entry_points={
        "console_scripts": [
            "imgtsero=imgtsero.__main__:main",
        ],
    },
    keywords="HLA IMGT immunogenetics bioinformatics",
    project_urls={
        "Bug Reports": "https://github.com/tavinathanson/imgtsero/issues",
        "Source": "https://github.com/tavinathanson/imgtsero",
    },
)