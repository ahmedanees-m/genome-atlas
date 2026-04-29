from setuptools import setup, find_packages

setup(
    name="genome-atlas",
    version="0.6.0",
    description="Programmable Genome-Writing Knowledge Graph (PEN-STACK Paper 1)",
    author="PEN-STACK Consortium",
    python_requires=">=3.10",
    packages=find_packages(exclude=["tests*", "notebooks*", "scripts*"]),
    install_requires=[
        "torch>=2.0",
        "torch-geometric>=2.3",
        "networkx>=3.0",
        "pandas>=2.0",
        "numpy>=1.24",
        "pyyaml>=6.0",
        "scikit-learn>=1.3",
        "click>=8.0",
    ],
    extras_require={
        "dev": ["pytest>=7.0", "pytest-cov>=4.0", "black>=23.0", "flake8>=6.0", "sphinx>=7.0"],
        "quantum": ["qiskit>=1.0", "qiskit-machine-learning>=0.7"],
        "all": [
            "pytest>=7.0", "pytest-cov>=4.0", "black>=23.0", "flake8>=6.0",
            "sphinx>=7.0", "qiskit>=1.0", "qiskit-machine-learning>=0.7",
        ],
    },
    entry_points={
        "console_scripts": ["genome-atlas=genome_atlas.cli:cli"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
