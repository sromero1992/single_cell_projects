"""
Setup script for scReservoir package.
"""

from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='scReservoir',
    version='0.1.0',
    author='scReservoir Contributors',
    author_email='contact@screservoir.org',
    description='Reservoir Computing Framework for Single-Cell RNA-seq GRN Inference',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/screservoir/screservoir',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.21',
        'scipy>=1.7',
        'scikit-learn>=1.0',
        'scanpy>=1.9',
        'anndata>=0.8',
        'matplotlib>=3.5',
        'seaborn>=0.11',
    ],
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=22.0',
            'flake8>=3.9',
        ],
    },
    include_package_data=True,
    keywords=['single-cell', 'RNA-seq', 'gene-regulatory-networks', 'reservoir-computing', 'echo-state-networks'],
    project_urls={
        'Bug Reports': 'https://github.com/screservoir/screservoir/issues',
        'Source': 'https://github.com/screservoir/screservoir',
        'Documentation': 'https://screservoir.readthedocs.io',
    },
)
