"""
Setup script for scReservoir.

Installs both the classical (Echo State Network) and quantum reservoir
computing packages for single-cell RNA-seq GRN inference.

Usage
-----
  pip install -e .                    # editable install (development)
  pip install -e ".[quantum]"         # include quantum extras
  pip install -e ".[dev]"             # include dev/test tools
  pip install -e ".[quantum,dev]"     # everything
"""

from setuptools import setup, find_packages
from pathlib import Path

long_description = Path('README.md').read_text(encoding='utf-8')

setup(
    name='scReservoir',
    version='0.2.0',
    author='scReservoir Contributors',
    author_email='contact@screservoir.org',
    description=(
        'Classical and Quantum Reservoir Computing for '
        'Single-Cell RNA-seq GRN Inference'
    ),
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/screservoir/screservoir',

    # -------------------------------------------------------------------------
    # Package discovery
    # The repo has two flat namespaces:
    #   code/classical/   — Echo State Network pipeline
    #   code/quantum/     — Quantum / NG-RC pipeline
    # Both are plain directories with __init__.py so Python finds them as
    # packages when you add code/ to sys.path (or install in editable mode).
    # -------------------------------------------------------------------------
    packages=find_packages(where='.'),
    package_dir={'': '.'},

    python_requires='>=3.8',

    # -------------------------------------------------------------------------
    # Core dependencies — needed for classical reservoir pipeline
    # -------------------------------------------------------------------------
    install_requires=[
        'numpy>=1.21',
        'scipy>=1.9',          # bumped: expm used in quantum_reservoir.py
        'scikit-learn>=1.0',
        'scanpy>=1.9',
        'anndata>=0.8',
        'matplotlib>=3.5',
        'seaborn>=0.12',
    ],

    # -------------------------------------------------------------------------
    # Optional extras
    # -------------------------------------------------------------------------
    extras_require={
        # quantum module uses scipy.linalg.expm (already in core) — no new
        # hard deps, but future versions may add qiskit / pennylane hooks
        'quantum': [
            'scipy>=1.9',      # expm, already required above
        ],
        'dev': [
            'pytest>=7.0',
            'pytest-cov>=4.0',
            'black>=23.0',
            'flake8>=6.0',
            'isort>=5.0',
        ],
    },

    include_package_data=True,

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Physics',
    ],

    keywords=[
        'single-cell', 'RNA-seq', 'scRNA-seq',
        'gene-regulatory-networks', 'GRN',
        'reservoir-computing', 'echo-state-networks',
        'quantum-computing', 'quantum-reservoir-computing',
        'pseudotime', 'developmental-biology',
        'Waddington-landscape', 'attractor',
    ],

    project_urls={
        'Bug Reports': 'https://github.com/screservoir/screservoir/issues',
        'Source':      'https://github.com/screservoir/screservoir',
        'Documentation': 'https://screservoir.readthedocs.io',
    },
)
