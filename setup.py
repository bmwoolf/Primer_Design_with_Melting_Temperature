from setuptools import setup, find_packages

setup(
    name='primer_design',
    version='0.1.0',
    description='PCR primer design tool with melting temperature constraints',
    author='Bradley Woolf',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.81',
        'numpy>=1.21.0',
    ],
    entry_points={
        'console_scripts': [
            'primer-design = primer_design.cli:main',
        ],
    },
    python_requires='>=3.7',
    include_package_data=True,
    license='MIT',
) 