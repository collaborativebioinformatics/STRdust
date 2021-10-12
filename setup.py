# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

exec(open('STRdust/version.py').read())

setup(
    name='STRdust',
    version=__version__,
    description='Discovery of STR lengths from long reads',
    long_description=open(path.join(here, "README.md")).read(),
    long_description_content_type="text/markdown",
    url='https://github.com/collaborativebioinformatics/STRdust',
    author='hackathon 2021 group2',
    author_email='decosterwouter@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='STR long reads',
    packages=find_packages(),
    python_requires='>=3',
    install_requires=['pysam'],
    package_data={'STRdust': []},
    package_dir={'STRdust': 'STRdust'},
    include_package_data=True,
    entry_points={
        'console_scripts': ['STRdust=STRdust.STRdust:main', ]}
)
