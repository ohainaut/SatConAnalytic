from setuptools import setup, find_packages

setup(
    name='SatConAnalytic',  # Name of your package
    version='2.0.4',   # Version number
    description='Effects of satellite constellations on astronomical observations',
    long_description=open('README.md').read(),  # Long description (from README)
    long_description_content_type='text/markdown',  # Format of README
    author='Olivier Hainaut',
    author_email='ohainaut@eso.org',
    url='https://github.com/ohainaut/SatConAnalytic',  # Project URL (if applicable)
    packages=find_packages(),  # Automatically find all packages
    include_package_data=True,  # Include files from MANIFEST.in
    install_requires=[  ],  # Add any dependencies here (e.g., ['numpy', 'requests'])
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',  # Minimum Python version
)