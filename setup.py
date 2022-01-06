from setuptools import setup, find_packages

# Get the documentation
with open("README.md", "r") as fh:
    long_description = fh.read()

# Get the requirements
with open('requirements_frozen.txt') as f: required_frozen = f.read().splitlines()
with open('requirements_web.txt')    as f: required_web    = f.read().splitlines()

CLASSIFIERS = [
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: Other/Proprietary License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Development Status :: 4 - Beta",
    "Programming Language :: Python :: 3.8",
]

setup(
    name="controlling_covid19_ttq",
    version="1.2.1",
    author="Cliff Kerr, Dina Mistry, Robyn Stuart, Daniel Klein, et al., on behalf of the IDM COVID-19 Response Team",
    author_email="info@covasim.org",
    description="Code for the 'Controlling COVID-19 via test-trace-quarantine' paper",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='http://ttq.covasim.org',
    keywords=["COVID-19", "SARS-CoV-2", "testing", "tracing", "quarantine", "covasim"],
    platforms=["OS Independent"],
    classifiers=CLASSIFIERS,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "statsmodels",
        "matplotlib",
        "seaborn",
        "sciris",
        "covasim",
    ],
    extras_require={
        "frozen": required_frozen,
        "web":    required_web,
    }
)
