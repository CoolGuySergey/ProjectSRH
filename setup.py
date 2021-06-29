from setuptools import setup



with open("README.md", "r") as fh:
    long_description = fh.read


    
setup(

    
    name = "SRHClusterMapper",
    version = "0.01",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    py_modules = [
        "functions",
        "SRHClusterMapper",
    ],
    package_dir = {"" : "src"},

    
    install_requires = [
        "scipy >= 1.6",
        "numpy >= 1.1",
        "seaborn >= 0.11",
        "matplotlib >= 3.4",
    ],
    extras_requires = {
        "dev": [
            "pytest == 6.2.4",
        ],
    },

    
    classifiers = [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
)
