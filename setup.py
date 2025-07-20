"""
Modern setup.py for CMTJ using pybind11 helpers
"""
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "_cmtj",  # Rename to avoid conflict with the Python package
        sorted([
            "python/cmtj.cpp",
        ]),
        include_dirs=[
            # Path to core headers
            "core",
            # Path to third party headers
            "third_party",
        ],
        cxx_std=17,
        language="c++",
    ),
]

setup(
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    packages=find_packages(include=["cmtj", "cmtj.*"]),
    package_data={
        "cmtj": ["py.typed", "**/*.pyi"],
    },
    include_package_data=True,
    zip_safe=False,
)
