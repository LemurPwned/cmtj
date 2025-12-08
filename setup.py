"""
Modern setup.py for CMTJ using pybind11 helpers
"""
import os
import sys
from pathlib import Path
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages

# Get version from setuptools_scm
try:
    from setuptools_scm import get_version
    version = get_version(root=".", relative_to=__file__)
except (ImportError, LookupError):
    # Fallback: try to read from _version.py if it exists
    try:
        from cmtj._version import version
    except ImportError:
        version = "dev"

# Handle platform-specific compilation issues
extra_compile_args = []
extra_link_args = []
define_macros = [("VERSION_INFO", f'"{version}"')]


if 'CXXFLAGS' not in os.environ:
    extra_compile_args.extend(['-O2']) # we avoid -ffastmath and native here
if 'LDFLAGS' not in os.environ:
    extra_link_args.extend(['-O2'])

if sys.platform == 'darwin':
    # macOS-specific flags
    extra_compile_args.extend([
        '-mmacosx-version-min=10.9',  # Ensure compatibility
        '-stdlib=libc++',  # Use libc++ on macOS
        '-fvisibility=hidden',  # Hide symbols by default
        '-Wno-unused-variable',  # Suppress warnings that might be errors
        '-Wno-deprecated-declarations',  # Handle deprecated API warnings
    ])
elif sys.platform == 'win32':
    # Windows-specific flags
    extra_compile_args.extend([
        '/EHsc',  # Enable exception handling
        '/bigobj',  # Allow large object files
        '/wd4244',  # Disable conversion warnings
        '/wd4267',  # Disable size_t conversion warnings
    ])
    define_macros.extend([
        ('WIN32_LEAN_AND_MEAN', None),
        ('NOMINMAX', None),  # Prevent Windows.h from defining min/max macros
        ('_USE_MATH_DEFINES', None),  # Enable M_PI and other math constants
    ])
else:
    # Linux/Unix flags
    extra_compile_args.extend([
        '-fvisibility=hidden',
        '-Wno-unused-variable',
    ])

# Define the extension module with cross-platform include paths
# Use absolute paths and ensure they exist
try:
    project_root = Path(__file__).parent.resolve()
except NameError:
    # Fallback when __file__ is not available (e.g., in some build environments)
    project_root = Path(".").resolve()

include_dirs = [
    str(project_root),  # Project root for relative includes
    str(project_root / "core"),  # Core headers
    str(project_root / "third_party"),  # Third party headers root
    str(project_root / "third_party" / "kissfft"),  # Specific kissfft path
]

# Verify all include directories exist
for inc_dir in include_dirs:
    if not Path(inc_dir).exists():
        print(f"Warning: Include directory does not exist: {inc_dir}")

ext_modules = [
    Pybind11Extension(
        "_cmtj",  # Internal C++ extension name
        sorted([
            "python/cmtj.cpp",  # Use relative path as required by setuptools
        ]),
        include_dirs=include_dirs,
        define_macros=define_macros,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
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
