from setuptools import find_namespace_packages, setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os

__version__ = '1.0.0'
"""
As per 
https://github.com/pybind/python_example
"""

# if sys.platform == 'darwin':
#     # use g++ instead of clang on Mac
#     os.environ["CXX"] = "g++"
#     os.environ['CC'] = "g++"


class get_pybind_include(object):
    """
    Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. 
    """
    def __str__(self):
        import pybind11
        return pybind11.get_include()


ext_modules = [
    Extension(
        'cmtj',
        # Sort input source files to ensure bit-for-bit reproducible builds
        # (https://github.com/pybind/python_example/pull/53)
        sorted([os.path.join('python', 'cmtj.cpp')]),
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
        ],
        libraries=[],
        library_dirs=['/usr/local/lib'],
        extra_compile_args=['-O3', '-v', '-shared'],
        language='c++'),
]


# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """
    Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    import os
    with tempfile.NamedTemporaryFile('w', suffix='.cpp', delete=False) as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        fname = f.name
    try:
        compiler.compile([fname], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
        return False
    finally:
        try:
            os.remove(fname)
        except OSError:
            pass
    return True


def cpp_flag(compiler):
    """
    Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++17', '-std=c++14', '-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag):
            return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


class BuildExt(build_ext):
    """
    A custom build extension for adding compiler-specific options.
    """
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }
    """
    TBD if below is a problem for some. Leaving JIC
    """

    # if sys.platform == 'darwin':
    # darwin_opts = [, '-mmacosx-version-min=10.7']
    # c_opts['unix'] += darwin_opts
    # l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')

        for ext in self.extensions:
            ext.define_macros = [
                ('VERSION_INFO',
                 '"{}"'.format(self.distribution.get_version()))
            ]
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        build_ext.build_extensions(self)


setup(
    name='cmtj',
    version=__version__,
    author='Jakub',
    keywords=['magnetics', 'physics', 'simulation'],
    author_email="mojsieju@agh.edu.pl",
    url='https://github.com/LemurPwned/cmtj',
    description='CMTJ - C Magnetic Tunnel Junctions.',
    long_description='Efficient library for simulating magnetic multilayers',
    ext_modules=ext_modules,
    include_package_data=True,
    namespace_packages=['cmtj'],
    package_data={'cmtj': ["py.typed", "__init__.pyi"]},
    setup_requires=['pybind11>=2.6.1'],
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)