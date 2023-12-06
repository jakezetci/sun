# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
from pybind11.setup_helpers import ParallelCompile

# Optional multithreaded build
ParallelCompile("NPY_NUM_BUILD_JOBS").install()

__version__ = "0.0.3"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/cpp_module/pull/53)

ext_modules = [
    Pybind11Extension(
        "cpp_module",
        ["src/main.cpp"],
        # Example: passing in the version to the compiled code
        define_macros=[("VERSION_INFO", __version__)],
    ),
]

setup(
    name="cpp_module",
    version=__version__,
    author="Sylvain Corlay",
    author_email="sylvain.corlay@gmail.com",
    url="https://github.com/pybind/cpp_module",
    description="A test project using pybind11",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
