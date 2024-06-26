# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(pybind11_DIR /Users/ktocz/AppData/Roaming/Python/Python312/site-packages/pybind11/share/cmake/pybind11)

#ignorowanie flagami
#add_compile_options(-Wno-ignored-attributes)

add_subdirectory(matplotplusplus)


# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
python_add_library(_core MODULE src/main.cpp WITH_SOABI)


target_link_libraries(_core PRIVATE pybind11::headers)

target_link_libraries(_core PUBLIC matplot)



# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION scikit_build_example)
