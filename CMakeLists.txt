# compas_cra._native: IPOPT bound into a Python extension module.
#
# IPOPT (and its MUMPS linear solver) are consumed as the static libraries produced by
# packaging/build_ipopt.sh; point IPOPT_PREFIX (env or cache var) at that stage tree.
# The result is one self-contained extension module inside the compas_cra package: no
# ipopt executable, no DLLs, no separate solver distribution.

cmake_minimum_required(VERSION 3.18...3.30)
project(compas_cra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT DEFINED IPOPT_PREFIX)
    if(DEFINED ENV{IPOPT_PREFIX})
        set(IPOPT_PREFIX "$ENV{IPOPT_PREFIX}")
    else()
        # the default location packaging/build_ipopt.sh stages into
        set(IPOPT_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/build/ipopt/stage")
    endif()
endif()
if(NOT EXISTS "${IPOPT_PREFIX}/include/coin-or/IpTNLP.hpp")
    message(FATAL_ERROR "No IPOPT build at IPOPT_PREFIX=${IPOPT_PREFIX}; run `invoke setup` (or packaging/build_ipopt.sh) first")
endif()
message(STATUS "IPOPT_PREFIX: ${IPOPT_PREFIX}")

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# nanobind is a build requirement; locate its CMake config through the interpreter
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_VARIABLE nanobind_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(nanobind CONFIG REQUIRED)

# IPOPT + MUMPS themselves are linked as the static archives from the stage tree. The
# Fortran/BLAS runtimes are linked dynamically instead: static Fortran runtimes do not
# link cleanly into shared objects (local symbols, non-PIC objects), and the ecosystem
# answer — what numpy/scipy wheels do — is to link the shared runtime and have the
# platform wheel-repair tool (delocate / auditwheel / delvewheel) graft it into the
# wheel afterwards.
# everything needed is derivable from the stage tree itself; a pkg-config executable
# is not reliably present on all builders (Windows in particular)
set(IPOPT_INCLUDE_DIRS "${IPOPT_PREFIX}/include/coin-or")
file(STRINGS "${IPOPT_PREFIX}/lib/pkgconfig/ipopt.pc" _ipopt_ver_line REGEX "^Version:")
string(REGEX REPLACE "^Version: *" "" IPOPT_VERSION "${_ipopt_ver_line}")
message(STATUS "IPOPT version: ${IPOPT_VERSION}")

nanobind_add_module(_core native/src/ipopt_nb.cpp)

target_include_directories(_core PRIVATE ${IPOPT_INCLUDE_DIRS})
target_link_libraries(_core PRIVATE
    "${IPOPT_PREFIX}/lib/libipopt${CMAKE_STATIC_LIBRARY_SUFFIX}"
    "${IPOPT_PREFIX}/lib/libcoinmumps${CMAKE_STATIC_LIBRARY_SUFFIX}")
target_compile_definitions(_core PRIVATE IPOPT_VERSION="${IPOPT_VERSION}")

# BLAS/LAPACK: Accelerate on macOS, OpenBLAS elsewhere (matching build_ipopt.sh)
if(APPLE)
    target_link_options(_core PRIVATE "-Wl,-framework,Accelerate")
else()
    target_link_libraries(_core PRIVATE openblas)
endif()

# the shared Fortran runtime, located through the Fortran compiler itself
if(NOT DEFINED FORTRAN_COMPILER)
    set(FORTRAN_COMPILER gfortran)
endif()
if(APPLE)
    set(_frt libgfortran.dylib)
elseif(WIN32)
    set(_frt libgfortran.dll.a)
else()
    set(_frt libgfortran.so)
endif()
execute_process(
    COMMAND ${FORTRAN_COMPILER} -print-file-name=${_frt}
    OUTPUT_VARIABLE _frt_path OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
if(NOT EXISTS "${_frt_path}")
    message(FATAL_ERROR "no shared Fortran runtime found (${FORTRAN_COMPILER} -print-file-name=${_frt} -> ${_frt_path})")
endif()
get_filename_component(_frt_dir "${_frt_path}" DIRECTORY)
message(STATUS "Fortran runtime: ${_frt_path}")
target_link_directories(_core PRIVATE "${_frt_dir}")
target_link_libraries(_core PRIVATE gfortran)
if(NOT WIN32)
    target_link_libraries(_core PRIVATE ${CMAKE_DL_LIBS})
endif()
if(WIN32)
    target_link_libraries(_core PRIVATE quadmath)
endif()

# extra link directories/flags for special setups (e.g. the msys2 lib dir on Windows)
if(DEFINED ENV{EXTRA_LINK_DIRS})
    target_link_directories(_core PRIVATE "$ENV{EXTRA_LINK_DIRS}")
endif()
if(DEFINED ENV{IPOPT_EXTRA_LINK})
    separate_arguments(_extra NATIVE_COMMAND "$ENV{IPOPT_EXTRA_LINK}")
    target_link_options(_core PRIVATE ${_extra})
endif()

# arm64 macOS kills processes that load unsigned code, so make sure the module keeps a
# valid ad-hoc signature whatever earlier steps did to it
if(APPLE)
    add_custom_command(TARGET _core POST_BUILD
        COMMAND codesign --force -s - $<TARGET_FILE:_core>
        COMMENT "ad-hoc signing $<TARGET_FILE:_core>")
endif()

install(TARGETS _core LIBRARY DESTINATION compas_cra/_native)
