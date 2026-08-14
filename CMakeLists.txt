cmake_minimum_required(VERSION 3.23)

# -----------------------------------------------------------------------------
# Compiler: Clang++
# -----------------------------------------------------------------------------
if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER clang++)
endif()

project(main_driver VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -----------------------------------------------------------------------------
# C++23
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# macOS / Apple Silicon (aarch64). Apple's own tooling calls this "arm64";
# CMAKE_OSX_ARCHITECTURES is the variable that actually controls it.
# -----------------------------------------------------------------------------
if(APPLE AND NOT CMAKE_OSX_ARCHITECTURES)
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architecture" FORCE)
endif()

# -----------------------------------------------------------------------------
# Build type: defaults to Release. Select explicitly with, e.g.:
#   cmake -S . -B build-debug   -DCMAKE_BUILD_TYPE=Debug   && cmake --build build-debug
#   cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release
# -----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# -----------------------------------------------------------------------------
# Dependencies
#   brew install fftw openblas zlog
#
# openblas is keg-only on Homebrew (macOS already ships Accelerate.framework
# as its BLAS/LAPACK, so brew won't symlink openblas over it), so its
# pkg-config file isn't on PKG_CONFIG_PATH by default -- add it explicitly.
# -----------------------------------------------------------------------------
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH /opt/homebrew)
    set(ENV{PKG_CONFIG_PATH} "/opt/homebrew/opt/openblas/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()

find_package(PkgConfig REQUIRED)

pkg_check_modules(FFTW3    REQUIRED IMPORTED_TARGET fftw3)    # double precision
pkg_check_modules(FFTW3F   REQUIRED IMPORTED_TARGET fftw3f)   # single precision (float)
pkg_check_modules(OPENBLAS REQUIRED IMPORTED_TARGET openblas)

find_package(Threads REQUIRED)

# zlog ships neither a pkg-config file nor a CMake package config (verified
# against its Homebrew formula), so it's located by hand. Override with
# -DZLOG_INCLUDE_DIR=... -DZLOG_LIBRARY=... if it lives somewhere unusual.
find_path(ZLOG_INCLUDE_DIR NAMES zlog.h)
find_library(ZLOG_LIBRARY NAMES zlog)

if(NOT ZLOG_INCLUDE_DIR OR NOT ZLOG_LIBRARY)
    message(FATAL_ERROR "zlog not found. Run `brew install zlog`, or pass -DZLOG_INCLUDE_DIR / -DZLOG_LIBRARY explicitly.")
endif()

add_library(zlog::zlog UNKNOWN IMPORTED)
set_target_properties(zlog::zlog PROPERTIES
    IMPORTED_LOCATION "${ZLOG_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${ZLOG_INCLUDE_DIR}"
)

# -----------------------------------------------------------------------------
# Target -- replace main.cpp with your actual sources
# -----------------------------------------------------------------------------
add_executable(${PROJECT_NAME}
    src/common/symbols_modulation/symbols_modulation.cpp
    main_driver.cpp
)

target_compile_options(${PROJECT_NAME} PRIVATE 
    -Wall -Wextra
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    PkgConfig::FFTW3
    PkgConfig::FFTW3F
    PkgConfig::OPENBLAS
    zlog::zlog
    Threads::Threads
)