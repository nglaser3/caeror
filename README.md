# Caeror

## Overview

Caeror is a library for ray tracing through constructive solid geometry (CSG) representations, targeting CPU and GPU backends.
It is intended for Monte Carlo particle transport and deterministic transport methods like the method of characteristics (MOC) and the random ray method (TRRM).

## Build

Requires a C++23-capable compiler, CMake 3.25, and a build tool such as Make or Ninja. BLT and RAJA are included as Git submodules.

Clone the repository and its dependencies:

```sh
git clone https://github.com/scram-lab/caeror.git
cd caeror
git submodule update --init --recursive
```

Configure and build:

```sh
cmake --workflow --preset release
```

The build is generated in `build/release`. Tests are enabled by default; run them with:

```sh
ctest --preset release
```
