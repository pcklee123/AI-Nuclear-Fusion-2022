# AI Nuclear Fusion 2022
## Table of Contents
1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
3. [Usage](#usage)

## Introduction
This project simulates Magnetized Liner Inertial Fusion (MagLIF). [Specifically, it simulates charged particles and their trajectories in β and e fields.]?

## Getting Started
### Prerequisites
- MSYS2
- GCC added to PATH in MSYS
    - In the root directory, run `export PATH=$PATH:/mingw64/bin`
<<<<<<< HEAD
    - this can also be added in windows "edit system environment variables" , "Path". Add C:\msys64\usr\bin and C:\msys64\....\bin

To install fftw3 recompiled with OMP enabled:
```
> wget https://www.fftw.org/fftw-3.3.10.tar.gz
> tar xvzf fftw-3.3.10.tar.gz
> cd fftw-3.3.10/
> ./configure --enable-threads --enable-openmp --enable-avx --enable-avx2 --enable-avx512 --enable-avx-128-fma --enable-float --with-our-malloc --enable-sse2
> make
> make install
```

To connect fftw3.h to VSCode, add the path of where fftw3.h is installed (fftw-3.3.10/api)

Add before fftw plans (what does this mean??)

`LIBS= -lm -lgsl -lOpenCL.dll -lomp.dll -lfftw3f -lfftw3f_omp`
=======
- Required Libraries
    - GSL
## Usage
>>>>>>> 25aed6a063a0bea6499efba7407cd9b89ecfc722
