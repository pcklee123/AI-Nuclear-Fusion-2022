# AI Nuclear Fusion 2022
## Table of Contents
1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
3. [Usage](#usage)
3. [Todo](#Todo)

## Introduction
This project aims to develop a Particle in Cell plasma code.
Original code in 2021 by Hilary,Yin Yue and Chloe, extensive improvements by Samuel,Ananth and Vishwa.

## Getting Started
### Prerequisites
- MSYS2
- to avoid confusion, use either of "mingw64" or "ucrt" and do not mix the two. The following examples make use pf ucrt  
- Tools 
    pacman -S base-devel cmake git mingw-w64-ucrt-x86_64-gcc paraview  
- Libs (Opencl, OpenMP, gsl, fftw, vtk)
    pacman -S mingw-w64-ucrt-x86_64-opencl-headers mingw-w64-ucrt-x86_64-opencl-clhpp mingw-w64-ucrt-x86_64-opencl-icd mingw-w64-ucrt-x86_64-openmp mingw-w64-ucrt-x86_64-gsl mingw-w64-ucrt-x86_64-fftw mingw-w64-ucrt-x86_64-vtk 

- GCC added to PATH in MSYS
    - In the root directory, run `export PATH=$PATH:/ucrt64/bin`
    - this can also be added in windows "edit system environment variables" , "Path". Add C:\msys64\usr\bin and C:\msys64\ucrt64\bin



## Usage

## Todo
- Setup particles for MagLIF
    - Two sets of plasmas
        - 1. Low density plasma cylinder
        - 2. High density thin cylindrical shell around plasma cylinder.
- Setup external fields. 
    - Electric field along the cylinder.
- add in "artificial viscosity" to simulate energy loss/gain
    - F=q(E+vxB)+rv w
    - here r is negative when there is energy loss and positive when there is energy gain

- to get more performance, you might want to recompile the libraries used. for example to install fftw3 recompiled with OMP enabled:
 > wget https://www.fftw.org/fftw-3.3.10.tar.gz

 > tar xvzf fftw-3.3.10.tar.gz

> cd fftw-3.3.10/

> ./configure --enable-threads --enable-openmp --enable-avx --enable-avx2 --enable-avx512 --enable-avx-128-fma --enable-float --with-our-malloc --enable-sse2

> make
 > make install
