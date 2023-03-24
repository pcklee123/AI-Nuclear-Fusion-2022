## Prerequisites
- MSYS2 Installed
- Add gcc to PATH in MSYS
    - In the root directory, run `export PATH=$PATH:/mingw64/bin`

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