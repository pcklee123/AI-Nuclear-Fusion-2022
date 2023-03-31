_DEPS = traj.h
_OBJ = utils.o TS3.o tnp.o generate.o sel_part_print.o save.o get_densityfields.o cl_code.o changedt.o calcEBV_FFT.o calcU.o calc_trilin_constants.o

IDIR = include

#https://stackoverflow.com/questions/14492436/g-optimization-beyond-o3-ofast
CC=g++
#ucrt64
CFLAGS= -I$(IDIR) -I /ucrt64/include/vtk -L /ucrt64/lib/vtk -fopenmp -fopenmp-simd -Ofast -march=native -malign-double -ftree-parallelize-loops=8 -std=c++2b
#mingw64
#CFLAGS= -I$(IDIR) -I /mingw64/include/vtk -L /ming64/lib/vtk -fopenmp -fopenmp-simd -Ofast -march=native -malign-double -ftree-parallelize-loops=8 -std=c++2b
#CFLAGS= -I$(IDIR) -fopenmp -fopenmp-simd -Ofast -march=native -malign-double -ftree-parallelize-loops=8 -std=c++2b

#LIBS= -lm -lgsl -lOpenCL.dll -lfftw3f -lomp.dll -lfftw3f_omp
LIBS= -lm -lgsl -lOpenCL.dll  -lgomp.dll -lfftw3f_omp -lfftw3f  \
#-lvtkCommonCore.dll  -lvtksys.dll -lvtkIOXML.dll -lvtkCommonDataModel.dll
#-lvtkIOLegacy.dll -lvtkCommonComputationalGeometry.dll -lvtkCommonSystem.dll
#-lvtkGraphics.dll -lvtkFiltersGeneral.dll -lvtkImagingCore.dll -lvtkFiltersGeneric.dll -lvtkIOCore.dll -lvtkIOImage.dll 
AFLAGS= -flto -funroll-loops -fno-signed-zeros -fno-trapping-math -D_GLIBCXX_PARALLEL -fgcse-sm -fgcse-las

#CC=clang++
#CFLAGS=-I$(IDIR) -fopenmp -fopenmp-simd -O3 -Ofast -mavx -mfma -ffast-math -ftree-vectorize -march=native -fomit-frame-pointer -malign-double -std=c++2b
#AFLAGS=
#LIBS=-lm -lgsl -lOpenCL -lfftw3f -lomp.dll

CFLAGS += $(AFLAGS)
CPUS ?= $(shell (nproc --all || sysctl -n hw.ncpu) 2>/dev/null || echo 1)
MAKEFLAGS += --jobs=$(CPUS)

ODIR=obj
LDIR=lib

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TS3: $(OBJ)
	$(CC) -v -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ TS3.exe