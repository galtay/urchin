FC = gfortran-4.4
CC = gcc 
FFLAGS = -g -fPIC
CFLAGS = -g -fPIC
LINKFLAGS = -I../../../hdf5/hdf5-1.8.6-gfortran-4.4/hdf5/lib -I../../../hdf5/hdf5-1.8.6-gfortran-4.4/hdf5/include/ -I/home/galtay/google-svn/sea-urchin/hdf5_wrapper/hdf5-1.8.6-gfortran-4.4/src/../lib -L/home/galtay/google-svn/sea-urchin/hdf5_wrapper/hdf5-1.8.6-gfortran-4.4/src/../lib -Wl,-rpath,/home/galtay/google-svn/sea-urchin/hdf5_wrapper/hdf5-1.8.6-gfortran-4.4/src/../lib -I/home/galtay/google-svn/sea-urchin/hdf5_wrapper/hdf5-1.8.6-gfortran-4.4/src/../lib -lhdfwrapper

.SUFFIXES:
.SUFFIXES: .f .c .o .f90

ALLFFLAGS = $(FFLAGS) $(LINKFLAGS)
ALLCFLAGS = $(CFLAGS)

.f.o:
	$(FC) $(ALLFFLAGS) -c $< -o $*.o
.f90.o:
	$(FC) $(ALLFFLAGS) -c $< -o $*.o
.c.o:
	$(CC) $(ALLCFLAGS) -c $< -o $*.o

all: test_write.exe test_read.exe test_extend.exe

test_write.exe: test_write.o
	$(FC)   test_write.o $(ALLFFLAGS) -o test_write.exe

test_read.exe: test_read.o
	$(FC)  test_read.o $(ALLFFLAGS) -o test_read.exe

test_extend.exe: test_extend.o
	$(FC)  test_extend.o $(ALLFFLAGS) -o test_extend.exe

clean:
	rm -f *.o *.mod
	rm -f test_read.exe test_write.exe test_extend.exe


