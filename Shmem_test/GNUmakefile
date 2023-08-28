#
# Makefile for building the SHMEM reader example.

BASELIBDIR = ${BASEDIR1}/Linux


NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

_Ishmem	      = -I../Lib_shmem1
_Lshmem	      = -L../Lib_shmem1 -lshmem1


.SUFFIXES: .o .F90 .x

all: shmem_reader.x

FC = mpiifort
FOPTS = -I. $(_Ishmem) $(_Inetcdf)

shmem_reader.x: shmem_reader.o
	$(FC) $(FOPTS) $^ -o $@ $(_Lnetcdf) $(_Lshmem)

clean:
	rm -f *.o *.mod core.* *~ 

distclean: clean
	rm -f *.x

.F90.o:
	$(FC) -c $(FOPTS) $< 
