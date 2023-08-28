#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

_Ibufr	      = -I${GEOSPATH}/Linux/include/NCEP_bufr_r4i4 
_Lbufr	      = -L${GEOSPATH}/Linux/lib -lNCEP_bufr_r4i4

_Ibasic	      = -I../Lib_basic1
_Lbasic	      = -L../Lib_basic1 -lbasic1

_Ishmem	      = -I../Lib_shmem1
_Lshmem	      = -L../Lib_shmem1 -lshmem1

FC	= ifort
FCFLAGS	= -O3 -I. $(_Ishmem) $(_Ibasic) $(_Inetcdf) $(_Ibufr)
mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lshmem) $(_Lbasic) $(_Lnetcdf) $(_Lbufr) 

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
                  write_raobs.f90

                  

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_GP		 = create_gp_obs.F90
OBJS_GP		 = $(addsuffix .o, $(basename $(SRCS_GP)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_GP) 

.SUFFIXES:
.SUFFIXES: .f90 .F90 .F .f .o

.f.o:
	$(FC) -c $(FCFLAGS) $<

.F.o:
	$(FC) -c $(FCFLAGS) $<

.f90.o:
	$(FC) -c $(FCFLAGS) $<

.F90.o:
	$(mpLD) -c $(FCFLAGS) $<


ALL	= create_gp_obs.x 
all: $(ALL)

# make excutable file for sim_rad.x
create_gp_obs.x: $(LIB) $(OBJS_GP) 
	$(mpLD) -o $@ $(OBJS_GP) $(_Lmine) $(LDFLAGS)

clean:
	rm -f *.o *.mod *.x *.a
lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

