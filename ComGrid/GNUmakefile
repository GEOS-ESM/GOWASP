#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

FC	= ifort
FCFLAGS	= -O3 -I. $(_Inetcdf) 

mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lnetcdf) 

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
		  read_nc4_file.f90	\
		  set_file_name.f90	\
		  m_sft.90      	\
		  m_write_nc4_file.f90 	\
		  m_write_ana_file.f90 	

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_CGRID	= common_grid.f90
SRCS_ANA	= common_grid_ana.f90
SRCS_BKG	= bkg_fields.f90
OBJS_CGRID	= $(addsuffix .o, $(basename $(SRCS_CGRID)))
OBJS_ANA	= $(addsuffix .o, $(basename $(SRCS_ANA)))
OBJS_BKG	= $(addsuffix .o, $(basename $(SRCS_BKG)))


SRCS_ALL	= $(SRCS_mylib) $(SRCS_CGRID) $(SRCS_ANA) $(SRCS_BKG) 

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


ALL	= common_grid.x common_grid_ana.x bkg_fields.x

all:	$(ALL)

# make excutable file for common_grid.x
common_grid.x: $(LIB) $(OBJS_CGRID) 
	$(mpLD)   -o $@ $(OBJS_CGRID) $(_Lmine) $(LDFLAGS) 

# make excutable file for common_grid_ana.x
common_grid_ana.x: $(LIB) $(OBJS_ANA) 
	$(mpLD)   -o $@ $(OBJS_ANA) $(_Lmine) $(LDFLAGS) 

# make excutable file for bkg_fields.x
bkg_fields.x: $(LIB) $(OBJS_BKG) 
	$(mpLD)   -o $@ $(OBJS_BKG) $(_Lmine) $(LDFLAGS) 

clean:
	rm -f *.o *.mod *.x

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

