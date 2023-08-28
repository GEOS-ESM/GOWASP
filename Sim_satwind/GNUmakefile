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

_Incdfods     = -I../Lib_ncdfods
_Lncdfods     = -L../Lib_ncdfods -lncdfods

FC	= ifort
FCFLAGS	= -O3 -I. $(_Ishmem) $(_Ibasic) $(_Incdfods) $(_Inetcdf) $(_Ibufr) 

#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Ishmem) $(_Ibasic) $(_Incdfods) $(_Inetcdf) $(_Ibufr) 


mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lshmem) $(_Lbasic) $(_Lncdfods) $(_Lnetcdf) $(_Lbufr)

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
		  m_kx_table.f90	\
		  m_satloc.f90		\
		  m_amv_view.f90	\
		  bufr_amv_write.f90

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_PARAMS	 = compute_params.F90
OBJS_PARAMS	 = $(addsuffix .o, $(basename $(SRCS_PARAMS)))

SRCS_AMVS	 = create_amv.F90
OBJS_AMVS	 = $(addsuffix .o, $(basename $(SRCS_AMVS)))

SRCS_SUM	 = sum_kx_table.f90
OBJS_SUM	 = $(addsuffix .o, $(basename $(SRCS_SUM)))

SRCS_DXMIN	 = compute_dx.f90
OBJS_DXMIN	 = $(addsuffix .o, $(basename $(SRCS_DXMIN)))

SRCS_TXT2B	 = satwind_txt2bufr.f90 
OBJS_TXT2B	 = $(addsuffix .o, $(basename $(SRCS_TXT2B)))

SRCS_ODS2T	 = satwind_ods2txt.f90 
OBJS_ODS2T	 = $(addsuffix .o, $(basename $(SRCS_ODS2T)))

SRCS_CNTAM	 = count_amv.f90 
OBJS_CNTAM	 = $(addsuffix .o, $(basename $(SRCS_CNTAM)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_PARAMS) $(SRCS_AMVS) $(SRCS_SUM) $(SRCS_DXMIN) $(SRCS_TXT2B) $(SRCS_ODS2T) $(SRCS_CNTAM) 

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

OBSFF90 =  m_amv_fields.o 

ALL	= compute_params.x create_amv.x sum_kx_table.x compute_dx.x satwind_txt2bufr.x satwind_ods2txt.x count_amv.x
all: $(ALL)

# make excutable file for sim_rad.x
compute_params.x: $(LIB) $(OBJS_PARAMS) $(OBSFF90) 
	$(mpLD) -o $@ $(OBJS_PARAMS) $(OBSFF90) $(_Lmine) $(LDFLAGS)

create_amv.x: $(LIB) $(OBJS_AMVS) $(OBSFF90) 
	$(mpLD) -o $@ $(OBJS_AMVS) $(OBSFF90) $(_Lmine) $(LDFLAGS)

sum_kx_table.x: $(LIB) $(OBJS_SUM) 
	$(FC) -o $@ $(OBJS_SUM) $(_Lmine) $(LDFLAGS)

compute_dx.x: $(OBJS_DXMIN) 
	$(FC) -O3 -o $@ $(OBJS_DXMIN) 

satwind_txt2bufr.x: $(LIB) $(OBJS_TXT2B) 
	$(FC) -O3 -o $@ $(OBJS_TXT2B) $(LDFLAGS)

satwind_ods2txt.x: $(LIB) $(OBJS_ODS2T) 
	$(mpLD) -O3 -o $@ $(OBJS_ODS2T) $(LDFLAGS)

count_amv.x: $(LIB) $(OBJS_CNTAM) 
	$(mpLD) -O3 -o $@ $(OBJS_CNTAM) $(_Lmine) $(LDFLAGS)

clean:
	rm -f *.o *.mod *.x *.a

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

