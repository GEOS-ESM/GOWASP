#!/usr/bin/make 

BASELIBDIR = ${BASEDIR1}/Linux
GEOSPATH = ${GEOSPATH1}

NETCDF_LIBDIR = ${BASELIBDIR}/lib
NETCDF_INCDIR = ${BASELIBDIR}/include/netcdf

_Inetcdf      = -I$(NETCDF_INCDIR) 
_Lnetcdf      = $(shell ${BASELIBDIR}/bin/nf-config --flibs)

_Ibufr	      = -I${GEOSPATH}/Linux/include/NCEP_bufr_r4i4 
_Lbufr	      = -L${GEOSPATH}/Linux/lib -lNCEP_bufr_r4i4

_Icrtm	      = -I${GEOSPATH}/src/NCEP_Shared/NCEP_crtm
_Lcrtm	      = -L${GEOSPATH}/Linux/lib -lNCEP_crtm

_Ibasic	      = -I../Lib_basic1
_Lbasic	      = -L../Lib_basic1 -lbasic1

_Ishmem	      = -I../Lib_shmem1
_Lshmem	      = -L../Lib_shmem1 -lshmem1

_Iobsrw	      = -I../Lib_obsrw1
_Lobsrw	      = -L../Lib_obsrw1 -lobsrw1

_Incdfods     = -I../Lib_ncdfods
_Lncdfods     = -L../Lib_ncdfods -lncdfods

FC	= ifort

FCFLAGS	= -O3 -I. $(_Iobsrw) $(_Ibasic) $(_Ibufr) $(_Incdfods) $(_Inetcdf) $(_Icrtm) $(_Ishmem) 
#FCFLAGS	= -O3 -I. $(_Ibasic) $(_Iobsrw) $(_Ibufr) $(_Icrtm) $(_Incdfods) $(_Inetcdf) $(_Ishmem) 


#FCFLAGS = -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0  -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -ftrapuv -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Incdfods) $(_Inetcdf) $(_Ibufr) $(_Icrtm) 

#FCFLAGS = -g -O0 -check all -traceback -I. $(_Iobsrw) $(_Ishmem) $(_Ibasic) $(_Incdfods) $(_Inetcdf) $(_Ibufr) $(_Icrtm) 

mpLD	= mpiifort
LD	= $(FC)
RM	= rm -f
AR	= ar cqs

LDFLAGS		= $(_Lobsrw) $(_Lbasic) $(_Lnetcdf) $(_Lbufr) $(_Lcrtm) 
LDFLAGS1	= $(_Lobsrw) $(_Lbasic) $(_Lnetcdf) $(_Lbufr) 
LDFLAGS2	= $(_Lobsrw) $(_Lshmem) $(_Lbasic) $(_Lnetcdf) 
LDFLAGS3	= $(_Lobsrw) $(_Lshmem) $(_Lbasic) $(_Lbufr) $(_Lcrtm) 
LDFLAGS4	= $(_Lobsrw) $(_Lbasic)
LDFLAGS5	= $(_Lbasic) $(_Lncdfods) $(_Lnetcdf)

LIB		= libmine.a
_Lmine		= -L. -lmine
_I		= -I.

SRCS_mylib	= \
		  compute_plevs.f90	\
		  prof_io_files.f90	\
                  m_copy_rad_obs.f90	\
                  m_crtm_interface.f90	\
		  m_obs_list.f90	\
		  m_prof_merge.f90	\
                  m_rad_index.f90	\
                  m_rad_prob.f90	\
		  m_rad_thin.f90	\
		  m_rad_thin_flds.f90	\
                  m_read_profiles.f90	\
		  m_set_unit_nums.f90	\
                  m_write_profiles.f90	\
                  m_realloc_bad.f90	\
                  m_realloc_get.f90	\
                  sub_prof_all_2d3d.f90 \
                  read_nc4_file.F90

OBJS_mylib	= $(addsuffix .o, $(basename $(SRCS_mylib)))

SRCS_RAD_BUFR	= create_rad_bufr.F90
SRCS_RAD_LIST	= create_rad_obs_list.f90
SRCS_RAD_PROFS	= create_rad_profs.F90
SRCS_RAD_REORD	= create_rad_reord.f90
SRCS_RAD_MERGE	= create_rad_merge.f90
SRCS_RAD_ODS	= rad_ods2txt.f90

OBJS_RAD_BUFR	= $(addsuffix .o, $(basename $(SRCS_RAD_BUFR)))
OBJS_RAD_LIST	= $(addsuffix .o, $(basename $(SRCS_RAD_LIST)))
OBJS_RAD_PROFS	= $(addsuffix .o, $(basename $(SRCS_RAD_PROFS)))
OBJS_RAD_REORD	= $(addsuffix .o, $(basename $(SRCS_RAD_REORD)))
OBJS_RAD_MERGE	= $(addsuffix .o, $(basename $(SRCS_RAD_MERGE)))
OBJS_RAD_ODS	= $(addsuffix .o, $(basename $(SRCS_RAD_ODS)))

SRCS_ALL	= $(SRCS_mylib) $(SRCS_RAD_PROFS) $(SRCS_RAD_BUFR) $(SRCS_RAD_LIST) $(SRCS_RAD_REORD) $(SRCS_RAD_MERGE) $(SRCS_RAD_ODS)

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


ALL	= create_rad_bufr.x create_rad_obs_list.x create_rad_profs.x create_rad_reord.x create_rad_merge.x rad_ods2txt.x

all:	$(ALL)

# make excutable file for sim_rad.x
create_rad_obs_list.x: $(LIB) $(OBJS_RAD_LIST) 
	$(mpLD)   -o $@  $(OBJS_RAD_LIST) $(_Lmine) $(LDFLAGS1)

create_rad_profs.x: $(LIB) $(OBJS_RAD_PROFS) 
	$(mpLD) -o $@ $(OBJS_RAD_PROFS) $(_Lmine) $(LDFLAGS2)  

create_rad_bufr.x: $(LIB) goss_mpeu.o m_mympi.o $(OBJS_RAD_BUFR) 
	$(mpLD) -o $@ $(OBJS_RAD_BUFR) goss_mpeu.o m_mympi.o $(_Lmine) $(LDFLAGS3)

create_rad_reord.x: $(LIB) $(OBJS_RAD_REORD) 
	$(LD) -o $@ $(OBJS_RAD_REORD) $(_Lmine) $(LDFLAGS4)

create_rad_merge.x: $(LIB) $(OBJS_RAD_MERGE) 
	$(LD) -o $@ $(OBJS_RAD_MERGE) $(_Lmine) $(LDFLAGS4)

rad_ods2txt.x: $(LIB) $(OBJS_RAD_ODS) 
	$(mpLD) -O3 -o $@ $(OBJS_RAD_ODS) $(LDFLAGS5)

clean:
	rm -f *.o *.mod *.x *.a

lib:	$(LIB)

$(LIB): $(OBJS_mylib)
	$(RM) $@
	$(AR) $@ $(OBJS_mylib)

include make.prereq

new_prereq: $(SRCS_ALL)
	@ for f in $(SRCS_ALL); do fdp $$f; done

