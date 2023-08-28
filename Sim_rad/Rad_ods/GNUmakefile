

#!/usr/bin/make
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 02Nov2011  Jing Guo	Created out of Ron Errico's GNUmakefile template.
#



# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
#  ESMADIR = ${NOBACKUP}/../nprive/GEOSadas
#  ESMADIR = /discover/nobackup/mkim1/progress_cvs/EnADAS-5_15_2p1/GEOSadas
#  ESMADIR = /discover/nobackup/rerrico/Tropics_Build/GEOSadas-5_14.OSSE.genrad
   ESMADIR = ${GOWASP_GSIPATH}
endif


# Make sure BASEDIR is defined
ifndef BASEDIR
#  BASEDIR = /discover/nobackup/projects/gmao/share/dao_ops/Baselibs/v3.2.0_build3
   BASEDIR = $basedir
#  BASEDIR = /discover/nobackup/projects/gmao/share/dao_ops/Baselibs/v4.0.2_build1
endif


default: esma_help

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS  := $(shell basename `pwd`)
LIB   :=
ifneq ( $(LIBSRCS),$(null) )
  LIB   := lib$(THIS).a
endif

esma_install install: all

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] $(APPS)

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] $(APPS)

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make [esma_]install    (builds and install under ESMADIR)"
	@echo "% make [esma_]clean      (removes deliverables: *.[aox], etc)"
	@echo "% make [esma_]distclean  (leaves in the same state as cvs co)"
	@echo "% make [esma_]help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"
	@echo "make's variables"
	@echo "         APPS = $(APPS)"
	@echo "      APPSRCS = $(APPSRCS)"
	@echo "      APPOBJS = $(APPOBJS)"
	@echo "          LIB = $(LIB)"
	@echo "      LIBSRCS = $(LIBSRCS)"
	@echo "      LIBOBJS = $(LIBOBJS)"
	@echo "         DEPS = $(DEPS)"
	@echo ""

more_details:
	@echo "make's variables"
	@echo "      INC_ODS = $(INC_ODS)"
	@echo "      LIB_ODS = $(LIB_ODS)"
	@echo "     INC_MPEU = $(INC_MPEU)"
	@echo "     LIB_MPEU = $(LIB_MPEU)"
	@echo "      INC_SDF = $(INC_SDF)"
	@echo "      LIB_SDF = $(LIB_SDF)"
	@echo "      LIB_SYS = $(LIB_SYS)"

INC_ODS  = $(ESMAINC)/GMAO_ods
INC_MPEU = $(ESMAINC)/GMAO_mpeu

LIB_ODS  = $(ESMALIB)/libGMAO_ods.a
LIB_MPEU = $(ESMALIB)/libGMAO_eu.a

MOD_DIRS  = . $(INC_ODS) $(INC_MPEU) $(INC_SDF) 

vpath % $(MOD_DIRS)

##--- user configuration starts
APPSRCS := rad_ods2txt.f90      



LIBSRCS := m_sat_info_table.f90
##--- user configuration ends

APPSRCS	:= $(wildcard $(APPSRCS))
APPOBJS	:= $(addsuffix .o, $(basename $(APPSRCS)))
APPS	:= $(addsuffix .x, $(basename $(APPSRCS)))

LIBOBJS := $(addsuffix .o, $(basename $(LIBSRCS)))

DEPS	:= $(addsuffix .d, $(basename $(LIBSRCS)))
DEPS	+= $(addsuffix .d, $(basename $(APPSRCS)))

FOPT = $(FOPT3)
FREAL = $(FREAL4)

MMOD_DIRS = $(M). $(M)$(INC_ODS) $(M)$(INC_MPEU) $(M)$(INC_SDF) 

USER_FFLAGS = $(BIG_ENDIAN) $(MMOD_DIRS) 
USER_LIBS   = $(LIB_ODS) $(LIB_MPEU) $(LIB_SDF) $(LIB_SYS) $(LIBOBJS) 

lib:	$(LIB)
ifneq ( $(LIB),$(null) )
$(LIB): $(LIBOBJS)
	$(AR) $(AR_FLAGS) $@ $(LIBOBJS)
endif

all : $(APPS)

%.x: $(LIB) %.o
	$(LD)  -$(LDFLAGS)  -o $@ $*.o $(LIB) $(USER_LIBS)

#                  --------------------
#                      Dependencies
#                  --------------------

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean, $(MAKECMDGOALS)),clean)
  ifneq ($(findstring  help, $(MAKECMDGOALS)), help)
    -include $(DEPS)
  endif
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
