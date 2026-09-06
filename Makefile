OSNF_DIR = osnf

.PHONY: all release debug osnf_code cleanall
CLEANDIRS = $(OSNF_DIR) ./

# LINUX or WIN64
#PLATFORM = WIN64   
PLATFORM = LINUX
DEBUG =#-gstabs -fbounds-check #-pg
MPI    =#-DMPI1
OPT    =-O3  #-fopenmp-simd -fopenmp 

# these three lines should be edited for your system. On systems 
# that do not have separate fortran and c libraries, set NETCDF_FOR and NETCDF_C
# to the same, and set NETCDF_LIB to -lnetcdf (i.e. without the extra f)
#NETCDF_FOR=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.4-mac/
#NETCDF_C=/Users/mccikpc2/Dropbox/programming/netcdf-4.4.1.1-mac/
NETCDF_FOR ?=
NETCDF_C ?=
NETCDF_LIB ?=-lnetcdff 

# The original NETCDF_FOR / NETCDF_C / NETCDFLIB / NETCDFMOD route is kept
# below for backward compatibility.  If none of those are set on Linux,
# nf-config is used automatically instead.
NF_CONFIG ?= nf-config
NETCDF_FFLAGS ?=
NETCDF_FLIBS ?=

#NETCDFLIB=-L/usr/lib/x86_64-linux-gnu/
#NETCDFMOD=/usr/include/

ifeq ($(strip $(PLATFORM)),$(strip LINUX))

# Use the original explicit paths when supplied.
ifneq ($(strip $(NETCDF_FOR)$(NETCDF_C)),)
NETCDFLIB ?= $(if $(strip $(NETCDF_FOR)),-L ${NETCDF_FOR}/lib/,) \
             $(if $(strip $(NETCDF_C)),-L ${NETCDF_C}/lib/,)
NETCDFMOD ?= $(if $(strip $(NETCDF_FOR)),${NETCDF_FOR}/include/,${NETCDF_C}/include/)
endif

# NETCDFLIB and NETCDFMOD can also still be supplied directly, as before.
ifneq ($(strip $(NETCDFMOD)),)
NETCDF_FFLAGS := -I ${NETCDFMOD}
endif
ifneq ($(strip $(NETCDFLIB)),)
NETCDF_FLIBS := ${NETCDFLIB} ${NETCDF_LIB}
endif

# Modern fallback: use nf-config only when no legacy flags were supplied.
ifeq ($(strip $(NETCDF_FFLAGS)),)
NETCDF_FFLAGS := $(shell $(NF_CONFIG) --fflags 2>/dev/null)
endif
ifeq ($(strip $(NETCDF_FLIBS)),)
NETCDF_FLIBS := $(shell $(NF_CONFIG) --flibs 2>/dev/null)
endif

FOR = gfortran -c  #-fno-underscoring 
FOR2 = gfortran  #-fno-underscoring 
AR = ar 
RANLIB = ranlib 
SIXTY_FOUR_F= 

else

# Keep the original Windows defaults.
NETCDFLIB ?= -L../netcdf-4.1.3-mingw/lib/
NETCDFMOD ?= ../netcdf-4.1.3-mingw/include/

ifeq ($(strip $(NETCDF_FFLAGS)),)
NETCDF_FFLAGS := -I ${NETCDFMOD}
endif
ifeq ($(strip $(NETCDF_FLIBS)),)
NETCDF_FLIBS := ${NETCDFLIB} ${NETCDF_LIB}
endif

FOR = x86_64-w64-mingw32-gfortran -c  #-fno-underscoring 
FOR2 = x86_64-w64-mingw32-gfortran  #-fno-underscoring 
AR = x86_64-w64-mingw32-ar
RANLIB = x86_64-w64-mingw32-ranlib
SIXTY_FOUR_F= -static -static-libgfortran #-m32 # note, uncomment the m32 for 32-bit model
endif



OBJ = o

FFLAGS = $(SIXTY_FOUR_F) $(OPT)  $(DEBUG) $(HELIOCENTRIC_FLAG) -o # heliocentric means that we keep the sun at (0,0,0)
FFLAGS2 = $(SIXTY_FOUR_F) $(OPT) $(DEBUG) -o # 
VAR_TYPE = 1 # 0 single, 1 double

# Backward-compatible compile-time frame switch.
# The historical/default behaviour is HELIOCENTRIC=1, which defines
# -Dheliocentric and keeps the Sun fixed at (0,0,0).
HELIOCENTRIC ?= 1
ifeq ($(strip $(HELIOCENTRIC)),1)
HELIOCENTRIC_FLAG = -Dheliocentric
else
HELIOCENTRIC_FLAG =
endif


main.exe	:  odelib.a  solar_system.$(OBJ) 
	$(FOR2) $(FFLAGS2)main.exe solar_system.$(OBJ)   \
	 -lm odelib.a ${NETCDF_FLIBS} $(DEBUG) -L$(OSNF_DIR) 

odelib.a	:   osnf_code
	$(AR) rc odelib.a $(OSNF_DIR)/numerics.$(OBJ) $(OSNF_DIR)/zeroin.$(OBJ) $(OSNF_DIR)/sfmin.$(OBJ) \
				$(OSNF_DIR)/fmin.$(OBJ) $(OSNF_DIR)/r1mach.$(OBJ) \
                $(OSNF_DIR)/d1mach.$(OBJ) $(OSNF_DIR)/dfsid1.$(OBJ) \
                $(OSNF_DIR)/poly_int.$(OBJ) $(OSNF_DIR)/find_pos.$(OBJ) \
                $(OSNF_DIR)/svode.$(OBJ) \
                $(OSNF_DIR)/slinpk.$(OBJ) $(OSNF_DIR)/vode.$(OBJ) \
                $(OSNF_DIR)/dlinpk.$(OBJ) $(OSNF_DIR)/vode_integrate.$(OBJ) \
                $(OSNF_DIR)/erfinv.$(OBJ) $(OSNF_DIR)/tridiagonal.$(OBJ) \
                $(OSNF_DIR)/hygfx.$(OBJ) $(OSNF_DIR)/random.$(OBJ)			

solar_system.$(OBJ)   : solar_system.f90 osnf_code
	$(FOR) -cpp solar_system.f90 -I$(OSNF_DIR) ${NETCDF_FFLAGS} -cpp $(FFLAGS)solar_system.$(OBJ)

osnf_code:
	$(MAKE) -C $(OSNF_DIR)

# Plain `make` retains the original optimised build. These are optional
# convenience targets only.
all: main.exe

release:
	$(MAKE) clean
	$(MAKE) OPT="-O3" DEBUG="" main.exe

debug:
	$(MAKE) clean
	$(MAKE) OPT="-O0" DEBUG="-g -fbounds-check -fbacktrace -Wall -Wextra" main.exe

clean: 
	rm -f *.exe  *.o *.mod *~ \
	odelib.a

cleanall:
	for i in $(CLEANDIRS); do \
		$(MAKE) -C $$i clean; \
	done
