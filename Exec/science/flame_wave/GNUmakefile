PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = Intel
FCOMP	   = Intel

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE

USE_MODELPARSER  = TRUE

USE_ROTATION  = TRUE
USE_DIFFUSION = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../..

ifdef MICROPHYSICS_HOME

  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_dir     := helmholtz

  # This sets the network directory in $(MICROPHYSICS_HOME)/networks
  Network_dir := triple_alpha_plus_cago

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Conductivity_dir := constant_opacity

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro