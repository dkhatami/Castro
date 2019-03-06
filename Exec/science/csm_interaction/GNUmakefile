PRECISION = DOUBLE
PROFILE   = FALSE

DEBUG	  = FALSE

DIM 	  = 1

COMP	  = gnu

USE_MPI   = TRUE
USE_OMP   = TRUE


USE_RAD   = TRUE
USE_NEUTRINO = FALSE

CASTRO_HOME  = ../../..
Opacity_dir := csm
EOS_dir		 := gamma_law_general

NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(NETWORK_DIR)/gammalaw.net
Bpack := ./Make.package
Blocs := .

include $(CASTRO_HOME)/Exec/Make.Castro
