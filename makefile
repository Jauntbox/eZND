#################################################################

# STEP 1: get the mesa makefile_header
# define MESA_DIR to be the path to the mesa directory

#MESA_DIR = /Users/Kevin/mesa_5271
MESA_DIR = /Users/Kevin/mesa_5819

include $(MESA_DIR)/utils/makefile_header


#LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -leos
#LOAD_OTHER = -L$(MESA_LIB_DIR) -leos -lchem -lkap -lnum $(LOAD_MATRIX) \
#   -lalert -linterp_1d -linterp_2d -lutils -lconst -lnet -lrates -lweak \
#   -lreaclib -lscreen

LOAD_OTHER = -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO)

#################################################################

# STEP 2: build the program

#PROG = run_profiler
#PROG_OBJS = wd_solver.o lane_emden_solver.o znd_profiler.o

PROG = run_new
PROG_OBJS = wd_solver.o lane_emden_solver.o znd.o

#PROG = run_le
#PROG_OBJS = lane_emden_solver.o lane_emden_driver.o

#PROG = run_wd
#PROG_OBJS = wd_solver.o wd_driver.o					

#PROG = test
#PROG_OBJS = array_test.o

PROG_DIR = .

$(PROG) : $(PROG_OBJS)
	$(FC) $(FCbasic) $(FCopenmp) -o$(PROG_DIR)/$(PROG) $(PROG_OBJS) $(LOAD_OTHER)

#################################################################

MY_FC_FLAGS = $(FCfree)
SRC_DIR = .

%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCbasic) $(MY_FC_FLAGS) -I$(MESA_INCLUDE_DIR) -c $<
