# This is the makefile for testing the mesa chem package 

#MESA_DIR = ../../..

#################################################################

# STEP 1: set up the standard info

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lchem
LOAD_OTHER = -L$(MESA_LIB_DIR) $(LOAD_MESA_NUMERICS)

#################################################################

# STEP 3: build the test program

TEST_OBJS = test_chem.o
TEST = tester

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)
	
#################################################################

%.o: $(TEST_SRC_DIR)/%.f
	$(TEST_COMPILE) $(FCfixed) $<

clean:
	-@rm -f *.o *.mod $(TEST_LIB) $(TEST_DIR)/$(TEST)
