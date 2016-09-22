################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include subdir_calibrate.mk

LDFLAGS:   -L/usr/local/lib -L/usr/lib 

all: staci_calibrate

staci_calibrate: $(OBJS) $(USER_OBJS)
	g++ -o staci_calibrate $(OBJS) $(USER_OBJS) -lga -lm -lumfpack

clean:
	-rm $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) staci_calibrate
