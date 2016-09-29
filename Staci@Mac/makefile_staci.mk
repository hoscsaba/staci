################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include subdir.mk

LDFLAGS:   -L/usr/local/lib -L/usr/lib

all: new_staci

new_staci: $(OBJS) $(USER_OBJS)
	g++ -o new_staci $(OBJS) $(USER_OBJS) -lga -lm -lumfpack

clean:
	-rm $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) new_staci
