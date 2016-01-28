################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# -include ../makefile.init

RM := rm -rf

CFLAGS = -Wall -O0 -pedantic

# All of the sources participating in the build are defined here
-include sources_staci.mk
-include subdir_staci.mk
-include objects_staci.mk

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: staci.exe

# Tool invocations
staci.exe: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ $(CFLAGS) -o"staci.exe" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) staci.exe
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

# -include ../makefile.targets
