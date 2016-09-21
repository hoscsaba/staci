################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

include /Users/hoscsaba/Utils/SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk

LDFLAGS:   -L/usr/local/lib -L/Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Lib	
CPPFLAGS:  -I/usr/local/include
LIBS = $(LAPACK) $(BLAS) $(XERBLA) $(LIB) $(LIB_WITH_CHOLMOD) $(LIB_WITH_PARTITION) $(CUBLAS_LIB) $(CUDART_LIB)
UMFPACK = /Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Lib/libumfpack.a /Users/hoscsaba/Utils/SuiteSparse/AMD/Lib/libamd.a \
    /Users/hoscsaba/Utils/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a 

INC:-I/Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Include -I/Users/hoscsaba/Utils/SuiteSparse/AMD/Include -I/Users/hoscsaba/Utils/SuiteSparse/SuiteSparse_config

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: new_staci staci_calibrate

staci_calibrate: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -I/Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Include\
	    -I/Users/hoscsaba/Utils/SuiteSparse/AMD/Include\
	    -I/Users/hoscsaba/Utils/SuiteSparse/SuiteSparse_config\
	    -o "staci_calibrate" $(OBJS) $(USER_OBJS) $(LIBS) $(UMFPACK)\
	    /Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Lib/libumfpack.a\
		/Users/hoscsaba/Utils/SuiteSparse/AMD/Lib/libamd.a\
		/Users/hoscsaba/Utils/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a\
		/Users/hoscsaba/Utils/SuiteSparse/CHOLMOD/Lib/libcholmod.a\
		/Users/hoscsaba/Utils/SuiteSparse/COLAMD/Lib/libcolamd.a\
	    -L/usr/local/lib\
	    -L/Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Lib\
	    $(LDFLAGS) 
	    -Lga -lga -lm 
	@echo 'Finished building target: $@'
	@echo ' '

# Tool invocations
new_staci: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -I/Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Include\
	    -I/Users/hoscsaba/Utils/SuiteSparse/AMD/Include\
	    -I/Users/hoscsaba/Utils/SuiteSparse/SuiteSparse_config\
	    -o "new_staci" $(OBJS) $(USER_OBJS) $(LIBS) $(UMFPACK)\
	    /Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Lib/libumfpack.a\
		/Users/hoscsaba/Utils/SuiteSparse/AMD/Lib/libamd.a\
		/Users/hoscsaba/Utils/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a\
		/Users/hoscsaba/Utils/SuiteSparse/CHOLMOD/Lib/libcholmod.a\
		/Users/hoscsaba/Utils/SuiteSparse/COLAMD/Lib/libcolamd.a\
	    -L/usr/local/lib\
	    -L/Users/hoscsaba/Utils/SuiteSparse/UMFPACK/Lib\
	    $(LDFLAGS) 
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) new_staci
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
