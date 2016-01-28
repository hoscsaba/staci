################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
staci_main.o 

CPP_SRCS += \
AnyOption.cpp \
Agelem.cpp \
AnyOption.cpp \
Cso.cpp \
Csomopont.cpp \
Csatorna.cpp \
BukoMutargy.cpp \
JelleggorbesFojtas.cpp \
KonstNyomas.cpp \
Staci.cpp \
Szivattyu.cpp \
Vegakna.cpp \
VisszacsapoSzelep.cpp \
data_io.cpp \
lubksb.cpp \
ludcmp.cpp \
xmlParser.cpp \
staci_main.cpp

OBJS += \
AnyOption.o \
Agelem.o \
AnyOption.o \
Cso.o \
Csomopont.o \
Csatorna.o \
BukoMutargy.o \
JelleggorbesFojtas.o \
KonstNyomas.o \
Staci.o \
Szivattyu.o \
Vegakna.o \
VisszacsapoSzelep.o \
data_io.o \
lubksb.o \
ludcmp.o \
xmlParser.o \
staci_main.o

CPP_DEPS += \
AnyOption.d \
Agelem.d \
AnyOption.d \
Cso.d \
Csomopont.d \
Csatorna.d \
BukoMutargy.d \
JelleggorbesFojtas.d \
KonstNyomas.d \
Staci.d \
Szivattyu.d \
Vegakna.d \
VisszacsapoSzelep.d \
data_io.d \
lubksb.d \
ludcmp.d \
xmlParser.d \
staci_main.d


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file AAAA: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -pedantic -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


