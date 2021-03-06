
#-include subdir_calibrate.mk

# Add inputs and outputs from these tool invocations to the build variables
O_SRCS += \
../Agelem.o \
../AnyOption.o \
../BukoMutargy.o \
../Csatorna.o \
../Cso.o \
../Csomopont.o \
../JelleggorbesFojtas.o \
../KonstNyomas.o \
../Staci.o \
../Staci_Runner.o \
../Szivattyu.o \
../Vegakna.o \
../VisszacsapoSzelep.o \
../bjg_tools.o \
../bjga_main.o \
../bjga_micro_ga.o \
../bjga_vizmu.o \
../data_io.o \
../lubksb.o \
../ludcmp.o \
../staci_split.o \
../xmlParser.o

CPP_SRCS += \
../Agelem.cpp \
../AnyOption.cpp \
../BukoMutargy.cpp \
../Csatorna.cpp \
../Cso.cpp \
../Csomopont.cpp \
../JelleggorbesFojtas.cpp \
../KonstNyomas.cpp \
../Staci.cpp \
../Szivattyu.cpp \
../Vegakna.cpp \
../VisszacsapoSzelep.cpp \
../data_io.cpp \
../lubksb.cpp \
../ludcmp.cpp \
../staci_split.cpp \
../xmlParser.cpp

OBJS += \
./Agelem.o \
./AnyOption.o \
./BukoMutargy.o \
./Csatorna.o \
./Cso.o \
./Csomopont.o \
./JelleggorbesFojtas.o \
./KonstNyomas.o \
./Staci.o \
./Szivattyu.o \
./Vegakna.o \
./VisszacsapoSzelep.o \
./data_io.o \
./lubksb.o \
./ludcmp.o \
./staci_split.o \
./xmlParser.o

CPP_DEPS += \
./Agelem.d \
./AnyOption.d \
./BukoMutargy.d \
./Csatorna.d \
./Cso.d \
./Csomopont.d \
./JelleggorbesFojtas.d \
./KonstNyomas.d \
./Staci.d \
./Szivattyu.d \
./Vegakna.d \
./VisszacsapoSzelep.d \
./data_io.d \
./lubksb.d \
./ludcmp.d \
./staci_split.d \
./xmlParser.d


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include -I/usr/local/include/eigen3 -O2 -g -Wall -c -fmessage-length=0 -pedantic -ansi -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@"  "$<"
	@echo 'Finished building: $<'
	@echo ' '


all: staci_split

staci_split: $(OBJS) $(USER_OBJS)
	g++ -o staci_split $(OBJS) $(USER_OBJS) -L/usr/local/lib -L/usr/lib -L/usr/local/lib -ligraph -lga -lm -lumfpack

clean:
	-rm $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) staci_split
