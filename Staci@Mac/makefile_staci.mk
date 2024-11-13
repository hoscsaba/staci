-include subdir.mk

LINK=   -lm -lumfpack
INC+=	-I/Users/hoscsaba/Documents/GitHub/SuiteSparse/UMFPACK/include\
-I/Users/hoscsaba/Documents/GitHub/SuiteSparse/SuiteSparse_config\
-I/Users/hoscsaba/Documents/GitHub/SuiteSparse/AMD/Include\
-L/Users/hoscsaba/Documents/GitHub/SuiteSparse/UMFPACK/build


all: new_staci

new_staci: $(OBJS) $(USER_OBJS)
	g++ -rpath /Users/hoscsaba/Documents/GitHub/SuiteSparse/UMFPACK/build $(INC) $(LINK) $(OBJS) $(USER_OBJS)  -o new_staci 


clean:
	-rm $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) new_staci
	