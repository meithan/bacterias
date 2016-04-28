# Compiler
COMPILER= g++

# User compiler flags
#USER_FLAGS = -g -Wall -pedantic
USER_FLAGS= -O3

# ==============================================================================

LDFLAGS= 
CCFLAGS= $(USER_FLAGS)
PROGRAMS= bacterias.exe

# ==============================================================================

bacterias : Foodsource.h Strain.h Bacterium.h bacterias.o 
	$(COMPILER) $(LDFLAGS) bacterias.o -o bacterias.exe

%.o : %.cpp
	$(COMPILER) $(CCFLAGS) -c -o $@ $<

%.h : ;

clean :
	rm -f *.o $(PROGRAMS)
