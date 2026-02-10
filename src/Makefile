# Compilers
CXX := g++
FXX := gfortran

# Flags
CXXFLAGS := -O3 -std=c++17 
FXXFLAGS := -O3 -fcheck=all -fbacktrace -I ../inc

# Target
TARGET := a.out

# Sources
CSRC := $(wildcard *.cpp)
FSRC := $(wildcard *.F)

# Object files
COBJ := $(CSRC:.cpp=.o)
FOBJ := $(FSRC:.F=.o)

# Default target
all: $(TARGET)

# Libraries
LIBS= /home/joesteer/lib/libocean.a /home/joesteer/slime/src/liblasta2.a -llapack -lblas

# Link step - use g++ as the linker and add Fortran runtime libraries
# Order: C++ objects, Fortran objects, external libs, Fortran runtime libs
$(TARGET): $(COBJ) $(FOBJ)
	$(CXX) $(COBJ) $(FOBJ) $(LIBS) -lgfortran -lquadmath -lm -o $(TARGET)

# Compile C++ sources
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile Fortran sources
%.o: %.F
	$(FXX) $(FXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(TARGET) *.o

.PHONY: all clean
