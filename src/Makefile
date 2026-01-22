# Compiler
CXX := g++

# Compiler flags
CXXFLAGS := -std=c++17 -Wall -Wextra

# Target executable name
TARGET := a.out

# Source files
SRC := $(wildcard *.cpp)

# Default target
all: $(TARGET)

# Link / build
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

# Clean
clean:
	rm -f $(TARGET)

.PHONY: all clean
