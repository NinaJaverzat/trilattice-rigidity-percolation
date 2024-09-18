# Compiler
CXX = g++-10

# Compiler flags
CXXFLAGS = -std=c++20 -Wall -Wextra -pedantic -O3

# List of source files
SRCS = main.cpp basic_functions.cpp rigidity_percolation.cpp pebble_game.cpp init.cpp plotting.cpp

# List of object files (derived from the source files)
OBJS = $(SRCS:.cpp=.o)

# Name of the executable
TARGET = rp.exe

# Default rule to build the program
all: $(TARGET)

# Rule to link the object files into the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile the source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to clean up the object files and the executable
clean:
	rm -f $(OBJS) $(TARGET)

clear:
	rm -f $(OBJS)

# Phony targets to avoid conflicts with files named 'all' or 'clean' or 'clean'
.PHONY: all clean clear
