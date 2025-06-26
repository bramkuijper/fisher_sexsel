# makefile for object-oriented simulation
# of Fisherian sexual selection

# name of executable
EXE=fisher_sexsel.exe

# compiler command
CXX=g++

# compiler flags
CXXFLAGS=-Wall -Wconversion -std=c++20 -O3

# required C++ files
CPP=individual.cpp fisher.cpp main.cpp

# required header files
HPP=individual.hpp fisher.hpp parameters.hpp

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean :
	rm -rf $(EXE)

