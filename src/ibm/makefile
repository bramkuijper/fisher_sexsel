# executable name
EXE=fisherian.exe

# c++ files needed to build
CPP=fisherian.cpp

# which compiler?
CXX=g++

# any flags that the compiler uses?
CXXFLAGS=-Wall -O3

# actually building the executable, this will run
# if you type 'make'
$(EXE) : $(CPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

# type 'make clean' in case you want to remove the exe
clean : 
	rm -rf $(EXE)
