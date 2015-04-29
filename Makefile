
NETCDF = /usr

CXX = g++

OPT_FLAGS = -O3
DEBUG_FLAGS = -g
CXX_FLAGS = -std=c++11

CXX_INCL = -I$(NETCDF)/include
CXX_LIBS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf

OBJ = netcdf_utils.o extract_max_level_cell.o

BIN = extract_max_level_cell

DEPS_HPP = netcdf_utils.h

all: CXX_FLAGS+=$(OPT_FLAGS)
all: $(BIN)

debug: CXX_FLAGS+=$(DEBUG_FLAGS)
debug: $(BIN)

clean:
	rm -f $(OBJ) $(BIN)


$(BIN): $(OBJ)
	$(CXX) -o $@ $^ $(CXX_FLAGS) $(CXX_LIBS)

%.o: %.cpp $(DEPS_HPP)
	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
