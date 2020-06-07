
program_NAME = run_mifas_pp

CXX      = g++

CXXFLAGS = -std=c++11 -O2
LIBS     = -lgsl -lgslcblas -lm
INCLUDE  =


SRC_DIR  = mifas_pp
BIN_DIR  = bin

SOURCES  = $(wildcard $(SRC_DIR)/*.cpp)
OBJS     = $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%.o,$(SOURCES))

 
all: $(program_NAME)


$(program_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS) $(INCLUDE)


$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/*.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<


.PHONY: clean

clean:
	rm -f $(program_NAME) $(OBJS)


