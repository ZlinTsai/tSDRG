CXX           := g++
CXXFLAGS      := -m64 -O3 -std=c++11
UNI10CXXFLAGS := -DUNI_CPU -DUNI_LAPACK
UNI10_ROOT    := /usr/local/uni10
ROOTS         := /usr/local/uni10
MAINSCRIPT    := $(code)
SRCDIRS       := ../Operator ../Hamiltonian ../MPO ../tSDRG_tools
SRC           := operator.cpp hamiltonian.cpp mpo.cpp tSDRG_tools.cpp measure.cpp

# Check variable
ifeq ($(wildcard $(code)),)
    $(error Enter code=XXX.cpp)
endif

ifeq ($(wildcard $(name)),)
    name := job.exe
endif

# The name of binary object.
EXU := $(MAINSCRIPT)

# Search the source files in under paths.
vpath %.cpp  $(SRCDIRS)

# The file names of source codes.
src := $(SRC)

COMMON_FLAGS_INCS += $(foreach root,$(ROOTS),-I$(root)/include)
COMMON_FLAGS_LIBS += $(foreach root,$(ROOTS),-L$(root)/lib)
CXXFLAGS += $(UNI10CXXFLAGS)
UNI10LDLIBRARY := -luni10

OBJ_DIR := objects
OBJECTS := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(src))

all: clean job.exe

job.exe: $(EXU) | $(OBJECTS)
	$(CXX) $(COMMON_FLAGS_INCS) $(COMMON_FLAGS_LIBS) $(CXXFLAGS) $^ -o $(name) $(OBJECTS) $(UNI10LDLIBRARY) 

$(OBJECTS): | $(OBJ_DIR)

$(OBJ_DIR):
	@mkdir -p $@

$(OBJ_DIR)/%.o: %.cpp
	$(CXX) $(COMMON_FLAGS_INCS) $(CXXFLAGS) -c $< -o $@

Clean:
	rm -rf job.exe $(OBJ_DIR)

clean:
	touch job.exe
	rm job.exe
