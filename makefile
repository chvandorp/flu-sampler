CC = g++ -std=c++0x
CC_FLAGS = -O3 -Wall -c
LD_FLAGS = -O3 -Wall
LIBS = -lgsl -lgslcblas -lm -lpthread
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

all: flu-sampler
 
clean:
	rm -f obj/*.o obj/*.d flu-sampler
 
flu-sampler: $(OBJ_FILES)
	$(CC) $(LD_FLAGS) $^ $(LIBS) -o $@
 
obj/%.o: src/%.cpp $^
	$(CC) $(CC_FLAGS) $< -o $@

$(OBJ_FILES): | obj

obj:
	mkdir -p obj
 
CC_FLAGS += -MMD
-include $(OBJ_FILES:.o=.d)
