# Declaration of variables

CC = g++  -std=c++17 -fconcepts

#CC = mpic++ 
LFLAGS          = $(OPTI) $(WARNINGS) -lgsl -lgslcblas -lboost_system -lboost_filesystem
GIT_COMMIT_VERSION := $(shell git rev-parse HEAD)
GIT_BRANCH	:= $(shell git rev-parse --abbrev-ref HEAD)
DEBUG           = -g -Wall -D DEBUGING=1
OPTI            = -O2 
CFLAGS          = $(OPTI) $(WARNINGS) $(DEBUG) -D__GIT_COMMIT__=\"$(GIT_COMMIT_VERSION)\" -D__GIT_BRANCH__=\"$(GIT_BRANCH)\" -I ./inc
   

# File names
INCDIR=./inc
SRCDIR=./src
OBJDIR=./lib

EXEC = bin/df_program
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

	 
# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(IA_OBJECTS) -o $(EXEC) $(LFLAGS)
	 
# To obtain object files
%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) -c $(CFLAGS) $< -o $@
# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS) $(IA_OBJECTS)

run:
	make;$(EXEC) starting.ini

# Statically analyse source code using 'cppcheck'
.PHONY: cppcheck
cppcheck:
	if type cppcheck ; then cppcheck . ; fi


