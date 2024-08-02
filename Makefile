SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio2211
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio2211/concert
CPLEXLIBDIR   = $(CPLEXDIR)/cplex/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CLNFLAGS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) -lilocplex -lconcert -lcplex -m64 -lm -lpthread

CC_DEBUG = -DDEBUG -g
CC_RELEASE = -DNDEBUG
CC_VALGRIND = -DNDEBUG -g -O0

COPT  = -m64 -O2 -fPIC -fexceptions $(CC_RELEASE) -DIL_STD -DLONG_MAX=0x7FFFFFFFL
GENERALINCDIR   = -I ./codes
CPLEXINCDIR   = -I $(CPLEXDIR)/cplex/include -I $(CPLEXDIR)/concert/include
CFLAGS = $(COPT) $(GENERALINCDIR) -std=c++17
CFLAGS2  = $(COPT) $(CPLEXINCDIR)

CC=ccache g++ -std=c++17

PROG_DIR=codes/src
PROG_BIN=codes/bin

MAIN_SRC=$(PROG_DIR)/main2.cpp

GENERAL_SRC=$(PROG_DIR)/general.cpp
GENERAL_H=$(PROG_DIR)/general.h
GENERAL_OBJ=$(PROG_BIN)/general.o

TIMER_SRC=$(PROG_DIR)/timer.cpp
TIMER_H=$(PROG_DIR)/timer.h
TIMER_OBJ=$(PROG_BIN)/timer.o

PROBLEM_SRC=$(PROG_DIR)/problem.cpp
PROBLEM_H=$(PROG_DIR)/problem.h
PROBLEM_OBJ=$(PROG_BIN)/problem.o

FEASIBILITY_PUMP_SRC=$(PROG_DIR)/feasibility_pump/feasibility_pump.cpp
FEASIBILITY_PUMP_H=$(PROG_DIR)/feasibility_pump/feasibility_pump.h
FEASIBILITY_PUMP_OBJ=$(PROG_BIN)/feasibility_pump.o

KERNEL_PUMP_SRC=$(PROG_DIR)/kernel_pump/kernel_pump.cpp
KERNEL_PUMP_H=$(PROG_DIR)/kernel_pump/kernel_pump.h
KERNEL_PUMP_OBJ=$(PROG_BIN)/kernel_pump.o

HEURISTIC_SOLUTION_SRC=$(PROG_DIR)/heuristic_solution.cpp
HEURISTIC_SOLUTION_H=$(PROG_DIR)/heuristic_solution.h
HEURISTIC_SOLUTION_OBJ=$(PROG_BIN)/heuristic_solution.o

SOLUTION_HPP=$(PROG_DIR)/solution.hpp
SOLUTION_OBJ=$(PROG_BIN)/solution.o

MATRIX_HPP=$(PROG_DIR)/matrix.hpp
MATRIX_OBJ=$(PROG_BIN)/matrix.o

general: $(GENERAL_SRC) $(GENERAL_H)
	$(CC) $(CFLAGS) -c $(GENERAL_SRC) -o $(GENERAL_OBJ)

timer: $(TIMER_SRC) $(TIMER_H)
	$(CC) $(CFLAGS) -c $(TIMER_SRC) -o $(TIMER_OBJ)

problem: $(PROBLEM_SRC) $(PROBLEM_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(PROBLEM_SRC) -o $(PROBLEM_OBJ)

feasibility_pump: $(FEASIBILITY_PUMP_SRC) $(FEASIBILITY_PUMP_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(FEASIBILITY_PUMP_SRC) -o $(FEASIBILITY_PUMP_OBJ)

kernel_pump: $(KERNEL_PUMP_SRC) $(KERNEL_PUMP_H)
	$(CC) $(CFLAGS) $(CFLAGS2) -c $(KERNEL_PUMP_SRC) -o $(KERNEL_PUMP_OBJ)

heuristic_solution: $(HEURISTIC_SOLUTION_SRC) $(HEURISTIC_SOLUTION_H)
	$(CC) $(CFLAGS) -c $(HEURISTIC_SOLUTION_SRC) -o $(HEURISTIC_SOLUTION_OBJ)

kp: general timer problem heuristic_solution feasibility_pump kernel_pump
	$(CC) $(CFLAGS) $(CFLAGS2) $(ARC_SRC) $(GENERAL_SRC) $(TIMER_SRC) $(PROBLEM_SRC) $(HEURISTIC_SOLUTION_SRC) $(FEASIBILITY_PUMP_SRC) $(KERNEL_PUMP_SRC) $(MAIN_SRC) -o $(PROG_BIN)/kp $(CLNFLAGS)

clean:
	rm $(PROG_BIN)/*
