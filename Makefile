CXX := g++
FLAGS := -march=native -std=c++11 -O3 -lm 
OPENMP_FLAG := -fopenmp

BIN_PATH := ./bin
SRC_PATH := ./src

SRCS := $(wildcard ${SRC_PATH}/*.cpp ${SRC_PATH}/*.c)

parallel: 
	mkdir ./bin/summary
	mkdir ./bin/helix
	mkdir ./bin/press
	@echo -n "Compiling code... "
	@${CXX} ${SRCS} ${FLAGS} ${OPENMP_FLAG} -o ${BIN_PATH}/shape
	@echo "Done"

