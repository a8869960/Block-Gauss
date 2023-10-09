FLAGS=-O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -c

all: main.o matrixInput.o gauss_func.o matrixOperations.o matrixOutput.o results.o otherFunctions.o inverseMatrix.o forGauss.o
	g++ main.o matrixInput.o gauss_func.o matrixOperations.o matrixOutput.o results.o otherFunctions.o inverseMatrix.o forGauss.o

main.o: main.cpp
	g++ $(FLAGS) main.cpp

matrixInput.o: matrixInput.cpp
	g++ $(FLAGS) matrixInput.cpp

gauss_func.o: gauss_func.cpp
	g++ $(FLAGS) gauss_func.cpp

matrixOperations.o: matrixOperations.cpp
	g++ $(FLAGS) matrixOperations.cpp

matrixOutput.o: matrixOutput.cpp
	g++ $(FLAGS) matrixOutput.cpp

results.o: results.cpp
	g++ $(FLAGS) results.cpp

otherFunctions.o: otherFunctions.cpp
	g++ $(FLAGS) otherFunctions.cpp

inverseMatrix.o: inverseMatrix.cpp
	g++ $(FLAGS) inverseMatrix.cpp

forGauss.o: forGauss.cpp
	g++ $(FLAGS) forGauss.cpp

clean:
	rm -f *.out *.o *.gch

