OBJECTS = main.o eigen_value_computation.o iterative_methods.o matrix.o cholesky.o func.o
CC = g++
CFLAGS = -Wall -Wextra -g



main: $(OBJECTS) 
	$(CC) -o main $(OBJECTS)

test: test.o $(OBJECTS)
	$(CC) -o test test.o matrix.o cholesky.o func.o

test.o: test.cpp
	$(CC) $(CFLAGS) -c test.cpp

iterative_methods.o: iterative_methods.cpp cholesky.h matrix.h
	$(CC) $(CFLAGS) -c iterative_methods.cpp

cholesky.o: cholesky.cpp matrix.h
	$(CC) $(CFLAGS) -c cholesky.cpp

eigen_value_computation.o: eigen_value_computation.cpp cholesky.h
	$(CC) $(CFLAGS) -c eigen_value_computation.cpp

matrix.o: matrix.cpp matrix.h 
	$(CC) $(CFLAGS) -c matrix.cpp

func.o: func.cpp 
	$(CC) $(CFLAGS) -c func.cpp

main.o: main.cpp cholesky.h func.h 
	$(CC) $(CFLAGS) -c main.cpp

clean:
	rm -rf  *.o
	rm -rf main
	rm -rf matrix
	rm -rf test