all: main.o matrix.o task.o functions.o
	g++ -pthread main.o matrix.o task.o functions.o

main.o: main.cpp matrix.h functions.h
	g++ -pthread -Wall -std=c++20 -O3 -c main.cpp

matrix.o: matrix.cpp matrix.h
	g++ -pthread -Wall -std=c++20 -O3 -c matrix.cpp

task.o: task.cpp task.h
	g++ -pthread -Wall -std=c++20 -O3 -c task.cpp

functions.o: functions.cpp functions.h
	g++ -pthread -Wall -std=c++20 -O3 -c functions.cpp

clean:
	rm -f *.o *.out
