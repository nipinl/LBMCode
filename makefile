
main: main.o D2Q4.o
	g++  main.o D2Q4.o -o main

main.o: main.cpp
	g++ -c  main.cpp

D2Q4.o: D2Q4.cpp
	g++ -c D2Q4.cpp

clean:
	rm *.o main
