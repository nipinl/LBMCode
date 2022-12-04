
main: main.o D2Q4.o D2Q5.o D2Q9.o
	g++  main.o D2Q4.o D2Q5.o D2Q9.o -o main

main.o: main.cpp
	g++ -c  main.cpp

D2Q4.o: D2Q4.cpp
	g++ -c D2Q4.cpp

D2Q5.o: D2Q5.cpp
	g++ -c D2Q5.cpp

D2Q9.o: D2Q9.cpp
	g++ -c D2Q9.cpp

clean:
	rm *.o main
