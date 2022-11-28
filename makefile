
d1q3: 2_d1q3.o D1Q3.o
	g++  2_d1q3.o D1Q3.o -o d1q3

2_d1q3.o: 2_d1q3.cpp
	g++ -c  2_d1q3.cpp

D1Q3.o: D1Q3.cpp
	g++ -c D1Q3.cpp

clean:
	rm *.o d1q3