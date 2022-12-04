#ifndef D2Q9_H_
#define D2Q9_H_
#include "D1Q3.h"
#include "D2Q4.h"
using namespace std;

class D2Q9: public D2Q4
{
protected:
	const double w0{4.0/9},w1{1.0/9},w5{1.0/36};
	const double w2{w1},w3{w1},w4{w1};
	const double w6{w5},w7{w5},w8{w5};

	const double cs2{1.0/3};

	double_ptr_2D f0,f5,f6,f7,f8;
public:
	D2Q9(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc);
	~D2Q9(){};
	void collide();
	void applyBc();
	void calculateT();	
	void circshift(double_ptr_2D arr, int xShift, int yShift);
	void stream();

};
#endif