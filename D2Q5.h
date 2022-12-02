#ifndef D2Q5_H_
#define D2Q5_H_

#include "D2Q4.h"
using namespace std;

class D2Q5: public D2Q4
{
protected:
	const double w0{2.0/6},w1{1.0/6};
	const double cs2{1.0/3};

	double_ptr_2D f0;
public:
	D2Q5(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc);
	~D2Q5(){};
	void collide();
	void applyBc();
	void calculateT();	
};
#endif