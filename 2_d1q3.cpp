/*
5.5.2 Heat Diffusion in an Infinite Slab Subjected
to a Constant Temperature
*/
#include <iostream>
#include<memory>
#include "D1Q3.h"
using namespace std;

const int N {100};
const int endTime {200};
int main()
{
	D1Q3 d1q3(N,endTime);
	d1q3.setLeftBC(bc(Dirichlet,1.0));
	d1q3.setRightBC(bc(Neumann,0.0));
	d1q3.solve();
	d1q3.write();

	return 0;
}

