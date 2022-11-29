
#include <iostream>
#include<memory>
#include "D1Q3.h"
using namespace std;

int main()
{
	/*
	A problem should have 
	1) a material(k, alpha) 
	2) solverSetting has N(Number of nodes) and endTime(end time) 
	3) boundary conditions. Two bcs for 1d problems; first bc is for left and other is for right.
	*/

//5.5.2 Heat Diffusion in an Infinite Slab Subjected to a Constant Temperature
/*
	Material m(0.25,0.5);
	solverSettings ss(100,200);
	bc lbc(Dirichlet,1);
	bc rbc(Neumann,0);
	D1Q3 d1q3(m,ss,lbc,rbc);
	d1q3.solve();
	d1q3.write();
	*/

//5.5.5 Heat Diffusion in an Infinite Slab Subjected to a Constant Heat flux
	
	Material m(20.0,0.25);
	solverSettings ss(100,200);
	bc lbc(Neumann,100);
	bc rbc(Neumann,0);
	D1Q3 d1q3(m,ss,lbc,rbc);
	d1q3.setUniformHeatSource(1);
	d1q3.solve();
	d1q3.write();
	

	return 0;
}

