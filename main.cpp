
#include <iostream>
#include<memory>
//#include "D1Q3.h"
#include "D2Q5.h"
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
	solverSettings ss(200,100);
	bc lbc(Dirichlet,1);
	bc rbc(Neumann,0);
	D1Q3 d1q3(m,ss,lbc,rbc);
	d1q3.solve();
	d1q3.write();
	*/

//5.5.5 Heat Diffusion in an Infinite Slab Subjected to a Constant Heat flux
	
	/* Material m(20.0,0.25);
	solverSettings ss(200,100);
	bc lbc(Neumann,100);
	bc rbc(Neumann,0);
	D1Q3 d1q3(m,ss,lbc,rbc);
	d1q3.setUniformHeatSource(1);
	d1q3.solve();
	d1q3.write(); */

	Material m(20.0,0.25);
	solverSettings ss(200,50,50);
	bc lbc(Dirichlet,1);
	bc rbc(Dirichlet,0);
	bc bbc(Dirichlet,1);
	bc tbc(Dirichlet,1);

	D2Q4 d2q4(m,ss,lbc,rbc,tbc,bbc);
	d2q4.initialize(2);
	d2q4.solve();
	d2q4.animate();

	/* D2Q5 d2q5(m,ss,lbc,rbc,tbc,bbc);
	d2q5.initialize(2);
	d2q5.solve();
	d2q5.animate(); */

	return 0;
}

