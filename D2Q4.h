#ifndef D2Q4_H_
#define D2Q4_H_
#include<memory>
#include <iostream>
#include<fstream>
#include <sys/stat.h>//to check if a file exist
#include "D1Q3.h"
using namespace std;
using double_ptr_1D = std::shared_ptr<double[]>;
using double_ptr_2D = std::shared_ptr<double_ptr_1D[]>;
//using double_ptr_3D = std::shared_ptr<double_ptr_2D[]>;

class D2Q4
{
protected:
	const double w0{0.25};//For clarity. All w = 0.25. 0.25 is used instead of w{i}
	int Nx{100},Ny{100};
	const double /* length{1.0}, */ cs2{1.0/2};

	double_ptr_2D T;
	double_ptr_2D f1,f2,f3,f4;
	std::shared_ptr<double[]>x;
	std::shared_ptr<double[]>y;
	double dt{1.0},dx{1.0},dy{1.0};
	int endTime{200};
	bc lbc, rbc, tbc, bbc;
	double k{0.25}, alpha{0.5};
	double uniformHeatSource{0.0};
	double omega{0},oneMinusOmega{1.0};
public:
	D2Q4(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc);
	~D2Q4(){};
	void initialize(double initialTemperature);
	void setUniformHeatSource(double uniformHeatSource){uniformHeatSource=uniformHeatSource;}
	void calculateT();
	void applyBc();
	void collide();
	void stream();
	void solve();
	void write();
	void animate();
	inline bool exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); //stat rerun 0 if file exists else -1
}
};
#endif