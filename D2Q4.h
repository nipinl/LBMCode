#ifndef D1Q3_H_
#define D1Q3_H_
#include<memory>
#include <iostream>
#include<fstream>
#include<vector>
#include <sys/stat.h>//to check if a file exist
using namespace std;
using double_ptr_1D = std::unique_ptr<double[]>;
using double_ptr_2D = std::unique_ptr<double_ptr_1D[]>;

double_ptr_2D dPtr( new double_ptr_1D[10] );

enum bcType {Dirichlet, Neumann, Mixed};
//Dirichlet: constant value of variablex
//Neumann: derivative of variable
//Mixed: derivative = a.variable + b

//Material gives thermal conductivity(k) and thermal diffusivity (alpha)
class Material{
	public:
	double k{0.5},alpha{0.25};
	Material(){};
	Material(double k, double alpha):k(k),alpha(alpha){}
	Material(const Material& m){
		k=m.k;
		alpha= m.alpha;
	}
};
//solverSetting has N(Number of Nx) and endTime(end time)
class solverSettings{
	public:
	int Nx{100}, endTime{200};
	solverSettings(){};
	solverSettings(int Nx, int endTime):Nx(Nx),endTime(endTime){}
	solverSettings(const solverSettings& ss){
		Nx=ss.Nx;
		endTime= ss.endTime;
	}
};

class bc
{
public:
	double val1{1.0},val2{0.0};
	bcType type{Dirichlet};
	bc():type(Dirichlet),val1(1){};
	bc(bcType bt, double val=0.0):type(bt), val1(val){};
	//copy constructor
	bc(const bc& other){*this = other;}
	//overloading assignment operator
	bc& operator=(const bc& other){
		type=other.type;
		val1=other.val1;
		val2=other.val2;
		return *this;
	}
	~bc(){};
};
class D2Q5
{
private:
	const double weights[5]= {2.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6};//For clarity. w0=2.0/6, w1=1.0/6 is enough
	int Nx{100},Ny{100};
	const double length{1}, cs2{1.0/3};

	
	std::unique_ptr<double[]>T;
	std::unique_ptr<double[]>x;
	std::unique_ptr<double[]>y;
	std::unique_ptr<double[]>f0;
	std::unique_ptr<double[]>f1;
	std::unique_ptr<double[]>f2;
	std::unique_ptr<double[]>f3;
	std::unique_ptr<double[]>f4;
	double dt{1.0},dx{1.0},dy{1.0};
	int endTime{200};
	bc lbc, rbc, tbc, bbc;
	double k{0.25}, alpha{0.5};
	double uniformHeatSource{0.0};
	double omega{0},oneMinusOmega{1.0};
public:
	D2Q5(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc);
	~D2Q5(){};
	void setUniformHeatSource(double uniformHeatSource){uniformHeatSource=uniformHeatSource;}
	void calculateT();
	void applyBc();
	void collide();
	void stream();
	void solve();
	void write();
	inline bool exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); //stat rerun 0 if file exists else -1
}
};
#endif