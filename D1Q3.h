#ifndef D1Q3_H_
#define D1Q3_H_
#include<memory>
#include <iostream>
#include<fstream>
#include<vector>
#include <sys/stat.h>//to check if a file exist
using namespace std;
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
//solverSetting has Nx(Number of nodes in x direction),Ny(you guessed it) and endTime(end time)
class solverSettings{
	public:
	int endTime{200}, Nx{100}, Ny{100};
	solverSettings(){};
	solverSettings(int endTime,int Nx):endTime(endTime),Nx(Nx){}
	solverSettings(int endTime,int Nx, int Ny):endTime(endTime),Nx(Nx),Ny(Ny){}
	solverSettings(const solverSettings& ss){
		endTime= ss.endTime;
		Nx=ss.Nx;
		Ny=ss.Ny;
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

class D1Q3
{
private:
	const double weights[3]= {4.0/6, 1.0/6, 1.0/6};
	int Nx{100};
	const double length{1}, cs2{1.0/3};
	std::unique_ptr<double[]>T;
	std::unique_ptr<double[]>x;
	std::unique_ptr<double[]>f0;
	std::unique_ptr<double[]>f1;
	std::unique_ptr<double[]>f2;
	double dt{1.0},dx{1.0};
	int endTime{200};
	bc lbc, rbc;
	double k{0.25}, alpha{0.5};
	double uniformHeatSource{0.0};
	double omega{0},oneMinusOmega{1.0};
public:
	D1Q3(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc);
	~D1Q3(){};
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