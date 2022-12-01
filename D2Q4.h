#ifndef D2Q4_H_
#define D2Q4_H_
#include<memory>
#include <iostream>
#include<fstream>
#include <sys/stat.h>//to check if a file exist
using namespace std;
using double_ptr_1D = std::unique_ptr<double[]>;
using double_ptr_2D = std::unique_ptr<double_ptr_1D[]>;
//using double_ptr_3D = std::unique_ptr<double_ptr_2D[]>;

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
	int endTime{200}, Nx{100}, Ny{100};
	solverSettings(){};
	solverSettings(int endTime,int Nx):Nx(Nx),endTime(endTime){}
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
class D2Q4
{
private:
	const double weights[4]= {0.25, 0.25, 0.25, 0.25};//For clarity. All w = 0.25. 0.25 is used instead of w{i}
	int Nx{100},Ny{100};
	const double /* length{1.0}, */ cs2{1.0/3};

	double_ptr_2D T;
	double_ptr_2D f1,f2,f3,f4;
	std::unique_ptr<double[]>x;
	std::unique_ptr<double[]>y;
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