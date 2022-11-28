#ifndef D1Q3_H_
#define D1Q3_H_
#include<memory>
#include <iostream>
#include<fstream>
#include<vector>
#include <sys/stat.h>//to check if a file exist
using namespace std;
enum bcType {Dirichlet, Neumann, Mixed};
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
class Q3
{
double f[3]={0.0};
public:
	Q3(){};
	Q3(double val0,double val1,double val2 ){
		f[0]=val0;f[1]=val1;f[2]=val2;
	}
	Q3(const Q3& other){
		*this = other;
	}
	void setf(double val, int i){f[i]=val;};
	double operator[](int i){
		double val{-1};
		if(i>=0 && i<3) {val= f[i];}
		return(val);
	}
	Q3 & operator=(const Q3& other){
		f[0]=other.f[0];f[1]=other.f[1];f[2]=other.f[2];
		return *this;
	}
};
class D1Q3
{
private:
	const double weights[3]= {4.0/6, 1.0/6, 1.0/6};
	int nodes{100};
	const double length{1}, cs2{1.0/3};
	std::unique_ptr<double[]>T;
	std::unique_ptr<double[]>x;
	std::unique_ptr<double[]>f0;
	std::unique_ptr<double[]>f1;
	std::unique_ptr<double[]>f2;
	double dt{1.0},dx{1.0};
	int endTime{200};
	bc lbc, rbc;
	Material material;
	double omega{0},oneMinusOmega{1.0};
public:
	D1Q3(int nodes, int endTime);
	~D1Q3(){};
	void setLeftBC(bc leftbc){lbc=leftbc;}
	void setRightBC(bc rightbc){rbc=rightbc;}
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