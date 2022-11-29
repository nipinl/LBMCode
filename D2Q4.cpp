/*
Number of di


*/


#include "D2Q5.h"
//constructor with arguments
D2Q5::D2Q5(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc)//default Material created if not supplied
	:Nx(ss.Nx),
	endTime(ss.endTime),
    k(m.k),
    alpha(m.alpha),
    lbc(lbc),
    rbc(rbc),
    tbc(tbc),
    bbc(bbc),
	x(new double[Nx]),
    y(new double[Ny])
	{
        double_ptr_2D T(new double_ptr_1D[Ny]);
        double_ptr_2D f0(new double_ptr_1D[Ny]);
        double_ptr_2D f1(new double_ptr_1D[Ny]);
        double_ptr_2D f2(new double_ptr_1D[Ny]);
        double_ptr_2D f3(new double_ptr_1D[Ny]);
        double_ptr_2D f4(new double_ptr_1D[Ny]);
        
        for(size_t i=0;i<Ny;i++){
            T[i] = double_ptr_1D(new double[Nx]);
            f0[i] = double_ptr_1D(new double[Nx]);
            f1[i] = double_ptr_1D(new double[Nx]);
            f2[i] = double_ptr_1D(new double[Nx]);
            f3[i] = double_ptr_1D(new double[Nx]);
            f4[i] = double_ptr_1D(new double[Nx]);
        }
        omega = 1.0 / (alpha / (dt*cs2) + 0.5);//Eqn 5.27
		oneMinusOmega = 1.0 - omega;
        for (int i = 0; i < Nx - 1; i++) {
		    x[i + 1] = x[i] + dx;
            y[i + 1] = y[i] + dy;
	    }
	};
    void D2Q5::collide(){
        double source = uniformHeatSource*alpha/k;
    for (int i = 0; i < Nx; i++)
		{
			double feq0 = weights[0]*T[i];
			double feq  = weights[1]*T[i];
			f0[i] = oneMinusOmega*f0[i] + omega*feq0 + source*weights[0];
			f1[i] = oneMinusOmega*f1[i] + omega*feq  + source*weights[1];//Eqn 5.21
			f2[i] = oneMinusOmega*f2[i] + omega*feq  + source*weights[1];
            f3[i] = oneMinusOmega*f3[i] + omega*feq  + source*weights[1];//Eqn 5.21
			f4[i] = oneMinusOmega*f4[i] + omega*feq  + source*weights[1];
		}
}

void D2Q5::stream(){
    for (int i = 1; i <= Nx-1; i++)
		{
			f1[Nx - i] = f1[Nx - i - 1];//f0 doesn't stream
			f2[i - 1] = f2[i];
        }
}
void D2Q5::applyBc(){
    //left boundary
    if(lbc.type==0){//Dirichlet. T(0)=TLeft
        double TLeft = lbc.val1;
        f1[0] = TLeft - (f0[0]+f2[0]);
    }
    if(lbc.type==1){//Neumann
        if(lbc.val1==0){// Adiabatic boundary condition, zero flux condition at left boundary
            f0[0]=f0[1];
            f1[0]=f1[1];
            f2[0]=f2[1];	
        }
        else{
            double q = lbc.val1;
            f1[0] = (f0[1]+f1[1]+f2[1]) + q/k - (f0[0]+f2[0]);
        }
    }
    //right boundary
    if(rbc.type==0){//Dirichlet. T(Nx-1)=TRight
        double TRight = rbc.val1;
        f2[Nx-1] = TRight - (f0[Nx-1]+f1[Nx-1]);
    }
    if(rbc.type==1){
        if(rbc.val1==0){//Adiabatic boundary condition, zero flux condition at right boundary
            f0[Nx-1]=f0[Nx-2];
            f1[Nx-1]=f1[Nx-2];
            f2[Nx-1]=f2[Nx-2];	
        }
        else{
            double q = rbc.val1;
             f2[Nx-1] = (f0[Nx-2]+f1[Nx-2]+f2[Nx-2]) + q/k - (f0[Nx-1]+f1[Nx-1]);
        }
    }
}
void D2Q5::calculateT(){
    for (int i = 0; i < Nx; i++)
    {
        T[i] = f0[i] + f1[i] + f2[i]+ f3[i] + f4[i];
    }
}
void D2Q5::solve(){
    for (int i = 0; i < endTime; i++)
    {
        collide();
        stream();
        applyBc();
        calculateT();
    }
    
}
void D2Q5::write(){
    std::ofstream outFile;
    if(exists("T")) remove("T");
	outFile.open("T");
	for (int i = 0; i < Nx; i++)
	{
		outFile << x[i] << "\t" << T[i] << endl;
	}
	outFile.close();
	system("./test.sh");
}


