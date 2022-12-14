#include "D1Q3.h"
//constructor with arguments
D1Q3::D1Q3(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc)//default Material created if not supplied
	:Nx(ss.Nx),
	endTime(ss.endTime),
    k(m.k),
    alpha(m.alpha),
    lbc(lbc),
    rbc(rbc),
    T(new double[Nx]),
	x(new double[Nx]),
	f0(new double[Nx]),
    f1(new double[Nx]),
    f2(new double[Nx])
	{
        omega = 1.0 / (alpha / (dt*cs2) + 0.5);//Eqn 5.27
		oneMinusOmega = 1.0 - omega;
        for (int i = 0; i < Nx - 1; i++) {
		    x[i + 1] = x[i] + dx;
	    }
	};
    void D1Q3::collide(){
        double source = uniformHeatSource*alpha/k;
    for (int i = 0; i < Nx; i++)
		{
			double feq0 = weights[0]*T[i];
			double feq  = weights[1]*T[i];
			f0[i] = oneMinusOmega*f0[i] + omega*feq0 + source*weights[0];
			f1[i] = oneMinusOmega*f1[i] + omega*feq  + source*weights[1];//Eqn 5.21
			f2[i] = oneMinusOmega*f2[i] + omega*feq  + source*weights[1];
		}
}

void D1Q3::stream(){
    for (int i = 1; i <= Nx-1; i++)
		{
			f1[Nx - i] = f1[Nx - i - 1];//f0 doesn't stream
			f2[i - 1] = f2[i];
        }
}
void D1Q3::applyBc(){
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
void D1Q3::calculateT(){
    for (int i = 0; i < Nx; i++)
    {
        T[i] = f0[i] + f1[i] + f2[i];
    }
}
void D1Q3::solve(){
    for (int i = 0; i < endTime; i++)
    {
        collide();
        stream();
        applyBc();
        calculateT();
    }
    
}
void D1Q3::write(){
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


