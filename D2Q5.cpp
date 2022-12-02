#include "D2Q5.h"
//constructor with arguments
D2Q5::D2Q5(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc)//default Material created if not supplied
	:D2Q4( m, ss, lbc, rbc, tbc, bbc),
    f0(new double_ptr_1D[Nx])
	{
        for(size_t i=0;i<Nx;i++){
            f0[i] = double_ptr_1D(new double[Ny]);
        }
    };
void D2Q5::collide(){
    double source = uniformHeatSource*alpha/k;
    double feq,feq0;
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            feq0 = T[i][j]*w0;
            feq = T[i][j]*w1;
            f0[i][j] = oneMinusOmega*f1[i][j] + omega*feq0  + source*0.25;
            f1[i][j] = oneMinusOmega*f1[i][j] + omega*feq  + source*0.25;//Eqn 5.21
            f2[i][j] = oneMinusOmega*f2[i][j] + omega*feq  + source*0.25;
            f3[i][j] = oneMinusOmega*f3[i][j] + omega*feq  + source*0.25;
            f4[i][j] = oneMinusOmega*f4[i][j] + omega*feq  + source*0.25;
        }
    }
}

//stream function is same as of D2Q4

void D2Q5::applyBc(){
    //left boundary
    if(lbc.type==0){//Dirichlet. T(0)=TLeft
        double TLeft = lbc.val1;
        for (int j = 0; j < Ny; j++){
            f1[0][j] = TLeft - (f0[0][j]+f2[0][j]+f3[0][j]+f4[0][j]);//5.86 on page #73     
        }
    }
    if(lbc.type==1){//Neumann
        if(lbc.val1==0){// Adiabatic boundary condition, zero flux condition at left boundary
            for (int j = 0; j < Ny; j++){
                f0[0][j] = f0[1][j];
                f1[0][j] = f1[1][j];//next line to 5.89 on page #73
                f2[0][j] = f2[1][j];
                f3[0][j] = f3[1][j];
                f4[0][j] = f4[1][j];            
            }	
        }
        else{
            double q = lbc.val1;
            for (int j = 0; j < Ny; j++){
                f1[0][j] = (f0[1][j]+f1[1][j]+f2[1][j]+f3[1][j]+f4[1][j]) + q/k - (f0[0][j]+f2[0][j]+f3[0][j]+f4[0][j]);//next line to 5.92.5.93           
            }
        }
    }
    //right boundary
    if(rbc.type==0){//Dirichlet. T[Nx-1][j]=TRight
        double TRight = rbc.val1;
        for (int j = 0; j < Ny; j++){
            f2[Nx-1][j] = TRight - (f0[Nx-1][j]+f1[Nx-1][j]+f3[Nx-1][j]+f4[Nx-1][j]);//5.86 on page #73     
        }
    }
    if(rbc.type==1){//Neumann
        if(rbc.val1==0){// Adiabatic boundary condition, zero flux condition at right boundary
            for (int j = 0; j < Ny; j++){
                f0[Nx-1][j] = f0[Nx-2][j];
                f1[Nx-1][j] = f1[Nx-2][j];//next line to 5.89 on page #73
                f2[Nx-1][j] = f2[Nx-2][j];
                f3[Nx-1][j] = f3[Nx-2][j];
                f4[Nx-1][j] = f4[Nx-2][j];            
            }	
        }
        else{
            double q = rbc.val1;
            for (int j = 0; j < Ny; j++){
                f2[Nx-1][j] = (f0[Nx-2][j]+f1[Nx-2][j]+f2[Nx-2][j]+f3[Nx-2][j]+f4[Nx-2][j]) + q/k - (f0[Nx-1][j]+f1[Nx-1][j]+f3[Nx-1][j]+f4[Nx-1][j]);//next line to 5.92.5.93           
            }
        }
    }

    //bottom boundary
    if(bbc.type==0){//Dirichlet. T[i][0]=TBottom
        double TBottom = bbc.val1;
        for (int i = 0; i < Nx; i++){
            f3[i][0] = TBottom - (f0[i][0]+f1[i][0]+f2[i][0]+f4[i][0]);//5.86 on page #73     
        }
    }
    if(bbc.type==1){//Neumann
        if(bbc.val1==0){// Adiabatic boundary condition, zero flux condition at bottom boundary
            for (int i = 0; i < Nx; i++){
                f0[i][0] = f0[i][1];
                f1[i][0] = f1[i][1];//next line to 5.89 on page #73
                f2[i][0] = f2[i][1];
                f3[i][0] = f3[i][1];
                f4[i][0] = f4[i][1];            
            }	
        }
        else{
            double q = bbc.val1;
            for (int i = 0; i < Nx; i++){
                f3[i][0] = (f0[i][1]+f1[i][1]+f2[i][1]+f3[i][1]+f4[i][1]) + q/k - (f1[i][0]+f2[i][0]+f4[i][0]);//next line to 5.92.5.93           
            }
        }
    }
    //top boundary
    if(tbc.type==0){//Dirichlet. T[i][0]=TTop
        double TTop = tbc.val1;
        for (int i = 0; i < Nx; i++){
            f4[i][Ny-1] = TTop - (f0[i][Ny-1]+f1[i][Ny-1]+f2[i][Ny-1]+f3[i][Ny-1]);//5.86 on page #73     
        }
    }
    if(tbc.type==1){//Neumann
        if(tbc.val1==0){// Adiabatic boundary condition, zero flux condition at bottom boundary
            for (int i = 0; i < Nx; i++){
                f0[i][Ny-1] = f0[i][Ny-2];
                f1[i][Ny-1] = f1[i][Ny-2];//next line to 5.89 on page #73
                f2[i][Ny-1] = f2[i][Ny-2];
                f3[i][Ny-1] = f3[i][Ny-2];
                f4[i][Ny-1] = f4[i][Ny-2];            
            }	
        }
        else{
            double q = tbc.val1;
            for (int i = 0; i < Nx; i++){
                f4[i][Ny-1] = (f0[i][Ny-2]+f1[i][Ny-2]+f2[i][Ny-2]+f3[i][Ny-2]+f4[i][Ny-2]) + q/k - (f0[i][Ny-1]+f1[i][Ny-1]+f2[i][Ny-1]+f3[i][Ny-1]);//next line to 5.92.5.93           
            }
        }
    }

}
void D2Q5::calculateT(){
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++)
        {
            T[i][j] = f0[i][j] +f1[i][j] + f2[i][j] + f3[i][j] + f4[i][j];
        }
    }
}

