#include "D2Q9.h"
//constructor with arguments
D2Q9::D2Q9(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc)//default Material created if not supplied
	:D2Q4( m, ss, lbc, rbc, tbc, bbc),
    f0(new double_ptr_1D[Nx]),
    f5(new double_ptr_1D[Nx]),
    f6(new double_ptr_1D[Nx]),
    f7(new double_ptr_1D[Nx]),
    f8(new double_ptr_1D[Nx])
	{
        for(size_t i=0;i<Nx;i++){
            f0[i] = double_ptr_1D(new double[Ny]);
            f5[i] = double_ptr_1D(new double[Ny]);
            f6[i] = double_ptr_1D(new double[Ny]);
            f7[i] = double_ptr_1D(new double[Ny]);
            f8[i] = double_ptr_1D(new double[Ny]);
        }
    };
void D2Q9::collide(){
    double source = uniformHeatSource*alpha/k;
    double feq,feq0;
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            feq0 = T[i][j]*w0;
            feq = T[i][j]*w1;

            f0[i][j] = oneMinusOmega*f1[i][j] + omega*feq0  + source*w0;//center
            f1[i][j] = oneMinusOmega*f1[i][j] + omega*feq  + source*w1;//along one axis
            f2[i][j] = oneMinusOmega*f2[i][j] + omega*feq  + source*w2;
            f3[i][j] = oneMinusOmega*f3[i][j] + omega*feq  + source*w3;
            f4[i][j] = oneMinusOmega*f4[i][j] + omega*feq  + source*w4;
            f5[i][j] = oneMinusOmega*f1[i][j] + omega*feq  + source*w5;//diagonal
            f6[i][j] = oneMinusOmega*f2[i][j] + omega*feq  + source*w6;
            f7[i][j] = oneMinusOmega*f3[i][j] + omega*feq  + source*w7;
            f8[i][j] = oneMinusOmega*f4[i][j] + omega*feq  + source*w8;
        }
    }
}

void D2Q9::stream(){
    //inner nodes: no check needed to stream
    for (int j = 1; j < Ny-1; j++){
        for (int i = 1; i < Nx-1; i++){
            circshift(f1,1 , 0);
            circshift(f2,0 , 1);
            circshift(f3,-1, 0);
            circshift(f4,0 ,-1);
            circshift(f5,1 , 1);
            circshift(f6,-1, 1);
            circshift(f7,-1,-1);
            circshift(f8,-1, 1);
        }
    }
}

void D2Q9::applyBc(){
    //left boundary
    if(lbc.type==0){//Dirichlet. T(0)=TLeft
        double TLeft = lbc.val1;
        for (int j = 0; j < Ny; j++){
            f1[0][j] = w1*TLeft + w3*TLeft - f3[0][j];
            f5[0][j] = w5*TLeft + w7*TLeft - f7[0][j];
            f8[0][j] = w8*TLeft + w6*TLeft - f6[0][j];
        }
    }
    if(lbc.type==1){//Neumann
        if(lbc.val1==0){// Adiabatic boundary condition, zero flux condition at left boundary
            for (int j = 0; j < Ny; j++){
                f1[0][j] = f1[1][j];
                f5[0][j] = f5[1][j];
                f8[0][j] = f8[1][j];           
            }	
        }
    }
    //right boundary
    if(rbc.type==0){//Dirichlet. T[Nx-1][j]=TRight
        double TRight = rbc.val1;
        for (int j = 0; j < Ny; j++){
            f3[Nx-1][j] = w1*TRight + w3*TRight - f1[0][j];
            f7[Nx-1][j] = w5*TRight + w7*TRight - f5[0][j];
            f6[Nx-1][j] = w8*TRight + w6*TRight - f8[0][j];
        }
    }
    if(rbc.type==1){//Neumann
        if(rbc.val1==0){// Adiabatic boundary condition, zero flux condition at right boundary
            for (int j = 0; j < Ny; j++){
                f3[Nx-1][j] = f3[Nx-2][j];
                f7[Nx-1][j] = f7[Nx-2][j];
                f6[Nx-1][j] = f6[Nx-2][j];             
            }	
        }
    }

    //bottom boundary
    if(bbc.type==0){//Dirichlet. T[i][0]=TBottom
        double TBottom = bbc.val1;
        for (int i = 0; i < Nx; i++){
            f2[i][0] = w4*TBottom + w2*TBottom - f4[i][0];
            f5[i][0] = w5*TBottom + w7*TBottom - f7[i][0];
            f6[i][0] = w8*TBottom + w6*TBottom - f8[i][0];
        }
    }
    if(bbc.type==1){//Neumann
        if(bbc.val1==0){// Adiabatic boundary condition, zero flux condition at bottom boundary
            for (int i = 0; i < Nx; i++){
                f2[i][0] = f2[i][1];
                f5[i][0] = f5[i][1];
                f6[i][0] = f6[i][1];           
            }	
        }
    }
    //top boundary
    if(tbc.type==0){//Dirichlet. T[i][0]=TTop
        double TTop = tbc.val1;
        for (int i = 0; i < Nx; i++){
            f4[i][Ny-1] = w4*TTop + w2*TTop - f2[i][Ny-1];
            f7[i][Ny-1] = w5*TTop + w7*TTop - f5[i][Ny-1];
            f8[i][Ny-1] = w8*TTop + w6*TTop - f6[i][Ny-1];

            
        }
    }
    if(tbc.type==1){//Neumann
        if(tbc.val1==0){// Adiabatic boundary condition, zero flux condition at bottom boundary
            for (int i = 0; i < Nx; i++){
                f4[i][Ny-1] = f4[i][Ny-2];
                f7[i][Ny-1] = f7[i][Ny-2];
                f8[i][Ny-1] = f8[i][Ny-2];            
            }	
        }
    }

}
void D2Q9::calculateT(){
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++)
        {
            T[i][j] = f0[i][j]+f1[i][j]+f2[i][j]+f3[i][j]+f4[i][j]+f5[i][j]+f6[i][j]+f7[i][j]+f8[i][j];
        }
    }
}

void D2Q9::circshift(double_ptr_2D arr, int xShift, int yShift)
{
    double_ptr_2D copyArr {new double_ptr_1D[Nx]};
    for(size_t i=0;i<Nx;i++){
            copyArr[i] = double_ptr_1D(new double[Ny]);
    }
    for (int j =0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            copyArr[j][i] = arr[j][i];
        }
    }
    int ii, jj;
    for (int j =0; j < Ny; j++){
        jj = (j + yShift) % Ny;
        if (jj <0) jj = jj + Ny;
        for (int i = 0; i < Nx; i++){
            ii = (i + xShift) % Nx;
            if (ii <0) ii = ii + Nx;
            arr[j] [i] = copyArr[jj][ii];
        }
    }
}

