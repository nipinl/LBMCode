/*
Number of di


*/


#include "D2Q4.h"
//constructor with arguments
D2Q4::D2Q4(const Material& m, const solverSettings& ss, bc& lbc, bc& rbc, bc& tbc, bc& bbc)//default Material created if not supplied
	:Nx(ss.Nx),
    Ny(ss.Ny),
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
        double_ptr_2D T(new double_ptr_1D[Nx]);
        double_ptr_2D f0(new double_ptr_1D[Nx]);
        double_ptr_2D f1(new double_ptr_1D[Nx]);
        double_ptr_2D f2(new double_ptr_1D[Nx]);
        double_ptr_2D f3(new double_ptr_1D[Nx]);
        double_ptr_2D f4(new double_ptr_1D[Nx]);
        
        for(size_t i=0;i<Nx;i++){
            T[i] = double_ptr_1D(new double[Ny]);
            f0[i] = double_ptr_1D(new double[Ny]);
            f1[i] = double_ptr_1D(new double[Ny]);
            f2[i] = double_ptr_1D(new double[Ny]);
            f3[i] = double_ptr_1D(new double[Ny]);
            f4[i] = double_ptr_1D(new double[Ny]);
        }
        omega = 1.0 / (alpha / (dt*cs2) + 0.5);//Eqn 5.27
		oneMinusOmega = 1.0 - omega;
        for (int i = 0; i < Nx - 1; i++) {
		    x[i + 1] = x[i] + dx;
	    }
        for (int i = 0; i < Ny - 1; i++) {
            y[i + 1] = y[i] + dy;
	    }
};
void D2Q4::collide(){
    double source = uniformHeatSource*alpha/k;
    double feq;
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            feq = T[i][j]*0.25;//All weights = 0.25, so using directly to save time
            f1[i][j] = oneMinusOmega*f1[i][j] + omega*feq  + source*0.25;//Eqn 5.21
            f2[i][j] = oneMinusOmega*f2[i][j] + omega*feq  + source*0.25;
            f3[i][j] = oneMinusOmega*f3[i][j] + omega*feq  + source*0.25;
            f4[i][j] = oneMinusOmega*f4[i][j] + omega*feq  + source*0.25;
        }
    }
}

void D2Q4::stream(){
    for (int j = 0; j <= Ny-1; j++){
        for (int i = 1; i <= Nx-1; i++){
            f1[Nx - i + 1][j] = f1[Nx - i][j];//f1: towards +ve x direction
            f2[i - 1][j] = f2[i][j];//f2: towards -ve x direction
        }
    }
    for (int i = 0; i <= Nx-1; i++){
        for (int j = 1; j <= Ny-1; j++){
            f3[i][Ny - j + 1] = f3[i][Ny - j];//f3: towards +ve y direction
            f4[i][j - 1] = f4[i][j];//f4: towards -ve y direction
        }
    }
}
void D2Q4::applyBc(){
    //left boundary
    if(lbc.type==0){//Dirichlet. T(0)=TLeft
        double TLeft = lbc.val1;
        for (int j = 0; j <= Ny-1; j++){
            f1[0][j] = TLeft - (f2[0][j]+f3[0][j]+f4[0][j]);//5.86 on page #73     
        }
    }
    if(lbc.type==1){//Neumann
        if(lbc.val1==0){// Adiabatic boundary condition, zero flux condition at left boundary
            for (int j = 0; j <= Ny-1; j++){
                f1[0][j] = f1[1][j];//next line to 5.89 on page #73
                f2[0][j] = f2[1][j];
                f3[0][j] = f3[1][j];
                f4[0][j] = f4[1][j];            
            }	
        }
        else{
            double q = lbc.val1;
            for (int j = 0; j <= Ny-1; j++){
                f1[0][j] = (f1[1][j]+f2[1][j]+f3[1][j]+f4[1][j]) + q/k - (f2[0][j]+f3[0][j]+f4[0][j]);//next line to 5.92.5.93           
            }
        }
    }
    //right boundary
    if(rbc.type==0){//Dirichlet. T[Nx-1][j]=TRight
        double TRight = rbc.val1;
        for (int j = 0; j <= Ny-1; j++){
            f2[Nx-1][j] = TRight - (f1[Nx-1][j]+f3[Nx-1][j]+f4[Nx-1][j]);//5.86 on page #73     
        }
    }
    if(rbc.type==1){//Neumann
        if(rbc.val1==0){// Adiabatic boundary condition, zero flux condition at right boundary
            for (int j = 0; j <= Ny-1; j++){
                f1[Nx-1][j] = f1[Nx-2][j];//next line to 5.89 on page #73
                f2[Nx-1][j] = f2[Nx-2][j];
                f3[Nx-1][j] = f3[Nx-2][j];
                f4[Nx-1][j] = f4[Nx-2][j];            
            }	
        }
        else{
            double q = rbc.val1;
            for (int j = 0; j <= Ny-1; j++){
                f2[Nx-1][j] = (f1[Nx-2][j]+f2[Nx-2][j]+f3[Nx-2][j]+f4[Nx-2][j]) + q/k - (f1[Nx-1][j]+f3[Nx-1][j]+f4[Nx-1][j]);//next line to 5.92.5.93           
            }
        }
    }

    //bottom boundary
    if(bbc.type==0){//Dirichlet. T[i][0]=TBottom
        double TBottom = bbc.val1;
        for (int i = 0; i <= Nx-1; i++){
            f3[i][0] = TBottom - (f1[i][0]+f2[i][0]+f4[i][0]);//5.86 on page #73     
        }
    }
    if(bbc.type==1){//Neumann
        if(bbc.val1==0){// Adiabatic boundary condition, zero flux condition at bottom boundary
            for (int i = 0; i <= Nx-1; i++){
                f1[i][0] = f1[i][1];//next line to 5.89 on page #73
                f2[i][0] = f2[i][1];
                f3[i][0] = f3[i][1];
                f4[i][0] = f4[i][1];            
            }	
        }
        else{
            double q = bbc.val1;
            for (int i = 0; i <= Nx-1; i++){
                f3[i][0] = (f1[i][1]+f2[i][1]+f3[i][1]+f4[i][1]) + q/k - (f1[i][0]+f2[i][0]+f4[i][0]);//next line to 5.92.5.93           
            }
        }
    }
    //top boundary
    if(tbc.type==0){//Dirichlet. T[i][0]=TTop
        double TTop = tbc.val1;
        for (int i = 0; i <= Nx-1; i++){
            f4[i][Ny-1] = TTop - (f1[i][Ny-1]+f2[i][Ny-1]+f3[i][Ny-1]);//5.86 on page #73     
        }
    }
    if(tbc.type==1){//Neumann
        if(tbc.val1==0){// Adiabatic boundary condition, zero flux condition at bottom boundary
            for (int i = 0; i <= Nx-1; i++){
                f1[i][Ny-1] = f1[i][Ny-2];//next line to 5.89 on page #73
                f2[i][Ny-1] = f2[i][Ny-2];
                f3[i][Ny-1] = f3[i][Ny-2];
                f4[i][Ny-1] = f4[i][Ny-2];            
            }	
        }
        else{
            double q = tbc.val1;
            for (int i = 0; i <= Nx-1; i++){
                f4[i][Ny-1] = (f1[i][Ny-2]+f2[i][Ny-2]+f3[i][Ny-2]+f4[i][Ny-2]) + q/k - (f1[i][Ny-1]+f2[i][Ny-1]+f3[i][Ny-1]);//next line to 5.92.5.93           
            }
        }
    }

}
void D2Q4::calculateT(){
    for (int i = 0; i < Nx; i++)
    {
        T[i] =  f1[i] + f2[i]+ f3[i] + f4[i];
    }
}
void D2Q4::solve(){
    for (int i = 0; i < endTime; i++)
    {
        collide();
        stream();
        applyBc();
        calculateT();
    }
    
}
void D2Q4::write(){
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


