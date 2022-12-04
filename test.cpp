#include<iostream>
#include<memory>
using namespace std;
using double_ptr_1D = std::unique_ptr<double[]>;
using double_ptr_2D = std::unique_ptr<double_ptr_1D[]>;
const int Nx{5},Ny{5};
void circshift(double arr[Nx][Ny], int xShift, int yShift)
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
    double val{0};
    for (int j =0; j < Ny; j++) 
    {
        jj = (j + yShift) % Ny;
        if (jj <0) jj = jj + Ny;
        for (int i = 0; i < Nx; i++) 
        {
            ii = (i + xShift) % Nx;
            if (ii <0) ii = ii + Nx;
            arr[j] [i] = copyArr[jj][ii];
        }
    }
}
void printArr(double arr[Nx][Ny])
{

    for (int j =0; j < Ny; j++) 
    {
        for (int i = 0; i < Nx; i++) 
        {
           cout<<arr[j][i]<<"   ";
        }
        cout<<endl;
    }
}
int main(){

    double a[5][5] = {  
                        {1, 2, 3, 4, 5},
                        {6, 7, 8, 9, 10},
                        {11, 12, 13, 14, 15},
                        {16, 17, 18, 19, 20},
                        {21, 22, 23, 24, 25}
                     };
    printArr(a);
    cout<<endl<<endl;
    circshift(a,0,1);
     printArr(a);
    
}


