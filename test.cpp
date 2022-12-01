#include<iostream>
#include<memory>
using namespace std;
using double_ptr_1D = std::unique_ptr<double[]>;
using double_ptr_2D = std::unique_ptr<double_ptr_1D[]>;
/* class Car
{
public:
int year{2000};
    Car(int y):year(y){};
    ~Car();
}; */
int main(){
    int m,n;
    double d{1.1};
    cout<<"Enter number of rows "<<endl;
    cin>>m;
    cout<<"Enter number of columns"<<endl;
    cin>>n;
    double_ptr_2D ptr1(new double_ptr_1D[m]);
    /* double_ptr_2D ptr1;
    double_ptr_2D ptr1(new double_ptr_1D[m]); */
    for(size_t i=0;i<m;i++){
        ptr1[i] = double_ptr_1D(new double[n]);
    }
    for(size_t i=0;i<m;i++){
        for(size_t j=0;j<n;j++){
        ptr1[i][j] = d;
        }
    }
    for(size_t i=0;i<m;i++){
        for(size_t j=0;j<n;j++){
        cout<<ptr1[i][j]<<" ";
        }
        cout<<endl;
    }
    return 0;
}


