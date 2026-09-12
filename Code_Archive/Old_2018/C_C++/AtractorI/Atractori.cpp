#include "helpers.h"
#include "Integratori.h"
#include<iostream>

using namespace std;

int main()
{
    double x,y,z,h,tMAX;
    int N;
    short Choose=1;

    cout<<"Variabile: N-nr. pasi, x0,y0,z0-primul set de conditii initiale,"<<'\n'<<"tMAX-intervalul de timpi (de la 0 la ceva): "<<'\n'<<'\n'<<"    ";

    cin>>N;
    cin>>x>>y>>z;
    cin>>tMAX;

    //cout<<"1-Lorentz, 0-Rossler "<<'\n'<<"  ";
    //cin>>Choose;

    h=(double)tMAX/N;

    if(Choose)
    RK4L(x,y,z,h,N);

    else
    RK4R(x,y,z,h,N);

    Distanta(h,N);
    extremelocale(h,N);

    return 0;
}
