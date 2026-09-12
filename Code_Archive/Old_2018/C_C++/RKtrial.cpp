#include<iostream>
#include<cmath>
using namespace std;
double functie(double x, double y)
{
    return 3*exp(-x)-0.4*y;
}
double Euler(double x, double y, double h, int N)
{
    for(int i=1;i<=N;i++)
    {
        y=y+h*functie(x,y);
        x+=h;
    }
    return y;
}
double RK4(double x, double y, double h, int N)
{
    double k1,k2,k3,k4;
    for(int i=1;i<=N;i++)
    {
        k1=functie(x,y);
        k2=functie(x+h/2,y+h*k1/2);
        k3=functie(x+h/2,y+h*k2/2);
        k4=functie(x+h,y+h*k3);
        y=y+(h/6)*(k1+2*k2+2*k3+k4);
        x+=h;
    }
    return y;
}
int main()
{
    double x0,y0,x,y,h,xM;
    int N;
    cout<<"Variabile: N, x0, xM, y0 "<<'\n';
    cin>>N;
    cin>>x0;
    cin>>xM;
    cin>>y0;
    h=(double)(xM-x0)/N;
    x=x0;
    y=y0;
    cout<<Euler(x,y,h,N)<<'\n';
    cout<<RK4(x,y,h,N)<<'\n';
    return 0;
}
