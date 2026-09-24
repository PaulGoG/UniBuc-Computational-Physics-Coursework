#include<iostream>
#include<fstream>
using namespace std;

ofstream g ("EnergiiExamen.txt");

int main()
{
    float x,y,x1,y1,dt,E,ct; //ct=Omega^2
    int N,Nmax;
    ct=1;
    cout<<"Introduceti numarul maxim de pasi";
    cin>>Nmax;
    //g<<"N"<<"         "<<"E"<<'\n';
    for(N=1;N<=Nmax;N++)
    {
    x1=1; y1=-1; //conditii initiale
    dt=(float)10/N;// t merge de la 0 la 10
    for(int i=1; i<=N; i++)
    {
        y=y1+dt*(-ct*x1);
        x=x1+dt*y1;
        y1=y;
        x1=x;
    }

    E=0.5*y*y+0.5*ct*x*x; // 0.5*y'+ct*0.5*y
    //g<<N<<"         "<<E<<'\n';
    g<<E<<'\n';
    }
    return 0;
}
