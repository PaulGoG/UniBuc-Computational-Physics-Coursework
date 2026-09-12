#include "helpers.h"
#include<cmath>
#include<fstream>
#include<iomanip>
#include<cmath>

using namespace std;

ifstream m ("C:\\Users\\Gogu\\Desktop\\IOFolder\\Atractor1.txt");
ifstream n ("C:\\Users\\Gogu\\Desktop\\IOFolder\\Atractor2.txt");
ofstream d ("C:\\Users\\Gogu\\Desktop\\IOFolder\\Distanta.txt");
ifstream p ("C:\\Users\\Gogu\\Desktop\\IOFolder\\Distanta.txt");
ofstream timp ("C:\\Users\\Gogu\\Desktop\\IOFolder\\TimpiUC.txt");

void Distanta(double h, int &N)
{
    double x1,y1,z1,x2,y2,z2,t=0;
    int Nskip=30/h;

    d.precision(15);
    for(int i=1;i<Nskip;i++)
    {
        m>>x1>>y1>>z1;
        n>>x2>>y2>>z2;
        t+=h;
    }
    for(int i=Nskip;i<=N-Nskip;i++)
    {
        m>>x1>>y1>>z1;
        n>>x2>>y2>>z2;
        d<<sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))<<' '<<t<<'\n';
        t+=h;
    }
    N-=Nskip;
}
void extremelocale(double h, int N)
{
    double T,D1,D2,D3,te=0,De,epsilon=0.0001;
    short low,high;
    p>>D1>>T;
    p>>D2>>T;
    p>>D3>>T;
    De=D1;
    N-=3;
    low=3;
    high=2;

    timp.precision(15);

    // Skip primul minim

    do
    {
        N--;
        if(D2>D1) low=1;
        if(D2<D1) low=0;
        if(D2>D3) high=1;
        if(D2<D3) high=0;
        if(low==high && abs(De-D2)/De>epsilon)
        {
            if(low==1) timp<<T-te-h<<' ';
            te=T-h;
            De=D2;
        }
        D1=D2;
        D2=D3;
        low=2;
        high=3;
    }while(te==0 && p>>D3 && p>>T);

    for(int i=1;i<=N;i++)
    {
        if(D2>D1) low=1;
        if(D2<D1) low=0;
        if(D2>D3) high=1;
        if(D2<D3) high=0;
        if(low==high && abs(De-D2)/De>epsilon)
        {
            timp<<T-te-h;
            te=T-h;
            if(low==1) timp<<' ';
                else timp<<'\n';
            De=D2;
        }
        if(D2!=D3) D1=D2;
        D2=D3;
        low=2;
        high=3;
        p>>D3;
        p>>T;
    }
}

