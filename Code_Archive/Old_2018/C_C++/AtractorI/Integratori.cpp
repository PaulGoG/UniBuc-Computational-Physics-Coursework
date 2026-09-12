#include "Integratori.h"
#include<fstream>
#include<iomanip>

#define sigma 10
#define beta 2.66
#define ro 28

#define a 0.2
#define b 0.2
#define c 5.7

using namespace std;

ofstream f ("C:\\Users\\Gogu\\Desktop\\IOFolder\\Atractor1.txt");
ofstream g ("C:\\Users\\Gogu\\Desktop\\IOFolder\\Atractor2.txt");

double functiexL(double x, double y)
{
    return sigma*(y-x);
}
double functieyL(double x, double y, double z)
{
    return x*(ro-z)-y;
}
double functiezL(double x, double y, double z)
{
    return x*y-beta*z;
}
void RK4L(double x, double y, double z, double h, int N)
{
    double k1[4], k3[4], k2[4], k4[4];
    double x0=x,y0=y,z0=z;

    f.precision(15);

    for(int i=1;i<=N;i++)
    {
        k1[1]=functiexL(x,y);
        k1[2]=functieyL(x,y,z);
        k1[3]=functiezL(x,y,z);

        k2[1]=functiexL(x+h*k1[1]/2,y+h*k1[2]/2);
        k2[2]=functieyL(x+h*k1[1]/2,y+h*k1[2]/2,z+h*k1[3]/2);
        k2[3]=functiezL(x+h*k1[1]/2,y+h*k1[2]/2,z+h*k1[3]/2);

        k3[1]=functiexL(x+h*k2[1]/2,y+h*k2[2]/2);
        k3[2]=functieyL(x+h*k2[1]/2,y+h*k2[2]/2,z+h*k2[3]/2);
        k3[3]=functiezL(x+h*k2[1]/2,y+h*k2[2]/2,z+h*k2[3]/2);

        k4[1]=functiexL(x+h*k3[1],y+h*k3[2]);
        k4[2]=functieyL(x+h*k3[1],y+h*k3[2],z+h*k3[3]);
        k4[3]=functiezL(x+h*k3[1],y+h*k3[2],z+h*k3[3]);

        x=x+(h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        y=y+(h/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        z=z+(h/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);
        f<<x<<' '<<y<<' '<<z<<'\n';
    }
    x=x0+0.0000001;
    y=y0+0.0000001;
    z=z0+0.0000001;
    for(int i=1;i<=N;i++)
    {
        k1[1]=functiexL(x,y);
        k1[2]=functieyL(x,y,z);
        k1[3]=functiezL(x,y,z);

        k2[1]=functiexL(x+h*k1[1]/2,y+h*k1[2]/2);
        k2[2]=functieyL(x+h*k1[1]/2,y+h*k1[2]/2,z+h*k1[3]/2);
        k2[3]=functiezL(x+h*k1[1]/2,y+h*k1[2]/2,z+h*k1[3]/2);

        k3[1]=functiexL(x+h*k2[1]/2,y+h*k2[2]/2);
        k3[2]=functieyL(x+h*k2[1]/2,y+h*k2[2]/2,z+h*k2[3]/2);
        k3[3]=functiezL(x+h*k2[1]/2,y+h*k2[2]/2,z+h*k2[3]/2);

        k4[1]=functiexL(x+h*k3[1],y+h*k3[2]);
        k4[2]=functieyL(x+h*k3[1],y+h*k3[2],z+h*k3[3]);
        k4[3]=functiezL(x+h*k3[1],y+h*k3[2],z+h*k3[3]);

        x=x+(h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        y=y+(h/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        z=z+(h/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);
        g<<x<<' '<<y<<' '<<z<<'\n';
    }
}
double functiexR(double y, double z)
{
    return -y-z;
}
double functieyR(double x, double y)
{
    return x+a*y;
}
double functiezR(double x, double z)
{
    return b+z*(x-c);
}
void RK4R(double x, double y, double z, double h, int N)
{
    double k1[4], k3[4], k2[4], k4[4];
    double x0=x,y0=y,z0=z;

    f.precision(15);

    for(int i=1;i<=N;i++)
    {
        k1[1]=functiexR(y,z);
        k1[2]=functieyR(x,y);
        k1[3]=functiezR(x,z);

        k2[1]=functiexR(y+h*k1[2]/2,z+h*k1[3]/2);
        k2[2]=functieyR(x+h*k1[1]/2,y+h*k1[2]/2);
        k2[3]=functiezR(x+h*k1[1]/2,z+h*k1[3]/2);

        k3[1]=functiexR(y+h*k2[2]/2,z+h*k2[3]/2);
        k3[2]=functieyR(x+h*k2[1]/2,y+h*k2[2]/2);
        k3[3]=functiezR(x+h*k2[1]/2,z+h*k2[3]/2);

        k4[1]=functiexR(y+h*k3[2],z+h*k3[3]);
        k4[2]=functieyR(x+h*k3[1],y+h*k3[2]);
        k4[3]=functiezR(x+h*k3[1],z+h*k3[3]);

        x=x+(h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        y=y+(h/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        z=z+(h/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);

        f<<x<<' '<<y<<' '<<z<<'\n';
    }
    x=x0+0.0000001;
    y=y0+0.0000001;
    z=z0+0.0000001;
    for(int i=1;i<=N;i++)
    {
        k1[1]=functiexR(y,z);
        k1[2]=functieyR(x,y);
        k1[3]=functiezR(x,z);

        k2[1]=functiexR(y+h*k1[2]/2,z+h*k1[3]/2);
        k2[2]=functieyR(x+h*k1[1]/2,y+h*k1[2]/2);
        k2[3]=functiezR(x+h*k1[1]/2,z+h*k1[3]/2);

        k3[1]=functiexR(y+h*k2[2]/2,z+h*k2[3]/2);
        k3[2]=functieyR(x+h*k2[1]/2,y+h*k2[2]/2);
        k3[3]=functiezR(x+h*k2[1]/2,z+h*k2[3]/2);

        k4[1]=functiexR(y+h*k3[2],z+h*k3[3]);
        k4[2]=functieyR(x+h*k3[1],y+h*k3[2]);
        k4[3]=functiezR(x+h*k3[1],z+h*k3[3]);

        x=x+(h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        y=y+(h/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        z=z+(h/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);

        g<<x<<' '<<y<<' '<<z<<'\n';
    }
}
