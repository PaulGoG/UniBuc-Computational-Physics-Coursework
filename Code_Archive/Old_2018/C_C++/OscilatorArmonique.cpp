#include<iostream>
#include<fstream>
#include<cmath>

#define omega 1

using namespace std;

ofstream f ("OscArm1.txt");
ofstream g ("OscArm2.txt");
ifstream m ("OscArm1.txt");
ifstream n ("OscArm2.txt");
ofstream d ("DistantaArm.txt");
ifstream p ("DistantaArm.txt");
ofstream timp ("TimpiUCArm.txt");


double functiex(double y)
{
    return y;
}
double functiey(double x)
{
    return -omega*omega*x;
}
void RK41(double x, double y,double h, int N)
{
    double k1[3], k3[3], k2[3], k4[3];

    for(int i=1;i<=N;i++)
    {
        k1[1]=functiex(y);
        k1[2]=functiey(x);

        k2[1]=functiex(y+h*k1[2]/2);
        k2[2]=functiey(x+h*k1[1]/2);

        k3[1]=functiex(y+h*k2[2]/2);
        k3[2]=functiey(x+h*k2[1]/2);

        k4[1]=functiex(y+h*k3[2]);
        k4[2]=functiey(x+h*k3[1]);

        x=x+(h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        y=y+(h/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        f<<x<<' '<<y<<'\n';
    }
}
void RK42(double x, double y,double h, int N)
{
    double k1[3], k3[3], k2[3], k4[3];

    for(int i=1;i<=N;i++)
    {
        k1[1]=functiex(y);
        k1[2]=functiey(x);

        k2[1]=functiex(y+h*k1[2]/2);
        k2[2]=functiey(x+h*k1[1]/2);

        k3[1]=functiex(y+h*k2[2]/2);
        k3[2]=functiey(x+h*k2[1]/2);

        k4[1]=functiex(y+h*k3[2]);
        k4[2]=functiey(x+h*k3[1]);

        x=x+(h/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        y=y+(h/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        g<<x<<' '<<y<<'\n';
    }
}
void Distanta(double h, int N)
{
    double x1,y1,x2,y2,t=0;
    for(int i=1;i<=N;i++)
    {
        m>>x1>>y1;
        n>>x2>>y2;
        d<<sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))<<' '<<t<<'\n';
        t+=h;
    }
}
void extremelocale(double h, int N)
{
    double T,D1,D2,D3,De,te=0,epsilon=0.0009;
    short low,high;
    p>>D1>>T;
    p>>D2>>T;
    p>>D3>>T;
    N-=3;
    De=D1;
    low=3;
    high=2;

    // Skip primul minim

    do
    {
        N--;
        if(D2>D1) low=1;
        if(D2<D1) low=0;
        if(D2>D3) high=1;
        if(D2<D3) high=0;
        if(low==high && abs(De-D2)>epsilon)
        {
            De=D2;
            if(low==1) timp<<T-te-h<<' ';
            te=T-h;
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
        if(low==high && abs(De-D2)>epsilon)
        {
            timp<<T-te-h;
            te=T-h;
            De=D2;
            if(low==1) timp<<' ';
                else timp<<'\n';
        }
        if(D2!=D3) D1=D2;
        D2=D3;
        low=2;
        high=3;
        p>>D3;
        p>>T;
    }
}
int main()
{
    double x,y,h,tMAX;
    int N;

    cout<<"Variabile: N-nr. pasi, x0,y0-primul set de conditii initiale,"<<'\n'<<"tMAX-intervalul de timpi (de la 0 la ceva): "<<'\n'<<'\n'<<"    ";

    cin>>N;
    cin>>x>>y;
    cin>>tMAX;
    h=(double)tMAX/N;
    RK41(x,y,h,N);
    cout<<'\n';
    cout<<"x0',y0'- al doilea set de conditii initiale ";
    cout<<'\n'<<'\n'<<"   ";
    cin>>x>>y;
    RK42(x,y,h,N);
    Distanta(h,N);
    extremelocale(h,N);
    return 0;
}
