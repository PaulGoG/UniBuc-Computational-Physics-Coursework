#include<iostream>
#include<fstream>
#include<ctime>
#include<cstdlib>
using namespace std;

//programul calculeaza energii pentru toate cele 3 metode in functie de N variabil si CI fixe in pt 2.
//si pt N fix si CI variabile calculeaza valori X,Y pentru reprezentare - spatiul fazelor pt 1.

ofstream g ("EnergieImplicit.txt");
ofstream h ("EnergieSIV1.txt");
ofstream p ("EnergieSIV2.txt");
ofstream j ("ConditiiInitialex.txt");
ofstream k ("ConditiiInitialey.txt");
ofstream m1 ("RezultatxImplicit.txt");
ofstream m2 ("RezultatyImplicit.txt");
ofstream n1 ("RezultatxSIV1.txt");
ofstream n2 ("RezultatySIV1.txt");
ofstream o1 ("RezultatxSIV2.txt");
ofstream o2 ("RezultatySIV2.txt");

void EulerImplicit(float &x, float &y, float dt, float ct, int N)
{
    float x1,y1;
    x1=x;
    y1=y;
    for(int i=1;i<=N;i++)
    {
        x=x1+dt*y1;
        y=y1+dt*(-ct*x1);
        x1=x;
        y1=y;
    }
}

void EulerImplicitEN(float &x, float &y, float dt, float ct, int N, float E)
{
    float x1,y1;
    x1=x;
    y1=y;
    for(int i=1;i<=N;i++)
    {
        E=0.5*y*y+0.5*ct*x*x;
        g<<E<<'\n';
        x=x1+dt*y1;
        y=y1+dt*(-ct*x1);
        x1=x;
        y1=y;
    }
     E=0.5*y*y+0.5*ct*x*x;
        g<<E<<'\n';
}


void EulerSemiImplicitV1(float &x, float &y, float dt, float ct, int N)
{
    float x1,y1;
    x1=x;
    y1=y;
    for(int i=1;i<=N;i++)
    {
        x=x1+dt*y1;
        y=y1+dt*(-ct*x);
        x1=x;
        y1=y;
    }
}

void EulerSemiImplicitV1EN(float &x, float &y, float dt, float ct, int N, float E)
{
    float x1,y1;
    x1=x;
    y1=y;
    for(int i=1;i<=N;i++)
    {
        E=0.5*y*y+0.5*ct*x*x;
        h<<E<<'\n';
        x=x1+dt*y1;
        y=y1+dt*(-ct*x);
        x1=x;
        y1=y;
    }
    E=0.5*y*y+0.5*ct*x*x;
        h<<E<<'\n';
}

void EulerSemiImplicitV2(float &x, float &y, float dt, float ct, int N)
{
    float x1,y1;
    x1=x;
    y1=y;
    for(int i=1;i<=N;i++)
    {
        y=y1+dt*(-ct*x1);
        x=x1+dt*y;
        x1=x;
        y1=y;
    }
}

void EulerSemiImplicitV2EN(float &x, float &y, float dt, float ct, int N, float E)
{
    float x1,y1;
    x1=x;
    y1=y;
    for(int i=1;i<=N;i++)
    {
        E=0.5*y*y+0.5*ct*x*x;
        p<<E<<'\n';
        y=y1+dt*(-ct*x1);
        x=x1+dt*y;
        x1=x;
        y1=y;
    }
     E=0.5*y*y+0.5*ct*x*x;
        p<<E<<'\n';

}


void generare(float &x, float&y, int s, int t)
{
    x=(rand() % (t - s + 1) + s);
    y=(rand() % (t - s + 1) + s);
}
int main()
{
    float x,y,x1,y1,ct,dt,E=0;
    int N,M,s,t;
        //  int Nmax;
        //  int t0,T;
        cout<<"Nr. de conditii initiale generate "<<'\n';
    cin>>M;
        cout<<"Nr. de pasi "<<'\n';
    cin>>N;
 //       cout<<"Nr. maxim de pasi"<<'\n';
 //  cin>>Nmax;
       // cin>>T; cin>>t0; dt=(float)(T-t0)/N;

    dt=(float)10/N; //rezolvam pt t intre 0 si 10
    ct=1;
    srand(time(NULL));

    s=0;    //intre ce valori intregi generam conditiile initiale
    t=100;

    for(int i=1;i<=M;i++)
    {
    generare(x1,y1,s,t);
    j<<x1<<'\n';
    k<<y1<<'\n';
    x=x1; y=y1;
        EulerImplicit(x,y,dt,ct,N);

    m1<<x<<'\n';
    m2<<y<<'\n';
    x=x1; y=y1;
        EulerSemiImplicitV1(x,y,dt,ct,N);

    n1<<x<<'\n';
    n2<<y<<'\n';
    x=x1; y=y1;
        EulerSemiImplicitV2(x,y,dt,ct,N);

    o1<<x<<'\n';
    o2<<y<<'\n';
    }

    generare(x1,y1,s,t);
    x=x1; y=y1;
    EulerImplicitEN(x,y,dt,ct,N,E);
    x=x1; y=y1;
    EulerSemiImplicitV1EN(x,y,dt,ct,N,E);
    x=x1; y=y1;
    EulerSemiImplicitV2EN(x,y,dt,ct,N,E);


/*    generare(x1,y1,s,t); // CI fixe de data asta si N variabil

 for(N=1;N<=Nmax;N++)
    {
        dt=(float)10/N;
        x=x1; y=y1;
            EulerImplicit(x,y,dt,ct,N);
        E=0.5*y*y+0.5*ct*x*x;
        g<<E<<'\n';
        x=x1; y=y1;
            EulerSemiImplicitV1(x,y,dt,ct,N);
        E=0.5*y*y+0.5*ct*x*x;
        h<<E<<'\n';
        x=x1; y=y1;
            EulerSemiImplicitV2(x,y,dt,ct,N);
        E=0.5*y*y+0.5*ct*x*x;
        p<<E<<'\n';

    }
*/ //practic am mai facut odata pb de la examen...
    return 0;
}
