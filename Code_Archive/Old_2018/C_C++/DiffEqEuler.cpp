#include<iostream>
#include<cmath>
using namespace std;
double function (double y)
{
    return y+sin(0.4*y);
}
int main()
{
    double x,y,step;
    int N;
    y=1;
    x=1; // xinitial=o
    //cout<<"Introduceti valoarea argumentului pentru care aflam functia"<<'\n';
    //cin>>x;
    cout<<"Introduceti numarul de diviziuni (N) ";
    cin>>N;
    cout<<'\n';
    step=x/N;
    for(int i=1;i<=N;i++)
    {
        y=y+function(y)*step;
    }
    cout<<"Valoarea functiei este: "<<y<<'\n';
    return 0;
    }
