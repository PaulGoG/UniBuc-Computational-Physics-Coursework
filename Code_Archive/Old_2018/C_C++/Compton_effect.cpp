#include<iostream>
#include<cmath>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<iomanip>
using namespace std;

#define PI 3.14159265
ofstream g ("Compton.txt");
ofstream h ("Tabel.txt");
double generareUnghi ()
        {
        return ((double)(rand()))/((double)(RAND_MAX))*180;
        }
double calculRecul (double x, double E, double E0)
        {
            return atan(((E*sin(x*PI/180))/(E0-E*cos(x*PI/180))))*180/PI;
        }
double calculE1 (double x, double E0)
        {
            return E0/(1+(E0/511)*(1-cos(x*PI/180)));
        }
double calculE2 (double x, double y, double E0)
        {
            return E0*(sin(y*PI/180)/sin((x+y)*PI/180));
        }
int main()
{
    double theta,phi,E1,E2,E,E0,Eprag,delta=0;
    int numara=0;
    cout<<"Introduceti energia initiala a fotonului incident in KeV: ";
    cin>>E0;
    cout<<"Introduceti energia de prag la care inceteaza interactiile de tip Compton in KeV: ";
    cin>>Eprag;
    cout<<'\n';
    if (E0<=Eprag) { cout<<"Energia initiala a fotonului nu este suficienta pentru a manifesta o interactie de tip Compton!";
                        g<<"Energia initiala a fotonului nu este suficienta pentru a manifesta o interactie de tip Compton!"; }
    else {
            srand((int)time(NULL));
                   h<<"-----------------------------------------------------------------------"<<'\n';
                   h<<"| Nr.ciocnire | Theta |  Phi  | Gamma |     E     |  Delta  | Abatere |"<<'\n';
                   h<<"|             |(grade)|(grade)|(grade)|   (KeV)   |  (KeV)  |   (%)   |"<<'\n';
                   h<<"-----------------------------------------------------------------------"<<'\n';
            while(E0>Eprag)
            {


                cout<<"Nr. ciocnire: "<<numara+1<<'.'<<'\n';
                   g<<"Nr. ciocnire: "<<numara+1<<'.'<<'\n';
                theta=generareUnghi(); { cout<<"Unghiul de deplasare a directiei fotonului este "<<(int)(theta)<<'.'<<(int)((theta)*100)%100<<" grade."<<'\n';
                                            g<<"Unghiul de deplasare a directiei fotonului este "<<(int)(theta)<<'.'<<(int)((theta)*100)%100<<" grade."<<'\n'; }
                E1= calculE1(theta,E0);
                phi= calculRecul(theta,E1,E0); { cout<<"Unghiul de deplasare a directiei electronului de recul este "<<(int)(phi)<<'.'<<(int)((phi)*100)%100<<" grade."<<'\n';
                                                    g<<"Unghiul de deplasare a directiei electronului de recul este "<<(int)(phi)<<'.'<<(int)((phi)*100)%100<<" grade."<<'\n'; }
                                               { cout<<"Unghiul dintre directiile fotonului si electronului de recul dupa ciocnire este "<<(int)(theta+phi)<<'.'<<(int)((theta+phi)*100)%100<<" grade."<<'\n';
                                                    g<<"Unghiul dintre directiile fotonului si electronului de recul dupa ciocnire este "<<(int)(theta+phi)<<'.'<<(int)((theta+phi)*100)%100<<" grade."<<'\n'; }
                E2= calculE2(theta,phi,E0);
                E=(E1+E2)/2; { cout<<"Energia fotonului dupa interactie este "<<(int)(E)<<'.'<<(int)((E)*100)%100<<" KeV."<<'\n';
                                  g<<"Energia fotonului dupa interactie este "<<(int)(E)<<'.'<<(int)((E)*100)%100<<" KeV."<<'\n'; }
                delta=delta+(E0-E);
                numara++;
                     { cout<<"Abaterea standard este: "<<(abs)((int)(100-(((abs)((E0-E)-(delta/numara)))/(delta/numara))*100)%100)<<'.'<<(abs)((int)(10000-(((abs)((E0-E)-(delta/numara)))/(delta/numara))*10000)%100)<<'%'<<'\n';
                          g<<"Abaterea standard este: "<<(abs)((int)(100-(((abs)((E0-E)-(delta/numara)))/(delta/numara))*100)%100)<<'.'<<(abs)((int)(10000-(((abs)((E0-E)-(delta/numara)))/(delta/numara))*10000)%100)<<'%'<<'\n'; }

                cout<<'\n'<<'\n';
                   g<<'\n'<<'\n';
                   h<<"|     "<<std::left<<setw(3)<<numara<<"     | "<<std::left<<setw(3)<<(int)(theta)<<'.'<<std::left<<setw(2)<<(int)((theta)*100)%100<<"| "<<std::left<<setw(3)<<(int)(phi)<<'.'<<std::left<<setw(2)<<(int)((phi)*100)%100<<"| "<<std::left<<setw(3)<<(int)(theta+phi)<<'.'<<std::left<<setw(2)<<(int)((theta+phi)*100)%100<<"|  "<<std::left<<setw(4)<<(int)(E)<<'.'<<std::left<<setw(2)<<(int)((E)*100)%100<<"  | "<<std::left<<setw(4)<<(int)(delta)<<'.'<<std::left<<setw(2)<<(int)((delta)*100)%100<<" |  "<<std::left<<setw(2)<<(abs)((int)(100-(((abs)((E0-E)-(delta/numara)))/(delta/numara))*100)%100)<<'.'<<std::left<<setw(2)<<(abs)((int)(10000-(((abs)((E0-E)-(delta/numara)))/(delta/numara))*10000)%100)<<"  |"<<'\n';
                    E0=E;
            }
            h<<"-----------------------------------------------------------------------"<<'\n';
            if(numara<10)
            { cout<<"S-au efectuat "<<numara<<" ciocniri pe parcursul simularii, iar diferenta medie intre energiile initiale si finale ale fiecarei ciocniri este "<<(int)(delta/numara)<<'.'<<(int)((delta/numara)*100)%100<<" KeV";
                 g<<"S-au efectuat "<<numara<<" ciocniri pe parcursul simularii, iar diferenta medie intre energiile initiale si finale ale fiecarei ciocniri este "<<(int)(delta/numara)<<'.'<<(int)((delta/numara)*100)%100<<" KeV"; }
            else if(numara>9)
            { cout<<"S-au efectuat "<<numara<<" de ciocniri pe parcursul simularii, iar diferenta medie intre energiile initiale si finale ale fiecarei ciocniri este "<<(int)(delta/numara)<<'.'<<(int)((delta/numara)*100)%100<<" KeV";
                 g<<"S-au efectuat "<<numara<<" de ciocniri pe parcursul simularii, iar diferenta medie intre energiile initiale si finale ale fiecarei ciocniri este "<<(int)(delta/numara)<<'.'<<(int)((delta/numara)*100)%100<<" KeV";}
         }
         return 0;
}
