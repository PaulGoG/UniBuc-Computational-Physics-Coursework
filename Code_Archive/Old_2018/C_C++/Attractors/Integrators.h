#ifndef INTEGRATORI_H
#define INTEGRATORI_H

double functiexL(double x, double y);
double functieyL(double x, double y, double z);
double functieyL(double x, double y, double z);
void RK4L(double x, double y, double z, double h, int N);

double functiexR(double y, double z);
double functieyR(double x, double y);
double functieyR(double x, double z);
void RK4R(double x, double y, double z, double h, int N);

#endif
