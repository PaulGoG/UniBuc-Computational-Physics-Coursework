#include<fstream>
#include<cmath>
#include<iostream>
#include<ctime>
#include<cstdlib>
using namespace std;

ifstream f("Icing.in");
ofstream g("Icing.out");



void citire_matrice(short **a, int n) //citeste configuratia initiala
{
	for (int i = 1; i<n+1; i++)
		for (int j = 1; j<n+1; j++)
			f >> a[i][j];
}
void afisare_matrice(short **a, int n) //afiseaza configuratia finala de energie minima
{
	for (int i = 1; i<n+1; i++)
	{
		for (int j = 1; j<n+1; j++)
			g << a[i][j] << ' ';
		g << '\n';
	}
}
int generator_pozitie(int m, int n) //genereaza o valoare aleatoare intreaga intre m si n
{
	return (int)(rand() % (n - m + 1) + m);
}
float generator_coeficient_agitatie()//genereala o valoare aleatoare reala intre 0 si 1
{
	return ((float)(rand())) / ((float)(RAND_MAX));
}
void bordare_matrice1(short **a, int n)// conditii periodice la limita
{
	for (int i = 1; i<n+1; i++)
	{
		a[i][0] = a[i][n];
		a[i][n+1] = a[i][1];
		a[0][i] = a[n][i];
		a[n+1][i] = a[1][i];
	}
}
void bordare_matrice0(short **a, int n)//bordare cu 0
{
	for (int i = 1; i<n+1; i++)
	{
		a[i][0] = 0;
		a[i][n+1] = 0;
		a[0][i] = 0;
		a[n+1][i] = 0;
	}
}
float deltaH(short **a, int n, int x, int y, int J) //calculeaza diferenta de energie intre 2 configuratii care difera prin termenul de coordonate x si y
{
	return 4 * a[x][y] * J*(a[x][y - 1] + a[x - 1][y] + a[x][y + 1] + a[x + 1][y]);
}
float H(short** a, int n, int J) //calculeaza energia configuratiei
{
	int S = 0;
	for (int i = 1; i<n+1; i++)
		for (int j = 1; j<n+1; j++)
			S += a[i][j] * (a[i][j - 1] + a[i - 1][j] + a[i][j + 1] + a[i + 1][j]);
	return (-1)*J*S;
}

int main()
{
	long double T;
	unsigned long long numarmare;
	float energie, coeficientAgitatie, step;
	int J;
	int x, y, ok;
	int n;
	f >> n;
	short** a;
	a = (short**)malloc((n+2) * sizeof(short*));
	for (int i = 0; i < n+2; i++) {
		a[i] = (short*)malloc(n * sizeof(short));
	}
	cout << "Introduceti 1 pentru conditii de periodicitate la limita si 0 pentru conditii de absorbtie (bordare cu 0-uri)" << '\n';
	cin >> ok;
	cout << "Introduceti un numar foarte mare care va determina de cate ori configuratia va suferi permutari" << '\n';
	cin >> numarmare;
	cout << "Introduceti rata de racire a temperaturii (termen supraunitar)" << '\n';
	cin >> step;
	cout << "Introduceti temperatura initiala (numar mare)" << '\n';
	cin >> T;
	cout << "Introduceti constanta J" << '\n';
	cin >> J;
	citire_matrice(a, n);
	srand(time(NULL));
	if (ok == 0) bordare_matrice0(a, n);
	else bordare_matrice1(a, n);
	cout << "Energia initiala a configuratiei este: " << H(a, n, J) << '\n';
	for (unsigned int i = 0; i <= numarmare; i++)
	{
		x = generator_pozitie(1, n );
		y = generator_pozitie(1, n );
		energie = deltaH(a, n, x, y, J);
		coeficientAgitatie = generator_coeficient_agitatie();
		if (exp(-energie / T)>coeficientAgitatie)
			a[x][y] = a[x][y] * (-1);
		T = T / step;
	}
	cout << "Energia minima a configuratiei este" << H(a, n, J) << '\n';

	afisare_matrice(a, n);
	for (int i = 0; i < n+2; i++) {
		free(a[i]);
	}
	free(a);
	return 0;
}
