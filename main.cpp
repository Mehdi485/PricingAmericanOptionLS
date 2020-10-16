#include "Header.hpp"

#define S0 36
#define STRIKE 110
#define N 1000
#define M 50
#define MATURITE 1
#define VOL 0.20
#define RATE 0.06
#define BASE 3 //nb de polynome pour la regression

using namespace std;

void info_prog();

int main()
{
	srand(time(NULL));
	
	double variance=0;
	double MC_erreur=0;
	clock_t start,inter, end;
	start = clock();
	info_prog();

	Matrix Spot(N,N);
	Spot.Mise_A_Zero();

	GenerSpot P(S0,MATURITE,VOL,RATE,N,M,Spot);
	
	inter = clock()-start;	
	cout << "le programme a mis " << " " << (((float)inter)/CLOCKS_PER_SEC) << " " << "pour creer la matrice des prix " << endl;
	cout << "Le programme va maintenir calculer le prix de l'option américaine, cela peut prendre un peu de temps" << endl;	
	
	cout << "En attendant que le programme s'exécute vous pouvez jetez un coup d'oeil au rapport de ce projet" << endl; cout << endl;
	LSmethod LS(BASE,STRIKE,P);
	Matrix EF(N,M );
	EF.Mise_A_Zero();
	EF = LS.Exer_mat(EF);
	double valuation = LS.Americain(EF, variance, MC_erreur);
	double europeen = LS.Europeen();
	end = clock() - start;
	cout << "Merci de votre patience " << endl; cout << endl;
	cout << "La valeur du put americain est " << " " << valuation << " de plus la variance est de " << variance << " et l'erreur de Monte Carlo est " << MC_erreur <<  endl;
	cout << "La valeur du put europeen est"	<< " " << europeen << endl; cout << endl;
	cout << "Le programme a mis " << (((float)end)/CLOCKS_PER_SEC) << "secondes pour évaluer le put Américain et le Put Européen " <<  endl; 
	return 0;

};





void info_prog()
{
	cout << "Pour rappel, voilà les constantes entrées dans le programme" << endl;

	cout << "le Prix de départ, la racine pour la simulation des prix de l'actif est " << " " << S0 << endl;
		
	cout << "le Strike du contrat est de " << " " << STRIKE << endl;
	
	cout << "la voloatilité entrée dans le programme est " << " " << VOL << endl;
	
	cout << "la maturité de l'actif entrée dans le programme est " << " " << MATURITE << endl;
	
	cout << "le taux sans risque entré dans le programme est " << " " << RATE << endl;
	
	cout << "le nombre de chemin simulé est de " << " " << N << endl;
	
	cout << "Avec le nombre de temps d'arret observé par l'actif suivant " << " " << M << endl;
	
	cout << "Voilà un résumé des conditions pour l'exécution de ce programme d'évalutation d'un Put Américain avec un seul sou-jacent" << endl; 
	
	cout << "le programme vient de débuter, selon les conditions il peut prendre plus ou moins de temps" << endl;
	cout << endl;
};
