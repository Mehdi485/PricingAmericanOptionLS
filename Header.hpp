#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<vector>

using namespace std;

class Matrix {
	public:
	
		/*constructeur, destructeur, constructeur par copie*/
		Matrix(int n, int m);
		Matrix (const Matrix & A);
		~Matrix();
		
		/*operator + methode afficher*/
		Matrix & operator=(const Matrix & M);
		Matrix & operator=(const vector<double> vec);
		Matrix operator*(const Matrix & M);
		vector<double> operator*(const vector<double> vec);
		Matrix & operator*=(const double a);
		double operator()(const int i, const int j)const {if(i>nblig || j>nbcol){return EXIT_FAILURE;} return v[i][j];};
		double & operator()(int i, int j) {return v[i][j];};
		friend ostream & operator<<(ostream & o, const Matrix & A);
		void afficher(int i, int j){cout << v[i][j];}
		
		/*methodes*/
		void Mise_A_Zero();
		Matrix Cholesky(const Matrix & M) const;
		Matrix & inverse(); // calcul de l'inverse par le pivot de Gauss 
		Matrix operator/(const double nb);
		Matrix transpose();

		/*accesseurs*/
		int get_nb_lignes(){return nblig;};
		int get_nb_colonnes(){return nbcol;};
		int get_nb_lignes()const {return nblig;};
		int get_nb_colonnes()const {return nbcol;};
	private:
		int nblig,nbcol;
		vector < vector<double> > v;
};


double unif_rand();
double norm_rand();
double LP(int R, double x);//Polynome de Laguerre avec un poids de exp(-x/2)
void info_prog();

class GenerSpot{
	public:
		GenerSpot(double S0, double T, double vol, double r, int nbsimilations, int nbstep, Matrix & npath);
	
		double S0, T, vol, r, strike; int nbsimulations, nbsteps;
		Matrix & npath;
};

class LSmethod{
	public:
		LSmethod(int nbBase, double K, GenerSpot & S):nbBase(nbBase),K(K), S(S){};
		Matrix regress(vector<double> X, vector<double> Y);
		void Exercice(Matrix & EF, const Matrix & CE, const Matrix & CC, int tps);
		Matrix BETA(const Matrix & CE,const Matrix & CC, int tps);
		void payoff(Matrix & CE, int tps);
		Matrix Exer_mat(Matrix & EF);
		double Americain(Matrix & EF, double & variance, double & MC_erreur);
		double Europeen();		
	private:
		int nbBase;
		GenerSpot & S;
		double K;
};
