#include "Header.hpp"

using namespace std;




double unif_rand() {return (rand()+0.5)/(RAND_MAX + 1.0);};

double norm_rand()
{
	double U=unif_rand();
	double V=2*M_PI*unif_rand();
	return (sqrt(-2*M_PI*log(U))*cos(V));
};

double LP(int R, double x)
{
	double LaguerrePolynome = 0;
	if(R==0) LaguerrePolynome = 1;
	else if (R==1) LaguerrePolynome = (1-x);
	else if (R==2) LaguerrePolynome = (1-2*x+0.5*x*x);
	return (LaguerrePolynome*exp(-x/2));
}

double max(double x, double y)
{
	if(x>y) return x;
	else return y;
}





GenerSpot::GenerSpot(double S0, double T, double vol, double r,int nbsimulations, int nbsteps, Matrix & npath)
			:S0(S0),T(T),vol(vol), r(r), nbsimulations(nbsimulations), nbsteps(nbsteps), npath(npath)
{
	nbsteps=nbsteps*T;
	double dt=1/(nbsteps);
	double drift = (r-0.5*vol*vol)*dt;
	double vsqrdt = vol*sqrt(dt);
	for(int i=0;i<nbsimulations;++i)
	{
		npath(i,0)=S0;
		for(int j=1;j<nbsteps;++j)
		{
			npath(i,j)= npath(i,j-1)*exp(drift + vsqrdt*norm_rand());
		}	
	
	}
};


Matrix LSmethod::regress(vector<double> X, vector<double> Y)
{
	Matrix beta(nbBase,1);
	beta.Mise_A_Zero();
	if(X.size()!=0)
	{	
		int dim = X.size();
		Matrix reg(dim, nbBase);
		Matrix t_reg(nbBase,dim);
		Matrix tregxreg(nbBase,nbBase);
		reg.Mise_A_Zero();
		t_reg.Mise_A_Zero();
		tregxreg(nbBase,nbBase);
		for(int i=0;i<dim;i++)
		{
			for(int r=0;r<nbBase;r++)
			{
				reg(i,r)= LP(nbBase-1-r, X[i]);		
			}
		}
		
		t_reg = reg.transpose();
		tregxreg = t_reg*reg;
		tregxreg.inverse();
		Y = t_reg*Y;
		beta = tregxreg * Y;
	}  
	return beta;
};


void LSmethod::payoff(Matrix & CE, int tps)
{

	for(int i=0;i<CE.get_nb_lignes();i++)
	{
		CE(i,tps) = max(K-S.npath(i,tps),0);
	}
};


Matrix LSmethod::BETA(const Matrix & CE, const Matrix & CC, int tps)
{
	Matrix beta(nbBase,1);
	
	vector<double> X; 
	vector<double> Y;
	double discount = exp((-S.r));
	
	for(int i = 0;i<CE.get_nb_lignes();i++)
	{
		if(CE(i,tps)>0)
		{
			X.push_back(S.npath(i,tps));
			Y.push_back(CC(i,tps+1)*discount);
		}
		
	}
	beta = regress(X,Y);
	return beta;
};

void LSmethod::Exercice(Matrix & EF, const Matrix & CE, const Matrix & CC, int tps)
{
	for(int i=0;i<(*this).S.nbsimulations;i++)
	{
		if(CE(i,tps)>CC(i,tps)) 
		{
			EF(i,tps) =1;
			for(int t = tps+1;t<(*this).S.nbsteps;t++) {EF(i,t)=0;}
		}
		else 
			EF(i,tps)=0;
	}
};

Matrix LSmethod::Exer_mat(Matrix & EF)
{
	int chemin = (*this).S.nbsimulations;
	int maturite = (*this).S.nbsteps;

	Matrix CC(chemin, maturite);
	Matrix CE(chemin, maturite);
	Matrix Beta(nbBase,1);
	
	CC.Mise_A_Zero();
	CE.Mise_A_Zero();
	Beta.Mise_A_Zero();
	double delta_t=((*this).S.T)/maturite;		
	double discount = exp(-(S.r)*delta_t);	
	
	for(int i=0;i<chemin; i++)
	{
		payoff(CE, maturite-1);
		CC(i,maturite-1) = CE(i,maturite-1);
		if (CE(i,maturite-1) >0) 
			EF(i,maturite-1) = 1;
		else EF(i,maturite-1) =0;
	}

	
	for(int tps = maturite-2; tps>=1; tps--)
	{
		payoff(CE,tps);
		
		for(int i=0;i<chemin;i++)
		{
			Beta = BETA(CE,CC,tps); 
			if(CE(i,tps)>0)
			{
				CC(i,tps) = Beta(0,0)*LP(2,S.npath(i,tps))+Beta(1,0)*LP(1,S.npath(i,tps))+Beta(2,0)*LP(0,S.npath(i,tps));
			}
		}
	
		Exercice(EF,CE, CC, tps);
	
		
		for(int i = 0;i<chemin;i++)
		{
			if(EF(i,tps)==0) CC(i,tps) = CC(i,tps)*discount;
			else if (EF(i,tps)==1) CC(i,tps) = CE(i,tps);
		}	
	}
	return EF;
};



double LSmethod::Americain(Matrix & EF, double & variance, double & MC_erreur)
{
	//cout << EF << endl;
	double value=0;
	double delta_t = S.T/S.nbsteps;
	vector<double> v;
	for(int i=1;i<(*this).S.nbsimulations;i++)
	{
		for(int j =0; j< EF.get_nb_colonnes();j++)
		{
			if(EF(i,j)==1)
				v.push_back(exp((-S.r)*(j)*delta_t)*S.npath(i,j));
		}
	}
	for (int i=0;i<v.size();i++)
	{
		value += v[i];	
	}
	value = value/ S.nbsimulations;
	for(int i=0;i<v.size();i++)
	{
		variance += pow(v[i]-value,2);
	}
	variance= variance/S.nbsimulations;
	MC_erreur = sqrt(variance/S.nbsimulations);
	return value;
};

double LSmethod::Europeen()
{
	double value=0;
	double delta_t=S.T/S.nbsteps;
	double discount = exp((-S.r)*(S.nbsteps)*delta_t);
	for(int i=0;i<S.nbsimulations;i++)
	{
		double payoff = K-S.npath(i,S.nbsteps-1);
		value = value + discount*(max( payoff ,0));
	}
	value = value /S.nbsimulations;
	return value;	
};


Matrix::Matrix(int n, int m): nblig(n),nbcol(m), v(n,vector<double>(m)){};

Matrix::~Matrix()
{
	//for(int i=0;i<nblig;i++) {v[i].clear();}
	v.clear();
};

void Matrix::Mise_A_Zero()
{
	for(int i=0;i<nblig;i++)
	{
		for(int j=0;j<nbcol;j++)
		{
			v[i][j]=0;
		}
	}
	
};

Matrix::Matrix(const Matrix & A)
{
	nblig = A.nblig;
	nbcol = A.nbcol;
	v.resize(nblig,vector<double>(nbcol));
	
	for(int i=0;i<nblig;i++)
	{
		for(int j=0;j<nbcol;j++)
		{
			v[i][j]=A.v[i][j];
		}
	}
};

Matrix Matrix::operator*(const Matrix & M)
{
	
	Matrix PROD(nblig, M.nbcol);
	for(int i=0;i<nblig;i++)
	{
		for(int j=0;j<M.nbcol;j++)
		{	
			for(int k=0;k<nbcol;k++)
			{
				PROD(i,j) = PROD(i,j) + (*this)(i,k)*M(k,j);
			}
		}
	}
	
	return PROD;
};

vector<double> Matrix::operator*(const vector<double> vec)
{
	vector<double> pro_vec(nblig,0);
	for(int i=0;i<nblig;i++)
	{
		for(int j=0;j<nbcol;j++)
		{
			pro_vec[i] = pro_vec[i] + ((*this)(i,j) * vec[j]);
		}
	}
	return pro_vec;
};

Matrix Matrix::transpose()
{

	Matrix TRANS(nbcol, nblig);
	for(int i=0;i<nbcol;i++)
	{
		for(int j=0; j<nblig;j++)
		{
			TRANS(i,j)=(*this)(j,i);
		}
	}

	return TRANS;
};

Matrix & Matrix::operator=(const Matrix & M)
{
	if(M.nblig != nblig || M.nbcol != nbcol)
	{
		v.resize(M.nblig, vector<double>(M.nbcol));
      	} 	
	for(int i=0;i<M.nblig;i++)
	{
		for(int j=0;j<M.nbcol;j++)
		{
			(*this)(i,j)=M(i,j);
		}
	}
	return (*this);
};

Matrix & Matrix::operator=(const vector<double> vec)
{
	if(nbcol!=1){cout << "erreur, revoir le code, la matrix a plus d'une colonne" << endl;}
	if(nblig!= vec.size() || nbcol!= 1)
	{
		//for(unsigned int i=0; i<nblig; i++) {v[i].clear();}
         	v.clear();
        	nblig = vec.size();
         	nbcol = 1;
         
         	v.resize(nblig,vector<double>(nbcol));
		//for(int i=0;i<nblig;i++){ v[i].resize(nbcol);}
      	} 	
	for(int i=0;i<nblig;i++)
	{
		(*this)(i,0)= vec[i];
	}

	return *this;
};

Matrix & Matrix::operator*=(const double a)
{
	for(int i=0;i<nblig;i++)
	{
		for(int j=0;j<nbcol;j++)
		{
			(*this)(i,j)*=a;
		}	
	}

	return *this;
};

Matrix Matrix::operator/(const double nb) 
{
   	if(nb!=0) 
   	{
      		for(int i=0; i<nblig; i++)
         		for(int j=0; j<nbcol; j++){(*this)(i,j)/=nb;}
   	}
  	else //exception division par 0
   	{
      		cout << "erreur division avec 0, code à revoir!!!" << endl;
   	}
   	return *this;
};

ostream & operator<<(ostream &o,const Matrix &M)
{
  	for(unsigned int i=0; i< M.nblig; i++)
   	{
      		for(unsigned int j=0; j<M.nbcol; j++)
        		 o << M(i,j) << " ";
      		o << endl;
   	}
   	return o;
};



Matrix Matrix::Cholesky(const Matrix & M)const 
{
	Matrix L(M.nblig,M.nbcol);
	double taille = M.nblig;
	double s=0;
	for(int j=0;j<taille;j++)
	{	
		s=0;
		for(int k=0;k<j;k++){s+=pow(M(j,k),2);}
		L(j,j)=sqrt(M(j,j)-s);
		for(int i=j;i<taille;i++)
		{
			s=0;
			for(int k=0;k<j;k++){s+=L(i,k)*L(j,k);}
			L(i,j) = (M(i,j)-s) / L(j,j); 			
		}
	}
	return L;
};


Matrix & Matrix::inverse()
{
	if(nblig != nbcol)
	{
		cout << "inversion imposssible d'une matrice nn carré" << endl;
		return *this;	
	}
	int n = nblig;
	Matrix MM(n,n);
	Matrix M1(n,n);
	Matrix Mi1(n,n);
	Matrix Minv(n,n);

	MM = (*this);
	M1 = MM;
	Mi1.Mise_A_Zero();
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			if(i==j)Mi1(i,j)=1;

	Minv = Mi1;
	for (int i=0;i<n;i++)
        {
             	for (int j=0;j<n;j++)
                {
                 	if (MM(i,i)==0)
                    	{
                    		cout << "inversion pivot GAUSS impossible : division par 0" << endl;
                    		return *this;
                    	}
                 	M1(i,j)=MM(i,j)/MM(i,i);
                 	Mi1(i,j)=Minv(i,j)/(MM(i,i));
                 }
             	 MM=M1;
            	 Minv=Mi1;
             	 for (int k=0;k<n;k++)// mise zéros de la colonne
             	 if (k!=i)
             	 for (int j=0;j<n;j++)
                 {
                 	M1(k,j)=MM(k,j)-MM(i,j)*MM(k,i);
                 	Mi1(k,j)=Minv(k,j)-Minv(i,j)*MM(k,i);
                 }
             	MM=M1;
             	Minv=Mi1;
      	 }
         *this=Minv;
         return *this;
};

