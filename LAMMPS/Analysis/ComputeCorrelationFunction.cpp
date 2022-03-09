#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <cmath>
#include <limits>
#include <unistd.h>

using namespace std;

////// GLOBAL VARIABLES
int natoms, lseq, s0, i00, i11;
double timestep;
vector <double> coordinates[3];
vector <int> tipologia;
ifstream leggi;
ofstream scrivi;

////// DECLARATION OF FUNCTIONS
//DUMP READER
int ReadDump(int);

//COMPUTATIONS
void cnt(double *, int);
double corr_func();

//VECTORIAL OPERATIONS
void rotate(double *, double *, double, double*);
void cross_product(double *, double *, double *);
double dot_product(double *, double *, int);
void versor(double *, int, double*);
double angle_vec(double *, double *);
double angle_vec_sign(double *, double *, double*);
void cp_vec(double *, double *, int);
void diff_vec(double *, double *, double *, int);
void print_vec(double *vettore, int vecsize);

///////  MAIN
int main(int argc, char*argv[])
{

	if (argc != 5+1)
        {
                cout<<"Inputs:"<<endl;
                cout<<"1 --> trajectory"<<endl;
                cout<<"2 --> s0 (# of bps separating two points)"<<endl;
                cout<<"3 --> start nucleotide"<<endl;
                cout<<"4 --> stop nucleotide"<<endl;
                cout<<"5 --> output file"<<endl;
                return 0;
        }
	string trajectory(argv[1]);
        s0 = atoi(argv[2]);
	i00 = atoi(argv[3]);
	i11 = atoi(argv[4]);
	string corr_file(argv[5]);

	cout << setprecision(10);
	
	//Import first frame
	int flag;
	if (access(trajectory.c_str(), F_OK ) != -1 )
        {
                leggi.open(trajectory.c_str(), ifstream::in);
                flag = ReadDump(1);
        }
	else
	{
		cout<<"Trajectory file not found! Exiting..."<<endl;
		return 0;
	}

	double lstep = 0;
        int nstep = 0;

	scrivi.open(corr_file.c_str(), ios::out);

	//Analyze trajectory
	while (flag > 0)
        {
		scrivi<<timestep;
	        lstep += corr_func();
                nstep ++;
                flag = ReadDump(0);
        }
        leggi.close();
	scrivi.close();

	lstep *= 1./nstep;
        char nome[1000];
	sprintf(nome, "%s_ds", corr_file.c_str());
	scrivi.open(nome, ios::out);
	scrivi<<lstep<<endl;
	scrivi.close();

	return 1;
}



///////////////// FUNCTIONS 

int ReadDump(int build)
{
	string comstr;
        getline(leggi, comstr);
	double comodo;
	int indice, tipo;
        leggi>>comodo;
        if (!leggi.eof())
        {
                timestep = comodo;
                getline(leggi, comstr);
                getline(leggi, comstr);
                //number of atoms
                leggi>>natoms;
		if (build == 1)
		{
			lseq = (natoms+2)/6;
                        tipologia.resize(natoms+1);
                        for (int k1=0; k1<3; k1++)
                        {
                                coordinates[k1].resize(natoms+1);
                        }
		}
                getline(leggi, comstr);
                getline(leggi, comstr);
                //size of box
		for (int ii = 1; ii<=6; ii++)
                        leggi>>comodo;
                getline(leggi, comstr);
                getline(leggi, comstr);
                //coordinates
                for (int i=1; i<=natoms; i++)
                {
                        leggi>>indice;
                        leggi>>tipo;
                        tipologia[indice] = tipo;
			leggi>>comodo;
                        coordinates[0][indice] = 0.1*comodo;
                        leggi>>comodo;
                        coordinates[1][indice] = 0.1*comodo;
                        leggi>>comodo;
                        coordinates[2][indice] = 0.1*comodo;

                }
                getline(leggi, comstr);
                return 1;
        }
        else
                return 0;
}

void cnt(double *vettore, int nuc)
{
	int iA = 3*nuc - 2;
	int iB = natoms - iA;
	for (int k1=0; k1<3; k1++)
	    	vettore[k1] = 0.5*(coordinates[k1][iA] + coordinates[k1][iB]);
	return;
}

double corr_func()
{
	double vec1[3], vec2[3];
        double lb, distbb = 0;
        int ibb = 0;
        //Precomputation: tangent vector in all relevant points and lb
        double tgvec[i11+1][3];
        for (int ibp1 = i00; ibp1<=i11; ibp1++)
        {
		cnt(vec1, ibp1-s0/2);
		cnt(vec2, ibp1+s0/2);
		diff_vec(vec2, vec1, vec1, 3);
                lb = sqrt(dot_product(vec1, vec1, 3));
                for (int k1=0; k1<3; k1++)
                        tgvec[ibp1][k1] = vec1[k1]/lb;
                distbb += lb;
                ibb++;
        }
        distbb *= 1./ibb;

        //Compute correlation
        double coseno;
	int i_coseno;
        for (int ds = 0; ds <= i11-i00; ds += s0)
	{
		coseno = 0;
		i_coseno = 0;
        	for (int ibp1 = i00; ibp1<=i11-ds; ibp1++)
	        {
        	        for (int k1=0; k1<3; k1++)
			{
                	        vec1[k1] = tgvec[ibp1][k1];
                	        vec2[k1] = tgvec[ibp1+ds][k1];
			}

                        coseno += dot_product(vec1, vec2, 3);
                        i_coseno ++;
                }
		scrivi<<" "<<coseno/i_coseno;
        }
	scrivi<<endl;

        return distbb;
}

void rotate(double *vettore, double *axis, double phi, double *vettoredest)
{
    //Use Rodrigues' formula
    double pscal = dot_product(axis, vettore, 3);
    double pvec[3];
    cross_product(axis, vettore, pvec);
    for (int k1=0; k1<3; k1++)
	    vettoredest[k1] = vettore[k1]*cos(phi) + pvec[k1]*sin(phi) + axis[k1]*pscal*(1-cos(phi));
    return;
}

void cross_product(double *vettore1, double *vettore2, double *vettoredest)
{
	vettoredest[0] = vettore1[1]*vettore2[2] - vettore1[2]*vettore2[1];
	vettoredest[1] = vettore1[2]*vettore2[0] - vettore1[0]*vettore2[2];
	vettoredest[2] = vettore1[0]*vettore2[1] - vettore1[1]*vettore2[0];
	return;
}

double dot_product(double *vettore1, double *vettore2, int sizevec)
{
	double ss = 0;
	for (int k1=0; k1<sizevec; k1++)
		ss += vettore1[k1]*vettore2[k1];
	return ss;
}

void versor(double *vettore, int sizevec, double *vettoredest)
{
	double nn = sqrt(dot_product(vettore, vettore, sizevec));
	for (int k1=0; k1<sizevec; k1++)
		vettoredest[k1] *= 1./nn;
}

double angle_vec(double *v1, double *v2)
{
    versor(v1,3,v1);
    versor(v2,3,v2);
    double sp = dot_product(v1,v2,3);
    if (sp > 1.0)
	    sp = 1.0;
    if (sp < -1.0)
	    sp = -1.0;
    return acos(sp);
}

double angle_vec_sign(double *v1, double *v2, double *asse)
{
	versor(v1,3,v1);
	versor(v2,3,v2);
	double theta = angle_vec(v1, v2);
	double cp[3];
	cross_product(v1, v2, cp);
	double checksegno = dot_product(cp, asse, 3);
	if (checksegno <0)
		theta *= -1;
	return theta;
}

void cp_vec(double *vettore_origine, double *vettore_destinazione, int sizevec)
{
	for (int k1=0; k1<sizevec; k1++)
		vettore_destinazione[k1] = vettore_origine[k1];
	return;
}

//destvec = vecA - vecB
void diff_vec(double *vecA, double *vecB, double *destvec, int sizevec)
{
	for (int k1=0; k1<sizevec; k1++)
		destvec[k1] = vecA[k1] - vecB[k1];
	return;
}

void print_vec(double *vettore, int vecsize)
{
        for (int i=0; i<vecsize; i++)
                cout<<vettore[i]<<" ";
        cout<<endl;
}
