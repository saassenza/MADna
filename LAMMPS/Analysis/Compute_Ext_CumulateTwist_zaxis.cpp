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
int natoms, lseq;
double timestep, haxis[3];
vector <double> coordinates[3];
vector <int> tipologia;
ifstream leggi;

////// DECLARATION OF FUNCTIONS
//DUMP READER
int ReadDump(double, double, int);

//COMPUTATIONS
void cnt_bases(double *, int);
void cnt(double *, int);
void bpvec(double *, int, int);
double compute_htwist(int);

//VECTORIAL OPERATIONS
void cross_product(double *, double *, double *);
double dot_product(double *, double *, int);
void versor(double *, int, double*);
double angle_vec(double *, double *);
double angle_vec_sign(double *, double *, double*);

///////  MAIN
int main(int argc, char*argv[])
{

	if (argc != 3+1)
	{
		cout<<"Inputs:"<<endl;
		cout<<"1 --> trajectory"<<endl;
		cout<<"2 --> start nucleotide"<<endl;
		cout<<"3 --> stop nucleotide"<<endl;
		cout<<"Output: timestep extension(nm) cumulate-twist(rad)"<<endl;
		return 0;
	}
	string trajectory(argv[1]);
	int i00 = atoi(argv[2]);
	int i11 = atoi(argv[3]);

	cout << setprecision(10);
	
	//Import first frame
	int flag;
	if (access(trajectory.c_str(), F_OK ) != -1 )
        {
                leggi.open(trajectory.c_str(), ifstream::in);
                flag = ReadDump(-1, 1e10, 1);
        }
	else
	{
		cout<<"Trajectory file not found! Exiting..."<<endl;
		return 0;
	}

	double cum_htwist, ext, vec1[3], vec2[3];

	haxis[0] = 0;
	haxis[1] = 0;
	haxis[2] = 1;

	//Analyze trajectory
	while (flag > 0)
        {
		//compute extension
		cnt_bases(vec1, i00);
		cnt_bases(vec2, i11);
		ext = 0.1*(vec2[2] - vec1[2]);

		//compute cumulate twist
		cum_htwist = 0;
                for (int iseq = i00; iseq < i11; iseq++)
                {
                        cum_htwist += compute_htwist(iseq);
                }

                cout<<timestep<<" "<<ext<<" "<<cum_htwist<<endl;
                flag = ReadDump(-1, 1e10, 0);
        }
        leggi.close();

	return 1;
}



///////////////// FUNCTIONS 

int ReadDump(double t00, double t11, int build)
{
	string comstr;
        getline(leggi, comstr);
	double comodo;
	int indice, tipo;
        leggi>>comodo;
        if (!leggi.eof())
        {
                if (comodo >= t00 && comodo <= t11)
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
                                coordinates[0][indice] = comodo;
                                leggi>>comodo;
                                coordinates[1][indice] = comodo;
                                leggi>>comodo;
                                coordinates[2][indice] = comodo;

                        }
                        getline(leggi, comstr);
                }
                else
                {
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        //number of atoms
                        leggi>>natoms;
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        for (int i=1; i<=natoms; i++)
                                getline(leggi, comstr);
                }
                return 1;
        }
        else
                return 0;
}

void cnt_bases(double *vettore, int nuc)
{
	int iA = 3*nuc - 1;
	int iB = natoms + 2 - iA;
	for (int k1=0; k1<3; k1++)
	    	vettore[k1] = 0.5*(coordinates[k1][iA] + coordinates[k1][iB]);
	return;
}

void cnt(double *vettore, int nuc)
{
	int iA = 3*nuc - 2;
	int iB = natoms - iA;
	for (int k1=0; k1<3; k1++)
	    	vettore[k1] = 0.5*(coordinates[k1][iA] + coordinates[k1][iB]);
	return;
}

void bpvec(double *vettore, int nuc)
{
	int iA = 3*nuc - 2;
	int iB = natoms - iA;
	for (int k1=0; k1<3; k1++)
	    	vettore[k1] = coordinates[k1][iB] - coordinates[k1][iA];
	return;
}

double compute_htwist(int nucinizio)
{
	double vec1[3], vec2[3];
	bpvec(vec1, nucinizio);	
	bpvec(vec2, nucinizio+1);
	double sp1 = dot_product(vec1, haxis, 3);
	double sp2 = dot_product(vec2, haxis, 3);
	for (int k1=0; k1<3; k1++)
	{
		vec1[k1] -= sp1*haxis[k1];
		vec2[k1] -= sp2*haxis[k1];
	}
	return angle_vec_sign(vec1, vec2, haxis);
}

double compute_hrise(int nucinizio)
{
    //h-rise is projection of dvec on helical axis
    double dvec[3], vec1[3], vec2[3];
    cnt(vec1, nucinizio);
    cnt(vec2, nucinizio+1);
    for (int k1=0; k1<3; k1++)
            dvec[k1] = vec2[k1] - vec1[k1];
    return dot_product(dvec, haxis, 3);
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
