#include<fstream>
#include<string>
#include<vector>
//#include "cnode_cedge.h"
using namespace std;

vector<double> GenerateDistribution(int*, int);

vector<double> GenerateDistribution(double*, int);

vector<double> AddVectors(vector<double>, vector<double>);

vector<pair<double,double> > AddVectors(vector<pair<double,double> >, vector<pair<double,double> >);

int Binom2(int);

double RandomNumberBetweenZeroAndOne();

bool Bernoulli(double Probability);

   int RandomElement(vector<int>);

string DoubleToString(double);

int SampleFromDegreeDistribution(char);

void Pause();

//int partition(cedge* , int , int );

//cedge quick_select(cedge* , int , int , int );

