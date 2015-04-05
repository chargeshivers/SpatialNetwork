#include<iostream>
#include<cmath>
#include<ctime>
//#include<fstream>
//#include<stdio.h>
//#include<stdlib.h>
#include<string>
#include <sstream>
//#include <boost/random.hpp>
#include <random>
#include<vector>
//#include"MiscellaneousFunctions.h"
using namespace std;
//#include"cnode_cedge.h"
vector<pair<double,double> > AddVectors(vector<pair<double,double> > Vector1, vector<pair<double,double> > Vector2)
{
      vector<pair<double,double> > Result;
      
      Result = Vector2;
      for(int i=0; i< Vector1.size();i++)
      {
	    Result[i].first += Vector1[i].first;
	    Result[i].second += Vector1[i].second;
      }
      return Result;
}

vector<double> AddVectors(vector<double> Vector1, vector<double> Vector2)
{
      vector<double> Result, SmallVector;
      if( Vector1.size() > Vector2.size() )
      {
	    Result = Vector1;
	    SmallVector = Vector2;
      }
      else
      {
	    Result = Vector2;
	    SmallVector = Vector1;
      }
      
      for(int i=0; i< SmallVector.size();i++ )
      {
	    Result[i] += SmallVector[i];
      }
      
      return Result;
}

vector<double> GenerateDistribution(int* _List, int SizeOfList)
{
      int Largest = 0;
      for(int i=0; i<SizeOfList; i++)
      {
	    if(_List[i] > Largest) {Largest = _List[i];}
      }
      
      vector<int> NumberOfEntriesOfValue (Largest+1,0);
      vector<double> FractionOfEntriesOfValue (Largest+1,0.0);
      
      for(int i=0; i<SizeOfList; i++)
      {
	    
	    NumberOfEntriesOfValue[_List[i]]++;
      }
      
      
      for(int j=0; j<Largest+1; j++)
      {
	    FractionOfEntriesOfValue[j] = double(NumberOfEntriesOfValue[j]) / double(SizeOfList) ;
      }
      return FractionOfEntriesOfValue;
}

vector<double> GenerateDistribution(double* _List, int SizeOfList)
{
      int _NewList[SizeOfList];
      
      for(int i=0; i<SizeOfList; i++ )
      {
	    _NewList[i] = floor(_List[i] + 0.5); //nearest integer
      }
      return GenerateDistribution(_NewList, SizeOfList);
}

bool IsExtrema(int xCurrent, int xLast, int xLastToLast  )
{
      if(( xLastToLast != xLast ) && ( (xCurrent > xLast )  != (xLast > xLastToLast )) ) return true; else return false;     
}

int Binom2(int x) 
{
      return x * (x-1) /2;
}

double RandomNumberBetweenZeroAndOne()
{ 
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
    }

bool Bernoulli(double Probability)
{
      if( RandomNumberBetweenZeroAndOne() <= Probability) return true; else  return false;
}

string DoubleToString(double x)
{
      string myString; 
      stringstream ss; 
      ss<< x;
      ss>>myString;
      return myString;
}

int RandomElement(vector<int> GivenVector)
{
      return GivenVector[ RandomNumberBetweenZeroAndOne() * GivenVector.size()];
}

int SampleFromDegreeDistribution(char NetworkType)
{
      if(NetworkType == 'P') {return  3 * pow( RandomNumberBetweenZeroAndOne(),-1.0/(2.5-1) );} //-1.0/(alpha-1) }
      if(NetworkType == 'R') {return 4;}
    return 0;
    
}

void Pause()
{
      char x;
    cout<<"\npausing....enter some key to continue....";
    cin>>x;
}