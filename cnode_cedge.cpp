#include<iostream>
#include<algorithm>
//#include<cmath>
#include"cnode_cedge.h"
#include "MiscellaneousFunctions.h"

int partition(cedge* input, int p, int r)
{
    //cedge pivot = input[r];
    double LengthOfPivot = input[r].Length;
    
    while ( p < r )
    {
        while ( input[p].Length < LengthOfPivot )
            p++;
        
        while ( input[r].Length > LengthOfPivot )
            r--;
        
        if ( input[p].Length == input[r].Length )
            p++;
        else if ( p < r ) {
            cedge tmp = input[p];
            input[p] = input[r];
            input[r] = tmp;
        }
    }
    
    return r;
}

cedge quick_select(cedge* input, int p, int r, int k)
{
    if ( p == r ) return input[p];
    int j = partition(input, p, r);
    int length = j - p + 1;
    if ( length == k ) return input[j];
    else if ( k < length ) return quick_select(input, p, j - 1, k);
    else  return quick_select(input, j + 1, r, k - length);
}


cedge::cedge(int a, int b, double givenLength) 
{
      SmallVertex = min(a,b);
      LargeVertex = max(a,b);
      Length = givenLength;
}

cedge::cedge(int a, int b) 
{
      SmallVertex = min(a,b);
      LargeVertex = max(a,b);
}

cedge::cedge()
{
}

bool cedge::operator==(cedge Other) 
{
      return ((SmallVertex == Other.SmallVertex)&&(LargeVertex==Other.LargeVertex));
}

void cedge::operator=(cedge Other)
{
      SmallVertex = Other.SmallVertex;
      LargeVertex = Other.LargeVertex;
      Length = Other.Length;
}

bool cedge::operator<(cedge Other)
{
      return (Length < Other.Length);
}

bool operator<(const cedge& lhs, const cedge& rhs )
{
    return lhs < rhs;
}


void cedge::Display()
{
      cout<<"["<<SmallVertex<<","<<LargeVertex<<"]";
}

cgraphstats::cgraphstats()
{
	EnergyDensity = 0.0;
	AverageDegree = 0.0;
	FractionOfVerticesInLargestComponent=0.0;
	ClusteringCoefficient=0.0;
	RouteFactor =0.0;
	AverageHopDistance =0.0;
	AverageRouteDistance = 0.0;
}

cgraphstats::cgraphstats(double givenEnergyDensity, double givenAverageDegree, double givenFractionOfVerticesInLargestComponent, double givenAverageHopDistance, double givenAverageRouteDistance,vector<double> givenFractionOfVerticesOfDegree, vector<double> givenFractionOfEdgeLengthsAround, vector<pair<double, double> > givenRobustness,double givenClusteringCoefficient, double givenRouteFactor)
{
	EnergyDensity = givenEnergyDensity;
	AverageDegree = givenAverageDegree;
	FractionOfVerticesInLargestComponent=givenFractionOfVerticesInLargestComponent;
	AverageHopDistance = givenAverageHopDistance;
	AverageRouteDistance =  givenAverageRouteDistance;
	ClusteringCoefficient= givenClusteringCoefficient;
	RouteFactor = givenRouteFactor;
	
	FractionOfVerticesOfDegree.assign( givenFractionOfVerticesOfDegree.begin(), givenFractionOfVerticesOfDegree.end());
	FractionOfEdgeLengthsAround.assign( givenFractionOfEdgeLengthsAround.begin(), givenFractionOfEdgeLengthsAround.end());
	Robustness.assign( givenRobustness.begin(), givenRobustness.end());
	
}


cgraphstats cgraphstats::operator+ (cgraphstats Given)
{
      return cgraphstats(EnergyDensity + Given.EnergyDensity, AverageDegree + Given.AverageDegree, FractionOfVerticesInLargestComponent + Given.FractionOfVerticesInLargestComponent, AverageHopDistance + Given.AverageHopDistance , AverageRouteDistance + Given.AverageRouteDistance, AddVectors(FractionOfVerticesOfDegree , Given.FractionOfVerticesOfDegree), AddVectors(FractionOfEdgeLengthsAround , Given.FractionOfEdgeLengthsAround), AddVectors( Robustness, Given.Robustness),ClusteringCoefficient + Given.ClusteringCoefficient, RouteFactor + Given.RouteFactor );
}


cgraphstats cgraphstats::operator/ (int factor)
{
      cgraphstats temp;
      temp.EnergyDensity = EnergyDensity / factor;
      temp.AverageDegree = AverageDegree / factor;
      temp.FractionOfVerticesInLargestComponent = FractionOfVerticesInLargestComponent / factor;
      temp.ClusteringCoefficient = ClusteringCoefficient / factor;
      temp.RouteFactor = RouteFactor / factor;
      
      temp.AverageHopDistance = AverageHopDistance / factor;
      temp.AverageRouteDistance = AverageRouteDistance / factor;
      
      temp.FractionOfVerticesOfDegree = FractionOfVerticesOfDegree;
      temp.FractionOfEdgeLengthsAround = FractionOfEdgeLengthsAround;
      
      for(int i=0; i< FractionOfVerticesOfDegree.size();i++ )
      {
	    temp.FractionOfVerticesOfDegree[i] = FractionOfVerticesOfDegree[i] / factor;	    
      }
      for(int i=0; i< FractionOfEdgeLengthsAround.size();i++ )
      {
	    temp.FractionOfEdgeLengthsAround[i] = FractionOfEdgeLengthsAround[i] / factor;
      }
    

      for(int i=0; i< Robustness.size();i++ )
      {
	    	    	
		temp.Robustness.push_back(make_pair(Robustness[i].first / factor ,Robustness[i].second / factor ) );	      
      }
      
      
      return temp;
}

void cgraphstats::Display()
{
      cout<<"\nEnergyDensity ="<<EnergyDensity<<"\nAverageDegree ="<<AverageDegree<<"\nFractionOfVerticesInLargestComponent = "<<FractionOfVerticesInLargestComponent;
      cout<<"\nDegree\tFraction Of Vertices";
      for(int i=0; i<FractionOfVerticesOfDegree.size();i++)
      {
	    cout<<"\n"<<i<<"\t"<<FractionOfVerticesOfDegree[i];
      }
}