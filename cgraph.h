//#include "HashDefineStuff.h"
//#include<iostream>
//#include<cmath>
//#include<ctime>
#include<fstream>
//#include<stdio.h>
//#include<stdlib.h>
#include<string>
//#include <sstream>
//#include <boost/random.hpp>
#include<vector>
#include"cnode_cedge.h"
#include<utility>

using namespace std;

class cgraph{
      protected:      
	    int n,m,Dimension;
	    
	    cnode* Vertices;
	    double Hamiltonian;
	    
	    double SpatialDistanceBetweenVertices(int, int);
	    
	    void CreateByConfigurationModel(char);
	    
	    void CreateByPercolationModel(double, int, double);
	    void CreateErdosRenyi(double);
	    void CreateConnectedGraph(double);
	    
	    //void CreateLattice();
	    //void CreateStar(int);
	    //void CreateTree(int);
	    
	    void SprinkleVerticesInSpace();
	    void ArrangeVerticesInGrid();
	    
	    void CreateEdge(int, int);
	    void BreakEdge(int, int);
	    void BreakEdge(cedge); 
	    void Rewire(int,int,int);	    
	    
	    bool AreVerticesConnected(int , int );
	    cedge RandomEdge();
	    int RandomVertexOutsideNeighborhoodOf(int);
	    vector<int> Neighbors(int);
	    
	    void operator=(cgraph); 
	    
      public:
	    double AverageDegree();
	    double FractionOfVerticesInLargestComponent();
	    double ClusteringCoefficient();
	    void ComponentSizeDistribution(int*);
	    bool EquilibriumReached(int, double*,bool,ofstream&);
	    double EnergyDensity();
	    virtual void Evolve(double);
	    //virtual void Evolve();
	    void Terminate();
	    //void WriteDegreeDistributionToFile(ofstream&);
	    int Degree(int);
	    vector<double> FractionOfVerticesOfDegree();
	    vector<double> FractionOfEdgeLengthsAround();
	    void Display();
	    void WriteToFile(ofstream&);
	    virtual void Initialize(int, int, double,char);
	    void InitializeForPercolation(int, int, double, double,char);
	    virtual void InitializeFromGivenData(string, string, double){}
	    void InitializeForInfiniteBeta(int, int, double, char);
	    virtual double AverageHopDistance() { return 0.0;}
	    virtual double AverageRouteDistance() { return 0.0;}
	    virtual double RouteFactor() { return 0.0;}
	    virtual bool IsGraphConnected() { return true;}
	    //virtual vector<pair<double,double> > Robustness(){ }
	    vector<pair<double,double> > Robustness();
};

class cCONNECTEDgraph: public cgraph {
      protected:
       	    bool IsConnectedGraphConnectedAfterRewire(int , int , int );
	    bool IsThereAPathBetweenVertices(int , int );
	    double AverageHopDistanceToAllOtherVertices(int);
	    double AverageRouteDistanceToAllOtherVertices(int);
	    double RouteFactorOfVertex(int );

      public:
	    void Evolve(double);
	    void Initialize(int, int, double, char);
	    void InitializeFromGivenData(string, string, double);
	    double AverageHopDistance();
	    double AverageRouteDistance();
	    double RouteFactor();
	    bool IsGraphConnected();
	    //vector<pair<double,double> > Robustness();
};