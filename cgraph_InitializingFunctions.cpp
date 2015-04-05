#include"cgraph.h"
#include"MiscellaneousFunctions.h"
//#include"cnode_cedge.h"
#include<vector>
#include<iostream>
#include<fstream>
#include<cmath>
#include<algorithm>
using namespace std;

void cgraph::operator=(cgraph GivenGraph)
{
      n= GivenGraph.n;
      m= GivenGraph.m;
      Dimension= GivenGraph.Dimension;
      Vertices = new cnode[n];
      for(int i=0;i<n;i++) 
      { 
	    Vertices[i].Location = GivenGraph.Vertices[i].Location;
	    Vertices[i].Neighbors = GivenGraph.Vertices[i].Neighbors;	    
      }
      Hamiltonian= GivenGraph.Hamiltonian;      
}

void cCONNECTEDgraph::Initialize(int givenN, int givenDimension, double givenAverageDegree, char VertexDistributionType)
{
      n=givenN;
      Dimension=givenDimension;
      Vertices = new cnode[n];
      for(int i=0;i<n;i++) { Vertices[i].Location = new double[givenDimension];}
      if(VertexDistributionType == 's')
      {
	    SprinkleVerticesInSpace();
      }
      else if(VertexDistributionType == 'l')
      {
	    ArrangeVerticesInGrid();		  
      }
      
      CreateConnectedGraph(givenAverageDegree);      
}

void cCONNECTEDgraph::InitializeFromGivenData(string VertexLocationFileName, string EdgeListFileName, double givenAverageDegree)
{
      ifstream VertexLocationFile(VertexLocationFileName.c_str() ), EdgeListFile(EdgeListFileName.c_str() );
      
      vector<double> xcoordinate, ycoordinate; 
      n =0;
      double x,y;
      while( VertexLocationFile>>x>>y)
      {
	    xcoordinate.push_back(x);ycoordinate.push_back(y);
	    //cout<<"\nx = "<<x<<"\ty = "<<y;
	    n++;//cout<<"\t n="<<n;
      }
      
      Vertices = new cnode[n];
      Dimension=2; for(int i=0;i<n;i++) { Vertices[i].Location = new double[Dimension];}
      
      for(int i=0; i<n; i++)
      {
		  Vertices[i].Location[0] = xcoordinate[i]; Vertices[i].Location[1] = ycoordinate[i];
      }
      
      //cout<<"\n n="<<n;
      Hamiltonian =0.0;m=0;
      int OneEnd, OtherEnd;
      if(givenAverageDegree == 0.0 ) //Edges are given in datafile
      {
	    while(EdgeListFile>>OneEnd>>OtherEnd )
	    {
		  //cout<<"\nreached Boston";
		  if( OneEnd < OtherEnd) {CreateEdge(OneEnd,OtherEnd);}//this statement assumes that edges are listed with the smaller label first; this is to avoid multiple edges
	    }
      }
      else
      {
	    CreateConnectedGraph(givenAverageDegree); 
      }
      EdgeListFile.close();
      VertexLocationFile.close();
      
      //cout<<"\nspatial distance between vertex pairs\n";
      //for(int i = 0; i < n ; i++)
      //{
	//    for(int j=i+1; j<n; j++)
	  //  {
		//  double xx = SpatialDistanceBetweenVertices(i,j);
		  // if(xx == 0) {cout<<"|( "<<i<<","<<j<<" )| = "<<xx; Pause();}
//	    }
  //    }
}

void cgraph::InitializeForInfiniteBeta(int givenN, int givenDimension, double givenAverageDegree, char VertexDistributionType)
{
      n=givenN;
      Dimension=givenDimension;
      Vertices = new cnode[n];
      for(int i=0;i<n;i++) { Vertices[i].Location = new double[givenDimension];}
      if(VertexDistributionType == 's')
      {
	    SprinkleVerticesInSpace();
      }
      else if(VertexDistributionType == 'l')
      {
	    ArrangeVerticesInGrid();		  
	    for(int i =0; i<n; i++)
	    {
		  for(int j=0; j<Dimension; j++)
		  {
			Vertices[i].Location[j] += (2.0*RandomNumberBetweenZeroAndOne()-1.0)*.0001;			
		  }
	    }
      }
      int NumberOfEdgesToBeCreated = givenAverageDegree * n /2;
      m = 0;
      Hamiltonian = 0.0;
	    
      cedge* _TemporaryListOfEdges; _TemporaryListOfEdges= new cedge[NumberOfEdgesToBeCreated];
	    
      int NumberOfEdgesAddedToList=0;
      int ii=0,jj=0;
	    
      for(; (NumberOfEdgesAddedToList < NumberOfEdgesToBeCreated)&&(ii<n); ii++)
      {
	    for(jj=ii+1; (NumberOfEdgesAddedToList < NumberOfEdgesToBeCreated)&&(jj<n); jj++)
	    {
		  _TemporaryListOfEdges[NumberOfEdgesAddedToList] = cedge(ii,jj,SpatialDistanceBetweenVertices(ii,jj));
		  NumberOfEdgesAddedToList++;
	    }
      }
	    
      make_heap(_TemporaryListOfEdges, _TemporaryListOfEdges+NumberOfEdgesToBeCreated);
	    
      for(; jj < n; jj++)
      {
	    cedge CurrentMaxEdge = _TemporaryListOfEdges[0], EdgeBeingConsidered(ii,jj,SpatialDistanceBetweenVertices(ii,jj));
		  
	    if( EdgeBeingConsidered < CurrentMaxEdge )
	    {
		  pop_heap(_TemporaryListOfEdges, _TemporaryListOfEdges+NumberOfEdgesToBeCreated); //_TemporaryListOfEdges.pop_back();
		  _TemporaryListOfEdges[NumberOfEdgesToBeCreated-1] = EdgeBeingConsidered; push_heap(_TemporaryListOfEdges, _TemporaryListOfEdges+NumberOfEdgesToBeCreated);
	    }		  
      }
	    
      for(; ii < n; ii++)
      {
	    for(jj=ii+1; jj<n; jj++)
	    {
		  cedge CurrentMaxEdge = _TemporaryListOfEdges[0], EdgeBeingConsidered(ii,jj,SpatialDistanceBetweenVertices(ii,jj));
			
		  if( EdgeBeingConsidered < CurrentMaxEdge )
		  {
			pop_heap(_TemporaryListOfEdges, _TemporaryListOfEdges+NumberOfEdgesToBeCreated); //_TemporaryListOfEdges.pop_back();
			_TemporaryListOfEdges[NumberOfEdgesToBeCreated-1] = EdgeBeingConsidered; push_heap(_TemporaryListOfEdges, _TemporaryListOfEdges+NumberOfEdgesToBeCreated);
		  }		  
	    }
      }
	    
      for(int i=0; i< NumberOfEdgesToBeCreated; i++)
      {
	    CreateEdge(_TemporaryListOfEdges[i].SmallVertex, _TemporaryListOfEdges[i].LargeVertex );
      }	    
      delete[] _TemporaryListOfEdges;
}

void cgraph::Initialize(int givenN, int givenDimension, double givenAverageDegree, char VertexDistributionType)
{
      n=givenN;
      Dimension=givenDimension;
      Vertices = new cnode[n];
      for(int i=0;i<n;i++) { Vertices[i].Location = new double[givenDimension];}
      if(VertexDistributionType == 's')
      {
	    SprinkleVerticesInSpace();
      }
      else if(VertexDistributionType == 'l')
      {
	    ArrangeVerticesInGrid();		  
      }
	    
      CreateErdosRenyi(givenAverageDegree);           
}

void cgraph::InitializeForPercolation(int givenN, int givenDimension, double Beta, double Kappa, char VertexDistributionType)
{
      n=givenN;
      Dimension=givenDimension;
      Vertices = new cnode[n];
      for(int i=0;i<n;i++) { Vertices[i].Location = new double[givenDimension];}
      if(VertexDistributionType == 's')
      {
	    SprinkleVerticesInSpace();
      }
      else if(VertexDistributionType == 'l')
      {
	    ArrangeVerticesInGrid();
      }
      CreateByPercolationModel(Beta, Dimension, Kappa);      
}

void cgraph::CreateErdosRenyi(double AverageDegree) 
{
      int NumberOfEdgesToBeCreated = AverageDegree * n /2;
      m = 0;
      Hamiltonian = 0.0;
      for(;m < NumberOfEdgesToBeCreated;) 
      { 
	    
	    int ChosenVertex = RandomNumberBetweenZeroAndOne()*n,  ChosenOtherVertex = RandomNumberBetweenZeroAndOne()*n;
	    CreateEdge(ChosenVertex, ChosenOtherVertex);
      }
}

void cgraph::CreateConnectedGraph(double AverageDegree) 
{
      int NumberOfEdgesToBeCreated = AverageDegree * n /2;
      m = 0;
      Hamiltonian = 0.0;
      
      for(int i=0; i<n-1; i++)
      {
	    CreateEdge(i,i+1);
      }
      
      for(;m < NumberOfEdgesToBeCreated;) 
      {	    
	    int ChosenVertex = RandomNumberBetweenZeroAndOne()*n,  ChosenOtherVertex = RandomNumberBetweenZeroAndOne()*n;
	    CreateEdge(ChosenVertex, ChosenOtherVertex);
      }
}

void cgraph::SprinkleVerticesInSpace()
{
      for(int i=0; i<n; i++)
      {
	    for(int j=0;j<Dimension;j++)
	    {
		  Vertices[i].Location[j] = RandomNumberBetweenZeroAndOne() * pow(double(n),1.0/double(Dimension));
	    }
      }
}

void cgraph::ArrangeVerticesInGrid()
{
      /*code for 2D only*/
      int LengthOfSide = sqrt(n);
      n = LengthOfSide * LengthOfSide;
      
      int CurrentVertex = 0;
      
      for(int Row= 0 ; Row<LengthOfSide; Row++)
      {
	    for(int Column=0;Column<LengthOfSide;Column++)
	    {
		  Vertices[CurrentVertex].Location[0] = Row + 0.5;
		  Vertices[CurrentVertex].Location[1] = Column + 0.5 ;
		  CurrentVertex++;
	    }
      }
}


void cgraph::CreateByPercolationModel(double Beta, int Dimension, double Kappa)
{
      m=0;
      Hamiltonian =0.0;
      for(int i=0; i<n; i++)
      {
	    for(int j=i+1; j<n;j++)
	    {
		  if(Bernoulli( Kappa / (Kappa + exp( Beta * SpatialDistanceBetweenVertices(i,j) ) ) ))
		  {
			CreateEdge(i,j);			
		  }		  
	    }
      }
}


void cgraph::Terminate() 
{
      delete[] Vertices;
}
