#include"cgraph.h"
#include"MiscellaneousFunctions.h"
#include<cmath>
#include<iostream>
#include<fstream>
#include<queue>
//#include <climits>


double cgraph::ClusteringCoefficient()
{
      int NumberOfWedges =0, NumberOfClosedWedges =0;
      //cout<<"\nNumberOfClosedWedges ="<<NumberOfClosedWedges<<"\tNumberOfWedges ="<<NumberOfWedges;
      for(int i=0; i<n; i++)
      {
	    int DegreeOfi = Degree(i);
	    NumberOfWedges += Binom2(DegreeOfi );	   
	    
	    for(int j=0; j< DegreeOfi; j++ )
	    {
		  for(int k=j+1; k< DegreeOfi; k++)
		  {
			if( AreVerticesConnected( Vertices[i].Neighbors[j],Vertices[i].Neighbors[k]) ) 
			{
			      NumberOfClosedWedges++;
			}
		  }	    
	    }
      }
      return double(NumberOfClosedWedges)/double(NumberOfWedges);
}

//Hop distance and related functions starts below..................................................

typedef pair<int,int> Vertex_HopDist_pair;

 vector<pair<double,double> > cgraph::Robustness()
 {
       vector<pair<double,double> > ReturnQuantity;
      double pmax = 1;
      int NumberOfpValues = 100;

             //cout<<"\nhere is the original graph.......";this->Display();

      ReturnQuantity.push_back(make_pair(this->AverageDegree(), 1) );
      for(int i= 0; i<NumberOfpValues; i++)
      {
	    double p = 0.0 + (pmax-0.0) * double(i+1)/double(NumberOfpValues);
	    
	    int NumberOfEdgesToBeDeleted = p*m;
	    
	    cgraph DamagedGraph;
	    DamagedGraph = *this;
	    
	    for(int j = 0; j< NumberOfEdgesToBeDeleted; j++)
	    {
		  DamagedGraph.BreakEdge(DamagedGraph.RandomEdge());
	    }
	    //now we have DamagedGraph with exactly (1-p) m edges	    
	    
	    ReturnQuantity.push_back(make_pair(DamagedGraph.AverageDegree(), DamagedGraph.FractionOfVerticesInLargestComponent()) );
	    
	    DamagedGraph.Terminate();
      }
      
      return ReturnQuantity; 
 }


class HopComparator
{
public:
      int operator() ( const Vertex_HopDist_pair& p1, const Vertex_HopDist_pair& p2)
      {
	    return p1.second > p2.second;
      }
};

double cCONNECTEDgraph::AverageHopDistance()
{
      double Sum =0.0;
      for(int i=0; i<n; i++)
      {
	    //cout<<"\nAverageHopDistances for vertex ["<< i<<" ] ="<<AverageHopDistanceToAllOtherVertices(i);
	    Sum += AverageHopDistanceToAllOtherVertices(i);
      }
      
      return Sum/double(n);      
}

double cCONNECTEDgraph::AverageHopDistanceToAllOtherVertices(int SourceVertex)
{
      int dist[n];     
      
      for(int i=0; i<n; i++)
      {
	    dist[i] = 999999;
	    //cout<<"\ndist["<<i<<"] = "<<dist[i];
      }
      dist[SourceVertex] =0;      
      
      priority_queue<Vertex_HopDist_pair,vector<Vertex_HopDist_pair>,HopComparator> PQ;
      
      PQ.push(Vertex_HopDist_pair(SourceVertex,dist[SourceVertex]));
      
      
      while(!PQ.empty())
      {
	    Vertex_HopDist_pair Current  = PQ.top(); PQ.pop();
	    int CurrentVertex = Current.first;
	    //int DistanceOfCurrentVertex = Current.second;
	    
	    //if(DistanceOfCurrentVertex <= dist[CurrentVertex]) //not sure if this is needed
	    //{
		  for(int i =0; i< Vertices[CurrentVertex].Neighbors.size(); i++  )
		  {
			int CurrentNeighbor = Vertices[CurrentVertex].Neighbors[i], DistanceFromCurrentVertexToNeighbor = 1;
			int DistanceOfNeighborThroughCurrentVertex = dist[CurrentVertex] + DistanceFromCurrentVertexToNeighbor;
			if( DistanceOfNeighborThroughCurrentVertex  < dist[CurrentNeighbor] )
			{
			      dist[CurrentNeighbor] = DistanceOfNeighborThroughCurrentVertex;
			      
			      PQ.push(Vertex_HopDist_pair(CurrentNeighbor,dist[CurrentNeighbor]));
			}
		  }
	    //}
      }
      
      int Sum =0;
      for(int i=0; i<n; i++)
      {
	    Sum += dist[i];	    
      }
      return double(Sum)/double(n-1);
}
//Hop distance and related functions ends above..................................................

//Route distance and related functions starts below..................................................

typedef pair<int,double> Vertex_RouteDist_pair;

class RouteComparator
{
public:
      int operator() ( const Vertex_RouteDist_pair& p1, const Vertex_RouteDist_pair& p2)
      {
	    return p1.second > p2.second;
      }
};

double cCONNECTEDgraph::RouteFactorOfVertex(int SourceVertex)
{
      double dist[n];     
      
      for(int i=0; i<n; i++)
      {
	    dist[i] = 999999;	    
      }
      dist[SourceVertex] =0.0;      
      
      priority_queue<Vertex_RouteDist_pair,vector<Vertex_RouteDist_pair>,RouteComparator> PQ;
      
      PQ.push(Vertex_RouteDist_pair(SourceVertex,dist[SourceVertex]));
      
      
      while(!PQ.empty())
      {
	    Vertex_RouteDist_pair Current  = PQ.top(); PQ.pop();
	    int CurrentVertex = Current.first;
	    double DistanceOfCurrentVertex = Current.second;
	    
	    if(DistanceOfCurrentVertex <= dist[CurrentVertex]) //not sure if this is needed
	    {
		  for(int i =0; i< Vertices[CurrentVertex].Neighbors.size(); i++  )
		  {
			int CurrentNeighbor = Vertices[CurrentVertex].Neighbors[i];
			double DistanceFromCurrentVertexToNeighbor = SpatialDistanceBetweenVertices(CurrentVertex, CurrentNeighbor);
			double DistanceOfNeighborThroughCurrentVertex = dist[CurrentVertex] + DistanceFromCurrentVertexToNeighbor;
			if( DistanceOfNeighborThroughCurrentVertex  < dist[CurrentNeighbor] )
			{
			      dist[CurrentNeighbor] = DistanceOfNeighborThroughCurrentVertex;
			      
			      PQ.push(Vertex_RouteDist_pair(CurrentNeighbor,dist[CurrentNeighbor]));
			}
		  }
	    }
      }
      
      double Sum =0.0;
      for(int i=0; i<n; i++)
      {
	    if(i != SourceVertex)
	    {
		  //cout<<"\ndist["<<i<<"] = "<<dist[i]<<"\t SpatialDistanceBetweenVertices("<< SourceVertex<<","<<i<<")"<<SpatialDistanceBetweenVertices(SourceVertex,i);
		  Sum += dist[i] / SpatialDistanceBetweenVertices(SourceVertex,i) -1;	    
		  //cout<<"\ndist[i] / SpatialDistanceBetweenVertices(SourceVertex,i) -1 = "<<dist[i] / SpatialDistanceBetweenVertices(SourceVertex,i) -1;
		  //cout<<"\nSum = "<<Sum;
		  if( Sum/double(n-1) > 9999) { cout<<"\nSourceVertex = "<<SourceVertex<<"\ti="<<i<<"\tdist[i] ="<<dist[i]<<"\tSpatialDistanceBetweenVertices(SourceVertex,i) = "<<SpatialDistanceBetweenVertices(SourceVertex,i)<<"\tdist[i] / SpatialDistanceBetweenVertices(SourceVertex,i) -1 = "<<dist[i] / SpatialDistanceBetweenVertices(SourceVertex,i) -1;Pause();}
	    }
      }
      if(Sum/double(n-1) > 9999) {cout<<"\nRouteFactorOfVertex "<<SourceVertex<<" = "<<Sum/double(n-1);Pause(); }
      
      return Sum/double(n-1);
}

double cCONNECTEDgraph::RouteFactor()
{	
      double Sum =0.0;
      
      for(int SourceVertex=0; SourceVertex <n; SourceVertex++)
      {
	    Sum += RouteFactorOfVertex(SourceVertex);
	    //cout<<"\nSum = "<<Sum; if(SourceVertex %100 ==0) Pause();
	    
      }
      cout<<"\nSum/double(n) = "<<Sum/double(n);
      return Sum/double(n);
}


double cCONNECTEDgraph::AverageRouteDistance()
{
      double Sum =0.0;
      for(int i=0; i<n; i++)
      {
	    Sum += AverageRouteDistanceToAllOtherVertices(i);
      }
      
      return Sum/double(n);      
}

double cCONNECTEDgraph::AverageRouteDistanceToAllOtherVertices(int SourceVertex)
{
      double dist[n];     
      
      for(int i=0; i<n; i++)
      {
	    dist[i] = 999999;	    
      }
      dist[SourceVertex] =0.0;      
      
      priority_queue<Vertex_RouteDist_pair,vector<Vertex_RouteDist_pair>,RouteComparator> PQ;
      
      PQ.push(Vertex_RouteDist_pair(SourceVertex,dist[SourceVertex]));
      
      
      while(!PQ.empty())
      {
	    Vertex_RouteDist_pair Current  = PQ.top(); PQ.pop();
	    int CurrentVertex = Current.first;
	   // double DistanceOfCurrentVertex = Current.second;
	    
	    //if(DistanceOfCurrentVertex <= dist[CurrentVertex]) //not sure if this is needed
	    //{
		  for(int i =0; i< Vertices[CurrentVertex].Neighbors.size(); i++  )
		  {
			int CurrentNeighbor = Vertices[CurrentVertex].Neighbors[i];
			double DistanceFromCurrentVertexToNeighbor = SpatialDistanceBetweenVertices(CurrentVertex, CurrentNeighbor);
			double DistanceOfNeighborThroughCurrentVertex = dist[CurrentVertex] + DistanceFromCurrentVertexToNeighbor;
			if( DistanceOfNeighborThroughCurrentVertex  < dist[CurrentNeighbor] )
			{
			      dist[CurrentNeighbor] = DistanceOfNeighborThroughCurrentVertex;
			      
			      PQ.push(Vertex_RouteDist_pair(CurrentNeighbor,dist[CurrentNeighbor]));
			}
		  }
	    //}
      }
      
      double Sum =0.0;
       //if(SourceVertex ==0) cout<<"\nThe following vertices are outside the giant component\n";
      for(int i=0; i<n; i++)
      {
	    Sum += dist[i];	
      }
      
      return Sum/double(n-1);
}
//Route distance and related functions ends above..................................................


bool cgraph::EquilibriumReached(int TimeStep, double* _EnergyDensityAtLastCheckPoint, bool RecordTimeSeries, ofstream& myFileT)
{
      //return true;
      if(TimeStep%(1000*m) ==0 )
      {
	    double delta = 0.005   /*0.005 if beta=0, 0.0001 if beta=10/*/;
	    
	    double EquilibriumTestingParameter = abs( (EnergyDensity() - *_EnergyDensityAtLastCheckPoint)/ *_EnergyDensityAtLastCheckPoint );
	    
	    *_EnergyDensityAtLastCheckPoint = EnergyDensity();
	    
	    if(RecordTimeSeries)
	    {
		  //double x = FractionOfVerticesInLargestComponent();
		  cout<<"\nTimeStep ="<<TimeStep<<"\tEnergyDensity = "<<EnergyDensity()/*<<"\tFraction of vertices in largest component = "<<x*/<<"\tEquilibriumTestingParameter = "<< EquilibriumTestingParameter; myFileT<<"\n"<<TimeStep<<"\t"<<EnergyDensity()/*<<"\t"<<x*/<<"\t"<< EquilibriumTestingParameter;
	    }
	    if( EquilibriumTestingParameter < delta  )
	    {
		  return true;
	    }
	    else
	    {		  
		  return false;
	    }
      }
      else
      {
	    return false;
      }
}

double cgraph::AverageDegree()
{
      return 2.0 * double(m) / double(n);
}

double cgraph::SpatialDistanceBetweenVertices(int ii, int jj)
{
      double sum = 0.0;
      for(int dimensionNumber=0; dimensionNumber < Dimension; dimensionNumber++)
      {
	    double x = Vertices[ii].Location[dimensionNumber] - Vertices[jj].Location[dimensionNumber];
	    sum += x * x;
      }
      return sqrt(sum);
}


cedge cgraph::RandomEdge()
{
      int PositionOfEdgeInDoubleListOfEdges = RandomNumberBetweenZeroAndOne()*2*m;
      int NumberOfEdgesParsed =0;
      for(int i=0;i<n;i++)
      {
	    for(int j=0; j<Vertices[i].Neighbors.size();j++)
	    {
		  if(PositionOfEdgeInDoubleListOfEdges == NumberOfEdgesParsed )
		  {
			return cedge(i,Vertices[i].Neighbors[j]);		 
		  }
		  NumberOfEdgesParsed++;		
	    }
      }
    return cedge();
}

int cgraph::RandomVertexOutsideNeighborhoodOf(int givenVertex)
{
      int ChosenVertex;
      for(;;)
      {
	    ChosenVertex = RandomNumberBetweenZeroAndOne()*n;
	    if( (ChosenVertex != givenVertex) && (!AreVerticesConnected(ChosenVertex, givenVertex) ) ) return ChosenVertex;
      }      
}

bool cgraph::AreVerticesConnected(int ii, int jj) 
{
      for(int i=0;i<Vertices[ii].Neighbors.size();i++) { 
	    if(Vertices[ii].Neighbors[i]==jj) return true;
      }
      return false;
}

int cgraph::Degree(int givenVertex)
{
      return static_cast<int>(Vertices[givenVertex].Neighbors.size());
}

double cgraph::EnergyDensity()
{
      return Hamiltonian/double(m);
}


double cgraph::FractionOfVerticesInLargestComponent()
{
      int SizeOfLargestComponent =1;

      bool visited[n];
      for(int i=0;i<n;i++)
      {
	    visited[i] = false;
      }
      
      for(int i=0;i<n;i++)
      {
	    if(!visited[i])
	    {
		  int SizeOfComponentContainingCurrentVertex = 1;
		  queue<int> q;
		  q.push(i);
		  visited[i] = true;
		  for(;!q.empty();)
		  {
			int w = q.front(); q.pop();
			for(int j=0; j<Vertices[w].Neighbors.size();j++)
			{      
			      if(!visited[Vertices[w].Neighbors[j]])
			      {
				    SizeOfComponentContainingCurrentVertex++;
				    visited[Vertices[w].Neighbors[j]] = true;
				    q.push(Vertices[w].Neighbors[j]);
			      }
			}
		  }
		  if(SizeOfComponentContainingCurrentVertex > SizeOfLargestComponent)
		  {
			SizeOfLargestComponent = SizeOfComponentContainingCurrentVertex;
		  }
	    }
      }
      return double(SizeOfLargestComponent)/double(n);
}


vector<double> cgraph::FractionOfVerticesOfDegree()
{
      int _ListOfVertexDegrees[n];
      
      for(int i=0; i<n; i++)
      {
	    _ListOfVertexDegrees[i] = Degree(i);
      }
      
      return GenerateDistribution(_ListOfVertexDegrees, n);
}

vector<double> cgraph::FractionOfEdgeLengthsAround()
{
      double _ListOfEdgeLengths[m];
      int NumberOfEdgesParsed =0;
      
      for(int i=0; i<n; i++)
      {
	    for(int j=0; j< Vertices[i].Neighbors.size(); j++ )
	    {
		  if( i < Vertices[i].Neighbors[j]  )
		  {
			_ListOfEdgeLengths[NumberOfEdgesParsed] = SpatialDistanceBetweenVertices(i,Vertices[i].Neighbors[j]);
			NumberOfEdgesParsed++;
		  }
	    }
      }
      return GenerateDistribution(_ListOfEdgeLengths, m);
}


void cgraph::ComponentSizeDistribution(int* NumberOfComponentsOfSize)
{
      int SizeOfLargestComponent =1;
      for(int i=0;i<n;i++)
      {
	    NumberOfComponentsOfSize[i] = 0;
      }

      bool visited[n];
      for(int i=0;i<n;i++)
      {
	    visited[i] = false;
      }
      
      for(int i=0;i<n;i++)
      {
	    if(!visited[i])
	    {
		  int SizeOfComponentContainingCurrentVertex = 1;
		  queue<int> q;
		  q.push(i);
		  visited[i] = true;
		  for(;!q.empty();)
		  {
			int w = q.front(); q.pop();
			for(int j=0; j<Vertices[w].Neighbors.size();j++)
			{      
			      if(!visited[Vertices[w].Neighbors[j]])
			      {
				    SizeOfComponentContainingCurrentVertex++;
				    visited[Vertices[w].Neighbors[j]] = true;
				    q.push(Vertices[w].Neighbors[j]);
			      }
			}
		  }
		  NumberOfComponentsOfSize[SizeOfComponentContainingCurrentVertex]++;
		  if(SizeOfComponentContainingCurrentVertex > SizeOfLargestComponent)
		  {
			SizeOfLargestComponent = SizeOfComponentContainingCurrentVertex;
		  }
	    }
      }
}

bool cCONNECTEDgraph::IsConnectedGraphConnectedAfterRewire(int Pivot, int OldNeighbor, int NewNeighbor)
{
      Rewire(Pivot, OldNeighbor, NewNeighbor);
      
      if(IsThereAPathBetweenVertices(Pivot,OldNeighbor)) { Rewire(Pivot, NewNeighbor, OldNeighbor); return true; } else {Rewire(Pivot, NewNeighbor, OldNeighbor); return false;}      
}



bool cCONNECTEDgraph::IsThereAPathBetweenVertices(int Source, int Destination)
{
      bool visited[n];
      for(int i=0;i<n;i++)
      {
	    visited[i] = false;
      }
		  queue<int> q;
		  q.push(Source);
		  visited[Source] = true;
		  for(;!q.empty();)
		  {
			int w = q.front(); q.pop();
			for(int j=0; j<Vertices[w].Neighbors.size();j++)
			{      
			      if(Vertices[w].Neighbors[j] == Destination ) return true;
			      
			      if(!visited[Vertices[w].Neighbors[j]])
			      {
				    visited[Vertices[w].Neighbors[j]] = true;
				    q.push(Vertices[w].Neighbors[j]);
			      }
			}
		  }
      return false;
}

bool cCONNECTEDgraph::IsGraphConnected()
{
      for(int i=0; i<n; i++)
      {
	    for(int j=i+1; j<n;j++)
	    {
		  if( !IsThereAPathBetweenVertices(i,j) ) {cout<<"\nThere is no path from vertex "<<i<<" to vertex "<<j<<"!!";return false;}
	    }
      }
      return true;
}
