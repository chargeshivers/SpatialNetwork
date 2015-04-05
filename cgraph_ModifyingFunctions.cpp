#include"cgraph.h"
#include"MiscellaneousFunctions.h"
#include<algorithm>
#include<iostream>
#include<cmath>

using namespace std;

void cgraph::Evolve(double beta)
{
      cedge ChosenEdge = RandomEdge();
      int ChosenPivot, ChosenOldNeighbor, ChosenNewNeighbor;
      if(Bernoulli(0.5))
      {
	    ChosenPivot = ChosenEdge.SmallVertex, ChosenOldNeighbor = ChosenEdge.LargeVertex; 
      }	    
      else
      {
	    ChosenPivot = ChosenEdge.LargeVertex, ChosenOldNeighbor = ChosenEdge.SmallVertex;
      }
      ChosenNewNeighbor = RandomVertexOutsideNeighborhoodOf(ChosenPivot);
      
      double OldEdgeLength = SpatialDistanceBetweenVertices(ChosenPivot, ChosenOldNeighbor), NewEdgeLength = SpatialDistanceBetweenVertices(ChosenPivot, ChosenNewNeighbor);
      if( NewEdgeLength < OldEdgeLength )
      {
	    Rewire(ChosenPivot,ChosenOldNeighbor,ChosenNewNeighbor);
      }
      else
      {
	    if( Bernoulli( exp( -beta*( NewEdgeLength - OldEdgeLength ) ) ) )
	    {
		  Rewire(ChosenPivot,ChosenOldNeighbor,ChosenNewNeighbor);
	    }
      }
}

void cCONNECTEDgraph::Evolve(double beta)
{
      cedge ChosenEdge = RandomEdge();
      int ChosenPivot, ChosenOldNeighbor, ChosenNewNeighbor;
      if(Bernoulli(0.5))
      {
	    ChosenPivot = ChosenEdge.SmallVertex, ChosenOldNeighbor = ChosenEdge.LargeVertex; 
      }	    
      else
      {
	    ChosenPivot = ChosenEdge.LargeVertex, ChosenOldNeighbor = ChosenEdge.SmallVertex;
      }
      ChosenNewNeighbor = RandomVertexOutsideNeighborhoodOf(ChosenPivot);
      
      double OldEdgeLength = SpatialDistanceBetweenVertices(ChosenPivot, ChosenOldNeighbor), NewEdgeLength = SpatialDistanceBetweenVertices(ChosenPivot, ChosenNewNeighbor);
      
      if( NewEdgeLength < OldEdgeLength )
      {
	    if(IsConnectedGraphConnectedAfterRewire(ChosenPivot,ChosenOldNeighbor,ChosenNewNeighbor)) 
	    {
		  Rewire(ChosenPivot,ChosenOldNeighbor,ChosenNewNeighbor);
	    }
      }
      else
      {
	    if( Bernoulli( exp( -beta*( NewEdgeLength - OldEdgeLength ) ) ) )
	    {
		  if(IsConnectedGraphConnectedAfterRewire(ChosenPivot,ChosenOldNeighbor,ChosenNewNeighbor)) 
		  {
			Rewire(ChosenPivot,ChosenOldNeighbor,ChosenNewNeighbor);
		  }
	    }
      }
}

void cgraph::Rewire(int Pivot, int OldNeighbor, int NewNeighbor)
{
      BreakEdge(Pivot, OldNeighbor);
      CreateEdge(Pivot,NewNeighbor);
}

void cgraph::BreakEdge(cedge GivenEdge) 
{
      BreakEdge(GivenEdge.SmallVertex, GivenEdge.LargeVertex);
}

void cgraph::BreakEdge(int ii, int jj) 
{
      if(ii == jj) 
      {	
	    cout<<"ERROR!, both vertices are same, so can't break edge";
      }
      else 
      {
	    if(!AreVerticesConnected(ii,jj))
	    {
		  //cout<<"ERROR!, "<<ii<<" and "<<jj<<" are not connected, so can't create break them";
	    }
	    
	    else 
	    {
		  for(int k=0; k<Vertices[ii].Neighbors.size(); k++ )
		  {
			if(Vertices[ii].Neighbors[k]==jj) 
			{
			      Vertices[ii].Neighbors.erase(Vertices[ii].Neighbors.begin()+k); 
			      break;
			}
		  }
		  
		  for(int k=0; k<Vertices[jj].Neighbors.size(); k++ )
		  {
			if(Vertices[jj].Neighbors[k]==ii) 
			{
			      Vertices[jj].Neighbors.erase(Vertices[jj].Neighbors.begin()+k); 
			      break;
			}
		  }		  
		  m--;
		  Hamiltonian -=  SpatialDistanceBetweenVertices(ii, jj);
	    }
      }
}

void cgraph::CreateEdge(int ii, int jj) 
{
      if(ii == jj) 
      {	
	    //cout<<"ERROR!, both vertices are same, so can't create edge";
      }
      else 
      {
	    if(false/*AreVerticesConnected(ii,jj)*/)
	    {
		  //cout<<"ERROR!, edge already exists between "<<ii<<" and "<<jj<<" , so can't create edge";
	    }
	    
	    else 
	    {
		  Vertices[ii].Neighbors.push_back(jj);
		  Vertices[jj].Neighbors.push_back(ii);
		  
		  m++;
		  Hamiltonian += SpatialDistanceBetweenVertices(ii, jj) ;
	    }
      }
}