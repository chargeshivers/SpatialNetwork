#include"cgraph.h"
#include<iostream>
#include"MiscellaneousFunctions.h"

void cgraph::WriteToFile(ofstream& Picture)
      {
	    //node x coordinate, y coordinate, label, neighbors...
	    
	    for(int i=0;i<n;i++) 
	    {
		  for(int j=0; j<Dimension; j++)
		  {
			Picture<<Vertices[i].Location[j]<<"\t";
		  }
		  Picture<<i+1<<"\t";
		  for(int j=0; j<Vertices[i].Neighbors.size();j++ )
		  {
			if( Vertices[i].Neighbors[j] > i ) 
			{
				Picture<<Vertices[i].Neighbors[j]+1<<"\t";
		 	} 
		}
		  Picture<<"\n";
	    }
      }

void cgraph::Display()
      { 
	    for(int i=0;i<n;i++) 
	    {
		  cout<<"\nVertices["<<i<<"].Location=(";
		  for(int j=0;j<Dimension;j++) 
		  {
			cout<<Vertices[i].Location[j]<<",";
		  }
		  cout<<")";
		  cout<<"\tVertices["<<i<<"].Neighbors=";
		  for(int j=0;j<Vertices[i].Neighbors.size();j++) 
		  {
			cout<<Vertices[i].Neighbors[j]<<",";
		  }
	    }
	    Pause();
      }
