#include<iostream>
#include<cmath>
//#include<ctime>
#include<fstream>
#include "cgraph.h"
#include<string>
#include"MiscellaneousFunctions.h"


cgraphstats SingleRun(int givenN,double givenDimension, double Kappa, double givenAverageDegree, double Beta, bool RecordTimeSeries,string Outputs, char VertexDistributionType, int TotalNumberOfParameterValues, string FileNameInputPart, char ModelType)
{    
      cgraph *_Graph;
      ofstream myFileT;

      cgraphstats ResultOfCurrentRun;

      if(ModelType == 's')
      {
	    _Graph = new cgraph;
	    if(Kappa ==0.0)
	    {
		  if(Beta > givenN) //beta = infinity
		  {
			// ***not sure why i have this condition***if(givenAverageDegree <0)
			if(givenAverageDegree >0)
			{
			      _Graph->InitializeForInfiniteBeta(givenN, givenDimension, givenAverageDegree, VertexDistributionType);
			}			      
		  } 
		  else 
		  {
			_Graph->Initialize(givenN, givenDimension, givenAverageDegree, VertexDistributionType);
			
			double EnergyDensityAtLastCheckPoint = _Graph->EnergyDensity();
			
			if(RecordTimeSeries)
			{
			      cout<<"\n Graph creation done......";
			      string TimeSeriesFileName = "IN"+FileNameInputPart+"OUTtimeseries";
			      myFileT.open(TimeSeriesFileName.c_str() );
			      
			      myFileT<<"n="<<givenN<<"D="<<givenDimension<<"mu="<<givenAverageDegree<<"beta="<<Beta;
			      myFileT<<"\nTimeStep\t EnergyDensity \tLargestComponent\t EquilibriumTestingParameter";     
			}
			
			int TimeStep =0;
			if(Beta != 0.0) 
			{ 
			      for(TimeStep=1;  !_Graph->EquilibriumReached(TimeStep,&EnergyDensityAtLastCheckPoint, RecordTimeSeries, myFileT);TimeStep++)
			      {
				    _Graph->Evolve(Beta);	    
			      }
			} 
		  }
	    }
	    else
	    {//percolation model
		  _Graph->InitializeForPercolation(givenN, givenDimension, Beta, Kappa, VertexDistributionType );
		  cout<<"\n_Graph->AverageDegree() ="<<_Graph->AverageDegree();	    
	    }      
      }
      else if(ModelType == 'o')
      {
	    _Graph = new cCONNECTEDgraph;
	    _Graph->Initialize(givenN, givenDimension, givenAverageDegree, VertexDistributionType);
      
	    double EnergyDensityAtLastCheckPoint = _Graph->EnergyDensity();
      
	    if(RecordTimeSeries)
	    {
		  cout<<"\n Graph creation done......";
		  string TimeSeriesFileName = "IN"+FileNameInputPart+"OUTtimeseries";
		  myFileT.open(TimeSeriesFileName.c_str() );
	    
		  myFileT<<"n="<<givenN<<"D="<<givenDimension<<"mu="<<givenAverageDegree<<"beta="<<Beta;
		  myFileT<<"\nTimeStep\t EnergyDensity \tLargestComponent\t EquilibriumTestingParameter";     
	    }
      
	    int TimeStep =0;
	    //if(Beta != 0.0) 
			
	    for(TimeStep=1;  !_Graph->EquilibriumReached(TimeStep,&EnergyDensityAtLastCheckPoint, RecordTimeSeries, myFileT);TimeStep++)
	    {
		  _Graph->Evolve(Beta);	    
	    }
      }

ResultOfCurrentRun.EnergyDensity = _Graph->EnergyDensity();
ResultOfCurrentRun.AverageDegree =_Graph->AverageDegree();

if(ModelType == 'o') {ResultOfCurrentRun.FractionOfVerticesInLargestComponent =1;}
if( Outputs.find('c') != std::string::npos ) {ResultOfCurrentRun.FractionOfVerticesInLargestComponent = _Graph->FractionOfVerticesInLargestComponent();}
if( Outputs.find('h') != std::string::npos ) { ResultOfCurrentRun.AverageHopDistance = _Graph->AverageHopDistance(); cout<<"\nHopDistance ="<< ResultOfCurrentRun.AverageHopDistance;}
if( Outputs.find('r') != std::string::npos ) { ResultOfCurrentRun.AverageRouteDistance = _Graph->AverageRouteDistance();}
if( Outputs.find('R') != std::string::npos ) { ResultOfCurrentRun.RouteFactor = _Graph->RouteFactor();}
if( Outputs.find('C') != std::string::npos ) { ResultOfCurrentRun.ClusteringCoefficient =_Graph->ClusteringCoefficient();}
if( Outputs.find('D') != std::string::npos ) {ResultOfCurrentRun.FractionOfVerticesOfDegree = _Graph->FractionOfVerticesOfDegree();}
if( Outputs.find('E') != std::string::npos ) { ResultOfCurrentRun.FractionOfEdgeLengthsAround = _Graph->FractionOfEdgeLengthsAround();}
if( Outputs.find('V') != std::string::npos ) {ResultOfCurrentRun.Robustness = _Graph->Robustness();}
if( Outputs.find('P') != std::string::npos )       
{
      ofstream myFilePicture;
      string PictureFileName = "IN"+FileNameInputPart+"OUTpicture";
      myFilePicture.open(PictureFileName.c_str());
      _Graph->WriteToFile(myFilePicture);
      myFilePicture.close();
}
_Graph->Terminate();
if(RecordTimeSeries) 
{
      myFileT.close();
}
return ResultOfCurrentRun;      
}