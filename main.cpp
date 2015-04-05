#include<iostream>
#include<cmath>
#include<ctime>
#include<fstream>
#include <sstream>
#include<cstdlib>
#include<string>
#include<vector>
#include"MiscellaneousFunctions.h"
#include"cnode_cedge.h"

int main(int argc, char* argv[])
{
      int CurrentParameterValueNumber = atoi(argv[1]);
      int TotalNumberOfParameterValues= atoi(argv[2]);
      double  StartOfParameterValueRange = atof(argv[3]);
      double EndOfParameterValueRange=atof(argv[4]);
      double Kappa = atof(argv[5]); 
      double AverageDegree = atof(argv[6]);
      double Beta = atof(argv[7]);
      char RecordTimeSeries = *argv[8];
      int NumberOfGs = atoi(argv[9]); 
      string Outputs = argv[10];
      int n=pow(10.0,atof(argv[11])); 
      int Dimension = atoi(argv[12]);
      char VertexDistributionType = *argv[13];
      char ModelType = *argv[14];


      
      
      double CurrentParameterValue;
      
      if(TotalNumberOfParameterValues > 1)
      {
	    CurrentParameterValue = StartOfParameterValueRange + (EndOfParameterValueRange -StartOfParameterValueRange) * double(CurrentParameterValueNumber) /double( TotalNumberOfParameterValues -1);
      }
      else
      {
	    CurrentParameterValue = StartOfParameterValueRange;
      }
      
      
      
      string FileNameOutputPart="";
      if(Kappa == 0.0) //rewiring model
      {
	    if( Beta < 0.0)
	    {
		  FileNameOutputPart += "b";
		  Beta = CurrentParameterValue;
	    }
	    else
	    {
		  FileNameOutputPart += "mu";
		  AverageDegree = CurrentParameterValue;
	    }	    
      }
      else // percolation model
      {
	    if( Beta < 0.0)
	    {
		  FileNameOutputPart += "b";
		  Beta = CurrentParameterValue;
	    }
	    else
	    {
		  FileNameOutputPart += "k";
		  Kappa = CurrentParameterValue;
	    }	    
      }
      
      /*this part for \beta-c plot   beta from 0 to 10,  */ 
      string FileNameForThisPurpose = "MsD2mu5i"; FileNameForThisPurpose += argv[1]; FileNameForThisPurpose += ".dat";
      ifstream infile(FileNameForThisPurpose.c_str());
      
      infile>>Kappa;
      infile.close();
      /************************/
      
      
      FileNameOutputPart +="ucmuC";
      
      if(ModelType == 'o')
      {
	    FileNameOutputPart +="hrR";
      }
      
      string FileNameInputPart="";
      FileNameInputPart += "M"; FileNameInputPart += argv[14];
      FileNameInputPart += "n"; FileNameInputPart += argv[11]; FileNameInputPart += "mu"; FileNameInputPart += DoubleToString(AverageDegree);FileNameInputPart += "D"; FileNameInputPart += argv[12]; FileNameInputPart += "b"; 
      if(Beta < n) {FileNameInputPart += DoubleToString(Beta);} else {FileNameInputPart += "inf";}
      FileNameInputPart += "k"; FileNameInputPart += DoubleToString(Kappa); FileNameInputPart += argv[13];
      
      string FileName = "IN"+ FileNameInputPart + "OUT" + FileNameOutputPart;     
      
      ofstream myFile; 
      myFile.open(FileName.c_str());     
      
      cgraphstats ResultOfCurrentRun, SumOfResultsOfIndividualRuns, AverageOfResultsOfIndividualRuns;
      
      string FileNameRunStatus = "IN"+ FileNameInputPart + "OUT" + "RunStatus";
      
      for(int i=0; i<NumberOfGs;i++)
      {	    
	    ofstream myFileRunStatus; 
	    myFileRunStatus.open(FileNameRunStatus.c_str());
	    myFileRunStatus<<"Now executing graph number "<<i+1<<" out of "<<NumberOfGs<<"\n";
	    myFileRunStatus.close();
      
	    ResultOfCurrentRun =  SingleRun(n,Dimension, Kappa, AverageDegree, Beta, RecordTimeSeries == 'y',Outputs,VertexDistributionType, TotalNumberOfParameterValues,FileNameInputPart, ModelType);
	    
	    SumOfResultsOfIndividualRuns = SumOfResultsOfIndividualRuns + ResultOfCurrentRun;	    
      }
      
      AverageOfResultsOfIndividualRuns = SumOfResultsOfIndividualRuns/NumberOfGs;
      
      myFile<<CurrentParameterValue<<"\t"<< AverageOfResultsOfIndividualRuns.EnergyDensity<<"\t"<<AverageOfResultsOfIndividualRuns.FractionOfVerticesInLargestComponent<<"\t"<<AverageOfResultsOfIndividualRuns.AverageDegree ;
      
      myFile<<"\t"<<AverageOfResultsOfIndividualRuns.ClusteringCoefficient;
      if(ModelType == 'o')
      {
	    myFile<<"\t"<<AverageOfResultsOfIndividualRuns.AverageHopDistance<<"\t"<<AverageOfResultsOfIndividualRuns.AverageRouteDistance<<"\t"<<AverageOfResultsOfIndividualRuns.RouteFactor;
      }
      
      
      if(TotalNumberOfParameterValues < 10)
      {
	    string FileNameDegreeDist = "IN"+ FileNameInputPart + "OUT" + "DegreeDist";
	    string FileNameEdgeLengthDist = "IN"+ FileNameInputPart + "OUT" + "EdgeLengthDist";
	    string FileNameRobustness = "IN"+ FileNameInputPart + "OUT" + "Robustness";
	    
	    ofstream myFileDegreeDist, myFileEdgeLengthDist, myFileRobustness; 
	    myFileDegreeDist.open(FileNameDegreeDist.c_str());
	    myFileEdgeLengthDist.open(FileNameEdgeLengthDist.c_str());
	    myFileRobustness.open(FileNameRobustness.c_str());
	    
	    for(int i=0; i<AverageOfResultsOfIndividualRuns.FractionOfVerticesOfDegree.size(); i++)
	    {
		  myFileDegreeDist<<i<<"\t"<<AverageOfResultsOfIndividualRuns.FractionOfVerticesOfDegree[i]<<"\n"; 
	    }
	    for(int i=0; i<AverageOfResultsOfIndividualRuns.FractionOfEdgeLengthsAround.size(); i++)
	    {
		  myFileEdgeLengthDist<<i<<"\t"<<AverageOfResultsOfIndividualRuns.FractionOfEdgeLengthsAround[i]<<"\n"; 
	    }
	    
	    for(int i=0; i<AverageOfResultsOfIndividualRuns.Robustness.size(); i++)
	    {
		  myFileRobustness<<AverageOfResultsOfIndividualRuns.Robustness[i].first<<"\t"<<AverageOfResultsOfIndividualRuns.Robustness[i].second<<"\n"; 
	    }
	    
	    
	    myFileDegreeDist.close();
	    myFileEdgeLengthDist.close();
	    myFileRobustness.close();
      }
      
      myFile<<'\n';
      
      myFile.close();
      
      cout<<'\n';
      return 0;
}		     
