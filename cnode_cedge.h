//#include<unordered_set>
#include<vector>
using namespace std;



class cnode{
public:
      double* Location;
      vector<int> Neighbors;
};


class cedge{
public:
      int SmallVertex;
      int LargeVertex;
      cedge();
	cedge(int, int, double );
	cedge(int, int);
      bool operator==(cedge);
      void operator=(cedge);
      bool operator<(cedge);
      void Display();
      double Length;
};

class cgraphstats{
public:
      cgraphstats();
      cgraphstats(double , double , double , double, double, vector<double> , vector<double> ,vector<pair<double, double> >,double,double);
      double EnergyDensity;
      double AverageDegree;
      double FractionOfVerticesInLargestComponent;
      double AverageHopDistance;
      double AverageRouteDistance;
      double ClusteringCoefficient;
      double RouteFactor;
      
      vector<double> FractionOfVerticesOfDegree;
      vector<double> FractionOfEdgeLengthsAround;
      vector<pair<double,double> > Robustness;
      cgraphstats operator+ (cgraphstats);
      cgraphstats operator/ (int);
      void Display();
};

bool operator<(const cedge& lhs, const cedge& rhs );

cgraphstats SingleRun(int, double,double , double , double , bool ,string, char,int,string, char);
cgraphstats SingleRun_data(string , string ,double , double , bool ,string , int , string );
int partition(cedge* , int , int );

cedge quick_select(cedge* , int , int , int );
