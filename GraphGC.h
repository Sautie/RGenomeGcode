#ifndef GRAPHGC_H
#define GRAPHGC_H
using namespace std;
#include <vector>
#include <list>

typedef struct{
   vector<  vector<int> > HsB;
   vector<int> HB;
}CodeN;

class GraphGC
{
    public:
        GraphGC(int n=20);
        GraphGC(const GraphGC& g);
        GraphGC(vector <double> &, vector <string> &, int, int, int, int, int);
        GraphGC& operator=(const GraphGC&);
        ~GraphGC();
       double getGW(int, int);
       double getPW(int, int);
       double getProd(int, int);
       //double* getWeight();
         bool STOP(char, char);
       bool METRP(int, int, vector< double >);
       bool transition(char, char);
       bool transversion(char, char);
       void clearD();
       void clearW();
       vector <double> allMeanSuppresors(vector <double>,vector< string >);
       double SumProd(int);
       int GCPart(vector<int>, vector<int>, int);
       int PartChanges(vector<int> );
       vector<vector<vector<double> > > hPart(int , int, vector< vector<double> > , vector< vector<double> > ,  vector<int> ,  vector<int> ,  vector<int> );
       vector<vector<vector<double> > > gcPart(int, vector<int> , vector< vector<double> >, vector< vector<int> > , vector<int>, vector< vector<int> >,  vector< vector<int> >,  vector< vector<int> > );
       vector<vector<vector<double> > > mePart(int, vector<int> , double, vector< vector<double> >, vector< vector<int> >, vector< vector<int> > ,  vector<double>);
       vector<vector<vector<double> > > vaPart(vector<vector<vector<double> > >, vector<int>, int, double, vector< vector<double> >, vector< vector<int> >, vector< vector<int> >,  vector< vector<double> >);
       vector <int> BChanges(const vector <string>&, int, int, int, int, int);
       vector<int> NumericRecode(string);
       vector< vector <int> > PairCompGC(vector<int>, vector<int>);
       vector <int> AANumbers(vector <int>, vector <int>);
       CodeN CodeNbd(string AAGCnum);
       double TNMean(int);
       double TNVar(int, double);
       double CUBound(double, double, double);
       double score(double, double, double);
       vector< double> insert20(vector< double >);
       vector< double> standarize(vector< double>);
       void AP1(int, bool[], int [], int[], int [], bool[]);
       bool* APN();
       void toADJL();
       void print();
       bool Bconnected(string, char, char, vector< float >, vector< string >);
       vector <double> Bprob(ofstream&, int, vector< float >, vector <double>, vector< string >, vector <double>, vector <double>, char, char, int);
       int contAA(string, char);
       vector <double>  Reassign1(vector < vector <double> >, ofstream &, string, vector< float >, vector< string >, bool);
       vector< double > permGenCodes(string, const vector <double> &, const vector <double> &, vector <double> &, const vector <string>&);
       vector< double > permGenCodes2(string, const vector <double> &, const vector <double> &,  const vector <string> &);
       vector< double > permGenCodes64(string, const vector <double>  &, const vector <double>  &,  const vector <string> &);
    private:
        int vertices=64;
        double **geneM;
        double **phenoM;
        list<int> *adL;
        double weights[3][2] ={{1,0.5},{0.5,0.1},{1,1}};
        string aa="ARN*DCQEGHILKMFPSTWYV";
        void setW(int, int, double);
        void setD(int, int, double);
        void setWeight(int, int, double);
        pair<int,char> AAnIdentify(int, string);
        vector<double> P20ToP64(string, vector<double>);
        double MeanSuppresor(vector <double>,vector< string >, string);
        void AP1(int, bool, int, int, int, bool);
};

#endif // GRAPHGC_H
