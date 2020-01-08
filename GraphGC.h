#ifndef GRAPHGC_H
#define GRAPHGC_H
using namespace std;
#include <vector>

typedef struct{
   vector<  vector<int> > HsB;
   vector<int> HB;
}CodeN;

class GraphGC
{
    public:
        GraphGC(int n=0);
        GraphGC(const GraphGC& g);
        GraphGC(const vector <double> &, const vector <string> &, int, int, int, int, int);
        GraphGC& operator=(const GraphGC&);
        ~GraphGC();
       double getGW(int, int);
       double getPW(int, int);
       double getProd(int, int);
       //double* getWeight();
       bool STOP(char, char);
       bool METRP(int, int, vector< float >&);
       bool transition(char, char);
       bool transversion(char, char);
       vector <double> allMeanSuppresors(vector <double>,vector< string >);
       double SumProd(int);
       int GCPart(vector<int>, vector<int>, int);
       vector <int> BChanges(const vector <string>&, int, int, int, int, int);
       vector<int> NumericRecode(string);
       vector< vector <int> > PairCompGC(vector<int>, vector<int>);
       vector <int> AANumbers(vector <int>, vector <int>);
       CodeN CodeNbd(vector<int> AAGCnum);
       double TNMean(int);
       double TNVar(int, double);
       double CUBound(double, double, double);
       double score(double, double, double);
       vector< double> insert20(vector< double >);
       vector< double> standarize(vector< double>);
       void print();

    private:
        int vertices;
        double **geneM;
        double **phenoM;
        double weights[3][2] ={{1,0.5},{0.5,0.1},{1,1}};
        string aa="ARN*DCQEGHILKMFPSTWYV";
        void setW(int, int, double);
        void setD(int, int, double);
        void setWeight(int, int, double);
        pair<int,char> AAnIdentify(int, string);
        vector<double> P20ToP64(string, vector<double>);
        double MeanSuppresor(vector <double>,vector< string >, string);
};

#endif // GRAPHGC_H
