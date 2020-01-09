#ifndef SYNCOD_H
#define SYNCOD_H
using namespace std;
#include <vector>
#include <string>

struct Kodon {
            string id;
           double  fr;
             };

class SynCod
{
    public:
        SynCod(double f=0, string i="aaa");
        SynCod(const SynCod&);
        SynCod& operator=(const SynCod& other);
        int  posCod(vector< Kodon >&, char, char, char);
        vector< Kodon > preInCodon(string);
        vector <vector <double>  >  posInCodon(vector< Kodon >&, vector <string> &);
        pair< vector  <  vector<float> >, vector <char> >  posInCodonAA(vector< Kodon >, vector <string>&, string);
        vector < double > inCG(const string file);

    private:
        Kodon Kod;

};

#endif // SYNCOD_H
