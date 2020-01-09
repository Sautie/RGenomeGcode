#include <fstream>
#include <iostream>
#include "SynCod.h"
using namespace std;

SynCod::SynCod(double f, string i)
{
   Kod.fr=f;
   Kod.id=i;
}

SynCod::SynCod(const SynCod& oth)
{
    Kod.fr = oth.Kod.fr;
    Kod.id = oth.Kod.id;//copy ctor
}

SynCod& SynCod::operator=(const SynCod& h)
{
    if (this == &h) return *this; // handle self assignment
      Kod.fr = h.Kod.fr;
      Kod.id = h.Kod.id;
    return *this;
}


int  SynCod::posCod(vector< Kodon > &codons, char a, char b, char c){
          string cod;
             string s1(1,a);
             string s2(1,b);
             string s3(1,c);
             cod=s1 + s2+s3;
             int p=0;
             while ((p<codons.size())&&(codons[p].id!=cod)) p++;
        return p;
       }


vector< Kodon > SynCod::preInCodon(string file) {
    vector< Kodon > codons;
    Kodon cod;
       ifstream ff(file.c_str());
        if (!ff.is_open())
         { cout << "Error opening file"; }
        string line;
       // cout <<file<<endl
        int pff=0, ppf, pp, pf;
        while(getline(ff,line))
        {
            pff=0;
          if (line.length()>0) {
            int c=0;
            while(c<4){
             ppf = line.find_first_of("ACGT",pff);
             pp = line.find_first_not_of("ACGT",ppf);
             cod.id=line.substr(ppf,pp-ppf);
             pf = line.find_first_of('(',pp);
             pff = line.find_first_of(')',pf);
             cod.fr=stof(line.substr(pf+1,pff-pf-1));
             codons.push_back(cod);
            c++;
            //cout<<stoi(line.substr(pf+1,pff-pf-1))<<endl;                               }
                                     }
                   }
                       }
        return codons;
    }

  vector <vector <double>  >  SynCod::posInCodon(vector< Kodon >& codons, vector <string> &GCodes) {
       vector <vector <double>  > outres;
       int h=0, hh;
        vector < string >::iterator itstr, itstr2;
        vector <double> Fcod, Faa, Rscv, FSyn;
        int tot=0;
        int t=0;
        while (t<codons.size()){
                tot=tot+codons[t].fr; t++;
               }

        while (h<GCodes[0].size()){
             string cod;
             double su=0, ns=0;
             int p=0;

             string s1(1,GCodes[2][h]);
             string s2(1,GCodes[3][h]);
             string s3(1,GCodes[4][h]);
             cod=s1 + s2+s3;

             while ((p<codons.size())&&(codons[p].id!=cod)) p++;
             if (p<codons.size()) {
                hh=0;su=0;ns=0;
                while(hh<GCodes[0].size()){
                   if (GCodes[0][h]==GCodes[0][hh]) {
                     int pp=0;
                     string s1(1,GCodes[2][hh]);
                     string s2(1,GCodes[3][hh]);
                     string s3(1,GCodes[4][hh]);
                     cod=s1 + s2+s3;

                     while ((pp<codons.size())&&(codons[pp].id!=cod)) pp++;
                     su=su+(codons[pp].fr);
                     ns=ns+1;
                          }
                    hh++;
                    }
                  //  cout<<endl;
                Rscv.push_back(double(codons[p].fr)/(double(su/ns)) );
                Fcod.push_back(double(codons[p].fr)/double(tot));
                FSyn.push_back(double(codons[p].fr)/double(su));
                Faa.push_back(double(su)/double(tot));

                      }
                ++h;
              }
          outres.push_back(Rscv);
          outres.push_back(Fcod); //1
          outres.push_back(FSyn);  //2
          outres.push_back(Faa);  //3
        return outres;
     }

pair< vector  <  vector<float> >, vector <char> >  SynCod::posInCodonAA(vector< Kodon > codons, vector <string> &GCodes, string file) {
      //freq relativa de codones
       int h=0, hh;
       vector  <  vector<float> > Prod;
        vector < string >::iterator itstr, itstr2;
        vector <float> Fcodons, aa, Fcod, FSyn;
        vector <char> alphAA;
         string ex, name1, name2, name3, dest2; //file: name of codon usage file
         dest2="CodBias";
         name1=dest2+file;
         dest2="aa";
         name2=dest2+file;
         dest2="CodFr";
         name3=dest2+file;

        int sum=0, c=0, f;
        h=0;
        while (h<codons.size()) { sum=sum+codons[h].fr; h++;} //
        h=0;
        while (h<GCodes[0].size()){
             string cod;
             float su=0;
             int p=0;f=0;
             if (h==0) { alphAA.push_back(GCodes[0][h]); c=1; }
             else {
                    int n=0;
                    while(n<alphAA.size()&&alphAA[n]!=GCodes[0][h]) n++;
                    if (n==alphAA.size()) {alphAA.push_back(GCodes[0][h]); c=1;}
                    else c=0;
                      }
             p=posCod(codons, GCodes[2][h],GCodes[3][h], GCodes[4][h]); //pos de codon en codons
             if (p<codons.size()) {
                hh=0;su=0; //ns=0;
                while(hh<GCodes[0].size()){
                   if (GCodes[0][h]==GCodes[0][hh]) {
                     int pp=0;
                     pp=posCod(codons, GCodes[2][hh],GCodes[3][hh], GCodes[4][hh]);
                     su=su+(codons[pp].fr);
                    // ns=ns+1;
                          }
                    hh++;
                    }

                if (c==1) {aa.push_back(float(su)/(float(sum))); }
                FSyn.push_back(float(codons[p].fr)/float(su));
                Fcod.push_back(float(codons[p].fr)/float(sum));
                      }
                ++h;
              }
         Prod.push_back(aa);
         Prod.push_back(FSyn);
         Prod.push_back(Fcod);
        return make_pair(Prod, alphAA);
     }

vector < double > SynCod::inCG(const string file)  {
    vector< string > codons;
    vector< int > Freqodons;
    ifstream ff(file.c_str());
      if (!ff.is_open())
        { cout << "Error opening file"; }
        string line;
        int pff=0, ppf, pp, pf;
        int sum=0;
        while(getline(ff,line))  //captura de datos
        {
            pff=0;
            if (line.length()>0) {
            int c=0;
            while(c<4){
                    //cout<<" "<<line.substr(ppf,pp-ppf)<<" ";
             ppf = line.find_first_of("ACGT",pff);
             pp = line.find_first_not_of("ACGT",ppf);
             codons.push_back(line.substr(ppf,pp-ppf));
             pf = line.find_first_of('(',pp);
             pff = line.find_first_of(')',pf);
             Freqodons.push_back(stoi(line.substr(pf+1,pff-pf-1)));
             sum=sum+stoi(line.substr(pf+1,pff-pf-1));
            c++;
                            }
                        }
                   }

        int c1=0, c2=0, c3=0;
        for (unsigned int r = 0; r<codons.size(); ++r) {
             string kodon=codons[r];
             if ((kodon[2]=='G')||(kodon[2]=='C')) c3=c3+Freqodons[r]; //pos 3 codon
             if ((kodon[1]=='G')||(kodon[1]=='C')) c2=c2+Freqodons[r]; //pos 2 codon
             if ((kodon[0]=='G')||(kodon[0]=='C')) c1=c1+Freqodons[r]; //pos 2 codon
                }
        vector < double > gc;
        gc.push_back(double(double(c1)/double(sum)));
        gc.push_back(double(double(c2)/double(sum)));
        gc.push_back(double(double(c3)/double(sum)));
      return gc;
        }
