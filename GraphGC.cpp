#include <iostream>
#include <cmath>
#include "GraphGC.h"
using namespace std;

GraphGC::GraphGC(int n):vertices(n)
{
    geneM = new double*[vertices];
    phenoM = new double*[vertices];
    for (int r = 0; r< vertices; r++) {
        geneM[r] = new double[vertices];
        phenoM[r] = new double[vertices];
        for (int c = 0; c < vertices; c++) {
            geneM[r][c] = 0;
            phenoM [r][c] = 0;
        }
    }
}
GraphGC::GraphGC(const GraphGC& g)
{
    vertices = g.vertices;
    if ( g.vertices )
    {
        geneM = new double*[vertices];
        phenoM= new double*[vertices];
        for (int i = 0; i < vertices; i++) {
            geneM[i] = new double[vertices];
            phenoM[i]= new double[vertices];
            }

        for (int i = 0; i < vertices; i++)
        {
            for (int j = 0; j < vertices; j++)
            {
                geneM[i][j] = g.geneM[i][j];
                phenoM[i][j] = g.phenoM[i][j];
            }
        }
    }
}


GraphGC& GraphGC::operator=(const GraphGC& g)
{
    if (this != &g)
    {
        for (int i = 0; i < vertices; i++){
                  delete[] geneM[i];
                  delete[] phenoM[i];
                   }
            delete[] geneM;
            delete[] phenoM;

        vertices = g.vertices;

         geneM = new double*[vertices];
        phenoM= new double*[vertices];
        for (int i = 0; i < vertices; i++) {
            geneM[i] = new double[vertices];
            phenoM[i]= new double[vertices];
            }

        for (int i = 0; i < vertices; i++)
        {
            for (int j = 0; j < vertices; j++)
            {
                geneM[i][j] = g.geneM[i][j];
                phenoM[i][j] = g.phenoM[i][j];
            }
        }

    }
    return *this;
}
 void GraphGC::setW(int i, int j, double w=0) {
                  geneM[i][j] = geneM[i][j]+w;
                  geneM[j][i] = geneM[j][i]+w;
    }
 void GraphGC::setD(int i, int j, double d=0) {
                  phenoM[i][j] = d;
                  phenoM[j][i] = d;
    }

 void GraphGC::setWeight(int i, int j, double w=0) {
                  weights[i][j]=w;
    }
//double* GraphGC::getWeight() {
//    return &weight;
//    }

double GraphGC::getGW(int i, int j) {
                  return geneM[i][j];
    }
double GraphGC::getPW(int i, int j) {
                  return phenoM[i][j];
    }

double GraphGC::getProd(int i, int j) {
                  return ((phenoM[i][j]*geneM[i][j])+(phenoM[j][i]*geneM[j][i]));
    }
double GraphGC::SumProd (int N){
  double s=0;
  for (int i = 0; i < vertices; i++)
   {
    for (int ii = i+1; ii < vertices; ii++)
     {
      s=s+this->getProd(i, ii);
         }
         }
    return (s/N);
}
int GraphGC::GCPart(vector<int> GCassign, vector<int> GCVassign, int iv=0){  //iv=0 GCassign: invariant set of nodes, iv=1 variable part
   int LC;
   for (int i = 0; i < vertices; i++)
   {
    for (int ii = 0; ii < vertices; ii++)
     {
        if (iv==0) {
          phenoM[GCassign[i]][ii]=0;
          //phenoM[ii][GCassign[i]]=0;
          LC=GCassign.size();
                }
          else {
           phenoM[GCVassign[i]][ii]=0;
           //phenoM[ii][GCVassign[i]]=0;
           LC=GCVassign.size();
                 }
         }
         }
      return LC;
    }
vector< vector<int> > GraphGC::PairCompGC(vector<int> AAGCnum1, vector<int> AAGCnum2){
  vector< vector<int> > GComp;
  vector< int> GCommon, GCdiff;
  for (int i = 0; i < vertices; i++)
   {
        if (AAGCnum1[i]==AAGCnum2[i])
            GCommon.push_back(i);
        else
            GCdiff.push_back(i);
       }
       GComp.push_back(GCommon);
       GComp.push_back(GCdiff);
   return GComp;
    }
bool GraphGC::STOP(char a, char b) {
       return ((a=='*')||(b=='*'));
    }
bool GraphGC::METRP(int c, int cc, vector< float >& Rscv){
       return ((Rscv[c]==1)||(Rscv[cc]==1));
    }
bool GraphGC::transition(char a, char b) {
       return (((a=='T')&&(b=='C'))||((b=='T')&&(a=='C'))||((a=='G')&&(b=='A'))||((a=='A')&&(b=='G')));
    }
bool GraphGC::transversion(char a, char b) {
      return (((a=='T')&&(b=='A'))||((b=='T')&&(a=='A'))||((a=='G')&&(b=='T'))||((a=='T')&&(b=='G'))||((a=='C')&&(b=='A'))||((b=='C')&&(a=='A'))||((a=='G')&&(b=='C'))||((a=='C')&&(b=='G')));
    }

pair<int,char> GraphGC::AAnIdentify(int ss, string alp){
    int s=0;
    //string aa="ARN*DCQEGHILKMFPSTWYV";
    string aa64=alp;
    while ((s<21)&&(aa[s]!=aa64[ss])) s++;
    if (s<21)
          return make_pair(s, aa[s]);
    else  return make_pair(s,'0');
    }
vector<int> GraphGC::NumericRecode(string aa64){
    pair<int,char> po;
    vector<int> AAGCnum;
   for (int i = 0; i < aa64.size(); i++)
     {
         po=AAnIdentify(i, aa64);
         AAGCnum.push_back(po.first);
          }
 return AAGCnum;
 }

vector<double> GraphGC::P20ToP64(string gc, vector<double> Paa20){
  vector<double> Paa64;
  pair<int,char> aaID;
  for (int p = 0;p < gc.length(); p++) {
      aaID=AAnIdentify(p, gc);
      Paa64.push_back(Paa20[aaID.first]);
           }
  return Paa64;
     }
CodeN GraphGC::CodeNbd(vector<int> AAGCnum){
   vector<  vector<int> > L (vertices, vector<int>(6, 0));
   vector<int> B;
   int l,a,q,c;
for (int i = 0; i < vertices; i++)
   {
     l=0; q=0;
    for (int ii = 0; ii < vertices; ii++)
     {
       if ((i!=ii)&&(AAGCnum[i]==AAGCnum[ii])){
            c=0; a=0; q=q+1; int j = 0;
            for (j = 0; j< vertices; j++) {
                  a=0;
                  for (int t = 0; t< vertices; t++) {
                       if ((phenoM[i][j]>0)&&(phenoM[ii][t]>0)){
                           if ((phenoM[i][j]==phenoM[ii][t])&&(AAGCnum[j]==AAGCnum[t])){
                                           a=1; }
                                       }
                                      }
                                if (a==0) {
                                    break;
                                }//it is enough that one vertex j connected to i has any vertex t specifying the same aa as that for the vertex j ....and/or with the same weight on the edge connected to the vertex ii
                              }
                              if (j==vertices){
                                L[i][l]=ii;
                                l=l+1;
                                    }
                             }
                          }
                          if (l==0) L[i][l]=-1;
                         if ((q>0)&&(l==q)) B.push_back(1);
                         else B.push_back(0);

            }
         CodeN Cneighb;
         Cneighb.HB=B;
         Cneighb.HsB=L;
       return Cneighb;
     }
double GraphGC::MeanSuppresor(vector <double> Paa,vector< string > GCodes, string Stopcodon)
 {
  double sum1=0, sum2=0, sum3=0;
  int n=0;
    for (int i = 0; i < GCodes[0].size(); i++)
     {
      if (GCodes[0][i]!='*'){
       if ((GCodes[2][i]!=Stopcodon[0])&&(GCodes[3][i]==Stopcodon[1])&&(GCodes[4][i]==Stopcodon[2])){
                     sum1= sum1+Paa[i]; n=n+1;
             }
       if ((GCodes[3][i]!=Stopcodon[1])&&(GCodes[2][i]==Stopcodon[0])&&(GCodes[4][i]==Stopcodon[2])){
                     sum2= sum2+Paa[i]; n=n+1;
             }
       if ((GCodes[4][i]!=Stopcodon[2])&&(GCodes[3][i]==Stopcodon[1])&&(GCodes[2][i]==Stopcodon[0])){
                     sum3= sum3+Paa[i]; n=n+1;
             }
             }
    }
double out=(sum1+sum2+sum3)/n;
    return out;
 }

vector <double>  GraphGC::allMeanSuppresors(vector <double> Paa,vector< string > GCodes){
  string cod;
   for (int i = 0; i < GCodes[0].size(); i++)
     {
         if (GCodes[0][i]=='*'){
             string s1(1,GCodes[2][i]);
             string s2(1,GCodes[3][i]);
             string s3(1,GCodes[4][i]);
             cod=s1+s2+s3;
             Paa[i]=MeanSuppresor(Paa,GCodes,cod);
          }   }
  return Paa;
 }


GraphGC::GraphGC(const vector <double>  &Paa, const vector <string> &GCodes, int vert, int gc, int p, int s, int v){
  /*  p=0 s=0 v=1  transv model
    p=0 s=1 v=0    trans model
    p=3 s=1 v=1    3 pos model
    p=2 s=1 v=1    2 pos model
    p=1 s=1 v=1    1 pos model
    gc=0  block model
    gc=1  sense model
    gc=2  codon model
    */
  vertices=vert;//for sgc 64 61 20
  double d, da1, da2;
  vector <double> Paa64;
  pair<int,char> po, po2;
  string alp=GCodes[0], gcode=GCodes[0];
  Paa64=P20ToP64(gcode, Paa);
  Paa64=allMeanSuppresors(Paa64,GCodes);
 for (int i = 0; i < GCodes[0].length(); i++)
   {
    po=AAnIdentify(i, alp);
    for (int ii = 0; ii < GCodes[0].size(); ii++)
     {
          po2=AAnIdentify(ii, alp);
           if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])) {
                  da1=Paa[po.first];
                  da2=Paa[po2.first];
                  d=((da1-da2)*(da1-da2));
                  setD(i,ii,d);
                     }
             if (gc==2) {
                  da1=Paa64[i];
                  da2=Paa64[ii];
                  d=((da1-da2)*(da1-da2));
                  setD(i,ii,d);
               }
              if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))){
                  da1=Paa64[i];
                  da2=Paa64[ii];
                  d=((da1-da2)*(da1-da2));
                  setD(i,ii,d);
               }

        //3rst pos
       if ((p==3)||(p==0)&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){

           if ((s==1)&&(transition(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[2][0]);
                if (gc==2)
                    setW(i, ii, weights[2][0]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[2][0]);
                }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[2][1]);
                if (gc==2)
                    setW(i, ii, weights[2][1]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[2][1]);
                }
          }
       if (((p==2)||(p==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {

           if ((s==1)&&(transition(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[1][0]);
                if (gc==2)
                    setW(i, ii, weights[1][0]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[1][0]);
            }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[1][1]);
                if (gc==2)
                    setW(i, ii, weights[1][1]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[1][1]);
          }

               }
         if (((p==1)||(p==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {

           if ((s==1)&&(transition(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[0][0]);
                if (gc==2)
                    setW(i, ii, weights[0][0]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[0][0]);
            }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[0][1]);
                if (gc==2)
                    setW(i, ii, weights[0][1]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[0][1]);
          }
                     }

    }

}
}

vector <int>  GraphGC::BChanges(const vector <string> &GCodes, int vert, int gc, int p, int s, int v){ //gc=0 codon model, gc=1 block sense models

  vertices=vert;//for sgc 64 61 20
  //vector <double> Paa64;
  vector <int> ChangesCount;
  //Paa64=P20ToP64(GCodes[0], Paa);
  //int trans0p3=0, trans1p3=0, trans0p2=0, trans1p2=0, trans0p1=0, trans1p1=0;
  //int tranv0p3=0, tranv1p3=0, tranv0p2=0, tranv1p2=0, tranv0p1=0, tranv1p1=0;
  for (int i = 0; i < 12; i++)
     ChangesCount.push_back(0);
 for (int i = 0; i < GCodes[0].length(); i++)
   {

    for (int ii = 0; ii < GCodes[0].size(); ii++)
     {

        //3rst pos
       if (((p==3)||(p==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){

           if ((s==1)&&(transition(GCodes[4][i], GCodes[4][ii]))){
                if (gc==0) {ChangesCount[0]=ChangesCount[0]+1; }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))){
                    ChangesCount[1]=ChangesCount[1]+1;}
                }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if (gc==0)
                    ChangesCount[2]=ChangesCount[2]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    ChangesCount[3]=ChangesCount[3]+1;
                }
          }
       if (((p==2)||(p==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {
           if ((s==1)&&(transition(GCodes[4][i], GCodes[4][ii]))){
                if (gc==0)
                    ChangesCount[4]=ChangesCount[4]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    ChangesCount[5]=ChangesCount[5]+1;
            }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if (gc==0)
                     ChangesCount[6]=ChangesCount[6]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    ChangesCount[7]=ChangesCount[7]+1;
                                       }
               }
         if (((p==1)||(p==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {

           if ((s==1)&&(transition(GCodes[4][i], GCodes[4][ii]))){
                if (gc==0)
                      ChangesCount[8]=ChangesCount[8]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                     ChangesCount[9]=ChangesCount[9]+1;
            }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if (gc==0)
                    ChangesCount[10]=ChangesCount[10]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    ChangesCount[11]=ChangesCount[11]+1;
          }
                     }

    }
 }
return ChangesCount;
}
vector <int> GraphGC::AANumbers(vector <int> AAGC20, vector <int> AAGC64){
  vector <int> nc;
  int s=0;
  for (int i = 0; i < vertices; i++)
   {
    s=0;
    for (int ii = 0; ii < vertices; ii++)
     {
       if (AAGC20[i]==AAGC64[ii]) s=s+1;
       }
       nc.push_back(s);
    }
  return nc;
   }
double GraphGC::TNMean(int N){
double Sg=0, Sp=0;
double nmean;
for (int i = 0; i < vertices; i++)
   {
    for (int ii = 0; ii < vertices; ii++)
     {int
            Sg=Sg+geneM[i][ii];
            Sp=Sp+phenoM[i][ii];
                }
                }
     nmean=(Sg*Sp)/(N*(vertices*(vertices-1)) ) ;
 return nmean;
}

double GraphGC::TNVar(int N, double nmean){
 double nvar,s4=0, sg4=0, su=0,sug=0,sud=0, sudg=0, sudo2=0, sudo=0, sudog2=0, sudog=0, sumdo=0, sumdo4=0, sumdog=0, sumdog4=0;

 for (int i = 0; i < vertices; i++)  //
     {
      for (int ii = 0; ii < vertices; ii++)
       {
          if (ii!=i){
              su=su+phenoM[i][ii];  //daa2
              sug=sug+geneM[i][ii];
                            }
                }
          }
    for (int i = 0;i < vertices; i++)
     {
      for (int ii = i; ii < vertices; ii++)
       {
          if (ii!=i) {                                          //daa4
              sud=sud+(phenoM[i][ii]*phenoM[i][ii]);
              sudg=sudg+(geneM[i][ii]*geneM[i][ii]);
                                }

                            }
                }
     for (int i = 0; i < vertices; i++)  //T3   falta probar semi en sud y tot en sumdo ya lo deje listo
     {
       sudo=0;sudo2=0;
       sudog=0;sudog2=0;
      for (int ii = 0; ii < vertices; ii++)
       {
          if (ii!=i)
            {                                       //daa4
              sudo2=sudo2+(phenoM[i][ii]*phenoM[i][ii]);
              sudo=sudo+phenoM[i][ii];
              sudog2=sudog2+(geneM[i][ii]*geneM[i][ii]);
              sudog=sudog+geneM[i][ii];
                            }
                              }
               sumdo=sumdo+(sudo*sudo)-sudo2;  //T3
               sumdo4=sumdo4+(sudo*sudo);
               sumdog=sumdog+(sudog*sudog)-sudog2;  //T3
               sumdog4=sumdog4+(sudog*sudog);
                }
   s4=(su*su)-(4*(sumdo4))+2*sud;
   sg4=(sug*sug)-(4*(sumdog4))+2*sudg;
   nvar=(( ((s4*sg4)/((vertices)*(vertices-1)*(vertices-2)*(vertices-3)))+ ((sumdo*sumdog)/((vertices)*(vertices-1)*(vertices-2))+((sud*sudg)/((vertices)*(vertices-1)))  ))/(N*N)) - (nmean*nmean);
   return nvar;
}

double GraphGC::CUBound(double nvar, double nmean, double rob) {
    double CB;
    CB=nvar/(nvar-((rob-nmean)*(rob-nmean)));
    return CB;
}

 double GraphGC::score(double nvar, double nmean, double rob) {
    double sc;
    sc=(rob-nmean)/nvar;
    return sc;
}
vector< double> GraphGC::insert20(vector< double > Paa) {
               Paa.insert(Paa.begin() + 3, 0);
               return Paa;
                }
vector< double> GraphGC::standarize(vector< double> Paa){
           double m, st, sumc=0, sum=0;
           for (int i = 0; i < Paa.size(); i++){
                        sum=sum+Paa[i];
                  }
         for (int i = 0; i < Paa.size(); i++){
                        sumc=sumc+(Paa[i]*Paa[i]);
                  }
        m=(sum/Paa.size());
        st=sqrt((sumc/Paa.size())-(m*m));
        for (int i = 0; i < Paa.size(); i++){
                        Paa[i]=(Paa[i]-m)/st;
                  }
     return Paa;
  }

void GraphGC::print() {
      cout <<"Genetic Matrix"<<endl;
      for (int i = 0; i < vertices; i++) {
                  cout << i << " : ";
                  for (int j = 0; j < vertices; j++)
                        cout << geneM[i][j] << " ";
                  cout << "\n";
      }
      cout <<"Pheno-Matrix"<<endl;
      for (int i = 0; i < vertices; i++) {
                  cout << i << " : ";
                  for (int j = 0; j < vertices; j++)
                        cout << phenoM[i][j] << " ";
                  cout << "\n";
      }
    }


GraphGC::~GraphGC()
{
    for (int i = 0; i < vertices; i++){
                  delete[] geneM[i];
                  delete[] phenoM[i];
                   }
            delete[] geneM;
            delete[] phenoM;
 }
