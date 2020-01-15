#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <set>
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
 void GraphGC::clearD() {
   for (int i = 0; i < vertices; i++){
    for (int ii = 0; ii < vertices; ii++) phenoM[i][ii] = 0;
     }
  }
void GraphGC::clearW() {
   for (int i = 0; i < vertices; i++){
    for (int ii = 0; ii < vertices; ii++) geneM[i][ii] = 0;
     }
  }
void GraphGC::toADJL()
 {
    for (int i = 0; i < vertices; i++){
      for (int ii = i+1; ii < vertices; ii++) {
        if (geneM[i][ii]>0){
            adL[i].push_back(ii);
            adL[ii].push_back(i);
                      }
        }
     }
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
bool GraphGC::METRP(int c, int cc, vector< double > Rscv){
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
    // string aa64=alp;
    while ((s<21)&&(aa[s]!=alp[ss])) s++;
    if (s<21)
          return make_pair(s, aa[s]);
    else  return make_pair(s,'0');
    }
vector<int> GraphGC::NumericRecode(string aa64){
    pair<int,char> po;
    vector<int> AAGCnum;
    for (int i = 0; i < aa64.size(); i++)
     {   po=AAnIdentify(i, aa64);
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

CodeN GraphGC::CodeNbd(string AAGC){
   vector<  vector<int> > L (vertices, vector<int>(6, 0));
   vector<int> B;
   int l,j,a,q,c;
for (int i = 0; i < vertices; i++)
   {
     q=0;
    for (int ii = 0; ii < vertices; ii++)
     {
       if ((i!=ii)&&(AAGC[i]==AAGC[ii])){
             l=0;a=0; q=q+1;
            for (j = 0; j< vertices; j++) {
                  a=0; c=0;
                  if (geneM[i][j]>0) {
                  for (int t = 0; t< vertices; t++) {
                       if (geneM[ii][t]>0){
                           if ((j!=t)&&(geneM[i][j]==geneM[ii][t])&&(AAGC[j]==AAGC[t])){
                                           a=1;
                                          // cout<<i;//cout<<j<<" "<<t<<" "<<vertices<<endl;}
                                           }
                                       }
                                      }
                                     // cout<<" "<<a<<endl;
                                if (a==0) {
                                    break;
                                   }
                                 }
                              }
                              cout<<j<<endl;
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


GraphGC::GraphGC(vector <double>  &Paa, vector <string> &GCodes, int vert, int gc, int p, int s, int v){
  /*  p=0 s=0 v=1  transv model
    p=0 s=1 v=0    trans model
    p=3 s=1 v=1    3 pos model
    p=2 s=1 v=1    2 pos model
    p=1 s=1 v=1    1 pos model         vert 64 gc 2 p=1 s=1 v=1
    p=0 s=1 v=1     tot model
    gc=0  block model
    gc=1  sense model
    gc=2  codon model
    */
  vertices=vert;//for sgc 64 61 20
  double d, da1, da2;
  vector <double> Paa64;
  pair<int,char> po, po2;
  string alp=GCodes[0], gcode=GCodes[0];
  Paa=insert20(Paa);
  Paa64=P20ToP64(alp, Paa);
  Paa64=allMeanSuppresors(Paa64,GCodes);

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

           if ((s==1)&&(transition(GCodes[3][i], GCodes[3][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[1][0]);
                if (gc==2)
                    setW(i, ii, weights[1][0]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[1][0]);
            }
          if ((v==1)&&(transversion(GCodes[3][i], GCodes[3][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[1][1]);
                if (gc==2)
                    setW(i, ii, weights[1][1]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[1][1]);
          }

               }
         if (((p==1)||(p==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {

           if ((s==1)&&(transition(GCodes[2][i], GCodes[2][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii]))
                    setW(i, ii, weights[0][0]);
                if (gc==2)
                    setW(i, ii, weights[0][0]);
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    setW(i, ii, weights[0][0]);
            }
          if ((v==1)&&(transversion(GCodes[2][i], GCodes[2][ii]))){

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
           if ((s==1)&&(transition(GCodes[3][i], GCodes[3][ii]))){
                if (gc==0)
                    ChangesCount[4]=ChangesCount[4]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    ChangesCount[5]=ChangesCount[5]+1;
            }
          if ((v==1)&&(transversion(GCodes[3][i], GCodes[3][ii]))){
                if (gc==0)
                     ChangesCount[6]=ChangesCount[6]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                    ChangesCount[7]=ChangesCount[7]+1;
                                       }
               }
         if (((p==1)||(p==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {

           if ((s==1)&&(transition(GCodes[2][i], GCodes[2][ii]))){
                if (gc==0)
                      ChangesCount[8]=ChangesCount[8]+1;
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii])))
                     ChangesCount[9]=ChangesCount[9]+1;
            }
          if ((v==1)&&(transversion(GCodes[2][i], GCodes[2][ii]))){
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
     {
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

void GraphGC::AP1(int j, bool vis[], int dis[], int l[], int par[], bool ap[]) {
   static int t=0;
   int child = 0;
   vis[j] = true;
   dis[j] = l[j] = ++t;
   list<int>::iterator i;
   for (i = adL[j].begin(); i != adL[j].end(); ++i) {
      int ii = *i;
      if (!vis[ii]) {
         child++;
         par[ii] = j;
         AP1(ii, vis, dis, l, par, ap);
         l[j] = (l[j] < l[ii]) ? l[j] : l[ii];
         if (par[j] == -1 && child> 1)
            ap[j] = true;
         if (par[j] != -1 && l[ii] >= dis[j])
            ap[j] = true;
      } else if (ii != par[j])
         l[j] = min(l[j], dis[ii]);
   }
}

bool* GraphGC::APN() {
   bool *vis = new bool[vertices];
   int *dis = new int[vertices];
   int *l = new int[vertices];
   int *par = new int[vertices];
   bool *ap = new bool[vertices];
   for (int i = 0; i < vertices; i++) {
      par[i] = -1;
      vis[i] = false;
      ap[i] = false;
   }
   for (int i = 0; i < vertices; i++)
      if (vis[i] == false)
         AP1(i, vis, dis, l, par, ap);
     return ap;
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

vector< double > GraphGC::permGenCodes(string alp, const vector <double>  &CodBias, const vector <double>  &Faas, vector <double>  &Paa, const vector <string> &GCodes) {

  int codonpos=0;
  pair<int,char> po, po2;
  double spaa=0,sp1=0, sp2=0, sp3=0, sumt1=0, sumt2=0,sumt3=0;
  double sumd1=0, sumd2=0,sumd3=0;

  int w, mt, nt=0, nt1=0, nt2=0, nt3=0;
  int nts1=0, nts2=0, nts3=0;
  int ntt=0, ntt1=0, ntt2=0, ntt3=0;
  int ntts1=0, ntts2=0, ntts3=0;

   double transd=0, transf=0;
   double tranvd=0,tranvf=0;
   int trans=0, tranv=0;
  int transs=0, tranvv=0;

  double dstop=0, sumt2stop=0, sumt1stop=0, sumt3stop=0, sp1stop=0, sp2stop=0,sp3stop=0;
  double spaaStop=0, transfStop=0, transdStop=0, tranvfStop=0, tranvdStop=0;

    set < int > aaconnected, aasconnected;
  int saa=0, sas=0;
  double d, da1, da2;
  vector <double> Paa64;
  Paa=insert20(Paa);
  Paa64=P20ToP64(GCodes[0], Paa);

 for (unsigned int i = 0; i < GCodes[0].size(); i++)
   {
    po=AAnIdentify(i, alp);
    for (unsigned int ii = 0; ii < GCodes[0].size(); ii++)
     {
              po2=AAnIdentify(ii, alp);

              da1=Paa[po.first];
              da2=Paa[po2.first];
              dstop=((da1-da2)*(da1-da2));

    if (GCodes[0][i]!=GCodes[0][ii]) {

       if (((codonpos==3)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){
                  if (transition(GCodes[4][i], GCodes[4][ii])){
                        if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *  model cg20
                        else {
                            d=((Paa[po.first]-Paa[po2.first])*(Paa[po.first]-Paa[po2.first]));
                            trans=trans+1;
                            } //model cg20
                       if (METRP(i, ii, CodBias)) mt=0;
                          else mt=1;

                       sumt3=sumt3+((CodBias[i])*d*mt);  //(this->weights[2][0])*
                       sp3=sp3+((weights[2][0])*(CodBias[i])*d*mt);
                       spaa=spaa+(weights[2][0])*Faas[i]*d*mt;
                        transf=transf+((CodBias[i])*d*mt); //******
                       transd=transd+((weights[2][0])*(CodBias[i])*d*mt); //******

                       sumt3stop=sumt3stop+((CodBias[i])*dstop*mt);  //(this->weights[2][0])*
                       sp3stop=sp3stop+((weights[2][0])*(CodBias[i])*dstop*mt);//
                       spaaStop=spaaStop+(weights[2][0])*Faas[i]*dstop*mt;
                        transfStop=transfStop+((CodBias[i])*dstop*mt); //******
                       transdStop=transdStop+((weights[2][0])*(CodBias[i])*dstop*mt); //*******

                       transs=transs+1;

                                           }
                  if (transversion(GCodes[4][i], GCodes[4][ii])){
                      if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                        else {
                                d=((Paa[po.first]-Paa[po2.first])*(Paa[po.first]-Paa[po2.first]));
                                tranv=tranv+1;
                                         }
                     if (METRP(i, ii, CodBias)) mt=0;
                          else mt=1;

                       sumt3=sumt3+((CodBias[i])*d*mt);
                       sp3=sp3+((weights[2][1])*(CodBias[i])*d*mt);
                        spaa=spaa+(weights[2][1])*Faas[i]*d*mt;
                        tranvf=tranvf+((CodBias[i])*d*mt); //******
                       tranvd=tranvd+((weights[2][1])*(CodBias[i])*d*mt); //******
                       tranvv=tranvv+1;

                     sumt3stop=sumt3stop+((CodBias[i])*dstop*mt);
                       sp3stop=sp3stop+((weights[2][1])*(CodBias[i])*dstop*mt);//
                       spaaStop=spaaStop+(weights[2][1])*Faas[i]*dstop*mt;
                        tranvfStop=tranvfStop+((CodBias[i])*dstop*mt); //******
                       tranvdStop=tranvdStop+((weights[2][1])*(CodBias[i])*dstop*mt); //****/
                                                      }
                      nts3=nts3+1;
                      if (!(STOP(GCodes[0][i], GCodes[0][ii]))) {
                        aaconnected.insert(po2.second);
                         nt3=nt3+1;
                                }
                      aasconnected.insert(po2.second);
              }
       if (((codonpos==2)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {

                 if (transition(GCodes[3][i], GCodes[3][ii])){
                     if (STOP(GCodes[0][i], GCodes[0][ii])) d=0;
                        else {
                            d=((Paa[po.first]-Paa[po2.first])*(Paa[po.first]-Paa[po2.first]));
                            trans=trans+1; }
                     if (METRP(i, ii, CodBias)) mt=0;
                          else mt=1;

                       sumt2=sumt2+((CodBias[i])*d*mt);
                       sp2=sp2+((weights[1][0])*(CodBias[i])*d*mt);
                       spaa=spaa+(weights[1][0])*Faas[i]*d*mt;
                       transf=transf+((CodBias[i])*d*mt); //******
                       transd=transd+((weights[1][0])*(CodBias[i])*d*mt); //******
                       transs=transs+1;

                       sumt2stop=sumt2stop+((CodBias[i])*dstop*mt);  //(this->weights[2][0])*
                       sp2stop=sp2stop+((weights[1][0])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[1][0])*Faas[i]*dstop*mt;
                        transfStop=transfStop+((CodBias[i])*dstop*mt); //******
                       transdStop=transdStop+((weights[1][0])*(CodBias[i])*dstop*mt); //******/
                               }
                 if (transversion(GCodes[3][i], GCodes[3][ii])){
                    if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                        else {
                                d=((Paa[po.first]-Paa[po2.first])*(Paa[po.first]-Paa[po2.first]));
                                tranv=tranv+1;
                                         }
                         if (METRP(i, ii, CodBias)) mt=0;
                          else mt=1;

                         sumt2=sumt2+((CodBias[i])*d*mt);
                         sp2=sp2+((weights[1][1])*(CodBias[i])*d*mt);
                           spaa=spaa+(weights[1][1])*Faas[i]*d*mt;
                        tranvf=tranvf+((CodBias[i])*d*mt); //******
                       tranvd=tranvd+((weights[1][1])*(CodBias[i])*d*mt); //******
                       tranvv=tranvv+1;

                      sumt2stop=sumt2stop+((CodBias[i])*dstop*mt);
                       sp2stop=sp2stop+((weights[1][1])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[1][1])*Faas[i]*dstop*mt;
                        tranvfStop=tranvfStop+((CodBias[i])*dstop*mt); //******
                       tranvdStop=tranvdStop+((weights[1][1])*(CodBias[i])*dstop*mt); //******/
                               }
                  nts2=nts2+1;
                  if (!(STOP(GCodes[0][i], GCodes[0][ii]))) {
                      aaconnected.insert(po2.second); //cg20
                    nt2=nt2+1;
                                      }
                 aasconnected.insert(po2.second); //cg21
             }
       if (((codonpos==1)||(codonpos==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {
                   if (transition(GCodes[2][i], GCodes[2][ii])){
                     if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                        else {
                            d=((Paa[po.first]-Paa[po2.first])*(Paa[po.first]-Paa[po2.first]));
                            trans=trans+1; }
                      if (METRP(i, ii, CodBias)) mt=0;
                          else mt=1;
                         sumt1=sumt1+((CodBias[i])*d*mt);
                       sp1=sp1+((weights[0][0])*(CodBias[i])*d*mt);
                        spaa=spaa+(weights[0][0])*Faas[i]*d*mt;
                        transf=transf+((CodBias[i])*d*mt); //******
                       transd=transd+((weights[0][0])*(CodBias[i])*d*mt); //******
                       transs=transs+1;

                      sumt1stop=sumt1stop+((CodBias[i])*dstop*mt);
                       sp1stop=sp1stop+((weights[0][0])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[0][0])*Faas[i]*dstop*mt;
                       transfStop=transfStop+((CodBias[i])*dstop*mt); //******
                       transdStop=transdStop+((weights[0][0])*(CodBias[i])*dstop*mt); //******/

                        }
                 if (transversion(GCodes[2][i], GCodes[2][ii])){
                     if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                        else {
                                d=((Paa[po.first]-Paa[po2.first])*(Paa[po.first]-Paa[po2.first]));
                                tranv=tranv+1;
                                         }
                      if (METRP(i, ii, CodBias)) mt=0;
                          else mt=1;
                       sumt1=sumt1+((CodBias[i])*d*mt);
                       sp1=sp1+((weights[0][1])*(CodBias[i])*d*mt);
                       spaa=spaa+(weights[0][1])*Faas[i]*d*mt;
                        tranvf=tranvf+((CodBias[i])*d*mt); //******
                       tranvd=tranvd+((weights[0][1])*(CodBias[i])*d*mt); //******
                       tranvv=tranvv+1;

                       sumt1stop=sumt1stop+((CodBias[i])*dstop*mt);
                       sp1stop=sp1stop+((weights[0][1])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[0][1])*Faas[i]*dstop*mt;
                       tranvfStop=tranvfStop+((CodBias[i])*dstop*mt); //******
                       tranvdStop=tranvdStop+((weights[0][1])*(CodBias[i])*dstop*mt); //*****/

                    }
                  nts1=nts1+1;
                    if (!(STOP(GCodes[0][i], GCodes[0][ii]))) {
                         aaconnected.insert(po2.second);
                        nt1=nt1+1;
                                          }
                       aasconnected.insert(po2.second);
                }
            }
      else if ((GCodes[0][i]==GCodes[0][ii])&&(i!=ii)) {
         if (METRP(i, ii, CodBias)) mt=0;
                          else mt=1;
       if (((codonpos==3)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){
                 ntts3=ntts3+1;  //cg21 cg19

                 if (!(STOP(GCodes[0][i], GCodes[0][ii]))){
                       if (transversion(GCodes[4][i], GCodes[4][ii])) tranv=tranv+1;
                       if (transition(GCodes[4][i], GCodes[4][ii]))  trans=trans+1;
                        ntt3=ntt3+1; //cg20 cg18
                        }//nor * //nor * and nor *
                        if (transversion(GCodes[4][i], GCodes[4][ii])){
                            tranvv=tranvv+1;
                            sumt3stop=sumt3stop+((CodBias[i])*dstop*mt);  //(this->weights[2][0])*
                            sp3stop=sp3stop+((weights[2][1])*(CodBias[i])*dstop*mt);//
                            spaaStop=spaaStop+(weights[2][1])*Faas[i]*dstop*mt;
                            tranvfStop=tranvfStop+((CodBias[i])*dstop*mt); //******
                             tranvdStop=tranvdStop+((weights[2][1])*(CodBias[i])*dstop*mt); //****/
                             }
                       if (transition(GCodes[4][i], GCodes[4][i])){
                            transs=transs+1;
                            sumt3stop=sumt3stop+((CodBias[i])*dstop*mt);
                            sp3stop=sp3stop+((weights[2][0])*(CodBias[i])*dstop*mt);//
                            spaaStop=spaaStop+(weights[2][0])*Faas[i]*dstop*mt;
                            transfStop=transfStop+((CodBias[i])*dstop*mt); //******
                            transdStop=transdStop+((weights[2][0])*(CodBias[i])*dstop*mt); //****
                       }
        }
       if (((codonpos==2)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {
                ntts2=ntts2+1;

                if (!(STOP(GCodes[0][i], GCodes[0][ii]))){
                       if (transversion(GCodes[3][i], GCodes[3][ii])) tranv=tranv+1;
                       if (transition(GCodes[3][i], GCodes[3][ii]))  trans=trans+1;
                        ntt2=ntt2+1;  }
                        //nor *  //nor * and nor *
                        if (transversion(GCodes[3][i], GCodes[3][ii])) {
                             tranvv=tranvv+1;
                                    sumt2stop=sumt2stop+((CodBias[i])*dstop*mt);
                       sp2stop=sp2stop+((weights[1][1])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[1][1])*Faas[i]*dstop*mt;
                        tranvfStop=tranvfStop+((CodBias[i])*dstop*mt); //******
                       tranvdStop=tranvdStop+((weights[1][1])*(CodBias[i])*dstop*mt);
                                }
                       if (transition(GCodes[3][i], GCodes[3][ii])) {
                        transs=transs+1;
                       sumt2stop=sumt2stop+((CodBias[i])*dstop*mt);
                       sp2stop=sp2stop+((weights[1][0])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[1][0])*Faas[i]*dstop*mt;
                        transfStop=transfStop+((CodBias[i])*dstop*mt); //******
                       transdStop=transdStop+((weights[1][0])*(CodBias[i])*dstop*mt); //******/
                           }
                    }
       if (((codonpos==1)||(codonpos==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {
                // cout<<po2.first<<" ";
                 ntts1=ntts1+1;
                 if (!(STOP(GCodes[0][i], GCodes[0][ii]))){
                     if (transversion(GCodes[2][i], GCodes[2][ii])) tranv=tranv+1;
                         if (transition(GCodes[2][i], GCodes[2][ii]))  trans=trans+1;
                        ntt1=ntt1+1;
                          }
                        //nor *  //nor * and nor *
                        if (transversion(GCodes[2][i], GCodes[2][ii]))  {
                            tranvv=tranvv+1;
                             sumt1stop=sumt1stop+((CodBias[i])*dstop*mt);
                       sp1stop=sp1stop+((weights[0][1])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[0][1])*Faas[i]*dstop*mt;
                       tranvfStop=tranvfStop+((CodBias[i])*dstop*mt); //******
                       tranvdStop=tranvdStop+((weights[0][1])*(CodBias[i])*dstop*mt); //*****/
                            }
                       if (transition(GCodes[2][i], GCodes[2][ii])) {
                        transs=transs+1;
                        sumt1stop=sumt1stop+((CodBias[i])*dstop*mt);
                       sp1stop=sp1stop+((weights[0][0])*(CodBias[i])*dstop*mt);//*
                       spaaStop=spaaStop+(weights[0][0])*Faas[i]*dstop*mt;
                       transfStop=transfStop+((CodBias[i])*dstop*mt); //******
                       transdStop=transdStop+((weights[0][0])*(CodBias[i])*dstop*mt); //******/
                       }
                }
            }
        }
          saa=saa+aaconnected.size();
          sas=sas+aasconnected.size();
          aaconnected.clear();
          aasconnected.clear();
      }
      vector< double > res;

      res.push_back((sumt1+sumt2+sumt3)/(nt1+nt2+nt3+ntt1+ntt2+ntt3));
      res.push_back((sp1+sp2+sp3)/(nt1+nt2+nt3+ntt1+ntt2+ntt3));
      res.push_back(spaa/(nt1+nt2+nt3+ntt1+ntt2+ntt3));
      res.push_back(sumt1/(nt1+ntt1));
      res.push_back(sumt2/(nt2+ntt2));
      res.push_back(sumt3/(nt3+ntt3));
      res.push_back(sp1/(nt1+ntt1));
      res.push_back(sp2/(nt2+ntt2));
      res.push_back(sp3/(nt3+ntt3));
      res.push_back(transf/trans);
      res.push_back(transd/trans);
      res.push_back(tranvf/tranv);
      res.push_back(tranvd/tranv);


      res.push_back(saa);
      res.push_back(nt1+nt2+nt3);
      res.push_back(nt1);
      res.push_back(nt2);
      res.push_back(nt3);
      res.push_back((nt1+nt2+nt3)+(ntt1+ntt2+ntt3));

       res.push_back(nt1+ntt1);
       res.push_back(nt2+ntt2);
       res.push_back(nt3+ntt3);

      res.push_back(sas);
      res.push_back(nts1+nts2+nts3);
      res.push_back(nts1);
      res.push_back(nts2);
      res.push_back(nts3);
      res.push_back((nts1+nts2+nts3)+(ntts1+ntts2+ntts3));

       res.push_back(nts1+ntts1);
       res.push_back(nts2+ntts2);
       res.push_back(nts3+ntts3);

       res.push_back(transd/trans);
       res.push_back(transf/trans);
       res.push_back(tranvd/tranv);
       res.push_back(tranvf/tranv);

   return res;
 }
 bool GraphGC::Bconnected(string gc, char aa, char ab, vector< float > Paa, vector< string > GCodes){
    pair<int,char> po, po2;
    bool connect=true;
    vector < int > resu;
     unsigned int g=0;
    while(g<gc.length()){
     po=AAnIdentify(g, gc);
     for (unsigned int ii = 0; ii < gc.length(); ii++)
      {
       if ((gc[g]==aa)&&(gc[ii]==ab)) {
         po2=AAnIdentify(ii, gc);
                if ((GCodes[2][g]==GCodes[2][ii])&&(GCodes[3][g]==GCodes[3][ii])&&(GCodes[4][g]!=GCodes[4][ii])){
                    return connect;
                                }
               if ((GCodes[2][g]==GCodes[2][ii])&&(GCodes[3][g]!=GCodes[3][ii])&&(GCodes[4][g]==GCodes[4][ii])){
                    return connect;
                             }
                if ((GCodes[2][g]!=GCodes[2][ii])&&(GCodes[3][g]==GCodes[3][ii])&&(GCodes[4][g]==GCodes[4][ii])){
                    return connect;
                         }
                    }
              ii++;
                 }
          g++;
       }
   connect=false;
return connect;
 }

 vector <double> GraphGC::Bprob(ofstream &InputFile, int h , vector< float > Paa, vector <double> vg, vector< string > GCodes, vector <double>  resy2, vector <double> resy, char aa, char ab, int g)  {
  bool selecton;
  string change;
  int L=9;
  vector < vector <double> > freqs;
    for (unsigned int ii = 0; ii < L; ii++)
     {
        if (g==1) selecton=(resy[ii]<resy2[ii]); //cantelli
        else selecton=(resy[ii]>resy2[ii]);     // scores
        if  (selecton) {  //
                vg[ii]++;

                       }
             }
  return vg;
       }

int GraphGC::contAA(string Gcode, char aa){
    int c=0;
    for (unsigned int i = 0; i <  Gcode.length(); i++){
             if (Gcode[i]==aa) c++;
            }
    return c;
 }

vector <double>  GraphGC::Reassign1(vector < vector <double> >GC201models, ofstream &InputFile, string name, vector< float > Paa, vector< string > gcodes, bool select) {

    vector< string > gcodes1(gcodes);
    string a="ARN*DCQEGHILKMFPSTWYV";
    vector < vector <double> >  res;
    vector <double> freqs, freqs1;
    string ma;
    if (select) ma="Mayorq2";
    else ma="Mayorq1";
      InputFile<<name<<"      "<<ma<<"      ";
    vector <double>  vg{0,0,0,0,0,0,0,0,0,0,0,0};   //prop of best codes
    vector <double>  resy20a, resy21Ca,resy20a0, resy20Ca0;
    string gc, gc1;
    int ff=0, g=0;
    gc=gcodes[0];
    //GC201models=GlobalMeanVar201(Paa,gcodes);
    int L=12;
    for (unsigned int i = 0; i < L; i++){
             resy20a0.push_back(GC201models[4][i]);
                     }
    g=0;
    bool selector=true;
    while(g<gc.length()){
         int gg=0;
         int cc;
         while(gg<a.length()){
           if (!select){
                cc=contAA(gc, gc[g]);
            if  ((a[gg]!=gc[g])&&(cc>1))
            {
               gc1=gc;
               gc1[g]=a[gg];
               gcodes[0]=gc1;
              // GC201models=GlobalMeanVar201(Paa,gcodes);
               for (unsigned int i = 0; i <L ; i++){
                       resy20a.push_back(GC201models[4][i]); //scores
                        }
                  if (ff==0){
                    freqs=Bprob(InputFile, g, Paa, vg, gcodes1, resy20a0, resy20a,gc[g], a[gg], 1);    //scores
                   // freqs1=Bprob(InputFile, g, Paa, vg, gcodes1, resy20Ca0, resy20Ca, gc[g], a[gg], 1); //cantelli
                                   }
                  else  {
                           vg=freqs; //scores
                           freqs=Bprob(InputFile, g, Paa, vg, gcodes1, resy20a0, resy20a,gc[g], a[gg], 1);

                    }
         resy20a.clear();
         ff++;
               }
                }
        else if (select) {
             cc=contAA(gc, gc[g]);
           if  ((a[gg]!=gc[g])&&(cc>2))  {  //((a[gg]!=gc[g])&&(gc[g]!='W')&&(gc[g]!='M')&&(gc[g]!='K')&&(gc[g]!='Q')&&(gc[g]!='D')&&(gc[g]!='N')&&(gc[g]!='C')&&(gc[g]!='H')&&(gc[g]!='Y')&&(gc[g]!='E')&&(gc[g]!='F')) {
               gc1=gc;
               gc1[g]=a[gg];  //replace(gc.begin(), gc.end(), gc[g], a[gg]);
               gcodes[0]=gc1;
               //GC201models=GlobalMeanVar201(Paa,gcodes);
               for (unsigned int i = 0; i <L ; i++){
                       resy20a.push_back(GC201models[4][i]); //scores
                       //resy20Ca.push_back(GC201models[4][i]); //cantelli's bound
                        }
                  if (ff==0){
                    freqs=Bprob(InputFile, g, Paa, vg, gcodes1, resy20a0, resy20a,gc[g], a[gg], 1);    //mayor que scores
                   // freqs1=Bprob(InputFile,g,  Paa, vg, gcodes1, resy20Ca0, resy20Ca, gc[g], a[gg], 1); //menor que  cantelli
                                   }
                  else  {                                                                                //scores
                           vg=freqs;
                           freqs=Bprob(InputFile, g, Paa, vg, gcodes1, resy20a0, resy20a,gc[g], a[gg], 1);
                                                                                 //cantelli
                           //for (unsigned int u = 0; u <L ; u++) vg[u]=freqs1[u];
                           //freqs1=Bprob(InputFile, g, Paa, vg, gcodes1, resy20Ca0, resy20Ca,gc[g], a[gg], 1);                                               }
         resy20a.clear();
         //resy20Ca.clear();
         ff++;   }
                }
                 }
            gg++;
              }
     g++;
     }
  InputFile<<"      "<<ff;
  vector <double> res3, res4;
  for (unsigned int ii = 0; ii < L; ii++)
     {
         InputFile<<"      "<<freqs[ii]; //<<"      "<<freqs1[ii]; //score cantelli
         res3.push_back(freqs[ii]/ff);
       //  res4.push_back(freqs1[ii]/ff);
         }
    InputFile<<endl;
    return res3;
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


 //vert 64 gc 2 p=1 s=1 v=1

