#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <set>
#include "GraphGC.h"
using namespace std;
//computing the Genetic robustness
GraphGC::GraphGC(int n):vertices(n)
{
    geneM = new double*[vertices];
    phenoM = new double*[vertices];
    for (int r = 0; r< vertices; r++) {
        geneM[r] = new double[vertices];    //genetic graph
        phenoM[r] = new double[vertices];    //phenotypic distance graph
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
//b weight setting for genetic graph
 void GraphGC::setWG(int i, int j, double w, double g) {
                  geneM[i][j] = geneM[i][j]+(w*g);
    }
 //weight setting for genetic graph
 void GraphGC::setW(int i, int j, double w=0) {
                  geneM[i][j] = geneM[i][j]+w;
                 // geneM[j][i] = geneM[j][i]+w;
    }
 //phenotypic distance setting
 void GraphGC::setD(int i, int j, double d=0) {
                  phenoM[i][j] = d;
                  //phenoM[j][i] = d;
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
//conversion function 
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
double** GraphGC::getGeneM() {
    return geneM;
    }
double** GraphGC::getPhenoM() {
    return phenoM;
    }

double GraphGC::getGW(int i, int j) {
                  return geneM[i][j];
    }
double GraphGC::getPW(int i, int j) {
                  return phenoM[i][j];
    }
//inner prod projecting
double GraphGC::getProd(int i, int j) {
                  return ((phenoM[i][j]*geneM[i][j]));
    }
    //(Frob prod)/N, N: number of single base changes
double GraphGC::SumProd (int N){
  double s=0;
  for (int i = 0; i < vertices; i++)
   {
	//for (int ii = i+1; ii < vertices; ii++)
    for (int ii = 0; ii < vertices; ii++)
     {
      s=s+this->getProd(i, ii);
         }
         }
    return (s/N);
}
// separating invariant from variable part
int GraphGC::GCPart(vector<int> GCassign, vector<int> GCVassign, int iv=0){  //iv=0 GCassign: invariant set of nodes, iv=1 variable part
   int LC;
   for (int i = 0; i < vertices; i++)
   {
    for (int ii = 0; ii < vertices; ii++)
     {
        if (iv==0) {
          phenoM[GCassign[i]][ii]=0;
          geneM[GCassign[i]][ii]=0;
          //phenoM[ii][GCassign[i]]=0;
          //geneM[ii][GCassign[i]]=0;
          LC=GCassign.size();
                }
          else {
           phenoM[GCVassign[i]][ii]=0;
           geneM[GCassign[i]][ii]=0;
           //phenoM[ii][GCVassign[i]]=0;
           //geneM[ii][GCassign[i]]=0;
           LC=GCVassign.size();
                 }
         }
         }
      return LC;
    }

double GraphGC::SumGeneM(){
   double s=0;
   for (int i = 0; i < vertices; i++)
   {
    for (int ii = 0; ii < vertices; ii++)
          s=s+geneM[i][ii];
         }
      return s;
    }
double GraphGC::SumPhenoM(){
   double s=0;
   for (int i = 0; i < vertices; i++)
   {
    for (int ii = 0; ii < vertices; ii++)
          s=s+phenoM[i][ii];
         }
      return s;
    }
//homogeneous and heterogenous codon block partitioning
vector<vector<vector<double> > >GraphGC::hPart(int N, int Nse, vector< vector<double> > mPaa, vector< vector<double> > mGenm,  vector<int> GCHassign,  vector<int> GCVassign,  vector<int> GCSassign){
//P2
vector<vector<vector<double> > > outs;
vector<vector<double> > se(mPaa.size(), vector<double>(mGenm[0].size(), 0));
vector<vector<double> > co(mPaa.size(), vector<double>(mGenm[0].size(), 0));
 for (int k = 0; k < mPaa.size(); k++) {
        double sh=0;
        for (int i = 0; i < GCHassign.size(); i++)
             for (int ii = 0; ii < vertices; ii++)
                     if ((geneM[GCHassign[i]][ii])>0)
                        sh=sh+ (geneM[GCHassign[i]][ii]*((mPaa[k][GCHassign[i]]-mPaa[k][ii])*(mPaa[k][GCHassign[i]]-mPaa[k][ii])));
                     for (int g = 0; g < mGenm.size(); g++){
                          int sg=0, ss=0;
                          for (int j = 0; j < GCVassign.size(); j++)
                            for (int jj = 0; jj < vertices; jj++)
                                  if ((geneM[GCVassign[j]][jj])>0)
                                      sg=sg+((mGenm[g][GCVassign[j]])*(geneM[GCVassign[j]][jj])*((mPaa[k][GCVassign[j]]-mPaa[k][jj])*(mPaa[k][GCVassign[j]]-mPaa[k][jj])));
                           se[k][g]=(sh+sg);
                           for (int j = 0; j < GCSassign.size(); j++)
                            for (int jj = 0; jj < vertices; jj++)
                                 if ((geneM[GCSassign[j]][jj])>0)
                                      ss=ss+((mGenm[g][GCSassign[j]])*(geneM[GCSassign[j]][jj])*((mPaa[k][GCSassign[j]]-mPaa[k][jj])*(mPaa[k][GCSassign[j]]-mPaa[k][jj])));
                           co[k][g]=(se[k][g]+ss)/N;
                           se[k][g]=se[k][g]/Nse;
                       }
                 }
     outs.push_back(se);
     outs.push_back(co);
    return outs;
    }
//GC p-comparing
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
//n of sb changes
int GraphGC::PartChanges(vector<int> AAssign){
  int changes=0;
  for (int i = 0; i < vertices; i++)
   {
      for (int j = 0; j < vertices; j++)
       {
        if (geneM[AAssign[i]][j]>0) changes=changes+1;
          }
       }
   return changes;
    }
// gc pairwise comparison partitioning
vector<vector<vector<double> > >GraphGC::gcPart(int N, vector<int> Nch, vector< vector<double> > mPaa, vector< vector<int> > codes, vector<int>  GCSgcassign, vector< vector<int> >GCBassign,  vector< vector<int> >GCPassign,  vector< vector<int> > GCSassign)
{
  vector<double> se(mPaa.size(), 0);
  vector<double> sco(mPaa.size(), 0);
  vector<vector<double> > scp(mPaa.size(), vector<double>(GCBassign.size(), 0));
  vector<vector<double> > sce(mPaa.size(), vector<double>(GCBassign.size(), 0));
  for (int k = 0; k < mPaa.size(); k++)
   {
       double s=0, sh=0, sp=0, sb=0, ss=0;
       for (int i = 0; i < vertices; i++)
             for (int ii = 0; ii < vertices; ii++)
                     if ((geneM[i][ii])>0)
                      s=s+ (geneM[i][ii]*((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]))); //for the sgc
        for (int i = 0; i < GCSgcassign.size(); i++)
             for (int ii = 0; ii <vertices ; ii++)
                  if ((geneM[GCSgcassign[i]][ii])>0)
                    sh=sh+ (geneM[GCSgcassign[i]][ii]*((mPaa[k][GCSgcassign[i]]-mPaa[k][ii])*(mPaa[k][GCSgcassign[i]]-mPaa[k][ii]))); //for the sgc
        se[k]=s/Nch[0];
        for (int c = 0; c < codes.size(); c++)
         {
                double sp=0, sb=0, ss=0;
                for (int f = 0; f < GCPassign.size(); f++)
                   for (int j = 0; j < vertices; j++)
                         if ((geneM[GCPassign[c][f]][j])>0)
                         sp=sp+(geneM[GCPassign[c][f]][j]*((mPaa[k][GCPassign[c][f]]-mPaa[k][j])*(mPaa[k][GCPassign[c][f]]-mPaa[k][j])));
                  for (int f = 0; f < GCBassign.size(); f++)
                   for (int ii = 0; ii < vertices; ii++)
                         if ((geneM[GCBassign[c][f]][ii])>0)
                         sb=sb+ (geneM[GCBassign[c][f]][ii]*((mPaa[k][GCBassign[c][f]]-mPaa[k][ii])*(mPaa[k][GCBassign[c][f]]-mPaa[k][ii])));
                 for (int f = 0; f < GCSassign.size(); f++)
                   for (int ii = 0; ii < vertices; ii++)
                         if ((geneM[GCSassign[c][f]][ii])>0)
                          ss=ss+ (geneM[GCSassign[c][f]][ii]*((mPaa[k][GCSassign[c][f]]-mPaa[k][ii])*(mPaa[k][GCSassign[c][f]]-mPaa[k][ii])));
                  sce[k][c]=((s+sh)-sp+sb-ss)/Nch[c];
                  scp[k][c]=((s+sh)-sp+sb)/N;
                                }
          sco[k]=(s+sh)/N;
     }
   for (int f = 0; f < mPaa.size(); f++) {
         scp[f][GCBassign.size()]=se[f];
         sce[f][GCBassign.size()]=sco[f]; }

   vector<vector<vector<double> > > outs;

   outs.push_back(scp);
   outs.push_back(sce);

   return outs;
}
// null mean hh partitioning

vector<vector<vector<double> > > GraphGC::mePart(int N, vector<int> Nbse, double Tco, vector< vector<double> > mPaa,  vector< vector<int> > nc, vector< vector<int> > nsc,  vector<double> Tbc)
{
  double sed;
  vector<double> sb(mPaa.size(), 0);
  vector<double>  Tse(nc.size(),0);
  vector<vector<double> > se(mPaa.size(), vector<double>(nc.size(), 0));
  vector<vector<double> > sco(mPaa.size(), vector<double>(nc.size(), 0));
   vector<vector<double> > mb(mPaa.size(), vector<double>(nc.size(), 0));
  vector<vector<double> > mse(mPaa.size(), vector<double>(nc.size(), 0));
  vector<vector<double> > mco(mPaa.size(), vector<double>(nc.size(), 0));
  //nsc vector per code of stoppositions, nc number of aas per cod

  for (int k = 0; k < mPaa.size(); k++)
   {
         for (int i = 0; i< 20; i++)
             for (int ii = 0; ii< 20; ii++) sb[k]=sb[k]+((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));

         for (int c = 0; c < nc.size(); c++)
         {
             for (int i = 0; i< 20; i++)
             for (int ii = 0; ii< 20; ii++) se[k][c]=se[k][c]+(nc[c][i]*nc[c][ii]*((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii])));

             for (int i = 0; i< nsc[c].size(); i++)
             for (int ii = 0; ii< 20; ii++)
                      sco[k][c]=se[k][c]+(nc[c][ii]*(mPaa[k][nsc[c][i]]-mPaa[k][ii])*(mPaa[k][nsc[c][i]]-mPaa[k][ii]));
              }
     }
  for (int c = 0; c < nc.size(); c++)
   {
       double sed=0;
       for (int i = 0; i< nsc[c].size(); i++)
             for (int ii = 0; ii< vertices; ii++)
                 sed=sed+geneM[nsc[c][i]][ii];
                 Tse[c]=  Tco-sed;
  for (int k = 0; k < mPaa.size(); k++)
   {
       mb[k][c]=(Tbc[c]*sb[k])/(Nbse[c]*380);
       mse[k][c]=((Tse[c])*(se[k][c]))/(Nbse[c]*(64-nsc[c].size())*(63-nsc[c].size()));
       mco[k][c]=(Tco*(sco[k][c]+se[k][c]))/(N*4032);
        }
            }
      vector<vector<vector<double> > > outs;
       outs.push_back(mb);
       outs.push_back(mse);
       outs.push_back(mco);
       return outs;
 }

//null mean GC pairwise comparison partitioining

 vector<vector<vector<double> > > GraphGC::meGPart(int N, int Nbse,vector< vector<double> > ng, vector< vector<double> > mPaa,  vector<int> nc, vector<int>nsc,  vector<double> Tbc)
{
  double Tse;
  vector<double> sb(20, 0);
  vector<vector<double> > se(20, 0));
  vector<vector<double> > sco(nsc.size(), 0);
   vector<vector<double> > mb(mPaa.size(), vector<double>(ng.size(), 0));
  vector<vector<double> > mse(mPaa.size(), vector<double>(ng.size(), 0));
  vector<vector<double> > mco(mPaa.size(), vector<double>(ng.size(), 0));
  //nsc vector per code of stoppositions, nc number of aas per cod

  for (int k = 0; k < mPaa.size(); k++)
   {
     for (int i = 0; i< 20; i++)
     for (int ii = 0; ii< 20; ii++) sb[k]=sb[k]+((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));

     for (int i = 0; i< 20; i++)
     for (int ii = 0; ii< 20; ii++) se[k]=se[k]+(nc[i]*nc[ii]*((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii])));

     for (int i = 0; i< nsc.size(); i++)
     for (int ii = 0; ii< 20; ii++) sco[k]=se[k]+(nc[ii]*(mPaa[k][nsc[i]]-mPaa[k][ii])*(mPaa[k][nsc[i]]-mPaa[k][ii]));

     }
  for (int g = 0; g< ng.size(); g++)
   {
       double sed=0;
       for (int i = 0; i< nsc.size(); i++)
             for (int ii = 0; ii< vertices; ii++)
                 sed=sed+(ng[g][nsc[i]]*geneM[nsc[i]][ii]);
       double ed=0;
       for (int i = 0; i< vertices; i++)
             for (int ii = 0; ii< vertices; ii++)
                 ed=ed+(ng[g][i]*geneM[i][ii]);
                 Tse= ed-sed;
  for (int k = 0; k < mPaa.size(); k++)
   {
       mb[k][g]=(Tbc[g]*sb[k])/(Nbse*380);
       mse[k][g]=((Tse)*(se[k]))/(Nbse*(64-nsc.size())*(63-nsc.size()));
       mco[k][g]=(ed*(sco[k]+se[k]))/(N*4032);
        }
            }
      vector<vector<vector<double> > > outs;
       outs.push_back(mb);
       outs.push_back(mse);
       outs.push_back(mco);
       return outs;
 }
 //null var GC/Genome hh partitioining
 vector<vector<vector<double> > > GraphGC::vaPart(vector<vector<vector<double> > > outs, vector<int> Nbse, int N, double Tco, vector< vector<double> > mPaa, vector< vector<int> > nc, vector< vector<int> > nsc,  vector< vector<double> > Tb)
 {
     double T4,s=0, ss=0, sw=0, ssw0=0,ssw=0, ssw2=0;
     vector<double>  pb(mPaa.size(),0);
     vector<double>  pb3(mPaa.size(),0);
     vector<double>  pb4(mPaa.size(),0);
     vector<vector<double> > pse(mPaa.size(), vector<double>(nc.size(), 0));
     vector<vector<double> > pse3(mPaa.size(), vector<double>(nc.size(), 0));
     vector<vector<double> > pse4(mPaa.size(), vector<double>(nc.size(), 0));
     vector<vector<double> > pco(mPaa.size(), vector<double>(nc.size(), 0));
     vector<vector<double> > pco3(mPaa.size(), vector<double>(nc.size(), 0));
     vector<vector<double> > pco4(mPaa.size(), vector<double>(nc.size(), 0));
     //vector<double>  sw(nc.size(),0);
    for (int i = 0; i< vertices; i++)
             for (int ii = 0; ii< vertices; ii++) s=s+geneM[i][ii];
    for (int i = 0; i< vertices; i++)
             for (int ii = 0; ii< vertices; ii++) ss=ss+(geneM[i][ii]*geneM[i][ii]);
    for (int i = 0; i< vertices; i++) {
             sw=0;
             for (int ii = 0; ii< vertices; ii++) sw=sw+geneM[i][ii];
                 ssw0=ssw0+(sw*sw);
                    }
    T4=(Tco*Tco)-(4*(ssw0))+(2*ss);
  for (int k = 0; k < mPaa.size(); k++)
   {
       ssw=0;ssw2=0;
        for (int i = 0; i< 20; i++)
             for (int ii = 0; ii< 20; ii++)  pb[k]=pb[k]+((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));

       for (int i = 0; i< 20; i++) {
             sw=0;
             for (int ii = 0; ii< 20; ii++) sw=sw+((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));
                 ssw=ssw+(sw*sw);
                 ssw2=ssw2+sw;
                    }
       pb3[k]=ssw-pb[k];
       pb4[k]=(ssw2*ssw2)-(4*pb3[k])+(2*pb[k]);
       for (int c = 0; c < nc.size(); c++)
       {
        double ssw3=0, sstp=0, sstp2=0, ssu=0, sstpa=0, sstpb=0;
           ssw=0;ssw2=0;
            for (int i = 0; i< 20; i++)
            for (int ii = 0; ii< 20; ii++)  pse[k][c]=pse[k][c]+(nc[c][i]*nc[c][ii])*((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));
            for (int i = 0; i< 20; i++)
            for (int ii = 0; ii< 20; ii++)  ssu=ssu+(nc[c][i]*nc[c][ii])*((mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));
            for (int i = 0; i< 20; i++) {
             sw=0;
            for (int ii = 0; ii< 20; ii++) {sw=sw+(nc[c][ii]*(mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));
                 ssw3=ssw3+(nc[c][i]*nc[c][ii]*(mPaa[k][i]-mPaa[k][ii])*(mPaa[k][i]-mPaa[k][ii]));
                                  }
                 ssw=ssw+(nc[c][i]*sw*sw);
                    }
             pse3[k][c]=ssw-pse[k][c];
             pse4[k][c]=(ssw3*ssw3)-(4*pse3[k][c])+(2*pse[k][c]);
             for (int i = 0; i< nsc[c].size(); i++)
                for (int ii = 0; ii< 20; ii++)
                      sstp=(nc[c][ii]*(mPaa[k][nsc[c][i]]-mPaa[k][ii])*(mPaa[k][nsc[c][i]]-mPaa[k][ii])*(mPaa[k][nsc[c][i]]-mPaa[k][ii])*(mPaa[k][nsc[c][i]]-mPaa[k][ii]));
              for (int i = 0; i< nsc[c].size(); i++)
                for (int ii = 0; ii< 20; ii++)
                      sstpb=sstpb+(nc[c][ii]*(mPaa[k][nsc[c][i]]-mPaa[k][ii])*(mPaa[k][nsc[c][i]]-mPaa[k][ii]));

              for (int i = 0; i< nsc[c].size(); i++) {
                sstpa=0;
                for (int ii = 0; ii< 20; ii++) sstpa=sstpa+ (nc[c][ii]*(mPaa[k][nsc[c][i]]-mPaa[k][ii])*(mPaa[k][nsc[c][i]]-mPaa[k][ii]));
                 sstp2=sstp2+(sstpa*sstpa);
                 }
               pco[k][c]=pse[k][c]+sstp;
               sstp2=sstp2-sstp;
               pco3[k][c]=pse3[k][c]+sstp2;
               pco4[k][c]=((sstpb+ssu)*(sstpb+ssu))-(4*pco3[k][c])+(2*pco[k][c]);
                    }
       }
      vector<double>  Tse(nc.size(),0);
      vector<double>  Tse3(nc.size(),0);
      vector<double>  Tse4(nc.size(),0);
      vector<double>  tb0(nc.size(),0);
      vector<double>  tb3(nc.size(),0);
      vector<double>  tb4(nc.size(),0);
      vector <vector <vector < double > > > ovar (3, vector < vector < double > > (mPaa.size(),vector<double> (nc.size(),0)));
   for (int c = 0; c < nc.size(); c++)
       {
            double tsd=0,stp=0,stp2=0, tss=0, tss1=0, tss2=0;
            for (int i = 0; i< nsc[c].size(); i++)
                for (int ii = 0; ii< vertices; ii++) stp=stp+geneM[nsc[c][i]][ii];

            for (int i = 0; i< nsc[c].size(); i++)
                for (int ii = 0; ii< vertices; ii++) tsd=tsd+(geneM[nsc[c][i]][ii]*geneM[nsc[c][i]][ii]);

               Tse[c]= Tco-tsd;
            for (int i = 0; i< nsc[c].size(); i++) {
                stp2=0;
                for (int ii = 0; ii< vertices; ii++) stp=stp+ (geneM[nsc[c][i]][ii]);
                 stp2=stp2+(stp*stp);
                 }
           Tse3[c]=ssw0-stp2;
           Tse4[c]=((s-stp)*(s-stp))-(4*Tse3[c])+(2*Tse[c]);
           for (int i = 0; i< Tb.size(); i++)
                for (int ii = 0; ii< Tb[i].size(); ii++) tb0[c]=tb0[c]+(Tb[i][ii]*Tb[i][ii]);
            for (int i = 0; i< Tb.size(); i++) {
                tss=0;
                for (int ii = 0; ii< Tb[i].size(); ii++) tss=tss+ (Tb[i][ii]);
                 tss2=tss2+(tss*tss);
                 tss1=tss1+(tss);
                 }
          tb3[c]=tss2-tb0[c];
          tb4[c]=(tss1*tss1)-(4*tb3[c])+(2*tb0[c]);
      for (int k = 0; k < mPaa.size(); k++)
       {
          ovar[0][k][c]=((((tb4[c]*pb4[k])/306) + ((4*tb3[c]*pb3[k])/18) + (2*tb0[c]*pb[k]))/(Nbse[c]*Nbse[c]*380))-(outs[0][k][c]*outs[0][k][c]);
          ovar[1][k][c]=((((Tse4[c]*pse4[k][c])/((61-nsc[c].size())*(62-nsc[c].size()))) + ((4*Tse3[c]*pse3[k][c])/(62-nsc[c].size()))+(2*Tse[c]*pse[k][c]))/((63-nsc[c].size())*(64-nsc[c].size())*Nbse[c]*Nbse[c]))-(outs[1][k][c]*outs[1][k][c]);
          ovar[2][k][c]=((((T4*pco4[k][c])/3782)+((4*ssw0*pco3[k][c])/62)+(2*ss*pco[k][c]))/(N*N*4032))-(outs[2][k][c]*outs[2][k][c]);
         }
               }
    return ovar;
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
                           if ((geneM[i][j]==geneM[ii][t])&&(AAGC[j]==AAGC[t])){
                                           a=1;
                                           }
                                       }
                                      }

                                if (a==0) {
                                    break;
                                   }
                                 }
                              }
                              //cout<<j<<endl;
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

GraphGC::GraphGC(vector <double>  &Paa, vector <string> &GCodes, int suppr, int ub, int vert, int gc, int p, int s, int v){
  /*  p=0 s=0 v=1  transv model
    p=0 s=1 v=0    trans model
    p=3 s=1 v=1    3 pos model
    p=2 s=1 v=1    2 pos model
    p=1 s=1 v=1    1 pos model         vert 64 gc 2 p=1 s=1 v=1
    p=0 s=1 v=1     tot model
    gc=0  block model
    gc=1  sense model
    gc=2  codon model
    ub=0  biased ub 1 unbiased
    */
  vertices=vert;//for sgc 64 61 20
  double d, da1, da2;
  vector <double> Paa64;
  pair<int,char> po, po2;
  string alp=GCodes[0], gcode=GCodes[0];
  Paa=insert20(Paa);
  if (suppr==0) Paa64=P20ToP64(alp, Paa);
  else {
      Paa64=P20ToP64(alp, Paa);
      Paa64=allMeanSuppresors(Paa64,GCodes);
                 }
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
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setW(i, ii, weights[2][0]);
                    else       setW(i, ii, 1);}
                if (gc==2){
                    if (ub==0) setW(i, ii, weights[2][0]);
                    else       setW(i, ii, 1);}
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                    if (ub==0) setW(i, ii, weights[2][0]);
                    else       setW(i, ii, 1);}
                }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])) {
                    if (ub==0) setW(i, ii, weights[2][1]);
                    else       setW(i, ii, 1);  }
                if (gc==2) {
                    if (ub==0) setW(i, ii, weights[2][1]);
                    else       setW(i, ii, 1); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                    if (ub==0) setW(i, ii, weights[2][1]);
                    else       setW(i, ii, 1); }
                }
          }
       if (((p==2)||(p==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {

           if ((s==1)&&(transition(GCodes[3][i], GCodes[3][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setW(i, ii, weights[1][0]);
                    else       setW(i, ii, 1); }
                if (gc==2){
                   if (ub==0) setW(i, ii, weights[1][0]);
                     else     setW(i, ii, 1); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                     if (ub==0) setW(i, ii, weights[1][0]);
                     else       setW(i, ii, 1); }
            }
          if ((v==1)&&(transversion(GCodes[3][i], GCodes[3][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])) {
                    if (ub==0) setW(i, ii, weights[1][1]);
                    else       setW(i, ii, 1); }
                if (gc==2) {
                   if (ub==0) setW(i, ii, weights[1][1]);
                   else       setW(i, ii, 1); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                   if (ub==0) setW(i, ii, weights[1][1]);
                   else       setW(i, ii, 1); }
                             }
               }
         if (((p==1)||(p==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {

           if ((s==1)&&(transition(GCodes[2][i], GCodes[2][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setW(i, ii, weights[0][0]);
                    else       setW(i, ii, 1); }
                if (gc==2){
                    if (ub==0) setW(i, ii, weights[0][0]);
                    else       setW(i, ii, 1); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))){
                    if (ub==0) setW(i, ii, weights[0][0]);
                    else       setW(i, ii, 1); }
            }
          if ((v==1)&&(transversion(GCodes[2][i], GCodes[2][ii]))){

                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setW(i, ii, weights[0][1]);
                    else       setW(i, ii, 1); }
                if (gc==2){
                    if (ub==0) setW(i, ii, weights[0][1]);
                    else       setW(i, ii, 1); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))){
                    if (ub==0) setW(i, ii, weights[0][1]);
                    else       setW(i, ii, 1); }
                                  }
                     }
    }

}
}
GraphGC::GraphGC(vector <double>  &codfreq, vector <double>  &Paa, vector <string> &GCodes, int suppr, int ub, int vert, int gc, int p, int s, int v){
  /*  p=0 s=0 v=1  transv model
    p=0 s=1 v=0    trans model
    p=3 s=1 v=1    3 pos model
    p=2 s=1 v=1    2 pos model
    p=1 s=1 v=1    1 pos model         vert 64 gc 2 p=1 s=1 v=1
    p=0 s=1 v=1     tot model
    gc=0  block model
    gc=1  sense model
    gc=2  codon model
    ub=0  biased ub 1 unbiased
    */
  vertices=vert;//for sgc 64 61 20
  double d, da1, da2;
  vector <double> Paa64;
  pair<int,char> po, po2;
  string alp=GCodes[0], gcode=GCodes[0];
  Paa=insert20(Paa);
  if (suppr==0) Paa64=P20ToP64(alp, Paa);
  else {
      Paa64=P20ToP64(alp, Paa);
      Paa64=allMeanSuppresors(Paa64,GCodes);
                 }
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
    for (int ii = 0; ii < GCodes[0].size(); ii++)
     {
        //3rst pos
       if ((p==3)||(p==0)&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){

           if ((s==1)&&(transition(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setWG(i, ii, weights[2][0], codfreq[i]);
                    else       setW(i, ii, codfreq[i]);}
                if (gc==2){
                    if (ub==0) setW(i, ii, weights[2][0], codfreq[i]);
                    else       setW(i, ii, codfreq[i]);}
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                    if (ub==0) setW(i, ii, weights[2][0], codfreq[i]);
                    else       setW(i, ii, codfreq[i]);}
                }
          if ((v==1)&&(transversion(GCodes[4][i], GCodes[4][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])) {
                    if (ub==0) setW(i, ii, weights[2][1], codfreq[i]);
                    else       setW(i, ii, codfreq[i]);  }
                if (gc==2) {
                    if (ub==0) setW(i, ii, weights[2][1], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                    if (ub==0) setW(i, ii, weights[2][1], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                }
          }
       if (((p==2)||(p==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {

           if ((s==1)&&(transition(GCodes[3][i], GCodes[3][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setW(i, ii, weights[1][0], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                if (gc==2){
                   if (ub==0) setW(i, ii, weights[1][0], codfreq[i]);
                     else     setW(i, ii, codfreq[i]); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                     if (ub==0) setW(i, ii, weights[1][0], codfreq[i]);
                     else       setW(i, ii, codfreq[i]); }
            }
          if ((v==1)&&(transversion(GCodes[3][i], GCodes[3][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])) {
                    if (ub==0) setW(i, ii, weights[1][1], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                if (gc==2) {
                   if (ub==0) setW(i, ii, weights[1][1], codfreq[i]);
                   else       setW(i, ii, codfreq[i]); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))) {
                   if (ub==0) setW(i, ii, weights[1][1], codfreq[i]);
                   else       setW(i, ii, codfreq[i]); }
                             }
               }
         if (((p==1)||(p==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {

           if ((s==1)&&(transition(GCodes[2][i], GCodes[2][ii]))){
                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setW(i, ii, weights[0][0], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                if (gc==2){
                    if (ub==0) setW(i, ii, weights[0][0], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))){
                    if (ub==0) setW(i, ii, weights[0][0], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
            }
          if ((v==1)&&(transversion(GCodes[2][i], GCodes[2][ii]))){

                if ((gc==0)&&(GCodes[0][i]!=GCodes[0][ii])){
                    if (ub==0) setW(i, ii, weights[0][1], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                if (gc==2){
                    if (ub==0) setW(i, ii, weights[0][1], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
                if ((gc==1)&&(!STOP(GCodes[0][i], GCodes[0][ii]))){
                    if (ub==0) setW(i, ii, weights[0][1], codfreq[i]);
                    else       setW(i, ii, codfreq[i]); }
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
//vector de conte
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

//permGenCodes
vector< double > GraphGC::permGenRobustness(string alp, const vector <double>  &CodBias, const vector <double>  &Faas, vector <double>  &Paa, const vector <string> &GCodes) {

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

                       sumt3=sumt3+((CodBias[i])*d*mt);  //
                       sp3=sp3+((weights[2][0])*(CodBias[i])*d*mt);
                       spaa=spaa+(weights[2][0])*Faas[i]*d*mt;
                        transf=transf+((CodBias[i])*d*mt); //******
                       transd=transd+((weights[2][0])*(CodBias[i])*d*mt); //******

                       sumt3stop=sumt3stop+((CodBias[i])*dstop*mt);  //
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
//permGenCodes64 
vector< double > GraphGC::permGenRobustness64(string alp, const vector <double>  &CodBias, const vector <double>  &Faas,  const vector <string> &GCodes) {
 //Codbias, faas, fcod para cg64
  int codonpos=0;
  pair<int,char> po, po2;

    int mt=0, nts1=0, nts2=0, nts3=0;
  int ntts1=0, ntts2=0, ntts3=0;

 double sumt1=0, sumt2=0, sumt3=0, sumw1=0, sumw2=0, sumw3=0;
  double sp1=0, sp2=0, sp3=0, spd1=0, spd2=0, spd3=0;
   double sumt1F=0, sumt2F=0, sumt3F=0, sumw1F=0, sumw2F=0, sumw3F=0;
  double sp1F=0, sp2F=0, sp3F=0, spd1F=0, spd2F=0, spd3F=0;

  double ttrans=0, ttranv=0;
 double  transd64=0, ttransd64=0, transdw64=0, ttransdw64=0,tranvd64=0, ttranvd64=0, tranvdw64=0, ttranvdw64=0;
 double  transd64F=0, ttransd64F=0, transdw64F=0, ttransdw64F=0     ,tranvd64F=0, ttranvd64F=0,tranvdw64F=0, ttranvdw64F=0;

  double sT3w1=0,sT3w2=0,sT3w3=0, sT3w=0, sT31=0, sT32=0,sT33=0, sT3=0, sT3wtrans=0,sT3trans=0, sT3wtranv=0, sT3tranv=0;
  double sT4w3=0,sT4w2=0,sT4w1=0, sT4w=0, sT41=0, sT42=0,sT43=0, sT4=0, sT4wtrans=0,sT4trans=0, sT4wtranv=0, sT4tranv=0;
  double d=1, spaa=0, spaa2=0, spaaF=0, spaa2F=0, sT4aa=0, sT3aa=0;

int ntt1=0, ntt2=0, ntt3=0;

 for (unsigned int i = 0; i < GCodes[0].size(); i++)
   {
       spd1F=0, spd2F=0, spd3F=0;sp1F=0, sp2F=0, sp3F=0;
    sumt1F=0, sumt2F=0,sumt3F=0, sumw3F=0, sumw2F=0, sumw1F=0;
    transd64F=0, ttransd64F=0,transd64F=0, ttransd64F=0,transdw64F=0, ttransdw64F=0, tranvdw64F=0, ttranvdw64F=0;
    spaaF=0, spaa2F=0;

    for (unsigned int ii = 0; ii < GCodes[0].size(); ii++)
     {
    //3rst pos
    if (ii!=i)  {
       if (((codonpos==3)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){
                  if (transition(GCodes[4][i], GCodes[4][ii])){

              //  if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                   //    else {d=1; ttrans=ttrans+1; ntt3=ntt3+1; }

                      if (CodBias[i]==1) mt=0; else
                        mt=1; d=1;
                        //global mean
                       sumt3=sumt3+((CodBias[i])*mt*d);  //cb cg20
                       sumw3=sumw3+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp3=sp3+((weights[2][0])*(CodBias[i])*mt*d); //w+cb cg20
                       spd3=spd3+(((weights[2][0])*(CodBias[i]))*((weights[2][0])*(CodBias[i])*mt)*d);// (w+cb)(w+cb) cg20

                       spaa=spaa+((weights[2][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2=spaa2+((weights[2][0])*(Faas[i])*(weights[2][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaaF=spaaF+((weights[2][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2F=spaa2F+((weights[2][0])*(Faas[i])*(weights[2][0])*(Faas[i])*mt*d); //w+cb cg20

                       transd64=transd64+((CodBias[i])*mt*d);  //****************
                       ttransd64=ttransd64+((CodBias[i])*(CodBias[i])*mt*d);
                       transdw64=transdw64+((weights[2][0])*(CodBias[i])*mt*d);  //****************
                       ttransdw64=ttransdw64+((weights[2][0])*(weights[2][0])*(CodBias[i])*(CodBias[i])*mt*d);

                       sumt3F=sumt3F+((CodBias[i])*mt*d);  //cb cg20
                       sumw3F=sumw3F+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp3F=sp3F+((weights[2][0])*(CodBias[i])*mt*d); //w+cb cg20
                       spd3F=spd3F+(((weights[2][0])*(CodBias[i]))*((weights[2][0])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       transd64F=transd64F+((CodBias[i])*mt*d);  //****************
                       ttransd64F=ttransd64F+((CodBias[i])*(CodBias[i])*mt*d);
                       transdw64F=transdw64F+((weights[2][0])*(CodBias[i])*mt*d);  //****************
                       ttransdw64F=ttransdw64F+((weights[2][0])*(weights[2][0])*(CodBias[i])*(CodBias[i])*mt*d);
                      // ttrans=ttrans+1;
                                           }
                  if (transversion(GCodes[4][i], GCodes[4][ii])){

                    //if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                     // else {d=1; ttranv=ttranv+1; ntt3=ntt3+1; }

                     if (CodBias[i]==1) mt=0; else
                        mt=1; d=1;
                        //global mean
                      sumt3=sumt3+((CodBias[i])*mt*d);  //cb cg20
                       sumw3=sumw3+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp3=sp3+((weights[2][1])*(CodBias[i])*mt*d); //w+cb cg20
                       spd3=spd3+(((weights[2][1])*(CodBias[i]))*((weights[2][1])*(CodBias[i])*mt)*d);// (w+cb)(w+cb) cg20

                        spaa=spaa+((weights[2][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2=spaa2+((weights[2][1])*(Faas[i])*(weights[2][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaaF=spaaF+((weights[2][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2F=spaa2F+((weights[2][1])*(Faas[i])*(weights[2][1])*(Faas[i])*mt*d); //w+cb cg20

                       tranvd64=tranvd64+((CodBias[i])*mt*d);  //****************
                       ttranvd64=ttranvd64+((CodBias[i])*(CodBias[i])*mt*d);
                       tranvdw64=tranvdw64+((weights[2][1])*(CodBias[i])*mt*d);  //****************
                       ttranvdw64=ttranvdw64+((weights[2][1])*(weights[2][1])*(CodBias[i])*(CodBias[i])*mt*d);

                       sumt3F=sumt3F+((CodBias[i])*mt*d);  //cb cg20
                       sumw3F=sumw3F+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp3F=sp3F+((weights[2][1])*(CodBias[i])*mt*d); //w+cb cg20
                       spd3F=spd3F+(((weights[2][1])*(CodBias[i]))*((weights[2][1])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       tranvd64F=tranvd64F+((CodBias[i])*mt*d);  //****************
                       ttranvd64F=ttranvd64F+((CodBias[i])*(CodBias[i])*mt*d);
                       tranvdw64F=tranvdw64F+((weights[2][1])*(CodBias[i])*mt*d);  //****************
                       ttranvdw64F=ttranvdw64F+((weights[2][1])*(weights[2][1])*(CodBias[i])*(CodBias[i])*mt*d);
                    //   ttranv=ttranv+1;
                               }
                                                        //nor *  //nor * and nor *
              }
       if (((codonpos==2)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {
                if (transition(GCodes[3][i], GCodes[3][ii])){
                     // if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                       // else {d=1; ttrans=ttrans+1; ntt2=ntt2+1; }

                        if (CodBias[i]==1) mt=0; else
                        mt=1; d=1;
                    //global mean
                       sumt2=sumt2+((CodBias[i])*mt*d);  //cb cg20
                       sumw2=sumw2+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp2=sp2+((weights[1][0])*(CodBias[i])*mt*d); //w+cb cg20
                       spd2=spd2+(((weights[1][0])*(CodBias[i]))*((weights[1][0])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       spaa=spaa+((weights[1][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2=spaa2+((weights[1][0])*(Faas[i])*(weights[1][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaaF=spaaF+((weights[1][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2F=spaa2F+((weights[1][0])*(Faas[i])*(weights[1][0])*(Faas[i])*mt*d); //w+cb cg20

                       transd64=transd64+((CodBias[i])*mt*d);  //****************
                       ttransd64=ttransd64+((CodBias[i])*(CodBias[i])*mt*d);
                       transdw64=transdw64+((weights[1][0])*(CodBias[i])*mt*d);  //****************
                       ttransdw64=ttransdw64+((weights[1][0])*(weights[1][0])*(CodBias[i])*(CodBias[i])*mt*d);

                       sumt2F=sumt2F+((CodBias[i])*mt*d);  //cb cg20
                       sumw2F=sumw2F+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp2F=sp2F+((weights[1][0])*(CodBias[i])*mt*d); //w+cb cg20
                       spd2F=spd2F+(((weights[1][0])*(CodBias[i]))*((weights[1][0])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       transd64F=transd64F+((CodBias[i])*mt*d);  //****************
                       ttransd64F=ttransd64F+((CodBias[i])*(CodBias[i])*mt*d);
                       transdw64F=transdw64F+((weights[1][0])*(CodBias[i])*mt*d);  //****************
                       ttransdw64F=ttransdw64F+((weights[1][0])*(weights[1][0])*(CodBias[i])*(CodBias[i])*mt*d);
                      // ttrans=ttrans+1;
                               }
                 if (transversion(GCodes[3][i], GCodes[3][ii])){
                      if (CodBias[i]==1) mt=0; else
                        mt=1; d=1;
                    //if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                      // else {d=1; ttranv=ttranv+1; ntt2=ntt2+1; }

                        //global mean
                        sumt2=sumt2+((CodBias[i])*mt*d);  //cb cg20
                       sumw2=sumw2+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp2=sp2+((weights[1][1])*(CodBias[i])*mt*d); //w+cb cg20
                       spd2=spd2+(((weights[1][1])*(CodBias[i]))*((weights[1][1])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                        spaa=spaa+((weights[1][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2=spaa2+((weights[1][1])*(Faas[i])*(weights[1][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaaF=spaaF+((weights[1][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2F=spaa2F+((weights[1][1])*(Faas[i])*(weights[1][1])*(Faas[i])*mt*d); //w+cb cg20

                       tranvd64=tranvd64+((CodBias[i])*mt*d);  //****************
                       ttranvd64=ttranvd64+((CodBias[i])*(CodBias[i])*mt*d);
                       tranvdw64=tranvdw64+((weights[1][1])*(CodBias[i])*mt*d);  //****************
                       ttranvdw64=ttranvdw64+((weights[1][1])*(weights[1][1])*(CodBias[i])*(CodBias[i])*mt*d);

                       sumt2F=sumt2F+((CodBias[i])*mt*d);  //cb cg20
                       sumw2F=sumw2F+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp2F=sp2F+((weights[1][1])*(CodBias[i])*mt*d); //w+cb cg20
                       spd2F=spd2F+(((weights[1][1])*(CodBias[i]))*((weights[1][1])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       tranvd64F=tranvd64F+((CodBias[i])*mt*d);  //****************
                       ttranvd64F=ttranvd64F+((CodBias[i])*(CodBias[i])*mt*d);
                       tranvdw64F=tranvdw64F+((weights[1][1])*(CodBias[i])*mt*d);  //****************
                       ttranvdw64F=ttranvdw64F+((weights[1][1])*(weights[1][1])*(CodBias[i])*(CodBias[i])*mt*d);
                    //   ttranv=ttranv+1;
                               }
                                 }
       if (((codonpos==1)||(codonpos==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {
                // cout<<po2.first<<" ";
                 //if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                   //     else {d=1; ttrans=ttrans+1; ntt1=ntt1+1; }

                  if (transition(GCodes[2][i], GCodes[2][ii])){
                      if (CodBias[i]==1) mt=0; else
                        mt=1; d=1;
                        //global mean
                       sumt1=sumt1+((CodBias[i])*mt*d);  //cb cg20
                       sumw1=sumw1+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp1=sp1+((weights[0][0])*(CodBias[i])*mt*d); //w+cb cg20
                       spd1=spd1+(((weights[0][0])*(CodBias[i]))*((weights[0][0])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                        spaa=spaa+((weights[0][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2=spaa2+((weights[0][0])*(Faas[i])*(weights[0][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaaF=spaaF+((weights[0][0])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2F=spaa2F+((weights[0][0])*(Faas[i])*(weights[0][0])*(Faas[i])*mt*d); //w+cb cg20

                       transd64=transd64+((CodBias[i])*mt*d);  //****************
                       ttransd64=ttransd64+((CodBias[i])*(CodBias[i])*mt*d);
                       transdw64=transdw64+((weights[0][0])*(CodBias[i])*mt*d);  //****************
                       ttransdw64=ttransdw64+((weights[0][0])*(weights[0][0])*(CodBias[i])*(CodBias[i])*mt*d);

                       sumt1F=sumt1F+((CodBias[i])*mt*d);  //cb cg20
                       sumw1F=sumw1F+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp1F=sp1F+((weights[0][0])*(CodBias[i])*mt*d); //w+cb cg20
                       spd1F=spd1F+(((weights[0][0])*(CodBias[i]))*((weights[0][0])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       transd64F=transd64F+((CodBias[i])*mt*d);  //****************
                       ttransd64F=ttransd64F+((CodBias[i])*(CodBias[i])*mt*d);
                       transdw64F=transdw64F+((weights[0][0])*(CodBias[i])*mt*d);  //****************
                       ttransdw64F=ttransdw64F+((weights[0][0])*(weights[0][0])*(CodBias[i])*(CodBias[i])*mt*d);
                    }
                 if (transversion(GCodes[2][i], GCodes[2][ii])){
                     // if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                     //  else {d=1; ttranv=ttranv+1; ntt1=ntt1+1; }
                       if (CodBias[i]==1) mt=0; else
                        mt=1; d=1;
                     //global mean
                     sumt1=sumt1+((CodBias[i])*mt*d);  //cb cg20
                       sumw1=sumw1+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp1=sp1+((weights[0][1])*(CodBias[i])*mt*d); //w+cb cg20
                       spd1=spd1+(((weights[0][1])*(CodBias[i]))*((weights[0][1])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       spaa=spaa+((weights[0][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2=spaa2+((weights[0][1])*(Faas[i])*(weights[0][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaaF=spaaF+((weights[0][1])*(Faas[i])*mt*d); //w+cb cg20
                       spaa2F=spaa2F+((weights[0][1])*(Faas[i])*(weights[0][1])*(Faas[i])*mt*d); //w+cb cg20

                       tranvd64=tranvd64+((CodBias[i])*mt*d);  //****************
                       ttranvd64=ttranvd64+((CodBias[i])*(CodBias[i])*mt*d);
                       tranvdw64=tranvdw64+((weights[0][1])*(CodBias[i])*mt*d);  //****************
                       ttranvdw64=ttranvdw64+((weights[0][1])*(weights[0][1])*(CodBias[i])*(CodBias[i])*mt*d);

                       sumt1F=sumt1F+((CodBias[i])*mt*d);  //cb cg20
                       sumw1F=sumw1F+((CodBias[i])*(CodBias[i])*mt*d); //w cg20
                       sp1F=sp1F+((weights[0][1])*(CodBias[i])*mt*d); //w+cb cg20
                       spd1F=spd1F+(((weights[0][1])*(CodBias[i]))*((weights[0][1])*(CodBias[i])*mt*d));// (w+cb)(w+cb) cg20

                       tranvd64F=tranvd64F+((CodBias[i])*mt*d);  //****************
                       ttranvd64F=ttranvd64F+((CodBias[i])*(CodBias[i])*mt*d);
                       tranvdw64F=tranvdw64F+((weights[0][1])*(CodBias[i])*mt*d);  //****************
                       ttranvdw64F=ttranvdw64F+((weights[0][1])*(weights[0][1])*(CodBias[i])*(CodBias[i])*mt*d);
                                         }

                     }
                 }
                 if (i!=ii){
                 if (((codonpos==3)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){
                     if (transition(GCodes[4][i], GCodes[4][ii]))   ttrans=ttrans+1;
                    if (transversion(GCodes[4][i], GCodes[4][ii])) ttranv=ttranv+1;
                    nts3=nts3+1;                                         //nor *  //nor * and nor *
                            }
       if (((codonpos==2)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {
                   if (transition(GCodes[3][i], GCodes[3][ii]))   ttrans=ttrans+1;
                    if (transversion(GCodes[3][i], GCodes[3][ii])) ttranv=ttranv+1;
                   nts2=nts2+1;
                }
       if (((codonpos==1)||(codonpos==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {
                        if (transition(GCodes[2][i], GCodes[2][ii]))  ttrans=ttrans+1;
                          if (transversion(GCodes[2][i], GCodes[2][ii])) ttranv=ttranv+1;
                           nts1=nts1+1;
                     }
                        }

             }
              sT4=sT4+((sumt3F+sumt2F+sumt1F)*(sumt3F+sumt2F+sumt1F));             //le term 2 dans T4
           sT41=sT41+(sumt1F*sumt1F);           //le term 2 dans T4
           sT42=sT42+(sumt2F*sumt2F);            //le term 2 dans T4
           sT43=sT43+(sumt3F*sumt3F);             //le term 2 dans T4

           sT4w=sT4w+((sp1F+sp2F+sp3F)*(sp1F+sp2F+sp3F));             //le term 2 dans T4
           sT4w1=sT4w1+(sp1F*sp1F);          //le term 2 dans T4
           sT4w2=sT4w2+(sp2F*sp2F);          //le term 2 dans T4
           sT4w3=sT4w3+(sp3F*sp3F);          //le term 2 dans T4

           sT4wtrans=sT4wtrans+(transdw64F*transdw64F); //le term T3
           sT4trans=sT4trans+(transd64F*transd64F); //le term T3
           sT4wtranv=sT4wtranv+(tranvdw64F*tranvdw64F); //le term T3
           sT4tranv=sT4tranv+(tranvd64F*tranvd64F); //le term T3
           sT4aa=sT4aa+(spaaF*spaaF);

            sT3aa=sT3aa+((spaaF*spaaF)-spaa2F);
           sT3w1=sT3w1+((sp1F*sp1F)-spd1F); //le term T3
           sT3w2=sT3w2+((sp2F*sp2F)-spd2F); //le term T3
           sT3w3=sT3w3+((sp3F*sp3F)-spd3F); //le term T3

           sT3w=sT3w+(((sp1F+sp2F+sp3F)*(sp1F+sp2F+sp3F))-(spd1F+spd2F+spd3F));
           sT31=sT31+((sumt1F*sumt1F)-sumw1F); //le term T3
           sT32=sT32+((sumt2F*sumt2F)-sumw2F); //le term T3
           sT33=sT33+((sumt3F*sumt3F)-sumw3F); //le term T3

           sT3=sT3+(((sumt3F+sumt2F+sumt1F)*(sumt3F+sumt2F+sumt1F))-(sumw3F+sumw2F+sumw1F));
           sT3wtrans=sT3wtrans+((transdw64F*transdw64F)-ttransdw64F); //le term T3
           sT3wtranv=sT3wtranv+((tranvdw64F*tranvdw64F)-ttranvdw64F); //le term T3
           sT3trans=sT3trans+((transd64F*transd64F)-ttransd64F); //le term T3
           sT3tranv=sT3tranv+((tranvd64F*tranvd64F)-ttranvd64F); //le term T3

      }
      vector< double > res;

      res.push_back((sumt2+sumt1+sumt3));
      res.push_back((sp2+sp1+sp3));    //media global  0   cb +w
        //cb
      res.push_back(spaa);                //faa+w
      res.push_back(sumt1);                //cb1
      res.push_back(sumt2);                //cb2
      res.push_back(sumt3);                 //cb3
      res.push_back(sp1);                   //cb+w1
      res.push_back(sp2);                   //cb+w2
      res.push_back(sp3);                   //cb+w3
      res.push_back(transd64);
      res.push_back(transdw64);
      res.push_back(tranvd64);
      res.push_back(tranvdw64);
     //12

      res.push_back(2*(sumw2+sumw1+sumw3));
       res.push_back(2*(spd2+spd1+spd3));          // cb+w * cb+w variance term 13
             //cb*cb
       res.push_back(2*spaa2);                      // faa*faa
       res.push_back(2*sumw1);                      //cb*cb1
       res.push_back(2*sumw2);                        //cb*cb2
       res.push_back(2*sumw3);                       //cb*cb3
       res.push_back(2*spd1);                        //
       res.push_back(2*spd2);
       res.push_back(2*spd3);
       res.push_back(2*ttransd64);
       res.push_back(2*ttransdw64);
        res.push_back(2*ttranvd64);
       res.push_back(2*ttranvdw64);  //25

      res.push_back(((sumt3+sumt2+sumt1)*(sumt3+sumt2+sumt1))+(2*(sumw3+sumw2+sumw1))-(4*sT4)); //          26    // pour le T4
      res.push_back(((sp3+sp2+sp1)*(sp3+sp2+sp1))+(2*(spd3+spd2+spd1))-(4*sT4w)); //       27
      res.push_back((spaa*spaa)+(2*spaa2)-(4*sT4aa));                             //   28
      res.push_back((sumt1*sumt1)+(2*sumw1)-(4*sT41));    //                                   29
      res.push_back((sumt2*sumt2)+(2*sumw2)-(4*sT42));   //                                    30
      res.push_back((sumt3*sumt3)+(2*sumw3)-(4*sT43));   //                                    31
      res.push_back((sp1*sp1)+(2*spd1)-(4*sT4w1));    //                                   32
      res.push_back((sp2*sp2)+(2*spd2)-(4*sT4w2));   //                                    33
      res.push_back((sp3*sp3)+(2*spd3)-(4*sT4w3));    //                                   34
      res.push_back((transd64*transd64)+(2*ttransd64)-(4*sT4trans));   //                                   35
      res.push_back((transdw64*transdw64)+(2*ttransdw64)-(4*sT4wtrans));  //                                   36
      res.push_back((tranvd64*tranvd64)+(2*ttranvd64)-(4*sT4tranv));   //                                   37
      res.push_back((tranvdw64*tranvdw64)+(2*ttranvdw64)-(4*sT4wtranv));   //  38


      res.push_back(4*sT3); //co, co*co cg20  -- daa  0              39       // pour le T3
      res.push_back(4*sT3w); //w+da cg20    -- daa                    40     //T3
      res.push_back(4*sT3aa);                        //41
      res.push_back(4*sT31);    //                                   42
      res.push_back(4*sT32);   //                                    43
      res.push_back(4*sT33);   //                                    44
      res.push_back(4*sT3w1);    //                                   45
      res.push_back(4*sT3w2);   //                                    46
      res.push_back(4*sT3w3);    //                                   47
      res.push_back(4*sT3trans);   //                                   48
      res.push_back(4*sT3wtrans);  //                                   49
      res.push_back(4*sT3tranv);   //                                   50
      res.push_back(4*sT3wtranv);   //                                  51

       res.push_back((nts1+nts2+nts3)); //(w+cb)*(w+cb) cg21 52
       res.push_back((nts1+nts2+nts3));  // co*co  53
      res.push_back((nts1+nts2+nts3));  //              54
      res.push_back((nts1)); //w*w cg21              55
      res.push_back((nts2));
      res.push_back((nts3));
      res.push_back((nts1));
      res.push_back((nts2));
      res.push_back((nts3)); //
      res.push_back(ttrans); //
      res.push_back(ttrans); //
      res.push_back(ttranv); //
      res.push_back(ttranv);  //64

    return res;
 }
//permGenCodes2 
 vector< double > GraphGC::permGenRobustness2(string alp, const vector <double>  &CodBias, const vector <double>  &Faas,  const vector <string> &GCodes) {
  //weight factors for global mean and variance
  int codonpos=0;
  pair<int,char> po, po2;

  int w, mt, nt=0, nt1=0, nt2=0, nt3=0;
  int ntt=0, ntt1=0, ntt2=0, ntt3=0;
  int nts1=0, nts2=0, nts3=0;
  int ntts1=0, ntts2=0, ntts3=0;

  double sumt1=0, sumt2=0, sumt3=0, sumw1=0, sumw2=0, sumw3=0;
  double sp1=0, sp2=0, sp3=0, spd1=0, spd2=0, spd3=0;
   double sumt1F=0, sumt2F=0, sumt3F=0, sumw1F=0, sumw2F=0, sumw3F=0;
  double sp1F=0, sp2F=0, sp3F=0, spd1F=0, spd2F=0, spd3F=0;

  double ttrans=0, ttranv=0;

  double sumc1=0, sumc2=0, sumc3=0, sumww1=0, sumww2=0, sumww3=0, trvwd20=0, trvd20=0, trswd20=0, trsd20=0, spao=0;

 double  transd20=0, ttransd20=0, transdw20=0, ttransdw20=0,tranvd20=0, ttranvd20=0, tranvdw20=0, ttranvdw20=0;
 double  transd20F=0, ttransd20F=0, transdw20F=0, ttransdw20F=0,tranvd20F=0, ttranvd20F=0, tranvdw20F=0, ttranvdw20F=0;
 double sT3w1=0,sT3w2=0,sT3w3=0, sT3w=0, sT31=0, sT32=0,sT33=0, sT3=0, sT3wtrans=0,sT3trans=0, sT3wtranv=0, sT3tranv=0;
 double sT4w3=0,sT4w2=0,sT4w1=0, sT4w=0, sT41=0, sT42=0,sT43=0, sT4=0, sT4wtrans=0,sT4trans=0, sT4wtranv=0, sT4tranv=0;

double d=0, spaa=0, spaa2=0, spaaF=0, spaa2F=0, sT4aa=0, sT3aa=0;

  //cout<<this->aa.length()<<endl;
 for (unsigned int i = 0; i < GCodes[0].size(); i++)
   {
    po=AAnIdentify(i, alp);
    spd1F=0, spd2F=0, spd3F=0; sp1F=0, sp2F=0, sp3F=0;
    sumt1F=0, sumt2F=0,sumt3F=0, sumw3F=0, sumw2F=0, sumw1F=0;
    transdw20F=0, ttransdw20F=0, tranvdw20F=0, ttranvdw20F=0;
    spaaF=0, spaa2F=0;
    for (unsigned int ii = 0; ii < GCodes[0].size(); ii++)
     {
    if (GCodes[0][i]!=GCodes[0][ii]) {
       po2=AAnIdentify(ii, alp);  //3rst pos
      // cout<<((Paa[po.first]-Paa[po2.first])*(Paa[po.first]-Paa[po2.first]))<<" ";
       if (((codonpos==3)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){
                  if (transition(GCodes[4][i], GCodes[4][ii])){
                      if (STOP(GCodes[0][i], GCodes[0][ii])) { d=0;   }//* or *
                         else {d=1; ttrans=ttrans+1; }
                    //d=1;
                     if (CodBias[i]==1) mt=0; else
                        mt=1;  // d=1;
                      // ttrans=ttrans+1;

                       sumt3=sumt3+((CodBias[i])*d*mt);  //cb cg20
                       sumw3=sumw3+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp3=sp3+((weights[2][0])*(CodBias[i])*d*mt); //w+cb cg20
                       spd3=spd3+(((weights[2][0])*(CodBias[i]))*((weights[2][0])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       spaa=spaa+((weights[2][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2=spaa2+((weights[2][0])*(Faas[i])*(weights[2][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaaF=spaaF+((weights[2][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2F=spaa2F+((weights[2][0])*(Faas[i])*(weights[2][0])*(Faas[i])*d*mt); //w+cb cg20

                       transd20=transd20+((CodBias[i])*d*mt);  //****************
                       ttransd20=ttransd20+((CodBias[i])*(CodBias[i])*d*mt);
                       transdw20=transdw20+((weights[2][0])*(CodBias[i])*d*mt);  //****************
                       ttransdw20=ttransdw20+((weights[2][0])*(weights[2][0])*(CodBias[i])*(CodBias[i])*d*mt);

                       sumt3F=sumt3F+((CodBias[i])*d*mt);  //cb cg20
                       sumw3F=sumw3F+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp3F=sp3F+((weights[2][0])*(CodBias[i])*d*mt); //w+cb cg20
                       spd3F=spd3F+(((weights[2][0])*(CodBias[i]))*((weights[2][0])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       transd20F=transd20F+((CodBias[i])*d*mt);  //****************
                       ttransd20F=ttransd20F+((CodBias[i])*(CodBias[i])*d*mt);
                       transdw20F=transdw20F+((weights[2][0])*(CodBias[i])*d*mt);  //****************
                       ttransdw20F=ttransdw20F+((weights[2][0])*(weights[2][0])*(CodBias[i])*(CodBias[i])*d*mt);

                        sumc3=sumc3+((CodBias[i])*d*mt); //co, co*co cg20   global mean
                       sumww3=sumww3+((CodBias[i])*(weights[2][0])*d*mt); //w cg20

                       trswd20=trswd20+((CodBias[i])*(weights[2][0])*d*mt);
                       trsd20=trsd20+((CodBias[i])*d*mt);
                       spao=spao+((weights[2][0])*(Faas[i])*d*mt);

                                           }
                  if (transversion(GCodes[4][i], GCodes[4][ii])){
                       if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                          else {d=1; ttranv=ttranv+1;  }
                       //d=1;
                         if (CodBias[i]==1) mt=0; else
                        mt=1;
                     //  ttranv=ttranv+1;
                       sumt3=sumt3+((CodBias[i])*d*mt);  //cb cg20
                       sumw3=sumw3+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp3=sp3+((weights[2][1])*(CodBias[i])*d*mt); //w+cb cg20

                       spd3=spd3+(((weights[2][1])*(CodBias[i]))*((weights[2][1])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20
                       spaa=spaa+((weights[2][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2=spaa2+((weights[2][1])*(Faas[i])*(weights[2][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaaF=spaaF+((weights[2][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2F=spaa2F+((weights[2][1])*(Faas[i])*(weights[2][1])*(Faas[i])*d*mt); //w+cb cg20

                       tranvd20=tranvd20+((CodBias[i])*d*mt);  //****************
                       ttranvd20=ttranvd20+((CodBias[i])*(CodBias[i])*d*mt);
                       tranvdw20=tranvdw20+((weights[2][1])*(CodBias[i])*d*mt);  //****************
                       ttranvdw20=ttranvdw20+((weights[2][1])*(weights[2][1])*(CodBias[i])*(CodBias[i])*d*mt);

                       sumt3F=sumt3F+((CodBias[i])*d*mt);  //cb cg20
                       sumw3F=sumw3F+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp3F=sp3F+((weights[2][1])*(CodBias[i])*d*mt); //w+cb cg20
                       spd3F=spd3F+(((weights[2][1])*(CodBias[i]))*((weights[2][1])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       tranvd20F=tranvd20F+((CodBias[i])*d*mt);  //****************
                       ttranvd20F=ttranvd20F+((CodBias[i])*(CodBias[i])*d*mt);
                       tranvdw20F=tranvdw20F+((weights[2][1])*(CodBias[i])*d*mt);  //****************
                       ttranvdw20F=ttranvdw20F+((weights[2][1])*(weights[2][1])*(CodBias[i])*(CodBias[i])*d*mt);

                       sumc3=sumc3+((CodBias[i])*d*mt); //co, co*co cg20   global mean
                       sumww3=sumww3+((CodBias[i])*(weights[2][1])*d*mt); //w cg20
                       trvwd20=trvwd20+((CodBias[i])*(weights[2][1])*d*mt);
                       trvd20=trvd20+((CodBias[i])*d*mt);
                       spao=spao+((weights[2][1])*(Faas[i])*d*mt);
                         }
                   nts3=nts3+1;
                    if (!(STOP(GCodes[0][i], GCodes[0][ii]))) {
                       nt3=nt3+1; }
                        //nor *  //nor * and nor *
              }
       if (((codonpos==2)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {
                 if (transition(GCodes[3][i], GCodes[3][ii])){
                     if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                       else  {d=1; ttrans=ttrans+1;}
                    //  ttrans=ttrans+1;
                        if (CodBias[i]==1) mt=0; else
                        mt=1;//d=1;

                       sumt2=sumt2+((CodBias[i])*d*mt);  //cb cg20
                       sumw2=sumw2+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp2=sp2+((weights[1][0])*(CodBias[i])*d*mt); //w+cb cg20
                       spd2=spd2+(((weights[1][0])*(CodBias[i]))*((weights[1][0])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       spaa=spaa+((weights[1][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2=spaa2+((weights[1][0])*(Faas[i])*(weights[1][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaaF=spaaF+((weights[1][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2F=spaa2F+((weights[1][0])*(Faas[i])*(weights[1][0])*(Faas[i])*d*mt); //w+cb cg20

                       transd20=transd20+((CodBias[i])*d*mt);  //****************
                       ttransd20=ttransd20+((CodBias[i])*(CodBias[i])*d*mt);
                       transdw20=transdw20+((weights[1][0])*(CodBias[i])*d*mt);  //****************
                       ttransdw20=ttransdw20+((weights[1][0])*(weights[1][0])*(CodBias[i])*(CodBias[i])*d*mt);

                       sumt2F=sumt2F+((CodBias[i])*d*mt);  //cb cg20
                       sumw2F=sumw2F+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp2F=sp2F+((weights[1][0])*(CodBias[i])*d*mt); //w+cb cg20
                       spd2F=spd2F+(((weights[1][0])*(CodBias[i]))*((weights[1][0])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       transd20F=transd20F+((CodBias[i])*d*mt);  //****************
                       ttransd20F=ttransd20F+((CodBias[i])*(CodBias[i])*d*mt);
                       transdw20F=transdw20F+((weights[1][0])*(CodBias[i])*d*mt);  //****************
                       ttransdw20F=ttransdw20F+((weights[1][0])*(weights[1][0])*(CodBias[i])*(CodBias[i])*d*mt);

                       sumc2=sumc2+((CodBias[i])*d*mt); //co, co*co cg20   global mean
                       sumww2=sumww2+((CodBias[i])*(weights[1][0])*d*mt); //w cg20
                       trswd20=trswd20+((CodBias[i])*(weights[1][0])*d*mt);
                       trsd20=trsd20+((CodBias[i])*d*mt);
                       spao=spao+((weights[1][0])*(Faas[i])*d*mt);
                               }
                 if (transversion(GCodes[3][i], GCodes[3][ii])){
                   if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                      else{ d=1; ttranv=ttranv+1; }
                    // ttranv=ttranv+1;
                       if (CodBias[i]==1) mt=0; else
                        mt=1; //d=1;

                       sumt2=sumt2+((CodBias[i])*d*mt);  //cb cg20
                       sumw2=sumw2+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp2=sp2+((weights[1][1])*(CodBias[i])*d*mt); //w+cb cg20
                       spd2=spd2+(((weights[1][1])*(CodBias[i]))*((weights[1][1])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                        spaa=spaa+((weights[1][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2=spaa2+((weights[1][1])*(Faas[i])*(weights[1][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaaF=spaaF+((weights[1][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2F=spaa2F+((weights[1][1])*(Faas[i])*(weights[1][1])*(Faas[i])*d*mt); //w+cb cg20

                       tranvd20=tranvd20+((CodBias[i])*d*mt);  //****************
                       ttranvd20=ttranvd20+((CodBias[i])*(CodBias[i])*d*mt);
                       tranvdw20=tranvdw20+((weights[1][1])*(CodBias[i])*d*mt);  //****************
                       ttranvdw20=ttranvdw20+((weights[1][1])*(weights[1][1])*(CodBias[i])*(CodBias[i])*d*mt);

                       sumt2F=sumt2F+((CodBias[i])*d*mt);  //cb cg20
                       sumw2F=sumw2F+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp2F=sp2F+((weights[1][1])*(CodBias[i])*d*mt); //w+cb cg20
                       spd2F=spd2F+(((weights[1][1])*(CodBias[i]))*((weights[1][1])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       tranvd20F=tranvd20F+((CodBias[i])*d*mt);  //****************
                       ttranvd20F=ttranvd20F+((CodBias[i])*(CodBias[i])*d*mt);
                       tranvdw20F=tranvdw20F+((weights[1][1])*(CodBias[i])*d*mt);  //****************
                       ttranvdw20F=ttranvdw20F+((weights[1][1])*(weights[1][1])*(CodBias[i])*(CodBias[i])*d*mt);

                       sumc2=sumc2+((CodBias[i])*d*mt); //co, co*co cg20   global mean
                       sumww2=sumww2+((CodBias[i])*(weights[1][1])*d*mt); //w cg20
                       trvwd20=trvwd20+((CodBias[i])*(weights[1][1])*d*mt);
                       trvd20=trvd20+((CodBias[i])*d*mt);
                       spao=spao+((weights[1][1])*(Faas[i])*d*mt);
                               }
                 nts2=nts2+1;
                 if (!(STOP(GCodes[0][i], GCodes[0][ii]))) {
                       nt2=nt2+1; }
                        //nor *  //nor * and nor *
             }
       if (((codonpos==1)||(codonpos==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {
                // cout<<po2.first<<" ";
                 if (transition(GCodes[2][i], GCodes[2][ii])){
                     if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                       else {ttrans=ttrans+1;d=1; }
                      // ttrans=ttrans+1;//d=1;
                         if (CodBias[i]==1) mt=0; else
                        mt=1; //cg20

                       sumt1=sumt1+((CodBias[i])*d*mt);  //cb cg20
                       sumw1=sumw1+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp1=sp1+((weights[0][0])*(CodBias[i])*d*mt); //w+cb cg20
                       spd1=spd1+(((weights[0][0])*(CodBias[i]))*((weights[0][0])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       transd20=transd20+((CodBias[i])*d*mt);  //****************
                       ttransd20=ttransd20+((CodBias[i])*(CodBias[i])*d*mt);
                       transdw20=transdw20+((weights[0][0])*(CodBias[i])*d*mt);  //****************
                       ttransdw20=ttransdw20+((weights[0][0])*(weights[0][0])*(CodBias[i])*(CodBias[i])*d*mt);

                        spaa=spaa+((weights[0][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2=spaa2+((weights[0][0])*(Faas[i])*(weights[0][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaaF=spaaF+((weights[0][0])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2F=spaa2F+((weights[0][0])*(Faas[i])*(weights[0][0])*(Faas[i])*d*mt); //w+cb cg20

                       sumt1F=sumt1F+((CodBias[i])*d*mt);  //cb cg20
                       sumw1F=sumw1F+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp1F=sp1F+((weights[0][0])*(CodBias[i])*d*mt); //w+cb cg20
                       spd1F=spd1F+(((weights[0][0])*(CodBias[i]))*((weights[0][0])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       transd20F=transd20F+((CodBias[i])*d*mt);  //****************
                       ttransd20F=ttransd20F+((CodBias[i])*(CodBias[i])*d*mt);
                       transdw20F=transdw20F+((weights[0][0])*(CodBias[i])*d*mt);  //****************
                       ttransdw20F=ttransdw20F+((weights[0][0])*(weights[0][0])*(CodBias[i])*(CodBias[i])*d*mt);

                         sumc1=sumc1+((CodBias[i])*d*mt); //co, co*co cg20   global mean
                       sumww1=sumww1+((CodBias[i])*(weights[0][0])*d*mt); //w cg20
                       trswd20=trswd20+((CodBias[i])*(weights[0][0])*d*mt);
                       trsd20=trsd20+((CodBias[i])*d*mt);
                       spao=spao+((weights[0][0])*(Faas[i])*d*mt);
                        }
                 if (transversion(GCodes[2][i], GCodes[2][ii])){
                    if (STOP(GCodes[0][i], GCodes[0][ii])) d=0; //* or *
                      else {d=1; ttranv=ttranv+1;}
                      //d=1; ttranv=ttranv+1;
                       if (CodBias[i]==1) mt=0; else
                        mt=1;

                       sumt1=sumt1+((CodBias[i])*d*mt);  //cb cg20
                        sp1=sp1+((weights[0][1])*(CodBias[i])*d*mt); //w+cb cg2
                       sumw1=sumw1+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       spd1=spd1+(((weights[0][1])*(CodBias[i]))*((weights[0][1])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       tranvd20=tranvd20+((CodBias[i])*d*mt);  //****************
                       ttranvd20=ttranvd20+((CodBias[i])*(CodBias[i])*d*mt);
                       tranvdw20=tranvdw20+((weights[0][1])*(CodBias[i])*d*mt);  //****************
                       ttranvdw20=ttranvdw20+((weights[0][1])*(weights[0][1])*(CodBias[i])*(CodBias[i])*d*mt);

                        spaa=spaa+((weights[0][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2=spaa2+((weights[0][1])*(Faas[i])*(weights[0][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaaF=spaaF+((weights[0][1])*(Faas[i])*d*mt); //w+cb cg20
                       spaa2F=spaa2F+((weights[0][1])*(Faas[i])*(weights[0][1])*(Faas[i])*d*mt); //w+cb cg20

                       sumt1F=sumt1F+((CodBias[i])*d*mt);  //cb cg20
                       sumw1F=sumw1F+((CodBias[i])*(CodBias[i])*d*mt); //w cg20
                       sp1F=sp1F+((weights[0][1])*(CodBias[i])*d*mt); //w+cb cg20
                       spd1F=spd1F+(((weights[0][1])*(CodBias[i]))*((weights[0][1])*(CodBias[i])*d*mt));// (w+cb)(w+cb) cg20

                       tranvd20F=tranvd20F+((CodBias[i])*d*mt);  //****************
                       ttranvd20F=ttranvd20F+((CodBias[i])*(CodBias[i])*d*mt);
                       tranvdw20F=tranvdw20F+((weights[0][1])*(CodBias[i])*d*mt);  //****************
                       ttranvdw20F=ttranvdw20F+((weights[0][1])*(weights[0][1])*(CodBias[i])*(CodBias[i])*d*mt);

                     sumc1=sumc1+((CodBias[i])*d*mt); //co, co*co cg20   global mean
                       sumww1=sumww1+((CodBias[i])*(weights[0][1])*d*mt); //w cg20
                       trvwd20=trvwd20+((CodBias[i])*(weights[0][1])*d*mt);
                       trvd20=trvd20+((CodBias[i])*d*mt);
                       spao=spao+((weights[0][1])*(Faas[i])*d*mt);
                    }
                  nts1=nts1+1;
                    if (!(STOP(GCodes[0][i], GCodes[0][ii]))) {
                        nt1=nt1+1; }
                         //nor * and nor *
                }
            }
      else if ((GCodes[0][i]==GCodes[0][ii])&&(i!=ii)) {
       po2=AAnIdentify(ii,alp);  //3rst pos
       if (((codonpos==3)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]!=GCodes[4][ii])){
                 ntts3=ntts3+1;
                if (!(STOP(GCodes[0][i], GCodes[0][ii]))){
                         if (transversion(GCodes[4][i], GCodes[4][ii])) ttranv=ttranv+1;
                       if (transition(GCodes[4][i], GCodes[4][ii]))  ttrans=ttrans+1;
                        ntt3=ntt3+1;   }
                        //nor * //nor * and nor *
        }
       if (((codonpos==2)||(codonpos==0))&&(GCodes[2][i]==GCodes[2][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[3][i]!=GCodes[3][ii])) {
              ntts2=ntts2+1;
              if (!(STOP(GCodes[0][i], GCodes[0][ii]))){
                         if (transversion(GCodes[3][i], GCodes[3][ii])) ttranv=ttranv+1;
                       if (transition(GCodes[3][i], GCodes[3][ii]))  ttrans=ttrans+1;
                        ntt2=ntt2+1;}
                        //nor *  //nor * and nor *
             }
       if (((codonpos==1)||(codonpos==0))&&(GCodes[3][i]==GCodes[3][ii])&&(GCodes[4][i]==GCodes[4][ii])&&(GCodes[2][i]!=GCodes[2][ii])) {
                // cout<<po2.first<<" ";
                  ntts1=ntts1+1;
                  if (!(STOP(GCodes[0][i], GCodes[0][ii]))){
                         if (transversion(GCodes[2][i], GCodes[2][ii])) ttranv=ttranv+1;
                       if (transition(GCodes[2][i], GCodes[2][ii]))  ttrans=ttrans+1;
                        ntt1=ntt1+1;
                        }
                    //nor *  //nor * and nor *
                }
            }
        }
           sT4=sT4+((sumt3F+sumt2F+sumt1F)*(sumt3F+sumt2F+sumt1F));             //le term 2 dans T4
           sT41=sT41+(sumt1F*sumt1F);           //le term 2 dans T4
           sT42=sT42+(sumt2F*sumt2F);            //le term 2 dans T4
           sT43=sT43+(sumt3F*sumt3F);             //le term 2 dans T4

           sT4w=sT4w+((sp1F+sp2F+sp3F)*(sp1F+sp2F+sp3F));             //le term 2 dans T4
           sT4w1=sT4w1+(sp1F*sp1F);          //le term 2 dans T4
           sT4w2=sT4w2+(sp2F*sp2F);          //le term 2 dans T4
           sT4w3=sT4w3+(sp3F*sp3F);          //le term 2 dans T4

           sT4wtrans=sT4wtrans+(transdw20F*transdw20F); //le term T3
           sT4trans=sT4trans+(transd20F*transd20F); //le term T3
           sT4wtranv=sT4wtranv+(tranvdw20F*tranvdw20F); //le term T3
           sT4tranv=sT4tranv+(tranvd20F*tranvd20F); //le term T3
           sT4aa=sT4aa+(spaaF*spaaF);

            sT3aa=sT3aa+((spaaF*spaaF)-spaa2F);
           sT3w1=sT3w1+((sp1F*sp1F)-spd1F); //le term T3
           sT3w2=sT3w2+((sp2F*sp2F)-spd2F); //le term T3
           sT3w3=sT3w3+((sp3F*sp3F)-spd3F); //le term T3
           sT3w=sT3w+(((sp1F+sp2F+sp3F)*(sp1F+sp2F+sp3F))-(spd1F+spd2F+spd3F));

           sT31=sT31+((sumt1F*sumt1F)-sumw1F); //le term T3
           sT32=sT32+((sumt2F*sumt2F)-sumw2F); //le term T3
           sT33=sT33+((sumt3F*sumt3F)-sumw3F); //le term T3
           sT3=sT3+(((sumt3F+sumt2F+sumt1F)*(sumt3F+sumt2F+sumt1F))-(sumw3F+sumw2F+sumw1F));

           sT3wtrans=sT3wtrans+((transdw20F*transdw20F)-ttransdw20F); //le term T3
           sT3wtranv=sT3wtranv+((tranvdw20F*tranvdw20F)-ttranvdw20F); //le term T3
           sT3trans=sT3trans+((transd20F*transd20F)-ttransd20F); //le term T3
           sT3tranv=sT3tranv+((tranvd20F*tranvd20F)-ttranvd20F); //le term T3

      }
      vector< double > res;

      res.push_back((sumt1+sumt2+sumt3));    //media global  0
      res.push_back((sp1+sp2+sp3));
      res.push_back(spaa);
      res.push_back(sumt1);
      res.push_back(sumt2);
      res.push_back(sumt3);
      res.push_back(sp1);
      res.push_back(sp2);
      res.push_back(sp3);
       res.push_back(transd20);
      res.push_back(transdw20);
      res.push_back(tranvd20);    //12
      res.push_back(ttranvdw20);
                        // variance term 13
       res.push_back(2*(sumw2+sumw1+sumw3)); //13
       res.push_back(2*(spd2+spd1+spd3));
       res.push_back(2*spaa2);
       res.push_back(2*sumw1);
       res.push_back(2*sumw2);
       res.push_back(2*sumw3);
       res.push_back(2*spd1);
       res.push_back(2*spd2);
       res.push_back(2*spd3);
       res.push_back(2*ttransd20);
       res.push_back(2*ttransdw20);
        res.push_back(2*ttranvd20);
       res.push_back(2*ttranvdw20);  //25

      res.push_back(((sumt3+sumt2+sumt1)*(sumt3+sumt2+sumt1))+(2*(sumw3+sumw2+sumw1))-(4*sT4)); //       26        // pour le T4
      res.push_back(((sp3+sp2+sp1)*(sp3+sp2+sp1))+(2*(spd3+spd2+spd1))-(4*sT4w)); //  27
      res.push_back((spaa*spaa)+(2*spaa2)-(4*sT4aa));                         //28
      res.push_back((sumt1*sumt1)+(2*sumw1)-(4*sT41));    //                         29
      res.push_back((sumt2*sumt2)+(2*sumw2)-(4*sT42));   //                                    30
      res.push_back((sumt3*sumt3)+(2*sumw3)-(4*sT43));   //                                    31
      res.push_back((sp1*sp1)+(2*spd1)-(4*sT4w1));    //                                   32
      res.push_back((sp2*sp2)+(2*spd2)-(4*sT4w2));   //                                    33
      res.push_back((sp3*sp3)+(2*spd3)-(4*sT4w3));    //                                   34
      res.push_back((transd20*transd20)+(2*ttransd20)-(4*sT4trans));   //                                   35
      res.push_back((transdw20*transdw20)+(2*ttransdw20)-(4*sT4wtrans));  //                                   36
      res.push_back((tranvd20*tranvd20)+(2*ttranvd20)-(4*sT4tranv));   //                                   37
      res.push_back((tranvdw20*tranvdw20)+(2*ttranvdw20)-(4*sT4wtranv));   //  38

      res.push_back(4*sT3); //co, co*co cg20  -- daa  0              39     // pour le T3
      res.push_back(4*sT3w); //w+da cg20    -- daa                    40    //T3
      res.push_back(4*sT3aa);                //                       41
      res.push_back(4*sT31);    //                                   42
      res.push_back(4*sT32);   //                                    43
      res.push_back(4*sT33);   //                                    44
      res.push_back(4*sT3w1);    //                                   45
      res.push_back(4*sT3w2);   //                                    46
      res.push_back(4*sT3w3);    //                                   47
      res.push_back(4*sT3trans);   //                                   48
      res.push_back(4*sT3wtrans);  //                                   49
      res.push_back(4*sT3tranv);   //                                   50
      res.push_back(4*sT3wtranv);   //                                  51

       res.push_back((nt1+nt2+nt3+ntt1+ntt2+ntt3)); //(w+cb)*(w+cb) 52
       res.push_back((nt1+nt2+nt3+ntt1+ntt2+ntt3));  // co*co
      res.push_back((nt1+nt2+nt3+ntt1+ntt2+ntt3));
      res.push_back((nt1+ntt1)); //w*w cg21
      res.push_back((nt2+ntt2));
      res.push_back((nt3+ntt3));
      res.push_back((nt1+ntt1));
      res.push_back((nt2+ntt2));
      res.push_back((nt3+ntt3)); //
      res.push_back(ttrans); //
      res.push_back(ttrans); //
      res.push_back(ttranv); //
      res.push_back(ttranv);  //64

    return res;
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


