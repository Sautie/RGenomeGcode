#include "AuxRandpp.h"
using namespace std;

vector<vector<int> >  randppp(vector<int> s, int m) {
     vector<vector<int> > randm;
     //vector<int> randv;
     for (int i = 0; i < m; i++)
     {
     unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
     shuffle(s.begin(), s.end(), default_random_engine(seed));
     randm.push_back(s);
           }
   return  randm;
}

vector<vector<int> >  randInvSwap(vector<vector<int> >  s, int m) {
  for (int y = 0; y < m; y++)
     {
   for (int i = 0; i < s.size(); i++)
     {
           srand((int)time(0));
           int r = (rand() % s[i].size()) ;
           int l = (rand() % s[i].size()) ;
           if (l>r)
               swap(l, r);
        int D=l+(r/2)+1;
           for (int j = l; j <D; j++)
              swap(s[i][j], s[i][D-j-1]);
                    }
       }
   return s;
}

long int  funcPos(vector<int> s){
    int p;
    long int ss=0;
    for (int i = 0; i < s.size(); i++)
     {
         p=1;
         for(int ii=0; ii<i; ii++)
            p *= 2;
        ss=ss+(s[i]*p);
         }
    return ss;
}

vector<vector<int> >  RandWOutRep(vector<int> s, int m) {
    vector<vector<int> > randm;
    set<long int>  rands;
    long int rs;
    int i=0;
     while (i<m)
     {
     unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
     shuffle(s.begin(), s.end(), default_random_engine(seed));
     rs=funcPos(s);
     if (rands.find(rs) == rands.end()) {
           randm.push_back(s);
           rands.insert(rs);
                            }
              i++;
           }
return randm;
   }
void EqGenerator(string dest1, vector<string> vect20, string dest2, vector<string> vect64)
{

    string file;
    ofstream outfile0(dest2.c_str());
    ofstream outfile1(dest1.c_str());
     outfile0<<"(mod0.lmer <- glmer(st ~ (1|P/C/G), data =score20, family = binomial))"<<endl;
     outfile1<<"(mod0.lmer <- glmer(st ~ (1|P/C/G), data =score64p3, family = binomial))"<<endl;
    int L=vect64.size();
    for (unsigned int p=100;p<L;p++)  {
          outfile0<<"(gm"<<to_string(p)<<"<- glmer(st ~ "<<vect20[p]<<" + (1|P/C/G), data = score20, family = binomial))"<<endl;
          outfile0<<"sqrt(diag(vcov(gm"<<to_string(p)<<")))"<<endl;
          outfile0<<"anova(mod0.lmer,gm"<<to_string(p)<<")"<<endl;
          outfile1<<"(gm"<<to_string(p)<<"<- glmer(st ~ "<<vect64[p]<<" + (1|P/C/G), data = score64p3, family = binomial))"<<endl;
          outfile1<<"sqrt(diag(vcov(gm"<<to_string(p)<<")))"<<endl;
          outfile1<<"anova(mod0.lmer,gm"<<to_string(p)<<")"<<endl;
         }

 }

void read_folder(const string& name, vector <string> & list)
{
    string pattern(name);
    pattern.append("\\*");
    WIN32_FIND_DATA data;
    HANDLE hFind;
    if ((hFind = FindFirstFile(pattern.c_str(), &data)) != INVALID_HANDLE_VALUE) {
        do {
            list.push_back(data.cFileName);
			//cout<<data.cFileName<<endl;
        } while (FindNextFile(hFind, &data) != 0);
        FindClose(hFind);
    }
}
vector< string > GenCodes(string GCname) {
vector< string > gcodes;
if (GCname=="sgc") {
    gcodes.push_back("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
    gcodes.push_back("---M------**--*----M---------------M----------------------------");
    gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
    gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
    gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
    }
else if (GCname=="vmc"){
      gcodes.push_back("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG");   //The Vertebrate Mitochondrial Code
      gcodes.push_back("----------**--------------------MMMM----------**---M------------");
      gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
      gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
      gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
         }
else if (GCname=="ymc") {
     gcodes.push_back("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG");  //The Yeast Mitochondrial Code
     gcodes.push_back("----------**----------------------MM---------------M------------");
     gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
     gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
     gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
           }
else if (GCname=="mmc") {
       gcodes.push_back("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
       gcodes.push_back("--MM------**-------M------------MMMM---------------M------------");
       gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
       gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
       gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
           }
else if (GCname=="imc") {
      gcodes.push_back("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"); //The Invertebrate Mitochondrial Code
      gcodes.push_back("---M------**--------------------MMMM---------------M------------");
      gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
      gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
      gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
    }
else if (GCname=="cnc") {
  gcodes.push_back("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); // The Ciliate, Dasycladacean and Hexamita Nuclear Code
  gcodes.push_back("--------------*--------------------M----------------------------");
  gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
  gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
  gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
}

else if (GCname=="emc") {
  gcodes.push_back("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"); //The Echinoderm and Flatworm Mitochondrial Code
  gcodes.push_back("----------**-----------------------M---------------M------------");
   gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
   gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
  gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
}

else if (GCname=="enc") {
  gcodes.push_back("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //The Euplotid Nuclear Code
  gcodes.push_back("----------**-----------------------M----------------------------");
  gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
  gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
  gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
}

else if (GCname=="aync") {
  gcodes.push_back("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //The Alternative Yeast Nuclear Code (transl_table=12)
  gcodes.push_back("----------**--*----M---------------M----------------------------");
  gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
  gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
  gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
}

else if (GCname=="amc") {
   gcodes.push_back("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"); //The Ascidian Mitochondrial Code
   gcodes.push_back("---M------**----------------------MM---------------M------------");
   gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
   gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
   gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
     }

 else if (GCname=="afmc") {
     gcodes.push_back("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"); //The Alternative Flatworm Mitochondrial Code (transl_table=14)
     gcodes.push_back("-----------*-----------------------M----------------------------");
     gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
     gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
     gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
      }
   else if (GCname=="cmc") {
     gcodes.push_back("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); // Chlorophycean Mitochondrial Code
     gcodes.push_back("----------*---*--------------------M----------------------------");
     gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
     gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
     gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
               }
   else if (GCname=="tmc") {
      gcodes.push_back("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"); //Trematode Mitochondrial Code (transl_table=21)
      gcodes.push_back("----------**-----------------------M---------------M------------");
      gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
      gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
      gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
        }
    else if (GCname=="smc") {
        gcodes.push_back("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //Scenedesmus obliquus Mitochondrial Code (transl_table=22)
        gcodes.push_back("------*---*---*--------------------M----------------------------");
        gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
        gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
        gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
      }
  else if (GCname=="pmc") {
    gcodes.push_back("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"); //Pterobranchia Mitochondrial Code (transl_table=24)
    gcodes.push_back("---M------**-------M---------------M---------------M------------");
    gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
    gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
    gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
          }
   else if (GCname=="cnc") {
        gcodes.push_back("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
        gcodes.push_back("---M------**-----------------------M---------------M------------");
        gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
        gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
        gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
         }
    else if (GCname=="ptnc") {
      gcodes.push_back("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //Pachysolen tannophilus Nuclear Code (transl_table=26)
      gcodes.push_back("----------**--*----M---------------M----------------------------");
      gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
      gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
      gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
             }
    else if (GCname=="knc") {
       gcodes.push_back("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //Karyorelict Nuclear Code
       gcodes.push_back("--------------*--------------------M----------------------------");
       gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
       gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
       gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
                       }
 else if (GCname=="conc") {
   gcodes.push_back("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //Condylostoma Nuclear Code (transl_table=28)
   gcodes.push_back("----------**--*--------------------M----------------------------");
   gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
   gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
   gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
                   }
else if (GCname=="mnc") {
  gcodes.push_back("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");//Mesodinium Nuclear Code (transl_table=29)
  gcodes.push_back("--------------*--------------------M----------------------------");
  gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
  gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
  gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
                }
else if (GCname=="pnc") {
  gcodes.push_back("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); // Peritrich Nuclear Code (transl_table=30)
  gcodes.push_back("--------------*--------------------M----------------------------");
  gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
  gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
  gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
        }
else if (GCname=="bnc") {
         gcodes.push_back("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //Blastocrithidia Nuclear Code (transl_table=31)
         gcodes.push_back("----------**-----------------------M----------------------------");
         gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
         gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
         gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
                  }
else if (GCname=="cmuc") {
    gcodes.push_back("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG");// Cephalodiscidae Mitochondrial UAA-Tyr Code (transl_table=33)
    gcodes.push_back("---M-------*-------M---------------M---------------M------------");
    gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
    gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
    gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
                                   }
 else if (GCname=="tmc") {
  gcodes.push_back("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"); //Thraustochytrium Mitochondrial Code (transl_table=23)
  gcodes.push_back("--*-------**--*-----------------M--M---------------M------------");
  gcodes.push_back("TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG");
  gcodes.push_back("TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG");
  gcodes.push_back("TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
                                 }
return gcodes;
}
	vector<double> 	AAprop(int j)	{
     vector< vector<double>  > AA{
                           {1.8,-4.5,-3.5,-3.5,2.5,-3.5,-3.5,-0.4,-3.2,4.5,3.8,-3.9,1.9,2.8,-1.6,-0.8,-0.7,-0.9,-1.3,4.2}, //  Hkyte  "Hkyte" "HSCTakano" "Hwoese" "Hpace" "HFauchere" "HWimley" "HKarplus" "HAbraham" "HBlack"
                            {9.8,7.3,3.6,4.9,3.0,2.4,4.4,0,11.9,17.2,17.0,10.5,11.9,23.0,15.0,2.6,6.9,24.2,17.2,15.3},      //// HSCTakano Side-chain contribution to protein stability  Takano-Yutani, 2001
                            {7.0,9.1,10.0,13.0,4.8,8.6,12.5,7.9,8.4,4.9,4.9,10.1,5.3,5.0,6.6,7.5,6.6,5.2,5.4,5.6},  // Hwoese
                            {0,0.21,0.65,0.69,0.68,0.39,0.40,1.00,0.66,0.41,0.21,0.26,0.24,0.54,3.16,0.50,0.66,0.49,0.53,0.61},  // Hpace
                            {0.42,-1.37,-0.82,-1.05,2.09,-0.30,-0.87,0.0,0.18,2.45,2.31,-1.35,1.67,2.43,0.98,-0.05,0.35,3.06,1.31,1.66},  // HFauchere Pliska 1983
                            {4.08,3.91,3.83,3.02,4.49,3.67,2.23,4.24,4.08,4.52,4.81,3.77,4.48,5.38,3.80,4.12,4.11,6.10,5.19,4.18},   // HWimley White
                            {2.15,2.23,1.05,1.13,1.20,1.65,1.73,1.18,2.45,3.88,4.10,3.05,3.43,3.46,3.10,1.40,2.25,4.11,2.81,3.38},  //HKarplus 1995
                            {0.440,-2.420,-1.320,-0.310,0.580,-0.710,-0.340,0.000,-0.010,2.460,2.460,-2.450,1.100,2.540,1.290,-0.840,-0.410,2.560,1.630,1.730},// HAbraham D.J., Leo A.J. 1987
                            {0.616,0.000,0.236,0.028,0.680,0.251,0.043,0.501,0.165,0.943,0.943,0.283,0.738,1.000,0.711,0.359,0.450,0.878,0.880,0.825},// HBlack Mould 1991
                            /////{300,-540,-170,-60,1320,-170,-150,-100,100,2400,1280,1500,2300,2700,1060,100,250,3200,1900,1300},//
                            {0.10,1.91,0.48,0.78,-1.42,0.95,0.83,0.33,-0.50,-1.13,-1.18,1.40,-1.59,-2.12,0.73,0.52,0.07,-0.51,-0.21,-1.27},// HGuy 1985
                            {0.20,-1.34,-0.69,-0.72,0.67,-0.74,-1.09,0.06,0.04,0.74,0.65,-2.0,0.71,0.67,-0.44,-0.34,-0.26,0.45,-0.22,0.61}, // Hmiller 1987
                            {-12.04,39.23,4.25,23.22,3.95,2.16,16.81,-7.85,6.28,-18.32,-17.79,9.71,-8.86,-21.98,5.82,-1.54,-4.15,-16.19,-1.51,-16.22},// HvHeijne energia libre de trasferencia Von Heijne  1979
                            {0.390,-3.950,-1.910,-3.810,0.250,-1.300,-2.910,0.000,-0.640,1.820,1.820,-2.770,0.960,2.270,0.990,-1.240,-1.000,2.130,1.470,1.300},// HRoseman 1988
                            {-0.48,-0.06,-0.87,-0.75,-0.32,-0.32,-0.71,0.0,-0.51,0.81,1.02,-0.09,0.81,1.03,2.03,0.05,-0.35,0.66,1.24,0.56}, // HLawsonCH 1984 CHP/AGUA
                            {5.330,4.180,3.710,3.590,7.930,3.870,3.650,4.480,5.100,8.830,8.470,2.950,8.950,9.030,3.870,4.090,4.490,7.660,5.890,7.630},// Hmiyazawa 1985
                            {12.28,11.49,11,10.97,14.93,11.28,11.19,12.01,12.84,14.77,14.10,10.80,14.33,13.43,11.19,11.26,11.65,12.95,13.29,15.07},// Hponnuswamy
                            {1.5,3,2,6,1,2,7.5,1,4,1,1,4.5,1,0,3,2,2,2,3,1}, // Hkuntz	 "Hmiyazawa"  "Hponnuswamy"	"Hkuntz" "HBull" "HMeek7" "HMeek2" "HCowan3" "HCowan7" "HTMtPilpel"	"HTMePilpel" "HTMcPilpel" "HTMiPilpel" "HTMtPilpel" "HAdamian" "Cadamian"
                            {-200,-120,80,-200,-450,160,-300,0,-120,2260,2460,-350,1470,2330,-980,-390,-520,2010,2240,1560}, //HBull Breese
                            {0.5,0.8,0.8,-8.2,-6.8,-4.8,-16.9,0,-3.5,13.9,8.8,0.1,4.8,13.2,6.1,1.2,2.7,14.9,6.1,2.7},//HMeek7 ph 7.4 1980
                            {-0.1,-4.5,-1.6,-2.8,-2.2,-2.5,-7.5,-0.5,0.8,11.8,10,-3.2,7.1,13.9,8,-3.7,1.5,18.1,8.2,3.3},// HMeek2 ph 2.1 1980
                            {0.420,-1.560,-1.030,-0.510,0.840,-0.960,-0.370,0.000,-2.280,1.810,1.800,-2.030,1.180,1.740,0.860,-0.640,-0.260,1.460,0.510,1.340},//HCowan3 R., Whittaker R.G. ph 3.4 1990
                            {0.350,-1.500,-0.990,-2.150,0.760,-0.930,-1.950,0.000,-0.650,1.830,1.800,-1.540,1.100,1.690,0.840,-0.630,-0.270,1.350,0.390,1.320},//HCowan7 R., Whittaker R.G. ph 7.5 1990
                            {0.02,0.18,-0.68,-0.87,0.27,-0.54,-0.86,-0.12,-0.34,0.12,0.19,0.25,-0.31,-0.11,-0.51,-0.21,-0.02,-0.12,-0.12,0.23},//HTMtPilpel Hidrofobicidad TM total Pilpel 1999
                            {0.03,-0.53,-0.69,-1.35,-0.53,-0.78,-1.01,-0.02,-0.21,0.34,0.11,-0.85,-0.28,-0.41,-0.16,-0.11,-0.03,0.25,-0.18,0.27},//HTMePilpel Hidrofobicidad TM extracelular Pilpel 1999
                            {0.09,-0.84,-1.11,-1.1,0.12,-0.83,-1.1,-0.05,-0.69,0.12,0.26,-0.56,-0.39,-0.16,-0.66,-0.22,-0.03,-0.65,-0.7,0.31},// HTMcPilpel Hidrofobicidad TM central Pilpel 1999
                            {-0.18,0.44,-0.75,-1.23,0.61,-0.39,-1.33,-0.33,0.21,0.09,0.17,0.66,-0.32,0.12,-0.75,-0.41,-0.43,0.29,0.26,0},//HTMiPilpel  Hidrofobicidad TM intracelular Pilpel 1999
                            {-0.09,0.4,-0.38,-0.73,0.5,-0.33,-0.7,-0.23,-0.15,0.13,0.09,0.46,-0.22,-0.07,-0.37,-0.21,0,0.16,0.23,0.11},// HTMtPilpel Hidrofobicidad TM both termini Pilpel 1999
                            {0.2,-0.23,0.03,-0.11,-0.09,0.3,0.06,0.44,0.12,-0.06,-0.1,-0.52,-0.19,-0.27,0.01,0.22,0.2,-0.53,-0.12,0.01},// HAdamian 2005
                            {0.12,1.03,1,0.73,-0.01,0.79,0.98,0.62,1.26,-0.38,-0.32,0.51,-0.01,-0.4,0.26,0.35,0.25,-0.3,0.04,-0.34},// CAdamian 2005 "MGromiha" "LCGromiha1120" "LCGromiha2130"	"LCGromiha3140" "LCGromiha4150"	"LCromihaM501" "LCromihaM502"
                            {2.11,1.94,1.84,1.8,1.88,2.03,2.09,1.53,1.98,1.77,2.19,1.96,2.27,1.98,1.32,1.57,1.57,1.9,1.67,1.63},// MGromiha   medium contacts globular prot. Gromiha Selvaraj 2004  78proimp
                            {0.56,0.58,0.43,0.38,1.16,0.49,0.40,0.59,0.56,0.92,0.71,0.44,0.71,0.74,0.55,0.58,0.77,0.89,0.81,0.92},// LCGromiha1120 Gromiha long range contacts 11-20    Gromiha Selvaraj 1999 312proimp
                            {0.62,0.46,0.50,0.44,0.99,0.43,0.37,0.61,0.52,0.79,0.61,0.38,0.61,0.63,0.55,0.48,0.49,0.44,0.60,0.82},// LCGromiha2130 Gromiha  long range contacts 21-30    Gromiha Selvaraj 1999 312proimp
                            {0.39,0.39,0.39,0.29,0.78,0.32,0.29,0.47,0.37,0.55,0.43,0.28,0.41,0.45,0.35,0.37,0.35,0.44,0.52,0.51},// LCGromiha3140 long range contacts 31-40    Gromiha Selvaraj 1999 312proimp
                            {0.31,0.22,0.24,0.23,0.51,0.20,0.22,0.34,0.32,0.41,0.31,0.19,0.30,0.30,0.29,0.26,0.29,0.31,0.37,0.44},// LCGromiha4150 long range contacts 41-50    Gromiha Selvaraj 1999 312proimp
                            {1.44,1.14,1.20,0.94,1.58,0.93,0.83,1.51,1.11,1.78,1.59,0.84,1.57,1.54,1.31,1.32,1.46,1.51,1.44,1.84},// LCromihaM50 long range contacts >50  1  Gromiha Selvaraj 1999 312proimp
                            {4.05,3.46,3.38,2.8,6.05,3.05,2.63,4.17,3.52,5.26,4.52,2.77,4.41,4.54,3.71,3.68,4.13,4.6,4.68,5.45},// LCGromihaM50 long range contacts >50 2  Gromiha Selvaraj 1999 312proimp
                            {0.717,0.814,0.851,0.921,0.668,0.849,0.963,0.843,0.754,0.632,0.681,0.912,0.685,0.599,0.850,0.840,0.758,0.626,0.615,0.619},// scFBt parametros de scala  factores b	"scFBt" "scFBRR" "scFBRF" "scFBFF"
                            {0.556,0.676,0.735,0.726,0.607,0.729,0.805,0.651,0.597,0.510,0.504,0.863,0.575,0.465,0.752,0.698,0.648,0.577,0.461,0.503},// scFBRR parametros de scala  factores b RR
                            {0.704,0.807,0.848,0.889,0.671,0.817,0.911,0.811,0.734,0.617,0.650,0.862,0.641,0.582,0.866,0.847,0.742,0.609,0.567,0.603},// scFBRF parametros de scala  factores b RF
                            {0.847,0.942,0.901,1.055,0.670,1.007,1.110,0.967,0.894,0.686,0.788,1.016,0.740,0.653,0.857,0.914,0.862,0.656,0.741,0.707},// scFBFF parametros de scala  factores b FF
                            {-3.93,2.74,0.48,-3.93,-2.23,-0.09,-3.93,-3.93,1.56,-3.93,-3.93,2.16,-2.58,-3.93,-3.93,0.42,-3.46,-1.96,-2.87,-3.93},//GaaDna_pr, g-aa  mandel 1998	 "GaaDna_pr" "AaaDna_pr" "TaaDna_pr" "CaaDna_pr" "IP_Pkeskin" "IP_PkeskinI"	"IP_PkeskinII" "IP_PkeskinIII"
                            {-3.93,0.34,1.93,-3.37,0.07,1.16,-1.24,-3.93,0.46,-3.93,-3.93,-0.08,-0.28,-3.93,-3.93,-0.68,-0.06,-3.93,-2.87,-3.93},//AaaDna_pr, a-aa  mandel 1998
                            {0.66,1.25,0.71,-3.93,-2.23,0.31,-3.93,-3.93,0.87,0.65,-0.94,0.21,0.42,-0.81,-0.30,-0.28,-0.06,-1.96,0.54,-0.17},//TaaDna_pr,  t-aa mandel 1998
                            {-3.72,-3.93,0.71,1.01,0.07,-3.09,0.55,-3.93,-0.23,-3.44,-3.93,-3.93,-0.28,-0.12,-3.29,-0.68,-1.16,-3.93,0.13,-3.57},//CaaDna_pr,  C-aa mandel 1998
                            {0.9,1.14,1.109,0.826,1.427,0.954,0.866,0.671,1.076,1.068,1.127,0.732,1.083,1.213,0.735,0.955,0.942,1.075,1.318,0.986},// IP_Pkeskin protein-protein interfaces 2004
                            {0.966,1.137,0.908,0.966,1.78,0.91,0.725,0.699,1.331,1.044,1.073,0.611,1.289,1.409,0.987,1.027,0.878,0.969,1.56,0.999},// IP_PkeskinI protein-protein interfaces I 2004
                            {0.761,0.754,1.215,0.602,1.338,1.1014,0.988,0.599,0.779,1.469,1.532,0.578,0.839,1.373,0.343,1.029,0.741,0.856,1.379,1.112},// IP_PkeskinII protein-protein interfaces II 2004
                            {0.897,1.276,1.261,0.775,1.153,0.951,0.982,0.644,0.9,0.999,1.039,0.887,0.982,0.956,0.626,0.864,1.080,1.26,1.102,0.865},//IP_PkeskinIII keskin protein-protein interfaces III 2004
                            {-0.5,3,0.2,3,-1,0.2,3,0,-0.5,-1.8,-1.8,3,-1.3,-2.5,0,0.3,-0.4,-3.4,-2.3,-1.5},// HHopp 1981 "HHopp" "HParker" "HWilson" "HCornette" "HEisenberg" "Hlevitt" "HZIMMERMAN" "HWolfenden"
                            {2.1,4.2,7.0,10.0,1.4,6.0,7.8,5.7,2.1,-8.0,-9.2,5.7,-4.2,-9.2,2.1,6.5,5.2,-10.0,-1.9,-3.7},// HParker 1986
                            {-0.300,-1.100,-0.200,-1.400,6.300,-0.200,0.000,1.200,-1.300,4.300,6.600,-3.600,2.500,7.500,2.200,-0.600,-2.200,7.900,7.100,5.900},// HWilson K.J., Honegger A., Stotzel R.P., Hughes G.J 1981
                            {-0.96,0.75,-1.94,-5.68,4.54,-5.30,-3.86,-1.28,-0.62,5.54,6.81,-5.62,4.76,5.06,-4.47,-1.92,-3.99,0.21,3.34,5.39},  //HCornette 1987
                            {0.25,-1.76,-0.64,-0.72,0.04,-0.69,-0.62,0.16,-0.40,0.73,0.53,-1.10,0.26,0.61,-0.07,-0.26,-0.18,0.37,0.02,0.54}, // HEisenberg (1982)
                            {-0.5,3.0,0.2,2.5,-1,0.2,2.5,0,-0.5,-1.8,-1.8,3,-1.3,-2.5,-1.4,0.3,-0.4,-3.4,-2.3,-1.5}, // Hlevitt 1976
                            {0.83,0.83,0.09,0.64,1.48,0,0.65,0.1,1.10,2.52,3.07,1.60,1.40,2.75,2.70,0.14,0.54,0.31,2.97,1.79},//HZIMMERMAN
                            {1.94,-19.92,-9.68,-10.95,-1.24,-9.38,-10.20,2.39,-10.27,2.15,2.28,-9.52,-1.48,-0.76,-3.68,-5.06,-4.88,-5.88,-6.11,1.99},  // HWolfenden, R.  Hydration potential 1981
                            {-2.61,-58.60,-32.86,-33.38,-10.23,-33.08,-34.09,-5.11,-28.38,1.25,0.80,-25.68,-4.70,-8.79,-2.91,-26.63,-19.35,-22.66,-38.36,-0.25}, // HOoi,Oobatake 1987
                            {-0.7,1.0,1.9,1.0,0.4,-0.6,5.6,-2.5,0.6,-3.5,-2.7,6.8,-0.9,1.4,-3.6,-0.8,-3.8,3.2,-0.6,-2.2},// HSamatey Hidrofobicidad Samatey  1995  "HOoi" "HSamatey" "HBeuming" "MLFBt" "MLFB2R" "MLFBRF" "MLFBFF"
                            {0.4,0.47,0.36,0.23,0.47,0.33,0.28,0.12,0,0.84,0.92,0.88,0.6,0.76,0.39,0.17,0.43,1,0.68,0.82},// HBeuming Hidrofobicidad structure SP Beuming 2004
                           // {-0.48,-0.06,-0.87,-0.75,-0.32,-0.32,-0.71,0.0,-0.51,0.81,1.02,-0.09,0.81,1.03,2.03,0.05,-0.35,0.66,1.24,0.56}, // HLawson 1984 CHP/AGUA
                            {-0.605,-0.448,-0.381,-0.279,-0.692,-0.368,-0.160,-0.537,-0.662,-0.682,-0.631,-0.043,-0.626,-0.719,-0.271,-0.424,-0.525,-0.727,-0.721,-0.669},  //MLFBt Mean Location parameters Gumbel ?Smith 2003
                            {-0.792,-0.638,-0.577,-0.584,-0.823,-0.563,-0.480,-0.760,-0.870,-0.889,-0.865,-0.243,-0.826,-0.934,-0.396,-0.641,-0.707,-0.886,-0.904,-0.847},  //MLFB2R Mean Location parameters 2R Gumbel ?Smith 2003
                            {-0.609,-0.412,-0.398,-0.285,-0.718,-0.362,-0.168,-0.588,-0.634,-0.693,-0.642,-0.065,-0.646,-0.737,-0.315,-0.412,-0.506,-0.754,-0.761,-0.673},  //MLFBRF  Mean Location parameters RF Gumbel ?Smith 2003
                            {-0.387,-0.314,-0.185,0.072,-0.593,-0.158,0.189,-0.241,-0.485,-0.538,-0.422,0.176,-0.485,-0.522,-0.120,-0.231,-0.390,-0.614,-0.503,-0.501},     //MLFBFF Mean Location parameters FF Gumbel ?Smith 2003
                            {0,0.06,0.6,0.59,0.6,0.32,0.34,1.1,0.62,0.35,0.19,0.15,0.21,0.47,2.72,0.52,0.57,0.47,0.47,0.46},// A Helix peptide Agadir 1995	AHAgadir   "AHAgadir" "HLawsonET" "BA_32" "AHBlaber" "AHChou" "AHChakra" "AHOneil" "AHDel"
                            {0.87,0.85,0.09,0.66,1.52,0.0,0.67,0.10,0.87,3.15,2.17,1.64,1.67,2.87,2.77,0.07,0.07,3.77,2.67,1.87}, //Lawson 1984 etOH/Agua "AHAgadir"  "HLawsonET" "BA_32" "AHBlaber" "AHChou" "AHChakra" "AHOneil" "AHDel" "AHMPRao" "AHpLiuH2O"
                            {0,0.14,0.66,0.71,1,0.48,0.55,0.91,0.78,0.81,0.35,0.19,0.31,0.69,4.08,0.41,0.79,0.98,0.82,0.88},  //BA-32	 BA_32
                            {0,0.19,0.57,0.54,0.54,0.16,0.43,0.96,0.39,0.12,0.04,0.23,0.1,0.37,3.46,0.43,0.42,0.38,0.24,0.33},// A Helix Blaber 1994  AHBlaber
                            {1.42,0.98,0.67,1.01,0.70,1.11,1.51,0.57,1.00,1.08,1.21,1.16,1.45,1.13,0.57,0.77,0.83,1.08,0.69,1.06},// A Helix Chou Fasman 1978 AHChou
                            {-0.258,-0.047,0.635,0.635,0.570,0.314,0.433,1.62,0.525,0.445,0.022,0.108,0.251,0.672,4,0.525,1.07,0.653,0.511,0.797},// A Helix peptide Chakrabartty 1990	AHChakra
                            {-0.77,-0.68,-0.07,-0.15,-0.23,-0.33,-0.27,0,-0.06,-0.23,-0.62,-0.65,-0.5,-0.41,3,-0.35,-0.11,-0.45,-0.17,-0.14},// A Helix ONeil DeGrado 1990	AHOneil
                            {1.489,1.224,0.772,0.924,0.966,1.164,1.504,0.51,1.003,1.003,1.236,1.172,1.363,1.195,0.492,0.739,0.785,1.09,0.787,0.990},// Helice a Deleage Roux 1987 AHDel
                            {1.360,0.150,0.330,0.110,1.270,0.330,0.250,1.090,0.680,1.440,1.470,0.090,1.420,1.570,0.540,0.970,1.080,1.000,0.830,1.370},// Rao M.J.K., Argos P 1986 AHMPRao
                            {3.72,2.84,2.82,2.66,3.28,2.87,2.54,3.44,2.9,3.88,3.84,2.65,3.67,3.78,1.70,2.99,3.27,3.20,3.33,3.82},// A Helix peptide en n-octanol Liu Deber 1998   AHpLiuOct
                            {2.06,1.35,0.59,1.47,0.91,1.19,2.25,0.34,0.79,1.54,1.80,1.12,1.33,1.12,0,0.79,0.9,1.13,1.16,1.03},// A Helix peptide en agua Liu Deber 1998	 AHpLiuH2O
                            {1.86,0.08,0.1,0.06,0.76,0.11,0.05,1.36,0.27,2.69,2.53,0.08,1.26,1.63,0.30,0.63,0.7,0.45,0.39,2.32},// A Helix peptide stats Liu Deber 1998  "AHpLiuStats" "AHRohl" "AHpYang" "AHpW" "AHRichardson" "AHWilliams" "AHTmGromiha" 	"AHpWarren" "AHMunoz" "AHPPMunoz"  "AHPhiplusBahar" "AHPhiMinusBahar" "AHPhiPlusMinBahar" "AHEaPhiPlusMinBahar"
                            {-0.27,-0.052,0.69,0.54,0.64,0.28,0.35,1.7,0.84,0.44,0.095,0.019,0.25,0.73,3.8,0.52,0.95,0.69,0.42,0.77},// A Helix peptide Rohl agua 1996	AHRohl
                            {-0.36,0.05,0.49,0.42,0.81,0.44,-0.15,0.69,0.175,0.02,0,0.225,0.165,0.75,2.2,0.4,0.525,0.75,1.585,0.37},// A Helix peptide yang 1997  AHpYang
                            {-0.04,-0.02,0.15,0.23,0.01,0.01,0.02,0.31,0.22,-0.08,0.08,0.04,-0.11,-0.05,1.01,0.16,0.12,-0.06,-0.01,0.03},// A Helix peptide Wojcik 1990 AHpW
                            {1.8,1.3,0.9,1,0.7,1.3,0.8,0.5,1,1.2,1.2,1.1,1.5,1.3,0.3,0.6,1,1.5,0.8,1.2},//H A Richardson  warren  AHRichardson
                            {1.41,1.21,0.76,0.99,0.66,1.27,1.59,0.43,1.05,1.09,1.34,1.23,1.3,1.16,0.34,0.57,0.76,1.02,0.74,0.98},//H A Williams  warren	 AHWilliams
                            {1.334,0.169,0.494,0.171,1.062,0.343,0.168,0.998,0.530,1.803,1.623,0.115,1.413,1.707,0.563,0.817,0.838,1.274,1.146,1.599},// TM proteins gromiha 1999 AHTmGromiha
                            {1.39,1.14,0.8,0.92,0.47,1.37,1.34,0.63,0.87,0.96,1.27,1.14,1.2,1.07,0.41,0.76,0.85,1.05,0.91,0.75},// a helice propensity thermophilic warren 1995	AHpWarren
                            {0.423,0.503,0.906,0.870,0.877,0.594,0.167,1.162,0.802,0.566,0.494,0.615,0.444,0.706,1.945,0.928,0.884,0.690,0.778,0.706},// HaMunoz-Serrano, 1994	AHMunoz
                            {0.617,0.753,1.089,0.932,1.107,0.770,0.675,1.361,1.034,0.876,0.740,0.784,0.736,0.968,1.780,0.969,1.053,0.910,1.009,0.939},// PP HaMunoz-Serrano, psifi 1994	 AHPPMunoz
                            {-1.65,-1.46,-1.07,-1.20,-1.01,-1.46,-1.6,-0.54,-1.1,-1.29,-1.47,-1.38,-1.53,-1.24,-1.17,-1.17,-1.01,-1.4,-1.07,-1.17},// Ha phi + Bahar 1997 1.pdf proteins AHPhiplusBahar
                            {-1.69,-1.49,-1.21,-1.35,-0.99,-1.53,-1.6,-0.68,-1.27,-1.15,-1.48,-1.42,-1.58,-1.32,0.11,-1.16,-1.05,-1.33,-1.09,-0.98},// Ha phi - Bahar 1997 1.pdf proteins AHPhiMinusBahar
                            {-0.52,-0.56,-0.74,-0.66,-1.01,-0.55,-0.43,-1.23,-0.72,-0.83,-0.61,-0.48,-0.57,-0.72,-0.86,-0.75,-0.87,-0.64,-0.83,-0.92},// Ha phi+ - Bahar 1997 1.pdf proteins AHPhiPlusMinBahar
                            {-3.85,-3.5,-3.08,-3.21,-3.01,-3.54,-3.64,-2.45,-3.09,-3.26,-3.55,-3.38,-3.68,-3.28,-1.92,-3.08,-2.91,-3.37,-3,-3.06},// Ha Ea phi+ - Bahar 1997 1.pdf proteins  AHEaPhiPlusMinBahar
                            {0.83,0.93,0.89,0.54,1.19,1.10,0.37,0.75,0.87,1.60,1.30,0.74,1.05,1.38,0.55,0.75,1.19,1.37,1.47,1.70},// Beta sheet Chou Fasman 1978  "BsChou" "BDeleage" "BmKim1" "BmKim2" "BesEswar" "BisEswar" "BexEswar" "Bmunoz" "BppMunoz"
                             //{1.420,0.980,0.670,1.010,0.700,1.110,1.510,0.570,1.000,1.080,1.210,1.160,1.450,1.130,0.570,0.770,0.830,1.080,0.690,1.060}    "ARNDCQEGHILKMFPSTWYV";
                           {0.709,0.920,0.604,0.541,1.191,0.840,0.567,0.657,0.863,1.799,1.261,0.721,1.210,1.393,0.354,0.928,1.221,1.306,1.266,1.965},// Beta Deleage Roux 1987	BDeleage
                            {0,0.45,-0.08,-0.94,0.52,0.23,0.01,-1.2,-0.02,1,0.51,0.27,0.72,0.86,-3,0.7,1.1,0.54,0.96,0.82},// Beta Minor Kim  feb 1994			  BmKim1
                            {0,-0.43,-0.24,-0.10,0.08,0.04,0.31,-0.85,-0.01,0.02,-0.24,-0.40,-0.02,0.16,-4,0.63,0.83,-0.17,0.11,0.17},// Beta Minor Kim  sep 1994	BmKim2
                            {0.73,1.00,0.70,0.61,1.12,0.81,0.81,0.45,1.02,1.42,1.12,0.85,1.07,1.29,1,0.94,1.36,1.47,1.41,1.68},// Beta edge strand Eswar 2003	BesEswar
                            {0.84,0.81,0.57,0.66,1.17,0.84,0.63,0.54,0.96,1.67,1.20,0.83,1.27,1.29,0.76,0.93,1.19,1.22,1.46,1.72},// Beta inner strand Eswar 2003  BisEswar
                            {0.74,1.02,0.76,0.78,1.27,0.94,0.92,0.37,0.92,1.18,1.01,0.98,1.01,1.18,2.27,0.96,1.23,0.90,1.01,1.25},// Beta Extended  Eswar 2003	  BexEswar
                            {1.080,0.976,1.197,1.266,0.733,1.050,1.085,1.104,0.906,0.583,0.789,1.026,0.812,0.685,1.412,0.987,0.784,0.755,0.665,0.546},// Beta Munoz-Serrano, 1994 Bmunoz
                            {0.978,0.784,0.915,1.038,0.573,0.863,0.962,1.405,0.724,0.502,0.766,0.841,0.729,0.585,2.613,0.784,0.569,0.671,0.560,0.444},// pp Beta psifi Munoz-Serrano, 1994	 BppMunoz
                            {1.40,1.23,1.61,1.89,1.14,1.33,1.42,2.06,1.25,1.02,1.33,1.34,1.12,1.07,3.90,1.20,0.99,1.10,0.98,0.87},// BetaR res Munoz-Serrano, 1994	  BRMunoz "BRMunoz" "BLinding" "CoilDeleage" "TurnDeleage" "TmAHMonne" "TmAHMonne2"
                            {-0.26154,-0.17659,0.22989,0.22763,-0.015152,-0.187677,-0.20469,0.43323,-0.0012174,-0.42224,-0.33793,-0.100092,-0.22590,-0.22557,0.55232,0.14288,0.0088780,-0.243375,-0.20751,-0.38618},// Linding 2003  3701  BLinding
                            {0.824,0.893,1.167,1.197,0.953,0.947,0.761,1.251,1.068,0.886,0.810,0.897,0.810,0.797,1.540,1.130,1.148,0.941,1.109,0.772},// coil Deleage Roux 1987		CoilDeleage
                            {0.788,0.912,1.572,1.197,0.965,0.997,1.149,1.860,0.970,0.240,0.670,1.302,0.436,0.624,1.415,1.316,0.739,0.546,0.795,0.387},// turn Deleage Roux 1987		TurnDeleage
                            {0.5,1.7,1.7,1.6,0.6,1.6,1.6,1.3,1.6,0.6,0.4,1.6,0.5,0.4,1.7,0.7,0.4,0.7,0.6,0.5},// turn membrane helix propensity Monné JMB 1999	 TmAHMonne
                            {0.4,1.5,1.6,1.5,0.7,1.4,1.3,1.1,1.4,0.5,0.3,1.4,0.5,0.3,1.6,0.9,0.7,0.9,0.9,0.4},// turn membrane helix propensity Monné JMB 1999 2  TmAHMonne2
                            {0.83,0.79,0.41,0.33,0.95,0.89,0.64,0.59,0.84,1.72,1.27,0.8,1.36,1.57,0.29,0.86,1.2,1.47,1.53,2.05},// BT Strand preference Daffner Argos 1994	BspDaffner  "BspDaffner" "BsiDaffner" "BsnDaffner" "TurnpDaffner" "LoopEswar"
                            {0.67,1.23,1.14,2,1,1,1,2.38,1.28,0.55,0.59,1.05,1,0.76,3.67,0.92,1.44,0.67,0.72,0.77},// BT Strand inside preference Daffner Argos 1994		BsiDaffner
                            {0.87,1.38,0.5,1.5,1.5,0.62,0.8,1.14,1.75,1,0.91,1.05,0.82,0.96,1.17,1.31,0.84,1.44,1.04,0.92},// BT Strand next to preference Daffner Argos 1994  BsnDaffner
                           {0.84,0.91,1.48,1.28,0.69,1,0.78,1.76,0.53,0.55,0.49,0.95,0.52,0.88,1.47,1.29,1.05,0.88,1.28,0.51},// Turn preference Daffner Argos 1994	   TpDaffner
                           {0.74,0.94,1.36,1.32,1.13,0.85,0.86,1.73,1.04,0.67,0.63,0.93,0.57,0.81,1.36,1.24,1.07,0.80,0.88,0.66},// Loop  Eswar 2003			LoopEswar
                           {0.06,-0.32,0.02,-0.07,-0.06,-0.45,-0.05,0.53,0.22,-0.29,-0.21,0.13,0.12,-0.98,0.43,-0.12,0.31,-0.28,-0.73,-0.24},//wrabl hilser low thermodynamic scale 2001   LtWrabl "LtWrabl" "MtWrabl" "GtWrabl"
                           {0.01,-0.19,0.08,0.25,-0.22,-0.09,0.15,-0.03,-0.09,0.08,0.01,0.1,-0.04,-0.4,-0.04,-0.08,-0.02,-1.18,-0.08,0.11},//wrabl hilser medium thermodynamic scale 2001  MtWrabl
                           {-0.07,0.37,-0.1,-0.25,0.23,0.37,-0.12,-1.17,-0.18,0.15,0.17,-0.27,-0.09,0.66,-0.7,0.17,-0.42,0.65,0.46,0.09},//wrabl hilser high thermodynamic scale 2001	    GtWrabl
                           {0.669,1.04,2.35,2.06,0.945,1.07,0.787,0.621,1.55,0.511,0.885,0.977,0.965,1.2,0.474,1.01,1.32,1.02,1.24,0.421},// Ha1 Shortle 2002 Ha1Shortle  "LtWrabl" "MtWrabl" "GtWrabl" "Ha1Shortle"  "Ha2Shortle" "Ha3Shortle" "L1Shortle" "L2Shortle"  "L3Shortle"
                           {1.11,1.03,1.09,1.16,0.878,1.15,1.25,0.482,1.05,0.561,1.07,1.13,1.06,0.772,1.38,1.43,1.09,1.04,0.781,0.572},// Ha2 Shortle 2002	  Ha2Shortle
                           {1.44,1.2,0.655,0.881,0.704,1.26,1.37,0.432,0.832,1.09,1.29,1.17,1.23,0.952,0.621,0.745,0.761,1.08,0.915,0.955},// Ha3 Shortle 2002 Ha3Shortle
                           {0.931,0.88,1.07,0.822,1.79,0.64,0.517,1.75,1.28,0.487,0.299,0.62,0.931,1.65,0.013,2.05,1.84,1.44,1.75,0.74},// L1 Shortle 2002     L1Shortle
                           {1.08,1.24,0.65,0.469,1.28,1.06,0.852,0.388,1.34,1.21,0.673,0.984,1.31,1.42,0.007,1.44,1.42,1.05,1.53,1.57},// L2 Shortle 2002     L2Shortle
                           {0.828,0.968,2.68,2.07,1.89,1.12,0.828,0.284,2.06,0.747,0.619,0.949,1.23,1.14,0.006,1.13,0.751,1.16,1.12,0.789},// L3 Shortle 2002  L3Shortle
                           {0.518,0.85,1.01,0.926,1.52,0.876,0.598,1.52,0.862,0.747,0.781,0.898,0.694,1.06,0.278,1.47,2.51,0.941,1.23,0.895},// m1 Shortle 2002	  m1Shortle "m1Shortle" "m2Shortle" "m3Shortle" "r1Shortle" "r2Shortle" "r3Shortle" "LHShortle"
                           {0.538,0.882,0.667,0.537,1.23,0.84,0.727,0.206,0.926,2.19,1.29,0.856,0.982,1.44,0.071,0.731,1.32,1.21,1.41,2.28},// m2 Shortle 2002	 m2Shortle
                           {0.529,0.921,2,2.01,1.46,0.847,0.697,0.184,1.64,1.59,1.03,0.857,0.96,1.5,0.0997,0.635,0.746,1.09,1.24,1.37},// m3 Shortle 2002	    m3Shortle
                           {0.762,0.578,0.899,1.25,1.13,0.609,0.484,1.73,0.629,0.237,0.521,0.61,0.524,0.465,4.38,2,1.48,0.542,0.521,0.275},// r1 Shortle 2002	r1Shortle
                           {0.996,0.846,0.791,1.05,1.21,0.783,0.866,0.297,0.849,0.873,0.945,0.949,0.842,0.829,3.82,1.02,0.792,1.02,0.848,0.887},// r2 Shortle 2002 r2Shortle
                           {0.744,0.812,2.13,2.51,1.38,0.778,0.881,0.395,1.05,0.795,0.879,0.904,0.815,0.841,2.4,0.725,0.442,0.828,0.669,0.73},// r3 Shortle 2002   r3Shortle
                           {0.307,0.582,2.53,1.17,0.5,0.655,0.473,7.07,1.04,0.049,0.198,0.785,0.33,0.465,0.006,0.475,0.152,0.277,0.402,0.065},// LH Shortle 2002   LHShortle
                           {0.288,0.256,0.266,0.321,0.226,0.176,0.204,10.9,0.312,0.0825,0.124,0.269,0.234,0.143,0.0205,0.436,0.192,0.15,0.147,0.0935},// phi Shortle 2002	phiShortle "phiShortle" "othShortle" "ratec-tReimer" "ratet-cReimer" "TurnSantai" "TurnSantai1" "TurnSantai2" "TurnSantai3"
                           {0.746,0.996,1.48,1.19,1.12,0.9,0.977,1.91,1.32,0.659,0.712,1.07,0.668,0.841,0.697,1.11,1.11,0.639,0.976,0.689},// other Shortle 2002 othShortle
                           {32,19,23,28,21,9,15,29,8,6,16,16,12,14,6,28,27,5,12,8},// rate constants c-t  Reimer 1998 45proimp his ph 8.0	 ratec-tReimer
                           {2.7,1.5,3.0,2.2,2,1.2,1.5,4.6,0.8,0.8,2.2,1.2,1.3,4.2,0.4,3.2,2.8,3.0,3.8,0.9},// rate constants t-c  Reimer 1998 45proimp his ph 8.0	 ratet-cReimer
                           {1.05,0.86,0.87,0.72,1.67,0.72,0.60,1.46,0.79,0.78,0.72,0.67,0.73,0.88,2.94,1.23,0.91,0.67,0.97,0.90},// i turn preferences Santa 2002	  TurnSantai
                           {1.04,1.06,0.89,1.86,0.44,1.08,1.17,0.31,0.75,0.69,0.79,1.49,0.55,0.55,2.23,0.98,1.16,0.53,0.66,0.76},// i+1 turn preferences Santa 2002	  TurnSantai1
                           {0.39,1.06,1.76,1.81,1.06,1.13,0.88,0.12,1.29,1.44,0.87,0.87,1.00,1.33,0.02,1.06,0.93,0.67,1.49,1.54},// i+2 turn preferences Santa 2002	  TurnSantai2
                           {0.77,0.63,0.57,0.84,0.67,0.87,0.58,0.65,0.54,1.16,0.75,1.11,0.64,0.63,4.30,0.90,1.07,0.87,0.71,1.35},// i+3 turn preferences Santa 2002	  TurnSantai3
                           {1.77,0.41,0.52,1.09,1.63,0.58,0.67,0.20,0.40,0.39,0.93,0.62,0.89,1.04,6.24,0.68,0.60,0.56,0.67,0.47},// Helice Polyproline  Eswar 2003	  PolypEswar "PolypEswar"  "polypIIStapley" "CubellisL3" "CubellisLp3" "CubellisLm3" "CubellisLm2" "CubellisLm1"
                           {0.84,1.09,0.71,0.88,0.82,1.24,0.81,0.27,0.55,0.69,0.91,1.06,0.86,0.70,5.06,0.7,1,0.61,0.57,0.82},// polyprolyne II Stapley 1999			  polypIIStapley  "Cubellisn1" "Cubellisn2" "Cubellismid" "Cubellisc2" "Cubellisc1"  "Cubellisc+1"  "Cubellisc+2"  "Cubellisc+3"
                           {0.86,0.86,0.73,0.71,0.9,0.91,0.9,0.45,0.76,1.17,0.95,1.05,1.06,0.73,3.84,0.92,1,0.72,0.44,0.8},// L=3 Cubellis 2005                       CubellisL3
                           {0.9,1.11,0.65,0.76,1.05,1.04,0.98,0.35,0.68,1.08,0.98,1.04,0.81,0.66,4.15,0.95,0.98,0.59,0.4,0.85},// L>3 Cubellis 2005                   CubellisLp3
                           {0.85,1.07,1.13,1.12,0.93,1.06,0.88,0.95,0.96,1.24,0.99,0.94,0.68,0.91,1.45,1.29,0.86,0.91,0.73,0.95},//n-3 cubellis 2005                  CubellisLm3
                           {0.86,0.94,1.47,1.07,1.02,1.01,0.87,1.26,1.10,1.15,0.71,0.94,0.78,0.87,1.49,1.18,1.20,1.09,0.56,0.77},//n-2 cubellis 2005                  CubellisLm2
                           {0.67,1.07,1.24,1.02,0.60,0.91,0.84,2.61,0.97,0.92,0.69,0.95,0.93,0.90,0.78,1.05,1.03,0.80,0.63,0.58},//n-1 cubellis 2005                  CubellisLm1
                           {0.93,1.25,0.72,0.54,1.21,1.16,1.00,0.41,0.88,1.49,1.13,1.14,1.12,1.07,2.05,0.84,1.21,0.84,0.53,1.09},//n1 cubellis 2005                   Cubellisn1
                           {0.92,1.20,0.34,0.63,1.16,1.06,1.28,0.36,0.54,1.14,0.81,1.35,0.49,0.59,4.85,0.68,0.80,0.51,0.48,0.88},//n2 cubellis 2005                   Cubellisn2
                           {0.85,0.87,0.54,0.57,1.02,1.10,0.88,0.32,0.66,1.21,1.03,1.26,1.30,0.56,5.16,0.66,0.86,0.77,0.34,0.85},//mid cubellis 2005                  Cubellismid
                           {1.01,0.91,0.73,0.59,0.93,0.92,0.93,0.29,0.74,1.22,1.32,0.93,1.32,0.60,4.71,0.57,0.68,0.96,0.32,0.96},//c2 cubellis 2005                   Cubellisc2
                           {0.87,0.80,0.98,1.14,0.71,0.74,0.89,0.57,0.75,0.78,0.74,0.80,0.59,0.59,4.68,1.55,1.15,0.41,0.44,0.53},//c1 cubellis 2005                   Cubellisc1
                           {0.79,0.78,1.06,1.03,0.72,0.77,0.98,2.52,0.62,0.83,0.69,0.80,0.60,0.72,2.33,1.05,0.93,0.65,0.40,0.56},//c+1 cubellis 2005                  Cubellisc+1
                           {0.92,0.74,1.16,1.40,0.85,1.16,1.48,0.93,1.13,0.98,0.63,0.78,0.64,0.93,0.72,1.55,1.36,1.18,0.71,0.77},//c+2 cubellis 2005                  Cubellisc+2
                           {0.82,0.90,0.97,1.22,0.85,1.14,1.35,1.04,0.96,1.27,0.78,0.77,0.68,0.96,1.10,1.08,1.25,1.37,0.72,0.95},//c+3 cubellis 2005                  Cubellisc+3
                           {0.94,0.47,1.60,1.40,0.91,1.11,1.17,1.45,1.29,1.19,1.10,0.48,1.13,1.54,1.46,1.15,1.09,1.70,1.55,1.22},//Persson Argos 1996      HPersson  "HPersson" "EntDoig" "EntAbsDoig" "EntAbagyan" "EntKoehl" "EntPickett" "EntLee" "EntHaAvbelj" "EntAvbelj" "DeltaEntAvbelj" "HaEntAvbelj" "SEAvbelj"
                           {0,-1.88,-1.03,-0.78,-0.85,-1.73,-1.46,0,-0.95,-0.76,-0.71,-1.89,-1.46,-0.62,-0.06,-1.11,-1.08,-0.99,-1.13,-0.43},   //Doig, Sternberg  1995 EntDoig
                           {0,-13.7,-5.6,-4.2,-3.9,-7.8,-6.4,0,-4.5,-2.6,-2.7,-9.3,-4.4,-3.1,0,-3.9,-3.9,-4.4,-5.6,-1.3},  //entropia absoluta Doig                    EntAbsDoig
                           {0,-2.13,-0.81,-0.61,-1.14,-2.02,-1.65,0,-0.99,-0.75,-0.75,-2.21,-1.53,-0.58,0,-1.19,-1.12,-0.97,-0.99,-0.5},// entropy abagyan totrov 1994  EntAbagyan
                           {0,-1.21,-0.75,-0.65,-0.63,-1.29,-1.31,0,-0.92,-0.94,-0.94,-1.63,-1.24,-0.65,-0.3,-0.43,-0.57,-1.14,-1.07,-0.62},// entropy Koehl Delarue 1994  EntKoehl
                         {0,-2.03,-1.57,-1.25,-0.55,-2.11,-1.81,0,-0.96,-0.89,-0.78,-1.94,-1.61,-0.58,0,-1.71,-1.63,-0.97,-0.98,-0.51},// entropy Pickett Sternberg 1993         EntPickett
                         {0,-2.13,-0.99,-0.6,-1.06,-1.51,-1.06,0,-1,-0.52,-0.49,-1.76,-1.37,-0.42,0,-1.1,-0.99,-0.82,-0.83,-0.04},// entropy Lee 1994                            EntLee
                         {0,-1.991,-1.436,-0.959,-0.535,-1.929,-1.547,0,-0.794,-0.481,-0.696,-1.849,-1.452,-0.409,0,-1.686,-1.363,-0.633,-0.858,-0.172},// entropy Ha Avbelj 1998 EntHaAvbelj
                         {0,-2.12,-1.708,-1.318,-0.572,-2.107,-1.763,0,-0.895,-0.926,-0.763,-1.973,-1.54,-0.544,0,-1.695,-1.618,-0.909,-1.019,-0.541},// entropy  Avbelj 1998     EntAvbelj
                         {0,0.13,0.27,0.36,0.04,0.18,0.22,0,0.1,0.45,0.07,0.12,0.09,0.13,0,0.01,0.25,0.28,0.16,0.37},// Delta T entropy  Avbelj 1998                             DeltaEntAvbelj
                         {0.54,0.58,0.67,0.56,0.54,0.59,0.55,0.74,0.62,0.45,0.53,0.60,0.57,0.56,0.41,0.54,0.47,0.52,0.55,0.46},//  main Ha entropy  Avbelj 1998                  HaEntAvbelj
                         {0.285,0.275,0.116,0.102,0.183,0.257,0.345,-0.103,0.176,0.404,0.367,0.244,0.358,0.331,-0.325,0.078,0.126,0.232,0.241,0.325},//  screening electrostatic   SEAvbelj 1998
                         {0.43,0.48,0.39,0.31,0.53,0.55,0.45,0.59,0.57,0.41,0.49,0.51,0.46,0.51,0.28,0.39,0.35,0.43,0.49,0.43},//  main B entropy  Avbelj 1998  "mBEntAvbelj" "EntCoilJha" "EntZamanUA" "EntZamanAA" "EntZamanGS" "EntZamanAmb" "AsaHelFan" "AsaSFan"
                         {0,0.174,0.283,0.275,-0.127,0.104,0.109,0.475,0.110,-0.138,0.007,0.156,-0.143,0.049,-0.675,0.174,0.134,-0.153,0.050,-0.046},// Backbone entropy Coil library Jha 2005 EntCoilJha
                         {0,-0.0635,-0.0515,-0.12495,-0.049,-0.0383,-0.1135,0.0455,-0.0955,-0.1565,-0.04995,-0.0455,-0.0515,-0.0645,-0.2505,-0.0235,-0.1005,-0.0095,-0.02,-0.1495},//Entropy Zaman 2003 OPLS-UA EntZamanUA
                         {0,0.011,-0.008,-0.852,-0.017,-0.01,-0.139,0.22,0.021,-0.674,-0.11,-0.111,-0.033,0.22,-0.473,-0.036,-0.164,-0.32,-0.065,-0.276},//Entropy Zaman 2003 OPLS-AA-01      EntZamanAA
                         {0,0.006,-0.053,-0.1485,0.0245,-0.0385,0.0004,0.045,0.0105,-0.181,-0.032,0.033,-0.09995,-0.2061,-0.398,0.025,-0.022,-0.0945,-0.125,0.067},//Entropy Zaman 2003 G-S-94  EntZamanGS
                         {0,-0.0415,0.0045,-0.251,0.06,-0.02,-0.041,0.498,0.1795,-0.0425,0.09,-0.0815,0.0255,0.3465,-0.2675,-0.024,-0.024,0.1015,-0.06765,-0.076},//Entropy Zaman 2003 AMBER 94  EntZamanAmb
                         {21.71,97.24,64.39,71.28,11.29,80.67,87.78,17.72,64.11,18.56,22.94,97.21,29.01,23.80,50.73,40.17,38.90,31.87,44.69,19.00},//Asa Helice Fan Jiang  2003                  AsaHelFan
                         {7.21,69.85,32.72,37.16,7.31,53.71,56.33,3.34,34.30,12.39,14.18,73.19,19.19,17.68,31.91,19.06,32.69,33.77,33.60,12.00},// Asa Strand Fan Jiang  2003                    AsaSFan
                         {39.07,122.22,79.66,78.03,21.90,97.47,116.04,36.27,69.09,37.03,43.27,119.09,58.46,43.35,63.52,55.33,64.53,78.55,70.70,38.38},// Asa Coil Fan Jiang  2003               AsaCoilFan "AsaCoilFan" "NcontSamanta" "AsaSamanta" "AsaWinLins" "AsaFtotLins" "AsaFcoilLins" "AsaFbetaLins" "AsaWinphoLins" "AsaWinphiLins"
                         {22.5,30.4,24.4,23.3,27.0,26.7,24.6,18.7,28.4,29.2,29.3,23.5,31.1,33.4,19.5,22.2,23.6,37.9,34.0,26.5},//numero de contactos medio samanta                              NcontSamanta
                         {24.1,30.6,30.0,32.2,9.9,32.1,33.5,29,25.4,15.6,16.5,39.8,18.7,16.4,31.9,27.9,27.3,18.0,19.1,16.3},// asa relativo samanta 2002                                        AsaSamanta
                         {111,250,166,160,157,194,187,86,191,173,179,212,201,208,135,125,144,249,227,149},  //asa win total lins 2003                                                           AsaWinLins
                         {12,35,36,39,3,38,44,22,24,4,5,48,7,6,36,28,26,10,14,5}, // asa folded total lins 2003                                                                                 AsaFtotLins
                         {39,42,48,52,9,50,56,34,33,17,15,57,25,15,45,51,40,16,19,21}, // asa folded coil lins 2003                                                                            AsaFcoilLins
                         {22,38,34,38,12,43,47,19,29,14,15,49,16,16,30,31,35,15,19,16}, // asa folded beta lins 2003                                                                           AsaFbetaLins
                         {66,65,27,30,119,45,49,29,79,137,140,101,163,170,103,36,65,174,130,112},  //asa win pho lins 2003                                                                     AsaWinphoLins
                         {45,187,140,131,38,152,141,56,111,35,39,114,38,37,33,89,79,76,96,36}, // asa win phi lins 2003                                                                        AsaWinphiLins
                         {10,31,29,36,1,34,41,27,25,3,4,47,6,5,39,28,26,8,10,4}, // asa folded pho lins 2003                                                                                   AsafphoLins   "AsafphoLins" "AsafphiLins" "AsafcoilphoLins" "AsafcoilphiLins" "AsafbetaphoLins" "Asafbetaphi" "LinkerRAGeorgetot"
                         {7,35,35,38,0,37,42,14,23,0,1,47,2,0,15,23,19,9,15,1}, // asa folded phi lins 2003             AsafphiLins
                         {41,44,48,56,6,52,62,46,34,15,13,58,23,13,47,55,42,16,17,19},//asa folded coil pho lins 2003                                                                          AsafcoilphoLins
                         {28,42,48,49,10,46,50,27,33,10,11,56,18,10,22,45,36,17,21,11},//asa folded coil phi lins 2003                                                                         AsafcoilphiLins
                         {18,32,30,41,4,34,38,17,28,8,9,46,12,10,33,33,34,11,14,10},//asa folded  beta pho lins 2003                                               AsafbetaphoLins
                         {19,38,33,36,15,42,47,13,28,16,18,49,19,18,11,22,25,24,25,18},//asa folded beta phi                                                        Asafbetaphi
                         {0.964,1.143,0.944,0.916,0.778,1.047,1.051,0.835,1.014,0.922,1.085,0.944,1.032,1.119,1.299,0.947,1.017,0.895,1,0.955},// linker propensities RA george all 2003       LinkerRAGeorgetot
                         {0.974,1.129,0.988,0.892,0.972,1.092,1.054,0.845,0.949,0.928,1.11,0.946,0.923,1.122,1.362,0.932,1.023,0.879,0.902,0.923},// linker propensities RA george 1 2003      LinkerRAGeorge1  "LinkerRAGeorge1" "LinkerRAGeorge2" "LinkerRAGeorge3" "LinkerRAGeosmall" "LinkerRAgeoMed" "LinkerRAlong" "LinkerRAgeorgHel" "LinkerRAgeorgNHel"
                         {0.938,1.137,0.902,0.857,0.6856,0.916,1.139,0.892,1.109,0.986,1,0.952,1.077,1.11,1.266,0.956,1.018,0.971,1.157,0.959 },// linker propensities RA george 2 2003        LinkerRAGeorge2
                         {1.042,1.069,0.828,0.97,0.5,1.111,0.992,0.743,1.034,0.852,1.193,0.979,0.998,0.981,1.332,0.984,0.992,0.96,1.12,1.001},// linker propensities RA george 3 2003          LinkerRAGeorge3
                         {1.065,1.131,0.762,0.836,1.015,0.861,0.736,1.022,0.973,1.189,1.192,0.478,1.369,1.368,1.241,1.097,0.822,1.017,0.836,1.14},// linker propensities RA george small 2003  LinkerRAGeosmall
                         {0.99,1.132,0.873,0.915,0.644,0.999,1.053,0.785,1.054,0.95,1.106,1.003,1.093,1.121,1.314,0.911,0.988,0.939,1.09,0.957},// linker propensities RA george medium 2003   LinkerRAgeoMed
                         {0.892,1.154,1.144,0.925,1.035,1.2,1.115,0.917,0.992,0.817,0.994,0.944,0.782,1.058,1.309,0.986,1.11,0.841,0.866,0.9},// linker propensities RA george long  2003      LinkerRAlong
                         {1.092,1.239,0.927,0.919,0.662,1.124,1.199,0.698,1.012,0.912,1.276,1.008,1.171,1.09,0.8,0.886,0.832,0.981,1.075,0.908},// linker propensities RA george helical 2003  LinkerRAgeorgHel
                         {0.843,1.038,0.956,0.906,0.896,0.968,0.9,0.978,1.05,0.946,0.885,0.893,0.878,1.151,1.818,1.003,1.189,0.852,0.945,0.999},// linker propensities RA george non helical 2003 LinkerRAgeorgNHel
                         {0,-0.1,-1.7,-1.6,-1.4,2.5,-0.7,-1.2,-0.7,-0.5,-0.7,0.1,-0.3,-0.9,-0.4,-1.2,-1.5,-0.6,-1.5,-0.1},// Ncap Doig 2002 NCapdoig   "NCapdoig" "N1Capdoig" "N2CapDoig" "PbioAkashi" "HbioAkashi" "tbioAkashi" "AHBuThompson"
                         {0,-0.1,-1.7,-1.6,-1.4,2.5,-0.7,-1.2,-0.7,-0.5,-0.7,0.1,-0.3,-0.9,-0.4,-1.2,-1.5,-0.6,-1.5,-0.1},// N1 cap Doig 2002 N1Capdoig
                         {0,0.79,1.69,-0.12,0.88,0.52,-0.39,3,0.77,0.58,0.48,0.88,0.67,0.93,3,0.65,0.45,0.83,3,0.49},// N2 cap Doig 2002  N2CapDoig
						   {1,10.7,3.3,1.3,7.3,3.7,2.7,2.3,20.3,4.3,2.7,4.3,9.7,13.3,3.7,2.3,3.3,27.7,13.3,2},// P biosintesis akashi 2002 PbioAkashi
						   {5.3,8.3,5.7,5.7,8.7,6.3,6.6,4.7,9,14,12.3,13,12.3,19.3,8.3,4.7,7.7,23.3,18.3,10.7},// H biosintesis akashi 2002 HbioAkashi
						   {11.7,27.3,14.7,12.7,24.7,16.3,15.3,11.7,38.3,32.3,27.3,30.3,34.3,52,20.3,11.7,18.7,74.3,50,23.3},// total biosintesis akashi 2002 tbioAkashi
						   {69,21,0,-22,-24,15,17,-52,30,49,74,-28,45,44,-76,-9,-21,44,21,34},// helice alfa buried thompson goldstein 1996  AHBuThompson aqui
						   {23,-16,3,-2,52,7,-25,1,22,79,40,-23,55,62,-39,33,35,72,65,97},// Beta buried thompson goldstein 1996             BBuThompson "BBuThompson" "TBuThompson" "cBuThompson" "AHexThompson" "BexThompson" "TexThompson"
						   {32,-17,25,18,40,-11,-10,114,40,-6,22,-70,39,53,77,13,52,78,-27,-23},// turn buried thompson goldstein 1996       TBuThompson
						   {0,2,63,68,79,21,-27,51,55,22,32,6,-10,37,89,68,59,5,34,-10},//coil buried thompson goldstein 1996                cBuThompson
						   {-34,3,-73,-40,-97,-18,44,-203,-61,-91,-133,50,-106,-183,-128,-73,-117,-323,-77,-184},//helice exposed thompson goldstein 1996  AHexThompson
						   {-108,63,-61,-64,-24,30,21,-105,-4,-180,-116,62,-94,-90,-41,-41,50,-63,29,-11},//Beta exposed thompson goldstein 1996    BexThompson
						   {-109,-71,-10,-4,-101,-67,-24,37,-81,-166,-159,-29,-102,-147,-2,-21,-86,-158,-106,-140},//turn exposed thompson goldstein 1996  TexThompson aqui
                         {-114,-4,6,15,-29,11,-42,-37,-87,-122,-108,-11,-20,-90,45,-28,-18,-114,-95,-139},//coil exposed thompson goldstein 1996       cexpThompson "cexpThompson" "AHaJiang" "AHcJiang" "AHBBJiang" "AHCBJiang" "AHhaBAJiang" "AHBAJiang"  "AHCBAJiang" "AHaBiaJiang" "AHBIAjiang" "AHCBIAJiang"
                         {1.2,1.18,0.86,0.87,0.81,1.07,1.23,0.59,1.17,1.08,1.18,0.98,1.19,1.07,0.45,0.8,0.85,1.04,0.93,1.13},// ha alfa proteins jiang 1997  AHaJiang
                         {0.72,0.75,1.23,1.21,1.25,0.88,0.65,1.65,0.77,0.88,0.73,1.03,0.72,0.88,1.87,1.3,1.18,0.94,1.05,0.73},// c alfa proteins jiang 1997  AHcJiang
                         {0.92,0.99,0.61,0.61,1.27,0.88,0.77,0.64,1.06,1.52,1.27,0.98,1.24,1.41,0.36,0.84,1.08,1.27,1.39,1.51},// B B proteins jiang 1997    AHBBJiang
                         {1.02,0.98,1.32,1.29,0.82,1.04,1.13,1.3,1,0.6,0.76,1,0.78,0.68,1.52,1.13,0.98,0.81,0.71,0.62},// C B proteins jiang 1997            AHCBJiang
                         {1.43,1.28,0.71,0.9,0.6,1.35,1.39,0.4,0.77,1.03,1.29,1.13,1.18,0.91,0.45,0.8,0.85,1.1,0.95,0.94},// Ha B+A proteins jiang 1997      AHhaBAJiang
                         {0.67,0.63,0.41,0.46,1.37,0.58,0.67,0.75,0.97,2.01,1.26,0.73,0.96,1.76,0.55,0.72,1.05,1.19,1.38,2.09},// B B+A proteins jiang 1997  AHBAJiang
                         {0.78,0.91,1.43,1.27,1.19,0.87,0.81,1.56,1.19,0.63,0.68,0.99,0.87,0.81,1.59,1.25,1.10,0.86,0.91,0.67},// C B+A proteins jiang 1997  AHCBAJiang
                         {1.62,1.41,0.82,0.97,0.93,1.51,1.52,0.33,0.76,0.95,1.56,1.01,1.61,0.80,0.31,0.66,0.67,0.70,0.83,0.95},// Ha BIA proteins jiang 1997 AHaBiaJiang
                         {0.67,0.82,0.55,0.62,1.02,0.91,0.82,0.71,0.79,1.74,1.07,0.87,1.06,1.76,0.38,0.82,1.32,1.47,1.27,1.78},// B BIA proteins jiang 1997  AHBIAjiang
                         {0.86,0.89,1.30,1.20,1.03,0.80,0.84,1.46,1.21,0.68,0.70,1.06,0.68,0.74,1.62,1.25,1.01,0.92,0.96,0.66},// C BIA proteins jiang 1997  AHCBIAJiang
                         {90.1,192.8,127.5,117.1,103.5,149.4,140.8,63.8,159.3,164.9,164.6,170,167.7,193.5,123.1,94.2,120,231.7,197.1,139.1},// volumen Harpaz 1994 tsai 1999 VHarpaz "VHarpaz" "VPontius" " VTsai" "RaSDCootes" "MWFasman" "IsoPSteinb"
		                 {91.5,196.1,138.3,135.2,102.4,156.4,154.6,67.5,162.3,162.6,163.4,162.5,165.9,198.8,123.4,102.0,126.0,237.2,209.8,138.4},// volumen Pontius 1996 tsai 1999 VPontius
		                 {90.0,194.0,124.7,117.3,103.3,149.4,142.2,64.9,160.0,163.9,164.0,167.3,167.0,191.9,122.9,95.4,121.5,228.2,197.0,139.0},// volumen BL tsai 1999       VTsai
		                 {0,2.63,1.59,1.57,0.91,1.86,1.84,0,1.83,1.68,1.64,2.27,1.86,2.03,1.31,0.71,1.33,2.57,2.37,1.37},// radio side chain Cootes 1998  RaSDCootes
		                 {89.09,174.20,132.12,133.10,121.15,146.15,147.13,75.07,155.16,131.17,131.17,146.19,149.21,165.19,115.13,105.09,119.12,204.24,181.19,117.15},// molecular weight fasman 1976 MWFasman
		                 {6.00,10.76,5.41,2.77,5.07,5.65,3.22,5.97,7.59,6.02,5.98,9.74,5.74,5.48,6.30,5.68,6.16,5.89,5.66,5.96},// isoelectric point  Alff Steinberger IsoPSteinb aqui
		                 //{89.09,174.20,132.12,133.10,121.15,146.15,147.13,75.07,155.16,131.17,131.17,146.19,149.21,165.19,115.13,105.09,119.12,204.24,181.19,117.15};// molecular weight fasman 1976
		                 {1.80,12.50,-5.60,5.05,-16.50,6.30,12.00,0.00,-38.50,12.40,-11.00,14.60,-10.00,-34.50,-86.20,-7.50,-28.00,-33.70,-10.00,5.63},// Optical rotation (Fasman, 1976) OptRFasman "OptRFasman" "RfMcMeekin" "PKNFasman" "PKCFasman" "MPFasman"
		                 {4.34,26.66,13.28,12.00,35.77,17.56,17.26,0.00,21.81,19.06,18.78,21.29,21.64,29.40,10.93,6.35,11.01,42.53,31.53,13.92},// refractividad (McMeekin et al., 1964) RfMcMeekin
		                 {9.69,8.99,8.80,9.60,8.35,9.13,9.67,9.78,9.17,9.68,9.60,9.18,9.21,9.18,10.64,9.21,9.10,9.44,9.11,9.62},//pK-N (Fasman, 1976) PKNFasman
		                 {2.34,1.82,2.02,1.88,1.92,2.17,2.10,2.35,1.82,2.36,2.36,2.16,2.28,2.16,1.95,2.19,2.09,2.43,2.20,2.32},//pK-C (Fasman, 1976) PKCFasman
		                 {297,238,236,270,178,185,249,290,277,284,337,224,283,284,222,228,253,282,344,293},//Melting point (Fasman, 1976) MPFasman
                         {11.50,14.28,12.82,11.68,13.46,14.45,13.57,3.40,13.69,21.40,21.40,15.71,16.25,19.80,17.43,9.47,15.77,21.67,18.03,21.57},//Bulkiness (Zimmerman et al., 1968)BZimmerman "BZimmerman" "XECollagPersikov" "YECollagPersikov" "xpECollagPersikov" "ypECollagPerskov"
                         {480,520,502,520,423,565,590,575,580,624,437,540,452,514,435,506,506,593,629,518},// X entalpia Collagen Persikov 2000 XECollagPersikov
                         {502,610,640,776,471,559,630,665,497,559,514,400,436,557,435,435,647,670,657,481},// Y entalpia Collagen Persikov 2000 YECollagPersikov
                         {11.1,2.8,2.1,4.9,0,2.9,13,1.6,1.6,2,7.8,3.6,0.9,3,32.9,4.9,1.8,0,0.5,2.6},// x p entalpia Collagen Persikov 2000  xpECollagPersikov
                         {10.6,11.4,2.1,4.8,0,6.9,2,0.7,0.5,2.1,1.7,9.0,9.0,0.2,34,4,4.2,0,0,4.3},// y p entalpia Collagen Persikov 2000    ypECollagPerskov aqui
                         {0.24,3.52,3.05,3.98,0.84,1.75,3.11,2.05,2.47,-3.89,-4.28,2.29,-2.85,-4.22,-1.66,2.39,0.75,-4.36,-2.54,-2.59},// ZZ1 Sandberg 1998 De genst ZZ1Sandberg "ZZ1Sandberg" "ZZ2Sandberg" "ZZ3Sandberg"
	                     {-2.32,2.50,1.62,0.93,-1.67,0.50,0.26,-4.06,1.95,-1.73,-1.30,0.89,-0.22,1.94,0.27,-1.07,-2.18,3.94,2.44,-2.64},// ZZ2 Sandberg 1998 De genst ZZ2Sandberg
                         {0.60,-3.50,1.04,1.93,3.71,-1.44,-0.11,0.36,0.26,-1.71,-1.49,-2.49,0.47,1.06,1.84,1.15,-1.12,0.59,0.43,-1.54},// ZZ3 Sandberg 1998 De genst  ZZ3Sandberg
                  };
                 vector<double> Vcodes;
				  for (int ii = 0; ii < AA[j].size(); ii++)
				       {
							Vcodes.push_back(AA[j][ii]);

								  }
					 return Vcodes;
             }

	string 	AApropName(int j)	{
	  vector< string > AAPnames {"Hkyte", "HSCTakano", "Hwoese", "Hpace", "HFauchere", "HWimley", "HKarplus", "HAbraham","HBlack",    // "HKrigbaum",
	                           "HGuy", "HMiller", "HvHeijne", "HRoseman", "HLawsonCH","Hmiyazawa", "Hponnuswamy",
	                           "Hkuntz", "HBull", "HMeek7", "HMeek2", "HCowan3","HCowan7", "HTMtPilpel","HTMePilpel","HTMcPilpel","HTMiPilpel","HTMtPilpel",
	                           "HAdamian","Cadamian", "MGromiha","LCGromiha1120","LCGromiha2130","LCGromiha3140","LCGromiha4150","LCromihaM501","LCromihaM502",
	                           "scFBt","scFBRR","scFBRF","scFBFF", "GaaDna_pr","AaaDna_pr","TaaDna_pr","CaaDna_pr","IP_Pkeskin","IP_PkeskinI","IP_PkeskinII",
	                           "IP_PkeskinIII","HHopp","HParker","HWilson","HCornette","HEisenberg","Hlevitt","HZIMMERMAN","HWolfenden","HOoi","HSamatey",
	                           "HBeuming","MLFBt","MLFB2R","MLFBRF","MLFBFF",  "AHAgadir",  "HLawsonET", "BA_32", "AHBlaber", "AHChou", "AHChakra","AHOneil","AHDel","AHMPRao", "AHpLiuOct", "AHpLiuH2O",
	                           "AHpLiuStats", "AHRohl", "AHpYang","AHpW","AHRichardson","AHWilliams","AHTmGromiha","AHpWarren","AHMunoz","AHPPMunoz","AHPhiplusBahar","AHPhiMinusBahar","AHPhiPlusMinBahar","AHEaPhiPlusMinBahar",
	                           "BsChou","BDeleage","BmKim1","BmKim2","BesEswar","BisEswar","BexEswar","Bmunoz","BppMunoz",    "BRMunoz","BLinding","CoilDeleage","TurnDeleage","TmAHMonne","TmAHMonne2",
	                             "BspDaffner","BsiDaffner", "BsnDaffner","TurnpDaffner","LoopEswar",     "LtWrabl","MtWrabl","GtWrabl","Ha1Shortle","Ha2Shortle","Ha3Shortle",
	                             "L1Shortle","L2Shortle","L3Shortle","m1Shortle","m2Shortle","m3Shortle","r1Shortle","r2Shortle","r3Shortle","LHShortle",
	                             "phiShortle","othShortle",     "ratec-tReimer","ratet-cReimer","TurnSantai","TurnSantai1","TurnSantai2","TurnSantai3",     "PolypEswar",
	                             "polypIIStapley","CubellisL3","CubellisLp3","CubellisLm3","CubellisLm2","CubellisLm1", "Cubellisn1","Cubellisn2","Cubellismid","Cubellisc2",
	                             "Cubellisc1","Cubellisc+1","Cubellisc+2","Cubellisc+3",   "HPersson","EntDoig","EntAbsDoig","EntAbagyan","EntKoehl","EntPickett",
	                             "EntLee","EntHaAvbelj","EntAvbelj","DeltaEntAvbelj","HaEntAvbelj","SEAvbelj",    "mBEntAvbelj","EntCoilJha","EntZamanUA","EntZamanAA","EntZamanGS",
	                             "EntZamanAmb","AsaHelFan","AsaSFan",       "AsaCoilFan","NcontSamanta","AsaSamanta","AsaWinLins","AsaFtotLins","AsaFcoilLins","AsaFbetaLins",
	                             "AsaWinphoLins","AsaWinphiLins",       "AsafphoLins","AsafphiLins","AsafcoilphoLins","AsafcoilphiLins","AsafbetaphoLins","Asafbetaphi",
	                             "LinkerRAGeorgetot","LinkerRAGeorge1","LinkerRAGeorge2","LinkerRAGeorge3","LinkerRAGeosmall","LinkerRAgeoMed","LinkerRAlong",
	                             "LinkerRAgeorgHel","LinkerRAgeorgNHel",     "NCapdoig","N1Capdoig","N2CapDoig","PbioAkashi","HbioAkashi","tbioAkashi","AHBuThompson",

	                             "BBuThompson","TBuThompson","cBuThompson","AHexThompson","BexThompson",   "TexThompson", "cexpThompson",
	                             "AHaJiang","AHcJiang","AHBBJiang","AHCBJiang","AHhaBAJiang","AHBAJiang","AHCBAJiang","AHaBiaJiang","AHBIAjiang","AHCBIAJiang",

	                             "VHarpaz","VPontius","VTsai","RaSDCootes","MWFasman","IsoPSteinb", "OptRFasman","RfMcMeekin","PKNFasman","PKCFasman","MPFasman",

	                             "BZimmerman","XECollagPersikov","YECollagPersikov","xpECollagPersikov","ypECollagPerskov","ZZ1Sandberg","ZZ2Sandberg","ZZ3Sandberg"};

			 return  AAPnames[j];
	   }

//int main()
//{
//   vector<vector<int> >  ms, ms1;
//   vector<int> sii{ 1, 2, 3, 4, 5,6,7,8,9 };
//   ms=randppp(sii, 10);
//   ms1=RandWOutRep(sii, 10);
//    for (int i = 0; i < ms1.size(); i++)
//     {
//    for (int ii = 0; ii < ms1[0].size(); ii++)
//         {
//        cout << ' ' << ms1[i][ii];
//    }
//    cout << endl;
//     }
//
//    return 0;
//}
