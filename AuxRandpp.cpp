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



