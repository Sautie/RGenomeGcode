#ifndef AUXRANDPP_H
#define AUXRANDPP_H
#include <bits/stdc++.h>
#include <Windows.h>
using namespace std;

vector<vector<int> >  randppp(vector<int> s, int m);
vector<vector<int> >  randInvSwap(vector<vector<int> >  s, int m);
long int  funcPos(vector<int> s);
vector<vector<int> >  RandWOutRep(vector<int> s, int m);
void EqGenerator(string dest1, vector<string> vect20, string dest2, vector<string> vect64);
void read_folder(const string& name, vector <string> & list);
vector< string > GenCodes(string GCname);
vector<double> 	AAprop(int j);
string 	AApropName(int j);


#endif