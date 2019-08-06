#ifndef _MEASURE_H_
#define _MEASURE_H_

#include <iostream>
#include <random>
#include <cmath>
#include <sstream>
#include <map>
#include <algorithm>
#include <uni10.hpp>
#include "../Operator/operator.h"
#include "../Hamiltonian/hamiltonian.h"
#include "../MPO/mpo.h"

using namespace std;

string Decision_tree(vector<int> Vs_loc, bool show);

vector<double> Schmidt_Value(uni10::UniTensor<double> top);

vector<double> Energy_Spectrum(vector<MPO>& MPO_chain, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_SxSx(int Si,int Sj,vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_SySy(int Si,int Sj,vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_SzSz(int Si,int Sj, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_StSt(int Si,int Sj,vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_String(int Si,int Sj, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_Sx(int Si, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_iSy(int Si, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_Sz(int Si, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Correlation_St(int Si, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Magnetization_Sz(vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

double Magnetization_Sx(vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc);

#endif