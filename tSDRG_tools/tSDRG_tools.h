#ifndef _tSDRG_TOOLS_H_
#define _tSDRG_TOOLS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <map>
#include <algorithm>
#include <uni10.hpp>
#include "../MPO/mpo.h"

using namespace std;

double Distribution_Random_Variable(int model, double var, double Jdis);

void Truncation(uni10::Matrix<double>& En, uni10::Matrix<double>& state, const int chi, int& m, bool& info);

double find_highest_gap(uni10::Matrix<double> En, const int chi, bool& info);

void tSDRG_OBC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool save_RG_info, bool& info);

void tSDRG_OBC_regular(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool save_RG_info, bool& info);

void tSDRG_PBC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool save_RG_info, bool& info);

void tSDRG_PBC_regular(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool save_RG_info, bool& info);

/* 
void tSDRG0(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info);

void tSDRG(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info);

void tSDRG2(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info);

void tSDRG_PBC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info);

void tSDRG2_PBC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info);

void tSDRG_PBC_layer(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info);

void tSDRG_PBC_layerlowest(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info);

void tSDRG_PBC_layer2(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info);

void tSDRG_PBC_layerBig(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info);

void tSDRG_PBC_layerBigC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info);

void tSDRG_PBC_layergg(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info);

void tSDRG0_PBC_test(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info);
 */
#endif