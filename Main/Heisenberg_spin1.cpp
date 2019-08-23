#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <uni10.hpp>
#include "../MPO/mpo.h"
#include "../tSDRG_tools/tSDRG_tools.h"
#include "../tSDRG_tools/measure_PBC.h"

using namespace std;

void tSDRG_XXZ(int L, int chi, int Pdis, double Jdis, string BC, double S, double Jz, double h, int Jseed)
{
    random_device rd;            // non-deterministic generator.
    mt19937 genRandom(rd() );    // use mersenne twister and seed is rd.
    mt19937 genFixed(Jseed);     // use mersenne twister and seed is fixed!

    uniform_real_distribution<double> Dist_J(nextafter(0.0, 1.0), 1.0); // probability distribution of J rand(0^+ ~ 1)

    /// create coupling list and MPO chain for OBC or PBC
    vector<double> J_list;
    vector<MPO> MPO_chain;
    if (BC == "PBC")
    {
        for(int i=0; i<L; i++)
        {
            double jvar = Dist_J(genFixed);
            jvar = Distribution_Random_Variable(Pdis, jvar, Jdis);
            J_list.push_back(jvar);
            MPO W("XXZ_PBC", 'm', S, J_list[i], Jz*J_list[i], h);
            MPO_chain.push_back(W);
        }
    }
    else if (BC == "OBC")
    {
        for(int i=0; i<L; i++)
        {
            double jvar = Dist_J(genFixed);
            jvar = Distribution_Random_Variable(Pdis, jvar, Jdis);

            if(i == 0)
            {
                J_list.push_back(jvar);     //random chain 
                MPO W("XXZ_OBC", 'l', S, J_list[i], Jz*J_list[i], h);
                MPO_chain.push_back(W);
            }
            else if (i > 0 && i < L-1)
            {
                J_list.push_back(jvar);     //random chain 
                MPO W("XXZ_OBC", 'm', S, J_list[i], Jz*J_list[i], h);
                MPO_chain.push_back(W);
            }
            else
            {
                MPO W("XXZ_OBC", 'r', S, 1.0, Jz, h);
                MPO_chain.push_back(W);
            }
        }
    }
    else
    {
        ostringstream err;
        err << "tSDRG support OBC and PBC only";
        throw runtime_error(err.str());
    }


    /// Save tree tensor network
    string file, dis, folder, file_name, file_name1, file_name2, file_nameS;
    if (Jdis < 1.0)
        dis = "0" + to_string( (int)(Jdis*10) );
    else if (Jdis == 1.0)
        dis = "10";
    else if (Jdis == 2.0)
        dis = "20";
    else if (Jdis == 3.0)
        dis = "30";
    else if (Jdis > 1.0)
        dis = to_string( (int)(Jdis*10) );

    /// create folder
    folder = "TTN/" + BC + "/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    string str = "mkdir -p " + folder;
    const char *mkdir = str.c_str();
    const int dir_err = system(mkdir);
    if (dir_err == -1)
    {
        cout << "Error creating directory!" << endl;
        exit(1);
    }

    vector<uni10::UniTensor<double> > w_up;      // w_up is isometry tensor = VT
    vector<int> w_loc;                           // location of each w
    vector<double> En = J_list;                  // J_list will earse to one, and return ground energy.
    bool info = 1;                               // True; if tSDRG can not find non-zero gap, info return 0, and stop this random seed.
    bool save_RG_info = 0;                       // save gaps at RG stage 
    tSDRG(MPO_chain, En, w_up, w_loc, chi, dis, Pdis, Jseed, save_RG_info, info);

    /// check info if can not RG_J
    if (info == 0)
    {
        cout << "random seed " + to_string(Jseed)  + " died (can not find non-zero gap) " << endl;
        return;
    }
    else
    {
        cout << "finish in " << folder << endl;
    }

    for (int i=0; i<10; i++)
        cout << En[i] << endl;

    //string top1 = Decision_tree(w_loc, true);

    /// create isometry of other part
    /*vector<uni10::UniTensor<double> > w_down;    // w_down
    uni10::UniTensor<double> kara;
    w_down.assign(L-1, kara);
    for (int i = 0; i < w_up.size(); i++)
    {
        w_down[i] = w_up[i];
        uni10::Permute(w_down[i], {-3, -1, 1}, 2, uni10::INPLACE);
    }

    /// save disorder coupling J list 
    vector<double> corr12;            // corr <S1><S2>
    file = folder + "/J_list.csv";
    ofstream fout(file);              // == fout.open(file);
    if (!fout)
    {
        ostringstream err;
        err << "Error: Fail to save (maybe need mkdir " << file << ")";
        throw runtime_error(err.str());
    }
    for (int i = 0; i < J_list.size(); i++)
    {
        fout << setprecision(16) << J_list[i] << endl;
        corr12.push_back(Correlation_St(i, w_up, w_down, w_loc) );
    }
    fout.flush();
    fout.close();*/
    
    /// save tensor
    /*file = folder;
    for (int i = 0; i < w_up.size(); i++)
    {
        file = folder + "/w_up" + to_string(i);
        w_up[i].Save(file);
    }*/

    /// save merge order
    /*file = folder + "/w_loc.csv";
    fout.open(file);
    if (!fout)
    {
        ostringstream err;
        err << "Error: Fail to save (maybe need mkdir " << file << ")";
        throw runtime_error(err.str());
    }
    for (int i = 0; i < w_loc.size(); i++)
    {
        fout << w_loc[i] << endl;
    }
    fout.flush();
    fout.close();

    /// save ground energy
    file = folder + "/energy.csv";
    fout.open(file);
    if (!fout)
    {
        ostringstream err;
        err << "Error: Fail to save (maybe need mkdir " << file << ")";
        throw runtime_error(err.str());
    }
    fout << "energy" << endl;
    for (int i=0; i<En.size(); i++)
    {
        fout << setprecision(16) << En[i] << endl;
    }
    fout.flush();
    fout.close();
    string top1 = Decision_tree(w_loc, true);
    
    /// correlation and string order parameter print
    double corr, corr1 ,corr2, Oz;
    file_name1 = "data/correlation/" + algo + "/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "_corr1.csv";
    file_name2 = "data/correlation/" + algo + "/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "_corr2.csv";
    file_nameS = "data/String/" + algo + "/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "_string.csv";
    ofstream fout1(file_name1);
    ofstream fout2(file_name2);
    ofstream foutS(file_nameS);
    fout1 << "x1,x2,corr" << endl;
    fout2 << "x1,x2,corr" << endl;
    foutS << "x1,x2,corr" << endl;
    int site1;
    int site2;
    int r;
    int loop = 0;
    for (site1=0; site1<1; site1 += 10)
    {
        int step = 2;
        for (site2=site1+4; site2<5; site2 += step)
        {
            r = site2 - site1;
            if (r >= 10 && r < 100)
                step = 4;
            else if (r >= 100)
                step = 10;

            if (r <= L/2)
            {
                cout << "TEST: " << site1 << " & " << site2 << endl;
                corr = Correlation_StSt(site1, site2, w_up, w_down, w_loc);
                corr1 = corr12[site1];
                corr2 = corr12[site2];

                fout1 << setprecision(16) << site1 << "," << site2 << "," << corr << endl;
                fout2 << setprecision(16) << site1 << "," << site2 << "," << corr - corr1*corr2 << endl;

                Oz = Correlation_String(site1, site2, w_up, w_down, w_loc);
                foutS << setprecision(16) << site1 << "," << site2 << "," << Oz << endl;
            }
        }
    }
    fout1.flush();
    fout2.flush();
    foutS.flush();
    fout1.close();
    fout2.close();
    foutS.close();
    
    /// Schmidt value = Entanglement spectrum
    string top = Decision_tree(w_loc, false);       // find location of top tensor
    vector<double> schmidt_value;
    schmidt_value = Schmidt_Value(w_up.back() );
    
    /// print out
    file_name = "data/Schmidt/" + algo + "/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed) + "_value" + top + ".csv";
    fout.open(file_name);
    if (!fout)
    {
        ostringstream err;
        err << "Error: Fail to save (maybe need mkdir " << file_name << ")";
        throw runtime_error(err.str());
    }
    fout << "value" << endl;
    fout.flush();
    for (int i = 0; i < schmidt_value.size(); i++)
    {
        fout << setprecision(16) << schmidt_value[i] << endl;
    }
    fout.flush();
    fout.close();*/
}

void errMsg(char *arg) 
{
    cerr << "Usage: " << arg << " [options]" << endl;
    cerr << "Need 8-parameter:" << endl;
    cerr << "./job.exe <system size> <keep state of RG procedure> <Prob distribution> <disorder> <algo> <seed1> <seed2>\n" << endl;
    cerr << "Example:" << endl;
    cerr << "./job.exe 32 8 10 0.1 PBC 1 1\n" << endl;
}

int main(int argc, char *argv[])
{
    int L;                      // system size
    int chi;                    // keep state of isometry
    string BC;                  // boundary condition
    int Pdis;                   // model of random variable disturbution
    double Jdis;                // J-coupling disorder strength
    int seed1;                  // random seed number in order to repeat data
    int seed2;                  // random seed number in order to repeat data
    double S      = 1.0;        // spin dimension
    double Jz     = 1.0;        // XXZ model
    double h      = 0.0;        // XXZ model

    if (argc == 8)
    {
        stringstream(argv[1]) >> L;
        stringstream(argv[2]) >> chi;
        stringstream(argv[3]) >> Pdis;
        stringstream(argv[4]) >> Jdis;
        BC = argv[5];
        stringstream(argv[6]) >> seed1;
        stringstream(argv[7]) >> seed2;
    }
    else
    {
        errMsg(argv[0]);
        return 1;
    }

    for (int Jseed=seed1; Jseed<=seed2; Jseed++)
    {
        tSDRG_XXZ(L, chi, Pdis, Jdis, BC, S, Jz, h, Jseed);
    }

    return 0;
}
