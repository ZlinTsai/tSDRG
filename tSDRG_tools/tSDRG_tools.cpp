#include "tSDRG_tools.h"

/// return random number by distribution
double Distribution_Random_Variable(int model, double var, double Jdis)
{
    double value;
    switch (model)
    {
        case 0:
        {
            /// non-disorder, all values are one.
            value = 1;
            break;
        }
        case 1:
        {
            /// box-like
            value = var;
            break;
        }
        case 10:
        {
            /// power-law
            value = pow(var, Jdis);
            break;
        }
        default:
        {
            ostringstream err;
            err << "Error: NO this model of disturbution (at tSDRG_tools/tSDRG_tools.cpp)\n";
            err << "Pdis: 0 is no disorder; 1 is box-like; 10 is power-law";
            throw runtime_error(err.str());
            break;
        }
    }
    
    return value;
}

/// isometry keep m-state , m = highest gap of below chi-state
void Truncation(uni10::Matrix<double>& En, uni10::Matrix<double>& state, const int chi, int& m, bool& info)
{
    int idx = chi - 1;  // c++ from 0 start
    m = chi;
    vector<double> gap;
    
    for (int i=0; i<En.col()-1; i++)
    {
        gap.push_back(En.At(i+1, i+1) - En.At(i, i) );  // At(i, j)
        //cout << setprecision(8) << gap[i] << endl;
    }
    
    while (gap[idx] <= pow(10, -12))
    {
        idx--;
        m--;
        if (idx < 0)
        {
            info = 0;
            return;
        }
    }
    
    uni10::Resize(En, m, m, uni10::INPLACE);
    uni10::Resize(state, m, state.col(), uni10::INPLACE); // NOTICE: Uni10 give me V_diagger
    
    //cout << En << state;
}

/// renormalization J
double find_highest_gap(uni10::Matrix<double> En, const int chi, bool& info)
{
    vector<double> gap;
    for (int i=0; i<En.col()-1; i++)
    {
        gap.push_back(En.At(i+1, i+1) - En.At(i, i) ); 
        //cout << setprecision(8) << gap[i] << endl;
    }
    //cout << "***********************" << endl;

    int idx;
    if (chi > gap.size() )
        idx = gap.size() - 1;
    else 
        idx = chi - 1;

    /// keep SU2 by discarding/keeping mutiplet
    while (gap[idx] <= pow(10, -12)) // <= pow(10, -20) /// check non-zero gap
    {
        idx--;

        if (idx < 0)
        {
            info = 0;
            return 0;
        }
    }
    double maxj = gap[idx];

    return maxj;
}

/// tSDRG
void tSDRG(vector<MPO>& MPO_chain, vector<double>& Q, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool save_RG_info, bool& info)
{
    /// loading Network file
    uni10::Network H12("../tSDRG_net/H12.net");

    /// find system size 
    int L = MPO_chain.size();

    /// check boundary condition, ex: model = XXZ_PBC, BC = model[-3] + model[-2] + model[-1] = PBC
    string model = MPO_chain[0].model;
    string BC = "";
    for (int i = 3; i >= 1; i--)
        BC += model.at(model.size() - i);

    //cout << BC << endl;

    /// create RG-layer
    vector<int> layer;
    layer.assign(L, 0);
    
    /// save RG infomation
    vector<int> RG_stage;
    vector<double> RG_lowest_gaps;
    vector<double> RG_highest_gaps;

    /// Contract any two site hamitonian in order to find highest gap
    uni10::UniTensor<double> H, H_last, H1, H2;             // H is hamitonian of site-1(= H1) and site-2(= H2)  
    H1 = MPO_chain[1].GetTensor('l');
    H2 = MPO_chain[2].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal local hamitonian
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H_block = H.GetBlock();
    uni10::EigH(H_block, En, state, uni10::INPLACE);
    double coeff = find_highest_gap(En, chi, info);

    /// transfor coupling to energy scale (highest gap); energy spactrum of S = 1 is -2J, -J, J, so that highest gaps = 2J
    for (int i=0; i<Q.size(); i++)
    {
        Q[i] = coeff * Q[i]; 
        //cout << Q[i] << endl;
    }
    
    while(MPO_chain.size() > 2)
    {
        /// find max gap = max{Q}
        auto Qmax = max_element(Q.begin(), Q.end() );      // max_element return iterators, not values. ( double Qmax = *Qmax; use *iterators to get values.
        int s1 = distance(Q.begin(), Qmax);                // return distance from 0 to Qmax. Note: distance(iterators, iterators) input is iter NOT value.
        int s2 = (s1 == MPO_chain.size()-1 ? 0 : s1 + 1);          // s1 = site1 and "coupling label", s2 = site2, for PBC s2+1 = 0 at last site.
        int chi_loc;                                       // chi of location (chi maybe cut at mutliplet, so chi will change).
        
        /// check 
        //cout << "MPO size = " << MPO_chain.size() << " , Qmax = " << *Qmax <<  ", MPO localtion = " << s1 << " and " << s2 << endl;
        //for (int i=0; i<Q.size(); i++)
            //cout << setprecision(16) << Q[i] << endl;
        
        /// Contract two site hamitonian
        H1 = MPO_chain[s1].GetTensor('l');
        H2 = MPO_chain[s2].GetTensor('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H);
        
        /// diagonal local hamitonian
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// save RG infomation
        if (save_RG_info)
        {
            layer[s1] = max(layer[s1], layer[s2]) + 1;
            layer.erase(layer.begin() + s2);

            double gap_lowest = abs(En.At(1, 1) - En.At(0, 0));
            RG_stage.push_back(layer[s1]);
            RG_highest_gaps.push_back(*Qmax);
            RG_lowest_gaps.push_back(gap_lowest);
        }

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        /// can not find non-zero gap, return
        if (info == 0)
            return;

        /// create isometry VT (V_Transpose)
        uni10::Bond leg_up(uni10::BD_IN, chi_loc);
        uni10::Bond leg_L(uni10::BD_OUT, H1.GetBlock().row() );
        uni10::Bond leg_R(uni10::BD_OUT, H2.GetBlock().col() );
        vector<uni10::Bond> bondVT;
        bondVT.push_back(leg_up); bondVT.push_back(leg_L); bondVT.push_back(leg_R);
        uni10::UniTensor<double> VT(bondVT);
        VT.PutBlock(state);
        VT.SetLabel({1, -1, -3});

        /// w_down = V
        uni10::UniTensor<double> V = VT;
        uni10::Permute(V, {-1, -3, 1}, 2, uni10::INPLACE);

        /// merge two site MPO, Note: this order is import for V and VT (check dim of tensor's legs)
        RG_MPO(MPO_chain[s1], MPO_chain[s2], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(s1);

        /// erase s2
        if (BC == "PBC")
        {
            MPO_chain.erase(MPO_chain.begin() + s2);
        }
        else if (BC == "OBC")
        {
            /// erase s2, except s2 is most left mpo 
            if (s2 == MPO_chain.size()-1)
                MPO_chain.erase(MPO_chain.begin() + s1);
            else
                MPO_chain.erase(MPO_chain.begin() + s2);
        }
        else
        {
            ostringstream err;
            err << "tSDRG support OBC and PBC only";
            throw runtime_error(err.str());
        }
        
        /// RG coupling into merge list, notice that eage mpo is OBC or PBC
        if (MPO_chain.size() != 2)
        {
            if (BC == "PBC")
            {
                if (s1 == 0)
                {
                    H1 = MPO_chain[Q.size()-2].GetTensor('l');   // -2 by erase MPO[s2] above line
                    H2 = MPO_chain[0         ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[Q.size()-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[0].GetTensor('l');
                    H2 = MPO_chain[1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[1] = find_highest_gap(En, chi, info);
                }
                else if (s1 == Q.size()-2)
                {
                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[s1  ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[s1].GetTensor('l');
                    H2 = MPO_chain[0 ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1+1] = find_highest_gap(En, chi, info);
                }
                else if (s1 == Q.size()-1)
                {
                    H1 = MPO_chain[s1-2].GetTensor('l');
                    H2 = MPO_chain[s1-1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[0   ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[0] = find_highest_gap(En, chi, info);
                }
                else
                {
                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[s1  ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[s1  ].GetTensor('l');
                    H2 = MPO_chain[s1+1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1+1] = find_highest_gap(En, chi, info);
                }
            }
        
            else if (BC == "OBC")
            {
                if (s1 != 0)
                {
                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[s1  ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);   // in mkl use uni10::EigHLazy()
                    Q[s1-1] = find_highest_gap(En, chi, info);
                }
                if (s1 != Q.size()-1)
                {
                    H1 = MPO_chain[s1  ].GetTensor('l');
                    H2 = MPO_chain[s1+1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);  // in mkl use uni10::EigHLazy()
                    Q[s1+1] = find_highest_gap(En, chi, info);
                }
            }
        }
        
        if (info == 0)
            return;
        
        /// erase coupling between s1 and s2
        Q.erase(Q.begin() + s1);
    }

    /// Contract last two site hamitonian
    if (BC == "PBC")
    {
        /// 1 and 0; use spcal tools GetTensorPBC for getting real term
        H1 = MPO_chain[1].GetTensorPBC('l');
        H2 = MPO_chain[0].GetTensorPBC('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H_last);
        //H_last.PrintDiagram();
        uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);
    }

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal
    if (BC == "PBC")
        H_block = H.GetBlock() + H_last.GetBlock();
    else if (BC == "OBC")
        H_block = H.GetBlock();
    uni10::EigH(H_block, En, state, uni10::INPLACE);
    uni10::Resize(state, 1, state.col(), uni10::INPLACE);

    /// create isometry VT = V_Transpose
    uni10::Bond leg_up(uni10::BD_IN, 1);
    uni10::Bond leg_L(uni10::BD_OUT, H1.GetBlock().row() );
    uni10::Bond leg_R(uni10::BD_OUT, H2.GetBlock().col() );
    vector<uni10::Bond> bondVT;
    bondVT.push_back(leg_up); bondVT.push_back(leg_L); bondVT.push_back(leg_R);
    uni10::UniTensor<double> VT(bondVT);
    VT.PutBlock(state);
    VT.SetLabel({1, -1, -3});

    /// uni10::permute to {up, R, L} (rotation 90) and save it
    uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
    VTs.push_back(VT);
    Vs_loc.push_back(0);

    /// return ground state energy save in queue Q
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );

    Q.clear();
    Q = energy;

    if (save_RG_info)
    {
        /// save last RG info
        double lowest_gap = energy[1] - energy[0];
        layer[0] = max(layer[0], layer[1]) + 1;
        layer.erase(layer.begin() + 1);
        RG_stage.push_back(layer[0]);
        RG_highest_gaps.push_back(lowest_gap);
        RG_lowest_gaps.push_back(lowest_gap);

        /// saving data into folder
        string file1, file2, folder;
        folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
        file1 = folder + "/RG_highest_gaps.csv";
        file2 = folder + "/RG_lowest_gaps.csv";

        ofstream fout1(file1);
        ofstream fout2(file2);
        fout1 << "layer,gap" << endl;
        fout2 << "layer,gap" << endl;

        for (int i = 0; i < RG_stage.size(); i++)
        {
            fout1 << setprecision(16) << RG_stage[i] << "," << RG_highest_gaps[i] << endl;
            fout2 << setprecision(16) << RG_stage[i] << "," << RG_lowest_gaps[i] << endl;
        }
    }
    
}

/// tSDRG with regular TTN as pyramid
void tSDRG_regular(vector<MPO>& MPO_chain, vector<double>& Q, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool save_RG_info, bool& info)
{
    /// loading Network file
    uni10::Network H12("../tSDRG_net/H12.net");

    /// find system size 
    int L = MPO_chain.size();

    /// check boundary condition, ex: model = XXZ_PBC, BC = model[-3] + model[-2] + model[-1] = PBC
    string model = MPO_chain[0].model;
    string BC = "";
    for (int i = 3; i >= 1; i--)
        BC += model.at(model.size() - i);

    //cout << BC << endl;

    // create regular TTN
    vector<int> TTN_Reg;
    int Q_idx = 0;
    if (L%4 == 0)
    {
        int temp  = 0;
        int layer = 1;
        int num_tensor = 0;
        for (int i=0; i<L-1; i++)
        {
            TTN_Reg.push_back(temp - num_tensor);
            temp += 2;
            num_tensor++;
            if (temp == L/layer)
            {
                temp = 0;
                num_tensor = 0;
                layer *= 2;
            }
        }
    }
    else
    {   
        ostringstream err;
        err << "No support for (L%4 != 0)";
        throw runtime_error(err.str());
    }

    /// create RG-layer
    vector<int> layer;
    layer.assign(L, 0);
    
    /// save RG infomation
    vector<int> RG_stage;
    vector<double> RG_lowest_gaps;
    vector<double> RG_highest_gaps;

    /// Contract any two site hamitonian in order to find highest gap
    uni10::UniTensor<double> H, H_last, H1, H2;             // H is hamitonian of site-1(= H1) and site-2(= H2)  
    H1 = MPO_chain[1].GetTensor('l');
    H2 = MPO_chain[2].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal local hamitonian
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H_block = H.GetBlock();
    uni10::EigH(H_block, En, state, uni10::INPLACE);
    double coeff = find_highest_gap(En, chi, info);

    /// transfor coupling to energy scale (highest gap); energy spactrum of S = 1 is -2J, -J, J, so that highest gaps = 2J
    for (int i=0; i<Q.size(); i++)
    {
        Q[i] = coeff * Q[i]; 
        //cout << Q[i] << endl;
    }
    
    while(MPO_chain.size() > 2)
    {
        /// find max gap = max{Q}
        auto Qmax = max_element(Q.begin(), Q.end() );      // max_element return iterators, not values. ( double Qmax = *Qmax; use *iterators to get values.
        int s1 = TTN_Reg[Q_idx];                           // return distance from 0 to Qmax. Note: distance(iterators, iterators) input is iter NOT value.
        int s2 = (s1 == MPO_chain.size()-1 ? 0 : s1 + 1);  // s1 = site1 and "coupling label", s2 = site2, for PBC s2+1 = 0 at last site.
        int chi_loc;                                       // chi of location (chi maybe cut at mutliplet, so chi will change).
        Q_idx++;
        
        /// check 
        //cout << "MPO size = " << MPO_chain.size() << " , Qmax = " << *Qmax <<  ", MPO localtion = " << s1 << " and " << s2 << endl;
        //for (int i=0; i<Q.size(); i++)
            //cout << setprecision(16) << Q[i] << endl;
        
        /// Contract two site hamitonian
        H1 = MPO_chain[s1].GetTensor('l');
        H2 = MPO_chain[s2].GetTensor('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H);
        
        /// diagonal local hamitonian
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// save RG infomation
        if (save_RG_info)
        {
            layer[s1] = max(layer[s1], layer[s2]) + 1;
            layer.erase(layer.begin() + s2);

            double gap_lowest = abs(En.At(1, 1) - En.At(0, 0));
            RG_stage.push_back(layer[s1]);
            RG_highest_gaps.push_back(*Qmax);
            RG_lowest_gaps.push_back(gap_lowest);
        }

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        /// can not find non-zero gap, return
        if (info == 0)
            return;

        /// create isometry VT (V_Transpose)
        uni10::Bond leg_up(uni10::BD_IN, chi_loc);
        uni10::Bond leg_L(uni10::BD_OUT, H1.GetBlock().row() );
        uni10::Bond leg_R(uni10::BD_OUT, H2.GetBlock().col() );
        vector<uni10::Bond> bondVT;
        bondVT.push_back(leg_up); bondVT.push_back(leg_L); bondVT.push_back(leg_R);
        uni10::UniTensor<double> VT(bondVT);
        VT.PutBlock(state);
        VT.SetLabel({1, -1, -3});

        /// w_down = V
        uni10::UniTensor<double> V = VT;
        uni10::Permute(V, {-1, -3, 1}, 2, uni10::INPLACE);

        /// merge two site MPO, Note: this order is import for V and VT (check dim of tensor's legs)
        RG_MPO(MPO_chain[s1], MPO_chain[s2], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(s1);

        /// erase s2
        if (BC == "PBC")
        {
            MPO_chain.erase(MPO_chain.begin() + s2);
        }
        else if (BC == "OBC")
        {
            /// erase s2, except s2 is most left mpo 
            if (s2 == MPO_chain.size()-1)
                MPO_chain.erase(MPO_chain.begin() + s1);
            else
                MPO_chain.erase(MPO_chain.begin() + s2);
        }
        else
        {
            ostringstream err;
            err << "tSDRG support OBC and PBC only";
            throw runtime_error(err.str());
        }
        
        /// RG coupling into merge list, notice that eage mpo is OBC or PBC
        if (MPO_chain.size() != 2)
        {
            if (BC == "PBC")
            {
                if (s1 == 0)
                {
                    H1 = MPO_chain[Q.size()-2].GetTensor('l');   // -2 by erase MPO[s2] above line
                    H2 = MPO_chain[0         ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[Q.size()-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[0].GetTensor('l');
                    H2 = MPO_chain[1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[1] = find_highest_gap(En, chi, info);
                }
                else if (s1 == Q.size()-2)
                {
                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[s1  ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[s1].GetTensor('l');
                    H2 = MPO_chain[0 ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1+1] = find_highest_gap(En, chi, info);
                }
                else if (s1 == Q.size()-1)
                {
                    H1 = MPO_chain[s1-2].GetTensor('l');
                    H2 = MPO_chain[s1-1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[0   ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[0] = find_highest_gap(En, chi, info);
                }
                else
                {
                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[s1  ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1-1] = find_highest_gap(En, chi, info);

                    H1 = MPO_chain[s1  ].GetTensor('l');
                    H2 = MPO_chain[s1+1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                    Q[s1+1] = find_highest_gap(En, chi, info);
                }
            }
        
            else if (BC == "OBC")
            {
                if (s1 != 0)
                {
                    H1 = MPO_chain[s1-1].GetTensor('l');
                    H2 = MPO_chain[s1  ].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);   // in mkl use uni10::EigHLazy()
                    Q[s1-1] = find_highest_gap(En, chi, info);
                }
                if (s1 != Q.size()-1)
                {
                    H1 = MPO_chain[s1  ].GetTensor('l');
                    H2 = MPO_chain[s1+1].GetTensor('r');
                    H12.PutTensor("H1", H1);
                    H12.PutTensor("H2", H2);
                    H12.Launch(H);
                    uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);  // in mkl use uni10::EigHLazy()
                    Q[s1+1] = find_highest_gap(En, chi, info);
                }
            }
        }
        
        if (info == 0)
            return;
        
        /// erase coupling between s1 and s2
        Q.erase(Q.begin() + s1);
    }

    /// Contract last two site hamitonian
    if (BC == "PBC")
    {
        /// 1 and 0; use spcal tools GetTensorPBC for getting real term
        H1 = MPO_chain[1].GetTensorPBC('l');
        H2 = MPO_chain[0].GetTensorPBC('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H_last);
        //H_last.PrintDiagram();
        uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);
    }

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal
    if (BC == "PBC")
        H_block = H.GetBlock() + H_last.GetBlock();
    else if (BC == "OBC")
        H_block = H.GetBlock();
    uni10::EigH(H_block, En, state, uni10::INPLACE);
    uni10::Resize(state, 1, state.col(), uni10::INPLACE);

    /// create isometry VT = V_Transpose
    uni10::Bond leg_up(uni10::BD_IN, 1);
    uni10::Bond leg_L(uni10::BD_OUT, H1.GetBlock().row() );
    uni10::Bond leg_R(uni10::BD_OUT, H2.GetBlock().col() );
    vector<uni10::Bond> bondVT;
    bondVT.push_back(leg_up); bondVT.push_back(leg_L); bondVT.push_back(leg_R);
    uni10::UniTensor<double> VT(bondVT);
    VT.PutBlock(state);
    VT.SetLabel({1, -1, -3});

    /// uni10::permute to {up, R, L} (rotation 90) and save it
    uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
    VTs.push_back(VT);
    Vs_loc.push_back(0);

    /// return ground state energy save in queue Q
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );

    Q.clear();
    Q = energy;

    if (save_RG_info)
    {
        /// save last RG info
        double lowest_gap = energy[1] - energy[0];
        layer[0] = max(layer[0], layer[1]) + 1;
        layer.erase(layer.begin() + 1);
        RG_stage.push_back(layer[0]);
        RG_highest_gaps.push_back(lowest_gap);
        RG_lowest_gaps.push_back(lowest_gap);

        /// saving data into folder
        string file1, file2, folder;
        folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
        file1 = folder + "/RG_highest_gaps.csv";
        file2 = folder + "/RG_lowest_gaps.csv";

        ofstream fout1(file1);
        ofstream fout2(file2);
        fout1 << "layer,gap" << endl;
        fout2 << "layer,gap" << endl;

        for (int i = 0; i < RG_stage.size(); i++)
        {
            fout1 << setprecision(16) << RG_stage[i] << "," << RG_highest_gaps[i] << endl;
            fout2 << setprecision(16) << RG_stage[i] << "," << RG_lowest_gaps[i] << endl;
        }
    }
    
}