#include "tSDRG_tools.h"

/// Set Tensor's bond
uni10::UniTensor<double> SetBond(int v_in, int p_in, int v_out, int p_out, int D, int d)
{
    uni10::Bond virt_in(uni10::BD_IN, D);
    uni10::Bond phys_in(uni10::BD_IN, d);
    uni10::Bond virt_out(uni10::BD_OUT, D);
    uni10::Bond phys_out(uni10::BD_OUT, d);

    vector<uni10::Bond> bond;

    for (int i=0; i!=v_in; i++)
        bond.push_back(virt_in);

    for (int i=0; i!=p_in; i++)
        bond.push_back(phys_in);

    for (int i=0; i!=v_out; i++)
        bond.push_back(virt_out);

    for (int i=0; i!=p_out; i++)
        bond.push_back(phys_out);

    uni10::UniTensor<double> tensor(bond);

    return tensor;
}

/// return random number by distribution
double Distribution_Random_Variable(int model, double var, double Jdis)
{
    double value;
    switch (model)
    {
        case 0:
        {
            /// no-disorder, all is one
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
double RG_J(uni10::Matrix<double> En, const int chi, bool& info)
{
    vector<double> gap;
    for (int i=0; i<En.col()-1; i++)
    {
        gap.push_back(En.At(i+1, i+1) - En.At(i, i) );  // At(i, j)
        //cout << setprecision(8) << gap[i] << endl;
    }

    int idx;
    if (chi > gap.size() )
        idx = gap.size() - 1;
    else 
        idx = chi - 1;

    /// keep SU2 by skip mutiplet
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

/// class: tSDRG for clean system
void tSDRG0(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info)
{
    int L = MPO_chain.size();
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");

    // check Clean system
    bool clean;
    int idx_clean = 0;
    vector<int> clean_j;
    if ( adjacent_find( J_list.begin(), J_list.end(), not_equal_to<double>() ) == J_list.end() )
    {
        clean = true;
        if (L%4 == 0)
        {
            int temp  = 0;
            int layer = 1;
            int num_tensor = 0;
            for (int i=0; i<L-1; i++)
            {
                clean_j.push_back(temp - num_tensor);
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
        
    }
    else
    {
        ostringstream err;
        err << "This is NOT clean system.";
        throw runtime_error(err.str());
    }

    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j;                                                  // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        if (clean)
        {
            j = clean_j[idx_clean];
            idx_clean++;
        }
        else
        {
            j = distance(J_list.begin(), jmax);
        }

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        H1 = MPO_chain[j  ].GetTensor('l');
        H2 = MPO_chain[j+1].GetTensor('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H);

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 
            
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

        /// RG two site MPO
        RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);

        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + j);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);
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

    /// w_down = V
    uni10::UniTensor<double> V = VT;
    uni10::Permute(V, {-1, -3, 1}, 2, uni10::INPLACE);

    /// uni10::permute to {up, R, L} (rotation 90) and save it
    uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
    VTs.push_back(VT);
    Vs_loc.push_back(0);
    
    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );
    J_list.clear();
    J_list = energy;
}

/// class: tSDRG 
void tSDRG(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info)
{
    int L = MPO_chain.size();
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        for (int i = 0; i < J_list.size(); i++)
        {
            J_list[i] = abs(J_list[i]);
        }
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j = distance(J_list.begin(), jmax);;                // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        H1 = MPO_chain[j  ].GetTensor('l');
        H2 = MPO_chain[j+1].GetTensor('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H);

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

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

        /// RG two site MPO
        RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);

        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + j);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (j != 0)
        {
            H1 = MPO_chain[j-1].GetTensor('l');
            H2 = MPO_chain[j  ].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);   // in mkl use uni10::EigHLazy()
            J_list[j-1] = RG_J(En, chi, info);
        }

        if (j != J_list.size()-1)
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);  // in mkl use uni10::EigHLazy()
            J_list[j+1] = RG_J(En, chi, info);
        }

        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);
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

    /// w_down = V
    uni10::UniTensor<double> V = VT;
    uni10::Permute(V, {-1, -3, 1}, 2, uni10::INPLACE);

    /// uni10::permute to {up, R, L} (rotation 90) and save it
    uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
    VTs.push_back(VT);
    Vs_loc.push_back(0);
    
    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );
    J_list.clear();
    J_list = energy;
}

/// tSDRG2 New criterion for select
void tSDRG2(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info)
{
    int L = MPO_chain.size();
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");

    /// energy of single site
    vector<double> enss;
    for (int i=0; i<L; i++)
    {
        uni10::UniTensor<double> H;
        H = MPO_chain[i].GetTensorSS();

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);

        enss.push_back(En.At(0, 0) );
    }

    vector<double> en0;
    for (int i=0; i<L-1; i++)
    {
        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        H1 = MPO_chain[i  ].GetTensor('l');
        H2 = MPO_chain[i+1].GetTensor('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H);

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);

        en0.push_back(enss[i] + enss[i+1] - En.At(0, 0) );
    }

    /// rho
    vector<double> rho_list;
    for (int i=0; i<L-1; i++)
    {
        if (i == 0)
        {
            rho_list.push_back(en0[i]/en0[i+1]);
        }
        else if (i == L-2)
        {
            rho_list.push_back(en0[i]/en0[i-1]);
        }
        else
        {
            double temp = max(en0[i-1], en0[i+1]);
            rho_list.push_back(en0[i]/temp);
        }
    }

    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto rhomax = max_element(rho_list.begin(), rho_list.end() ); // max_element return iterators, not values. 
        int j = distance(rho_list.begin(), rhomax);                   // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                                  // chi of location (chi maybe cut at mutliplet, so chi will change)

        //for (int i=0; i<rho_list.size(); i++)
            //cout << setprecision(16) << rho_list[i] << endl;
        
        //cout << setprecision(8) << "Jmax = " << *rhomax << ", r = " << j << endl;
        //cout << "*** mpo size = " << MPO_chain.size() << " ***" << endl;

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        H1 = MPO_chain[j  ].GetTensor('l');
        H2 = MPO_chain[j+1].GetTensor('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H);
        
        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 
            
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

        /// RG two site MPO
        RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);

        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == rho_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + j);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// erase energy of single site j+1
        enss.erase(enss.begin() + j + 1);

        /// update enss
        H = MPO_chain[j].GetTensorSS();
        uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
        enss[j] = En.At(0, 0);

        /// RG coupling into merge list, notice eage mpo
        if (j != 0)
        {
            H1 = MPO_chain[j-1].GetTensor('l');
            H2 = MPO_chain[j  ].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);

            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
            en0[j-1] = enss[j-1] + enss[j] - En.At(0, 0);
        }

        if (j != rho_list.size()-1)
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);

            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
            en0[j+1] = enss[j] + enss[j+1] - En.At(0, 0);
        }

        /// update rho list
        double temp;
        if (j == 0)
        {
            if (rho_list.size() > 3)
            {
                /// update j+1
                temp = en0[j+2];
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = max(en0[j+1], en0[j+3]);
                rho_list[j+2] = en0[j+2]/temp;
            }
            else
            {
                /// update j+1
                temp = en0[j+2];
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = en0[j+1];
                rho_list[j+2] = en0[j+2]/temp;
            }
        }
        else if (j == 1)
        {
            if (rho_list.size() > 4)
            {
                /// update j-1
                temp = en0[j+1];
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = max(en0[j-1], en0[j+2]);
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = max(en0[j+1], en0[j+3]);
                rho_list[j+2] = en0[j+2]/temp;
            }
            else
            {
                /// update j-1
                temp = en0[j+1];
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = max(en0[j-1], en0[j+2]);
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = en0[j+1];
                rho_list[j+2] = en0[j+2]/temp;
            }
        }
        else if (j == rho_list.size()-2)
        {
            if (rho_list.size() > 4)
            {
                /// update j-2
                temp = max(en0[j-3], en0[j-1]);
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = max(en0[j-2], en0[j+1]);
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = en0[j-1];
                rho_list[j+1] = en0[j+1]/temp;
            }
            else
            {
                /// update j-2
                temp = en0[j-1];
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = max(en0[j-2], en0[j+1]);
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = en0[j-1];
                rho_list[j+1] = en0[j+1]/temp;
            }
        }
        else if (j == rho_list.size()-1)
        {
            if (rho_list.size() > 3)
            {
                /// update j-2
                temp = max(en0[j-3], en0[j-1]);
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = en0[j-2];
                rho_list[j-1] = en0[j-1]/temp;
            }
            else
            {
                /// update j-2
                temp = en0[j-1];
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = en0[j-2];
                rho_list[j-1] = en0[j-1]/temp;
            }
        }
        else
        {
            if (rho_list.size() == 5)
            { 
                /// update j-2
                temp = en0[j-1];
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = max(en0[j-2], en0[j+1]);
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = max(en0[j-1], en0[j+2]);
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = en0[j+1];
                rho_list[j+2] = en0[j+2]/temp;
            }
            else if (j == 2)
            {
                /// update j-2
                temp = en0[j-1];
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = max(en0[j-2], en0[j+1]);
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = max(en0[j-1], en0[j+2]);
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = max(en0[j+1], en0[j+3]);
                rho_list[j+2] = en0[j+2]/temp;
            }
            else if (j == rho_list.size()-3)
            {
                /// update j-2
                temp = max(en0[j-3], en0[j-1]);
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = max(en0[j-2], en0[j+1]);
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = max(en0[j-1], en0[j+2]);
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = en0[j+1];
                rho_list[j+2] = en0[j+2]/temp;
            }
            else
            {
                /// update j-2
                temp = max(en0[j-3], en0[j-1]);
                rho_list[j-2] = en0[j-2]/temp;
            
                /// update j-1
                temp = max(en0[j-2], en0[j+1]);
                rho_list[j-1] = en0[j-1]/temp;

                /// update j+1
                temp = max(en0[j-1], en0[j+2]);
                rho_list[j+1] = en0[j+1]/temp;

                /// update j+2
                temp = max(en0[j+1], en0[j+3]);
                rho_list[j+2] = en0[j+2]/temp;
            }
        }

        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        en0.erase(en0.begin() + j);
        rho_list.erase(rho_list.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);
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

    /// w_down = V
    uni10::UniTensor<double> V = VT;
    uni10::Permute(V, {-1, -3, 1}, 2, uni10::INPLACE);

    /// uni10::permute to {up, R, L} (rotation 90) and save it
    uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
    VTs.push_back(VT);
    Vs_loc.push_back(0);
    
    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );

    J_list.clear();
    J_list = energy;
}

/// class: tSDRG 
void tSDRG_PBC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info)
{
    int L = MPO_chain.size();
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");

    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j = distance(J_list.begin(), jmax);                 // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        /*cout << "Init MPO = " << MPO_chain.size() << endl;
        for (int i=0; i<J_list.size(); i++)
            cout << setprecision(16) << J_list[i] << endl;*/

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = RG_J(En, chi, info);
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );
    J_list.clear();
    J_list = energy;
}

/// tSDRG2 New criterion for select
void tSDRG2_PBC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, bool& info)
{
    int L = MPO_chain.size();
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");

    /// energy of single site
    vector<double> enss;
    for (int i=0; i<L; i++)
    {
        uni10::UniTensor<double> H;
        H = MPO_chain[i].GetTensorSS();

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);

        enss.push_back(En.At(0, 0) );
    }

    vector<double> en0;
    for (int i=0; i<L; i++)
    {
        /// PBC postion
        int loc1, loc2;    
        loc1 = i;
        if (i == L-1)
            loc2 = 0;
        else    
            loc2 = i+1;

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        H1 = MPO_chain[loc1].GetTensor('l');
        H2 = MPO_chain[loc2].GetTensor('r');
        H12.PutTensor("H1", H1);
        H12.PutTensor("H2", H2);
        H12.Launch(H);

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);

        en0.push_back(enss[i] + enss[i+1] - En.At(0, 0) );
    }

    /// rho
    vector<double> rho_list;
    for (int i=0; i<L; i++)
    {
        /// PBC postion
        int locL, locM, locR;    
        locM = i;
        if (i == 0)
        {
            locL = L-1;
            locR = i+1;
        }
        else if (i == L-1)
        {
            locL = i-1;
            locR = 0;
        }
        else    
        {
            locL = i-1;
            locR = i+1;
        }

        double temp = max(en0[locL], en0[locR]);
        rho_list.push_back(en0[locM]/temp);
    }

    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto rhomax = max_element(rho_list.begin(), rho_list.end() ); // max_element return iterators, not values. 
        int j = distance(rho_list.begin(), rhomax);                   // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                                  // chi of location (chi maybe cut at mutliplet, so chi will change)

        //for (int i=0; i<rho_list.size(); i++)
            //cout << setprecision(16) << rho_list[i] << endl;
        
        //cout << setprecision(8) << "Jmax = " << *rhomax << ", r = " << j << endl;
        //cout << "*** mpo size = " << MPO_chain.size() << " ***" << endl;

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == rho_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }
        
        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 
        
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
        
        /// RG two site MPO
        if (j == rho_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);
        
        /// erase j+1, except j+1 is most left mpo and update enss
        if (j == rho_list.size()-1)
        {
            MPO_chain.erase(MPO_chain.begin() + 0);
            enss.erase(enss.begin() + 0);

            H = MPO_chain[j-1].GetTensorSS();
            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
            enss[j-1] = En.At(0, 0);
        }
        else
        {
            MPO_chain.erase(MPO_chain.begin() + j + 1);
            enss.erase(enss.begin() + j + 1);

            H = MPO_chain[j].GetTensorSS();
            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
            enss[j] = En.At(0, 0);
        }
        
        /// RG coupling into merge list, notice eage mpo
        if (j != 0)
        {
            /// PBC postion
            int loc1, loc2;
            if (j != rho_list.size()-1)   
            {
                loc1 = j-1;
                loc2 = j;
            }   
            else
            {
                loc1 = j-2;
                loc2 = j-1;
            }

            H1 = MPO_chain[loc1].GetTensor('l');
            H2 = MPO_chain[loc2].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);

            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);

            if (j != rho_list.size()-1)
                en0[loc1] = enss[loc1] + enss[loc2] - En.At(0, 0);
            else
                en0[loc2] = enss[loc1] + enss[loc2] - En.At(0, 0);
        }
        else
        {
            /// PBC postion
            int loc1, loc2;    
            loc1 = MPO_chain.size()-1;
            loc2 = 0;

            H1 = MPO_chain[loc1].GetTensor('l');
            H2 = MPO_chain[loc2].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);

            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
            en0[en0.size()-1] = enss[loc1] + enss[loc2] - En.At(0, 0);
        }

        if (j != rho_list.size()-1)
        {
            /// PBC postion
            int loc1, loc2;
            if (j != rho_list.size()-2)   
            {
                loc1 = j;
                loc2 = j+1;
            }   
            else
            {
                loc1 = j;
                loc2 = 0;
            }

            H1 = MPO_chain[loc1].GetTensor('l');
            H2 = MPO_chain[loc2].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);

            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);

            if (j != rho_list.size()-2)
                en0[loc2] = enss[loc1] + enss[loc2] - En.At(0, 0);
            else
                en0[en0.size()-1] = enss[loc1] + enss[loc2] - En.At(0, 0);
        }
        else
        {
            /// PBC postion
            int loc1, loc2;    
            loc1 = MPO_chain.size()-1;
            loc2 = 0;

            H1 = MPO_chain[loc1].GetTensor('l');
            H2 = MPO_chain[loc2].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);

            uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
            en0[0] = enss[loc1] + enss[loc2] - En.At(0, 0);
        }

        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        en0.erase(en0.begin() + j);
        rho_list.erase(rho_list.begin() + j);

        /// update rho list
        for (int i=0; i<rho_list.size(); i++)
        {
            /// PBC postion
            int locL, locM, locR;    
            locM = i;
            if (i == 0)
            {
                locL = rho_list.size()-1;
                locR = i+1;
            }
            else if (i == rho_list.size()-1)
            {
                locL = i-1;
                locR = 0;
                //cout << "TEST " << locL << endl;
            }
            else    
            {
                locL = i-1;
                locR = i+1;
            }

            double temp = max(en0[locL], en0[locR]);
            rho_list[i] = en0[locM]/temp;
        }

    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// w_down = V
    uni10::UniTensor<double> V = VT;
    uni10::Permute(V, {-1, -3, 1}, 2, uni10::INPLACE);

    /// uni10::permute to {up, R, L} (rotation 90) and save it
    uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
    VTs.push_back(VT);
    Vs_loc.push_back(0);
    
    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );

    J_list.clear();
    J_list = energy;
}

/// class: tSDRG 
void tSDRG_PBC_layer(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info)
{
    int L = MPO_chain.size();
    vector<int> layer;
    layer.assign(L, 0);

    string file, folder;
    folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    file = folder + "/layer.csv";
    ofstream fout(file);
    fout << "layer,gap" << endl;
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j = distance(J_list.begin(), jmax);                 // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        //cout << "Init MPO = " << MPO_chain.size() << endl;
        //for (int i=0; i<J_list.size(); i++)
            //cout << setprecision(16) << J_list[i] << endl;
        
        
        if (j == J_list.size()-1)
        {
            layer[j] = max(layer[j], layer[0]) + 1;
            layer.erase(layer.begin() + 0);
        }
        else
        {
            layer[j] = max(layer[j], layer[j+1]) + 1;
            layer.erase(layer.begin() + j + 1);
        }
            
        fout << setprecision(16) << layer[j] << "," << *jmax << endl;

        //gap_value.push_back(J_list[i]);
        //layer.push_back(gap_value);

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = RG_J(En, chi, info);
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );
    J_list.clear();
    J_list = energy;

    layer[0] = max(layer[0], layer[1]) + 1;
    layer.erase(layer.begin() + 1);
    fout << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    
}

/// class: tSDRG 
void tSDRG_PBC_layerlowest(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info)
{
    int L = MPO_chain.size();
    vector<int> layer;
    layer.assign(L, 0);
    vector<double> J_list_lowest = J_list;

    string file, folder;
    folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    file = folder + "/layerlowest.csv";
    ofstream fout(file);
    fout << "layer,gap" << endl;
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j = distance(J_list.begin(), jmax);                 // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        //cout << "Init MPO = " << MPO_chain.size() << endl;
        //for (int i=0; i<J_list.size(); i++)
            //cout << setprecision(16) << J_list[i] << endl;
        
        
        if (j == J_list.size()-1)
        {
            layer[j] = max(layer[j], layer[0]) + 1;
            layer.erase(layer.begin() + 0);
        }
        else
        {
            layer[j] = max(layer[j], layer[j+1]) + 1;
            layer.erase(layer.begin() + j + 1);
        }
            
        
        fout << setprecision(16) << layer[j] << "," << J_list_lowest[j] << endl;

        //gap_value.push_back(J_list[i]);
        //layer.push_back(gap_value);

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = RG_J(En, chi, info);
                J_list_lowest[J_list.size()-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = RG_J(En, chi, info);
                J_list_lowest[0] = En.At(1, 1) - En.At(0, 0);
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        J_list_lowest.erase(J_list_lowest.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );
    J_list.clear();
    J_list = energy;

    layer[0] = max(layer[0], layer[1]) + 1;
    layer.erase(layer.begin() + 1);
    fout << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    
}

/// class: tSDRG 
void tSDRG_PBC_layer2(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info)
{
    int L = MPO_chain.size();

    vector<int> layer;
    layer.assign(L, 0);

    vector<double> J_list_lowest = J_list;

    string file1, file2, folder;
    folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    file1 = folder + "/layer.csv";
    file2 = folder + "/layerlowest.csv";
    ofstream fout1(file1);
    ofstream fout2(file2);
    fout1 << "layer,gap" << endl;
    fout2 << "layer,gap" << endl;
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j = distance(J_list.begin(), jmax);                 // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        //cout << "Init MPO = " << MPO_chain.size() << endl;
        //for (int i=0; i<J_list.size(); i++)
            //cout << setprecision(16) << J_list[i] << endl;
        
        
        if (j == J_list.size()-1)
        {
            layer[j] = max(layer[j], layer[0]) + 1;
            layer.erase(layer.begin() + 0);
        }
        else
        {
            layer[j] = max(layer[j], layer[j+1]) + 1;
            layer.erase(layer.begin() + j + 1);
        }
            
        fout1 << setprecision(16) << layer[j] << "," << *jmax << endl;
        fout2 << setprecision(16) << layer[j] << "," << J_list_lowest[j] << endl;

        //gap_value.push_back(J_list[i]);
        //layer.push_back(gap_value);

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = RG_J(En, chi, info);
                J_list_lowest[J_list.size()-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = RG_J(En, chi, info);
                J_list_lowest[0] = En.At(1, 1) - En.At(0, 0);
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        J_list_lowest.erase(J_list_lowest.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );
    J_list.clear();
    J_list = energy;

    layer[0] = max(layer[0], layer[1]) + 1;
    layer.erase(layer.begin() + 1);
    fout1 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    fout2 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    
}

/// class: tSDRG 
void tSDRG_PBC_layerBig(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info)
{
    int L = MPO_chain.size();

    vector<int> layer;
    layer.assign(L, 0);

    vector<double> J_list_lowest = J_list;

    string file1, file2, folder;
    folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    file1 = folder + "/layerBig.csv";
    file2 = folder + "/layerlowestBig.csv";
    ofstream fout1(file1);
    ofstream fout2(file2);
    fout1 << "layer,gap" << endl;
    fout2 << "layer,gap" << endl;
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() );    // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        double jmin = *min_element(J_list.begin(), J_list.end() ); // *max_element return double
        int j = distance(J_list.begin(), jmax);                    // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                               // chi of location (chi maybe cut at mutliplet, so chi will change)


        //cout << "Init MPO = " << MPO_chain.size() << endl;
        //for (int i=0; i<J_list.size(); i++)
            //cout << setprecision(16) << J_list[i] << endl;
        
        
        if (j == J_list.size()-1)
        {
            layer[j] = max(layer[j], layer[0]) + 1;
            layer.erase(layer.begin() + 0);
        }
        else
        {
            layer[j] = max(layer[j], layer[j+1]) + 1;
            layer.erase(layer.begin() + j + 1);
        }
            
        fout1 << setprecision(16) << layer[j] << "," << *jmax << endl;
        fout2 << setprecision(16) << layer[j] << "," << J_list_lowest[j] << endl;

        //gap_value.push_back(J_list[i]);
        //layer.push_back(gap_value);

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        H = (1.0/jmin) * H;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = jmin * RG_J(En, chi, info);
                J_list_lowest[J_list.size()-1] = jmin * (En.At(1, 1) - En.At(0, 0));

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = jmin * RG_J(En, chi, info);
                J_list_lowest[j+1] = jmin * (En.At(1, 1) - En.At(0, 0));
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = jmin * RG_J(En, chi, info);
                J_list_lowest[j-1] = jmin * (En.At(1, 1) - En.At(0, 0));

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE);
                J_list[j+1] = jmin * RG_J(En, chi, info);
                J_list_lowest[j+1] = jmin * (En.At(1, 1) - En.At(0, 0));
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = jmin * RG_J(En, chi, info);
                J_list_lowest[j-1] = jmin * (En.At(1, 1) - En.At(0, 0));

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = jmin * RG_J(En, chi, info);
                J_list_lowest[0] = jmin * (En.At(1, 1) - En.At(0, 0));
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = jmin * RG_J(En, chi, info);
                J_list_lowest[j-1] = jmin * (En.At(1, 1) - En.At(0, 0));

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                H = (1.0/jmin) * H;
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = jmin * RG_J(En, chi, info);
                J_list_lowest[j+1] = jmin * (En.At(1, 1) - En.At(0, 0));
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        J_list_lowest.erase(J_list_lowest.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H       = (100.0) * H;
    H_last  = (100.0) * H_last;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
    {
        energy.push_back( (1.0/100.0) * En.At(i, i) );
    }  
    J_list.clear();
    J_list = energy;

    layer[0] = max(layer[0], layer[1]) + 1;
    layer.erase(layer.begin() + 1);
    fout1 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    fout2 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    
}

/// class: tSDRG 
void tSDRG_PBC_layerBigC(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info)
{
    int L = MPO_chain.size();
    double coeff = 1.0;

    vector<int> layer;
    layer.assign(L, 0);

    // transfor coupling to energy scale (highest gap); S = 1 highest gap = 2; S = 0.5, highest gap = 1 
    for (int i=0; i<J_list.size(); i++)
        J_list[i] = 2 * J_list[i];     

    vector<double> J_list_lowest = J_list;

    string file1, file2, folder;
    folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    file1 = folder + "/layerBigC.csv";
    file2 = folder + "/layerlowestBigC.csv";
    ofstream fout1(file1);
    ofstream fout2(file2);
    fout1 << "layer,gap" << endl;
    fout2 << "layer,gap" << endl;
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j = distance(J_list.begin(), jmax);                 // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        //cout << "Init MPO = " << MPO_chain.size() << endl;
        //for (int i=0; i<J_list.size(); i++)
            //cout << setprecision(16) << J_list[i] << endl;
        
        if (j == J_list.size()-1)
        {
            layer[j] = max(layer[j], layer[0]) + 1;
            layer.erase(layer.begin() + 0);
        }
        else
        {
            layer[j] = max(layer[j], layer[j+1]) + 1;
            layer.erase(layer.begin() + j + 1);
        }
            
        fout1 << setprecision(16) << layer[j] << "," << *jmax << endl;
        fout2 << setprecision(16) << layer[j] << "," << J_list_lowest[j] << endl;

        //gap_value.push_back(J_list[i]);
        //layer.push_back(gap_value);

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        //H = coeff * H;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = RG_J(En, chi, info);
                J_list_lowest[J_list.size()-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = RG_J(En, chi, info);
                J_list_lowest[0] = En.At(1, 1) - En.At(0, 0);
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        J_list_lowest.erase(J_list_lowest.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    H       = (100.0) * H;
    H_last  = (100.0) * H_last;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
    {
        energy.push_back( (1.0/100.0) * En.At(i, i) );
    }
    J_list.clear();
    J_list = energy;

    layer[0] = max(layer[0], layer[1]) + 1;
    layer.erase(layer.begin() + 1);
    fout1 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    fout2 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    
}

/// class: tSDRG 
void tSDRG_PBC_layergg(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info)
{
    int L = MPO_chain.size();
    double coeff = 1.0;

    vector<int> layer;
    layer.assign(L, 0);

    // transfor coupling to energy scale (highest gap); S = 1 highest gap = 2; S = 0.5, highest gap = 1 
    for (int i=0; i<J_list.size(); i++)
    {
        J_list[i] = 2 * J_list[i]; 
        cout << J_list[i] << endl;
    }

    string file1, file2, folder;
    folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    file1 = folder + "/layergg.csv";
    file2 = folder + "/layerlowestgg.csv";
    ofstream fout1(file1);
    ofstream fout2(file2);
    fout1 << "layer,gap" << endl;
    fout2 << "layer,gap" << endl;
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j = distance(J_list.begin(), jmax);                 // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        //cout << "Init MPO = " << MPO_chain.size() << endl;
        //for (int i=0; i<J_list.size(); i++)
            //cout << setprecision(16) << J_list[i] << endl;
        
        if (j == J_list.size()-1)
        {
            layer[j] = max(layer[j], layer[0]) + 1;
            layer.erase(layer.begin() + 0);
        }
        else
        {
            layer[j] = max(layer[j], layer[j+1]) + 1;
            layer.erase(layer.begin() + j + 1);
        }
            
        fout1 << setprecision(16) << layer[j] << "," << *jmax << endl;

        //gap_value.push_back(J_list[i]);
        //layer.push_back(gap_value);

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        //H = coeff * H;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);


        /// return ground state energy save in J list
        double gap_lowest = 0.0;
        gap_lowest = abs(En.At(1, 1) - En.At(0, 0));
        for (int n=0; n<En.col(); n++) 
            cout << En.At(n, n) << endl;
        cout << "************" << endl;
        fout2 << setprecision(16) << layer[j] << "," << gap_lowest << endl;

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = RG_J(En, chi, info);
                //J_list_lowest[J_list.size()-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                //J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                //J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                //J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                //J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = RG_J(En, chi, info);
                //J_list_lowest[0] = En.At(1, 1) - En.At(0, 0);
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                //J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                //J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        //J_list_lowest.erase(J_list_lowest.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    //H       = (100.0) * H;
    //H_last  = (100.0) * H_last;
    cout << H.GetBlock() << endl;
    cout << H_last.GetBlock() << endl;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
    {
        energy.push_back( En.At(i, i) );
        //cout << energy[i] << endl;
    }
    J_list.clear();
    J_list = energy;

    layer[0] = max(layer[0], layer[1]) + 1;
    layer.erase(layer.begin() + 1);
    fout1 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    fout2 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    
}

/// class: tSDRG 
void tSDRG0_PBC_test(vector<MPO>& MPO_chain, vector<double>& J_list, vector<uni10::UniTensor<double> >& VTs, vector<int>& Vs_loc, const int chi, string dis, const int Pdis, const int Jseed, bool& info)
{
    int L = MPO_chain.size();
    double coeff = 1.0;

    // check Clean system
    bool clean;
    int idx_clean = 0;
    vector<int> clean_j;
    clean = true;
    if (L%4 == 0)
    {
        int temp  = 0;
        int layer = 1;
        int num_tensor = 0;
        for (int i=0; i<L-1; i++)
        {
            clean_j.push_back(temp - num_tensor);
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

    vector<int> layer;
    layer.assign(L, 0);

    // transfor coupling to energy scale (highest gap); S = 1 highest gap = 2; S = 0.5, highest gap = 1 
    for (int i=0; i<J_list.size(); i++)
    {
        J_list[i] = 2 * J_list[i]; 
        cout << J_list[i] << endl;
    }

    string file1, file2, folder;
    folder = "TTN/algo/Jdis" + dis + "/L" + to_string(L) + "_P" + to_string(Pdis) + "_m" + to_string(chi) + "_" + to_string(Jseed);
    file1 = folder + "/layergg.csv";
    file2 = folder + "/layerlowestgg.csv";
    ofstream fout1(file1);
    ofstream fout2(file2);
    fout1 << "layer,gap" << endl;
    fout2 << "layer,gap" << endl;
    
    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");
    
    while(MPO_chain.size() > 2)
    {
        /// find max J
        auto jmax = max_element(J_list.begin(), J_list.end() ); // max_element return iterators, not values. ( double Jmax = *jmax; use *iterators to get values.
        int j;                                                  // return distance from 0 to jmax. Note: distance(iterators, iterators) NOT value.
        int chi_loc;                                            // chi of location (chi maybe cut at mutliplet, so chi will change)
        
        //cout << "Init MPO = " << MPO_chain.size() << endl;
        //for (int i=0; i<J_list.size(); i++)
            //cout << setprecision(16) << J_list[i] << endl;

        j = clean_j[idx_clean];
        idx_clean++;
        
        if (j == J_list.size()-1)
        {
            layer[j] = max(layer[j], layer[0]) + 1;
            layer.erase(layer.begin() + 0);
        }
        else
        {
            layer[j] = max(layer[j], layer[j+1]) + 1;
            layer.erase(layer.begin() + j + 1);
        }
            
        fout1 << setprecision(16) << layer[j] << "," << *jmax << endl;

        //gap_value.push_back(J_list[i]);
        //layer.push_back(gap_value);

        /// Contract two site hamitonian
        uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
        if (j == J_list.size()-1)
        {
            H1 = MPO_chain[j].GetTensor('l');
            H2 = MPO_chain[0].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);   
        }
        else
        {
            H1 = MPO_chain[j  ].GetTensor('l');
            H2 = MPO_chain[j+1].GetTensor('r');
            H12.PutTensor("H1", H1);
            H12.PutTensor("H2", H2);
            H12.Launch(H);
        }

        /// diagonal
        uni10::Matrix<double> En;                               // eigen energy
        uni10::Matrix<double> state;                            // eigen state
        uni10::Matrix<double> H_block;
        //H = coeff * H;
        H_block = H.GetBlock();
        uni10::EigH(H_block, En, state, uni10::INPLACE);


        /// return ground state energy save in J list
        double gap_lowest = 0.0;
        gap_lowest = abs(En.At(1, 1) - En.At(0, 0));
        for (int n=0; n<En.col(); n++) 
            cout << En.At(n, n) << endl;
        //cout << "************" << endl;
        fout2 << setprecision(16) << layer[j] << "," << gap_lowest << endl;

        /// truncation
        chi_loc = H.GetBlock().col();
        if (chi_loc > chi)
            Truncation(En, state, chi, chi_loc, info); 

        if (info == 0)
            return;
        //cout << chi_loc << H1.GetBlock().row() << H2.GetBlock().col() << endl;
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

        /// RG two site MPO
        if (j == J_list.size()-1)
            RG_MPO(MPO_chain[j], MPO_chain[0], VT, V);
        else
            RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);
        
        /// uni10::permute to {up, R, L} (rotation 90 in order to return normal form) and save it
        uni10::Permute(VT, {1, -3, -1}, 1, uni10::INPLACE);
        VTs.push_back(VT);
        Vs_loc.push_back(j);

        /// erase j+1, except j+1 is most left mpo 
        if (j == J_list.size()-1)
            MPO_chain.erase(MPO_chain.begin() + 0);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

        /// RG coupling into merge list, notice eage mpo
        if (MPO_chain.size() != 2)
        {
            if (j == 0)
            {
                H1 = MPO_chain[J_list.size()-2].GetTensor('l');
                H2 = MPO_chain[j              ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[J_list.size()-1] = RG_J(En, chi, info);
                //J_list_lowest[J_list.size()-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                //J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-2)
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                //J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                //J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
            else if (j == J_list.size()-1)
            {
                H1 = MPO_chain[j-2].GetTensor('l');
                H2 = MPO_chain[j-1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                //J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[0].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[0] = RG_J(En, chi, info);
                //J_list_lowest[0] = En.At(1, 1) - En.At(0, 0);
            }
            else
            {
                H1 = MPO_chain[j-1].GetTensor('l');
                H2 = MPO_chain[j  ].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j-1] = RG_J(En, chi, info);
                //J_list_lowest[j-1] = En.At(1, 1) - En.At(0, 0);

                H1 = MPO_chain[j  ].GetTensor('l');
                H2 = MPO_chain[j+1].GetTensor('r');
                H12.PutTensor("H1", H1);
                H12.PutTensor("H2", H2);
                H12.Launch(H);
                uni10::EigHLazy(H.GetBlock(), En, state, uni10::INPLACE); 
                J_list[j+1] = RG_J(En, chi, info);
                //J_list_lowest[j+1] = En.At(1, 1) - En.At(0, 0);
            }
        }
        
        if (info == 0)
            return;
        
        /// erase J list at j (merge coupling)
        J_list.erase(J_list.begin() + j);
        //J_list_lowest.erase(J_list_lowest.begin() + j);
    }

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H_last, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)

    /// 1 and 0
    H1 = MPO_chain[1].GetTensorPBC('l');
    H2 = MPO_chain[0].GetTensorPBC('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H_last);
    //H_last.PrintDiagram();
    uni10::Permute(H_last, {-3, -1, -4, -2}, 2, uni10::INPLACE);

    /// 0 and 1 , use H1 and H2 to find isometry two legs of origin picture
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);
    
    /// diagonal
    uni10::Matrix<double> En;                               // eigen energy
    uni10::Matrix<double> state;                            // eigen state
    uni10::Matrix<double> H_block;
    //H       = (100.0) * H;
    //H_last  = (100.0) * H_last;
    //cout << H.GetBlock() << endl;
    //cout << H_last.GetBlock() << endl;
    H_block = H.GetBlock() + H_last.GetBlock();
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

    /// return ground state energy save in J list
    vector<double> energy;
    for (int i=0; i<En.col(); i++)
    {
        energy.push_back( En.At(i, i) );
        //cout << energy[i] << endl;
    }
    J_list.clear();
    J_list = energy;

    layer[0] = max(layer[0], layer[1]) + 1;
    layer.erase(layer.begin() + 1);
    fout1 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    fout2 << setprecision(16) << layer[0] << "," << energy[1] - energy[0] << endl;
    
}