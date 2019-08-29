#include "measure.h"

/// return location of top-tensor Decision tree
string Decision_tree(vector<int> Vs_loc, bool show)
{
    map<int, string> tree;
    int L = Vs_loc.size() + 1;

    vector<string> mpo;
    for(int i=0; i<L; i++)
        mpo.push_back(to_string(i) );

    string t1, t2, order, right_block;
    int loc1, loc2;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpo.size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;

        t1 = mpo[loc1];
        t2 = mpo[loc2];

        if (i != Vs_loc.size()-1)
            order = "(" + t1 + ", " + t2 + ")";
        else
        {
            order = "{" + t1 + " | " + t2 + "}";
            right_block = t2;
        }
        tree[i] = order;

        mpo[loc1] = order;
        mpo.erase(mpo.begin() + loc2);
    }

    if (show)
    {
        auto it = tree.begin();
        while(it != tree.end() )
        {
            cout << it->second << endl;
            it++;
        }
    }

    char chars[] = "(),";
    for (unsigned int i = 0; i < strlen(chars); ++i)
    {
        right_block.erase(remove(right_block.begin(), right_block.end(), chars[i] ), right_block.end() );
    }

    istringstream iss(right_block);
    string top;
    iss >> top;
    
    return top;
}

/// return Schmidt value of cut top-tensor 
vector<double> Schmidt_Value(uni10::UniTensor<double> top)
{
    uni10::Permute(top, top.label(), 2, uni10::INPLACE);
    vector<uni10::Matrix<double>> top_svd = uni10::Svd(top.GetBlock());
    int s = top_svd[1].col();

    vector<double> value;
    for (int i=0; i<s; i++)
        value.push_back(top_svd[1].At(i, i) );
    
    return value;
}

/// return ground state energy spectrum
vector<double> Energy_Spectrum(vector<MPO>& MPO_chain, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    int L = MPO_chain.size();

    int j;
    for (int i=0; i<Vs_loc.size()-1; i++)
    {
        int j = Vs_loc[i];

        uni10::UniTensor<double> VT, V;
        VT = VTs[i];
        V = Vs[i];
        uni10::Permute(VT, {1, -1, -3}, 1, uni10::INPLACE);
        uni10::Permute(V, {-1, -3, 1}, 2, uni10::INPLACE);

        /// RG two site MPO
        RG_MPO(MPO_chain[j], MPO_chain[j+1], VT, V);

        /// erase j+1, except j+1 is most left mpo 
        if (j == MPO_chain.size()-2)
            MPO_chain.erase(MPO_chain.begin() + j);
        else
            MPO_chain.erase(MPO_chain.begin() + j + 1);

    }

    /// Network
    uni10::Network H12("../tSDRG_net/H12.net");

    /// Contract last two site hamitonian
    uni10::UniTensor<double> H, H1, H2;                     // H is hamitonian of site1(= H1) and site2(= H2)
    H1 = MPO_chain[0].GetTensor('l');
    H2 = MPO_chain[1].GetTensor('r');
    H12.PutTensor("H1", H1);
    H12.PutTensor("H2", H2);
    H12.Launch(H);

    /// diagonal
    uni10::Matrix<double> En;                        
    uni10::Matrix<double> state;                    
    uni10::EigH(H.GetBlock(), En, state, uni10::INPLACE);

    vector<double> energy;
    for (int i=0; i<En.col(); i++)
        energy.push_back(En.At(i, i) );

    return energy;
}

/// Sz-Sz Correlation
double Correlation_SzSz(int site1, int site2, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double corr = 0;
    int L = Vs_loc.size() + 1;
    double spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> kara;
    Sz.PutBlock(matSz(spin));

    vector<uni10::UniTensor<double> > mpo;
    mpo.assign(L, kara);
    
    mpo[site1] = Sz;
    mpo[site2] = Sz;
    
    int loc1, loc2;
    uni10::UniTensor<double> temp, tempL, tempR;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpo.size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;
        
        if (mpo[loc1].BondNum() != 0 || mpo[loc2].BondNum() != 0)
        {
            if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() != 0)
            {
                tempL = mpo[loc1];
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1});
                tempL.SetLabel({-1, -3});
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() == 0)
            {
                tempL = mpo[loc1];

                VTs[i].SetLabel({1, -2, -1}); 
                tempL.SetLabel({-1, -3});
                Vs[i].SetLabel({-2, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() == 0 && mpo[loc2].BondNum() != 0)
            {
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1}); 
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -1, 2});

                temp = uni10::Contract(VTs[i], tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else
            {
                ostringstream err;
                err << "Error: SzSz Correlation fail contraction.";
                throw runtime_error(err.str());
            }

            mpo[loc1] = temp;
            mpo.erase(mpo.begin() + loc2);
        }
        else
        {
            mpo.erase(mpo.begin() + loc2);
        }
    }

    corr = temp[0];
    return corr;
}

/// Sx-Sx Correlation
double Correlation_SxSx(int site1, int site2, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double corr = 0;
    int L = Vs_loc.size() + 1;
    double spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sx(bonds);
    uni10::UniTensor<double> kara;
    Sx.PutBlock(matSx(spin));

    vector<uni10::UniTensor<double> > mpo;
    mpo.assign(L, kara);
    
    mpo[site1] = Sx;
    mpo[site2] = Sx;

    int loc1, loc2;
    uni10::UniTensor<double> temp, tempL, tempR;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpo.size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;
        
        if (mpo[loc1].BondNum() != 0 || mpo[loc2].BondNum() != 0)
        {
            if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() != 0)
            {
                tempL = mpo[loc1];
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1});
                tempL.SetLabel({-1, -3});
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() == 0)
            {

                tempL = mpo[loc1];

                VTs[i].SetLabel({1, -2, -1}); 
                tempL.SetLabel({-1, -3});
                Vs[i].SetLabel({-2, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() == 0 && mpo[loc2].BondNum() != 0)
            {
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1}); 
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -1, 2});

                temp = uni10::Contract(VTs[i], tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else
            {
                ostringstream err;
                err << "Error: SzSz Correlation fail contraction.";
                throw runtime_error(err.str());
            }

            mpo[loc1] = temp;
            mpo.erase(mpo.begin() + loc2);
        }
        else
        {
            mpo.erase(mpo.begin() + loc2);
        }
    }

    corr = temp[0];
    return corr;
}

/// Stot-Stot Correlation Stot = 0.5*(SmSp + SpSm) + SzSz
double Correlation_StSt(int site1, int site2, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double corr = 0;
    int L = Vs_loc.size() + 1;
    double spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sm(bonds);
    uni10::UniTensor<double> Sp(bonds);
    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> kara;
    Sm.PutBlock(matSm(spin));
    Sp.PutBlock(matSp(spin));
    Sz.PutBlock(matSz(spin));

    vector<uni10::UniTensor<double> > mpo;
    mpo.assign(L, kara);

    vector<vector<uni10::UniTensor<double> > > mpos;
    mpos.assign(3, mpo);

    mpos[0][site1] = 0.5*Sp;
    mpos[0][site2] = Sm;

    mpos[1][site1] = 0.5*Sm;
    mpos[1][site2] = Sp;

    mpos[2][site1]  = Sz;
    mpos[2][site2]  = Sz;

    vector<uni10::UniTensor<double> > temp;
    temp.assign(3, kara);
    uni10::UniTensor<double> tempL, tempR;
    int loc1, loc2;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpos[2].size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;
    
        if (mpos[2][loc1].BondNum() != 0 || mpos[2][loc2].BondNum() != 0)
        {
            if (mpos[2][loc1].BondNum() != 0 && mpos[2][loc2].BondNum() != 0)
            {
                for (int S=0; S<3; S++)
                {
                    tempL = mpos[S][loc1];
                    tempR = mpos[S][loc2];

                    VTs[i].SetLabel({1, -2, -1});
                    tempL.SetLabel({-1, -3});
                    tempR.SetLabel({-2, -4});
                    Vs[i].SetLabel({-4, -3, 2});
                    //cout << loc2 << " <- loc2 \n" << VTs[i] << tempL << tempR << Vs[i] << endl;

                    temp[S] = uni10::Contract(VTs[i], tempL);
                    temp[S] = uni10::Contract(temp[S], tempR);
                    temp[S] = uni10::Contract(temp[S], Vs[i]);
                }
            }
            else if (mpos[2][loc1].BondNum() != 0 && mpos[2][loc2].BondNum() == 0)
            {
                for (int S=0; S<3; S++)
                {
                    tempL = mpos[S][loc1];

                    VTs[i].SetLabel({1, -2, -1}); 
                    tempL.SetLabel({-1, -3});
                    Vs[i].SetLabel({-2, -3, 2});

                    temp[S] = uni10::Contract(VTs[i], tempL);
                    temp[S] = uni10::Contract(temp[S], Vs[i]);
                }
            }
            else if (mpos[2][loc1].BondNum() == 0 && mpos[2][loc2].BondNum() != 0)
            {
                for (int S=0; S<3; S++)
                {
                    tempR = mpos[S][loc2];

                    VTs[i].SetLabel({1, -2, -1}); 
                    tempR.SetLabel({-2, -4});
                    Vs[i].SetLabel({-4, -1, 2});

                    temp[S] = uni10::Contract(VTs[i], tempR);
                    temp[S] = uni10::Contract(temp[S], Vs[i]);
                }
            }
            else
            {
                ostringstream err;
                err << "Error: StSt Correlation fail contraction.";
                throw runtime_error(err.str());
            }

            for (int S=0; S<3; S++)
            {
                mpos[S][loc1] = temp[S];
                mpos[S].erase(mpos[S].begin() + loc2);
            }
        }
        else
        {
            for (int S=0; S<3; S++)
            {
                mpos[S].erase(mpos[S].begin() + loc2);
            }
        }
    }

    for (int S=0; S<3; S++)
    {
        corr += temp[S][0];
    }


    return corr;
}

/// <Sz>>
double Correlation_Sz(int site1, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double corr = 0;
    int L = Vs_loc.size() + 1;
    double spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> kara;
    Sz.PutBlock(matSz(spin));

    vector<uni10::UniTensor<double> > mpo;
    mpo.assign(L, kara);
    
    mpo[site1] = Sz;
    
    int loc1, loc2;
    uni10::UniTensor<double> temp, tempL, tempR;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpo.size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;
        
        if (mpo[loc1].BondNum() != 0 || mpo[loc2].BondNum() != 0)
        {
            if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() != 0)
            {
                tempL = mpo[loc1];
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1});
                tempL.SetLabel({-1, -3});
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() == 0)
            {

                tempL = mpo[loc1];

                VTs[i].SetLabel({1, -2, -1}); 
                tempL.SetLabel({-1, -3});
                Vs[i].SetLabel({-2, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() == 0 && mpo[loc2].BondNum() != 0)
            {
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1}); 
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -1, 2});

                temp = uni10::Contract(VTs[i], tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else
            {
                ostringstream err;
                err << "Error: SzSz Correlation fail contraction.";
                throw runtime_error(err.str());
            }

            mpo[loc1] = temp;
            mpo.erase(mpo.begin() + loc2);
        }
        else
        {
            mpo.erase(mpo.begin() + loc2);
        }
    }

    corr = temp[0];
    return corr;
}

/// <Stot> 
double Correlation_St(int site1, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double corr = 0;
    int L = Vs_loc.size() + 1;
    float spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sx(bonds);
    uni10::UniTensor<double> iSy(bonds);
    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> kara;
    Sx.PutBlock(matSx(spin));
    iSy.PutBlock(matiSy(spin));
    Sz.PutBlock(matSz(spin));

    for (int S=0; S<=2; S++)
    {
        vector<uni10::UniTensor<double> > mpo;
        uni10::UniTensor<double> kara;
        mpo.assign(L, kara);
    
        if (S == 0)
        {
            mpo[site1] = Sx;
        }
        else if (S == 1)
        {
            mpo[site1] = iSy;
        }
        else
        {
            mpo[site1] = Sz;
        }

        int loc1, loc2;
        uni10::UniTensor<double> temp, tempL, tempR;
        for (int i=0; i<Vs_loc.size(); i++)
        {
            loc1 = Vs_loc[i];
            if (loc1 == mpo.size() - 1) 
                loc2 = 0;
            else
                loc2 = loc1 + 1;
        
            if (mpo[loc1].BondNum() != 0 || mpo[loc2].BondNum() != 0)
            {
                if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() == 0)
                {
                    tempL = mpo[loc1];

                    VTs[i].SetLabel({1, -2, -1}); 
                    tempL.SetLabel({-1, -3});
                    Vs[i].SetLabel({-2, -3, 2});

                    temp = uni10::Contract(VTs[i], tempL);
                    temp = uni10::Contract(temp, Vs[i]);
                }
                else if (mpo[loc1].BondNum() == 0 && mpo[loc2].BondNum() != 0)
                {
                    tempR = mpo[loc2];

                    VTs[i].SetLabel({1, -2, -1}); 
                    tempR.SetLabel({-2, -4});
                    Vs[i].SetLabel({-4, -1, 2});


                    temp = uni10::Contract(VTs[i], tempR);
                    temp = uni10::Contract(temp, Vs[i]);
                }
                else
                {
                    ostringstream err;
                    err << "Error: StSt Correlation fail contraction.";
                    throw runtime_error(err.str());
                }
                mpo[loc1] = temp;
                mpo.erase(mpo.begin() + loc2);
            }
            else
            {
                mpo.erase(mpo.begin() + loc2);
            }
            
        }

        if (S == 1) // Sy = iSy/i  
        {
            if (temp[0] > pow(10, -8))
            {
                ostringstream err;
                err << "Error: StSt Correlation fail contraction.";
                throw runtime_error(err.str());
            }
        }
        else
        {
            corr += temp[0];
        }
    }

    return corr;
}

/// String order parameter
double Correlation_String(int site1, int site2, vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double corr = 0;
    int L = Vs_loc.size() + 1;
    double spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> kara;
    Sz.PutBlock(matSz(spin));

    vector<uni10::UniTensor<double> > mpo;
    mpo.assign(L, kara);
    
    mpo[site1] = Sz;
    mpo[site2] = Sz;

    uni10_complex128 t(0, M_PI);
    uni10::UniTensor<uni10_complex128> Sz_kC(bonds);
    Sz_kC.PutBlock(uni10::ExpH(t, matSz(spin) ));
    uni10::UniTensor<double> Sz_k = Sz_kC;

    /// for spin-1
    /*uni10_double64 mat_elem[] = {\
       -1.0, 0.0, 0.0,\
        0.0, 1.0, 0.0,\
        0.0, 0.0,-1.0};
    uni10::UniTensor<double> Sz_k(bonds);
    Sz_k.PutBlock(uni10::Matrix<uni10_double64>(3, 3, mat_elem) );
    //cout << Sz_k << endl;*/

    for (int i=site1+1; i<site2; i++)
        mpo[i] = Sz_k;

    int loc1, loc2;
    uni10::UniTensor<double> temp, tempL, tempR;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpo.size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;
        
        if (mpo[loc1].BondNum() != 0 || mpo[loc2].BondNum() != 0)
        {
            if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() != 0)
            {
                tempL = mpo[loc1];
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1});
                tempL.SetLabel({-1, -3});
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() != 0 && mpo[loc2].BondNum() == 0)
            {

                tempL = mpo[loc1];

                VTs[i].SetLabel({1, -2, -1}); 
                tempL.SetLabel({-1, -3});
                Vs[i].SetLabel({-2, -3, 2});

                temp = uni10::Contract(VTs[i], tempL);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else if (mpo[loc1].BondNum() == 0 && mpo[loc2].BondNum() != 0)
            {
                tempR = mpo[loc2];

                VTs[i].SetLabel({1, -2, -1}); 
                tempR.SetLabel({-2, -4});
                Vs[i].SetLabel({-4, -1, 2});

                temp = uni10::Contract(VTs[i], tempR);
                temp = uni10::Contract(temp, Vs[i]);
            }
            else
            {
                ostringstream err;
                err << "Error: SzSz Correlation fail contraction.";
                throw runtime_error(err.str());
            }

            mpo[loc1] = temp;
            mpo.erase(mpo.begin() + loc2);
        }
        else
        {
            mpo.erase(mpo.begin() + loc2);
        }
    }

    corr = -temp[0];
    return corr;
}

/// Magnetization of z-axis
double Magnetization_Sz(vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double M = 0;
    int L = Vs_loc.size() + 1;
    double spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sz(bonds);
    Sz.PutBlock(matSz(spin));

    vector<uni10::UniTensor<double> > mpo;
    mpo.assign(L, Sz);
    
    int loc1, loc2; 
    uni10::UniTensor<double> temp, tempL, tempR;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpo.size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;

        tempL = mpo[loc1];
        tempR = mpo[loc2];

        VTs[i].SetLabel({1, -2, -1});
        tempL.SetLabel({-1, -3});
        tempR.SetLabel({-2, -4});
        Vs[i].SetLabel({-4, -3, 2});
                
        temp = uni10::Contract(VTs[i], tempL);
        temp = uni10::Contract(temp, tempR);
        temp = uni10::Contract(temp, Vs[i]);
        mpo[loc1] = temp;
        mpo.erase(mpo.begin() + loc2);
    }

    M = temp[0];
    
    return M;
}

/// Magnetization of x-axis
double Magnetization_Sx(vector<uni10::UniTensor<double> > VTs, vector<uni10::UniTensor<double> > Vs, vector<int> Vs_loc)
{
    double M = 0;
    int L = Vs_loc.size() + 1;
    double spin = (sqrt(Vs[0].GetBlock().row() ) - 1)/2;                   // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)

    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sx(bonds);
    Sx.PutBlock(matSx(spin));

    vector<uni10::UniTensor<double> > mpo;
    mpo.assign(L, Sx);
    
    int loc1, loc2;
    uni10::UniTensor<double> temp, tempL, tempR;
    for (int i=0; i<Vs_loc.size(); i++)
    {
        loc1 = Vs_loc[i];
        if (loc1 == mpo.size() - 1) 
            loc2 = 0;
        else
            loc2 = loc1 + 1;

        tempL = mpo[loc1];
        tempR = mpo[loc2];

        VTs[i].SetLabel({1, -2, -1});
        tempL.SetLabel({-1, -3});
        tempR.SetLabel({-2, -4});
        Vs[i].SetLabel({-4, -3, 2});
                
        temp = uni10::Contract(VTs[i], tempL);
        temp = uni10::Contract(temp, tempR);
        temp = uni10::Contract(temp, Vs[i]);
        mpo[loc1] = temp;
        mpo.erase(mpo.begin() + loc2);
    }

    M = temp[0];
    
    return M;
}