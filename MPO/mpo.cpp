#include "mpo.h"

MPO::MPO(string model_name, char loc, float spin, ...)
{
    model = model_name;
    vector<double> para;
    va_list vl;
    va_start(vl, spin);
    if (model_name == "XXZ_OBC")
    {
        for (int i=0; i<3; i++)
            para.push_back(va_arg(vl, double) );

        MPO_XXZ_OBC(loc, spin, para[0], para[1], para[2]);     // para : Jx, Jz, h
    }
    else if (model_name == "XXZ_PBC")
    {
        for (int i=0; i<3; i++)
            para.push_back(va_arg(vl, double) );
            
        MPO_XXZ_PBC(loc, spin, para[0], para[1], para[2]);    // para : Jx, Jz, h
    }
    else if (model_name == "Ising")
    {
        for (int i=0; i<3; i++)
            para.push_back(va_arg(vl, double) );

        MPO_Ising(loc, spin, para[0], para[1], para[2]);      // para : J, hx, hz
    }
    else
    {
        ostringstream err;
        err << "Error: No support this model.";
        throw runtime_error(err.str());
    }
    va_end(vl);
}

/// Declare MPO of Ising model (OBC)
void MPO::MPO_Ising(char loc, float spin, float J, float hx, float hz) 
{
    mpo_loc  = loc;
    phys_dim = (int)(2*spin) + 1 ;
    virt_dim = 3;
    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> Sx(bonds);
    uni10::UniTensor<double> Id(bonds);
    uni10::UniTensor<double> Zero(bonds);

    Sz.PutBlock(2*matSz(spin)); // matSz = {1/2 0; 0 -1/2}
    Sx.PutBlock(2*matSx(spin));
    Id.Identity();
    Zero.Zeros();
    
    switch (loc)
    {
        case 'l':
            mpo_l.assign(virt_dim, Zero);
            mpo_l[0] = -hx*Sx - hz*Sz;
            mpo_l[1] = J*Sz;
            mpo_l[2] = Id;

            break;

        case 'm':
            mpo_l.assign(virt_dim, Zero);
            mpo_l[0] = -hx*Sx - hz*Sz;
            mpo_l[1] = J*Sz;
            mpo_l[2] = Id;

            mpo_r.assign(virt_dim, Zero);
            mpo_r[0] = Id;
            mpo_r[1] = Sz;
            mpo_r[2] = -hx*Sx - hz*Sz;

            mpo_m.assign(virt_dim*virt_dim, Zero);
            mpo_m[0] = mpo_r[0];
            mpo_m[3] = mpo_r[1];
            mpo_m[6] = mpo_r[2];
            mpo_m[7] = mpo_l[1];
            mpo_m[8] = mpo_l[2];

            break;

        case 'r':
            mpo_r.assign(virt_dim, Zero);
            mpo_r[0] = Id;
            mpo_r[1] = Sz;
            mpo_r[2] = -hx*Sx - hz*Sz;

            break;
    }
    
}

/// Declare MPO of XXZ model (OBC)
void MPO::MPO_XXZ_OBC(char loc, float spin, float Jx, float Jz, float h) 
{
    mpo_loc  = loc;
    phys_dim = (int)(2*spin) + 1 ;
    virt_dim = 5;
    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sp(bonds);
    uni10::UniTensor<double> Sm(bonds);
    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> Id(bonds);
    uni10::UniTensor<double> Zero(bonds);

    Sp.PutBlock(matSp(spin));
    Sm.PutBlock(matSm(spin));
    Sz.PutBlock(matSz(spin));
    Id.Identity();
    Zero.Zeros();
    
    switch (loc)
    {
        case 'l':
            mpo_l.assign(virt_dim, Zero);
            mpo_l[0] = -h*Sz;
            mpo_l[1] = 0.5*Jx*Sm;
            mpo_l[2] = 0.5*Jx*Sp;
            mpo_l[3] = Jz*Sz;
            mpo_l[4] = Id;

            break;

        case 'm':
            mpo_l.assign(virt_dim, Zero);
            mpo_l[0] = -h*Sz;
            mpo_l[1] = 0.5*Jx*Sm;
            mpo_l[2] = 0.5*Jx*Sp;
            mpo_l[3] = Jz*Sz;
            mpo_l[4] = Id;

            mpo_r.assign(virt_dim, Zero);
            mpo_r[0] = Id;
            mpo_r[1] = Sp;
            mpo_r[2] = Sm;
            mpo_r[3] = Sz;
            mpo_r[4] = -h*Sz;

            mpo_m.assign(virt_dim*virt_dim, Zero);
            mpo_m[0]  = mpo_r[0];
            mpo_m[5]  = mpo_r[1];
            mpo_m[10] = mpo_r[2];
            mpo_m[15] = mpo_r[3];
            mpo_m[20] = mpo_r[4];
            mpo_m[21] = mpo_l[1];
            mpo_m[22] = mpo_l[2];
            mpo_m[23] = mpo_l[3];
            mpo_m[24] = mpo_l[4];

            break;

        case 'r':
            mpo_r.assign(virt_dim, Zero);
            mpo_r[0] = Id;
            mpo_r[1] = Sp;
            mpo_r[2] = Sm;
            mpo_r[3] = Sz;
            mpo_r[4] = -h*Sz;

            break;
    }
    
}

void MPO::MPO_XXZ_PBC(char loc, float spin, float Jx, float Jz, float h) 
{
    mpo_loc  = loc;
    phys_dim = (int)(2*spin) + 1 ;
    virt_dim = 5;
    uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
    uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
    vector<uni10::Bond> bonds = {bdi, bdo};

    uni10::UniTensor<double> Sp(bonds);
    uni10::UniTensor<double> Sm(bonds);
    uni10::UniTensor<double> Sz(bonds);
    uni10::UniTensor<double> Id(bonds);
    uni10::UniTensor<double> Zero(bonds);

    Sp.PutBlock(matSp(spin));
    Sm.PutBlock(matSm(spin));
    Sz.PutBlock(matSz(spin));
    Id.Identity();
    Zero.Zeros();
    
    switch (loc)
    {
        case 'l':
            mpo_l.assign(virt_dim, Zero);
            mpo_l[0] = -h*Sz;
            mpo_l[1] = 0.5*Jx*Sm;
            mpo_l[2] = 0.5*Jx*Sp;
            mpo_l[3] = Jz*Sz;
            mpo_l[4] = Id;

            break;

        case 'm':
            mpo_l.assign(virt_dim, Zero);
            mpo_l[0] = -h*Sz;
            mpo_l[1] = 0.5*Jx*Sm;
            mpo_l[2] = 0.5*Jx*Sp;
            mpo_l[3] = Jz*Sz;
            mpo_l[4] = Id;

            mpo_r.assign(virt_dim, Zero);
            mpo_r[0] = Id;
            mpo_r[1] = Sp;
            mpo_r[2] = Sm;
            mpo_r[3] = Sz;
            mpo_r[4] = -h*Sz;

            mpo_m.assign(virt_dim*virt_dim, Zero);
            mpo_m[0]  = mpo_r[0];
            mpo_m[5]  = mpo_r[1];
            mpo_m[10] = mpo_r[2];
            mpo_m[15] = mpo_r[3];
            mpo_m[20] = mpo_r[4];
            mpo_m[21] = mpo_l[1];
            mpo_m[22] = mpo_l[2];
            mpo_m[23] = mpo_l[3];
            mpo_m[24] = mpo_l[4];

            break;

        case 'r':
            mpo_r.assign(virt_dim, Zero);
            mpo_r[0] = Id;
            mpo_r[1] = Sp;
            mpo_r[2] = Sm;
            mpo_r[3] = Sz;
            mpo_r[4] = -h*Sz;

            break;
    }
    
}

/// create Tensor of MPO 
void MPO::Launch(uni10::UniTensor<double> location, vector<uni10::UniTensor<double> > mpo_frame, vector<int> v, uni10::UniTensor<double>& mpo)
{
    for (auto it = v.begin(); it != v.end(); it++)
    {
        location.Zeros();
        uni10::Matrix<double> temp = location.GetBlock();
        temp[*it] = 1.0;
        location.PutBlock(temp);
        if (mpo.BondNum() == 0)
            mpo = uni10::Otimes(location, mpo_frame[*it]);
        else
            mpo +=  uni10::Otimes(location, mpo_frame[*it]);
    }
}

/// return UniTensor of MPO (loc)
uni10::UniTensor<double> MPO::GetTensor()
{
    char loc = mpo_loc;
    return GetTensor(loc);
}

/// return MPO
uni10::UniTensor<double> MPO::GetTensor(char loc)      // loc = mpo_loc ==> invalid use of non-static data member
{
    uni10::UniTensor<double> mpo;
    vector<int> v;                                     // given non-zero vlaue location

    /// check edge
    if ( (mpo_loc == 'l' && loc != 'l' )  || (mpo_loc == 'r' && loc != 'r')  )
    {
        ostringstream err;
        err << "Error: NOTICE edge!!!";
        throw runtime_error(err.str());
    }

    switch (loc)
    {
        case 'l':
            {                                          // important """{ }""" avoid to redeclarate
                /// Declare bonds and location tensor  (5) x (2, 2) = (2, 10)
                vector<uni10::Bond> bonds;
                bonds.push_back(uni10::Bond(uni10::BD_OUT, virt_dim) );
                uni10::UniTensor<double> location(bonds);

                for (int i=0; i<virt_dim; i++)
                    v.push_back(i);
        
                Launch(location, mpo_l, v, mpo);
            }
            break;

        case 'm':
            {
                /// Declare bonds and location tensor  (5, 5) x (2, 2) = (10, 10)
                vector<uni10::Bond> bonds;
                bonds.push_back(uni10::Bond(uni10::BD_IN, virt_dim) );
                bonds.push_back(uni10::Bond(uni10::BD_OUT, virt_dim) );
                uni10::UniTensor<double> location(bonds);

                for (int i=0; i<virt_dim; i++)
                    v.push_back(virt_dim*i);

                for (int i=1; i<virt_dim; i++)
                    v.push_back(virt_dim*(virt_dim-1) + i);

                Launch(location, mpo_m, v, mpo);
            }
            break;

        case 'r':
            {
                /// Declare bonds and location tensor  (5) x (2, 2) = (10, 2)
                vector<uni10::Bond> bonds;
                bonds.push_back(uni10::Bond(uni10::BD_IN, virt_dim) );
                uni10::UniTensor<double> location(bonds);

                for (int i=0; i<virt_dim; i++)
                    v.push_back(i);
        
                Launch(location, mpo_r, v, mpo);
            }
            break;
    }

    return mpo;
}

/// return MPO of PBC need last coupling
uni10::UniTensor<double> MPO::GetTensorPBC(char loc)      // Q? loc = mpo_loc ==> invalid use of non-static data member
{
    uni10::UniTensor<double> mpo;
    vector<int> v;                                     // given non-zero vlaue location
    /// check edge
    if ( (mpo_loc == 'l' && loc != 'l' )  || (mpo_loc == 'r' && loc != 'r')  )
    {
        ostringstream err;
        err << "Error: NOTICE edge!!!";
        throw runtime_error(err.str());
    }

    switch (loc)
    {
        case 'l':
            {                                          // important """{ }""" avoid to redeclarate
                /// Declare bonds and location tensor  (5) x (2, 2) = (2, 10)
                vector<uni10::Bond> bonds;
                bonds.push_back(uni10::Bond(uni10::BD_OUT, virt_dim) );
                uni10::UniTensor<double> location(bonds);

                for (int i=1; i<virt_dim-1; i++)
                    v.push_back(i);
        
                Launch(location, mpo_l, v, mpo);
            }
            break;

        case 'm':
            {
                ostringstream err;
                err << "Error: PBC last part need loc='l' or 'r'";
                throw runtime_error(err.str());
            }
            break;

        case 'r':
            {
                /// Declare bonds and location tensor  (5) x (2, 2) = (10, 2)
                vector<uni10::Bond> bonds;
                bonds.push_back(uni10::Bond(uni10::BD_IN, virt_dim) );
                uni10::UniTensor<double> location(bonds);

                for (int i=1; i<virt_dim-1; i++)
                    v.push_back(i);
        
                Launch(location, mpo_r, v, mpo);
            }
            break;
    }

    return mpo;
}

/// return MPO of single site = SS
uni10::UniTensor<double> MPO::GetTensorSS()
{
    uni10::UniTensor<double> mpo;

    switch (mpo_loc)
    {
        case 'l':
            {    
                mpo = mpo_l[0];
            }
            break;

        case 'm':
            {
                int single = virt_dim * (virt_dim - 1);
                mpo = mpo_m[single];
            }
            break;

        case 'r':
            {
                mpo = mpo_r.back();
            }
            break;
    }

    return mpo;
}

/// renormalize group of two MPO
void RG_MPO(MPO &W1, MPO& W2, uni10::UniTensor<double> VT, uni10::UniTensor<double> V)
{
    int row = W1.virt_dim;
    VT.SetLabel({-1, 0, 1});
    V.SetLabel({2, 3, -2});

    if (W1.mpo_loc == 'l')
    {
        vector<uni10::UniTensor<double> > temp_l;

        for (int i=0; i<row; i++)
        {   
            uni10::UniTensor<double> temp;
            for (int j=0; j<row; j++)
            {
                if (temp.BondNum() == 0)
                    temp = uni10::Otimes(W1.mpo_l[j], W2.mpo_m[i + j*row]);
                else
                    temp += uni10::Otimes(W1.mpo_l[j], W2.mpo_m[i + j*row]);
            }
            temp = uni10::Contract(VT, temp);
            temp = uni10::Contract(temp, V);
            temp_l.push_back(temp);
        }

        W1.mpo_l = temp_l;
    }
    else if (W2.mpo_loc == 'r')
    {
        vector<uni10::UniTensor<double> > temp_r;

        for (int i=0; i<row; i++)
        {   
            uni10::UniTensor<double> temp;
            for (int j=0; j<row; j++)
            {
                if (temp.BondNum() == 0)
                    temp = uni10::Otimes(W1.mpo_m[i*row + j], W2.mpo_r[j]);
                else
                    temp += uni10::Otimes(W1.mpo_m[i*row + j], W2.mpo_r[j]);
            }
            temp = uni10::Contract(VT, temp);
            temp = uni10::Contract(temp, V);
            temp_r.push_back(temp);
        }

        W2.mpo_r = temp_r;
    }
    else // W1.mpo_loc == 'm' && W2.mpo_loc == 'm'
    {
        vector<uni10::UniTensor<double> > temp_l;
	    vector<uni10::UniTensor<double> > temp_m;
	    vector<uni10::UniTensor<double> > temp_r;

        uni10::Bond bdi(uni10::BD_IN, V.GetBlock().col() );           // row * col; Vs = (3*3 x chi); VTs = (chi x 3*3)
        uni10::Bond bdo(uni10::BD_OUT, V.GetBlock().col() );
        vector<uni10::Bond> bonds = {bdi, bdo}; 
        uni10::UniTensor<double> Zero(bonds);
        Zero.Zeros();
        temp_m.assign(row*row, Zero);

        for (int i=0; i<row; i++)
        {   
            uni10::UniTensor<double> temp;
            for (int j=0; j<row; j++)
            {
                if (temp.BondNum() == 0)
                    temp = uni10::Otimes(W1.mpo_m[i*row + j], W2.mpo_r[j]);
                else
                    temp += uni10::Otimes(W1.mpo_m[i*row + j], W2.mpo_r[j]);
            }

            temp = uni10::Contract(VT, temp);
            temp = uni10::Contract(temp, V);
            temp_r.push_back(temp);
            temp_m[i*row] = temp;
        }

        for (int i=0; i<row; i++)
        {   
            uni10::UniTensor<double> temp;
            for (int j=0; j<row; j++)
            {
                if (temp.BondNum() == 0)
                    temp = uni10::Otimes(W1.mpo_l[j], W2.mpo_m[i + j*row]);
                else
                    temp += uni10::Otimes(W1.mpo_l[j], W2.mpo_m[i + j*row]);
            }
            temp = uni10::Contract(VT, temp);
            temp = uni10::Contract(temp, V);
            temp_l.push_back(temp);
            temp_m[row*(row-1) + i] = temp;
        }

        W1.mpo_l = temp_l;
        W1.mpo_m = temp_m;
        W1.mpo_r = temp_r;
    } 
}
