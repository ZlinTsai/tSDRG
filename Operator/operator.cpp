#include "operator.h"
#include <cmath> 
#include <assert.h> 
#include <exception>

using namespace std;

const uni10_int32 CURRENT_SUPPORTED_SPIN_DIM = 5;

void spin_check(uni10_float32 spin){

  if(!(spin > 0 && floor(2 * spin) == 2 * spin)){
      std::ostringstream err;
      err<<"The given spin is not half integer.";
      throw std::runtime_error(err.str());
  }

}

uni10::Matrix<uni10_double64> matSp(uni10_float32 spin){

  spin_check(spin);

  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }

  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0, 1,\
      0, 0};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
      0.0, 0.0, 1.0,\
      0.0, 0.0, 0.0};
    return sqrt(2) * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
      0.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,     0.0,     0.0, sqrt(6), 0.0,\
      0.0,     0.0,     0.0,     0.0, 2.0,\
      0.0,     0.0,     0.0,     0.0, 0.0,};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }

  return uni10::Matrix<uni10_double64>();

}

uni10::Matrix<uni10_double64> matSm(uni10_float32 spin){

  spin_check(spin);
  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0, 0,\
      1, 0};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 0.0, 0.0,\
      1.0, 0.0, 0.0,\
      0.0, 1.0, 0.0};
    return sqrt(2) * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     0.0,     0.0,     0.0, 0.0,\
      2.0,     0.0,     0.0,     0.0, 0.0,\
      0.0, sqrt(6),     0.0,     0.0, 0.0,\
      0.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,     0.0,     0.0,     2.0, 0.0 };
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }

  return uni10::Matrix<uni10_double64>();

}

uni10::Matrix<uni10_double64> matSx(uni10_float32 spin){

  spin_check(spin);
  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0,   0.5,\
      0.5, 0  };
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
      1.0, 0.0, 1.0,\
      0.0, 1.0, 0.0};
    return (1.0 / sqrt(2)) * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
      2.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0, sqrt(6),     0.0, sqrt(6), 0.0,\
      0.0,     0.0, sqrt(6),     0.0, 2.0,\
      0.0,     0.0,     0.0,     2.0, 0.0 };
    return 0.5 * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
	
  return uni10::Matrix<uni10_double64>();

}

uni10::Matrix<uni10_double64> matiSy(uni10_float32 spin){

  spin_check(spin);
  int dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0,   0.5,\
     -0.5, 0};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
     -1.0, 0.0, 1.0,\
      0.0,-1.0, 0.0};
    return (1.0 / sqrt(2)) * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
     -2.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,-sqrt(6),     0.0, sqrt(6), 0.0,\
      0.0,     0.0,-sqrt(6),     0.0, 2.0,\
      0.0,     0.0,     0.0,    -2.0, 0.0,};
    return 0.5 * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }

  return uni10::Matrix<uni10_double64>();

}

uni10::Matrix<uni10_double64> matSz(uni10_float32 spin){

  spin_check(spin);
  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0.5,  0,\
      0,   -0.5  };
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      1.0, 0.0,  0.0,\
      0.0, 0.0,  0.0,\
      0.0, 0.0, -1.0};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      2.0, 0.0, 0.0, 0.0, 0.0,\
      0.0, 1.0, 0.0, 0.0, 0.0,\
      0.0, 0.0, 0.0, 0.0, 0.0,\
      0.0, 0.0, 0.0,-1.0, 0.0,\
      0.0, 0.0, 0.0, 0.0,-2.0,};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }

  return uni10::Matrix<uni10_double64>();

}

uni10::UniTensor<uni10_double64> periodicSx(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> Sx = matSx(spin);  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, Sx.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, Sx.col()));
  uni10::UniTensor<uni10_double64> perSx(bondI);
  perSx.PutBlock(Sx);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.Identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    perSx = Otimes(perSx, Id);

  uni10::UniTensor<uni10_double64> tmp = perSx;
  vector<uni10_int32> labels = perSx.label();
  vector<uni10_int32> per_labels(labels.size());
  uni10_int32 bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(uni10_int32 l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    perSx += Permute(tmp, per_labels, bondNum / 2);
    tmp.SetLabel(labels);
  }

  return perSx;

}

/*
uni10::UniTensor periodicSy(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix Sy = matSy(spin);  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, Sy.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, Sy.col()));
  uni10::UniTensor perSy(bondI);
  perSy.PutBlock(Sy);
  uni10::UniTensor Id(uni10::CTYPE, bondI);
  Id.Identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    perSy = Otimes(perSy, Id);

  uni10::UniTensor tmp = perSy;
  vector<int> labels = perSy.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    perSy += tmp.Permute(uni10::CTYPE, per_labels, bondNum / 2);
    tmp.SetLabel(labels);
  }

  return perSy;

}
*/

uni10::UniTensor<uni10_double64> periodicSz(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> Sz = matSz(spin);  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, Sz.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, Sz.col()));
  uni10::UniTensor<uni10_double64> perSz(bondI);
  perSz.PutBlock(Sz);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.Identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    perSz = Otimes(perSz, Id);
	
  uni10::UniTensor<uni10_double64> tmp = perSz;
  vector<int> labels = perSz.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    perSz += Permute(tmp, per_labels, bondNum / 2);
    tmp.SetLabel(labels);
  }

  return perSz;

}

uni10::UniTensor<uni10_double64> periodicsqSz(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> sqSz = uni10::Dot(matSz(spin), matSz(spin));  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, sqSz.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, sqSz.col()));
  uni10::UniTensor<uni10_double64> persqSz(bondI);
  persqSz.PutBlock(sqSz);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.Identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    persqSz = Otimes(persqSz, Id);

  uni10::UniTensor<uni10_double64> tmp = persqSz;
  vector<int> labels = persqSz.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();

  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    persqSz += Permute(tmp, per_labels, bondNum / 2);
    tmp.SetLabel(labels);
  }

  return persqSz;

}

uni10::UniTensor<uni10_double64> periodicsqSx(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> sqSx = uni10::Dot(matSx(spin), matSx(spin));  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, sqSx.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, sqSx.col()));
  uni10::UniTensor<uni10_double64> persqSx(bondI);
  persqSx.PutBlock(sqSx);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.Identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    persqSx = Otimes(persqSx, Id);

  uni10::UniTensor<uni10_double64> tmp = persqSx;
  vector<uni10_int32> labels = persqSx.label();
  vector<uni10_int32> per_labels(labels.size());
  uni10_int32 bondNum = per_labels.size();

  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    persqSx += Permute(tmp, per_labels, bondNum / 2);
    tmp.SetLabel(labels);
  }

  return persqSx;

}

uni10::UniTensor<uni10_double64> periodicsqSy(uni10_uint64 siteNum, uni10_float32 spin){

  spin_check(spin);
  int dim = spin * 2 + 1;
  uni10::Matrix<uni10_double64> matSy_Rpart;
  uni10::Matrix<uni10_double64> sqSy;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0, 1,\
     -1, 0};
    matSy_Rpart = uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
    sqSy = (-1.0/4.0) * uni10::Dot(matSy_Rpart, matSy_Rpart);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
     -1.0, 0.0, 1.0,\
      0.0,-1.0, 0.0};
    matSy_Rpart = uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
    sqSy = (-1.0/2.0) * uni10::Dot(matSy_Rpart, matSy_Rpart);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
     -2.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,-sqrt(6),     0.0, sqrt(6), 0.0,\
      0.0,     0.0,-sqrt(6),     0.0, 2.0,\
      0.0,     0.0,     0.0,    -2.0, 0.0,};
    matSy_Rpart = uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
    sqSy = (-1.0/4.0) * uni10::Dot(matSy_Rpart, matSy_Rpart);
  }

  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, sqSy.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, sqSy.col()));
  uni10::UniTensor<uni10_double64> persqSy(bondI);
  persqSy.PutBlock(sqSy);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.Identity();
  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    persqSy = Otimes(persqSy, Id);
  uni10::UniTensor<uni10_double64> tmp = persqSy;
  vector<int> labels = persqSy.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    persqSy += Permute(tmp, per_labels, bondNum / 2);
    tmp.SetLabel(labels);
  }
  return persqSy;
}


