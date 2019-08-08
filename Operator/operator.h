#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include "uni10.hpp"

void spin_check(uni10_float32 spin);

uni10::Matrix<double> matSx(uni10_float32 spin = 0.5);

uni10::Matrix<double> matiSy(uni10_float32 spin = 0.5);

uni10::Matrix<double> matSp(uni10_float32 spin = 0.5);

uni10::Matrix<double> matSm(uni10_float32 spin = 0.5);

uni10::Matrix<double> matSz(uni10_float32 spin = 0.5);

uni10::UniTensor<double> periodicSx(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<double> periodicSy(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<double> periodicSz(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<double> periodicsqSx(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<double> periodicsqSy(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<double> periodicsqSz(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

#endif
