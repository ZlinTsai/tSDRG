#ifndef _MPO_H_
#define _MPO_H_

#include <string>
#include <cstdarg>
#include <exception>
#include "../Operator/operator.h"
#include "../Hamiltonian/hamiltonian.h"

using namespace std;

/// Declare Matrix product operator
class MPO
{
	private:
		char mpo_loc;
		int phys_dim;
		int virt_dim;
		vector<uni10::UniTensor<double> > mpo_l;
		vector<uni10::UniTensor<double> > mpo_m;
		vector<uni10::UniTensor<double> > mpo_r;
		void Launch(uni10::UniTensor<double>, vector<uni10::UniTensor<double> >, vector<int>, uni10::UniTensor<double>&);

		/// MPO-model 
		void MPO_XXZ_OBC(char, float, float, float, float);     // XXZ-model,   loc, spin, Jx, Jz, h (OBC)
		void MPO_XXZ_PBC(char, float, float, float, float);     // XXZ-model,   loc, spin, Jx, Jz, h (PBC)
		void MPO_Ising(char, float, float, float, float);		// Ising-model, loc, spin, J, hx, hz (OBC)
	
	public:
		MPO(string, char, float, ...);

		string model;
		uni10::UniTensor<double> GetTensor();
		uni10::UniTensor<double> GetTensor(char);
		uni10::UniTensor<double> GetTensorPBC(char);
		uni10::UniTensor<double> GetTensorSS();
	
	friend void RG_MPO(MPO&, MPO&, uni10::UniTensor<double>, uni10::UniTensor<double>);
};

#endif
