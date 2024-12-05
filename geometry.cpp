#include "geometry.hpp"
#include <set>

namespace symmetry {
	std::vector<int> perm_from_cycles(
			const std::vector<std::vector<int>>& cycles){
		std::vector<int> res;
		int max_idx = 0;
		
		std::set<int> used_indices;
		for (auto& c : cycles){
			for (auto& i : c){
				if (i > max_idx) max_idx = i;
				if (i < 0) throw std::logic_error("Only positive indices allowed");
				auto res = used_indices.insert(i);
				if(res.second == false) {
					throw std::logic_error("An index was used twice in a cycle-spec");
				}
			}
		}
		std::vector<int> retval;
		retval.resize(max_idx+1);
		// fill with the identity permutation
		for (int J=0; J<max_idx; J++){
			retval[J] = J;
		}

		for (auto& c : cycles){
			for (int J=0; J<c.size(); J++){
				retval[c[J]] = c[(J+1)%c.size()];
			}
		}

		return retval;
	}



	std::vector<int> C3(){
		return perm_from_cycles(C3_cycles);
	}
	std::vector<int> sigma(){
		return perm_from_cycles(sigma_cycles);
	}
	std::vector<int> I(){
		return perm_from_cycles(I_cycles);
	}


};
