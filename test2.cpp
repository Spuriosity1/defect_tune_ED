#include "geometry.hpp"
#include <xdiag/all.hpp>
#include <xdiag/operators/opsum.hpp>

using namespace xdiag;


const int N_SITES = 7;

/* Set up spins in opposing pyramids:
 *      1
 *     2 3 
 *
 *      0
 *
 *     6 5
 *      4
 *
 *  Top View:
 *
 *       1
 *  6         5
 *       0
 *  2         3
 *       4
 *
 */
const std::vector<std::pair<int, int>> cluster_bonds = {
  {1,2},
  {2,3},
  {3,1},
  {0,1},
  {0,2},
  {0,3},
  {0,4},
  {0,5},
  {0,6},
  {4,5},
  {5,6},
  {6,4}
};


int main(int argc, char* argv[]) try {
 
  OpSum bonds;

  // Generate all the bonds
  for (auto& [i1, i2] : cluster_bonds) {
    bonds += Op("HB", "J", {i1, i2});
  }

  bonds["J"] = 1.0; // Sets bond strength of bondspec "J"

  // generate the symmetries
    //             {0, 1, 2, 3, 4, 5, 6}
  Permutation Id = {0, 1, 2, 3, 4, 5, 6};
  Permutation C3 = {0, 2, 3, 1, 5, 6, 4};
  Permutation C2  = {0, 1, 3, 2, 4, 6, 5};
  Permutation I = {0, 4, 5, 6, 1, 2, 3};

  set_verbosity(3);

  // Create the symmetry group
  auto group = PermutationGroup({Id,    C3,    C3*C3,    I,    I*C3, I*C3*C3,
                                 C2, C2*C3, C2*C3*C3, C2*I, C2*I*C3, C2*I*C3*C3});
   
  XDIAG_SHOW(group.size());


  Representation irrep(symmetry::irreps_D3d[argv[1]]);


  XDIAG_SHOW(irrep);

  auto block = Spinhalf(N_SITES, group, irrep);
  double e0 = eigval0(bonds, block);
  Log("e0: {}", e0);

} catch (Error e) {
  error_trace(e);
}
