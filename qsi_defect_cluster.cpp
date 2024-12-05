#include <complex>
#include <fstream>
#include <nlohmann/json_fwd.hpp>
#include <sstream>
#include <xdiag/all.hpp>
#include "geometry.hpp"
#include <nlohmann/json.hpp>

using namespace xdiag;
using namespace arma;
using json=nlohmann::json;

const int N_SITES = 24;

// converts from 1-based to zero based
// this is a trick to remove the nonexistent defect site from the calculation
int convert_index(int id){
  assert(id-1 >= 0);
  assert(id-1 < N_SITES);
  return id-1;
}

std::string pprint(const std::vector<int>& v){
  std::stringstream s;
  s << "[ ";
  for (auto& val : v){
    s << val << " ";
  }
  s << "]";
  return s.str();
}


const char* field_proj_label[4] = {"h0", "h1", "h2", "h3"}; 

const static std::vector<ivec3> pyro_pos = {
	{1,1,1},
	{1,-1,-1},
	{-1,1,-1},
	{-1,-1,1}
};


int main(int argc, char* argv[]) try { 

	if (argc < 7){
		cout <<"USAGE: "<<argv[0]<<" <out_folder> <Jpm/Jzz> <hx> <hy> <hz> <lanczos_dim> [num_kept_states = 4]\n";
		return 1;
	}

  json out;

	set_verbosity(1);// set verbosity for monitoring progress
	
  std::string out_folder(argv[1]);
    double Jpm = atof(argv[2]);
	arma::vec3 bfield = {atof(argv[3]), atof(argv[4]), atof(argv[5])};
	int lanczos_dim = atoi(argv[6]); 
	int num_kept_states = 4;
	if (argc >= 8) num_kept_states = atoi(argv[7]);

    out["Jpm"] = Jpm;
    out["B"] = bfield;
    out["lanczos_dim"] = lanczos_dim;


  OpSum ops;


	for (int mu=0; mu<4; mu++){
		ops[field_proj_label[mu]] = (double) arma::dot(bfield,pyro_pos[mu])/sqrt(3);
		cout << "bfield.e_mu "<<mu<<" "<<ops[field_proj_label[mu]]<<"\n";
	}

  // Generate the XXZ ops
  for (auto& [i1, i2] : pyro_bonds_internal) {
    if (i1 == 0 || i2 == 0) continue; // delete ops connected to the defect
    ops += Op("ISING", "Jy", {convert_index(i1), convert_index(i2)});
    ops += Op("EXCHANGE", "Jpm_neg", {convert_index(i1), convert_index(i2)});
  }

  for (int J=0; J<N_SITES; J++){
    switch (sublattice[J+1]) {
      case 0: ops += Op("Sx", "h0", J); break;
      case 1: ops += Op("Sx", "h1", J); break;
      case 2: ops += Op("Sx", "h2", J); break;
      case 3: ops += Op("Sx", "h3", J); break;
    }
  }

  ops["Jy"] = 1.; // Sets bond strength of the Ising coupling
  ops["Jpm_neg"] = Jpm; // Sets bond strength of the XXZ hopping (note sign)
                           
  ops["h0"] = arma::dot(bfield, arma::vec3({ 1, 1, 1})) / sqrt(3);
  ops["h1"] = arma::dot(bfield, arma::vec3({ 1,-1,-1})) / sqrt(3);
  ops["h2"] = arma::dot(bfield, arma::vec3({-1, 1,-1})) / sqrt(3);
  ops["h3"] = arma::dot(bfield, arma::vec3({-1,-1, 1})) / sqrt(3);

  auto C3 = Permutation(symmetry::C3());
  auto sigma = Permutation(symmetry::sigma());
  auto I = Permutation(symmetry::I());
  auto ident = I*I;

  XDIAG_SHOW(C3);
  XDIAG_SHOW(sigma);
  XDIAG_SHOW(I);
  
  auto group = PermutationGroup({
      ident, C3, C3*C3,
      });
  std::vector<Representation> irreps = {
    generated_irrep(C3, 1),
    generated_irrep(C3, std::polar(2*pi/3)),
    generated_irrep(C3, std::polar(-2*pi/3))
  };


  set_verbosity(2);             // set verbosity for monitoring progress

  std::cout << ("Diagonalizing H\n");
  std::cout << "ready to go" <<std::endl;

  std::vector<arma::vec> E;
  for (int i=0; i<3; i++){
      std::cout << ("Diagonalising H, block ")<<i<<"\n";

    Spinhalf block(24, 0, group, irreps[i]);

	std::vector<std::string> statev = {"Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn","Up","Dn"};
	auto init_state = product(block, statev);

	eigs_lanczos_result_t lanczos_res = eigs_lanczos(ops, block, 
			//init_state, 
			lanczos_dim,
			/* precision */ 1e-14,
			/* max iter */ 10000,
			/* force complex */ false
			);
	    cerr << "iter criterion -> "<< lanczos_res.criterion << endl;
      cout << "Eigenvalues: " << lanczos_res.eigenvalues;
	    E.push_back( lanczos_res.eigenvalues);
  }

  out["energies"] = E;

  std::stringstream label;
  label << "%jpm=" << ops["Jpm"] << "%B=" << bfield[0]<<","<<bfield[1]<<","<<bfield[2];

	std::ofstream file(out_folder + "/out_pyro16_"+label.str()+".json");
  file << out;

} catch (Error e) {
  error_trace(e);
}
