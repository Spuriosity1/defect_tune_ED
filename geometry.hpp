#pragma once
#include <complex>
#include <stdexcept>
#include <vector>
#include <map>
#include <xdiag/common.hpp>

/*   Bottom Layer
 *              
 *         ()--------[]
 *        /            \  11
 *       /              \     10
 *     []                ()--------[]
 *       \              /            \ 
 *    12  \      8     / 7            \
 *         ()--------[]                ()
 *   13   /            \              /
 *       /           9  \            /  
 *     []                ()--------[]
 *       \              /    15
 *        \            / 14
 *         ()--------[]
 *                
 * 
 * [] denotes link to top layer
 *
 *  Links
 *
 *                   [] 2
 *                      
 *                       
 *     [] 3                        [] 1
 *                                    
 *                                     
 *                   [] 0                 
 *                                     
 *                                    
 *     [] 4                        [] 6
 *                       
 *                      
 *                   [] 5
 *
 *
 *  Top Layer
 *
 *
 *                    {}--------()
 *                20 /            \
 *          21      /              \
 *      {}--------()                {}
 *     /            \              /  
 *    /          17  \     16     /  19 
 *  ()                {}--------()    
 *    \              /            \ 24
 *     \            / 18           \  
 *      {}--------()                {}
 *          22      \              /
 *               23  \            /
 *                    {}--------()
 *
 *
 *
 * Sublattices
 *
 *
 *       [] 0  -> 0
 *      \  -> 1
 *      /  -> 2
 *      -- -> 3
 */

const std::vector<int> sublattice = {
//  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 
    0,0,0,0,0,0,0,2,3,1, 3, 1, 1, 2, 2, 3, 3, 1, 2, 2, 2, 3, 3, 1, 1
};

const std::vector<std::pair<unsigned, unsigned>> pyro_bonds_internal = {
	// bottom layer grouped by tetrahedron
	{0,7},{0,8},{0,9},{7,8},{8,9},{9,7},
	{7,10},{10,11},{11,7},
	{8,12},{12,13},{13,8},
	{9,14},{14,15},{15,9},
	// middle layer
	{10,1},{1,19},
	{11,2},{2,20},
	{12,3},{3,21},
	{13,4},{4,22},
	{14,5},{5,23},
	{15,6},{6,24},
	// top layer grouped by tetrahedron
	{0,16},{0,17},{0,18},{16,17},{17,18},{18,16},
	{16,19},{19,24},{24,16},
	{20,21},{21,17},{17,20},
	{22,23},{23,18},{18,22}
};


namespace symmetry {
	std::vector<int> perm_from_cycles(
			const std::vector< std::vector<int>>& cycles);

	// CCW rotation by 2pi/3
	const std::vector< std::vector<int>> C3_cycles = {
		{1, 3, 5}, {2, 4, 6}, {7,8,9}, {10, 12, 14}, {11, 13, 15}, 
		{16, 17, 18}, {19, 21, 23}, {20, 22, 24}
	};

	// mirror plane in the 8-0-16 plane
	const std::vector< std::vector<int>> sigma_cycles = {
		{1,6}, {3,4}, {2,5}, {7,9}, {10,15}, {11, 14}, {12, 13},
		{19, 24}, {17,18}, {20,23}, {21,22}
	};


	// inversion
	const std::vector< std::vector<int>> I_cycles = {
		{1,4}, {2,5}, {3,6}, {7, 18}, {8, 16}, {9,17},
		{10,22}, {11, 23}, {12, 24}, {13, 19}, {14, 20}, {15, 21}
	};

	std::vector<int> C3();
	std::vector<int> sigma();
	std::vector<int> I();

	static std::map< std::string, std::vector<xdiag::complex> > irreps_D3d = {
		{"A1g", { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
		{"A2g", { 1, 1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1}},
		{"Eg",  { 2,-1,-1, 0, 0, 0, 2,-1,-1, 0, 0, 0}},
		{"A1u", { 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1}},
		{"A2u", { 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1}},
		{"Eu",  { 2,-1,-1, 0, 0, 0,-2, 1, 1, 0, 0, 0}}
	};

};
