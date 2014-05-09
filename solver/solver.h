
#ifndef _ALGEBRAIC_SOLVER_H
#define _ALGEBRAIC_SOLVER_H

/*
The MIT License (MIT)

Copyright (c) 2014 Charles J. Quarra

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <vector>
#include <poly.h>
#include <tensor_merge_map.h>
#include <iostream>
#include <cmath>



template <typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY>
struct AlgebraicIterativeSolver
{

	BigDataPoly< ValueType, Order , RankX , RankY> system;
	BigDataPoly< ValueType, Order-1, RankX, RankY+RankX > system_diff;
	const typename MD<ValueType, RankY>::RegArray y_goal;
	typename MD<ValueType, RankX>::RegArray* current_x_test;
	typename MD<ValueType, RankY>::RegArray y_test;
	typename MD<ValueType, RankX>::RegArray delta_x;
	typename MD<ValueType, RankX>::RegArray delta_x_2;

	template <typename Loader>
	AlgebraicIterativeSolver(Loader& loader, const typename MD<ValueType, RankY>::RegArray& _y_g) : y_goal(_y_g), system()
	{
		system.load_data(loader);
		compute_system_diff();
	}

	void setCurrentXTest(typename MD<ValueType, RankX>::RegArray& x_test)
	{
		current_x_test = & x_test;
	}

	
	void compute_system_diff()
	{
		DerivePoly< ValueType, Order , RankX , RankY > dp(system);
		dp.derivative(system_diff);
	}

// we compute the form phi = ( Y_test( X_test ) - Y_goal )^2
// and update X_test_j with alpha * 2 * ( Y_test - Y_goal)_i diff( Y_test_i , X_test_j )
// with alpha some small coefficient.

/*
	delta phi = 2 * ( Y_test(X_test) - Y_goal)_i * diff( Y_i, X_j )(X_test) delta_X_j + 
				2 * [ diff( Y_i , X_k )(X_test) * diff( Y_i , X_j )(X_test) + ( Y_test(X_test) - Y_goal)_i * diff^2( Y_i , X_j , X_k )(X_test) ] * delta_X_j * delta_X_k 


	assume that Y_test - Y_goal = dY_0, with all dY_i = 0 with i>0. It means that the delta_X_j components 


*/
	void compute_y_test_delta()
	{
		//std::cout << " compute_y_test_delta before: x = ["<< (*current_x_test)( _Index({0}) ) << "," << (*current_x_test)( _Index({1}) ) << "," << (*current_x_test)( _Index({2}) ) << "]" << std::endl;
		system.eval( *current_x_test , y_test);
		//std::cout << " compute_y_test_delta after: x = ["<< (*current_x_test)( _Index({0}) ) << "," << (*current_x_test)( _Index({1}) ) << "," << (*current_x_test)( _Index({2}) ) << "]" << std::endl;
		y_test -= y_goal;
		//std::cout << "compute_y_test_delta(): system eval solution: y = ["<< y_test( _Index({0}) ) << "," << y_test( _Index({1}) ) << "," << y_test( _Index({2}) ) << "] , dot^2 =  " << pow( y_test( _Index({0}) ) , 2 ) + pow( y_test( _Index({1}) ) , 2 ) + pow( y_test( _Index({2}) ) , 2) << std::endl;
		//std::cout << "compute_y_test_delta(): system goal: y_goal = ["<< y_goal( _Index({0}) ) << "," << y_goal( _Index({1}) ) << "," << y_goal( _Index({2}) ) << "] " << std::endl;

		//std::cout << "compute_y_test_delta(): system eval solution minus goal: y( x = " << (*current_x_test)( _Index({0}) ) << "," << (*current_x_test)( _Index({1}) ) << "," << (*current_x_test)( _Index({2}) ) << " ) = ["<< y_test( _Index({0}) ) << "," << y_test( _Index({1}) ) << "," << y_test( _Index({2}) ) << "]  " << std::endl;
	}

	ValueType compute_phi() const
	{
		ContractedIndexValues_t c_idx;
		for (int i=0; i<RankY; i++)
		{
			c_idx[i] = i;
		}
		ValueType local_phi;
		TensorMergeMap<ValueType, RankY, ValueType, RankY, ValueType, 0> dotProduct(y_test, c_idx, y_test, local_phi);
		//std::cout << "compute_phi: system eval solution: y = ["<< y_test( _Index({0}) ) << "," << y_test( _Index({1}) ) << "," << y_test( _Index({2}) ) << "]" << std::endl;
		//std::cout << "compute_phi: local_phi: " << local_phi << " dot^2 =  " << pow( y_test( _Index({0}) ) , 2 ) + pow( y_test( _Index({1}) ) , 2 ) + pow( y_test( _Index({2}) ) , 2) << std::endl;
		return local_phi;
	}

	void compute_delta_X()
	{
		delta_x.clear();
		typename MD<ValueType, RankY + RankX>::RegArray y_grad_ij(0.0);
		system_diff.eval( *current_x_test, y_grad_ij);
		//std::cout << "compute_delta_X(): y_grad_ij: " << std::endl;
		//y_grad_ij.inspect();
		ContractedIndexValues_t c_idx;
		for (int i=0; i< RankY; i++)
		{
			c_idx[i] = i;
		}
		/*
		{
			delta_x.clear();
			typename MD<ValueType, RankY>::RegArray y_blah(0.0);
			y_blah.assign( _Index({0}), 1.0);
			TensorMergeMap< ValueType, RankY, ValueType, RankY + RankX, ValueType, RankX> dot( y_blah, c_idx, y_grad_ij, delta_x);
			std::cout << " delta_X WOULD BE::..... delta_x = ["<< delta_x( _Index({0}) ) << "," << delta_x( _Index({1}) ) << "," << delta_x( _Index({2}) ) << "]  " << std::endl;
		}*/
		TensorMergeMap< ValueType, RankY, ValueType, RankY + RankX, ValueType, RankX> dot( y_test, c_idx, y_grad_ij, delta_x);
	}

	ValueType loop(ValueType alpha)
	{
		if (current_x_test == 0)
			throw -1; // x_test must be set
		
		compute_y_test_delta();
		double phi = compute_phi();
		compute_delta_X();
		ContractedIndexValues_t c_idx;
		TensorMergeMap< ValueType, 0, ValueType, RankX, ValueType, RankX> scalarMult( 2.0*alpha, c_idx, delta_x, delta_x_2);
		//std::cout << "delta_x_increment(): delta_x_2 = ["<< delta_x_2( _Index({0}) ) << "," << delta_x_2( _Index({1}) ) << "," << delta_x_2( _Index({2}) ) << "]  " << std::endl;
			
		(*current_x_test) -= delta_x_2;

		return phi;
	}

	void line_search(std::map< ValueType, ValueType >& values , ValueType alpha_from, ValueType alpha_to, int elements)
	{
		if (current_x_test == 0)
			throw -1; // x_test must be set

		
		ValueType elem_factor = pow( (alpha_to / alpha_from) , 1.0/elements );
		ValueType alpha = alpha_from;
		compute_y_test_delta();
		//values[ 0.0 ] = compute_phi();
		compute_delta_X();
		ContractedIndexValues_t c_idx;
		while (fabs(alpha) < fabs(alpha_to))
		{
			TensorMergeMap< ValueType, 0, ValueType, RankX, ValueType, RankX> scalarMult( 2*alpha, c_idx, delta_x, delta_x_2);
			//std::cout << "delta_x_increment(): delta_x_2 = ["<< delta_x_2( _Index({0}) ) << "," << delta_x_2( _Index({1}) ) << "," << delta_x_2( _Index({2}) ) << "]  " << std::endl;
			(*current_x_test) -= delta_x_2;
			compute_y_test_delta();
			values[ alpha ] = compute_phi();
			backtrack();
			alpha *= elem_factor;
		}
	}

	void line_search_between(std::map< ValueType, ValueType >& values , typename MD<ValueType, RankX>::RegArray& from, typename MD<ValueType, RankX>::RegArray& to, int elements)
	{
		typename MD<ValueType, RankX>::RegArray* store = current_x_test;
		typename MD<ValueType, RankX>::RegArray current;
		typename MD<ValueType, RankX>::RegArray dist;
		current = to;
		current -= from;
		ContractedIndexValues_t c_idx;
		TensorMergeMap< ValueType, 0, ValueType, RankX, ValueType, RankX> scalarMult( 1.0/elements, c_idx, current, dist);
		current = from;
		/*
		std::cout << " line_search_between: to = ["<< to( _Index({0}) ) << "," << to( _Index({1}) ) << "," << to( _Index({2}) ) << "]" << std::endl;
		std::cout << " line_search_between: from = ["<< from( _Index({0}) ) << "," << from( _Index({1}) ) << "," << from( _Index({2}) ) << "]" << std::endl;
		std::cout << " line_search_between: dist = ["<< dist( _Index({0}) ) << "," << dist( _Index({1}) ) << "," << dist( _Index({2}) ) << "]" << std::endl;
		*/
		setCurrentXTest( current );
		for (int i=0; i<= elements; i++)
		{
			compute_y_test_delta();
			ValueType local_phi = compute_phi();
			values[ (float)i / elements ] = local_phi;
			//std::cout << " line_search_between: x = ["<< current( _Index({0}) ) << "," << current( _Index({1}) ) << "," << current( _Index({2}) ) << "] phi = " << local_phi << std::endl;
			current += dist;
		}
		current_x_test = store;
	}

	ValueType line_search_update()
	{
		std::map< ValueType , ValueType > line_search_results;
		ValueType  lowerP = 0.000000001;
		ValueType factor = 10.0;
		unsigned int retries = 0;
		retry:
		line_search_results.clear();
		line_search( line_search_results, factor*lowerP , factor , 8 );
		//line_search( line_search_results, -factor*lowerP , -factor , 100 );
		ValueType best_phi = 9e30;
		ValueType best_alpha = 0.0;
		for (auto it = line_search_results.begin(), it_e = line_search_results.end(); it != it_e ; it++)
		{
			//std:: cout << " line search update: alpha = " << it->first << " -> phi = " << it->second << std::endl;
			if (it->second < best_phi)
			{
				best_phi = it->second;
				best_alpha = it->first;
			}
		}
		if (best_phi < 9e30)
		{
			ContractedIndexValues_t c_idx;
			TensorMergeMap< ValueType, 0, ValueType, RankX, ValueType, RankX> scalarMult( 2*best_alpha, c_idx, delta_x, delta_x_2);
			(*current_x_test) -= delta_x_2;
			compute_y_test_delta();
			return best_phi;
		}
		
		factor *= 0.1;
		retries++;
		if (retries > 3)
		{
			std::cout << " unable to find a reasonable point. " << std::endl;
			typename MD<ValueType, RankY + RankX>::RegArray y_grad_ij(0.0);
			system_diff.eval( *current_x_test, y_grad_ij);
			//std::cout << "y_grad_ij: " << std::endl;
			//y_grad_ij.inspect();
			for (auto it = line_search_results.begin(), it_e = line_search_results.end(); it != it_e ; it++)
			{
				std:: cout << " line search update: alpha = " << it->first << " -> phi = " << it->second << std::endl;
			}
			return best_phi;
		}
		goto retry;
		
	}

	void backtrack()
	{
		(*current_x_test) += delta_x_2;
	}

};


#endif
