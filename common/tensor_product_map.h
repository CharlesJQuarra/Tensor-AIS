
#ifndef _TENSOR_MERGE_MAP_H_
#define _TENSOR_MERGE_MAP_H_

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

#include <sparse_arr.h>
#include <vector>


template<typename ValueTypeA, unsigned int RankA, typename ValueTypeB, unsigned int RankB, typename ValueTypeP, unsigned int RankP>
struct TensorProductMap
{

	typedef std::pair< unsigned int , unsigned int > index_contraction_t;
	typedef typename MD<ValueTypeA, RankA>::RegArray RegArrayA_t;
	typedef typename MD<ValueTypeB, RankB>::RegArray RegArrayB_t;

	typedef typename RegArrayA_t::Iterator IteratorA_t;
	typedef typename RegArrayA_t::DetailIterator DetailIteratorA_t;
	typedef typename RegArrayB_t::Iterator IteratorB_t;
	typedef typename RegArrayB_t::DetailIterator DetailIteratorB_t;

	/*
		- get list of contractions in A and B, sorted by rank order
		- recurse on A, if rank is in list of A contractions, append iterator to a pair. If A is of rank Na, then a field in the k-column, will be
		- an iterator in MD< N - k >::RegArray::DetailIterator
		- but these are determined at compile-time, so I'm not able to construct a pair with contraction iterators
		- but assuming we have some way to refer to those iterators:
		- for every pair of contracted indices
			- take iterators of pair.A and pair.B
			- for each free index in A, with it.A at start, each free index in B, it.B at start
				- construct index expression
				- compute map
				- loop over all it.A == it.B
			- loop over all contracted indices recursively. With it2.A++ == it2.B++, recurse again over all it.A == it.B
				

		something like this would have to be 

			//
			node_map_recurse()
				if Rank belongs to some it.A with it.A == it.B of contracted indices
					//we need to set common value it.A == it.B according to current label
					
				else 
					// we assume this is a free index, build the index expression

			vector< [contracted_pair_0.A, contracted_pair_0.B], [contracted_pair_1.A, contracted_pair_1.B] ... >
	    
			vector< common value of pair N >
			for All free indices in A, contracted indices in vector< common value of pair N> 
				for All free indices in B, contracted indices in vector< common value of pair N>
					map(A, B)

			for a sparse multidimensional array, the set of common values depends on the preceeding indices, both free and contracted.
			What we do is that we begin the contracted indices for the first item in the pair at the beginning of the array for current preceeding indices.
			when we reach the the dimension of the second item in the pair, we try to find a matching index, or return with no match
				
			
	*/

	TensorProductMap(const RegArrayA_t& tA, const ContractedIndexValues_t& contracted_indices, const RegArrayB_t& tB, typename MD<ValueTypeP, RankP>::RegArray& tP)
	{
		//assert that RankA + RankB == RankP - 2*contracted_indices.size()
		tP.clear();
		ContractedIndexValues_t contracted_index_values;
		auto ff_A = [&](ValueTypeA vA, const int* baseA, const unsigned int index_of_last_update_A)->void
		{
			auto indices_A_total = _IndexFromRange(baseA, RankA);
			std::deque< unsigned int > indices_A;

			//skip contracted indices
			for (unsigned int iA = 0, iAe = indices_A_total.size(); iA < iAe; iA++)
			{
				if (contracted_indices.find(iA) == contracted_indices.end())
				{
					indices_A.push_back( indices_A_total[iA] );
				}
			}

			auto ff_B = [&](ValueTypeB vB, const int* baseB, const unsigned int index_of_last_update_B)->void
			{
				auto indices_B_total = _IndexFromRange(baseB, RankB);
				std::deque< unsigned int > indices_B;

				//skip contracted indices
				for (unsigned int iB = 0, iBe = indices_B_total.size(); iB < iBe; iB++)
				{
					if (contracted_index_values.find(iB) == contracted_index_values.end())
					{
						indices_B.push_back( indices_B_total[iB] );
					}
				}
				auto indices = indices_A;
				indices.insert( indices.end() , indices_B.begin(), indices_B.end() );
				auto indices_read = indices;
				ValueTypeP v = tP( std::move(indices_read) );
				v += vA * vB;
				/*
				std::cout << " tP[";
				for (int iii=0 ; iii < indices.size(); iii++)
					std::cout << indices[iii] << ",";
				std::cout << "] = " << v << std::endl;*/
				tP.assign( std::move(indices), v );
			};

			tB.skip_to_contract_map(contracted_index_values, ff_B);
		};
		
		tA.contract_index_map(contracted_indices, contracted_index_values, ff_A);
	};
};


template<typename ValueTypeA, unsigned int RankA, typename ValueTypeB, unsigned int RankB, typename ValueTypeP>
struct TensorProductMap<ValueTypeA, RankA, ValueTypeB, RankB, ValueTypeP, 0>
{

	typedef std::pair< unsigned int , unsigned int > index_contraction_t;
	typedef typename MD<ValueTypeA, RankA>::RegArray RegArrayA_t;
	typedef typename MD<ValueTypeB, RankB>::RegArray RegArrayB_t;

	typedef typename RegArrayA_t::Iterator IteratorA_t;
	typedef typename RegArrayA_t::DetailIterator DetailIteratorA_t;
	typedef typename RegArrayB_t::Iterator IteratorB_t;
	typedef typename RegArrayB_t::DetailIterator DetailIteratorB_t;

	TensorProductMap(const RegArrayA_t& tA, const ContractedIndexValues_t& contracted_indices, const RegArrayB_t& tB, ValueTypeP& tP)
	{
		//assert that RankA + RankB == 2*contracted_indices.size()
		tP = 0.0;
		ContractedIndexValues_t contracted_index_values;
		auto ff_A = [&](ValueTypeA vA, const int* baseA, const unsigned int index_of_last_update_A)->void
		{

			auto ff_B = [&](ValueTypeB vB, const int* baseB, const unsigned int index_of_last_update_B)->void
			{
				//ValueTypeP v = tP( std::move(indices_read) );
				tP += vA * vB;
				//tP.assign( std::move(indices), v );
			};

			tB.skip_to_contract_map(contracted_index_values, ff_B);
		};
		
		tA.contract_index_map(contracted_indices, contracted_index_values, ff_A);
	};
};




template<typename ValueTypeA, typename ValueTypeB, unsigned int RankB, typename ValueTypeP>
struct TensorProductMap<ValueTypeA, 0, ValueTypeB, RankB, ValueTypeP, RankB>
{

	typedef typename MD<ValueTypeB, RankB>::RegArray RegArrayB_t;

	typedef typename RegArrayB_t::Iterator IteratorB_t;
	typedef typename RegArrayB_t::DetailIterator DetailIteratorB_t;


	TensorProductMap(const ValueTypeA& tA, const ContractedIndexValues_t& contracted_indices, const RegArrayB_t& tB, typename MD<ValueTypeP, RankB>::RegArray& tP)
	{
		//assert that RankA + RankB == RankP - 2*contracted_indices.size()
		tP.clear();
		auto ff_B = [&](float vB, const int* baseB, const unsigned int index_of_last_update_B)->void
		{
			auto indices_B = _IndexFromRange(baseB, RankB);
			tP.assign( std::move(indices_B), tA*vB );
		};

		tB.map(ff_B);
	};
};


#endif


