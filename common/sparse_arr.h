
#ifndef _SPARSE_ARR_H_
#define _SPARSE_ARR_H_

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

#include <boost/unordered_map.hpp>
#include <map>
#include <set>
#include <deque>
#include <iostream>
#include <functional>
#include <utility>

template <typename StoreType>
struct SparseArray
{
	//typedef std::map< unsigned int, StoreType > map_t;
	typedef boost::unordered_map< unsigned int, StoreType > map_t;
	typedef typename map_t::const_iterator iterator;
	map_t map;
};

std::deque<unsigned int> _Index(std::initializer_list<unsigned int>&& vs) {
	std::deque<unsigned int> dq;
	for (auto i = vs.begin(); i != vs.end(); i++)
	{
      dq.push_back(*i);
  }
	return dq;
};

std::deque<unsigned int> _IndexFromRange(const int* base, unsigned int size) {
	std::deque<unsigned int> dq;
	for (int i = 0; i < size; i++)
	{
      dq.push_back(base[i]);
  }
	return dq;
};

std::deque<unsigned int> _IndexFromRangeDebug(const int* base, unsigned int size) {
	std::deque<unsigned int> dq;
	std::cout << "IndexFromRangeDebug: size " << size << std::endl;
	for (int i = 0; i < size; i++)
	{
      dq.push_back(base[i]);
			std::cout << "IndexFromRangeDebug: base[ " << i << " ] = " << base[i] << std::endl;
  }
	for (int i = 0; i < size; i++)
	{
			std::cout << "IndexFromRangeDebug: dq[ " << i << " ] = " << dq[i] << std::endl;
  }
	return dq;
};

typedef std::map< unsigned int, unsigned int > ContractedIndexValues_t;

template < typename ValueType, unsigned int Rank>
struct MD
{
	static constexpr unsigned int rank = Rank;

	struct RegArray
	{
		ValueType DefValue;

		RegArray(ValueType _d) : DefValue(_d), dataset() {}
		RegArray() : DefValue(0.0), dataset() {}

		typedef typename MD< ValueType , Rank-1 >::RegArray RegArray_t;
		SparseArray< RegArray_t > dataset;

		struct DetailIterator
		{
			typename SparseArray< RegArray_t >::iterator current_iterator;
			typename RegArray_t::DetailIterator sub_iterator;
		};

		struct Iterator
		{
			int Indices[Rank];
			unsigned int index_of_last_update;
			DetailIterator detail_iterator;

			std::deque<unsigned int> getIndex() const
			{
				std::deque<unsigned int> dq;
				for (unsigned int i = 0; i < Rank; i++)
				{
						dq.push_back(Indices[i]);
				}
				return dq;
			}
		};

		void clear()
		{
			dataset.map.clear();
		}

		void begin(Iterator& b) const
		{
			b.index_of_last_update = 0;
			begin_detail( b.detail_iterator , b.Indices , 0 );
		}

		bool begin_detail(DetailIterator& iterator, int* b, unsigned int depth) const
		{
			//int* b = iterator.Indices + depth;
			iterator.current_iterator = dataset.map.begin();
			if (iterator.current_iterator != dataset.map.end())
			{
				b[0] = iterator.current_iterator->first;
				return iterator.current_iterator->second.begin_detail( iterator.sub_iterator, b+1 , depth+1 );
			}
			else
				return false;
		}

		bool end(Iterator& e) const
		{
			return end_detail(e.detail_iterator);
		}

		bool end_detail(DetailIterator& e) const
		{
			return (e.current_iterator == dataset.map.end()) || (e.current_iterator->second.end_detail(e.sub_iterator));
		}

		bool next(Iterator& iter) const
		{
			iter.index_of_last_update = Rank;
			return next_detail( iter.detail_iterator, iter.Indices, iter.index_of_last_update , 0 );
		}

		void inspect() const
		{
			std::cout << " inspect : Rank: " << Rank << std::endl;
			auto fd = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
			{
				std::cout << " inspect : iterator index: it[";
				for (int i=0; i< Rank; i++)
				{
					std::cout << base[i] << ",";
				}
				std::cout  << "] = " << v << std::endl;
			};

			map(fd);
		}
		
		void map(std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) > f) const
		{
			Iterator it;
			it.index_of_last_update = 0;
			map_recurse(it.detail_iterator, it.Indices, 0, it.index_of_last_update, f);
		}

		void map_recurse(DetailIterator& iterator, int* base, unsigned int depth, unsigned int& index_of_last_update, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) >& f) const
		{
			auto it_e = dataset.map.end();
			for( iterator.current_iterator = dataset.map.begin(); iterator.current_iterator != it_e; iterator.current_iterator++)
			{
				//std::cout << " node map_recurse: depth -> " << depth << " index of last update " << index_of_last_update << std::endl;
				base[depth] = iterator.current_iterator->first;
				iterator.current_iterator->second.map_recurse( iterator.sub_iterator, base, depth+1, index_of_last_update, f);
				index_of_last_update = depth;
			}
		}

		void skip_to_contract_map(const ContractedIndexValues_t& contracted_index_values, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) > f) const
		{
			Iterator it;
			it.index_of_last_update = 0;
			map_recurse_and_skip_to(contracted_index_values, it.detail_iterator, it.Indices, 0, it.index_of_last_update, f);
		}

		void contract_index_map(const ContractedIndexValues_t& contracted_indices, ContractedIndexValues_t& contracted_index_values, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) > f) const
		{
			Iterator it;
			it.index_of_last_update = 0;
			map_recurse_contract_indices(contracted_indices, contracted_index_values, it.detail_iterator, it.Indices, 0, it.index_of_last_update, f);
		}
		
		void map_recurse_and_skip_to(const ContractedIndexValues_t& contracted_index_values, DetailIterator& iterator, int* base, unsigned int depth, unsigned int& index_of_last_update, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) >& f) const
		{
			auto it_e = dataset.map.end();
			auto it_contracted = contracted_index_values.find(depth);
			if (it_contracted != contracted_index_values.end())
			{
				//it is a contracted index
				unsigned int contracted_value = it_contracted->second;
				iterator.current_iterator = dataset.map.find(contracted_value);

				if (iterator.current_iterator != it_e)
				{
					base[depth] = contracted_value;
					iterator.current_iterator->second.map_recurse_and_skip_to( contracted_index_values, iterator.sub_iterator, base, depth+1, index_of_last_update, f);
					index_of_last_update = depth;
				}
				// if iterator is above contracted_value, then we skip further recursion
			}
			else
			{
				//it is a free index
				for( iterator.current_iterator = dataset.map.begin(); iterator.current_iterator != it_e; iterator.current_iterator++)
				{
					//std::cout << " node map_recurse: depth -> " << depth << " index of last update " << index_of_last_update << std::endl;
					base[depth] = iterator.current_iterator->first;
					iterator.current_iterator->second.map_recurse_and_skip_to( contracted_index_values, iterator.sub_iterator, base, depth+1, index_of_last_update, f);
					index_of_last_update = depth;
				}
			}
		}

		void map_recurse_contract_indices(const ContractedIndexValues_t& contracted_indices, ContractedIndexValues_t& contracted_index_values, DetailIterator& iterator, int* base, unsigned int depth, unsigned int& index_of_last_update, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) >& f) const
		{
			auto it_e = dataset.map.end();
			auto it_must_contract = contracted_indices.find(depth);
			bool must_be_contracted = (it_must_contract != contracted_indices.end());
			for( iterator.current_iterator = dataset.map.begin(); iterator.current_iterator != it_e; iterator.current_iterator++)
			{
				//std::cout << " node map_recurse: depth -> " << depth << " index of last update " << index_of_last_update << std::endl;
				if (must_be_contracted)
				{
					// the value must be added to the contracted values
					contracted_index_values[ it_must_contract->second ] = iterator.current_iterator->first;
				}
				base[depth] = iterator.current_iterator->first;
				//f( iterator.current_iterator->second , base , index_of_last_update );
				iterator.current_iterator->second.map_recurse_contract_indices( contracted_indices, contracted_index_values, iterator.sub_iterator, base, depth+1, index_of_last_update, f);
				index_of_last_update = depth;
			}
		}

		bool next_detail(DetailIterator& iterator, int* iter, unsigned int& index_of_last_update, unsigned int depth) const
		{
			//int* iter = iterator.Indices + depth;
			if (iterator.current_iterator->second.next_detail( iterator.sub_iterator , iter+1, index_of_last_update, depth+1 ))
			{
				return true;
			}
			do
			{
				iterator.current_iterator++;
				if (iterator.current_iterator == dataset.map.end())
				{
					iter[0] = -1;
					return false;
				}
				iter[0] = iterator.current_iterator->first;
				if (index_of_last_update > depth)
					index_of_last_update = depth;
			} while (! iterator.current_iterator->second.begin_detail( iterator.sub_iterator , iter+1 , depth+1 ));
			// - o -
		}

		
		
		/*
		template <int N, int... Rest>
		ValueType entry() 
		{
			auto it = dataset.map.find(N);
			if (it == dataset.map.end())
				return DefValue;
			return it->second.entry<Rest...>();
		}*/

		ValueType operator()(std::deque< unsigned int >&& idxs) const
		{
			unsigned int ii = idxs.front();
			idxs.pop_front();
			auto it = dataset.map.find(ii);
			if (it == dataset.map.end())
				return DefValue;
			return it->second(std::move(idxs));
		}

		bool empty() const
		{
			return dataset.map.empty();
		}

		RegArray& remove_entry(std::deque< unsigned int >&& idxs)
		{
			unsigned int ii = idxs.front();
			idxs.pop_front();
			auto it = dataset.map.find(ii);
			if (it == dataset.map.end())
			{
				return *this;
			}
			it->second.remove_entry(std::move(idxs));
			if (it->second.empty())
				dataset.map.erase(it);
			return *this;
		}

		RegArray& assign(std::deque< unsigned int >&& idxs, ValueType v)
		{
			if (v == DefValue)
				return remove_entry( std::move(idxs) );
			return assign_detail(std::move(idxs), v);
		}

		RegArray& assign_detail(std::deque< unsigned int >&& idxs, ValueType v)
		{
			if (v == DefValue)
				return remove_entry( std::move(idxs) );
			unsigned int ii = idxs.front();
			idxs.pop_front();
			auto it = dataset.map.find(ii);
			if (it == dataset.map.end())
			{
				RegArray_t arr(DefValue);
				dataset.map[ii] = arr;
				dataset.map[ii].assign_detail(std::move(idxs), v);
				return *this;
			}
			it->second.assign_detail(std::move(idxs), v);
			return *this;
		}

		RegArray& operator +=( const RegArray& other)
		{
			bool y_value_has_been_loaded = false;
			float current_y;
			std::deque< unsigned int > current_y_idx;
			std::deque< unsigned int > temp;
			auto ff = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
			{
				if (y_value_has_been_loaded)
				{
					this->assign( std::move(current_y_idx) , current_y );
				}
				current_y_idx = _IndexFromRange(base, Rank);
				temp = current_y_idx;
				current_y =	other( std::move(temp) );
				y_value_has_been_loaded = true;
				
				current_y += v;
			};

			this->map(ff);

			if (y_value_has_been_loaded)
			{
				this->assign( std::move(current_y_idx) , current_y );
			}
			return *this;
		}

		RegArray& operator -=( const RegArray& other)
		{
			bool y_value_has_been_loaded = false;
			float current_y;
			std::deque< unsigned int > current_y_idx;
			std::deque< unsigned int > temp;
			auto ff = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
			{
				if (y_value_has_been_loaded)
				{
					this->assign( std::move(current_y_idx) , current_y );
				}
				current_y_idx = _IndexFromRange(base, Rank);
				temp = current_y_idx;
				current_y =	v - other( std::move(temp) );
				y_value_has_been_loaded = true;
			};

			this->map(ff);

			if (y_value_has_been_loaded)
			{
				this->assign( std::move(current_y_idx) , current_y );
			}
			return *this;
		}

		typedef std::pair< unsigned int , unsigned int > index_contraction;
		template<typename ValueTypeA, unsigned int RankA, typename ValueTypeB, unsigned int RankB, typename ValueTypeP, unsigned int RankP>
		friend class TensorMergeMap;

		template <unsigned int RankOther, unsigned int RankProduct>
		void mult( const typename MD<ValueType, RankOther>::RegArray& other, typename MD<ValueType, RankProduct>::RegArray& product, const std::vector< index_contraction >& v_idx_contract) const
		{
			//assert RankProduct == Rank + RankOther - 2*v_idx_contract.size()
		}
		
	};
};






//--------------------------


template <typename ValueType>
struct MD<ValueType, 1>
{

	static constexpr unsigned int rank = 1;

	struct RegArray
	{

		ValueType DefValue;

		RegArray(ValueType _d) : DefValue(_d), value_dataset() {}
		RegArray() : DefValue(0.0), value_dataset() {}

		SparseArray< ValueType > value_dataset;
	
		ValueType operator()(std::deque< unsigned int >&& idxs) const
		{
			unsigned int ii = idxs.front();
			idxs.pop_front();
			auto it = value_dataset.map.find(ii);
			if (it == value_dataset.map.end())
				return DefValue;
			return it->second;
		}

		void inspect() const
		{
			std::cout << " inspect : Rank: " << 1 << std::endl;
			auto fd = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
			{
				std::cout << " inspect : iterator index: it[";
				for (int i=0; i< 1; i++)
				{
					std::cout << base[i] << ",";
				}
				std::cout  << "] = " << v << std::endl;
			};

			map(fd);
		}

		bool empty() const
		{
			return value_dataset.map.empty();
		}

		void clear()
		{
			value_dataset.map.clear();
		}

		RegArray& remove_entry(std::deque< unsigned int >&& idxs)
		{
			unsigned int ii = idxs.front();
			idxs.pop_front();
			auto it = value_dataset.map.find(ii);
			if (it == value_dataset.map.end())
			{
				return *this;
			}
			value_dataset.map.erase(it);
			return *this;
		}

		RegArray& assign(std::deque< unsigned int >&& idxs, ValueType v)
		{
			if (v == DefValue)
				return remove_entry( std::move(idxs) );
			return assign_detail(std::move(idxs), v);
		}

		RegArray& assign_detail(std::deque< unsigned int >&& idxs, ValueType v)
		{
			unsigned int ii = idxs.front();
			idxs.pop_front();
			auto it = value_dataset.map.find(ii);
			if (it == value_dataset.map.end())
			{
				value_dataset.map[ii] = v;
				return *this;
			}
			it->second = v;
			return *this;
		}

		// Iterable
		struct DetailIterator
		{
			typename SparseArray< ValueType >::iterator current_iterator;
		};

		struct Iterator
		{
			int Indices[1];
			DetailIterator detail_iterator;
			unsigned int index_of_last_update;

			std::deque<unsigned int>&& getIndex() const
			{
				std::deque<unsigned int> dq;
				dq.push_back(Indices[0]);
				return std::forward< std::deque<unsigned int> >(dq);
			}
		};

		void begin(Iterator& b) const
		{
			b.index_of_last_update = 1;
			begin_detail( b.detail_iterator , b.Indices , 0 );
		}

		bool begin_detail(DetailIterator& iterator, int* b, unsigned int depth) const
		{
			//int* b = iterator.Indices + depth;
			iterator.current_iterator = value_dataset.map.begin();
			if (iterator.current_iterator != value_dataset.map.end())
			{
				b[0] = iterator.current_iterator->first;
				std::cout << " leaf begin detail: depth -> " << depth << " current iterator -> " << b[0] << std::endl;
				return true;
			}
			else
				return false;
		}

		bool end(Iterator& e) const
		{
			return end_detail(e.detail_iterator);
		}

		bool end_detail(DetailIterator& e) const
		{
			return (e.current_iterator == value_dataset.map.end());
		}

		void map(std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) > f) const
		{
			Iterator it;
			it.index_of_last_update = 0;
			map_recurse(it.detail_iterator, it.Indices, 0, it.index_of_last_update, f);
		}

		void map_recurse(DetailIterator& iterator, int* base, unsigned int depth, unsigned int& index_of_last_update, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) >& f) const
		{
			auto it_e = value_dataset.map.end();
			for( iterator.current_iterator = value_dataset.map.begin(); iterator.current_iterator != it_e; iterator.current_iterator++)
			{
				//iterator.current_iterator->second.map_recurse( iterator.sub_iterator, base, depth+1, index_of_last_update, f);
				//std::cout << " leaf map_recurse: depth -> " << depth << " index of last update " << index_of_last_update << std::endl;
				base[depth] = iterator.current_iterator->first;
				f( iterator.current_iterator->second , base , index_of_last_update );
				index_of_last_update = depth;
			}
		}

		void skip_to_contract_map(const ContractedIndexValues_t& contracted_index_values, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) > f) const
		{
			Iterator it;
			it.index_of_last_update = 0;
			map_recurse_and_skip_to(contracted_index_values, it.detail_iterator, it.Indices, 0, it.index_of_last_update, f);
		}

		void contract_index_map(const ContractedIndexValues_t& contracted_indices, ContractedIndexValues_t& contracted_index_values, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) > f) const
		{
			Iterator it;
			it.index_of_last_update = 0;
			map_recurse_contract_indices(contracted_indices, contracted_index_values, it.detail_iterator, it.Indices, 0, it.index_of_last_update, f);
		}

		void map_recurse_and_skip_to(const ContractedIndexValues_t& contracted_index_values, DetailIterator& iterator, int* base, unsigned int depth, unsigned int& index_of_last_update, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) >& f) const
		{
			auto it_e = value_dataset.map.end();
			auto it_contracted = contracted_index_values.find(depth);
			if (it_contracted != contracted_index_values.end())
			{
				//it is a contracted index
				unsigned int contracted_value = it_contracted->second;
				iterator.current_iterator = value_dataset.map.find(contracted_value);

				if (iterator.current_iterator != it_e)
				{
					base[depth] = contracted_value;
					f( iterator.current_iterator->second , base , index_of_last_update );
					index_of_last_update = depth;
				}
			}
			else
			{
				//it is a free index
				for( iterator.current_iterator = value_dataset.map.begin(); iterator.current_iterator != it_e; iterator.current_iterator++)
				{
					//std::cout << " node map_recurse: depth -> " << depth << " index of last update " << index_of_last_update << std::endl;
					base[depth] = iterator.current_iterator->first;
					f( iterator.current_iterator->second , base , index_of_last_update );
					index_of_last_update = depth;
				}
			}
		}

		void map_recurse_contract_indices(const ContractedIndexValues_t& contracted_indices, ContractedIndexValues_t& contracted_index_values, DetailIterator& iterator, int* base, unsigned int depth, unsigned int& index_of_last_update, std::function< void(ValueType , const int* /*Indices array*/, const unsigned int /*index of last update*/) >& f) const
		{
			auto it_e = value_dataset.map.end();
			auto it_must_contract = contracted_indices.find(depth);
			for( iterator.current_iterator = value_dataset.map.begin(); iterator.current_iterator != it_e; iterator.current_iterator++)
			{
				//std::cout << " node map_recurse: depth -> " << depth << " index of last update " << index_of_last_update << std::endl;
				if (it_must_contract != contracted_indices.end())
				{
					// the value must be added to the contracted values
					contracted_index_values[ it_must_contract->second ] = iterator.current_iterator->first;
				}
				base[depth] = iterator.current_iterator->first;
				f( iterator.current_iterator->second , base , index_of_last_update );
				index_of_last_update = depth;
			}
		}

		bool next(Iterator& iter) const
		{
			iter.index_of_last_update = 1;
			return next_detail( iter.detail_iterator, iter.Indices, iter.index_of_last_update , 0 );
		}

		bool next_detail(DetailIterator& iterator, int* iter, unsigned int& index_of_last_update, unsigned int depth) const
		{
			//int* iter = iterator.Indices + depth;
			iterator.current_iterator++;
			if (iterator.current_iterator == value_dataset.map.end())
			{
				iter[0] = -1;
				return false;
			}
			iter[0] = iterator.current_iterator->first;
			if (index_of_last_update > depth)
				index_of_last_update = depth;
			std::cout << " leaf next detail: depth -> " << depth << " current iterator -> " << iter[0] << std::endl;
			return true;
			// - o -
		}

		template<typename ValueTypeA, unsigned int RankA, typename ValueTypeB, unsigned int RankB, typename ValueTypeP, unsigned int RankP>
		friend class TensorMergeMap;

		RegArray& operator +=( const RegArray& other)
		{
			bool y_value_has_been_loaded = false;
			ValueType current_y;
			std::deque< unsigned int > current_y_idx;
			std::deque< unsigned int > temp;
			auto ff = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
			{
				if (y_value_has_been_loaded)
				{
					this->assign( std::move(current_y_idx) , current_y );
				}
				current_y_idx = _IndexFromRange(base, 1);
				temp = current_y_idx;
				current_y =	other( std::move(temp) );
				y_value_has_been_loaded = true;
				
				current_y += v;
			};

			this->map(ff);

			if (y_value_has_been_loaded)
			{
				this->assign( std::move(current_y_idx) , current_y );
			}
			return *this;
		}

		RegArray& operator -=( const RegArray& other)
		{
			bool y_value_has_been_loaded = false;
			ValueType current_y;
			std::deque< unsigned int > current_y_idx;
			std::deque< unsigned int > temp;
			auto ff = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
			{
				current_y_idx = _IndexFromRange(base, 1);
				temp = current_y_idx;
				current_y =	v - other( std::move(temp) );
				this->assign( std::move(current_y_idx) , current_y );
			};

			this->map(ff);
			return *this;
		}


	};
};





#endif
