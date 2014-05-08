

#ifndef _POLY_H_
#define _POLY_H_

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




template<typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY>
struct BigDataPoly
{
	typename MD<ValueType, RankY + Order*RankX >::RegArray arr;
	BigDataPoly< ValueType, Order-1 , RankX , RankY > sub_poly;

	void eval( const typename MD<ValueType, RankX>::RegArray& x_test, typename MD<ValueType, RankY>::RegArray& y_test) const;
	void derivative_coeff(BigDataPoly< ValueType, Order-1, RankX, RankY+RankX >& diff) const;

	template <typename Loader>
	void load_data(Loader& loader)
	{
		arr.clear();
		sub_poly.load_data(loader.subloader);
		loader.load_array(arr);
	}

	void inspect()
	{
		arr.inspect();
		sub_poly.inspect();
	}

	static constexpr unsigned int ord = Order;
	static constexpr unsigned int rank_x = RankX;
	static constexpr unsigned int rank_y = RankY;
};


template<typename ValueType, unsigned int RankX, unsigned int RankY>
struct BigDataPoly<ValueType, 0, RankX, RankY>
{
	typename MD<ValueType, RankY>::RegArray arr;
	static constexpr unsigned int ord = 0;
	static constexpr unsigned int rank_x = RankX;
	static constexpr unsigned int rank_y = RankY;

	void inspect()
	{
		arr.inspect();
	}

	void eval( const typename MD<ValueType, RankX>::RegArray& x_test, typename MD<ValueType, RankY>::RegArray& y_test) const
	{
		y_test.clear();
		//y_test += arr;
		auto ff = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
		{
			ValueType current_y =	y_test( _IndexFromRange(base, RankY) );
			current_y += v;
			y_test.assign( _IndexFromRange(base, RankY) , current_y );
		};
		arr.map(ff);
	}

	template <typename Loader>
	void load_data(Loader& loader)
	{
		arr.clear();
		loader.load_array(arr);
	}

};

template<typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY>
void BigDataPoly<ValueType, Order, RankX, RankY>::eval( const typename MD<ValueType, RankX>::RegArray& x_test, typename MD<ValueType, RankY>::RegArray& y_test) const
{
	sub_poly.eval( x_test, y_test );
	/*
	auto displ_y = [&](float v, const int* base, const unsigned int index_of_last_update)->void
	{
		std::cout << "displ_Y lambda: arr value: y_test = " << v << " Order -> " << Order << std::endl;
	};
	y_test.map( displ_y );*/

	// y_test_{i_free} += arr_{{i_free} i0 i1 ... iOrder} * x_{i0} * x_{i1} * ... * x_{iOrder}
	bool y_value_has_been_loaded = false;
	ValueType current_y;
	std::deque< unsigned int > current_y_idx;
	std::deque< unsigned int > temp;
	std::vector< ValueType > current_x_powers;
	auto ff = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
	{
		//std::cout << "EVAL lambda: arr value: " << v << " Order -> " << Order << std::endl;
		int x_index_base = (int)index_of_last_update - (int)RankY;
		//std::cout << "EVAL lambda: x_index_base.. " << x_index_base << std::endl;
		if (x_index_base < 0)
		{
			//y value must be loaded. If a value is already loaded, it should be stored before loading a new one.
			if (y_value_has_been_loaded)
			{
				y_test.assign( std::move(current_y_idx) , current_y );
			}
			current_y_idx = _IndexFromRange(base, RankY);
			temp = current_y_idx;
			current_y =	y_test( std::move(temp) );
			y_value_has_been_loaded = true;
			//std::cout << "EVAL lambda: loading y_value " << current_y << std::endl;
		}
		// x updates
		int x_power_must_be_updated = (x_index_base / (int)Order);
		unsigned int new_power_index = x_power_must_be_updated < 0 ? 0 : x_power_must_be_updated;
		/*std::cout << "EVAL lambda: x_power index.. " << x_power_must_be_updated << " , " << new_power_index << std::endl;
		for (int ii = 0; ii < current_x_powers.size() ; ii++)
		{
			std::cout << "EVAL lambda: x_power index array[ " << ii << " ] = " << current_x_powers[ii] << std::endl;
		}*/
		current_x_powers.resize(new_power_index);
		ValueType prev = (new_power_index > 0) ? current_x_powers[ new_power_index - 1 ] : ValueType(1.0);
		std::deque< unsigned int > current_x_idx;
		for (unsigned int i=new_power_index; i < (Order - 1) ; i++)
		{
			current_x_idx = _IndexFromRange(base + RankY + RankX*i, RankX);
			prev *= x_test( std::move(current_x_idx) );
			current_x_powers.push_back( prev );
		}
		current_x_idx = _IndexFromRange(base + RankY + RankX*(Order - 1), RankX);
		prev *= x_test( std::move( current_x_idx ) );
		
		current_y += v * prev;
	};

	arr.map(ff);

	if (y_value_has_been_loaded)
	{
		y_test.assign( std::move(current_y_idx) , current_y );
	}
  // - 0 -
}

template<typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY>
void BigDataPoly<ValueType, Order, RankX, RankY>::derivative_coeff(BigDataPoly<ValueType, Order-1, RankX, RankY+RankX >& diff) const
{

	/* Diff( y_{i_free}, x_j ) =
 		 Diff( arr_{{i_free} i0 i1 ... iOrder} * x_{i0} * x_{i1} * ... * x_{iOrder} , x_j)
	 = arr_{{i_free} j i1 ... iOrder} * x_{i1} * ... * x_{iOrder}
		+  arr_{{i_free} i0 j ... iOrder} * x_{i0} * ... * x_{iOrder}
		+ .... + arr_{{i_free} i0 i1 ... j} * x_{i0} * x_{i1} * ... * x_{iOrder-1}

		Diff_{{i_free} j i0 i1 ... iOrder-1} = arr{{i_free} j i0 i1 ... iOrder-1} + arr{{i_free} i0 j i1 ... i_{Order -1} } + .... + arr{{i_free} i0 i1 ... j}
	*/
	diff.arr.clear();
	auto ff = [&](ValueType v, const int* base, const unsigned int index_of_last_update)->void
	{
		std::deque< unsigned int > indices = _IndexFromRangeDebug(base, RankY + Order*RankX);
		auto free_indices = _IndexFromRangeDebug(base, RankY);
		auto contract_poly_indices = _IndexFromRangeDebug(base + RankY, Order*RankX);
		auto write_indices = indices;
		
		std::cout << " BigDataPoly<Order = "<< Order << ", RankX = "<< RankX << ", RankY = "<< RankY <<">::derivative_coeff. Before reordering: "<< std::endl;
		for (int ii=0; ii < indices.size(); ii++)
		{
			std :: cout << indices[ii] << ",";
		}
		std::cout << std::endl;
		std::cout << " free indices: ";
		for (int ii=0; ii < free_indices.size(); ii++)
		{
			std :: cout << free_indices[ii] << ",";
		}
		std::cout << std::endl;
		diff.arr.inspect();
		
		float diff_v = diff.arr( std::move(indices) ) + v;
		diff.arr.assign( std::move(write_indices), diff_v );
		diff.arr.inspect();
		for (unsigned int ii = RankX; ii < Order*RankX; ii++)
		{
			indices.clear();
			indices.insert(indices.end(), free_indices.begin(), free_indices.end());
			auto contract_shifted_indices = contract_poly_indices;
			std::cout << " contract BEFORE shifted indices: ";
			for (int ii=0; ii < contract_shifted_indices.size(); ii++)
			{
				std :: cout << contract_shifted_indices[ii] << ",";
			}
			std::cout << std::endl;
			for (unsigned int shX = 0; shX < RankX; shX++)
			{
				unsigned int index_swap = contract_shifted_indices[ii + shX];
				contract_shifted_indices.erase( contract_shifted_indices.begin() + ii + shX);
				contract_shifted_indices.insert( contract_shifted_indices.begin() + shX , index_swap );
			}
			indices.insert(indices.end(), contract_shifted_indices.begin(), contract_shifted_indices.end());
			write_indices = indices;
			
			std::cout << " BigDataPoly<Order = "<< Order << ", RankX = "<< RankX << ", RankY = "<< RankY <<">::derivative_coeff. AFTER reordering: "<< std::endl;
			for (int kk=0; kk < indices.size(); kk++)
			{
				std :: cout << indices[kk] << ",";
			}
			std::cout << std::endl;
			std::cout << " contract AFTER shifted indices: ";
			for (int ii=0; ii < contract_shifted_indices.size(); ii++)
			{
				std :: cout << contract_shifted_indices[ii] << ",";
			}
			std::cout << std::endl;
			diff.arr.inspect();
			float diff_v = diff.arr( std::move(indices) );
			diff_v += v;
			diff.arr.assign(std::move(write_indices), diff_v);
			diff.arr.inspect();
		}
	};
	arr.map(ff);
}


template<typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY>
struct DerivePoly
{
	BigDataPoly<ValueType, Order, RankX, RankY>& bdp;
	DerivePoly(BigDataPoly<ValueType, Order, RankX, RankY>& _bdp) : bdp(_bdp) {}
	
	void derivative(BigDataPoly< ValueType, Order-1, RankX, RankY+RankX >& diff) const
	{
		bdp.derivative_coeff(diff);
		DerivePoly< ValueType, Order-1, RankX, RankY > sub_derive_poly(bdp.sub_poly);
		sub_derive_poly.derivative( diff.sub_poly );
	}
};

template<typename ValueType, unsigned int RankX, unsigned int RankY>
struct DerivePoly< ValueType, 1, RankX, RankY>
{
	BigDataPoly<ValueType, 1, RankX, RankY>& bdp;
	DerivePoly(BigDataPoly<ValueType, 1, RankX, RankY>& _bdp) : bdp(_bdp) {}
	
	void derivative(BigDataPoly< ValueType, 0, RankX, RankY+RankX >& diff) const
	{
		bdp.derivative_coeff(diff);
	}
};


#endif
