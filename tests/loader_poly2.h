
#ifndef _LOADER_POLY2_H
#define _LOADER_POLY2_H

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

template <typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY, bool MakeNonLinear= false>
struct Poly2LoaderExample
{
	void load_array( typename MD<ValueType, RankY + Order*RankX >::RegArray& arr)
	{
	}

	Poly2LoaderExample< ValueType, Order-1 , RankX , RankY> subloader;
};




template <typename ValueType>
struct Poly2LoaderExample< ValueType, 0 , 1 , 1>
{
	void load_array( typename MD<ValueType, 1 >::RegArray& arr)
	{
		arr.assign( _Index({ 0}) , ValueType(-5.0));
		arr.assign( _Index({ 1}) , ValueType(5.0));
		arr.assign( _Index({ 2}) , ValueType(2.0));
	}

};

template <typename ValueType>
struct Poly2LoaderExample< ValueType, 1 , 1 , 1>
{
	void load_array( typename MD<ValueType, 2 >::RegArray& arr)
	{
		arr.assign( _Index({ 0, 0}) , ValueType(1.0));
		arr.assign( _Index({ 0, 1}) , ValueType(0.5));
		arr.assign( _Index({ 0, 2}) , ValueType(8.0));
	
		arr.assign( _Index({ 1, 0}) , ValueType(0.5));
		arr.assign( _Index({ 1, 1}) , ValueType(-3.0));
		arr.assign( _Index({ 1, 2}) , ValueType(5.0));

		arr.assign( _Index({ 2, 0}) , ValueType(8.0));
		arr.assign( _Index({ 2, 1}) , ValueType(5.0));
		arr.assign( _Index({ 2, 2}) , ValueType(3.2));
	}

	Poly2LoaderExample< ValueType, 0 , 1 , 1> subloader;
};

template <typename ValueType, bool MakeNonLinear>
struct Poly2LoaderExample< ValueType, 2 , 1 , 1, MakeNonLinear>
{
	void load_array( typename MD<ValueType, 3 >::RegArray& arr)
	{
		if (MakeNonLinear)
		{
			arr.assign( _Index({ 0, 0, 1}) , ValueType(0.1));

			arr.assign( _Index({ 1, 1, 2}) , ValueType(0.1));

			arr.assign( _Index({ 2, 2, 1}) , ValueType(0.1));
		}
	}

	Poly2LoaderExample< ValueType, 1 , 1 , 1> subloader;
};


#endif
