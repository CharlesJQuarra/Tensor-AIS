
#ifndef _LOADER_POLY_SYM_H_
#define _LOADER_POLY_SYM_H_

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

#include "symbolicc++.h"

template <typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY>
struct PolySymLoaderExample
{
	void load_array( typename MD<ValueType, RankY + Order*RankX >::RegArray& arr)
	{
	}

	PolySymLoaderExample< ValueType, Order-1 , RankX , RankY> subloader;
};


template <>
struct PolySymLoaderExample< Symbolic, 0 , 1 , 1>
{
	void load_array( typename MD<Symbolic, 1 >::RegArray& arr)
	{
		arr.assign( _Index({ 0}) , Symbolic("c_0_0"));
		arr.assign( _Index({ 1}) , Symbolic("c_0_1"));
	}

};

template <>
struct PolySymLoaderExample< Symbolic, 1 , 1 , 1>
{
	void load_array( typename MD<Symbolic, 2 >::RegArray& arr)
	{
		arr.assign( _Index({ 0, 0}) , Symbolic("c_1_00"));
		arr.assign( _Index({ 0, 1}) , Symbolic("c_1_01"));
		arr.assign( _Index({ 1, 0}) , Symbolic("c_1_10"));
		arr.assign( _Index({ 1, 1}) , Symbolic("c_1_11"));
	}

	PolySymLoaderExample< Symbolic, 0 , 1 , 1> subloader;
};

template <>
struct PolySymLoaderExample< Symbolic, 2 , 1 , 1>
{
	void load_array( typename MD<Symbolic, 3 >::RegArray& arr)
	{
		arr.assign( _Index({ 0, 0, 0}) , Symbolic("c_2_000"));
		arr.assign( _Index({ 0, 0, 1}) , Symbolic("c_2_001"));
		arr.assign( _Index({ 0, 1, 0}) , Symbolic("c_2_010"));
		arr.assign( _Index({ 0, 1, 1}) , Symbolic("c_2_011"));

		arr.assign( _Index({ 1, 0, 0}) , Symbolic("c_2_100"));
		arr.assign( _Index({ 1, 0, 1}) , Symbolic("c_2_101"));
		arr.assign( _Index({ 1, 1, 0}) , Symbolic("c_2_110"));
		arr.assign( _Index({ 1, 1, 1}) , Symbolic("c_2_111"));
	}

	PolySymLoaderExample< Symbolic, 1 , 1 , 1> subloader;
};

#endif
