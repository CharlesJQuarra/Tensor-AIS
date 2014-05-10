
#ifndef _LOADER_POLY1_H
#define _LOADER_POLY1_H

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

template <typename ValueType, unsigned int Order, unsigned int RankX, unsigned int RankY, bool SingleVariable>
struct Poly1LoaderExample
{
	void load_array( typename MD<ValueType, RankY + Order*RankX >::RegArray& arr)
	{
	}

	Poly1LoaderExample< ValueType, Order-1 , RankX , RankY , SingleVariable > subloader;
};



//single variable case: (x-1)(x-3)(x+3) = x**2 - 4x + 3 (x+3) = x**3 - 4 x**2 + 3x + 3 x**2 - 12 x + 9 = x**3 - x**2 - 9 x + 9



//multiple variable case: a root in (x0r, x1r, x2r) = (-3.0 , 5.0 , 1.5)
// test multinomial:  (x0 - x0r)(x1 - x1r)(x2 - x2r)(x0 + - 2 x1 + x2 + 5.0) =
/*
	( x0 x1 - x0r x1 - x1r x0 + x0r x1r ) ( x2 - x2r ) = x0 x1 x2 - x0r x1 x2 - x1r x0 x2 + x0r x1r x2 - x2r x0 x1 + x0r x2r x1 + x1r x2r x0 - x0r x1r x2r
  ... * (x0 + -2 x1 + x2 + 5.0) = 
		x0^2 x1 x2 - x0r x0 x1 x2 - x1r x0^2 x2 + x0r x1r x0 x2 - x2r x0^2 x1 + x0r x2r x0 x1 + x1r x2r x0^2 - x0r x1r x2r x0
	- 2 x0 x1^2 x2 + 2 x0r x1^2 x2 + 2 x1r x0 x1 x2 - 2 x0r x1r x1 x2 + 2 x2r x0 x1^2 - 2 x0r x2r x1^2 - 2 x1r x2r x0 x1 + 2 x0r x1r x2r x1
	+ x0 x1 x2^2 - x0r x1 x2^2 - x1r x0 x2^2 + x0r x1r x2^2 - x2r x0 x1 x2 + x0r x2r x1 x2 + x1r x2r x0 x2 - x0r x1r x2r x2
	+ 5 x0 x1 x2 - 5 x0r x1 x2 - 5 x1r x0 x2 + 5 x0r x1r x2 - 5 x2r x0 x1 + 5 x0r x2r x1 + 5 x1r x2r x0 - 5 x0r x1r x2r

	4-order:
		x0^2 x1 x2 - 2 x0 x1^2 x2 + x0 x1 x2^2
	3-order:
	(2 x1r -x0r - x2r + 5) x0 x1 x2 + (-x1r) x0^2 x2 + (-x2r) x0^2 x1 + (2 x0r) x1^2 x2 + (2 x2r) x0 x1^2 + (-x0r) x1 x2^2 + (-x1r) x0 x2^2 
	2-order:
	( x0r + x2r - 5 ) x1r  x0 x2 + (x0r - 2 x1r - 5) x2r x0 x1 + (x1r x2r) x0^2 + (x2r - 2 x1r - 5) x0r x1 x2 + (-2 x0r x2r) x1^2 + (x0r x1r) x2^2
	1-order:
	( 5 - x0r ) x1r x2r x0 + (2 x1r + 5) x0r x2r x1 + (5 - x2r) x0r x1r x2
	0-order:
	-5 x0r x1r x2r



static constexpr float x0r = -3.0;
static constexpr float x1r = 5.0;
static constexpr float x2r = 1.5;*/

//0.609692,3.91482,1.5
/*
static constexpr float x0r = 0.609692;
static constexpr float x1r = 3.91482;
static constexpr float x2r = 1.5;
*/
//4.79883,6.78763,3.77642
static constexpr double x0r = 4.79883;
static constexpr double x1r = 6.78763;
static constexpr double x2r = 3.77642;

double poly_test_eval(double x0, double x1, double x2)
{
	return (x0- x0r)*(x1- x1r)*(x2-x2r)*(x0 - 2.0*x1 + x2 + 5.0);
};

double poly_test_eval_expanded(double x0, double x1, double x2)
{
	return pow( x0 , 2 ) * x1 * x2 - 2.0 * x0 * pow( x1 , 2 ) * x2 + x0 * x1 * pow( x2 , 2 )
		+ ( 2.0* x1r - x0r - x2r + 5.0)* x0*x1*x2 - x1r * pow( x0, 2)* x2 - x2r * pow( x0 , 2)*x1 + 2.0*x0r* pow(x1, 2)*x2 + 2.0*x2r*x0*pow(x1,2) - x0r*x1*pow(x2, 2) - x1r*x0*pow(x2, 2)
		+( x0r + x2r - 5.0)*x1r*x0*x2 + (x0r - 2.0*x1r - 5.0)*x2r*x0*x1 + x1r*x2r*pow(x0,2) + (x2r - 2.0*x1r- 5.0)*x0r*x1*x2 - 2.0*x0r*x2r*pow(x1, 2) + x0r*x1r*pow(x2, 2)
		+ (5.0 - x0r)*x1r*x2r*x0 + (2.0*x1r + 5.0)*x0r*x2r*x1 + (5.0-x2r)*x0r*x1r*x2
		- 5.0 *x0r*x1r*x2r;
};


template <typename ValueType, bool SingleVariable>
struct Poly1LoaderExample< ValueType, 0 , 1 , 1 , SingleVariable>
{
	void load_array( typename MD<ValueType, 1 >::RegArray& arr)
	{

		if (!SingleVariable)
		{
			arr.assign( _Index({ 0}) , ValueType(-5.0 * x0r * x1r* x2r));
		} else
			arr.assign( _Index({ 0}) , ValueType(3.0));
	}

};

template <typename ValueType, bool SingleVariable>
struct Poly1LoaderExample< ValueType, 1 , 1 , 1 , SingleVariable>
{
	void load_array( typename MD<ValueType, 2 >::RegArray& arr)
	{
		/*
			( 5 - x0r ) x1r x2r x0 + (2 x1r + 5) x0r x2r x1 + (5 - x2r) x0r x1r x2
		*/
		if (!SingleVariable)
		{
			arr.assign( _Index({ 0, 0}) , ValueType(x1r * x2r * (5.0 - x0r)));
			arr.assign( _Index({ 0, 1}) , ValueType(x0r * x2r * (2.0 * x1r + 5.0)));
			arr.assign( _Index({ 0, 2}) , ValueType(x0r * x1r * (5.0 - x2r)));
		}
		else
			arr.assign( _Index({ 0, 0}) , ValueType(-9.0));
	}

	Poly1LoaderExample< ValueType, 0 , 1 , 1 , SingleVariable> subloader;
};

template <typename ValueType, bool SingleVariable>
struct Poly1LoaderExample< ValueType, 2 , 1 , 1 , SingleVariable>
{
	void load_array( typename MD<ValueType, 3 >::RegArray& arr)
	{
		/*
	( x0r + x2r - 5 ) x1r  x0 x2 + (x0r - 2 x1r - 5) x2r x0 x1 + (x1r x2r) x0^2 + (x2r - 2 x1r - 5) x0r x1 x2 + (-2 x0r x2r) x1^2 + (x0r x1r) x2^2
	2-order:(0, 0, 0) = x1r x2r
					(0, 0, 1) = x2r ( x0r - 2 x1r  - 5)
					(0, 0, 2) = x1r ( x0r + x2r - 5)
					(0, 1, 1) = - 2 x0r x2r
					(0, 1, 2) = x0r (x2r - 2 x1r - 5)
					(0, 2, 2) = x0r x1r
		*/
		if (!SingleVariable)
		{
			arr.assign( _Index({ 0, 0, 0}) , ValueType(x1r * x2r));
			arr.assign( _Index({ 0, 0, 1}) , ValueType(x2r * ( x0r - 2.0 * x1r - 5.0)));
			arr.assign( _Index({ 0, 0, 2}) , ValueType(x1r * ( x0r + x2r - 5.0)));



			arr.assign( _Index({ 0, 1, 1}) , ValueType(-2.0 * x0r * x2r));
			arr.assign( _Index({ 0, 1, 2}) , ValueType(x0r * (x2r - 2.0*x1r - 5.0)));
			arr.assign( _Index({ 0, 2, 2}) , ValueType(x0r * x1r));
		}
		else
			arr.assign( _Index({ 0, 0, 0}) , ValueType(-1.0));
	}

	Poly1LoaderExample< ValueType, 1 , 1 , 1 , SingleVariable> subloader;
};

template <typename ValueType, bool SingleVariable>
struct Poly1LoaderExample< ValueType, 3 , 1 , 1 , SingleVariable>
{
	void load_array( typename MD<ValueType, 4 >::RegArray& arr)
	{
			/*
	(2 x1r -x0r - x2r + 5) x0 x1 x2 + (-x1r) x0^2 x2 + (-x2r) x0^2 x1 + (2 x0r) x1^2 x2 + (2 x2r) x0 x1^2 + (-x0r) x1 x2^2 + (-x1r) x0 x2^2
	3-order:(0, 0, 0, 1) = - x2r
					(0, 0, 0, 2) = - x1r
					(0, 0, 1, 2) = 2 x1r - x0r - x2r + 5
					(0, 0, 1, 1) = 2 x2r
					(0, 0, 2, 2) = - x1r
					(0, 1, 1, 2) = 2 x0r
					(0, 1, 2, 2) = - x0r
			*/
			if (!SingleVariable)
			{

				arr.assign( _Index({ 0, 0, 0, 1 }) , ValueType(-x2r) );
				arr.assign( _Index({ 0, 0, 0, 2 }) , ValueType(-x1r) );
				arr.assign( _Index({ 0, 0, 1, 2 }) , ValueType(2.0*x1r - x0r - x2r + 5.0));
				arr.assign( _Index({ 0, 0, 1, 1 }) , ValueType(2.0*x2r));
				arr.assign( _Index({ 0, 0, 2, 2 }) , ValueType(-x1r));
				arr.assign( _Index({ 0, 1, 1, 2 }) , ValueType(2.0*x0r));
				arr.assign( _Index({ 0, 1, 2, 2 }) , ValueType(-x0r));
			}
			else
				arr.assign( _Index({ 0, 0, 0, 0 }) , ValueType(1.0));
	}

	Poly1LoaderExample< ValueType, 2 , 1 , 1 , SingleVariable> subloader;
};



template <typename ValueType, bool SingleVariable>
struct Poly1LoaderExample< ValueType, 4 , 1 , 1 , SingleVariable>
{
	void load_array( typename MD<ValueType, 5 >::RegArray& arr)
	{
			/*
4-order:(0, 0, 0, 1, 2) = 1
					(0, 0, 1, 1, 2) = -2
					(0, 0, 1, 2, 2) = 1
*/
			if (!SingleVariable)
			{
				arr.assign( _Index({ 0, 0, 0, 1, 2 }) , ValueType(1.0) );
				arr.assign( _Index({ 0, 0, 1, 1, 2 }) , ValueType(-2.0) );
				arr.assign( _Index({ 0, 0, 1, 2, 2 }) , ValueType(1.0) );
			}
	}

	Poly1LoaderExample< ValueType, 3 , 1 , 1 , SingleVariable> subloader;
};


#endif
