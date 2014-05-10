
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

#include <gtest/gtest.h>

#include <sparse_arr.h>
#include <utility>
#include <iostream>
#include <tensor_product_map.h>
#include "symbolicc++.h"



class SparseMDArrayTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.


  SparseMDArrayTest() {
    // You can do set-up work for each test here.
  }

  virtual ~SparseMDArrayTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects declared here can be used by all tests in the test case for Foo.
};

TEST_F(SparseMDArrayTest, MatrixVectorProduct) {
	MD<Symbolic, 2>::RegArray M_arr(0.0);
	MD<Symbolic, 1>::RegArray s_arr(0.0);

	Symbolic s0 = Symbolic("s_0");
	Symbolic s1 = Symbolic("s_1");

	s_arr.assign( _Index({0}), s0);
	s_arr.assign( _Index({1}), s1);

	Symbolic M00 = Symbolic("M_00");
	Symbolic M01 = Symbolic("M_01");
	Symbolic M10 = Symbolic("M_10");
	Symbolic M11 = Symbolic("M_11");
	M_arr.assign( _Index({0,0}), M00);
	M_arr.assign( _Index({0,1}), M01);
	M_arr.assign( _Index({1,0}), M10);
	M_arr.assign( _Index({1,1}), M11);

	MD<Symbolic, 1>::RegArray t_arr(0.0);
	ContractedIndexValues_t c_idx;
	c_idx[1] = 0;
	TensorProductMap<Symbolic, 2, Symbolic, 1, Symbolic, 1> MatrixVectorProdSym( M_arr, c_idx, s_arr, t_arr);

	EXPECT_TRUE( t_arr( _Index({0}) ) == M00*s0 + M01*s1 );
	EXPECT_TRUE( t_arr( _Index({1}) ) == M10*s0 + M11*s1 );
}

TEST_F(SparseMDArrayTest, DotProduct) {
	MD<Symbolic, 1>::RegArray u_arr(0.0);
	MD<Symbolic, 1>::RegArray v_arr(0.0);

	Symbolic u0 = Symbolic("u_0");
	Symbolic u1 = Symbolic("u_1");
	
	Symbolic v0 = Symbolic("v_0");
	Symbolic v1 = Symbolic("v_1");

	u_arr.assign( _Index({0}), u0);
	u_arr.assign( _Index({1}), u1);

	v_arr.assign( _Index({0}), v0);
	v_arr.assign( _Index({1}), v1);



	Symbolic dotProductResult;
	ContractedIndexValues_t c_idx;
	c_idx[0] = 0;
	TensorProductMap<Symbolic, 1, Symbolic, 1, Symbolic, 0> dotProductSym( u_arr, c_idx, v_arr, dotProductResult);

	EXPECT_TRUE( dotProductResult == u0*v0 + u1*v1 );
}


TEST_F(SparseMDArrayTest, MatrixMatrixProduct) {
	MD<Symbolic, 2>::RegArray M_arr(0.0);
	MD<Symbolic, 2>::RegArray N_arr(0.0);
	MD<Symbolic, 2>::RegArray P_arr(0.0);


	Symbolic M00 = Symbolic("M_00");
	Symbolic M01 = Symbolic("M_01");
	Symbolic M10 = Symbolic("M_10");
	Symbolic M11 = Symbolic("M_11");

	M_arr.assign( _Index({0,0}), M00);
	M_arr.assign( _Index({0,1}), M01);
	M_arr.assign( _Index({1,0}), M10);
	M_arr.assign( _Index({1,1}), M11);

	Symbolic N00 = Symbolic("N_00");
	Symbolic N01 = Symbolic("N_01");
	Symbolic N10 = Symbolic("N_10");
	Symbolic N11 = Symbolic("N_11");

	N_arr.assign( _Index({0,0}), N00);
	N_arr.assign( _Index({0,1}), N01);
	N_arr.assign( _Index({1,0}), N10);
	N_arr.assign( _Index({1,1}), N11);

	ContractedIndexValues_t c_idx;
	c_idx[1] = 0;
	TensorProductMap<Symbolic, 2, Symbolic, 2, Symbolic, 2> MatrixMatrixProdSym( M_arr, c_idx, N_arr, P_arr);

	EXPECT_TRUE( P_arr( _Index({0,0}) ) == M00*N00 + M01*N10 );
	EXPECT_TRUE( P_arr( _Index({0,1}) ) == M00*N01 + M01*N11 );
	EXPECT_TRUE( P_arr( _Index({1,0}) ) == M10*N00 + M11*N10 );
	EXPECT_TRUE( P_arr( _Index({1,1}) ) == M10*N01 + M11*N11 );
}


TEST_F(SparseMDArrayTest, TensorMatrixContraction) {
	MD<Symbolic, 3>::RegArray M_arr(0.0);
	MD<Symbolic, 2>::RegArray N_arr(0.0);
	MD<Symbolic, 1>::RegArray P_arr(0.0);


	Symbolic M000 = Symbolic("M_000");
	Symbolic M001 = Symbolic("M_001");
	Symbolic M010 = Symbolic("M_010");
	Symbolic M011 = Symbolic("M_011");
	Symbolic M100 = Symbolic("M_100");
	Symbolic M101 = Symbolic("M_101");
	Symbolic M110 = Symbolic("M_110");
	Symbolic M111 = Symbolic("M_111");

	M_arr.assign( _Index({0,0,0}), M000);
	M_arr.assign( _Index({0,0,1}), M001);
	M_arr.assign( _Index({0,1,0}), M010);
	M_arr.assign( _Index({0,1,1}), M011);
	M_arr.assign( _Index({1,0,0}), M100);
	M_arr.assign( _Index({1,0,1}), M101);
	M_arr.assign( _Index({1,1,0}), M110);
	M_arr.assign( _Index({1,1,1}), M111);

	Symbolic N00 = Symbolic("N_00");
	Symbolic N01 = Symbolic("N_01");
	Symbolic N10 = Symbolic("N_10");
	Symbolic N11 = Symbolic("N_11");

	N_arr.assign( _Index({0,0}), N00);
	N_arr.assign( _Index({0,1}), N01);
	N_arr.assign( _Index({1,0}), N10);
	N_arr.assign( _Index({1,1}), N11);

	// the contraction is M_ijk N_kj, the convention is that the ContractedIndexValues_t array has the position of the contracted index in the first tensor as the key, and the position of the contracted
	// index in the second tensor as the value. So in this case, the contracted index j is at position 1 in tensor M, and at position 1 in tensor N. Which means that we need to add an entry 1 => 1 in the map.
	// Likewise, the repeated index k is at position 2 on M tensor, and at position 0 on N tensor, so we need to add the entry 2 => 0 in the map.
	ContractedIndexValues_t c_idx;
	c_idx[1] = 1;
	c_idx[2] = 0;
	TensorProductMap<Symbolic, 3, Symbolic, 2, Symbolic, 1> MatrixMatrixProdSym1( M_arr, c_idx, N_arr, P_arr);

	EXPECT_TRUE( P_arr( _Index({0}) ) == M000*N00 + M001*N10 + M010*N01 + M011*N11 );
	EXPECT_TRUE( P_arr( _Index({1}) ) == M100*N00 + M101*N10 + M110*N01 + M111*N11 );

	// now we'll calculate M_ijk N_jk
	c_idx[1] = 0;
	c_idx[2] = 1;
	TensorProductMap<Symbolic, 3, Symbolic, 2, Symbolic, 1> MatrixMatrixProdSym2( M_arr, c_idx, N_arr, P_arr);

	EXPECT_TRUE( P_arr( _Index({0}) ) == M000*N00 + M001*N01 + M010*N10 + M011*N11 );
	EXPECT_TRUE( P_arr( _Index({1}) ) == M100*N00 + M101*N01 + M110*N10 + M111*N11 );
}


int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
};
