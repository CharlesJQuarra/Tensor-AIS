

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

#include <solver.h>
#include <poly.h>
#include <iostream>
#include <loader_poly1.h>
#include <loader_poly2.h>


typedef long double real_t;
//typedef DelayedFloat<double> real_t;

class AlgebraicIterativeSolverTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

	MD< real_t , 1 >::RegArray y_goal;
	MD< real_t , 1 >::RegArray x_test;
	real_t phi;

  AlgebraicIterativeSolverTest() : phi(9e30), y_goal(), x_test() {
    // You can do set-up work for each test here.
  }

  virtual ~AlgebraicIterativeSolverTest() {
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

TEST_F(AlgebraicIterativeSolverTest, LinearSolution) {
	x_test.assign( _Index({0}) , 9.1 );
	x_test.assign( _Index({1}) , -3.1 );
	x_test.assign( _Index({2}) , -10.2 );

  Poly2LoaderExample< real_t, 2 , 1 , 1 , false> linearSymmetricLoader;
	AlgebraicIterativeSolver< real_t, 2 , 1 , 1 > ais_2_linear(linearSymmetricLoader, y_goal);
	ais_2_linear.setCurrentXTest(x_test);
	int iterations = 500;
	do
	{
		phi = ais_2_linear.line_search_update();
	} while (phi > 1e-10 && (iterations-- > 0));
	std::cout << "phi = " << phi << std::endl;
	std::cout << " test data solution: x = ["<< x_test( _Index({0}) ) << "," << x_test( _Index({1}) ) << "," << x_test( _Index({2}) ) << "]  " << std::endl;
	EXPECT_TRUE(phi < 1e-10);
	// x = [-2.12614,2.53306,0.732442]
	EXPECT_TRUE( fabs( x_test( _Index({0}) ) + 2.12614 ) < 1e-5 );
	EXPECT_TRUE( fabs( x_test( _Index({1}) ) - 2.53306 ) < 1e-5 );
	EXPECT_TRUE( fabs( x_test( _Index({2}) ) - 0.732442 ) < 1e-5 );
}

TEST_F(AlgebraicIterativeSolverTest, NonlinearSolution) {
	x_test.assign( _Index({0}) , 9.1 );
	x_test.assign( _Index({1}) , -3.1 );
	x_test.assign( _Index({2}) , -10.2 );

  Poly2LoaderExample< real_t, 2 , 1 , 1 , true> nonlinearSymmetricLoader;
	AlgebraicIterativeSolver< real_t, 2 , 1 , 1 > ais_2_nonlinear(nonlinearSymmetricLoader, y_goal);
	ais_2_nonlinear.setCurrentXTest(x_test);
	int iterations = 500;
	do
	{
		phi = ais_2_nonlinear.line_search_update();
	} while (phi > 1e-10 && (iterations-- > 0));
	std::cout << "phi = " << phi << std::endl;
	std::cout << " test data solution: x = ["<< x_test( _Index({0}) ) << "," << x_test( _Index({1}) ) << "," << x_test( _Index({2}) ) << "]  " << std::endl;
	EXPECT_TRUE(phi < 1e-10);
	// x = [-2.31105,2.72661,0.822234]
	EXPECT_TRUE( fabs( x_test( _Index({0}) ) + 2.31105 ) < 1e-5 );
	EXPECT_TRUE( fabs( x_test( _Index({1}) ) - 2.72661 ) < 1e-5 );
	EXPECT_TRUE( fabs( x_test( _Index({2}) ) - 0.822234 ) < 1e-5 );
}


int main(int argc, char **argv)
{

	::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
};


