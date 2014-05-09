import os

env=Environment()
os.system('echo "check pkg-config: PKG_CONFIG_PATH=$PKG_CONFIG_PATH"')
env["ENV"]["PKG_CONFIG_PATH"] = os.environ.get("PKG_CONFIG_PATH")
env.ParseConfig('pkg-config --cflags --libs boost-1.52')
env.ParseConfig('pkg-config --cflags --libs symbolic')
env.ParseConfig('pkg-config --cflags --libs gtest')
env.MergeFlags('-lpthread')

solver_source = "solver/"
tests = "tests/"

cpppath = ['common/']
tests_paths = ['common/', solver_source, tests]

tests_paths.extend(env['CPPPATH'])


ccflags = ['-std=c++11']
ccflags.extend(env['CCFLAGS'])



testDir = Dir(tests)

sparseArrTest = testDir.glob("test_sparse_arr.cpp", True, True, False)
bigDataPolyTest = testDir.glob("test_bigdatapoly.cpp", True, True, False)
bigDataPolyDiffTest = testDir.glob("test_bigdatapoly_diff.cpp", True, True, False)
iterativeSolverTest = testDir.glob("test_alg_iterative_solver.cpp", True, True, False)


env.Program('sparse_arr', sparseArrTest, CPPPATH=tests_paths, CCFLAGS=ccflags)
env.Program('bigdatapoly', bigDataPolyTest, CPPPATH=tests_paths, CCFLAGS=ccflags)
env.Program('bigdatapoly_diff', bigDataPolyDiffTest, CPPPATH=tests_paths, CCFLAGS=ccflags)
env.Program('iterative_solver', iterativeSolverTest, CPPPATH=tests_paths, CCFLAGS=ccflags)




