#include <boost/python.hpp>
//#include <boost/python/numpy.hpp>
//#include <numpy/numpyconfig.h>
//#include <numpy/ndarrayobject.h>
//#include <iostream>

namespace p = boost::python;
namespace np = boost::python::numpy;

p::tuple give_back_matrices(np::ndarray A, np::ndarray B)
{
  return p::make_tuple(A, B);
}

void hello_world(void)
{
  std::cout << "Hello World Boost Python!\n";
}

BOOST_PYTHON_MODULE(boosttest)
{
  p::def("give_back_matrices", give_back_matrices, "");
  p::def("hello_world", hello_world, "");
}
