#include <boost/python.hpp>

namespace p = boost::python;

char const* hello_world()
{
    return "Hello World\n";
}

BOOST_PYTHON_MODULE(GapsHelloWorld)
{
    Py_Initialize();
    p::def("hello_world", hello_world, "");
}