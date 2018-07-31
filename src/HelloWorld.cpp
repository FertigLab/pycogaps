#include <boost/python.hpp>
#include "Rpackage/src/GapsRunner.h"
#include "Rpackage/src/math/Random.h"

namespace p = boost::python;

char const* hello_world()
{
    GapsRunner runner("data/GIST.tsv", false, 3, false,
        std::vector<unsigned>(1));

    gaps::random::setSeed(123);
    runner.recordSeed(123);
    runner.setMaxIterations(1000);
    runner.setSparsity(0.01f, 0.01f, false);
    runner.setMaxGibbsMass(100.f, 100.f);
    runner.setMaxThreads(1);
    runner.setPrintMessages(true);
    runner.setOutputFrequency(250);
    runner.setCheckpointOutFile("gaps_checkpoint.out");
    runner.setCheckpointInterval(0);

    runner.run();

    return "Hello World\n";
}

BOOST_PYTHON_MODULE(CogapsPy)
{
    Py_Initialize();
    p::def("hello_world", hello_world, "");
}