// Enable functions with up to 20 parameters to be called.  Default of 15 is insufficient
#define BOOST_PYTHON_MAX_ARITY 20
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "DotSpheres.h"

using namespace boost::python;
using namespace molprobity::probe;

/// @brief Helper function to wrap the Point internal array so we can read its elements.
///
/// If you want to change values in a Point, construct a new one using the 3-parameter
/// constructor.
boost::python::tuple wrap_vec3_array(Point const& d) {
  boost::python::list a;
  for (int i = 0; i < d.size(); ++i) {
    a.append(d.elems[i]);
  }
  return boost::python::tuple(a);
}

BOOST_PYTHON_MODULE(mmtbx_probe_ext)
{
  // Describe and name compound classes that we need access to.
  class_<Point>("Point", init<double, double, double>())
    .def(init<>())
    .def("size", &Point::size)
    .add_property("elems", wrap_vec3_array)
  ;

  // Describe vectors we need access to
  typedef std::vector<Point> PointList;
  class_<PointList>("PointList")
    .def(vector_indexing_suite<PointList>())
  ;

  // Export the classes we define
  class_<DotSphere>("DotSphere", init<double, double>())
    .def(init<>())
    .def("dots", &DotSphere::dots, return_internal_reference<>())
    .def("radius", &DotSphere::radius)
    .def("density", &DotSphere::density)
    .def("test", &DotSphere::test)
  ;

  class_<DotSphereCache>("DotSphereCache", init<double>())
    .def("get_sphere", &DotSphereCache::get_sphere, return_internal_reference<>())
    .def("size", &DotSphereCache::size)
    .def("test", &DotSphereCache::test)
  ;

  // Export the vector indexing of objects that we'll use vectors for.

  // Export the global functions
  def("DotSpheres_test", DotSpheres_test, "Test all classes defined in DotSpheres.h.");
}

