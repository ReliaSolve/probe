// Enable functions with up to 20 parameters to be called.  Default of 15 is insufficient
#define BOOST_PYTHON_MAX_ARITY 20
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "DotSpheres.h"
#include "SpatialQuery.h"
#include "Scoring.h"

using namespace boost::python;
using namespace molprobity::probe;

/// @brief Helper function to wrap the Point internal array so we can read its elements.
///
/// If you want to change values in a Point, construct a new one using the 3-parameter
/// constructor.
boost::python::tuple wrap_vec3_array(Point const& d) {
  boost::python::list a;
  for (size_t i = 0; i < d.size(); ++i) {
    a.append(d.elems[i]);
  }
  return boost::python::tuple(a);
}

BOOST_PYTHON_MODULE(mmtbx_probe_ext)
{
  // Describe and name compound classes that we need access to.

  // This class is already exposed as its more basic vec3<double> class, but if we
  // don't include this we get an error saying that no class is registered...
  class_<Point>("Point", init<double, double, double>())
    .def(init<>())
    .def("size", &Point::size)
    .add_property("elems", wrap_vec3_array)
  ;

  class_<ContactResult>("ContactResult", init<>())
    .add_property("closestContact", &ContactResult::closestContact)
    .add_property("distAboveSurface", &ContactResult::distAboveSurface)
    ;

  class_<ExtraAtomInfo>("ExtraAtomInfo", init<>())
    .add_property("vdwRadius", &ExtraAtomInfo::getVdwRadius, &ExtraAtomInfo::setVdwRadius)
    .add_property("isAcceptor", &ExtraAtomInfo::getIsAcceptor, &ExtraAtomInfo::setIsAcceptor)
    .add_property("isDonor", &ExtraAtomInfo::getIsDonor, &ExtraAtomInfo::setIsDonor)
    .add_property("isDummyHydrogen", &ExtraAtomInfo::getIsDummyHydrogen, &ExtraAtomInfo::setIsDummyHydrogen)
    ;

  // These are read-only methods.  They cannot be filled in by Python.
  class_<AtomVsAtomDotScorer::ScoreDotsResult>("ScoreDotsResult", init<>())
    .add_property("valid", &AtomVsAtomDotScorer::ScoreDotsResult::valid)
    .add_property("bumpSubScore", &AtomVsAtomDotScorer::ScoreDotsResult::bumpSubScore)
    .add_property("hBondSubScore", &AtomVsAtomDotScorer::ScoreDotsResult::hBondSubScore)
    .add_property("attractSubScore", &AtomVsAtomDotScorer::ScoreDotsResult::attractSubScore)
    .add_property("hasBadBump", &AtomVsAtomDotScorer::ScoreDotsResult::hasBadBump)
    .def("totalScore", &AtomVsAtomDotScorer::ScoreDotsResult::totalScore)
    ;

  // Export the vector indexing of objects that we'll use vectors for.
  // NOTE: Everything that is using scitbx::af::shared "flex" arrays is
  // automatically wrapped for us in ways that let them be used as standard
  // Python iterators so we don't need to add the wrapping.
  //typedef scitbx::af::shared<Point> PointList;
  //class_<PointList>("PointList");

  // Export the classes we define
  class_<DotSphere>("DotSphere", init<double, double>())
    .def(init<>())
    .def("dots", &DotSphere::dots, return_internal_reference<>())
    .def("radius", &DotSphere::radius)
    .def("density", &DotSphere::density)
    .def("test", &DotSphere::test)
    .def("XXX", &DotSphere::XXX)
  ;

  class_<DotSphereCache>("DotSphereCache", init<double>())
    .def("get_sphere", &DotSphereCache::get_sphere, return_internal_reference<>())
    .def("size", &DotSphereCache::size)
    .def("test", &DotSphereCache::test)
  ;

  class_<SpatialQuery>("SpatialQuery", init<Point, Point, Point>())
    .def(init<scitbx::af::shared<iotbx::pdb::hierarchy::atom> const>())
    .def("add", &SpatialQuery::add)
    .def("remove", &SpatialQuery::remove)
    .def("neighbors", &SpatialQuery::neighbors)
    .def("test", &SpatialQuery::test)
  ;

  class_<AtomVsAtomDotScorer>("AtomVsAtomDotScorer", init<scitbx::af::shared<ExtraAtomInfo>,
        optional<double, double, double, double, double, double> >())
    .def("score_dots", &AtomVsAtomDotScorer::score_dots)
    .def("test", &AtomVsAtomDotScorer::test)
    ;

  // Export the global functions
  def("closest_contact", closest_contact, "Point of closest contact and distance for dot on atom.");
  def("atom_charge", atom_charge, "Integer charge on an atom given its string charge description.");
  def("dot2srcCenter", dot2srcCenter, "Distance from dot to point on the source surface closest to target.");
  def("kissEdge2bullsEye", kissEdge2bullsEye, ".");
  def("annularDots", annularDots, ".");

  def("DotSpheres_test", DotSpheres_test, "Test all classes defined in DotSpheres.h.");
  def("SpatialQuery_test", SpatialQuery_test, "Test all classes defined in SpatialQuery.h.");
  def("Scoring_test", Scoring_test, "Test all classes defined in Scoring.h.");

}

