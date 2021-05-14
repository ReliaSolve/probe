// Enable functions with up to 20 parameters to be called.  Default of 15 is insufficient
#define BOOST_PYTHON_MAX_ARITY 20
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <scitbx/boost_python/container_conversions.h>
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
  // Describe and name compound classes that we need access to beyond those that are
  // already defined for us by scitbx arrays that are defined elsewhere.

  class_<ContactResult>("ContactResult", init<>())
    .add_property("closestContact", &ContactResult::closestContact)
    .add_property("distAboveSurface", &ContactResult::distAboveSurface)
    ;

  class_<ExtraAtomInfo>("ExtraAtomInfo", init< optional<double, bool, bool, bool> >())
    .add_property("vdwRadius", &ExtraAtomInfo::getVdwRadius, &ExtraAtomInfo::setVdwRadius)
    .add_property("isAcceptor", &ExtraAtomInfo::getIsAcceptor, &ExtraAtomInfo::setIsAcceptor)
    .add_property("isDonor", &ExtraAtomInfo::getIsDonor, &ExtraAtomInfo::setIsDonor)
    .add_property("isDummyHydrogen", &ExtraAtomInfo::getIsDummyHydrogen, &ExtraAtomInfo::setIsDummyHydrogen)
    ;
  // Define the flex array wrapping for this class because we take it as a parameter.
  scitbx::boost_python::container_conversions::tuple_mapping_variable_capacity<
    scitbx::af::shared<ExtraAtomInfo> >();

  enum_<DotScorer::InteractionType>("InteractionType")
    .value("None", DotScorer::InteractionType::None)
    .value("Clash", DotScorer::InteractionType::Clash)
    .value("NearContact", DotScorer::InteractionType::NearContact)
    .value("HydrogenBond", DotScorer::InteractionType::HydrogenBond)
    ;

  class_<DotScorer::CheckDotResult>("CheckDotResult", init<>())
    .add_property("overlapType", &DotScorer::CheckDotResult::overlapType)
    .add_property("cause", &DotScorer::CheckDotResult::cause)
    .add_property("overlap", &DotScorer::CheckDotResult::overlap)
    .add_property("gap", &DotScorer::CheckDotResult::gap)
    .add_property("annular", &DotScorer::CheckDotResult::annular)
    ;

  class_<DotScorer::ScoreDotsResult>("ScoreDotsResult", init<>())
    .add_property("valid", &DotScorer::ScoreDotsResult::valid)
    .add_property("bumpSubScore", &DotScorer::ScoreDotsResult::bumpSubScore)
    .add_property("hBondSubScore", &DotScorer::ScoreDotsResult::hBondSubScore)
    .add_property("attractSubScore", &DotScorer::ScoreDotsResult::attractSubScore)
    .add_property("hasBadBump", &DotScorer::ScoreDotsResult::hasBadBump)
    .def("totalScore", &DotScorer::ScoreDotsResult::totalScore)
    ;

  class_<DotSphere>("DotSphere", init<double, double>())
    .def(init<>())
    .def("dots", &DotSphere::dotsCopyForPythonWrapping)
    .def("radius", &DotSphere::radius)
    .def("density", &DotSphere::density)
    .def("test", &DotSphere::test)
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

  class_<DotScorer>("DotScorer", init<scitbx::af::shared<ExtraAtomInfo>,
        optional<double, double, double, double, double, double> >())
    .def("check_dot", &DotScorer::check_dot)
    .def("score_dots", &DotScorer::score_dots)
    .def("test", &DotScorer::test)
    ;

  // Export the vector indexing of objects that we'll use vectors for.
  // NOTE: Everything that is using scitbx::af::shared "flex" arrays is
  // automatically wrapped for us in ways that let them be used as standard
  // Python iterators so we don't need to add the wrapping.  We only need
  // to describe any std::vector values that we use.

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

