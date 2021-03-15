// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#include <cmath>
#include <scitbx/constants.h>
#include <algorithm>
#include "Scoring.h"

namespace molprobity {
  namespace probe {


ContactResult closest_contact(Point dot, Point atom, double atom_radius)
{
  ContactResult ret;
  Point diff = dot - atom;

  // if the dot is at the center of the atom, then pick an arbitrary point
  // on the surface and return it.
  if (diff.length_sq() == 0) {
    ret.dist_above_surface = -atom_radius;
    ret.closest_contact = atom + Point(atom_radius, 0, 0);

  // Otherwise, find the point of closest approach and its distance and return
  // them.
  } else {
    double len = diff.length();
    ret.dist_above_surface = len - atom_radius;
    ret.closest_contact = atom + diff * atom_radius / len;
  }
  return ret;
}

AtomVsAtomDotScorer::ScoreDotsResult AtomVsAtomDotScorer::score_dots(
  iotbx::pdb::hierarchy::atom sourceAtom, iotbx::pdb::hierarchy::atom targetAtom,
  std::vector<Point> dots, double density, bool onlyBumps)
{

}

//===========================================================================================================
// Testing code below here

/// @brief Returns true of the two floating-point numbers are nearly equal
static bool closeTo(double a, double b) {
  return fabs(a - b) < 1e-10;
}

std::string Scoring_test()
{
  std::string ret;
  ContactResult res;

  // Test the AtomVsAtomDotScorer class
  ret = AtomVsAtomDotScorer::test();
  if (ret.size() > 0) {
    return std::string("Scoring_test: failed: ") + ret;
  }

  // Test that the distance from a location to itself is negative radius, and that
  // the projection is the radius away from the origin.
  Point ones(1, 1, 1);
  res = closest_contact(ones, ones, 1);
  if (!closeTo(res.dist_above_surface, -1)) {
    return "Scoring_test: closest_contact() returned bad distance for same point";
  }
  if (!closeTo((ones - res.closest_contact).length(), 1)) {
    return "Scoring_test: closest_contact() returned bad closest contact for same point";
  }

  // Test that the projection distance is the same for all points on the cardinal directions
  // and diagonal directions around a sphere.
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      for (int z = -1; z <= 1; z++) {
        if (x != 0 || y != 0 || z != 0) {
          Point test(x, y, z);
          test += ones;
          double rad = 0.5;
          res = closest_contact(test, ones, rad);
          if (!closeTo(rad, (res.closest_contact - ones).length())) {
            std::cout << "XXX at " << res.closest_contact[0] << ", " << res.closest_contact[1] << ", " << res.closest_contact[2] << std::endl;
            return "Scoring_test: closest_contact() returned bad closest contact for surrounding point";
          }
        }
      }
    }
  }

  // Test that the distance scales as expected as we move away from the center of an atom
  for (int x = 0; x < 100; x++) {
    Point test(x, 0, 0);
    test += ones;
    double rad = 25;
    res = closest_contact(test, ones, rad);
    if (!closeTo(x - rad, res.dist_above_surface)) {
      return "Scoring_test: closest_contact() returned bad distance for point";
    }
  }

  // All tests passed.
  return "";
}


} // end namespace probe
} // end namespace molprobity 
