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

int atom_charge(iotbx::pdb::hierarchy::atom const& atom)
{
  // Get the tidy version of the charge string stripped of any extra characters.
  // Detangle the fact that the string claims to have the potential to be optional.
  // This will make the string either blank or consisting of two characters, [12][+-],
  // with the first indicating the magnitude and the second indicating the sign.
  int ret = 0;
  std::string chStr;
  auto opt = atom.charge_tidy(true);
  if (opt) { chStr = opt.get(); }

  if (chStr.size() > 0) {
    ret = chStr[0] - '0';
    if (chStr[1] == '-') { ret *= -1; }
  }

  return ret;
}

AtomVsAtomDotScorer::ScoreDotsResult AtomVsAtomDotScorer::score_dots(
  iotbx::pdb::hierarchy::atom sourceAtom, double minOccupancy,
  SpatialQuery spatialQuery, double nearbyRadius, double probeRadius,
  std::vector<iotbx::pdb::hierarchy::atom> exclude,
  std::vector<Point> dots, double density, bool onlyBumps)
{
  // This method is based on AtomPositions::atomScore() from Reduce.
  // It is passed only the dots that it should score rather than excluding them
  // internally like that function does.

  // Default return has 0 value for all subscores.
  ScoreDotsResult ret;

  // Check for invalid parameters
  if (density <= 0) {
    return ret;
  }

  // Make sure that the occupancy of the atom is high enough to score it.
  if (std::abs(sourceAtom.data->occ) < minOccupancy) {
    return ret;
  }
  /// @todo Check the conformation to make sure this atom is visible in it

  // Find the neighboring atoms that are potentially interacting.
  // The nonzero minimum distance prevents us from selecting the source atom.
  std::vector<iotbx::pdb::hierarchy::atom> neighbors = spatialQuery.neighbors(sourceAtom.data->xyz, 0.001, nearbyRadius);

  // Select only those atoms actually interacting: that have sufficient occupancy, for whom the
  // gap between the Van Der Waals surfaces is less than the probe radius, and which
  // are not in the excluded list.
  /// @todo Check the conformation to make sure they are both visible to each other.
  std::vector<iotbx::pdb::hierarchy::atom> interacting;
  unsigned sourceID = sourceAtom.data->i_seq;
  ExtraAtomInfo const& sourceExtra = m_extraInfo[sourceID];
  for (iotbx::pdb::hierarchy::atom const &a : neighbors) {
    unsigned aID = a.data->i_seq;
    if ((sourceID >= m_extraInfo.size()) || (aID >= m_extraInfo.size())) {
      return ret;
    }
    ExtraAtomInfo const& aExtra = m_extraInfo[aID];
    double nonBondedDistance = sourceExtra.vdwRadius + aExtra.vdwRadius;
    if ((std::abs(a.data->occ) >= minOccupancy)
          && ((a.data->xyz - sourceAtom.data->xyz).length() <= nonBondedDistance + probeRadius)
      @todo replace this find with a scan looking for common data pointer get() outputs due to no == being defined
          && (std::find(exclude.begin(), exclude.end(), a) == exclude.end())
        ) {
      interacting.push_back(a);
    }
  }

  /// @todo Useful tidbits.  Remove when done.
  bool h = sourceAtom.element_is_hydrogen();
  /// @todo end useful tidbits

  // Run through all of the dots and determine whether and how to score each.
  for (Point const& d : dots) {
    // Find the world-space location of the dot by adding it to the location of the source atom.
    // The probe location is in the same direction as d from the source but is further away by the
    // probe radius.
    Point q = sourceAtom.data->xyz + d;
    Point probLoc = sourceAtom.data->xyz + d.normalize() * (d.length() + probeRadius);

    double minGap = 1e10;                 ///< Nearest atom that we found
    bool isHydrogenBond = false;          ///< Are we looking at a hydrogen bond to our neighbor?
    bool tooCloseHydrogenBond = false;    ///< Are we too close to be a hydrogen bond?
    double hydgrogenBondMinDist = 1e10;   ///< Nearest hydrogen bond we found
    bool keepDot = false;                 ///< Did we find a neighbor and we're not in a bonded atom?

    iotbx::pdb::hierarchy::atom const *cause = nullptr;

    // Look through each atom to find the one with the smallest gap, which is the one that
    // the dot would interact with.
    for (iotbx::pdb::hierarchy::atom const& b : interacting) {
      ExtraAtomInfo const& bExtra = m_extraInfo[b.data->i_seq];
      Point locb = b.data->xyz;
      double vdwb = bExtra.vdwRadius;

      // See if we are too far away to interact, bail if so.
      double squareDist = (probLoc - locb).length_sq();
      double pRadPlusVdwb = vdwb + probeRadius;
      if (squareDist > pRadPlusVdwb * pRadPlusVdwb) {
        continue;
      }

      // See if we are within the probe radius past the edge.  If so, we're in contention
      // to be the nearest atom.
      double dist = sqrt(squareDist);
      double probeGap = dist - pRadPlusVdwb;
      double gap = dist - vdwb;
      if (probeGap < 0) {

        // See if we replace the currently-closest atom.
        if (gap < minGap) {

          // Figure out what kind of interaction this is based on the atom types and
          // charge status of the two atoms.
          int chargeSource = atom_charge(sourceAtom);
          int chargeB = atom_charge(b);

          bool bothCharged = (chargeSource != 0) && (chargeB != 0);
          bool chargeComplement = bothCharged && (chargeSource * chargeB < 0);

          // See if one of the atoms is a hydrogen donor and the other can accept hydrogen bonds.
          bool couldHBond = (sourceExtra.isDonor && bExtra.isAcceptor)
            || (sourceExtra.isAcceptor && bExtra.isDonor);
          if (couldHBond && ((!bothCharged) || chargeComplement)) {
            isHydrogenBond = true;
            hydgrogenBondMinDist = bothCharged ? m_minChargedHydrogenBondGap : m_minRegularHydrogenBondGap;
            tooCloseHydrogenBond = (gap < -hydgrogenBondMinDist);
          } else {
            // If this is a dummy hydrogen, then we skip it, it can only be a hydrogen-bond partner
            if (bExtra.isDummyHydrogen) { continue; }
            // This is not a hydrogen bond.
            isHydrogenBond = tooCloseHydrogenBond = false;
          }

          // Record which atom is the closest and mark the dot to be kept because we found an
          // atom that is close enough.
          cause = &b;
          keepDot = true;
          minGap = gap;
        }
      }
    }

    // If the dot was close enough to some non-bonded atom, check to see if it should be removed
    // from consideration because it is also inside an excluded atom.
    if (keepDot) {
      /// @todo Check to make sure they are in interacting alternate conformations.
      for (iotbx::pdb::hierarchy::atom const& e : exclude) {
        double vdwe = m_extraInfo[e.data->i_seq].vdwRadius;
        if ((q - e.data->xyz).length_sq() < vdwe * vdwe) {
          keepDot = false;
          break;
        }
      }
    }
    if (!keepDot) { continue; }

    // We've gotten this far, so we score the dot and add it to the relevant subscore.

    /// @todo


    /// @todo Make sure that we add to the attraction subscore in the overlapType == 0 case
  }

  ret.bumpSubScore /= density;
  ret.hBondSubScore /= density;
  ret.attractSubScore /= density;
  ret.valid = true;
  return ret;
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
