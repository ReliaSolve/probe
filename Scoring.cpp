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
    ret.distAboveSurface = -atom_radius;
    ret.closestContact = atom + Point(atom_radius, 0, 0);

  // Otherwise, find the point of closest approach and its distance and return
  // them.
  } else {
    double len = diff.length();
    ret.distAboveSurface = len - atom_radius;
    ret.closestContact = atom + diff * atom_radius / len;
  }
  return ret;
}

double dot2srcCenter(const Point& dot, const Point& srcLoc, double srcVDWRad, const Point& targLoc) {
  // The vector from the source pointing towards the target that is the radius of the source atom
  Point src2targVec = (targLoc - srcLoc).normalize() * srcVDWRad;
  // The point on the surface of the source atom that is closest to the target
  Point srcSurfacePoint = src2targVec + srcLoc;
  // The distance from the dot to the point on the source surface that is closest to the target
  return (srcSurfacePoint - dot).length();
}

double kissEdge2bullsEye(double ra, double rb, double rp) {
  return 2 * ra * sqrt(rb * rp / ((ra + rb) * (ra + rp)));
}
bool annularDots(const Point& dot, const Point& srcLoc, double srcVDWRad,
  const Point& targLoc, double targVDWRad, double probeRadius) {
  return dot2srcCenter(dot, srcLoc, srcVDWRad, targLoc) > kissEdge2bullsEye(srcVDWRad, targVDWRad, probeRadius);
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
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude,
  scitbx::af::shared<Point> dots, double density, bool onlyBumps)
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
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> neighbors = 
    spatialQuery.neighbors(sourceAtom.data->xyz, 0.001, nearbyRadius);

  // Select only those atoms actually interacting: that have sufficient occupancy, for whom the
  // gap between the Van Der Waals surfaces is less than the probe radius, and which
  // are not in the excluded list.
  /// @todo Check the conformation to make sure they are both visible to each other.
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> interacting;
  unsigned sourceID = sourceAtom.data->i_seq;
  ExtraAtomInfo const& sourceExtra = m_extraInfo[sourceID];
  for (iotbx::pdb::hierarchy::atom const &a : neighbors) {
    unsigned aID = a.data->i_seq;
    if ((sourceID >= m_extraInfo.size()) || (aID >= m_extraInfo.size())) {
      return ret;
    }
    ExtraAtomInfo const& aExtra = m_extraInfo[aID];
    double nonBondedDistance = sourceExtra.getVdwRadius() + aExtra.getVdwRadius();
    bool excluded = false;
    for (iotbx::pdb::hierarchy::atom const& e : exclude) {
      if (e.data.get() == a.data.get()) {
        excluded = true;
        break;
      }
    }
    if ((std::abs(a.data->occ) >= minOccupancy)
          && ((a.data->xyz - sourceAtom.data->xyz).length() <= nonBondedDistance + probeRadius)
          && (!excluded)
        ) {
      interacting.push_back(a);
    }
  }

  /// @todo Useful tidbits.  Remove when done.
  // bool h = sourceAtom.element_is_hydrogen();
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
      double vdwb = bExtra.getVdwRadius();

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
          bool couldHBond = (sourceExtra.getIsDonor() && bExtra.getIsAcceptor())
            || (sourceExtra.getIsAcceptor() && bExtra.getIsDonor());
          if (couldHBond && ((!bothCharged) || chargeComplement)) {
            isHydrogenBond = true;
            hydgrogenBondMinDist = bothCharged ? m_minChargedHydrogenBondGap : m_minRegularHydrogenBondGap;
            tooCloseHydrogenBond = (gap < -hydgrogenBondMinDist);
          } else {
            // If this is a dummy hydrogen, then we skip it, it can only be a hydrogen-bond partner
            if (bExtra.getIsDummyHydrogen()) { continue; }
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
        double vdwe = m_extraInfo[e.data->i_seq].getVdwRadius();
        if ((q - e.data->xyz).length_sq() < vdwe * vdwe) {
          keepDot = false;
          break;
        }
      }
    }
    if (!keepDot) { continue; }

    // We've gotten this far, so we score the dot and add it to the relevant subscore.
    int overlapType = 0;
    double overlap = 0;

    // Determine the overlap type and amount of overlap.
    if (minGap > 0) {
      overlap = 0;
      overlapType = 0;
    } else if (isHydrogenBond && tooCloseHydrogenBond) {
      minGap += hydgrogenBondMinDist;
      overlap = -0.5 * minGap;
      overlapType = -1;
    } else if (isHydrogenBond) {
      overlap = -0.5 * minGap;
      overlapType = 1;
    } else if (minGap < 0) {
      overlap = -0.5 * minGap;
      overlapType = -1;
    }

    // Compute the score for the dot based on the overlap type and amount of overlap.
    // Assign it to the appropriate subscore.
    double dotScore = 0;
    switch (overlapType) {

    case -1:  // Clash
      ret.bumpSubScore += -m_bumpWeight * overlap;
      // See if we should flag this atom as having a bad bump
      if (minGap < -m_badBumpBondGap) {
        ret.hasBadBump = true;
      }
      break;

    case 0:   // Contact dot
      if ((!onlyBumps) && !annularDots(q, sourceAtom.data->xyz, sourceExtra.getVdwRadius(),
          cause->data->xyz, m_extraInfo[cause->data->i_seq].getVdwRadius(), probeRadius)) {

        double scaledGap = minGap / m_gapWeight;
        ret.attractSubScore += exp(-scaledGap * scaledGap);
      }
      break;

    case 1:   // Hydrogen bond
      if (!onlyBumps) {
        ret.hBondSubScore += m_hBondWeight * overlap;
      } else {  // In this case, we treat it as a bump
        ret.bumpSubScore += -m_bumpWeight * overlap;
      }
      
    default:
      // This should never happen.  Returns with ret invalid to indicate an error.
      return ret;
    }
  }

  // Normalize the score by the density so that it does not depend on the number of dots
  // that were constructed for the atom.
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

/// @brief holds parameters needed to initialize an atom and its associated extra information.
class AtomInfo {
public:
  Point loc;
  double occ;
  std::string charge;
  double radius;
  bool isAcceptor;
  bool isDonor;
  bool isDummyHydrogen;
};
static std::vector<AtomInfo> testAtoms = {
  { {0,0,0}, 1, "", 1.0, false, true, false}      // Uncharged hydrogen at the origin
};

std::string AtomVsAtomDotScorer::test()
{
  // Construct a model with some atoms in it so that we can use it for our tests.
  // Construct atoms with locations and extra info pulled from the data table above
  iotbx::pdb::hierarchy::model m;
  iotbx::pdb::hierarchy::chain c;
  iotbx::pdb::hierarchy::residue_group rg;
  iotbx::pdb::hierarchy::atom_group ag;
  size_t numAtoms = 10;
  Coord spacing = 5;
  for (int x = 0; x < numAtoms; x++) {
    for (int y = 0; y < numAtoms; y++) {
      for (int z = 0; z < numAtoms; z++) {
        Point v(x * spacing, y * spacing, z * spacing);
        iotbx::pdb::hierarchy::atom a(v, v);
        ag.append_atom(a);
      }
    }
  }
  rg.append_atom_group(ag);
  c.append_residue_group(rg);
  m.append_chain(c);

  // Sweep various atom types from far away to near and make sure their interaction
  // curves match what is expected.

  /// @todo
  return "";
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
  if (!closeTo(res.distAboveSurface, -1)) {
    return "Scoring_test: closest_contact() returned bad distance for same point";
  }
  if (!closeTo((ones - res.closestContact).length(), 1)) {
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
          if (!closeTo(rad, (res.closestContact - ones).length())) {
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
    if (!closeTo(x - rad, res.distAboveSurface)) {
      return "Scoring_test: closest_contact() returned bad distance for point";
    }
  }

  // All tests passed.
  return "";
}


} // end namespace probe
} // end namespace molprobity 
