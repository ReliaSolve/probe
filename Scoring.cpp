// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#include <cmath>
#include <scitbx/constants.h>
#include <algorithm>
#include "Scoring.h"
#include "DotSpheres.h"

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
  SpatialQuery &spatialQuery, double nearbyRadius, double probeRadius,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &exclude,
  scitbx::af::shared<Point> const &dots, double density, bool onlyBumps)
{
  // This method is based on AtomPositions::atomScore() from Reduce.
  // It is passed only the dots that it should score rather than excluding them
  // internally like that function does.

  // If we're looking for other contacts besides bumps, we need to add twice the
  // probe radius to the neighbor list to ensure that we include atoms where the
  // probe spans from the source to the target.
  if (!onlyBumps) {
    nearbyRadius += 2 * probeRadius;
  }

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

  // Find the neighboring atoms that are potentially interacting.
  // The nonzero minimum distance prevents us from selecting the source atom.
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> neighbors = 
    spatialQuery.neighbors(sourceAtom.data->xyz, 0.001, nearbyRadius);

  // Select only those atoms actually interacting: that have sufficient occupancy, for whom the
  // gap between the Van Der Waals surfaces is less than the probe radius, and which
  // are not in the excluded list.
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

  // Run through all of the dots and determine whether and how to score each.
  for (Point const& d : dots) {
    // Find the world-space location of the dot by adding it to the location of the source atom.
    // The probe location is in the same direction as d from the source but is further away by the
    // probe radius.
    Point absoluteDotLocation = sourceAtom.data->xyz + d;
    Point probLoc = sourceAtom.data->xyz + d.normalize() * (d.length() + probeRadius);

    double minGap = 1e10;                 ///< Nearest atom that we found
    bool isHydrogenBond = false;          ///< Are we looking at a hydrogen bond to our neighbor?
    bool tooCloseHydrogenBond = false;    ///< Are we too close to be a hydrogen bond?
    double hydrogenBondMinDist = 0;       ///< Hydrogen bond minimum distance based on the atom types (will be set below).
    bool keepDot = false;                 ///< Did we find a neighbor and we're not in a bonded atom?

    iotbx::pdb::hierarchy::atom const *cause = nullptr;

    // Look through each atom to find the one with the smallest gap, which is the one that
    // the dot would interact with.
    for (iotbx::pdb::hierarchy::atom const& b : interacting) {
      ExtraAtomInfo const& bExtra = m_extraInfo[b.data->i_seq];
      Point locb = b.data->xyz;
      double vdwb = bExtra.getVdwRadius();

      // See if we are too far away to interact, bail if so.
      double squareProbeDist = (probLoc - locb).length_sq();
      double pRadPlusVdwb = vdwb + probeRadius;
      if (squareProbeDist > pRadPlusVdwb * pRadPlusVdwb) {
        continue;
      }

      // At this point, we are within the probe radius past the edge, so we're in contention
      // to be the nearest atom.  Find the distance from the dot rather than from the probe
      // to check our actual interaction behavior.
      double dist = (absoluteDotLocation - locb).length();
      double gap = dist - vdwb;

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
          hydrogenBondMinDist = bothCharged ? m_minChargedHydrogenBondGap : m_minRegularHydrogenBondGap;
          tooCloseHydrogenBond = (gap < -hydrogenBondMinDist);
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

    // If the dot was close enough to some non-bonded atom, check to see if it should be removed
    // from consideration because it is also inside an excluded atom.
    if (keepDot) {
      for (iotbx::pdb::hierarchy::atom const& e : exclude) {
        double vdwe = m_extraInfo[e.data->i_seq].getVdwRadius();
        if ((absoluteDotLocation - e.data->xyz).length_sq() < vdwe * vdwe) {
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
    if (minGap >= 0) {
      overlap = 0;
      overlapType = 0;
    } else if (isHydrogenBond) {
      overlap = -0.5 * minGap;
      if (tooCloseHydrogenBond) {
        minGap += hydrogenBondMinDist;
        overlapType = -1;
      } else {
        overlapType = 1;
      }
    } else {  // minGap < 0 and not a hydrogen bond
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
      if ((!onlyBumps) && !annularDots(absoluteDotLocation, sourceAtom.data->xyz, sourceExtra.getVdwRadius(),
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
      break;
      
    default:
      // This should never happen.  Returns with ret invalid to indicate an error.
      std::cerr << "AtomVsAtomDotScorer::score_dots(): Internal error: Unrecognized overlap type: " << overlapType << std::endl;
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

std::string AtomVsAtomDotScorer::test()
{
  // Construct test cases with all combinations of charges and extra information, holding the
  // radii of the neighbor atom and probe atom constant.  Do this in combination with adding or
  // not adding an excluded atom that completely covers the neighbor atom.
  // Run tests against all of these cases to ensure that the bahavior is as expected in each case.
  // This tests the case logic within the code.
  // The target radius has to be large enough to get a bad bump even for hydrogen bond cases.
  double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;

  // Construct the dot sphere to be used.
  DotSphere ds(sourceRad, 200);

  for (std::string targetCharge : {"--", "-", "", "+", "++"}) {
  for (std::string sourceCharge : {"--", "-", "", "+", "++"}) {
   for (bool targetAccept : { false, true}) {
    for (bool sourceAccept : { false, true}) {
     for (bool targetDonor : { false, true}) {
      for (bool sourceDonor : { false, true}) {
       for (bool targetDummy : { false, true}) {
         for (bool onlyBumps : { false, true}) {
           for (bool excludeAtom : { false, true}) {

             //================================================================
             // Test the scoring for various cases to ensure that they all behave as expected
             unsigned int atomSeq = 0;

             // Construct and fill the SpatialQuery information
             // with a vector of a single target atom, including its extra info looked up by
             // its i_seq value.
             iotbx::pdb::hierarchy::atom a;
             a.set_charge(targetCharge.c_str());
             a.set_xyz({ 0,0,0 });
             a.set_occ(1);
             a.data->i_seq = atomSeq++;
             scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
             atoms.push_back(a);
             SpatialQuery sq(atoms);
             ExtraAtomInfo e(targetRad, targetAccept, targetDonor, targetDummy);
             scitbx::af::shared<ExtraAtomInfo> infos;
             infos.push_back(e);

             // Construct the source atom, including its extra info looked up by
             // its i_seq value.
             iotbx::pdb::hierarchy::atom source;
             source.set_charge(sourceCharge.c_str());
             source.set_occ(1);
             source.data->i_seq = atomSeq++;
             ExtraAtomInfo se(sourceRad, sourceAccept, sourceDonor, false);
             infos.push_back(se);

             // Construct the scorer to be used.
             AtomVsAtomDotScorer as(infos);

             // Determine our hydrogen-bond state
             bool compatibleCharge = atom_charge(source) * atom_charge(a) <= 0;
             bool compatible = (sourceDonor && targetAccept) || (sourceAccept && targetDonor);
             bool hBond = compatibleCharge && compatible;

             // If we have an excluded atom, we should always get no values or bumping.
             // Skip the remainder of the tests in this case
             scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;
             if (excludeAtom) {
                // Describe the extra atom to the system, including its extra info looked up by
                // its i_seq value.
                iotbx::pdb::hierarchy::atom ea;
                ea.set_xyz({ 0,0,0 });
                ea.set_occ(1);
                ea.data->i_seq = atomSeq++;
                ExtraAtomInfo ex(targetRad + 0.2, targetAccept, targetDonor, targetDummy);
                infos.push_back(ex);
                exclude.push_back(ea);

                // Even when we have a close clash, we should get no response.
                source.set_xyz({ sourceRad,0,0 });
                ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                  probeRad, exclude, ds.dots(), ds.density(), onlyBumps);
                if (!res.valid) {
                  return "AtomVsAtomDotScorer::test(): Could not score dots for excluded-atom case";
                }
                if ((res.totalScore() != 0) || res.hasBadBump) {
                  return "AtomVsAtomDotScorer::test(): Got unexpected result for excluded-atom case";
                }

                // Skip the rest of the tests for this case.
                continue;
             }

             // If we have a dummy hydrogen and we cannot be a hydrogen-bond pair,
             // we should always get no bumping.
             if (targetDummy) {
               if (!hBond) {
                 // Even when we have a close clash, we should get no response.
                 source.set_xyz({ sourceRad,0,0 });
                 ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                   probeRad, exclude, ds.dots(), ds.density(), onlyBumps);
                 if (!res.valid) {
                   return "AtomVsAtomDotScorer::test(): Could not score dots for dummy hydrogen case";
                 }
                 if ((res.bumpSubScore != 0) || res.hasBadBump) {
                   return "AtomVsAtomDotScorer::test(): Got unexpected result for dummy hydrogen case";
                 }
                 // Skip the rest of the tests for this case.
                 continue;
               }
             }

             // When we get so close that the source atom radius touches the center of the target,
             // we should get bad bumps in all cases.
             {
               source.set_xyz({ sourceRad,0,0 });
               ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                 probeRad, exclude, ds.dots(), ds.density(), onlyBumps);
               if (!res.valid) {
                 return "AtomVsAtomDotScorer::test(): Could not score dots for bad-bump case";
               }
               if (!res.hasBadBump) {
                 return "AtomVsAtomDotScorer::test(): Got no bad bump for bad-bump case case";
               }
             }

             // When we are only checking for bumps, we should get no interaction when the
             // atoms are not touching.  Otherwise, slight interaction.
             {
               source.set_xyz({ sourceRad + targetRad + 0.001,0,0 });
               ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                 probeRad, exclude, ds.dots(), ds.density(), onlyBumps);
               if (!res.valid) {
                 return "AtomVsAtomDotScorer::test(): Could not score dots for bump-only test case";
               }
               if (onlyBumps) {
                 if (res.totalScore() != 0) {
                   return "AtomVsAtomDotScorer::test(): Got value when not expected for bump-only test case";
                 }
               } else {
                 if (res.totalScore() == 0) {
                   return "AtomVsAtomDotScorer::test(): Got no value when one expected for non-bump-only test case";
                 }
               }
             }

             // When we are only checking for bumps, even hydrogen bonds should be counted as bumps.
             if (onlyBumps) {
               source.set_xyz({ sourceRad + targetRad - 0.1,0,0 });
               ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                 probeRad, exclude, ds.dots(), ds.density(), onlyBumps);
               if (!res.valid) {
                 return "AtomVsAtomDotScorer::test(): Could not score dots for bump-only hydrogen-bond test case";
               }
               if (res.hBondSubScore != 0) {
                 return "AtomVsAtomDotScorer::test(): Got unexpected hydrogen bond score for bump-only test case";
               }
               if (res.bumpSubScore >= 0) {
                 return "AtomVsAtomDotScorer::test(): Got unexpected bump score for bump-only test case";
               }
             }

           }
         }
       }
      }
     }
    }
   }
  }
  }

  // Sweep an atom from just touching to far away and make sure the attract
  // curve is monotonically decreasing to 0.
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a single target atom, including its extra info looked up by
    // its i_seq value.
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz({ 0,0,0 });
    a.set_occ(1);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct the source atom, including its extra info looked up by
    // its i_seq value.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    infos.push_back(se);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct the scorer to be used.
    AtomVsAtomDotScorer as(infos);

    // Sweep the source atom
    double lastAttract = 1e10;
    bool foundNonzero = false;
    for (double gap = 0; gap < 10; gap += 0.1) {
      source.set_xyz({ targetRad + sourceRad + gap, 0, 0 });
      ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density());
      if (!res.valid) {
        return "AtomVsAtomDotScorer::test(): Could not score dots for swept-distance case";
      }
      if ((res.attractSubScore != res.totalScore()) || (res.attractSubScore > lastAttract)) {
        return "AtomVsAtomDotScorer::test(): Non-monotonic scores for swept-distance case";
      }
      lastAttract = res.attractSubScore;
      if (lastAttract != 0) { foundNonzero = true; }
    }
    if (!foundNonzero) {
      return "AtomVsAtomDotScorer::test(): No nonzero scores for swept-distance case";
    }
    if (lastAttract != 0) {
      return "AtomVsAtomDotScorer::test(): Non-empty last score for swept-distance case";
    }
  }

  // Test the setting of weights for the various subscores.
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a two target atoms, including extra info looked up by
    // their i_seq values.  One is an acceptor and the other is not.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    scitbx::af::shared<ExtraAtomInfo> infos;

    { // The first is a bond acceptor at the origin
      iotbx::pdb::hierarchy::atom a;
      a.set_xyz({ 0,0,0 });
      a.set_occ(1);
      a.data->i_seq = atomSeq++;
      atoms.push_back(a);
      ExtraAtomInfo e(targetRad, true);
      infos.push_back(e);
    }

    { // The second is not a bond acceptor and is a bit less than the source atom away from the edge of the first.
      iotbx::pdb::hierarchy::atom a;
      a.set_xyz({ 2 * targetRad + 2*(sourceRad*0.8),0,0 });
      a.set_occ(1);
      a.data->i_seq = atomSeq++;
      atoms.push_back(a);
      ExtraAtomInfo e(targetRad, false);
      infos.push_back(e);
    }
    SpatialQuery sq(atoms);

    // Construct the source atom, including its extra info looked up by
    // its i_seq value.  It is a hydrogen donor and is located halfway
    // between the two atoms.
    iotbx::pdb::hierarchy::atom source;
    source.set_xyz({ targetRad + sourceRad * 0.8 ,0,0 });
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad, false, true);
    infos.push_back(se);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Results holds a vector of six values, three for the gap, bump, and bond weights used in
    // each run and three for the attract, bump, and hBond scores measured with these weights.
    std::vector< std::vector<double> > results;

    // Run tests with a hydrogen source molecule that can H-bond with one of two target atoms so that
    // we get results in all three of the score results when we compare it.
    for (double wGap = 0.25; wGap < 10; wGap += 3) {
      for (double wBump = 0.25; wBump < 10; wBump += 3) {
        for (double wBond = 0.25; wBond < 10; wBond += 3) {

          // Construct the scorer to be used.
          AtomVsAtomDotScorer as(infos,wGap, wBump, wBond);

          // Get the dot results and store all of the params and results in a vector
          ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
            probeRad, exclude, ds.dots(), ds.density());
          if (!res.valid) {
            return "AtomVsAtomDotScorer::test(): Could not score dots for weight scaling cases";
          }
          std::vector<double> row = { wGap, wBump, wBond, res.attractSubScore, res.bumpSubScore, res.hBondSubScore };
          results.push_back(row);
        }
      }
    }

    // Ensure that the ratio of the scores matches the ratio of the weights for each entry between
    // the first row and all other rows.
    for (size_t i = 1; i < results.size(); i++) {
      // The attraction score is not a linear weighting of the
      // dot scores, so we can only check those for equality or inequality.
      bool paramClose = closeTo(1, results[0][0] / results[i][0]);
      bool resultClose = closeTo(1, results[0][3 + 0] / results[i][3 + 0]);
      if (paramClose != resultClose) {
        return "AtomVsAtomDotScorer::test(): Inconsistent gap ratio for weight scaling case";
      }

      // Check the other two scores for linearity.
      for (size_t j = 1; j < 3; j++) {
        double paramRatio = results[0][j] / results[i][j];
        double resultRatio = results[0][3 + j] / results[i][3 + j];
        if (!closeTo(paramRatio, resultRatio)) {
          return "AtomVsAtomDotScorer::test(): Unequal ratio for weight scaling case";
        }
      }
    }

  }

  // Test the setting of bond-gap distances.
  {
    double badBondGap = 0.2;
    double minRegularHydrogenBondGap = badBondGap + 0.2;
    double minChargedHydrogenBondGap = minRegularHydrogenBondGap + 0.2;

    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a single target atom, including its extra info looked up by
    // its i_seq value.  It will be an acceptor and it will have a negative charge
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz({ 0,0,0 });
    a.set_occ(1);
    a.set_charge("-");
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad, true);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct a source atom, including its extra info looked up by
    // its i_seq value.  This will be a hydrogen but not a donor to check
    // for the standard bad-bump result.
    // Test it against both sides of the bad-bump line to see if it responds correctly.
    {
      iotbx::pdb::hierarchy::atom source;
      source.set_occ(1);
      source.data->i_seq = atomSeq++;
      ExtraAtomInfo se(sourceRad);
      infos.push_back(se);
      ScoreDotsResult res;

      // Construct the scorer to be used with the specified bond gaps.
      AtomVsAtomDotScorer as(infos, 0.25, 10.0, 4.0, minRegularHydrogenBondGap, minChargedHydrogenBondGap, badBondGap);

      // Check the source atom against outside and inside the gap
      source.set_xyz({ targetRad + sourceRad - badBondGap + 0.1, 0, 0 });
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density());
      if (!res.valid) {
        return "AtomVsAtomDotScorer::test(): Could not score dots for badBondGap setting case";
      }
      if (res.hasBadBump) {
        return "AtomVsAtomDotScorer::test(): Bad bump found when not expected for badBondGap setting case";
      }

      source.set_xyz({ targetRad + sourceRad - badBondGap - 0.1, 0, 0 });
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density());
      if (!res.valid) {
        return "AtomVsAtomDotScorer::test(): Could not score dots for badBondGap setting case";
      }
      if (!res.hasBadBump) {
        return "AtomVsAtomDotScorer::test(): Bad bump not found when expected for badBondGap setting case";
      }
    }

    // Construct a source atom, including its extra info looked up by
    // its i_seq value.  This will be an uncharged hydrogen donor to check
    // for the non-charged hydrogen bad-bump result.
    // Test it against both sides of the bad-bump line to see if it responds correctly.
    {
      iotbx::pdb::hierarchy::atom source;
      source.set_occ(1);
      source.data->i_seq = atomSeq++;
      ExtraAtomInfo se(sourceRad,false, true);
      infos.push_back(se);
      ScoreDotsResult res;

      // Construct the scorer to be used with the specified bond gaps.
      AtomVsAtomDotScorer as(infos, 0.25, 10.0, 4.0, minRegularHydrogenBondGap, minChargedHydrogenBondGap, badBondGap);

      // Check the source atom against outside and inside the gap
      source.set_xyz({ targetRad + sourceRad - minRegularHydrogenBondGap - badBondGap + 0.1, 0, 0 });
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density());
      if (!res.valid) {
        return "AtomVsAtomDotScorer::test(): Could not score dots for minRegularHydrogenBondGap setting case";
      }
      if (res.hasBadBump) {
        return "AtomVsAtomDotScorer::test(): Bad bump found when not expected for minRegularHydrogenBondGap setting case";
      }

      source.set_xyz({ targetRad + sourceRad - minRegularHydrogenBondGap - badBondGap - 0.1, 0, 0 });
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density());
      if (!res.valid) {
        return "AtomVsAtomDotScorer::test(): Could not score dots for minRegularHydrogenBondGap setting case";
      }
      if (!res.hasBadBump) {
        return "AtomVsAtomDotScorer::test(): Bad bump not found when expected for minRegularHydrogenBondGap setting case";
      }
    }

    // Construct a source atom, including its extra info looked up by
    // its i_seq value.  This will be a charged hydrogen donor to check
    // for the charged hydrogen bad-bump result.
    // Test it against both sides of the bad-bump line to see if it responds correctly.
    {
      iotbx::pdb::hierarchy::atom source;
      source.set_occ(1);
      source.set_charge("+");
      source.data->i_seq = atomSeq++;
      ExtraAtomInfo se(sourceRad, false, true);
      infos.push_back(se);
      ScoreDotsResult res;

      // Construct the scorer to be used with the specified bond gaps.
      AtomVsAtomDotScorer as(infos, 0.25, 10.0, 4.0, minRegularHydrogenBondGap, minChargedHydrogenBondGap, badBondGap);

      // Check the source atom against outside and inside the gap
      source.set_xyz({ targetRad + sourceRad - minChargedHydrogenBondGap - badBondGap + 0.1, 0, 0 });
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density());
      if (!res.valid) {
        return "AtomVsAtomDotScorer::test(): Could not score dots for minChargedHydrogenBondGap setting case";
      }
      if (res.hasBadBump) {
        return "AtomVsAtomDotScorer::test(): Bad bump found when not expected for minChargedHydrogenBondGap setting case";
      }

      source.set_xyz({ targetRad + sourceRad - minChargedHydrogenBondGap - badBondGap - 0.1, 0, 0 });
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density());
      if (!res.valid) {
        return "AtomVsAtomDotScorer::test(): Could not score dots for minChargedHydrogenBondGap setting case";
      }
      if (!res.hasBadBump) {
        return "AtomVsAtomDotScorer::test(): Bad bump not found when expected for minChargedHydrogenBondGap setting case";
      }
    }

  }

  // Test the control of occupancy level
  /// @todo

  /// @todo
  return "";
}

std::string Scoring_test()
{
  std::string ret;
  ContactResult res;

  // Test the atom-charge code.
  std::vector<std::string> charges = { "--", "-", "", "+", "++", "+2", "-1", "0" };
  std::vector<int> expectedCharge = { -2, -1, 0, 1, 2, 2, -1, 0 };
  for (size_t i = 0; i < charges.size(); i++) {
    iotbx::pdb::hierarchy::atom a;
    a.set_charge(charges[i].c_str());
    if (atom_charge(a) != expectedCharge[i]) {
      return "Scoring_test: atom_charge() failed";
    }
  }


  // Test the AtomVsAtomDotScorer class
  ret = AtomVsAtomDotScorer::test();
  if (ret.size() > 0) {
    return std::string("Scoring_test: AtomVsAtomDotScorer failed: ") + ret;
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
