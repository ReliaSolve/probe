// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#pragma once

#include "Common.h"
#include <cctbx/geometry_restraints/nonbonded.h>
#include <iotbx/pdb/hierarchy.h>

namespace molprobity {
  namespace probe {

    /// @brief Structure to hold the results from a call to closest_contact()
    typedef struct ContactResult_ {
      Point   closest_contact;      ///< The point on the radius of the tested sphere closest to the dot
      double  dist_above_surface;   ///< Distance that the dot is above the tested sphere (negative for inside)
    } ContactResult;

    /// @brief Find the point of closest contact and distance above/below atomic surface
    /// @param [in] dot The dot whose location is to be projected onto the atom
    /// @param [in] atom The center of the atom
    /// @param [in] atom_radius The radius of the atom
    /// @return The point of closest contact on the atom's surface to the dot and the
    ///       signed distance that the dot is above that surface, with negative indicating
    ///       that it is inside the atom's surface.
    ContactResult closest_contact(Point dot, Point atom, double atom_radius);

    //=====================================================================================================

    /// @brief Class to handle scoring vectors of dots given two atoms.  Used by Reduce to compute energies.
    class AtomVsAtomDotScorer {
    public:
      /// @brief Constructor stores the non-bonded parameters to be used to determine atom features.
      /// @param [in] params Non-bonded parameters computed by the Geometry Restraints Manager.
      /// @param [in] gapWeight Factor to apply to gap between atoms
      /// @param [in] bumpWeight Factor to apply when atoms are in bumping overlap
      /// @param [in] hBondWeight Factor to apply to hydrogen-bond overlaps
      AtomVsAtomDotScorer(cctbx::geometry_restraints::nonbonded_params const& params
        , double gapWeight = 0.25
        , double bumpWeight = 10.0
        , double hBondWeight = 4.0
      ) : m_params(params), m_gapWeight(gapWeight), m_bumpWeight(bumpWeight), m_hBondWeight(hBondWeight) {};

      /// @brief Structure to hold the results from a call to score_dots()
      typedef struct ScoreDotsResult_ {
        double  bumpSubScore;       ///< Portion of the score due to bumps
        double  hBondSubScore;      ///< Portion of the score due to hydrogen bonds
        double  totalScore;         ///< Sum of the bumpSubScore and the hBondSubScore, the total score
      } ScoreDotsResult;

      /// @brief Determine the bump and hydrogen-bond subscores for a set of dots between two atoms
      /// @param [in] sourceAtom Atom that the dots are surrounding
      /// @param [in] targetAtom Atom to be tested against for bumping or hydrogen bonding.
      /// @param [in] dots Vector of dots to compare.  Each is added to the sourceAtom origin.
      /// @param [in] density Density of the dots on the probe sphere, used to normalize results
      /// @param [in] onlyBumps If true, ignore near touches and count even hydrogen bonds as bumps
      /// @return Normalized sum of scores, also broken down by hydrogen bond vs. bump scores.
      ScoreDotsResult score_dots(iotbx::pdb::hierarchy::atom sourceAtom, iotbx::pdb::hierarchy::atom targetAtom,
        std::vector<Point> dots, double density, bool onlyBumps = false);

      //===========================================================================
      // Seldom-used methods below here.

      /// @brief Test method to verify that the class is behaving as intended.
      /// @return Empty string on success, string telling what went wrong on failure.
      static std::string test();

    protected:
      /// Parameters stored from constructor.
      cctbx::geometry_restraints::nonbonded_params const& m_params;
      double m_gapWeight;
      double m_bumpWeight;
      double m_hBondWeight;
    };

    /// @todo Figure out what all of the things needed by Probe are.

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string Scoring_test();

  }
}
