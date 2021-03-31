// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#pragma once

#include "SpatialQuery.h"
#include <cctbx/geometry_restraints/nonbonded.h>
#include <iotbx/pdb/hierarchy.h>

namespace molprobity {
  namespace probe {

    //=====================================================================================================
    // Helper functions outside the class.

    /// @brief Structure to hold the results from a call to closest_contact()
    class ContactResult {
    public:
      Point   closestContact;      ///< The point on the radius of the tested sphere closest to the dot
      double  distAboveSurface;   ///< Distance that the dot is above the tested sphere (negative for inside)
    };

    /// @brief Find the point of closest contact and distance above/below atomic surface
    /// @param [in] dot The dot whose location is to be projected onto the atom
    /// @param [in] atom The center of the atom
    /// @param [in] atom_radius The radius of the atom
    /// @return The point of closest contact on the atom's surface to the dot and the
    ///       signed distance that the dot is above that surface, with negative indicating
    ///       that it is inside the atom's surface.
    ContactResult closest_contact(Point dot, Point atom, double atom_radius);

    /// @brief Return the signed-integer charge of the atom.
    /// @param atom Atom whose charge is to be determined.
    /// return Integer representing the charge: -2, -1, 0, 1, 2
    int atom_charge(iotbx::pdb::hierarchy::atom const& atom);

    // functions used to restrict annular rings of good dots around clashes

    /// @brief The distance from a dot to the point on the source surface that is closest to the target
    /// @param [in] dot The dot being considered
    /// @param [in] srcLoc The center of the source atom
    /// @param [in] srcVDWRad The van Der Waals radius of the source atom in Angstroms
    /// @param [in] targLoc The center of the target atom
    /// @return The distance from a dot to the point on the source surface that is closest to the target
    double dot2srcCenter(const Point& dot, const Point& srcLoc, double srcVDWRad, const Point& targLoc);

    /// @brief @todo Figure out what this is computing.
    double kissEdge2bullsEye(double ra, double rb, double rp);

    /// @brief @todo Figure out what this is computing.
    bool annularDots(const Point& dot, const Point& srcLoc, double srcVDWRad,
      const Point& targLoc, double targVDWRad, double probeRadius);


    //=====================================================================================================
    /// @brief Class to hold data values for an atom beyond those present in the hierarchy::atom class itself
    // that are needed by the Probe calculations.

    class ExtraAtomInfo {
    public:
      /// @brief Get and set methods
      double  getVdwRadius() const { return m_vdwRadius; }
      void setVdwRadius(double val) { m_vdwRadius = val; }

      bool getIsAcceptor() const { return m_isAcceptor; }
      void setIsAcceptor(bool val) { m_isAcceptor = val; }
      bool getIsDonor() const { return m_isDonor; }
      void setIsDonor(bool val) { m_isDonor = val; }
      bool getIsDummyHydrogen() const { return m_isDummyHydrogen; }
      void setIsDummyHydrogen(bool val) { m_isDummyHydrogen = val; }

      /// @brief == operator is required so that we can wrap the standard vector operators in Boost::Python
      bool operator ==(ExtraAtomInfo const& o) {
        return ((getVdwRadius() == o.getVdwRadius())
          && (getIsAcceptor() == o.getIsAcceptor())
          && (getIsDonor() == o.getIsDonor())
          && (getIsDummyHydrogen() == o.getIsDummyHydrogen()));
      }

    protected:
      double m_vdwRadius = 0;          ///< van Der Waals radius of the atom
      bool m_isAcceptor = false;       ///< Does this accept hydrogen bonds (aromatic carbon, nitrogen acceptor,
                                       ///  oxygen, sulfur, fluorine, chlorine, bromine, or iodine?
      bool m_isDonor = false;          ///< Is this a donor hydrogen (from polar, aromatic polar, or water)?
      bool m_isDummyHydrogen = false;  ///< These are inserted on Oxygens that are waters to provide
                                       ///  bonds that can go in any direction.
    };

    //=====================================================================================================

    /// @brief Class to handle scoring vectors of dots given two atoms.  Used by Reduce to compute energies.
    class AtomVsAtomDotScorer {
    public:
      /// @brief Constructor stores the non-bonded parameters to be used to determine atom features.
      /// @param [in] extraInfo Vector of extra information pertaining to each atom in the structure
      ///         that is being examined (the one that was used to construct the SpatialQuery
      ///         structure).  Warning: The i_sel values from the atoms in the structure and the atoms
      ///         passed as parameters are used to look up directly in this vector so they must not
      ///         have changed (due to structure modification) since the extaInfo vector or the
      ///         SpatialQuery structure were generated.
      /// @param [in] gapWeight Factor to apply to gap between atoms
      /// @param [in] bumpWeight Factor to apply when atoms are in bumping overlap
      /// @param [in] hBondWeight Factor to apply to hydrogen-bond overlaps
      /// @param [in] minRegularHydrogenBondGap How much overlap can there be between a hydrogen
      ///             and the atom it is hydrogen-bonded to before we call it a clash when the atoms
      ///             are not both charged.
      /// @param [in] minChargedHydrogenBondGap How much overlap can there be between a hydrogen
      ///             and the atom it is hydrogen-bonded to before we call it a clash when the
      ///             atoms are both charged.
      /// @param [in] badBumpBondGap Dots that are closer than this will cause bad bump to be flagged.
      AtomVsAtomDotScorer(scitbx::af::shared<ExtraAtomInfo> extraInfo
        , double gapWeight = 0.25
        , double bumpWeight = 10.0
        , double hBondWeight = 4.0
        , double minRegularHydrogenBondGap = 0.6
        , double minChargedHydrogenBondGap = 0.8
        , double badBumpBondGap = 0.4
      ) : m_extraInfo(extraInfo)
        , m_gapWeight(gapWeight), m_bumpWeight(bumpWeight), m_hBondWeight(hBondWeight)
        , m_minRegularHydrogenBondGap(minRegularHydrogenBondGap)
        , m_minChargedHydrogenBondGap(minChargedHydrogenBondGap)
        , m_badBumpBondGap(badBumpBondGap) {};

      /// @brief Structure to hold the results from a call to score_dots()
      class ScoreDotsResult {
      public:
        bool    valid = false;          ///< False if this information has not yet been computed.
        double  bumpSubScore = 0;       ///< Portion of the score due to bumps
        double  hBondSubScore = 0;      ///< Portion of the score due to hydrogen bonds
        double  attractSubScore = 0;    ///< Portion of the score due to non-bumping attraction
        bool    hasBadBump = false;     ///< Did this atom have a bad bump for any of its dots?
        /// @brief Sum of all of the sub-scores, the total score
        double  totalScore() const { return bumpSubScore + hBondSubScore + attractSubScore; }
      };

      /// @brief Determine the bump and hydrogen-bond subscores for a vector of dots on an atom
      /// @param [in] sourceAtom Atom that the dots are surrounding.
      /// @param [in] minOccupancy The minimum occupancy of the atom to be scored.
      /// @param [in] spatialQuery Structure to ask for neighbors of the atom.
      /// @param [in] nearbyRadius Maximum distance that an atom can be away and still be a neighbor
      /// @param [in] probeRadius Radius of the probe rolled between the two potentially-contacting atoms
      /// @param [in] excluded Atoms that are to be excluded from contact.
      /// @param [in] dots Vector of dots to compare.  Each is added to the sourceAtom origin.
      /// @param [in] density Density of the dots on the probe sphere, used to normalize results.
      /// @param [in] onlyBumps If true, ignore near touches and count even hydrogen bonds as bumps.
      /// @return Normalized sum of scores, also broken down by hydrogen bond vs. bump scores.
      ScoreDotsResult score_dots(iotbx::pdb::hierarchy::atom sourceAtom, double minOccupancy,
        SpatialQuery spatialQuery, double nearbyRadius, double probeRadius,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude,
        scitbx::af::shared<Point> dots, double density, bool onlyBumps = false);

      //===========================================================================
      // Seldom-used methods below here.

      /// @brief Test method to verify that the class is behaving as intended.
      /// @return Empty string on success, string telling what went wrong on failure.
      static std::string test();

    protected:
      /// Parameters stored from constructor.
      scitbx::af::shared<ExtraAtomInfo> m_extraInfo;
      double m_gapWeight;
      double m_bumpWeight;
      double m_hBondWeight;
      double m_minRegularHydrogenBondGap;
      double m_minChargedHydrogenBondGap;
      double m_badBumpBondGap;
    };

    /// @todo Figure out what all of the things needed by Probe are.

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string Scoring_test();

  }
}
