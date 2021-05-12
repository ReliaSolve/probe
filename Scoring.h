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
      Point   closestContact;     ///< The point on the radius of the tested sphere closest to the dot
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

    /// @brief This is the distance from the point on the source atom closest to the target to the edge
    ///        of contact when the source and target are in optimal contact?
    ///        @todo Figure out for sure what this is computing.
    /// @param [in] ra Radius of the source atom
    /// @param [in] rb Radius of the target atom
    /// @param [in] rp Probe radius
    double kissEdge2bullsEye(double ra, double rb, double rp);

    /// @brief A dot is annular if it is further from the center of contact than edge of the overlap
    /// region is at optimum contact.
    /// 
    /// This checks to make sure that dots that would not have contributed to a good score at optimium
    /// contact are not considered to contribute to a good score when the atoms are overlapping or far from
    /// each other?
    /// @todo Figure out for sure what this is computing.
    bool annularDots(const Point& dot, const Point& srcLoc, double srcVDWRad,
      const Point& targLoc, double targVDWRad, double probeRadius);


    //=====================================================================================================
    /// @brief Class to hold data values for an atom beyond those present in the hierarchy::atom class itself
    // that are needed by the Probe calculations.  These must be filled in by the client, perhaps using data
    // from the mmtbx.monomer_library.server.ener_lib() function to get a library and looking things up
    // in it based on getting a monomer lib query mon_lib_query() from reduce_hydrogen, then calling its
    // atom_dict() method and then looking up the atom by name in that dictionary to get its type_energy
    // value.

    class ExtraAtomInfo {
    public:
      /// @brief Constructor with default parameters
      ExtraAtomInfo(double vdwRadius = 0, bool isAcceptor = false, bool isDonor = false,
        bool isDummyHydrogen = false)
        : m_vdwRadius(vdwRadius), m_isAcceptor(isAcceptor), m_isDonor(isDonor)
        , m_isDummyHydrogen(isDummyHydrogen) {}

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

    /// @brief Class to handle scoring dots given sets of atoms.
    class DotScorer {
    public:
      /// @brief Constructor stores the non-bonded parameters to be used to determine atom features.
      /// @param [in] extraInfo Vector of extra information pertaining to each atom in the structure
      ///         that is being examined (the one that was used to construct the SpatialQuery
      ///         structure).  Warning: The i_sel values from the atoms in the structure and the atoms
      ///         passed as parameters are used to look up directly in this vector so they must not
      ///         have changed (due to structure modification) since the extaInfo vector or the
      ///         SpatialQuery structure were generated.
      /// @param [in] gapScale Scale factor to apply to gap between atoms (gap is divided by this)
      /// @param [in] bumpWeight Factor to apply when atoms are in bumping overlap
      /// @param [in] hBondWeight Factor to apply to hydrogen-bond overlaps
      /// @param [in] maxRegularHydrogenOverlap How much overlap can there be between a hydrogen
      ///             and the atom it is hydrogen-bonded to before we call it a clash when the atoms
      ///             are not both charged.  It must go badBumpOverlap beyond this before we call
      ///             it a bad clash.
      /// @param [in] maxChargedHydrogenOverlap How much overlap can there be between a hydrogen
      ///             and the atom it is hydrogen-bonded to before we call it a clash when the
      ///             atoms are both charged.  It must go badBumpOverlap beyond this before we call
      ///             it a bad clash.
      /// @param [in] badBumpOverlap Atoms that overlap more than this will cause bad bump to be flagged.
      DotScorer(scitbx::af::shared<ExtraAtomInfo> extraInfo
        , double gapScale = 0.25
        , double bumpWeight = 10.0
        , double hBondWeight = 4.0
        , double maxRegularHydrogenOverlap = 0.6
        , double maxChargedHydrogenOverlap = 0.8
        , double badBumpOverlap = 0.4
      ) : m_extraInfo(extraInfo)
        , m_gapScale(gapScale), m_bumpWeight(bumpWeight), m_hBondWeight(hBondWeight)
        , m_maxRegularHydrogenOverlap(maxRegularHydrogenOverlap)
        , m_maxChargedHydrogenOverlap(maxChargedHydrogenOverlap)
        , m_badBumpOverlap(badBumpOverlap) {};

      /// @brief Structure to hold the results from a call to check_dot()
      class CheckDotResult {
      public:
        int     overlapType = -2;           ///< -2 for no result, -1 for clash, 0 for near contact, 1 for hydrogen bond
        iotbx::pdb::hierarchy::atom cause;  ///< Cause of the overlap, if overlapType != -2
        double  overlap = 0;                ///< Amount of overlap if there is a clash
        double  gap = 1e100;                ///< Gap distance (overlap may only be a fraction of this).
        bool    annular = false;            ///< Was this an annular dot?
      };

      /// @brief Score an individual dot against a specific set of interacting atoms unless within an excluded atom
      /// @param [in] sourceAtom Atom that the dot is offset with respect to.
      /// @param [in] dotOffset Offset from the center of the source atom to the dot
      /// @param [in] probeRadius Radius of the probe rolled between the two potentially-contacting atoms
      ///             If this is < 0, an invalid result will be returned.
      /// @param [in] interacting The atoms that are to be checked because they are close enough to sourceAtom
      ///             to interact with it.
      /// @param [in] excluded Atoms that are to be excluded from contact, for example this could be a list
      ///             of atoms bonded to sourceAtom.  If the dot is inside an excluded atom, it will not be
      ///             considered even if it is overlapping with an interacting atom.
      CheckDotResult check_dot(iotbx::pdb::hierarchy::atom sourceAtom, ExtraAtomInfo const& sourceExtra,
        Point const& dotOffset, double probeRadius,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& interacting,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& exclude);

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
      /// @param [in] minOccupancy The minimum occupancy of the atom to be scored (applies to both
      ///             source and target; if either is below this, it will not be scored.
      ///             If the source is below this, the return value will be marked valid but will
      ///             have 0 in all of its interactions.
      /// @param [in] spatialQuery Structure to ask for neighbors of the atom.  This must contain
      ///             only atoms that are to be considered; those that are in the same conformation
      ///             or in all conformations.
      /// @param [in] nearbyRadius Maximum distance that an atom can be away and still be a neighbor.
      ///             This should NOT include consideration of the probe radius, which will be added
      ///             inside this function when it is required (depending on whether onlyBumps is
      ///             set).
      /// @param [in] probeRadius Radius of the probe rolled between the two potentially-contacting atoms
      ///             If this is < 0, an invalid result will be returned.
      /// @param [in] excluded Atoms that are to be excluded from contact, for example this could be a list
      ///             of atoms bonded to sourceAtom.
      /// @param [in] dots Vector of dots to compare.  Each is added to the sourceAtom origin.
      /// @param [in] density Density of the dots on the probe sphere, used to normalize results.
      ///             If this is <= 0, an invalid result will be returned.
      /// @param [in] onlyBumps If true, ignore near touches and count even hydrogen bonds as bumps.
      /// @return Normalized sum of scores, also broken down by hydrogen bond vs. bump scores.
      ScoreDotsResult score_dots(iotbx::pdb::hierarchy::atom sourceAtom, double minOccupancy,
        SpatialQuery &spatialQuery, double nearbyRadius, double probeRadius,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &exclude,
        scitbx::af::shared<Point> const &dots, double density, bool onlyBumps = false);

      //===========================================================================
      // Seldom-used methods below here.

      /// @brief Test method to verify that the class is behaving as intended.
      /// @return Empty string on success, string telling what went wrong on failure.
      static std::string test();

    protected:
      /// Parameters stored from constructor.
      scitbx::af::shared<ExtraAtomInfo> m_extraInfo;
      double m_gapScale;
      double m_bumpWeight;
      double m_hBondWeight;
      double m_maxRegularHydrogenOverlap;
      double m_maxChargedHydrogenOverlap;
      double m_badBumpOverlap;
    };

    /// @todo Figure out what all of the things needed by Probe (as opposed to Reduce) are.

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string Scoring_test();

  }
}
