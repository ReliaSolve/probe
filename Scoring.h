// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#pragma once

#include "Common.h"

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

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string Scoring_test();
  }
}
