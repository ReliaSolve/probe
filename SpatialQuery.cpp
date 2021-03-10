// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#include <cmath>
#include <scitbx/constants.h>
#include <algorithm>
#include "SpatialQuery.h"

namespace molprobity {
  namespace probe {

void SpatialQuery::initialize(Point lowerBounds, Point upperBounds, Point binSize)
{
  // Set the bin size
  m_binSize = binSize;

  // Determine the proper location and size for the grid.
  for (size_t i = 0; i < 3; i++) {

    // Upper bounds must be at or above lower bounds
    double val = lowerBounds[i];
    m_lowerBounds[i] = val;
    if (upperBounds[i] < val) { upperBounds[i] = val; }

    // If the bin size is less than or equal to zero, set it to 1.
    if (binSize[i] <= 0) { binSize[i] = 1; }
    m_binSize[i] = binSize[i];

    // Find the grid size on each axis, which must be large enough to span the whole
    // range and must have at least one entry on each axis.
    m_gridSize[i] = static_cast<size_t>(ceil((upperBounds[i] - lowerBounds[i]) / m_binSize[i]));
    if (m_gridSize[i] == 0) { m_gridSize[i] = 1; }
  }

  // Allocate the space for the grid.  It will start with empty vectors for each element
  m_grid.resize(m_gridSize[0] * m_gridSize[1] * m_gridSize[2]);
}

SpatialQuery::SpatialQuery(Point lowerBounds, Point upperBounds, Point binSize)
{
  initialize(lowerBounds, upperBounds, binSize);
}

SpatialQuery::SpatialQuery(iotbx::pdb::hierarchy::model m)
{
  // Compute the parameters needed and initialize the grid
  Point lowerBounds(1e10, 1e10, 1e10);
  Point upperBounds(-1e10, -1e10, -1e10);
  for (auto chain : m.chains()) {
    for (auto residue_group : chain.residue_groups()) {
      for (auto atom_group : residue_group.atom_groups()) {
        for (iotbx::pdb::hierarchy::atom const& a : atom_group.atoms()) {
          for (size_t i = 0; i < 3; i++) {
            Point loc = a.data->xyz;
            if (loc[i] < lowerBounds[i]) { lowerBounds[i] = loc[i]; }
            if (loc[i] > upperBounds[i]) { upperBounds[i] = loc[i]; }
          }
        }
      }
    }
  }

  // Make sure that we don't have more than 50 bins on a given axis
  Point binSize(3,3,3);
  for (size_t i = 0; i < 3; i++) {
    Coord minSize = (upperBounds[i] - lowerBounds[i]) / 50;
    if (binSize[i] < minSize) { binSize[i] = minSize; }
  }
  initialize(lowerBounds, upperBounds, binSize);

  // Add all of the atoms in the model
  for (auto chain : m.chains()) {
    for (auto residue_group : chain.residue_groups()) {
      for (auto atom_group : residue_group.atom_groups()) {
        for (iotbx::pdb::hierarchy::atom const& a : atom_group.atoms()) {
          add(a);
        }
      }
    }
  }
}

bool SpatialQuery::add(iotbx::pdb::hierarchy::atom a)
{
  // Look up the index of the grid element for the location of the atom.
  // Attempt to insert the atom at that location.
  // Return whether the insertion was a success.
  return m_grid[grid_index(a.data->xyz)].insert(a).second;
}

bool SpatialQuery::remove(iotbx::pdb::hierarchy::atom a)
{
  // Look up the index of the grid element for the location of the atom.
  // Attempt to remove the atom at that location.
  // Return whether the insertion was a success (a 1 means a success).
  return m_grid[grid_index(a.data->xyz)].erase(a) == 1;
}

std::vector<iotbx::pdb::hierarchy::atom> SpatialQuery::neighbors(
  Point const& p, double min_distance, double max_distance)
{
  // Find the range of bins that are within the maximum distance of the point.
  // Handle cases where the point is at any location relative to the grid: inside,
  // outside to left or right on one or more axes.
  // We overestimate here, counting all bins within range in X, Y, or Z.
  std::array<size_t, 3> lowerIndex, upperIndex;
  for (size_t i = 0; i < 3; i++) {
    Coord lower = p[i] - max_distance;
    if (lower < m_lowerBounds[i]) { lowerIndex[i] = 0; }
    else { lowerIndex[i] = static_cast<size_t>(floor((p[i] - m_lowerBounds[i]) / m_binSize[i])); }
    if (lowerIndex[i] >= m_gridSize[i]) { lowerIndex[i] = m_gridSize[i] - 1; }

    Coord upper = p[i] + max_distance;
    if (upper < m_lowerBounds[i]) { upperIndex[i] = 0; }
    else { upperIndex[i] = static_cast<size_t>(floor((p[i] - m_lowerBounds[i]) / m_binSize[i])); }
    if (upperIndex[i] >= m_gridSize[i]) { upperIndex[i] = m_gridSize[i] - 1; }
  }

  // Look through each atom in each potential bin and add it to the result if
  // it is within the distance range.  We square the distances ranges and compare
  // against these to avoid having to do square roots.
  double min2 = min_distance * min_distance;
  double max2 = max_distance * max_distance;
  std::vector<iotbx::pdb::hierarchy::atom> ret;
  for (size_t x = lowerIndex[0]; x <= upperIndex[0]; x++) {
    for (size_t y = lowerIndex[1]; y <= upperIndex[1]; y++) {
      for (size_t z = lowerIndex[2]; z <= upperIndex[2]; z++) {
        size_t index = x + m_gridSize[0] * (y + m_gridSize[1] * (z));
        for (iotbx::pdb::hierarchy::atom const &a : m_grid[index]) {
          double dist2 = (a.data->xyz - p).length_sq();
          if (dist2 >= min2 && dist2 <= max2) {
            ret.push_back(a);
          }
        }
      }
    }
  }

  return ret;
}

std::string SpatialQuery::test()
{
  // Test creation of a grid and see if it is the size that we expect it to be
  /// @todo

  // Test adding and removing elements, making sure we cannot re-insert or
  // re-remove the same element
  /// @todo

  // Test checking for neighbors in the case where they are in our bin.
  /// @todo

  // Test checking for neighbors when they are in multiple bins.
  /// @todo

  // Test checking for neighbors when the point is outside of the grid on any
  // side.
  /// @todo


  // Test the hierarchy-based constructor
  /// @todo

  // All tests passed.
  return "";
}

std::string SpatialQuery_test()
{
  std::string ret;

  /// Test SpatialQuery class
  ret = SpatialQuery::test();
  if (!ret.empty()) {
    return std::string("molprobity::probe::SpatialQuery_test(): failed: ") + ret;
  }

  // All tests passed.
  return "";
}


} // end namespace probe
} // end namespace molprobity 
