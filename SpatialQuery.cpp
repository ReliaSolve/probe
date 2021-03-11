// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#include <cmath>
#include <scitbx/constants.h>
#include <algorithm>
#include "SpatialQuery.h"

namespace molprobity {
  namespace probe {

const Coord molprobity::probe::SpatialQuery::DEFAULT_BIN_SIZE = 3;  ///< Default size of a grid bin in X, Y, and Z.

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
  Point binSize(DEFAULT_BIN_SIZE, DEFAULT_BIN_SIZE, DEFAULT_BIN_SIZE);
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
    else { lowerIndex[i] = static_cast<size_t>(floor((lower - m_lowerBounds[i]) / m_binSize[i])); }
    if (lowerIndex[i] >= m_gridSize[i]) { lowerIndex[i] = m_gridSize[i] - 1; }

    Coord upper = p[i] + max_distance;
    if (upper < m_lowerBounds[i]) { upperIndex[i] = 0; }
    else { upperIndex[i] = static_cast<size_t>(floor((upper - m_lowerBounds[i]) / m_binSize[i])); }
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
          if ((dist2 >= min2) && (dist2 <= max2)) {
            ret.push_back(a);
          }
        }
      }
    }
  }

  return ret;
}

//============================================================================================
// Test data and code below here.

std::string SpatialQuery::test()
{
  // Construct a model with some atoms in it so that we can use it for our tests.
  // Make it a bunch of Hydrogen atoms in a regular 3D grid
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

  // Test creation of a grid and see if it is the size that we expect it to be
  {
    Point lower(-10, -10, -10);
    Point upper( 11,  11,  11);
    Point binSize(5, 5, 5);
    SpatialQuery q(lower, upper, binSize);
    if (q.m_lowerBounds != lower) {
      return "molprobity::probe::SpatialQuery::test(): Unexpected lower bound on standard construction";
    }
    if (q.m_gridSize != std::array<size_t, 3>({ 5, 5, 5 })) {
      return "molprobity::probe::SpatialQuery::test(): Unexpected grid size on standard construction";
    }
  }

  // Test creation of a grid with upper below lower
  {
    Point lower(-10, -10, -10);
    Point upper(-11, -11, -11);
    Point binSize(5, 5, 5);
    SpatialQuery q(lower, upper, binSize);
    if (q.m_lowerBounds != lower) {
      return "molprobity::probe::SpatialQuery::test(): Unexpected lower bound on inverted construction";
    }
    if (q.m_gridSize != std::array<size_t, 3>({ 1, 1, 1 })) {
      return "molprobity::probe::SpatialQuery::test(): Unexpected grid size on inverted construction";
    }
  }

  // Test adding and removing elements, making sure we cannot re-insert or
  // re-remove the same element
  {
    Point lower(-10, -10, -10);
    Point upper(11, 11, 11);
    Point binSize(5, 5, 5);
    SpatialQuery q(lower, upper, binSize);

    iotbx::pdb::hierarchy::atom a(lower, lower);
    if (!q.add(a)) {
      return "molprobity::probe::SpatialQuery::test(): Could not add atom";
    }
    if (q.add(a)) {
      return "molprobity::probe::SpatialQuery::test(): Could double add atom";
    }
    if (!q.remove(a)) {
      return "molprobity::probe::SpatialQuery::test(): Could not remove atom";
    }
    if (q.remove(a)) {
      return "molprobity::probe::SpatialQuery::test(): Could double remove atom";
    }
  }

  // Test checking for neighbors in the case where they are in our bin.
  {
    Point lower(-10, -10, -10);
    Point upper(11, 11, 11);
    Point binSize(5, 5, 5);
    SpatialQuery q(lower, upper, binSize);
    atom_less lessThan; ///< Comparison object to see if one atom is less than another

    iotbx::pdb::hierarchy::atom a(lower, lower), b(upper, upper);
    q.add(a);
    q.add(b);

    std::vector<iotbx::pdb::hierarchy::atom> n1 = q.neighbors(lower, 0, 1);
    if ((n1.size() != 1) || (lessThan(n1[0],a)) || (lessThan(a, n1[0]))) {
      return "molprobity::probe::SpatialQuery::test(): Did not find expected lower neighbor";
    }

    std::vector<iotbx::pdb::hierarchy::atom> n2 = q.neighbors(upper, 0, 1);
    if ((n2.size() != 1) || (lessThan(n2[0], b)) || (lessThan(b, n2[0]))) {
      return "molprobity::probe::SpatialQuery::test(): Did not find expected upper neighbor";
    }

    // Check for the case where we are right on top of them to be sure they are excluded if
    // our minimum distance is greater than 0.
    std::vector<iotbx::pdb::hierarchy::atom> n3 = q.neighbors(lower, 0.1, 1);
    if (n3.size() != 0) {
      return "molprobity::probe::SpatialQuery::test(): Found lower neighbor when not expected";
    }
  }

  // Test the hierarchy-based constructor, which will be used for later tests as well.
  {
    SpatialQuery q(m);
    if (q.m_lowerBounds != Point(0, 0, 0)) {
      return "molprobity::probe::SpatialQuery::test(): Unexpected lower bound on model-based construction";
    }
    size_t expectedSize = static_cast<size_t>((numAtoms - 1) * spacing / DEFAULT_BIN_SIZE);
    if (q.m_gridSize != std::array<size_t, 3>({ expectedSize, expectedSize, expectedSize })) {
      return "molprobity::probe::SpatialQuery::test(): Unexpected grid size on model-based construction";
    }

    // Test checking for neighbors when they are in all bins.
    // Test checking for neighbors when the point is inside the grid and when it is
    // outside of the grid on any side.
    for (int x = -90; x <= 110; x += 100) {
      for (int y = -90; y <= 110; y += 100) {
        for (int z = -90; z <= 110; z += 100) {
          Point p(x, y, z);
          std::vector<iotbx::pdb::hierarchy::atom> n = q.neighbors(p, 0, 2000);
          if (n.size() != numAtoms * numAtoms * numAtoms) {
            return std::string("molprobity::probe::SpatialQuery::test(): Wrong number of neighbors in full-grid test at ")
              + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z)
              + " (found " + std::to_string(n.size()) + ", expected " + std::to_string(numAtoms * numAtoms * numAtoms) + ")";
          }
        }
      }
    }

    // Test checking for neighbors when they are in multiple bins.  Check at one corner at slightly over the
    // requested spacing so we should get us and our three neighbors.
    {
      Point p(0, 0, 0);
      std::vector<iotbx::pdb::hierarchy::atom> n = q.neighbors(p, 0, spacing+0.1);
      if (n.size() != 4) {
        return "molprobity::probe::SpatialQuery::test(): Wrong number of neighbors in multi-element test";
      }
    }
  }

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
