// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#include <cmath>
#include <scitbx/constants.h>
#include <algorithm>
#include "DotSpheres.h"

namespace molprobity {
namespace probe {

DotSphere::DotSphere(double radius, double density)
  : m_rad(radius), m_dens(density)
{
  // Clamp our inputs to non-negative
  if (m_rad < 0) { m_rad = 0; }
  if (m_dens < 0) { m_dens = 0; }

  // If we have a zero radius or density, we're done
  if (m_rad == 0 || m_dens == 0) { return; }

  // Estimate the number of dots to be placed given the radius and density.
  // Do an overestimate.  This code is pulled from DotDph.cpp in Reduce.
  size_t num_dots = static_cast<size_t>(ceil(4.0 * scitbx::constants::pi * density * (radius * radius)));

  // Generate the dots, spreading them across the sphere.
  double offset = 0.2;
  double ang, cosang, sinang, phi, theta, xy0, x0, y0, z0;
  int nequator, nvert, nhoriz;
  bool odd = true;

  nequator = static_cast<int>(ceil(sqrt(num_dots * scitbx::constants::pi)));

  ang = 5.0 * scitbx::constants::pi / 360.0;
  cosang = cos(ang);
  sinang = sin(ang);

  nvert = nequator / 2;
  for (int j = 0; j <= nvert; j++) {
    phi = (scitbx::constants::pi * j) / nvert;
    z0 = cos(phi) * radius;
    xy0 = sin(phi) * radius;

    nhoriz = static_cast<int>(floor(nequator * sin(phi)));
    if (nhoriz < 1) { nhoriz = 1; }
    for (int k = 0; k < nhoriz; k++) {
      if (odd) { theta = (2.0 * scitbx::constants::pi * k + offset) / nhoriz; }
      else { theta = (2.0 * scitbx::constants::pi * k) / nhoriz; }
      x0 = cos(theta) * xy0;
      y0 = sin(theta) * xy0;

      m_vec.push_back(molprobity::probe::Dot(x0, y0 * cosang - z0 * sinang, y0 * sinang + z0 * cosang));
    }
    odd = !odd;
  }
}

std::string DotSphere::test()
{
  // Test creation with negative density and/or radius
  {
    DotSphere d(-1, 10);
    if (d.radius() != 0 || d.dots().size() != 0) {
      return "molprobity::probe::DotSphere::test(): Construction with negative radius failed";
    }
  }
  {
    DotSphere d(5, -3);
    if (d.density() != 0 || d.dots().size() != 0) {
      return "molprobity::probe::DotSphere::test(): Construction with negative density failed";
    }
  }
  {
    DotSphere d(-5, -3);
    if (d.radius() != 0 || d.density() != 0 || d.dots().size() != 0) {
      return "molprobity::probe::DotSphere::test(): Construction with negative radius and density failed";
    }
  }

  // Test creating with extremely small density to be sure we get a single dot
  {
    DotSphere d(1, 1e-10);
    if (d.radius() != 1 || d.density() != 1e-10 || d.dots().size() > 2) {
      return "molprobity::probe::DotSphere::test(): Construction with small density failed";
    }
  }

  // Test creating with a different densities and see whether the number of dots is close.
  {
    for (double den = 1; den < 128; den *= 1.3) {
      DotSphere d(1, den);
      double expected = ceil(4.0 * scitbx::constants::pi * den * (1 * 1));
      long found = static_cast<long>(d.dots().size());
      double diff = abs(found - expected);
      if (diff > std::max(3.0, expected*0.05)) {
        return std::string("molprobity::probe::DotSphere::test(): Construction with density")
          +std::to_string(den) + " failed: "
          + "found " + std::to_string(found) + " dots but expected " + std::to_string(expected);
      }
    }
  }

  // Test creating with reasonable density to be sure we get dots in all octants
  {
    DotSphere d(1, 5);
    bool px = false, mx = false, py = false, my = false, pz = false, mz = false;
    const std::vector<Dot>& dots = d.dots();
    for (size_t i = 0; i < dots.size(); i++) {
      Dot dot = dots[i];
      if (dot[0] > 0) { px = true; }
      if (dot[0] < 0) { mx = true; }
      if (dot[1] > 0) { py = true; }
      if (dot[1] < 0) { my = true; }
      if (dot[2] > 0) { pz = true; }
      if (dot[2] < 0) { mz = true; }
    }
    if (!px || !mx || !py || !my || !pz || !mz) {
      return "molprobity::probe::DotSphere::test(): Construction with reasoneble density "
        "did not have dots in all octants";
    }
  }

  // All tests passed.
  return "";
}

const DotSphere& DotSphereCache::get_sphere(double radius)
{
  std::map<double, DotSphere>::const_iterator ret = m_spheres.find(radius);
  if (ret == m_spheres.end()) {
    // We don't have a sphere with this radius -- create one and insert it
    std::pair<std::map<double, DotSphere>::iterator, bool> iRet = 
      m_spheres.emplace(radius, DotSphere(radius, m_dens));
    ret = iRet.first;
  }
  return ret->second;
}

std::string DotSphereCache::test()
{
  // Object to use for our tests
  DotSphereCache dsc(10);

  // Test creation of a single sphere
  const DotSphere& sp1 = dsc.get_sphere(1.0);
  if (dsc.size() != 1) {
    return "molprobity::probe::DotSphereCache::test(): Single sphere creation failed";
  }

  // Ask for another sphere of the same size and make sure we get the same one.
  const DotSphere& sp2 = dsc.get_sphere(1.0);
  if (dsc.size() != 1 || sp1 != sp2) {
    return "molprobity::probe::DotSphereCache::test(): Identical sphere creation failed";
  }

  // Ask a sphere of a differnt size and make sure we get a different one.
  const DotSphere& sp3 = dsc.get_sphere(2.0);
  if (dsc.size() != 2 || sp1 == sp3) {
    return "molprobity::probe::DotSphereCache::test(): Unique sphere creation failed";
  }

  // All tests passed.
  return "";
}

std::string DotSpheres_test()
{
  std::string ret;

  /// Test DotSphere class
  ret = DotSphere::test();
  if (!ret.empty()) {
    return std::string("molprobity::probe::DotSpheres_test(): failed: ") + ret;
  }

  /// Test DotSphereCache class
  ret = DotSphereCache::test();
  if (!ret.empty()) {
    return std::string("molprobity::probe::DotSpheres_test(): failed: ") + ret;
  }

  // All tests passed.
  return "";
}


} // end namespace probe
} // end namespace molprobity 
