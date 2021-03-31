// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#pragma once

#include <vector>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>

namespace molprobity {
  namespace probe {

    /// The type to be used for a single coordinate in space
    typedef double Coord;

    /// @brief A location in space
    ///
    /// This is a class rather than a typedef because the typedef confused
    /// Boost Python when we tried provide an element accessor to it; it claimed
    /// that we had already registered a to-Python converter, but when we did
    /// not add one it failed because none had been defined...
    class Point : public scitbx::vec3<Coord> {
    public:
      Point(scitbx::vec3<Coord> const &c) : scitbx::vec3<Coord>(c[0],c[1],c[2]) { }
      Point(Coord x, Coord y, Coord z) : scitbx::vec3<Coord>(x,y,z) { }
      Point() : scitbx::vec3<Coord>() { }
    };

  }
}
