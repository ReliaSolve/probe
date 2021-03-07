// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#include "DotSpheres.h"
#include <iostream>

int main(int argc, const char* argv[])
{
  // Test DotSpheres module.
  std::string ret = molprobity::probe::DotSpheres_test();
  if (!ret.empty()) {
    std::cerr << "Error: " << ret << std::endl;
    return 1;
  }

  std::cout << "Success!" << std::endl;
  return 0;
}
