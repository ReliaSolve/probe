##################################################################################
# This is a test program to validate that the Python wrapping of Probe worked.
#

import sys
import mmtbx_probe_ext as probe

def RunProbeTests():

  # Make sure we can get at the DotSphere objects and their methods
  print('Making DotSphere from cache and getting its dots')
  cache = probe.DotSphereCache(10)
  sphere1 = cache.get_sphere(1)
  dots = sphere1.dots()
  print(' Found',len(dots),'dots')
  print('First dot is at',dots[0].elems[0],dots[0].elems[1],dots[0].elems[2])

  print('Testing DotSphere objects')
  ret = probe.DotSpheres_test()
  if len(ret) > 0:
    return ret

  print('Testing SpatialQuery objects')
  ret = probe.SpatialQuery_test()
  if len(ret) > 0:
    return ret

  return ret

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  ret = RunProbeTests()
  if len(ret) == 0:
    print('Success!')
  sys.exit(ret)
