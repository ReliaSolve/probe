##################################################################################
# This is a test program to validate that the Python wrapping of Probe worked.
#

import sys
import mmtbx_probe_ext as probe

def RunProbe():
  ret = 0
  print('Testing Probe objects')
  probe.DotSpheres_test()

  print('Making DotSphere from cache and getting its dots')
  cache = probe.DotSphereCache(10)
  sphere1 = cache.get_sphere(1)
  dots = sphere1.dots()
  print(' Found',len(dots),'dots')
  if len(dots) > 0:
    print('First dot is',dots[0])
    print('First dot is at',dots[0][0],dots[0][1],dots[0][2])

  return ret

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be the name of a PDB file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    fileName = sys.argv[i]

  sys.exit(RunProbe())
