##################################################################################
# This is a test program to validate that the Python wrapping of Probe worked.
#

import sys
import mmtbx_probe_ext as probe
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
import mmtbx

def RunProbeTests():

  #========================================================================
  # Make sure we can get at the DotSphere objects and their methods
  print('Making DotSphere from cache and getting its dots')
  cache = probe.DotSphereCache(10)
  sphere1 = cache.get_sphere(1)
  dots = sphere1.dots()
  print(' Found',len(dots),'dots')
  print('First dot is at',dots[0].elems[0],dots[0].elems[1],dots[0].elems[2])

  #========================================================================
  # Make sure we can fill in an ExtraAtomInfoList and pass it to scoring
  # Generate an example data model with a small molecule in it
  print('Filling in extra atom information needed for probe score')
  mmm=map_model_manager()         #   get an initialized instance of the map_model_manager
  mmm.generate_map()              #   get a model from a generated small library model and calculate a map for it
  model = mmm.model()             #   get the model

  # Fill in an ExtraAtomInfoList with an entry for each atom in the hierarchy.
  # We first find the largest i_seq sequence number in the model and reserve that
  # many entries so we will always be able to fill in the entry for an atom.
  atoms = model.get_atoms()
  print(' Found',len(atoms),'atoms')
  maxI = atoms[0].i_seq
  for a in atoms:
    if a.i_seq > maxI:
      maxI = a.i_seq
  extra = []
  for i in range(maxI+1):
    extra.append(probe.ExtraAtomInfo())

  # Traverse the hierarchy and look up the extra data to be filled in.
  mon_lib_srv = model.get_mon_lib_srv()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  ph = model.get_hierarchy()
  for m in ph.models():
    for chain in m.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          md, ani = mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
                residue_name=ag.resname, atom_names=ag.atoms().extract_name())
          atom_dict = md.atom_dict()

          for a in ag.atoms():
            te = atom_dict[a.name.strip()].type_energy
            extra[a.i_seq].vdwRadius = ener_lib.lib_atom[te].vdw_radius
            hb_type = ener_lib.lib_atom[te].hb_type
            if hb_type == "A":
              extra[a.i_seq].isAcceptor = True
            if hb_type == "D":
              extra[a.i_seq].isDonor = True

  # @todo

  #========================================================================
  # Call the test functions for all of our files.

  print('Testing DotSphere objects')
  ret = probe.DotSpheres_test()
  assert (len(ret) == 0)

  print('Testing SpatialQuery objects')
  ret = probe.SpatialQuery_test()
  assert (len(ret) == 0)

  print('Testing Scoring objects')
  ret = probe.Scoring_test()
  assert (len(ret) == 0)

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

  assert (len(ret) == 0)
