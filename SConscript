import libtbx.load_env

Import("env_base", "env_etc")

env = env_base.Clone(
  LIBS=env_etc.libm)
if (env_etc.compiler != "win32_cl"):
  env.Replace(LINK=env_base["CC"])
  env.Append(CXXFLAGS=[
    "-std=c++11"]
  )
env.StaticLibrary(
  target=["#probe/lib/probelib"],
  source=[
    "Scoring.cpp",
    "DotSpheres.cpp",
    "SpatialQuery.cpp"])

if (not env_etc.no_boost_python):
  Import("env_cctbx_boost_python_ext")
  env_bpl = env_cctbx_boost_python_ext.Clone()
  env_bpl.Append(
    SHCXXFLAGS=[
      "-std=c++11"])
  if (env_etc.compiler != "win32_cl"):
    env.Append(SHCXXFLAGS=[
      "-std=c++11"]
    )
  env_bpl.Append(LIBPATH=["#probe/lib"])
  env_bpl.Prepend(LIBS=["probelib"])
  env_bpl.SharedLibrary(
    target="#lib/mmtbx_probe_ext",
    source=["probe_bpl.cpp"
            ])