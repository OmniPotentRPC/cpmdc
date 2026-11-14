/* LAMMPS plugin entry point that registers pair_style cpmdc.
 *
 *   plugin load cpmdcplugin.so
 *
 * The plugin is built against the headers of the LAMMPS it is loaded
 * into and links libcpmdc through its public C header only.
 */
#if __has_include("lammpsplugin.h")
#include "lammpsplugin.h"
#else
/* The C plugin interface of LAMMPS' src/PLUGIN/lammpsplugin.h, for header
 * sets that do not install it. */
extern "C" {
typedef void *(lammpsplugin_factory1)(void *);
typedef void *(lammpsplugin_factory2)(void *, int, char **);
typedef struct {
  const char *version;
  const char *style;
  const char *name;
  const char *info;
  const char *author;
  union {
    lammpsplugin_factory1 *v1;
    lammpsplugin_factory2 *v2;
  } creator;
  void *handle;
} lammpsplugin_t;
typedef void (*lammpsplugin_regfunc)(lammpsplugin_t *, void *);
void lammpsplugin_init(void *, void *, void *);
}
#endif

#include "lammps.h"
#if __has_include("version.h")
#include "version.h"
#endif
#include "pair_cpmdc.h"

using namespace LAMMPS_NS;

static Pair *cpmdc_creator(LAMMPS *lmp)
{
  return new PairCPMDC(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  auto register_plugin = (lammpsplugin_regfunc) regfunc;

#ifdef LAMMPS_VERSION
  plugin.version = LAMMPS_VERSION;
#else
  // Header sets without version.h: the loader compares against itself.
  plugin.version = static_cast<LAMMPS *>(lmp)->version;
#endif
  plugin.style = "pair";
  plugin.name = "cpmdc";
  plugin.info = "cpmdc pair style: energy, forces, and virial from libcpmdc";
  plugin.author = "Rohit Goswami (rgoswami@ieee.org)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &cpmdc_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
