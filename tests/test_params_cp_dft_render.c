#include "cpmdc_params.h"
#include "cpmdc_features.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned char *slurp(const char *p, size_t *n) {
  FILE *f = fopen(p, "rb");
  if (!f) return NULL;
  fseek(f, 0, SEEK_END);
  long L = ftell(f);
  rewind(f);
  unsigned char *b = malloc((size_t)L);
  *n = fread(b, 1, (size_t)L, f);
  fclose(f);
  return b;
}

static int render_deck_file(const char *bin, char *deck, size_t deck_size) {
  size_t sz = 0;
  unsigned char *msg = slurp(bin, &sz);
  if (!msg) {
    fprintf(stderr, "missing %s\n", bin);
    return -1;
  }
  struct capn arena;
  CPMDParams_ptr root;
  int rc = cpmdc_params_root(msg, sz, &arena, &root);
  if (rc == 0) {
    rc = cpmdc_params_render_input_deck(root, deck, deck_size);
    cpmdc_params_release(&arena);
  }
  free(msg);
  return rc;
}

static int check_deck(const char *bin, const char **need, int nneed) {
  char deck[CPMDC_BLOCKS];
  if (render_deck_file(bin, deck, sizeof(deck)) != 0)
    return -1;
  for (int i = 0; i < nneed; ++i) {
    if (!strstr(deck, need[i])) {
      fprintf(stderr, "missing token [%s] in deck from %s\n", need[i], bin);
      fprintf(stderr, "--- deck ---\n%s\n", deck);
      return -1;
    }
  }
  return 0;
}

struct CatalogCoverage {
  const char *feature_id;
  const char *field_name;
  const char *render_token;
};

static int starts_with(const char *s, const char *prefix) {
  return strncmp(s, prefix, strlen(prefix)) == 0;
}

static const struct CatalogCoverage cpmd_coverage[] = {
    {"catalog.cpmd.OPTIMIZE_WAVEFUNCTION", "optimizeWavefunction",
     " OPTIMIZE WAVEFUNCTION\n"},
    {"catalog.cpmd.OPTIMIZE_GEOMETRY", "optimizeGeometry",
     " OPTIMIZE GEOMETRY\n"},
    {"catalog.cpmd.MOLECULAR_DYNAMICS", "molecularDynamics",
     " MOLECULAR DYNAMICS\n"},
    {"catalog.cpmd.MOLECULAR_DYNAMICS_CP", "molecularDynamicsCp",
     " MOLECULAR DYNAMICS CP\n"},
    {"catalog.cpmd.MOLECULAR_DYNAMICS_BO", "molecularDynamicsBo",
     " MOLECULAR DYNAMICS BO\n"},
    {"catalog.cpmd.MOLECULAR_DYNAMICS_EH", "molecularDynamicsEh",
     " MOLECULAR DYNAMICS EH\n"},
    {"catalog.cpmd.MOLECULAR_DYNAMICS_PT", "molecularDynamicsPt",
     " MOLECULAR DYNAMICS PT\n"},
    {"catalog.cpmd.MOLECULAR_DYNAMICS_CLASSICAL",
     "molecularDynamicsClassical", " MOLECULAR DYNAMICS CLASSICAL\n"},
    {"catalog.cpmd.MOLECULAR_DYNAMICS_FILE", "molecularDynamicsFile",
     " MOLECULAR DYNAMICS FILE\n"},
    {"catalog.cpmd.EMASS", "electronMass", " EMASS\n"},
    {"catalog.cpmd.TIMESTEP", "timestep", " TIMESTEP\n"},
    {"catalog.cpmd.MAXSTEP", "maxStep", " MAXSTEP\n"},
    {"catalog.cpmd.MAXITER", "maxIter", " MAXITER\n"},
    {"catalog.cpmd.CONVERGENCE_ORBITALS", "convergenceOrbitals",
     " CONVERGENCE ORBITALS\n"},
    {"catalog.cpmd.CONVERGENCE_GEOMETRY", "convergenceGeometry",
     " CONVERGENCE GEOMETRY\n"},
    {"catalog.cpmd.NOSE", "nose", " NOSE\n"},
    {"catalog.cpmd.NOSE_IONS", "noseIons", " NOSE IONS\n"},
    {"catalog.cpmd.NOSE_ELECTRONS", "noseElectrons",
     " NOSE ELECTRONS\n"},
    {"catalog.cpmd.BERENDSEN", "berendsen", " BERENDSEN\n"},
    {"catalog.cpmd.LANGEVIN", "langevin", " LANGEVIN\n"},
    {"catalog.cpmd.ANNEALING", "annealing", " ANNEALING\n"},
    {"catalog.cpmd.QUENCH", "quench", " QUENCH\n"},
    {"catalog.cpmd.RATTLE", "rattle", " RATTLE\n"},
    {"catalog.cpmd.SHAKE", "shake", " SHAKE\n"},
    {"catalog.cpmd.CONSTRAINT", "constraint", " CONSTRAINT\n"},
    {"catalog.cpmd.TROTTER", "trotter", " TROTTER\n"},
    {"catalog.cpmd.RESTART", "restart", " RESTART\n"},
    {"catalog.cpmd.RESTART_WAVEFUNCTION", "restartWavefunction",
     " RESTART WAVEFUNCTION\n"},
    {"catalog.cpmd.TRAJECTORY", "trajectory", " TRAJECTORY\n"},
    {"catalog.cpmd.PRINT", "printOptions", " PRINT\n"},
    {"catalog.cpmd.STORE", "storeOptions", " STORE\n"},
    {"catalog.cpmd.CENTER_MOLECULE_OFF", "centerMoleculeOff",
     " CENTER MOLECULE OFF\n"},
    {"catalog.cpmd.CENTER_MOLECULE_ON", "centerMoleculeOn",
     " CENTER MOLECULE ON\n"},
    {"catalog.cpmd.DIIS", "diis", " DIIS\n"},
    {"catalog.cpmd.ODIIS", "odiis", " ODIIS\n"},
    {"catalog.cpmd.PCG", "pcg", " PCG\n"},
    {"catalog.cpmd.DIAGONALIZATION", "diagonalization",
     " DIAGONALIZATION\n"},
    {"catalog.cpmd.FREE_ENERGY", "freeEnergy", " FREE-ENERGY\n"},
    {"catalog.cpmd.INTERFACE", "interface", " INTERFACE\n"},
    {"catalog.cpmd.QMMM", "qmmm", " QMMM\n"},
    {"catalog.cpmd.BICANONICAL_ENSEMBLE", "bicanonicalEnsemble",
     " BICANONICAL ENSEMBLE\n"},
    {"catalog.cpmd.CDFT", "cdft", " CDFT\n"},
    {"catalog.cpmd.PROPERTIES", "properties", " PROPERTIES\n"},
};

static const struct CatalogCoverage dft_scalar_coverage[] = {
    {"catalog.dft.LSD", "lsd", " LSD\n"},
    {"catalog.dft.GC_CUTOFF", "gcCutoff", " GC-CUTOFF\n"},
    {"catalog.dft.OLDCODE", "oldCode", " OLDCODE\n"},
    {"catalog.dft.NEWCODE", "newCode", " NEWCODE\n"},
    {"catalog.dft.XC_DRIVER", "xcDriver", " XC_DRIVER\n"},
    {"catalog.dft.LIBXC", "libxc", " LIBXC\n"},
    {"catalog.dft.LR_KERNEL", "lrKernel", " LR KERNEL\n"},
    {"catalog.dft.REFUNCT", "refunct", " REFUNCT\n"},
    {"catalog.dft.MTS_HIGH_FUNC", "mtsHighFunc", " MTS_HIGH_FUNC\n"},
    {"catalog.dft.MTS_LOW_FUNC", "mtsLowFunc", " MTS_LOW_FUNC\n"},
    {"catalog.dft.HFX", "hfx", " HFX\n"},
    {"catalog.dft.HFX_SCREENING", "hfxScreening", " HFX-SCREENING\n"},
    {"catalog.dft.HUBBARD", "hubbard", " HUBBARD\n"},
    {"catalog.dft.CORRELATION", "correlation", " CORRELATION\n"},
    {"catalog.dft.EXCHANGE", "exchange", " EXCHANGE\n"},
    {"catalog.dft.BECKE88", "becke88", " BECKE88\n"},
    {"catalog.dft.ALPHA", "alpha", " ALPHA\n"},
    {"catalog.dft.BETA", "beta", " BETA\n"},
};

static const struct CatalogCoverage *find_coverage(
    const struct CatalogCoverage *coverage, size_t ncoverage,
    const char *feature_id) {
  for (size_t i = 0; i < ncoverage; ++i) {
    if (strcmp(coverage[i].feature_id, feature_id) == 0)
      return &coverage[i];
  }
  return NULL;
}

static int decks_have_token(char **decks, int ndecks, const char *token) {
  for (int i = 0; i < ndecks; ++i) {
    if (strstr(decks[i], token))
      return 1;
  }
  return 0;
}

static int coverage_entry_is_valid(const struct CatalogCoverage *entry,
                                   char **decks, int ndecks) {
  const CPMDCFeatureEntry *feature = cpmdc_feature_find(entry->feature_id);
  if (!feature) {
    fprintf(stderr, "missing catalog feature [%s]\n", entry->feature_id);
    return 0;
  }
  if (feature->kind != CPMDC_FEATURE_KEYWORD) {
    fprintf(stderr, "catalog feature [%s] is not keyword kind\n",
            entry->feature_id);
    return 0;
  }
  if (!decks_have_token(decks, ndecks, entry->render_token)) {
    fprintf(stderr, "catalog feature [%s] mapped to field [%s] lacks rendered token [%s]\n",
            entry->feature_id, entry->field_name, entry->render_token);
    return 0;
  }
  return 1;
}

static int catalog_keyword_is_covered(const char *feature_id) {
  if (find_coverage(cpmd_coverage,
                    sizeof(cpmd_coverage) / sizeof(cpmd_coverage[0]),
                    feature_id))
    return 1;
  if (find_coverage(dft_scalar_coverage,
                    sizeof(dft_scalar_coverage) /
                        sizeof(dft_scalar_coverage[0]),
                    feature_id))
    return 1;
  if (strcmp(feature_id, "catalog.dft.FUNCTIONAL") == 0 ||
      starts_with(feature_id, "catalog.dft.FUNCTIONAL_"))
    return 1;
  return 0;
}

static int check_catalog_coverage(char **decks, int ndecks) {
  for (size_t i = 0; i < sizeof(cpmd_coverage) / sizeof(cpmd_coverage[0]);
       ++i) {
    if (!coverage_entry_is_valid(&cpmd_coverage[i], decks, ndecks))
      return -1;
  }
  for (size_t i = 0;
       i < sizeof(dft_scalar_coverage) / sizeof(dft_scalar_coverage[0]); ++i) {
    if (!coverage_entry_is_valid(&dft_scalar_coverage[i], decks, ndecks))
      return -1;
  }
  if (!decks_have_token(decks, ndecks, " FUNCTIONAL ")) {
    fprintf(stderr, "catalog.dft.FUNCTIONAL lacks rendered token\n");
    return -1;
  }
  const CPMDCFeatureEntry *features = cpmdc_feature_table();
  size_t nfeatures = cpmdc_feature_count();
  for (size_t i = 0; i < nfeatures; ++i) {
    const char *id = features[i].feature_id;
    if (features[i].kind != CPMDC_FEATURE_KEYWORD)
      continue;
    if (!starts_with(id, "catalog.cpmd.") && !starts_with(id, "catalog.dft."))
      continue;
    if (!catalog_keyword_is_covered(id)) {
      fprintf(stderr, "catalog keyword [%s] has no typed coverage mapping\n",
              id);
      return -1;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 8) {
    fprintf(stderr,
            "usage: %s cp_md.bin dft_func.bin cpmd_geometry.bin dft_scalars.bin cpmd_dynamics.bin cpmd_misc.bin long_tail_sections.bin\n",
            argv[0]);
    return 2;
  }
  /* inventory must resolve CP and DFT catalog ids */
  if (!cpmdc_feature_find("catalog.section.PIMD") ||
      !cpmdc_feature_find("catalog.section.VDW") ||
      !cpmdc_feature_find("catalog.cpmd.EMASS") ||
      !cpmdc_feature_find("catalog.dft.FUNCTIONAL_BLYP") ||
      !cpmdc_feature_find("catalog.dft.FUNCTIONAL_PBE0")) {
    fprintf(stderr, "inventory find failed for CP/DFT catalog ids\n");
    return 1;
  }
  const char *cp_need[] = {
      "&CPMD", "MOLECULAR DYNAMICS", "EMASS", "400.0", "NOSE IONS",
      "ANNEALING IONS", "TIMESTEP", "RESTART WAVEFUNCTION", "TRAJECTORY",
      "&SYSTEM", "SCALE", "0.5", "CHARGE", "1", "&DFT", "FUNCTIONAL BLYP",
  };
  const char *dft_need[] = {
      "&DFT", "FUNCTIONAL PBE0", "LSD", "GC-CUTOFF", "&PIMD", "TROTTER",
      "&VDW", "VDW CORRECTION",
  };
  const char *geometry_need[] = {
      "&CPMD", "OPTIMIZE GEOMETRY", "CONVERGENCE ORBITALS", "2e-06",
      "CONVERGENCE GEOMETRY", "0.0001", "MAXSTEP", "50", "MAXITER", "12",
      "TIMESTEP", "3", "EMASS", "450", "&SYSTEM", "CUTOFF", "80", "&DFT",
      "FUNCTIONAL PBE",
  };
  const char *dft_scalars_need[] = {
      "&DFT", "FUNCTIONAL PBE0", "LSD", "GC-CUTOFF", "1e-08",
      "XC_DRIVER", "LIBXC", "GGA_X_PBE GGA_C_PBE", "LR KERNEL", "REFUNCT",
      "MTS_HIGH_FUNC", "PBE0", "MTS_LOW_FUNC", "PBE", "HFX",
      "HFX-SCREENING", "0.2", "HUBBARD", "U 1 4.0", "ALPHA", "0.25",
      "BETA", "0.75", "OLDCODE", "NEWCODE", "CORRELATION", "LYP",
      "EXCHANGE", "B88", "BECKE88",
  };
  const char *cpmd_dynamics_need[] = {
      "&CPMD", "MOLECULAR DYNAMICS", "MOLECULAR DYNAMICS CP",
      "MOLECULAR DYNAMICS BO", "MOLECULAR DYNAMICS EH",
      "MOLECULAR DYNAMICS PT", "MOLECULAR DYNAMICS CLASSICAL",
      "MOLECULAR DYNAMICS FILE", "TRAJECTORY.in", "NOSE", "NOSE IONS",
      "NOSE ELECTRONS", "BERENDSEN", "300 100", "LANGEVIN", "ANNEALING",
      "IONS 300 50", "QUENCH", "RATTLE", "SHAKE", "CONSTRAINT", "FIX COM",
      "TROTTER", "8", "RESTART",
  };
  const char *cpmd_misc_need[] = {
      "&CPMD", "PRINT", "FORCES ON", "STORE", "WAVEFUNCTION",
      "CENTER MOLECULE OFF", "CENTER MOLECULE ON", "DIIS", "ODIIS", "PCG",
      "DIAGONALIZATION", "FREE-ENERGY", "INTERFACE", "QMMM",
      "BICANONICAL ENSEMBLE", "CDFT", "PROPERTIES",
  };
  const char *long_tail_sections_need[] = {
      "&ATOM", "ATOM SYMBOL", "H", "&BASIS", "BASIS SET", "DZVP",
      "&CLAS", "CLASSICAL FORCEFIELD", "demo.ff", "&EAM", "EAM POTENTIAL",
      "Cu.eam", "&EXTE", "EXTERNAL FIELD", "1.0", "&HARDNESS",
      "CHEMICAL HARDNESS", "1.25", "&INFO", "PRINT LEVEL", "LOW",
      "&LINRES", "DAVIDSON PARAMETER", "12", "&MOLSTATES", "STATE FILE",
      "states.dat", "&MTS", "MTS FACTOR", "4", "&NLCC", "NLCC ON",
      "&PATH", "REACTION COORDINATE", "path.dat", "&PIMD",
      "TROTTER FACTOR", "16", "&POTENTIAL", "POTENTIAL FILE", "bias.dat",
      "&PROP", "DIPOLE MOMENT", "&PTDDFT", "PROP_TSTEP", "0.25", "&RESP",
      "RESP CHARGES", "ON", "&TDDFT", "LR-TDDFT", "SINGLET", "&VDW",
      "VDW PARAMETERS", "GRIMME", "&VECTORS", "VECTOR FILE", "wf.vec",
      "&WAVEFUNCTION", "WAVEFUNCTION FILE", "RESTART.1",
  };
  if (check_deck(argv[1], cp_need, (int)(sizeof(cp_need) / sizeof(cp_need[0]))) != 0)
    return 1;
  if (check_deck(argv[2], dft_need, (int)(sizeof(dft_need) / sizeof(dft_need[0]))) != 0)
    return 1;
  if (check_deck(argv[3], geometry_need,
                 (int)(sizeof(geometry_need) / sizeof(geometry_need[0]))) != 0)
    return 1;
  if (check_deck(argv[4], dft_scalars_need,
                 (int)(sizeof(dft_scalars_need) /
                       sizeof(dft_scalars_need[0]))) != 0)
    return 1;
  if (check_deck(argv[5], cpmd_dynamics_need,
                 (int)(sizeof(cpmd_dynamics_need) /
                       sizeof(cpmd_dynamics_need[0]))) != 0)
    return 1;
  if (check_deck(argv[6], cpmd_misc_need,
                 (int)(sizeof(cpmd_misc_need) /
                       sizeof(cpmd_misc_need[0]))) != 0)
    return 1;
  if (check_deck(argv[7], long_tail_sections_need,
                 (int)(sizeof(long_tail_sections_need) /
                       sizeof(long_tail_sections_need[0]))) != 0)
    return 1;
  char *decks[7] = {0};
  for (int i = 0; i < 7; ++i) {
    decks[i] = malloc(CPMDC_BLOCKS);
    if (!decks[i])
      return 1;
    if (render_deck_file(argv[i + 1], decks[i], CPMDC_BLOCKS) != 0) {
      free(decks[i]);
      return 1;
    }
  }
  int coverage_rc = check_catalog_coverage(decks, 7);
  for (int i = 0; i < 7; ++i)
    free(decks[i]);
  if (coverage_rc != 0)
    return 1;
  printf("OK cp_md, dft_multi, cpmd_geometry, dft_scalars, cpmd_dynamics, cpmd_misc, and long_tail_sections render + inventory finds\n");
  return 0;
}
