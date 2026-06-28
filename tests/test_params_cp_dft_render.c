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

static int check_render_fails(const char *bin) {
  char deck[CPMDC_BLOCKS];
  if (render_deck_file(bin, deck, sizeof(deck)) == 0) {
    fprintf(stderr, "expected render failure from %s\n", bin);
    fprintf(stderr, "--- deck ---\n%s\n", deck);
    return -1;
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
     " OPTIMIZE GEOMETRY"},
    {"catalog.cpmd.OPTIMIZE_GEOMETRY_OPTIONS", "optimizeGeometryOptions",
     " OPTIMIZE GEOMETRY "},
    {"catalog.cpmd.OPTIMIZE_GEOMETRY_SAMPLE", "optimizeGeometrySample", "3"},
    {"catalog.cpmd.OPTIMIZE_COMBINED", "optimizeCombinedOptions",
     " OPTIMIZE COMBINED "},
    {"catalog.cpmd.OPTIMIZE_COMBINED_SAMPLE", "optimizeCombinedSample", "4"},
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
    {"catalog.cpmd.MOLECULAR_DYNAMICS_FILE_OPTIONS",
     "molecularDynamicsFileOptions", " MOLECULAR DYNAMICS FILE "},
    {"catalog.cpmd.MOLECULAR_DYNAMICS_BD", "molecularDynamicsBdTrajectories",
     " MOLECULAR DYNAMICS BD\n"},
    {"catalog.cpmd.PARRINELLO_RAHMAN", "parrinelloRahmanOptions",
     " PARRINELLO-RAHMAN "},
    {"catalog.cpmd.CHEBY", "cheby", " CHEBY\n"},
    {"catalog.cpmd.CAYLEY", "cayley", " CAYLEY\n"},
    {"catalog.cpmd.RUNGE_KUTTA", "rungeKutta", " RUNGE-KUTTA\n"},
    {"catalog.cpmd.FORCEMATCH", "forceMatch", " FORCEMATCH\n"},
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
    {"catalog.cpmd.GLOCALIZATION_PARAMETERS", "glocalizationParameters",
     " GLOCALIZATION PARAMETERS\n"},
    {"catalog.cpmd.GLOCALIZATION_OPTIMIZATION", "glocalizationOptimization",
     " GLOCALIZATION OPTIMIZATION "},
    {"catalog.cpmd.GFUNCTIONAL_TYPE", "gfunctionalType",
     " GFUNCTIONAL TYPE "},
    {"catalog.cpmd.SPREAD_RSPACE", "spreadRspace", " SPREAD RSPACE="},
    {"catalog.cpmd.G_UNITARITY_OPTIONS", "gUnitarityOptions", " PIPPO "},
    {"catalog.cpmd.STEP_TUNING", "stepTuning", " STEP TUNING\n"},
    {"catalog.cpmd.G_ANTISYM", "gAntisym", " G_ANTISYM"},
    {"catalog.cpmd.G_ANTISYM_PENALTY", "gAntisymPenalty",
     " G_ANTISYM PENALTY"},
    {"catalog.cpmd.G_KICK", "gKick", " G_KICK\n"},
    {"catalog.cpmd.G_COMPLEX", "gComplex", " G_COMPLEX\n"},
    {"catalog.cpmd.G_REAL", "gReal", " G_REAL\n"},
    {"catalog.cpmd.READ_MATRIX", "readMatrix", " READ MATRIX\n"},
    {"catalog.cpmd.G_STEP_TUNE", "gStepTune", " G_STEP TUNE\n"},
    {"catalog.cpmd.GLOC_WFNOUT", "glocWfnOutOptions", " GLOC WFNOUT "},
    {"catalog.cpmd.NO_GEO_CHECK", "noGeoCheck", " NO_GEO_CHECK\n"},
    {"catalog.cpmd.BROKEN", "brokenSymmetry", " BROKEN\n"},
    {"catalog.cpmd.DISTRIBUTED_LINALG", "distributedLinalg",
     " DISTRIBUTED LINALG "},
    {"catalog.cpmd.LINALG_NEWORTHO", "linalgNewOrtho",
     " LINALG NEWORTHO "},
    {"catalog.cpmd.DISORTHO_BSIZE", "disorthoBlockSize",
     " DISORTHO_BSIZE\n"},
    {"catalog.cpmd.BLOCKSIZE_STATES", "statesBlockSize",
     " BLOCKSIZE STATES\n"},
    {"catalog.cpmd.ALLTOALL", "allToAllPrecision", " ALLTOALL "},
    {"catalog.cpmd.GSHELL", "gshell", " GSHELL\n"},
    {"catalog.cpmd.LOCAL_POTENTIAL", "localPotential",
     " LOCAL POTENTIAL\n"},
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
    {"catalog.cpmd.CDFT", "cdft", " CDFT"},
    {"catalog.cpmd.CDFT_OPTIONS", "cdftOptions", " CDFT "},
    {"catalog.cpmd.CDFT_PAYLOAD", "cdftPayload", "1.0 2.0 0.3 0.4 50"},
    {"catalog.cpmd.CDFT_HDA_PAYLOAD", "cdftHdaPayload", "3.0 4.0"},
    {"catalog.cpmd.VGFACTOR", "vgfactor", " VGFACTOR\n"},
    {"catalog.cpmd.VMIRROR", "vMirror", " VMIRROR\n"},
    {"catalog.cpmd.COMBINE_SYSTEMS", "combineSystemsOptions",
     " COMBINE SYSTEMS "},
    {"catalog.cpmd.COMBINE_SYSTEMS_PAYLOAD", "combineSystemsPayload",
     "1 2 3 4 4 0.25"},
    {"catalog.cpmd.COMBINE_SYSTEMS_SAB_PAYLOAD",
     "combineSystemsSabPayload", "5 6"},
    {"catalog.cpmd.KSHAM", "kshamOptions", " KSHAM "},
    {"catalog.cpmd.KSHAM_PAYLOAD", "kshamPayload", "7 8"},
    {"catalog.cpmd.CZONES", "czonesSet", " CZONES SET\n"},
    {"catalog.cpmd.WOUT", "woutOptions", " WOUT "},
    {"catalog.cpmd.WOUT_PAYLOAD", "woutPayload", "0.25 10"},
    {"catalog.cpmd.XFMQC", "xfmqcTrajectories", " XFMQC\n"},
    {"catalog.cpmd.DEBUG", "debugOptions", " DEBUG "},
    {"catalog.cpmd.KOHN_SHAM_ENERGIES", "kohnShamEnergiesOptions",
     " KOHN-SHAM ENERGIES "},
    {"catalog.cpmd.KOHN_SHAM_ENERGIES_COUNT", "kohnShamEnergiesCount", "6"},
    {"catalog.cpmd.SURFACE_HOPPING", "surfaceHoppingOptions",
     " SURFACE HOPPING "},
    {"catalog.cpmd.ROKS", "roksOptions", " ROKS "},
    {"catalog.cpmd.ROKS_EXPERT_PAYLOAD", "roksExpertPayload",
     "1 2 3 4 5 6 7 8 9 10 11 12"},
    {"catalog.cpmd.PATH_SAMPLING", "pathSampling", " PATH SAMPLING\n"},
    {"catalog.cpmd.FIXRHO_UPWFN", "fixrhoUpwfnOptions", " FIXRHO UPWFN "},
    {"catalog.cpmd.FIXRHO_VECT", "fixrhoVectors", " FIXRHO VECT\n"},
    {"catalog.cpmd.FIXRHO_LOOP", "fixrhoLoop", " FIXRHO LOOP\n"},
    {"catalog.cpmd.FIXRHO_WFTOL", "fixrhoWftol", " FIXRHO WFTOL\n"},
    {"catalog.cpmd.BOGOLIUBOV_CORRECTION", "bogoliubovCorrection",
     " BOGOLIUBOV CORRECTION "},
    {"catalog.cpmd.VIBRATIONAL_ANALYSIS", "vibrationalAnalysisOptions",
     " VIBRATIONAL ANALYSIS "},
    {"catalog.cpmd.VIBRATIONAL_ANALYSIS_SAMPLE",
     "vibrationalAnalysisSample", "9"},
    {"catalog.cpmd.VIBRATIONAL_ANALYSIS_MODE", "vibrationalAnalysisMode",
     "MODE=3"},
    {"catalog.cpmd.ELECTRONIC_SPECTRA", "electronicSpectra",
     " ELECTRONIC SPECTRA\n"},
    {"catalog.cpmd.SPIN_ORBIT_COUPLING", "spinOrbitCouplingStates",
     " SPIN-ORBIT COUPLING\n"},
    {"catalog.cpmd.PROPAGATION_SPECTRA", "propagationSpectra",
     " PROPAGATION SPECTRA\n"},
    {"catalog.cpmd.PROPAGATION_DISTRUB", "propagationDistrub",
     " PROPAGATION DISTRUB\n"},
    {"catalog.cpmd.GAUGEPULSE", "gaugePulse", " GAUGEPULSE\n"},
    {"catalog.cpmd.GAUGEFIELD", "gaugeFieldFrequency", " GAUGEFIELD\n"},
    {"catalog.cpmd.NACV", "nacv", " NACV\n"},
    {"catalog.cpmd.ORBITAL_HARDNESS", "orbitalHardnessOptions",
     " ORBITAL HARDNESS "},
    {"catalog.cpmd.PATH_INTEGRAL", "pathIntegral", " PATH INTEGRAL\n"},
    {"catalog.cpmd.PATH_MINIMIZATION", "pathMinimization",
     " PATH MINIMIZATION\n"},
    {"catalog.cpmd.LANGEVIN_OPTIONS", "langevinOptions", " LANGEVIN "},
    {"catalog.cpmd.LANGEVIN_PARAMETER", "langevinParameter", "0.02"},
    {"catalog.cpmd.QMMMEASY", "qmmmEasy", " QMMMEASY\n"},
    {"catalog.cpmd.INTERFACE_OPTIONS", "interfaceOptions", " INTERFACE "},
    {"catalog.cpmd.TROTTER_FACTOR", "trotterFactorCount",
     " TROTTER FACTOR=2\n"},
    {"catalog.cpmd.TROTTER_FACTOR_PAYLOAD", "trotterFactorPayload",
     "0.1 0.2"},
    {"catalog.cpmd.LINEAR_RESPONSE", "linearResponse",
     " LINEAR RESPONSE\n"},
    {"catalog.cpmd.HARMONIC_REFERENCE", "harmonicReference",
     " HARMONIC REFERENCE "},
    {"catalog.cpmd.SCALED_MASSES", "scaledMasses", " SCALED MASSES "},
    {"catalog.cpmd.TDDFT", "tddft", " TDDFT\n"},
    {"catalog.cpmd.SSIC", "ssic", " SSIC\n"},
    {"catalog.cpmd.NONORTHOGONAL_ORBITALS", "nonorthogonalOrbitalsOptions",
     " NONORTHOGONAL ORBITALS "},
    {"catalog.cpmd.NONORTHOGONAL_ORBITALS_LIMIT",
     "nonorthogonalOrbitalsLimit", "1e-06"},
    {"catalog.cpmd.LANCZOS_DIAGONALIZATION", "lanczosDiagonalizationOptions",
     " LANCZOS DIAGONALIZATION "},
    {"catalog.cpmd.LANCZOS_PARAMETERS", "lanczosParametersCount",
     " LANCZOS PARAMETERS N=2\n"},
    {"catalog.cpmd.LANCZOS_PARAMETERS_PAYLOAD", "lanczosParametersPayload",
     "4 16 2 1e-08"},
    {"catalog.cpmd.DAVIDSON_DIAGONALIZATION", "davidsonDiagonalization",
     " DAVIDSON DIAGONALIZATION\n"},
    {"catalog.cpmd.DAVIDSON_PARAMETERS", "davidsonParameters",
     " DAVIDSON PARAMETERS\n"},
    {"catalog.cpmd.ALEXANDER_MIXING", "alexanderMixing",
     " ALEXANDER MIXING\n"},
    {"catalog.cpmd.ANDERSON_MIXING", "andersonMixingCount",
     " ANDERSON MIXING "},
    {"catalog.cpmd.ANDERSON_MIXING_PAYLOAD", "andersonMixingPayload",
     "0.2 0.1"},
    {"catalog.cpmd.BROYDEN_MIXING", "broydenMixingOptions",
     " BROYDEN MIXING "},
    {"catalog.cpmd.BROYDEN_MIXING_PAYLOAD", "broydenMixingPayload",
     "MIX=0.3"},
    {"catalog.cpmd.DIIS_MIXING", "diisMixingCount", " DIIS MIXING N=2\n"},
    {"catalog.cpmd.DIIS_MIXING_PAYLOAD", "diisMixingPayload", "0.1 4"},
    {"catalog.cpmd.MOVERHO", "moverhoMixing", " MOVERHO\n"},
    {"catalog.cpmd.EXTRAPOLATE_WFN", "extrapolateWfnOptions",
     " EXTRAPOLATE WFN "},
    {"catalog.cpmd.EXTRAPOLATE_WFN_ORDER", "extrapolateWfnOrder", "5"},
    {"catalog.cpmd.EXTRAPOLATE_CONSTRAINT", "extrapolateConstraintOrder",
     " EXTRAPOLATE CONSTRAINT\n"},
    {"catalog.cpmd.TSDE", "tsdeOptions", " TSDE "},
    {"catalog.cpmd.TSDP", "tsdpOptions", " TSDP "},
    {"catalog.cpmd.TCGP", "tcgp", " TCGP\n"},
    {"catalog.cpmd.TSDC", "tsdc", " TSDC\n"},
    {"catalog.cpmd.STEEPEST_DESCENT", "steepestDescentOptions",
     " STEEPEST DESCENT "},
    {"catalog.cpmd.CONJUGATE_GRADIENT", "conjugateGradientOptions",
     " CONJUGATE GRADIENT "},
    {"catalog.cpmd.ODIIS_OPTIONS", "odiisOptions", " ODIIS "},
    {"catalog.cpmd.ODIIS_VECTORS", "odiisVectors", "7"},
    {"catalog.cpmd.HAMILTONIAN_CUTOFF", "hamiltonianCutoff",
     " HAMILTONIAN CUTOFF\n"},
    {"catalog.cpmd.GDIIS", "gdiisVectors", " GDIIS\n"},
    {"catalog.cpmd.PROPERTIES", "properties", " PROPERTIES\n"},
    {"catalog.cpmd.VDW_CORRECTION", "vdwCorrection", " VDW CORRECTION"},
    {"catalog.cpmd.VDW_WANNIER", "vdwWannier", " VDW WANNIER"},
    {"catalog.cpmd.DCACP", "dcacp", " DCACP\n"},
    {"catalog.cpmd.ISOLATED_MOLECULE", "isolatedMolecule",
     " ISOLATED MOLECULE\n"},
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
  if (argc < 19) {
    fprintf(stderr,
            "usage: %s cp_md.bin dft_func.bin cpmd_geometry.bin dft_scalars.bin cpmd_dynamics.bin cpmd_misc.bin long_tail_sections.bin vdw_controls.bin system_controls.bin system_kpoints_monkhorst.bin system_kpoints_bands.bin system_cdft_acceptor_wmult.bin system_couplings_prod.bin system_couplings_linres.bin system_couplings_lists.bin system_cell_qualifiers.bin system_cell_vectors.bin bad_occupation.bin\n",
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
      "&SYSTEM", "SCALE S=0.5", "CHARGE", "1", "&DFT", "FUNCTIONAL BLYP",
  };
  const char *dft_need[] = {
      "&DFT", "FUNCTIONAL PBE0", "LSD", "GC-CUTOFF", "&PIMD", "TROTTER",
      "&VDW", "VDW CORRECTION",
  };
  const char *geometry_need[] = {
      "&CPMD", "OPTIMIZE GEOMETRY CLASSICAL XYZ SAMPLE", "3",
      "OPTIMIZE COMBINED XYZ SAMPLE", "4", "CONVERGENCE ORBITALS", "2e-06",
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
      "MOLECULAR DYNAMICS FILE", "TRAJECTORY.in",
      "MOLECULAR DYNAMICS FILE XYZ NSKIP=2 NSAMPLE=5",
      "MOLECULAR DYNAMICS BD", "4",
      "PARRINELLO-RAHMAN NPT SHOCK TOLKINC=0.01", "CHEBY", "CAYLEY",
      "RUNGE-KUTTA", "FORCEMATCH", "MAXRUNTIME", "3600",
      "TIMESTEP ELECTRONS", "2.5", "TIMESTEP IONS", "5", "CMASS", "1500",
      "CONVERGENCE CELL", "0.0002", "CONVERGENCE ADAPT", "0.0003",
      "CONVERGENCE ENERGY", "0.0004", "CONVERGENCE CALFOR", "0.0005",
      "CONVERGENCE RELAX", "12", "CONVERGENCE RHOFIX", "0.0006",
      "CONVERGENCE INITIAL", "0.0007", "CONVERGENCE CONSTRAINT",
      "0.01 0.02",
      "TEMPERATURE ELECTRON", "0.125", "TEMPERATURE RAMP", "300 400 25",
      "RESCALE OLD VELOCITIES", "REVERSE VELOCITIES", "SUBTRACT COMVEL",
      "20", "SUBTRACT ROTVEL", "40", "PRNGSEED", "12345",
      "TEMPCONTROL IONS", "300 0.01", "TEMPCONTROL ELECTRONS",
      "0.5 0.001", "TEMPCONTROL CELL", "1.5 0.02", "BERENDSEN IONS",
      "300 100", "BERENDSEN ELECTRONS", "0.5 50", "BERENDSEN CELL",
      "1.5 75", "NOSE", "NOSE IONS", "300 1200", "NOSE ELECTRONS",
      "0.5 900", "NOSE CELL", "1.5 800", "NOSE PARAMETERS",
      "3 4 5 6.0 7 8", "BERENDSEN", "LANGEVIN", "ANNEALING",
      "ANNEALING IONS", "0.98", "ANNEALING ELECTRONS", "0.97",
      "ANNEALING CELL", "0.96", "IONS 300 50", "DAMPING IONS", "0.11",
      "DAMPING ELECTRONS", "0.12", "DAMPING CELL", "0.13",
      "HESSIAN SCHLEGEL PARTIAL", "PROJECT FULL", "STRESS TENSOR VIRIAL",
      "15", "CLASSTRESS", "20", "QUENCH", "RATTLE", "SHAKE",
      "CONSTRAINT", "FIX COM", "TROTTER", "8", "RESTART",
  };
  const char *cpmd_misc_need[] = {
      "&CPMD", "PRINT", "FORCES ON", "STORE", "WAVEFUNCTION",
      "STORE WAVEFUNCTION DENSITY", "25 SC=5", "STORE OFF POTENTIAL",
      "RESTFILE SAMPLE", "3", "10 20 30",
      "TRAJECTORY XYZ BINARY FORCES RANGE SAMPLE", "1 4", "7",
      "MOVIE SAMPLE", "9", "MOVIE OFF",
      "ENERGYBANDS", "EXTERNAL POTENTIAL ADD",
      "ELECTROSTATIC POTENTIAL SAMPLE=4", "DIPOLE DYNAMICS SAMPLE WANNIER",
      "6", "RHOOUT BANDS SAMPLE=5", "2", "1 3", "ELF PARAMETER",
      "0.8 0.01",
      "WANNIER PARAMETER", "0.1 1e-06 0.5 200",
      "WANNIER OPTIMIZATION JACOBI", "WANNIER TYPE RESTA",
      "WANNIER REFERENCE", "0.1 0.2 0.3", "WANNIER SERIAL",
      "WANNIER DOS", "WANNIER MOLECULAR", "WANNIER WFNOUT PARTIAL DENSITY",
      "1 3",
      "COMPRESS WRITE32", "MEMORY CHECK SMALL", "REAL SPACE WFN KEEP SIZE",
      "2048", "SPLINE POINTS QFUNCTION INIT RANGE", "64", "FINITE DIFFERENCES",
      "0.01 COORD=1 2 3 RADIUS=4", "TASKGROUPS CARTESIAN", "2",
      "DISTRIBUTE FNL ON", "FILEPATH", "restart/path", "BENCHMARK",
      "1 2 3 4 5 6 7 8 9 10", "MIRROR", "SHIFT POTENTIAL", "0.25",
      "GLOCALIZATION PARAMETERS", "0.2 1e-05 0.01 150 -1",
      "GLOCALIZATION OPTIMIZATION SD INTEGRATION",
      "GFUNCTIONAL TYPE ZETANUMBER", "SPREAD RSPACE=0.7 0.2",
      "PIPPO LAGRANGE=1e-08 40", "STEP TUNING", "G_ANTISYM PENALTY",
      "G_KICK", "G_COMPLEX", "G_REAL", "READ MATRIX", "G_STEP TUNE",
      "GLOC WFNOUT REAL PARTIAL", "2 5",
      "NO_GEO_CHECK", "BROKEN", "DISTRIBUTED LINALG ON",
      "LINALG NEWORTHO OFF", "DISORTHO_BSIZE", "32", "BLOCKSIZE STATES",
      "64", "ALLTOALL SINGLE", "GSHELL", "LOCAL POTENTIAL",
      "CENTER MOLECULE OFF", "CENTER MOLECULE ON", "DIIS", "ODIIS", "PCG",
      "DIAGONALIZATION", "FREE-ENERGY", "INTERFACE", "QMMM",
      "BICANONICAL ENSEMBLE",
      "CDFT ALL NEWTON HDA PROJECT AUTO NOCCOR RESWF PHIOUT",
      "1.0 2.0 0.3 0.4 50", "3.0 4.0", "VGFACTOR", "0.2 0.1",
      "VMIRROR", "COMBINE SYSTEMS REFLECT SAB NONORTH",
      "1 2 3 4 4 0.25", "5 6", "KSHAM MATRIX ROUT", "7 8",
      "CZONES SET", "0.1 0.2 0.3", "0.4 0.5", "WOUT SLICE",
      "0.25 10", "XFMQC", "11", "PROPERTIES",
      "DEBUG FORCES FILE IO NOACC HUBBARD",
      "KOHN-SHAM ENERGIES NOWAVEFUNCTION", "6",
      "SURFACE HOPPING NEW RESTART SHIFT=0.25", "ROKS TRIPLET LOCALIZED",
      "1 2 3 4 5 6 7 8 9 10 11 12", "PATH SAMPLING",
      "FIXRHO UPWFN PCG", "FIXRHO VECT", "5", "FIXRHO LOOP", "2 9",
      "FIXRHO WFTOL", "1e-05", "BOGOLIUBOV CORRECTION OFF",
      "VIBRATIONAL ANALYSIS LR IN ACLIMAX GAUSS SAMPLE MODE=3", "9",
      "ELECTRONIC SPECTRA", "SPIN-ORBIT COUPLING", "2 4",
      "PROPAGATION SPECTRA", "PROPAGATION DISTRUB", "GAUGEPULSE",
      "GAUGEFIELD", "0.75", "NACV", "ORBITAL HARDNESS LR",
      "PATH INTEGRAL", "PATH MINIMIZATION",
      "LANGEVIN SMART MOVECM CENTROIDOFF", "0.02", "QMMMEASY",
      "INTERFACE GMX ESP PCGFIRST", "TROTTER FACTOR=2", "0.1 0.2",
      "0.3 0.4", "LINEAR RESPONSE", "HARMONIC REFERENCE OFF",
      "SCALED MASSES OFF", "TDDFT", "SSIC", "0.5",
      "NONORTHOGONAL ORBITALS ON", "1e-06",
      "LANCZOS DIAGONALIZATION ALL NEW PERIODIC=4",
      "LANCZOS PARAMETERS N=2", "4 16 2 1e-08", "0.5 2e-08",
      "DAVIDSON DIAGONALIZATION", "DAVIDSON PARAMETERS", "8 1e-06 4",
      "ALEXANDER MIXING", "0.15", "ANDERSON MIXING G-SPACE N=2",
      "0.2 0.1", "0.5 0.05", "BROYDEN MIXING ADAPT",
      "MIX=0.3 CUTBROY=5 W02=0.01 NFRBROY=4 RESET=8 KERMIX=0.7",
      "DIIS MIXING N=2", "0.1 4", "0.4 8", "MOVERHO", "0.25",
      "EXTRAPOLATE WFN STORE ASPC CSTEPS=3", "5",
      "EXTRAPOLATE CONSTRAINT", "2",
      "TSDE NOPREC", "TSDP LINE", "TCGP", "TSDC",
      "STEEPEST DESCENT ELECTRONS NOPREC IONS LINE CELL",
      "CONJUGATE GRADIENT ELECTRONS NOPREC MINIMIZE IONS",
      "ODIIS NOPREC NO_RESET=6", "7", "HAMILTONIAN CUTOFF", "0.001",
      "GDIIS", "8",
  };
  const char *vdw_controls_need[] = {
      "&CPMD",
      "VDW CORRECTION ON",
      "VDW WANNIER OFF",
      "DCACP",
      "&VDW",
      "EMPIRICAL CORRECTION",
      "VDW PARAMETERS",
      "ALL DFT-D2",
      "S6GRIMME",
      "PBE",
      "CUTOFF",
      "0.01",
      "CELL",
      "1 1 0",
      "DFT-D3 PARAMETERS THREEBODY NUMGRAD",
      "94.868329 40.0",
      "DFT-D3 FUNCTIONAL VERSION=4 TZ",
      "PBE0",
      "DFT-D3 NO_CONTRIBUTION",
      "1 2",
      "END EMPIRICAL CORRECTION",
      "CUSTOM VDW NOTE",
      "kept",
  };
  const char *system_controls_need[] = {
      "&SYSTEM", "SYMMETRY", "0", "ANGSTROM", "CELL", "14 1 1 0 0 0",
      "REFERENCE CELL", "15 1 1 0 0 0", "CLASSICAL CELL",
      "16 1.1 1.2 0 0 0", "ISOTROPIC CELL", "ZFLEXIBLE CELL", "CUTOFF",
      "CHECK SYMMETRY", "2e-05", "90", "NOSPHERICAL CUTOFF", "HFX CUTOFF",
      "35 140", "DENSITY CUTOFF", "360", "DENSITY CUTOFF NUMBER", "64",
      "DUAL", "5", "CONSTANT CUTOFF", "0.8 0.2 75", "BOX WALLS", "6.5",
      "CHARGE", "1", "NSUP", "2", "STATES", "4", "OCCUPATION FIXED",
      "2 1 1 0", "EXTERNAL FIELD", "0.01 0.02 0.03", "WCUT 0.45",
      "WGAUSS 2", "0.1234", "0.5678", "COUPLINGS FD=0.01",
      "DONOR 2 WMULT", "1 3", "2 4", "ACCEPTOR 2 HDAS", "2 4", "1 3",
      "POINT GROUP MOLECULE DELTA=0.0001", "DNH 2",
      "LOW SPIN EXCITATION ROKS PENALTY LSETS", "0.75",
      "LSE PARAMETERS", "2 1", "MODIFIED GOEDECKER PARAMETERS",
      "-0.25 0.25", "ENERGY PROFILE",
      "POISSON SOLVER HOCKNEY PARAMETER", "1.2", "MESH", "24 24 32",
      "KPOINTS SCALED ONLYDIAG BLOCK=2 ALL NOSWAP", "0 0 0 0.5",
      "0.5 0 0 0.5",
      "SCALE CARTESIAN SX=2 SY=3 SZ=4", "DOUBLE GRID ON", "PRESSURE", "12.5",
      "STRESS TENSOR", "1 0.1 0.2 0.1 2 0.3 0.2 0.3 3",
      "SYMMETRIZE COORDINATES", "TESR", "4", "ISOLATED MOLECULE",
      "CENTER MOLECULE ON", "CENTER MOLECULE OFF", "SURFACE XY", "POLYMER",
      "CLUSTER", "SHOCK VELOCITY", "4500",
  };
  const char *system_monkhorst_need[] = {
      "&SYSTEM", "KPOINTS MONKHORST-PACK SYMMETRIZED FULL KDP",
      "2 3 4 SHIFT=0.5 0.25 0",
  };
  const char *system_kpoint_bands_need[] = {
      "&SYSTEM", "KPOINTS BANDS SCALED", "4 0 0 0 0.5 0.5 0",
      "0 0 0 0 0 0 0",
  };
  const char *system_cdft_acceptor_wmult_need[] = {
      "&SYSTEM", "ACCEPTOR 2 WMULT", "2 4", "5 6",
  };
  const char *system_couplings_prod_need[] = {
      "&SYSTEM", "COUPLINGS PROD=0.02",
  };
  const char *system_couplings_linres_need[] = {
      "&SYSTEM", "COUPLINGS LINRES TOL=1e-05 NVECT=8 SPECIFY THRESHOLDS BRUTE FORCE",
      "0.1 0.2", "0.3 0.4", "0.5 0.6",
  };
  const char *system_couplings_lists_need[] = {
      "&SYSTEM", "COUPLINGS NSURF=2", "1 2 0.5", "2 3 0.75",
      "COUPLINGS NAT=3", "1 4 7",
  };
  const char *system_cell_qualifiers_need[] = {
      "&SYSTEM", "CELL ABSOLUTE DEGREE", "10 11 12 90 100 110",
      "REFERENCE CELL ABSOLUTE DEGREE", "13 14 15 80 85 95",
      "CLASSICAL CELL ABSOLUTE DEGREE", "16 17 18 70 75 85",
  };
  const char *system_cell_vectors_need[] = {
      "&SYSTEM", "CELL VECTORS", "1 0 0 0 2 0 0 0 3",
      "REFERENCE CELL VECTORS", "2 0 0 0 3 0 0 0 4",
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
  if (check_deck(argv[8], vdw_controls_need,
                 (int)(sizeof(vdw_controls_need) /
                       sizeof(vdw_controls_need[0]))) != 0)
    return 1;
  if (check_deck(argv[9], system_controls_need,
                 (int)(sizeof(system_controls_need) /
                       sizeof(system_controls_need[0]))) != 0)
    return 1;
  if (check_deck(argv[10], system_monkhorst_need,
                 (int)(sizeof(system_monkhorst_need) /
                       sizeof(system_monkhorst_need[0]))) != 0)
    return 1;
  if (check_deck(argv[11], system_kpoint_bands_need,
                 (int)(sizeof(system_kpoint_bands_need) /
                       sizeof(system_kpoint_bands_need[0]))) != 0)
    return 1;
  if (check_deck(argv[12], system_cdft_acceptor_wmult_need,
                 (int)(sizeof(system_cdft_acceptor_wmult_need) /
                       sizeof(system_cdft_acceptor_wmult_need[0]))) != 0)
    return 1;
  if (check_deck(argv[13], system_couplings_prod_need,
                 (int)(sizeof(system_couplings_prod_need) /
                       sizeof(system_couplings_prod_need[0]))) != 0)
    return 1;
  if (check_deck(argv[14], system_couplings_linres_need,
                 (int)(sizeof(system_couplings_linres_need) /
                       sizeof(system_couplings_linres_need[0]))) != 0)
    return 1;
  if (check_deck(argv[15], system_couplings_lists_need,
                 (int)(sizeof(system_couplings_lists_need) /
                       sizeof(system_couplings_lists_need[0]))) != 0)
    return 1;
  if (check_deck(argv[16], system_cell_qualifiers_need,
                 (int)(sizeof(system_cell_qualifiers_need) /
                       sizeof(system_cell_qualifiers_need[0]))) != 0)
    return 1;
  if (check_deck(argv[17], system_cell_vectors_need,
                 (int)(sizeof(system_cell_vectors_need) /
                       sizeof(system_cell_vectors_need[0]))) != 0)
    return 1;
  if (check_render_fails(argv[18]) != 0)
    return 1;
  char *decks[17] = {0};
  for (int i = 0; i < 17; ++i) {
    decks[i] = malloc(CPMDC_BLOCKS);
    if (!decks[i])
      return 1;
    if (render_deck_file(argv[i + 1], decks[i], CPMDC_BLOCKS) != 0) {
      free(decks[i]);
      return 1;
    }
  }
  int coverage_rc = check_catalog_coverage(decks, 11);
  for (int i = 0; i < 17; ++i)
    free(decks[i]);
  if (coverage_rc != 0)
    return 1;
  printf("OK cp_md, dft_multi, cpmd_geometry, dft_scalars, cpmd_dynamics, cpmd_misc, long_tail_sections, vdw_controls, system_controls, system_monkhorst, system_kpoint_bands, system_cdft_acceptor_wmult, system_couplings_prod, system_couplings_linres, system_couplings_lists, system_cell_qualifiers, and system_cell_vectors render + inventory finds\n");
  return 0;
}
