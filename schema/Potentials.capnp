# @brief RPC schema for distributed potential evaluations (cpmdc / rgpot).
#
# Shared ForceInput / PotentialResult carrier matches nwchemc and rgpot-core so
# the same Cap'n Proto messages work for direct in-process "socket" calls and
# for remote Potential.calculate RPCs.

@0xc8a4f2e19b3d7051;

# @struct ForceInput
# @brief Input configuration for a potential energy evaluation.
# @field lengthUnit Unit expression for positions/box (default "angstrom").
# @field energyUnit Unit expression for energy output (default "eV").
# Unit strings are parsed by cpmdc unit helpers (same names as rgpot).
struct ForceInput {
  pos        @0 :List(Float64); # Flat array of atomic coordinates [natoms * 3].
  atmnrs     @1 :List(Int32);   # Array of atomic numbers [natoms].
  box        @2 :List(Float64); # Simulation cell vectors [9] (row-major 3x3).
  lengthUnit @3 :Text = "angstrom"; # Unit for positions and box vectors.
  energyUnit @4 :Text = "eV";       # Unit for energy and forces output.
}

# @struct PotentialResult
# @brief Results returned from a potential energy evaluation.
struct PotentialResult {
  energy @0 :Float64;       # Calculated potential energy.
  forces @1 :List(Float64); # Flat array of atomic forces [natoms * 3].
}

# --- CPMD structured input (section-oriented, mirrors CPMD &SECTION decks) ---

struct CPMDDirective {
  keyword @0 :Text;       # Line keyword inside a section, e.g. "CONVERGENCE ORBITALS".
  args    @1 :List(Text); # Tokenized arguments on following lines or same line.
}

struct CPMDGenericSection {
  name       @0 :Text;                 # Section name without ampersands, e.g. "CPMD".
  directives @1 :List(CPMDDirective);  # Structured section body.
}

struct CPMDSetDirective {
  key   @0 :Text; # Dotted SECTION.KEYWORD, e.g. "CPMD.PRINT FORCES ON".
  value @1 :Text; # Optional value emitted on the following indented line.
}

struct CPMDKPoint {
  coordinates @0 :List(Float64); # Explicit KPOINTS coordinates kx, ky, kz.
  weight      @1 :Float64 = 1.0; # KPOINTS integration weight.
}

struct CPMDKPointBand {
  points @0 :Int32 = 0;      # Number of interpolated KPOINTS BANDS points.
  start  @1 :List(Float64);  # Initial KPOINTS BANDS vector.
  end    @2 :List(Float64);  # Final KPOINTS BANDS vector.
}

struct CPMDSystemSection {
  symmetry        @0 :Int32 = 0;          # SYMMETRY code.
  angstrom        @1 :Bool = true;        # Emit ANGSTROM for CELL / ATOMS.
  cell            @2 :List(Float64);      # CELL vectors or a,b,c,alpha,beta,gamma (6 or 9).
  cutOffRy        @3 :Float64 = 70.0;     # CUTOFF in Rydberg (CPMD default unit).
  scale           @4 :Float64 = 0.0;      # Optional SCALE factor; 0 => omit.
  charge          @5 :Int32 = 0;          # CHARGE of the system.
  multiplicity    @6 :Int32 = 1;          # Spin multiplicity 2S+1 (maps to LSD/spin when >1).
  directives      @7 :List(CPMDDirective);
  densityCutOffRy @8 :Float64 = 0.0;      # DENSITY CUTOFF; 0 => omit.
  poissonSolver   @9 :Text;               # POISSON SOLVER argument, e.g. HOCKNEY.
  poissonParameter @10 :Float64 = 0.0;    # Optional POISSON ... PARAMETER value.
  surface         @11 :Text;              # SURFACE direction, e.g. XY/YZ/ZX.
  referenceCell   @12 :List(Float64);      # REFERENCE CELL values (6) or vectors (9).
  classicalCell   @13 :List(Float64);      # CLASSICAL CELL values (6).
  isotropicCell   @14 :Bool = false;       # ISOTROPIC CELL.
  zFlexibleCell   @15 :Bool = false;       # ZFLEXIBLE CELL.
  densityCutoffNumber @16 :Int32 = 0;      # DENSITY CUTOFF NUMBER; 0 => omit.
  dual            @17 :Float64 = 0.0;      # DUAL density cutoff factor; 0 => omit.
  constantCutoff  @18 :List(Float64);      # CONSTANT CUTOFF akin, skin, eckin (3).
  mesh            @19 :List(Int32);        # MESH nr1, nr2, nr3.
  scaleCartesian  @20 :Bool = false;       # SCALE CARTESIAN; combines with scale as S=.
  doubleGrid      @21 :Text;               # DOUBLE GRID argument, e.g. ON/OFF.
  symmetrizeCoordinates @22 :Bool = false; # SYMMETRIZE COORDINATES.
  tesr            @23 :Int32 = 0;          # TESR value; 0 => omit.
  polymer         @24 :Bool = false;       # POLYMER.
  cluster         @25 :Bool = false;       # CLUSTER.
  cutoffShape     @26 :Text;               # CUTOFF qualifier, e.g. SPHERICAL/NOSPHERICAL.
  hfxCutoff       @27 :List(Float64);       # HFX CUTOFF hfxwfe, hfxdee (2).
  boxWalls        @28 :Float64 = 0.0;       # BOX WALLS skin; 0 => omit.
  nSup            @29 :Int32 = 0;           # NSUP alpha spin states; 0 => omit.
  states          @30 :Int32 = 0;           # STATES count; 0 => omit.
  occupation      @31 :List(Float64);       # OCCUPATION values; length should match states.
  occupationFixed @32 :Bool = false;        # OCCUPATION FIXED.
  externalField   @33 :List(Float64);       # EXTERNAL FIELD vector (3).
  pressure        @34 :Float64 = 0.0;       # PRESSURE; 0 => omit.
  stressTensor    @35 :List(Float64);       # STRESS TENSOR values (9).
  shockVelocity   @36 :Float64 = 0.0;       # SHOCK VELOCITY; 0 => omit.
  checkSymmetryPrecision @37 :Float64 = 0.0; # CHECK SYMMETRY precision; 0 => omit.
  checkSymmetryOff @38 :Bool = false;       # CHECK SYMMETRY OFF.
  wCut            @39 :Float64 = 0.0;       # WCUT; 0 => omit.
  wGauss          @40 :List(Float64);       # WGAUSS sigma values.
  lowSpinExcitation @41 :Text;              # LOW SPIN EXCITATION options, e.g. ROKS.
  lowSpinExcitationLsets @42 :Bool = false; # LOW SPIN EXCITATION LSETS.
  lseParameters   @43 :List(Float64);       # LSE PARAMETERS lsea, lseb (2).
  modifiedGoedecker @44 :Bool = false;      # MODIFIED GOEDECKER.
  modifiedGoedeckerParameters @45 :List(Float64); # MODIFIED GOEDECKER PARAMETERS lambda_ab, lambda_ba (2).
  energyProfile   @46 :Bool = false;        # ENERGY PROFILE.
  pointGroup       @47 :Text;                # POINT GROUP selector line, e.g. AUTO, NAME=C2v, or 26.
  pointGroupDelta  @48 :Float64 = 0.0;       # POINT GROUP DELTA accuracy; 0 => omit.
  pointGroupMolecule @49 :Bool = false;      # POINT GROUP MOLECULE.
  scaleX         @50 :Float64 = 0.0;          # SCALE SX= value; 0 => omit.
  scaleY         @51 :Float64 = 0.0;          # SCALE SY= value; 0 => omit.
  scaleZ         @52 :Float64 = 0.0;          # SCALE SZ= value; 0 => omit.
  kpoints        @53 :List(CPMDKPoint);       # KPOINTS explicit weighted points.
  kpointsScaled  @54 :Bool = false;           # KPOINTS SCALED.
  kpointsOnlyDiagonal @55 :Bool = false;      # KPOINTS ONLYDIAG.
  kpointsMonkhorstPack @56 :List(Int32);      # KPOINTS MONKHORST-PACK nk1, nk2, nk3.
  kpointsMonkhorstSymmetrized @57 :Bool = false; # KPOINTS MONKHORST-PACK SYMMETRIZED.
  kpointsMonkhorstFull @58 :Bool = false;     # KPOINTS MONKHORST-PACK FULL.
  kpointsMonkhorstKdp @59 :Bool = false;      # KPOINTS MONKHORST-PACK KDP.
  kpointsMonkhorstShift @60 :List(Float64);   # KPOINTS MONKHORST-PACK mesh-line SHIFT vector.
  kpointBands    @61 :List(CPMDKPointBand);   # KPOINTS BANDS segments.
  kpointsBlock   @62 :Int32 = 0;              # KPOINTS BLOCK=n; 0 => omit.
  kpointsBlockAll @63 :Bool = false;          # KPOINTS BLOCK ALL.
  kpointsBlockCalculated @64 :Bool = false;   # KPOINTS BLOCK CALCULATED.
  kpointsBlockNoSwap @65 :Bool = false;       # KPOINTS BLOCK NOSWAP.
}

struct CPMDCpmdSection {
  optimizeWavefunction @0 :Bool = true;   # OPTIMIZE WAVEFUNCTION.
  molecularDynamics    @1 :Bool = false;  # MOLECULAR DYNAMICS.
  convergenceOrbitals  @2 :Float64 = 1.0e-6;
  maxStep              @3 :Int32 = 0;     # MAXSTEP; 0 => omit (CPMD default).
  timestep             @4 :Float64 = 0.0; # TIMESTEP; 0 => omit.
  restartWavefunction  @5 :Bool = false;  # RESTART WAVEFUNCTION.
  trajectory           @6 :Bool = false;  # TRAJECTORY.
  directives           @7 :List(CPMDDirective);
  optimizeGeometry     @8 :Bool = false;  # OPTIMIZE GEOMETRY.
  maxIter              @9 :Int32 = 0;     # MAXITER; 0 => omit.
  convergenceGeometry  @10 :Float64 = 0.0; # CONVERGENCE GEOMETRY; 0 => omit.
  electronMass         @11 :Float64 = 0.0; # EMASS; 0 => omit.
  molecularDynamicsCp        @12 :Bool = false; # MOLECULAR DYNAMICS CP.
  molecularDynamicsBo        @13 :Bool = false; # MOLECULAR DYNAMICS BO.
  molecularDynamicsEh        @14 :Bool = false; # MOLECULAR DYNAMICS EH.
  molecularDynamicsPt        @15 :Bool = false; # MOLECULAR DYNAMICS PT.
  molecularDynamicsClassical @16 :Bool = false; # MOLECULAR DYNAMICS CLASSICAL.
  molecularDynamicsFile      @17 :Text;         # MOLECULAR DYNAMICS FILE.
  nose                 @18 :Bool = false; # NOSE.
  noseIons             @19 :Bool = false; # NOSE IONS.
  noseElectrons        @20 :Bool = false; # NOSE ELECTRONS.
  berendsen            @21 :Text;         # BERENDSEN.
  langevin             @22 :Bool = false; # LANGEVIN.
  annealing            @23 :Text;         # ANNEALING.
  quench               @24 :Bool = false; # QUENCH.
  rattle               @25 :Bool = false; # RATTLE.
  shake                @26 :Bool = false; # SHAKE.
  constraint           @27 :Text;         # CONSTRAINT.
  trotter              @28 :Text;         # TROTTER.
  restart              @29 :Bool = false; # RESTART.
  printOptions         @30 :Text;         # PRINT.
  storeOptions         @31 :Text;         # STORE.
  centerMoleculeOff    @32 :Bool = false; # CENTER MOLECULE OFF.
  centerMoleculeOn     @33 :Bool = false; # CENTER MOLECULE ON.
  diis                 @34 :Bool = false; # DIIS.
  odiis                @35 :Bool = false; # ODIIS.
  pcg                  @36 :Bool = false; # PCG.
  diagonalization      @37 :Bool = false; # DIAGONALIZATION.
  freeEnergy           @38 :Bool = false; # FREE-ENERGY.
  interface            @39 :Bool = false; # INTERFACE.
  qmmm                 @40 :Bool = false; # QMMM.
  bicanonicalEnsemble  @41 :Bool = false; # BICANONICAL ENSEMBLE.
  cdft                 @42 :Bool = false; # CDFT.
  properties           @43 :Bool = false; # PROPERTIES.
  vdwCorrection        @44 :Text;         # VDW CORRECTION argument, e.g. ON/OFF.
  vdwWannier           @45 :Text;         # VDW WANNIER argument, e.g. ON/OFF.
  dcacp                @46 :Bool = false; # DCACP.
  isolatedMolecule     @47 :Bool = false; # ISOLATED MOLECULE.
}

struct CPMDDftSection {
  functional    @0 :Text = "BLYP";        # FUNCTIONAL name (LDA, BLYP, PBE, ...).
  lsd           @1 :Bool = false;         # LSD (spin polarized).
  directives    @2 :List(CPMDDirective);
  gcCutoff      @3 :Float64 = 0.0;        # GC-CUTOFF; 0 => omit.
  xcDriver      @4 :Text;                 # XC_DRIVER.
  libxc         @5 :Text;                 # LIBXC.
  lrKernel      @6 :Text;                 # LR KERNEL.
  refunct       @7 :Text;                 # REFUNCT.
  mtsHighFunc   @8 :Text;                 # MTS_HIGH_FUNC.
  mtsLowFunc    @9 :Text;                 # MTS_LOW_FUNC.
  hfx           @10 :Bool = false;        # HFX.
  hfxScreening  @11 :Text;                # HFX-SCREENING.
  hubbard       @12 :Text;                # HUBBARD.
  alpha         @13 :Float64 = 0.0;       # ALPHA; 0 => omit.
  beta          @14 :Float64 = 0.0;       # BETA; 0 => omit.
  oldCode       @15 :Bool = false;        # OLDCODE.
  newCode       @16 :Bool = false;        # NEWCODE.
  correlation   @17 :Text;                # CORRELATION.
  exchange      @18 :Text;                # EXCHANGE.
  becke88       @19 :Bool = false;        # BECKE88.
}

struct CPMDAtomsPseudopotential {
  element @0 :Text;  # Element symbol, e.g. "H".
  path    @1 :Text;  # Pseudopotential file path or library token.
  lmax    @2 :Int32 = -1; # Optional LMAX; -1 => omit.
}

struct CPMDAtomsSection {
  pseudopotentials @0 :List(CPMDAtomsPseudopotential);
  # Explicit Cartesian coordinates are supplied per step via ForceInput; this
  # section only carries PP and fixed &ATOMS directives for deck rendering.
  directives       @1 :List(CPMDDirective);
}

struct CPMDDirectiveSection {
  directives @0 :List(CPMDDirective); # Keyword/value lines inside a named section.
  subsections @1 :List(CPMDGenericSection); # Nested NAME ... END NAME blocks.
}

enum CPMDSectionKind {
  generic @0;
  system  @1;
  cpmd    @2;
  dft     @3;
  atoms   @4;
  set     @5;
  raw     @6;
  atom    @7;
  basis   @8;
  clas    @9;
  eam     @10;
  exte    @11;
  hardness @12;
  info    @13;
  linres  @14;
  molstates @15;
  mts     @16;
  nlcc    @17;
  path    @18;
  pimd    @19;
  potential @20;
  prop    @21;
  ptddft  @22;
  resp    @23;
  tddft   @24;
  vdw     @25;
  vectors @26;
  wavefunction @27;
}

struct CPMDInputSection {
  union {
    generic @0 :CPMDGenericSection;
    system  @1 :CPMDSystemSection;
    cpmd    @2 :CPMDCpmdSection;
    dft     @3 :CPMDDftSection;
    atoms   @4 :CPMDAtomsSection;
    set     @5 :CPMDSetDirective;
    raw     @6 :Text; # Full &SECTION ... &END block text.
    atom    @7 :CPMDDirectiveSection;
    basis   @8 :CPMDDirectiveSection;
    clas    @9 :CPMDDirectiveSection;
    eam     @10 :CPMDDirectiveSection;
    exte    @11 :CPMDDirectiveSection;
    hardness @12 :CPMDDirectiveSection;
    info    @13 :CPMDDirectiveSection;
    linres  @14 :CPMDDirectiveSection;
    molstates @15 :CPMDDirectiveSection;
    mts     @16 :CPMDDirectiveSection;
    nlcc    @17 :CPMDDirectiveSection;
    path    @18 :CPMDDirectiveSection;
    pimd    @19 :CPMDDirectiveSection;
    potential @20 :CPMDDirectiveSection;
    prop    @21 :CPMDDirectiveSection;
    ptddft  @22 :CPMDDirectiveSection;
    resp    @23 :CPMDDirectiveSection;
    tddft   @24 :CPMDDirectiveSection;
    vdw     @25 :CPMDDirectiveSection;
    vectors @26 :CPMDDirectiveSection;
    wavefunction @27 :CPMDDirectiveSection;
  }
}

# @struct CPMDParams
# @brief CPMD-specific knobs (backend arm inside PotentialConfig / rgpot params).
#
# Not a standalone config language — only used when the active potential is CPMD.
# Same fields in/out via Cap'n Proto. Geometry for each evaluation is ForceInput.
struct CPMDParams {
  functional   @0 :Text = "BLYP";     # Default DFT functional (also in dft section).
  cutOffRy     @1 :Float64 = 70.0;    # Plane-wave cutoff in Rydberg.
  charge       @2 :Int32 = 0;
  multiplicity @3 :Int32 = 1;         # 2S+1; >1 enables LSD when rendering defaults.
  task         @4 :Text = "gradient"; # energy | gradient | md | optimize (frontend hint).
  title        @5 :Text = "";         # Optional comment header in rendered deck.
  memoryMb     @6 :UInt32 = 0;        # Frontend memory hint; 0 => environment defaults.
  scratchDir   @7 :Text = "";         # Fallback &CPMD FILEPATH directory.
  permanentDir @8 :Text = "";         # Preferred &CPMD FILEPATH restart directory.
  cpmdRoot     @9 :Text = "";         # OpenCPMD source/build tree; empty => env CPMD_ROOT.
  enginePath   @10 :Text = "";        # Frontend engine selection hint.
  inputBlocks  @11 :List(Text);       # Raw &SECTION blocks applied before generated ones.
  inputSections @12 :List(CPMDInputSection); # Structured CPMD input sections.
  # Long-tail options: use CPMDInputSection.raw, inputBlocks, or set.
}

# @struct PotentialConfig
# @brief rgpot / OmniPotentRPC user parameters (extensible, Cap'n Proto only).
#
# Tagged union: exactly one backend's options (or none). Geometry stays on
# ForceInput; this struct is method/backend setup only.
struct PotentialConfig {
  union {
    none   @0 :Void;
    nwchem @1 :Void;      # Reserved ordinal aligned with nwchemc; not implemented here.
    cpmd   @2 :CPMDParams; # OpenCPMD / cpmdc backend.
  }
}

# @interface Potential
# @brief RPC interface for remote calculations (same contract as nwchemc/rgpot).
interface Potential {
  calculate @0 (fip :ForceInput) -> (result :PotentialResult);
  configure @1 (config :PotentialConfig) -> (ok :Bool, message :Text);
}
