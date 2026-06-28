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

struct CPMDCouplingSurface {
  stateI      @0 :Int32 = 0;      # First coupled KS state index.
  stateJ      @1 :Int32 = 0;      # Second coupled KS state index.
  coefficient @2 :Float64 = 0.0;  # Coupling coefficient for NSURF.
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
  lowSpinExcitationPenalty @66 :Float64 = 0.0; # LOW SPIN EXCITATION PENALTY; 0 => omit.
  cdftDonorAtoms @67 :List(Int32);             # DONOR atom indexes.
  cdftDonorWeights @68 :List(Int32);           # DONOR WMULT integer weights.
  cdftAcceptorAtoms @69 :List(Int32);          # ACCEPTOR atom indexes.
  cdftAcceptorHdasDonors @70 :List(Int32);     # ACCEPTOR HDAS donor atom indexes.
  cdftAcceptorWeights @71 :List(Int32);        # ACCEPTOR WMULT integer weights.
  couplingsFiniteDifference @72 :Bool = false; # COUPLINGS FD.
  couplingsFiniteDifferenceDisplacement @73 :Float64 = 0.0; # COUPLINGS FD=eps; 0 => omit value.
  couplingsProductDisplacement @74 :Float64 = 0.0; # COUPLINGS PROD=eps; 0 => omit.
  couplingsLinres @75 :Bool = false;           # COUPLINGS LINRES.
  couplingsLinresTolerance @76 :Float64 = 0.0; # COUPLINGS LINRES TOL=; 0 => omit.
  couplingsLinresNvects @77 :Int32 = 0;        # COUPLINGS LINRES NVECT=; 0 => omit.
  couplingsLinresSpecify @78 :Bool = false;    # COUPLINGS LINRES NVECT=... SPECIFY.
  couplingsLinresBruteForce @79 :Bool = false; # COUPLINGS LINRES BRUTE FORCE.
  couplingsLinresThresholds @80 :List(Float64); # COUPLINGS LINRES THRESHOLDS low/med/high pairs.
  couplingsSurfaces @81 :List(CPMDCouplingSurface); # COUPLINGS NSURF surface triples.
  couplingsFiniteDifferenceAtoms @82 :List(Int32); # COUPLINGS NAT atom indexes.
  cellAbsolute @83 :Bool = false;          # CELL ABSOLUTE.
  cellDegree @84 :Bool = false;            # CELL DEGREE.
  referenceCellAbsolute @85 :Bool = false; # REFERENCE CELL ABSOLUTE.
  referenceCellDegree @86 :Bool = false;   # REFERENCE CELL DEGREE.
  classicalCellAbsolute @87 :Bool = false; # CLASSICAL CELL ABSOLUTE.
  classicalCellDegree @88 :Bool = false;   # CLASSICAL CELL DEGREE.
  cellVectors @89 :Bool = false;            # CELL VECTORS.
  referenceCellVectors @90 :Bool = false;   # REFERENCE CELL VECTORS.
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
  maxRuntime           @48 :Float64 = 0.0; # MAXRUNTIME; 0 => omit.
  timestepElectrons    @49 :Float64 = 0.0; # TIMESTEP ELECTRONS; 0 => omit.
  timestepIons         @50 :Float64 = 0.0; # TIMESTEP IONS; 0 => omit.
  cellMass             @51 :Float64 = 0.0; # CMASS; 0 => omit.
  temperatureElectron  @52 :Float64 = 0.0; # TEMPERATURE ELECTRON; 0 => omit.
  temperature          @53 :Float64 = 0.0; # TEMPERATURE; 0 => omit.
  temperatureRamp      @54 :Bool = false; # TEMPERATURE RAMP.
  temperatureRampTime  @55 :Float64 = 0.0; # TEMPERATURE RAMP second value.
  temperatureRampRate  @56 :Float64 = 0.0; # TEMPERATURE RAMP third value.
  rescaleOldVelocities @57 :Bool = false; # RESCALE OLD VELOCITIES.
  reverseVelocities    @58 :Bool = false; # REVERSE VELOCITIES.
  subtractComVelocity  @59 :Int32 = 0; # SUBTRACT COMVEL; 0 => omit.
  subtractRotVelocity  @60 :Int32 = 0; # SUBTRACT ROTVEL; 0 => omit.
  prngSeed             @61 :Int32 = 0; # PRNGSEED; 0 => omit.
  tempControlIons      @62 :Text; # TEMPCONTROL IONS payload.
  tempControlElectrons @63 :Text; # TEMPCONTROL ELECTRONS payload.
  tempControlCell      @64 :Text; # TEMPCONTROL CELL payload.
  berendsenIons        @65 :Text; # BERENDSEN IONS payload.
  berendsenElectrons   @66 :Text; # BERENDSEN ELECTRONS payload.
  berendsenCell        @67 :Text; # BERENDSEN CELL payload.
  noseIonsThermostat      @68 :Text; # NOSE IONS payload.
  noseElectronsThermostat @69 :Text; # NOSE ELECTRONS payload.
  noseCellThermostat      @70 :Text; # NOSE CELL payload.
  noseParameters          @71 :Text; # NOSE PARAMETERS payload.
  convergenceCell         @72 :Float64 = 0.0; # CONVERGENCE CELL; 0 => omit.
  convergenceAdapt        @73 :Float64 = 0.0; # CONVERGENCE ADAPT; 0 => omit.
  convergenceEnergy       @74 :Float64 = 0.0; # CONVERGENCE ENERGY; 0 => omit.
  convergenceCalfor       @75 :Float64 = 0.0; # CONVERGENCE CALFOR; 0 => omit.
  convergenceRelax        @76 :Int32 = 0; # CONVERGENCE RELAX; 0 => omit.
  convergenceRhofix       @77 :Float64 = 0.0; # CONVERGENCE RHOFIX; 0 => omit.
  convergenceInitial      @78 :Float64 = 0.0; # CONVERGENCE INITIAL; 0 => omit.
  convergenceConstraint   @79 :Text; # CONVERGENCE CONSTRAINT payload.
  annealingIons           @80 :Float64 = 0.0; # ANNEALING IONS; 0 => omit.
  annealingElectrons      @81 :Float64 = 0.0; # ANNEALING ELECTRONS; 0 => omit.
  annealingCell           @82 :Float64 = 0.0; # ANNEALING CELL; 0 => omit.
  dampingIons             @83 :Float64 = 0.0; # DAMPING IONS; 0 => omit.
  dampingElectrons        @84 :Float64 = 0.0; # DAMPING ELECTRONS; 0 => omit.
  dampingCell             @85 :Float64 = 0.0; # DAMPING CELL; 0 => omit.
  hessian                 @86 :Text; # HESSIAN inline options.
  project                 @87 :Text; # PROJECT inline option.
  stressTensorSample      @88 :Int32 = 0; # STRESS TENSOR sample; 0 => omit.
  stressTensorVirial      @89 :Bool = false; # STRESS TENSOR VIRIAL.
  classStressSample       @90 :Int32 = 0; # CLASSTRESS sample; 0 => omit.
  storeSelection          @91 :Text; # STORE inline selection.
  storeInterval           @92 :Int32 = 0; # STORE interval; 0 => omit.
  storeSelfConsistentInterval @93 :Int32 = 0; # STORE SC interval; 0 => omit.
  storeOffSelection       @94 :Text; # STORE OFF inline selection.
  restFileCount           @95 :Int32 = 0; # RESTFILE count; 0 => omit.
  restFileSample          @96 :Text; # RESTFILE SAMPLE payload.
  trajectoryOptions       @97 :Text; # TRAJECTORY inline options.
  trajectorySample        @98 :Int32 = 0; # TRAJECTORY SAMPLE value; 0 => omit.
  trajectoryRange         @99 :Text; # TRAJECTORY RANGE payload.
  movieSample             @100 :Int32 = 0; # MOVIE SAMPLE value; 0 => omit.
  movieOff                @101 :Bool = false; # MOVIE OFF.
  energyBands             @102 :Bool = false; # ENERGYBANDS.
  externalPotential       @103 :Bool = false; # EXTERNAL POTENTIAL.
  externalPotentialAdd    @104 :Bool = false; # EXTERNAL POTENTIAL ADD.
  electrostaticPotential  @105 :Bool = false; # ELECTROSTATIC POTENTIAL.
  electrostaticPotentialSample @106 :Int32 = 0; # ELECTROSTATIC POTENTIAL SAMPLE.
  dipoleDynamicsSample    @107 :Int32 = 0; # DIPOLE DYNAMICS SAMPLE.
  dipoleDynamicsWannier   @108 :Bool = false; # DIPOLE DYNAMICS WANNIER.
  rhoOut                  @109 :Bool = false; # RHOOUT.
  rhoOutSample            @110 :Int32 = 0; # RHOOUT SAMPLE.
  rhoOutBandsCount        @111 :Int32 = 0; # RHOOUT BANDS count.
  rhoOutBands             @112 :Text; # RHOOUT BANDS payload.
  elf                     @113 :Bool = false; # ELF.
  elfParameters           @114 :Text; # ELF PARAMETER payload.
  wannierParameters       @115 :Text; # WANNIER PARAMETER payload.
  wannierOptimization     @116 :Text; # WANNIER OPTIMIZATION option.
  wannierType             @117 :Text; # WANNIER TYPE option.
  wannierReference        @118 :Text; # WANNIER REFERENCE payload.
  wannierSerial           @119 :Bool = false; # WANNIER SERIAL.
  wannierDos              @120 :Bool = false; # WANNIER DOS.
  wannierMolecular        @121 :Bool = false; # WANNIER MOLECULAR.
  wannierWfnOutOptions    @122 :Text; # WANNIER WFNOUT inline options.
  wannierWfnOutPayload    @123 :Text; # WANNIER WFNOUT payload.
  compress                @124 :Text; # COMPRESS inline option.
  memory                  @125 :Text; # MEMORY inline options.
  realSpaceWfnKeep        @126 :Bool = false; # REAL SPACE WFN KEEP.
  realSpaceWfnSize        @127 :Float64 = 0.0; # REAL SPACE WFN SIZE; 0 => omit.
  splineOptions           @128 :Text; # SPLINE inline options.
  splinePoints            @129 :Int32 = 0; # SPLINE POINTS payload; 0 => omit.
  splineRange             @130 :Float64 = 0.0; # SPLINE RANGE payload; 0 => omit.
  finiteDifferences       @131 :Text; # FINITE DIFFERENCES payload.
  taskGroups              @132 :Text; # TASKGROUPS inline option.
  taskGroupsCount         @133 :Int32 = 0; # TASKGROUPS count; 0 => omit.
  distributeFnl           @134 :Text; # DISTRIBUTE FNL ON/OFF.
  filePath                @135 :Text; # FILEPATH payload.
  benchmark               @136 :Text; # BENCHMARK payload.
  mirror                  @137 :Bool = false; # MIRROR.
  shiftPotential          @138 :Text; # SHIFT POTENTIAL payload.
  glocalizationParameters @139 :Text; # GLOCALIZATION PARAMETERS payload.
  glocalizationOptimization @140 :Text; # GLOCALIZATION OPTIMIZATION option.
  gfunctionalType         @141 :Text; # GFUNCTIONAL TYPE option.
  spreadRspace            @142 :Text; # SPREAD RSPACE payload after '='.
  gUnitarityOptions       @143 :Text; # PIPPO unitarity options.
  stepTuning              @144 :Bool = false; # STEP TUNING.
  gAntisym                @145 :Bool = false; # G_ANTISYM.
  gAntisymPenalty         @146 :Bool = false; # G_ANTISYM PENALTY.
  gKick                   @147 :Bool = false; # G_KICK.
  gComplex                @148 :Bool = false; # G_COMPLEX.
  gReal                   @149 :Bool = false; # G_REAL.
  readMatrix              @150 :Bool = false; # READ MATRIX.
  gStepTune               @151 :Bool = false; # G_STEP TUNE.
  glocWfnOutOptions       @152 :Text; # GLOC WFNOUT inline options.
  glocWfnOutPayload       @153 :Text; # GLOC WFNOUT payload line.
  noGeoCheck              @154 :Bool = false; # NO_GEO_CHECK.
  brokenSymmetry          @155 :Bool = false; # BROKEN symmetry.
  distributedLinalg       @156 :Text; # DISTRIBUTED LINALG ON/OFF.
  linalgNewOrtho          @157 :Text; # LINALG NEWORTHO ON/OFF.
  disorthoBlockSize       @158 :Int32 = 0; # DISORTHO_BSIZE; 0 => omit.
  statesBlockSize         @159 :Int32 = 0; # BLOCKSIZE STATES; 0 => omit.
  allToAllPrecision       @160 :Text; # ALLTOALL precision option.
  gshell                  @161 :Bool = false; # GSHELL.
  localPotential          @162 :Bool = false; # LOCAL POTENTIAL.
  cdftOptions             @163 :Text; # CDFT inline options.
  cdftPayload             @164 :Text; # CDFT payload line.
  cdftHdaPayload          @165 :Text; # CDFT HDA payload line.
  vgfactor                @166 :Text; # VGFACTOR payload line.
  vMirror                 @167 :Bool = false; # VMIRROR.
  combineSystemsOptions   @168 :Text; # COMBINE SYSTEMS inline options.
  combineSystemsPayload   @169 :Text; # COMBINE SYSTEMS payload line.
  combineSystemsSabPayload @170 :Text; # COMBINE SYSTEMS SAB payload line.
  kshamOptions            @171 :Text; # KSHAM inline options.
  kshamPayload            @172 :Text; # KSHAM payload line.
  czonesSet               @173 :Text; # CZONES SET payload lines.
  woutOptions             @174 :Text; # WOUT inline options.
  woutPayload             @175 :Text; # WOUT payload line.
  xfmqcTrajectories       @176 :Int32 = 0; # XFMQC trajectory count; 0 => omit.
  molecularDynamicsFileOptions @177 :Text; # MOLECULAR DYNAMICS FILE inline options.
  molecularDynamicsBdTrajectories @178 :Int32 = 0; # MOLECULAR DYNAMICS BD payload.
  parrinelloRahmanOptions @179 :Text; # PARRINELLO-RAHMAN inline options.
  optimizeGeometryOptions @180 :Text; # OPTIMIZE GEOMETRY inline options.
  optimizeGeometrySample  @181 :Int32 = 0; # OPTIMIZE GEOMETRY SAMPLE payload.
  optimizeCombinedOptions @182 :Text; # OPTIMIZE COMBINED inline options.
  optimizeCombinedSample  @183 :Int32 = 0; # OPTIMIZE COMBINED SAMPLE payload.
  cheby                   @184 :Bool = false; # CHEBY.
  cayley                  @185 :Bool = false; # CAYLEY.
  rungeKutta              @186 :Bool = false; # RUNGE-KUTTA.
  forceMatch              @187 :Bool = false; # FORCEMATCH.
  debugOptions            @188 :Text; # DEBUG inline options.
  kohnShamEnergiesOptions @189 :Text; # KOHN-SHAM ENERGIES inline options.
  kohnShamEnergiesCount   @190 :Int32 = 0; # KOHN-SHAM ENERGIES payload count.
  surfaceHoppingOptions   @191 :Text; # SURFACE HOPPING inline options.
  roksOptions             @192 :Text; # ROKS inline options.
  roksExpertPayload       @193 :Text; # ROKS EXPERT payload line.
  pathSampling            @194 :Bool = false; # PATH SAMPLING.
  fixrhoUpwfnOptions      @195 :Text; # FIXRHO UPWFN inline options.
  fixrhoVectors           @196 :Int32 = 0; # FIXRHO VECT payload; 0 => omit.
  fixrhoLoop              @197 :Text; # FIXRHO LOOP min/max payload.
  fixrhoWftol             @198 :Float64 = 0.0; # FIXRHO WFTOL; 0 => omit.
  bogoliubovCorrection    @199 :Text; # BOGOLIUBOV CORRECTION inline option.
  vibrationalAnalysisOptions @200 :Text; # VIBRATIONAL ANALYSIS inline options.
  vibrationalAnalysisSample @201 :Int32 = 0; # VIBRATIONAL ANALYSIS SAMPLE payload.
  vibrationalAnalysisMode @202 :Int32 = 0; # VIBRATIONAL ANALYSIS MODE= value.
  electronicSpectra       @203 :Bool = false; # ELECTRONIC SPECTRA.
  spinOrbitCouplingStates @204 :List(Int32); # SPIN-ORBIT COUPLING state pair.
  propagationSpectra      @205 :Bool = false; # PROPAGATION SPECTRA.
  propagationDistrub      @206 :Bool = false; # PROPAGATION DISTRUB.
  gaugePulse              @207 :Bool = false; # GAUGEPULSE.
  gaugeFieldFrequency     @208 :Float64 = 0.0; # GAUGEFIELD payload; 0 => omit.
  nacv                    @209 :Bool = false; # NACV.
  orbitalHardnessOptions  @210 :Text; # ORBITAL HARDNESS inline options.
  pathIntegral            @211 :Bool = false; # PATH INTEGRAL.
  pathMinimization        @212 :Bool = false; # PATH MINIMIZATION.
  langevinOptions         @213 :Text; # LANGEVIN inline options.
  langevinParameter       @214 :Text; # LANGEVIN payload line.
  qmmmEasy                @215 :Bool = false; # QMMMEASY.
  interfaceOptions        @216 :Text; # INTERFACE inline options.
  trotterFactorCount      @217 :Int32 = 0; # TROTTER FACTOR= count; 0 => omit.
  trotterFactorPayload    @218 :Text; # TROTTER FACTOR payload lines.
  linearResponse          @219 :Bool = false; # LINEAR RESPONSE.
  harmonicReference       @220 :Text; # HARMONIC REFERENCE option.
  scaledMasses            @221 :Text; # SCALED MASSES option.
  tddft                   @222 :Bool = false; # TDDFT.
  ssic                    @223 :Float64 = 0.0; # SSIC payload; 0 => omit.
  nonorthogonalOrbitalsOptions @224 :Text; # NONORTHOGONAL ORBITALS option.
  nonorthogonalOrbitalsLimit @225 :Float64 = 0.0; # NONORTHOGONAL ORBITALS payload.
  lanczosDiagonalizationOptions @226 :Text; # LANCZOS DIAGONALIZATION inline options.
  lanczosParametersCount  @227 :Int32 = 0; # LANCZOS PARAMETERS N= count.
  lanczosParametersPayload @228 :Text; # LANCZOS PARAMETERS payload lines.
  davidsonDiagonalization @229 :Bool = false; # DAVIDSON DIAGONALIZATION.
  davidsonParameters      @230 :Text; # DAVIDSON PARAMETERS payload line.
  alexanderMixing         @231 :Float64 = 0.0; # ALEXANDER MIXING payload; 0 => omit.
  andersonMixingGspace    @232 :Bool = false; # ANDERSON MIXING G-SPACE.
  andersonMixingCount     @233 :Int32 = 0; # ANDERSON MIXING N= count.
  andersonMixingPayload   @234 :Text; # ANDERSON MIXING payload lines.
  broydenMixingOptions    @235 :Text; # BROYDEN MIXING inline options.
  broydenMixingPayload    @236 :Text; # BROYDEN MIXING payload line.
  diisMixingCount         @237 :Int32 = 0; # DIIS MIXING N= count.
  diisMixingPayload       @238 :Text; # DIIS MIXING payload lines.
  moverhoMixing           @239 :Float64 = 0.0; # MOVERHO payload; 0 => omit.
  extrapolateWfnOptions   @240 :Text; # EXTRAPOLATE WFN inline options.
  extrapolateWfnOrder     @241 :Int32 = 0; # EXTRAPOLATE WFN payload; 0 => omit.
  extrapolateConstraintOrder @242 :Int32 = 0; # EXTRAPOLATE CONSTRAINT payload.
  tsdeOptions             @243 :Text; # TSDE inline options.
  tsdpOptions             @244 :Text; # TSDP inline options.
  tcgp                    @245 :Bool = false; # TCGP.
  tsdc                    @246 :Bool = false; # TSDC.
  steepestDescentOptions  @247 :Text; # STEEPEST DESCENT inline options.
  conjugateGradientOptions @248 :Text; # CONJUGATE GRADIENT inline options.
  odiisOptions            @249 :Text; # ODIIS inline options.
  odiisVectors            @250 :Int32 = 0; # ODIIS payload; 0 => omit.
  hamiltonianCutoff       @251 :Float64 = 0.0; # HAMILTONIAN CUTOFF payload.
  gdiisVectors            @252 :Int32 = 0; # GDIIS payload; 0 => omit.
  lbfgsOptions            @253 :Text; # LBFGS inline options.
  lbfgsPayload            @254 :Text; # LBFGS payload line.
  prfoOptions             @255 :Text; # PRFO inline options.
  prfoPayload             @256 :Text; # PRFO payload line.
  hesscore                @257 :Bool = false; # HESSCORE.
  bfgs                    @258 :Bool = false; # BFGS.
  rfoOrder                @259 :Int32 = 0; # RFO ORDER= value; 0 => omit.
  inrParametersCount      @260 :Int32 = 0; # INR PARAMETERS N= count.
  inrParametersPayload    @261 :Text; # INR PARAMETERS payload lines.
  implicitNewtonOptions   @262 :Text; # IMPLICIT NEWTON inline options.
  implicitNewtonMaxIter   @263 :Int32 = 0; # IMPLICIT NEWTON payload; 0 => omit.
  mixsd                   @264 :Float64 = 0.0; # MIXSD payload; 0 => omit.
  mixdiis                 @265 :Float64 = 0.0; # MIXDIIS payload; 0 => omit.
  restartOptions          @266 :Text; # RESTART inline options.
  intFileOptions          @267 :Text; # INTFILE inline options.
  intFileName             @268 :Text; # INTFILE filename payload.
  initializeWavefunctionOptions @269 :Text; # INITIALIZE WAVEFUNCTION inline options.
  rattleParameters        @270 :Text; # RATTLE payload line.
  orthogonalizationOptions @271 :Text; # ORTHOGONALIZATION inline options.
  quenchOptions           @272 :Text; # QUENCH inline options.
  randomizeOptions        @273 :Text; # RANDOMIZE inline options.
  randomizeAmplitude      @274 :Float64 = 0.0; # RANDOMIZE payload; 0 => omit.
  useMts                  @275 :Bool = false; # USE_MTS.
  nabdyZmax               @276 :Int32 = 0; # NABDY_ZMAX payload; 0 => omit.
  nabdySoft               @277 :Float64 = 0.0; # NABDY_SOFT payload; 0 => omit.
  nabdyRedistributeAmplitude @278 :Bool = false; # NABDY_REDISTR_AMPLI.
  nabdyScaleP             @279 :Bool = false; # NABDY_SCALEP.
  nabdyThermo             @280 :Text; # NABDY_THERMO payload line.
  noseIonsOptions         @281 :Text; # NOSE IONS inline options.
  useInStream             @282 :Bool = false; # USE_IN_STREAM.
  useOutStream            @283 :Bool = false; # USE_OUT_STREAM.
  useCublas               @284 :Bool = false; # USE_CUBLAS.
  useCufft                @285 :Bool = false; # USE_CUFFT.
  blasNStreamsPerDevice   @286 :Int32 = 0; # BLAS_N_STREAMS_PER_DEVICE payload.
  blasNDevicesPerTask     @287 :Int32 = 0; # BLAS_N_DEVICES_PER_TASK payload.
  fftNStreamsPerDevice    @288 :Int32 = 0; # FFT_N_STREAMS_PER_DEVICE payload.
  fftNDevicesPerTask      @289 :Int32 = 0; # FFT_N_DEVICES_PER_TASK payload.
  useMpiIo                @290 :Bool = false; # USE_MPI_IO.
  traceOptions            @291 :Text; # TRACE inline options.
  traceProcedure          @292 :Text; # TRACE_PROCEDURE payload line.
  traceMaxDepth           @293 :Int32 = 0; # TRACE_MAX_DEPTH payload.
  traceMaxCalls           @294 :Int32 = 0; # TRACE_MAX_CALLS payload.
  structureOptions        @295 :Text; # STRUCTURE inline options.
  structureSelection      @296 :Text; # STRUCTURE SELECT payload lines.
  wannierNproc            @297 :Int32 = 0; # WANNIER NPROC payload.
  wannierRelocalizeInScf  @298 :Bool = false; # WANNIER RELOCALIZE_IN_SCF.
  wannierRecomputeDipoleMatricesEvery @299 :Int32 = 0; # WANNIER RECOMPUTE_DIPOLE_MATRICES_EVERY payload.
  wannierRelocalizeEvery  @300 :Int32 = 0; # WANNIER RELOCALIZE_EVERY payload.
  paraUseMpiInPlace       @301 :Bool = false; # PARA_USE_MPI_IN_PLACE.
  paraBuffSize            @302 :Int32 = 0; # PARA_BUFF_SIZE payload.
  paraStackBuffSize       @303 :Int32 = 0; # PARA_STACK_BUFF_SIZE payload.
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
