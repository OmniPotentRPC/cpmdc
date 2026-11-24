/* pair_style cpmdc: cpmdc as a LAMMPS pair style, on the public C ABI.
 *
 *   pair_style cpmdc params.bin [width 0.30]
 *   pair_coeff * * O H NULL
 *
 * Element entries are QM and NULL entries are Gaussian MM charges taken from
 * atom->q, with width in Angstrom. With no NULL entry the style evaluates the
 * whole system with cpmdc_session_calculate_au() once per force call. With
 * one, the QM atoms go to cpmdc_session_calculate_au() and the MM atoms to
 * cpmdc_session_set_external_charges(), the MM forces coming back from
 * cpmdc_session_last_external_charge_forces(). Electrostatic embedding is
 * isolated only, so the virial -dE/d(strain) is zero. It requests no neighbor list: every rank sends the
 * unwrapped positions of the atoms it owns, by atom ID, to world rank 0,
 * which owns the session, evaluates, and broadcasts the forces. The step is
 * ordered by increasing atom ID, so sparse IDs left by delete_atoms work.
 * Rank 0 adds the energy and the virial -dE/d(strain) to its tallies; LAMMPS
 * sums those over the ranks.
 *
 * The energy and the virial are totals of the whole system. A density
 * functional energy is not a sum of atom terms, and the virial is a cell
 * derivative, not sum r f, so neither has a unique per-atom split. A compute
 * that asks for either (pe/atom, stress/atom, centroid/stress/atom) is refused in init_style(),
 * before the first evaluation; per-atom forces are available.
 *
 * Units metal, real, and electron are converted; any other style is an
 * error. The cell is the restricted triclinic box with one lattice vector
 * per row: a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz). A fully
 * periodic box sends that cell; a box without periodic axes sends none, and
 * the CPMDParams deck sets the boundary conditions.
 */
#include "pair_cpmdc.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include "cpmdc.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <numeric>

using namespace LAMMPS_NS;

namespace {

/* Same values as CPMDC_BOHR_TO_ANGSTROM, cpmdc.HARTREE_TO_EV, and
 * cpmdc.HARTREE_TO_KCAL_PER_MOL (CODATA 2018 Hartree times the Avogadro
 * constant over 4184 J/kcal). */
constexpr double kBohrToAngstrom = 0.529177210903;
constexpr double kHartreeToEv = 27.211386245988;
constexpr double kHartreeToKcalPerMol = 4.3597447222071e-18 * 6.02214076e23 / 4184.0;

/* The library refuses a narrower Gaussian; refusing it here names the setting. */
constexpr double kMinExternalWidthAng = 1.0e-12;

constexpr int kMaxZ = 118;
const char *const kSymbols[kMaxZ] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",
    "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
    "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

/* An element symbol or an atomic number; 0 when it is neither. */
int element_number(const char *text)
{
  char *end = nullptr;
  long z = std::strtol(text, &end, 10);
  if (end != text && *end == '\0') return (z >= 1 && z <= kMaxZ) ? (int) z : 0;
  for (int i = 0; i < kMaxZ; ++i)
    if (std::strcmp(kSymbols[i], text) == 0) return i + 1;
  return 0;
}

bool unit_factors(const char *style, double &length_per_bohr, double &energy_per_hartree)
{
  if (std::strcmp(style, "metal") == 0) {
    length_per_bohr = kBohrToAngstrom;
    energy_per_hartree = kHartreeToEv;
  } else if (std::strcmp(style, "real") == 0) {
    length_per_bohr = kBohrToAngstrom;
    energy_per_hartree = kHartreeToKcalPerMol;
  } else if (std::strcmp(style, "electron") == 0) {
    length_per_bohr = 1.0;
    energy_per_hartree = 1.0;
  } else {
    return false;
  }
  return true;
}

}    // namespace

PairCPMDC::PairCPMDC(LAMMPS *lmp) :
    Pair(lmp), session(nullptr), length_per_bohr(1.0), energy_per_hartree(1.0),
    external_width_ang(0.30), external_active(false), adopted(0)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  no_virial_fdotr_compute = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
}

PairCPMDC::~PairCPMDC()
{
  close_session();
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

void PairCPMDC::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n, n, "pair:cutsq");
}

void PairCPMDC::close_session()
{
  if (session) cpmdc_session_destroy(session);
  session = nullptr;
  composition.clear();
  external_active = false;
}

/* Rank 0 only. A rank of a multi-rank run evaluates alone, so it first
 * makes MPI_COMM_SELF its calculator communicator. A library without MPI
 * refuses that and needs no communicator; an OpenCPMD build must accept. */
int PairCPMDC::open_session(std::string &why)
{
  if (comm->nprocs > 1 && !adopted) {
    MPI_Comm self = MPI_COMM_SELF;
    if (cpmdc_adopt_calculator_comm(&self, sizeof(self), 1) != 0 && cpmdc_available()) {
      const char *text = cpmdc_last_error();
      why = text && *text ? text : "cpmdc_adopt_calculator_comm failed";
      return -1;
    }
    adopted = 1;
  }
  session = cpmdc_session_create(params.data(), params.size());
  if (!session) {
    const char *text = cpmdc_last_error();
    why = text && *text ? text : "cpmdc_session_create failed";
    return -1;
  }
  return 0;
}

void PairCPMDC::settings(int narg, char **arg)
{
  if (narg != 1 && narg != 3)
    error->all(FLERR, "Illegal pair_style cpmdc command: use pair_style cpmdc <CPMDParams file> [width <Angstrom>]");
  external_width_ang = 0.30;
  if (narg == 3) {
    if (std::strcmp(arg[1], "width") != 0)
      error->all(FLERR, "pair_style cpmdc optional setting is width <Angstrom>");
    char *end = nullptr;
    external_width_ang = std::strtod(arg[2], &end);
    if (!end || *end != '\0' || !std::isfinite(external_width_ang) ||
        external_width_ang < kMinExternalWidthAng)
      error->all(FLERR, "pair_style cpmdc width must be finite and at least 1e-12 Angstrom");
  }
  if (!unit_factors(update->unit_style, length_per_bohr, energy_per_hartree))
    error->all(FLERR, "pair_style cpmdc supports units metal, real, and electron; this input uses units {}",
               update->unit_style);
  close_session();
  params_path = arg[0];
  params.clear();
  int ok = 1;
  if (comm->me == 0) {
    std::ifstream in(params_path, std::ios::binary);
    if (in) params.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    ok = !params.empty();
  }
  MPI_Bcast(&ok, 1, MPI_INT, 0, world);
  if (!ok) error->all(FLERR, "pair_style cpmdc cannot read a CPMDParams message from {}", params_path);
}

void PairCPMDC::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  const int ntypes = atom->ntypes;
  if (narg != 2 + ntypes || std::strcmp(arg[0], "*") != 0 || std::strcmp(arg[1], "*") != 0)
    error->all(FLERR, "pair_coeff for pair_style cpmdc is pair_coeff * * followed by one element per atom type ({} types)",
               ntypes);
  type_z.assign(ntypes + 1, 0);
  for (int i = 1; i <= ntypes; i++) {
    type_z[i] = std::strcmp(arg[i + 1], "NULL") == 0 ? 0 : element_number(arg[i + 1]);
    if (!type_z[i] && std::strcmp(arg[i + 1], "NULL") != 0)
      error->all(FLERR, "pair_coeff cpmdc: {} is not NULL, an element symbol, or an atomic number 1..118", arg[i + 1]);
  }
  if (std::none_of(type_z.begin() + 1, type_z.end(), [](int z) { return z != 0; }))
    error->all(FLERR, "pair_coeff cpmdc needs at least one QM atom type");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++) setflag[i][j] = 1;
}

void PairCPMDC::init_style()
{
  if (!atom->tag_enable) error->all(FLERR, "pair_style cpmdc needs atom IDs");
  const int periodic = domain->xperiodic + domain->yperiodic + domain->zperiodic;
  if (periodic != 0 && periodic != 3)
    error->all(FLERR, "pair_style cpmdc takes a fully periodic box or one without periodic axes");
  if (atom->natoms > MAXSMALLINT / 3) error->all(FLERR, "pair_style cpmdc: too many atoms");
  const bool have_mm =
      std::any_of(type_z.begin() + 1, type_z.end(), [](int z) { return z == 0; });
  if (have_mm && periodic)
    error->all(FLERR, "pair_style cpmdc: periodic electrostatic embedding is unsupported");
  if (have_mm && !atom->q_flag)
    error->all(FLERR, "pair_style cpmdc: MM atom types need an atom style with charges");
  for (int i = 0; i < modify->ncompute; i++) {
    const Compute *request = modify->compute[i];
    if (request->peatomflag || request->pressatomflag)
      error->all(FLERR, "pair_style cpmdc computes no per-atom energy or virial (compute {})",
                 request->id);
  }
  // No neighbor list: the whole system goes to cpmdc.
}

double PairCPMDC::init_one(int, int)
{
  return 0.0;
}

void PairCPMDC::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  if (eflag_atom || vflag_atom)
    error->all(FLERR, "pair_style cpmdc computes no per-atom energy or virial");

  const int natoms = (int) atom->natoms;
  const int nlocal = atom->nlocal;
  const tagint *tag = atom->tag;
  const int *type = atom->type;
  const imageint *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  const double *q = atom->q;

  std::vector<long long> local_tags(nlocal);
  std::vector<double> local_positions(3 * (size_t) nlocal);
  std::vector<double> local_charges(nlocal, 0.0);
  std::vector<int> local_owner(nlocal);
  int bad = 0;
  for (int i = 0; i < nlocal; i++) {
    if (tag[i] < 1) {
      bad = 1;
      continue;
    }
    local_tags[i] = static_cast<long long>(tag[i]);
    double xu[3];
    domain->unmap(x[i], image[i], xu);
    for (int d = 0; d < 3; d++)
      local_positions[3 * i + d] = xu[d] / length_per_bohr;
    if (q) local_charges[i] = q[i];
    local_owner[i] = type[i];
  }
  int any_bad = 0;
  MPI_Allreduce(&bad, &any_bad, 1, MPI_INT, MPI_MAX, world);
  if (any_bad) error->all(FLERR, "pair_style cpmdc needs positive atom IDs");

  std::vector<int> counts, offsets, coordinate_counts, coordinate_offsets;
  if (comm->me == 0) counts.resize(comm->nprocs);
  MPI_Gather(&nlocal, 1, MPI_INT, comm->me == 0 ? counts.data() : nullptr, 1, MPI_INT, 0,
             world);
  if (comm->me == 0) {
    offsets.resize(comm->nprocs);
    coordinate_counts.resize(comm->nprocs);
    coordinate_offsets.resize(comm->nprocs);
    int total = 0;
    for (int rank = 0; rank < comm->nprocs; ++rank) {
      offsets[rank] = total;
      coordinate_counts[rank] = 3 * counts[rank];
      coordinate_offsets[rank] = 3 * total;
      total += counts[rank];
    }
    bad = total != natoms;
    tags.resize(natoms);
    positions.resize(3 * (size_t) natoms);
    owner.resize(natoms);
    charges.resize(natoms);
  }
  MPI_Bcast(&bad, 1, MPI_INT, 0, world);
  if (bad) error->all(FLERR, "pair_style cpmdc received an incomplete atom set");
  MPI_Gatherv(local_tags.data(), nlocal, MPI_LONG_LONG, comm->me == 0 ? tags.data() : nullptr,
              comm->me == 0 ? counts.data() : nullptr, comm->me == 0 ? offsets.data() : nullptr,
              MPI_LONG_LONG, 0, world);
  MPI_Gatherv(local_positions.data(), 3 * nlocal, MPI_DOUBLE,
              comm->me == 0 ? positions.data() : nullptr,
              comm->me == 0 ? coordinate_counts.data() : nullptr,
              comm->me == 0 ? coordinate_offsets.data() : nullptr, MPI_DOUBLE, 0, world);
  MPI_Gatherv(local_charges.data(), nlocal, MPI_DOUBLE, comm->me == 0 ? charges.data() : nullptr,
              comm->me == 0 ? counts.data() : nullptr, comm->me == 0 ? offsets.data() : nullptr,
              MPI_DOUBLE, 0, world);
  MPI_Gatherv(local_owner.data(), nlocal, MPI_INT, comm->me == 0 ? owner.data() : nullptr,
              comm->me == 0 ? counts.data() : nullptr, comm->me == 0 ? offsets.data() : nullptr,
              MPI_INT, 0, world);

  // status, energy, and the row-major virial, from rank 0
  double head[11] = {0.0};
  char why[512] = {0};
  forces.assign(3 * (size_t) natoms, 0.0);
  if (comm->me == 0) {
    std::string problem;
    std::vector<int> order(natoms);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [this](int left, int right) {
      return tags[left] < tags[right];
    });
    for (int k = 1; k < natoms; ++k)
      if (tags[order[k - 1]] == tags[order[k]]) {
        problem = "pair_style cpmdc needs unique atom IDs";
        break;
      }
    std::vector<long long> ordered_tags(natoms);
    std::vector<double> ordered_positions(3 * (size_t) natoms);
    std::vector<double> ordered_charges(natoms);
    std::vector<int> ordered_owner(natoms);
    for (int k = 0; k < natoms; ++k) {
      const int source = order[k];
      ordered_tags[k] = tags[source];
      ordered_owner[k] = owner[source];
      ordered_charges[k] = charges[source];
      for (int d = 0; d < 3; ++d)
        ordered_positions[3 * k + d] = positions[3 * source + d];
    }
    tags.swap(ordered_tags);
    positions.swap(ordered_positions);
    charges.swap(ordered_charges);
    owner.swap(ordered_owner);
    std::vector<int> qm_index, mm_index;
    numbers.clear();
    for (int k = 0; k < natoms && problem.empty(); k++) {
      if (owner[k] < 1 || owner[k] > atom->ntypes)
        problem = "pair_style cpmdc received an invalid atom type";
      else if (type_z[owner[k]] == 0)
        mm_index.push_back(k);
      else {
        qm_index.push_back(k);
        numbers.push_back(type_z[owner[k]]);
      }
    }
    if (problem.empty() && qm_index.empty()) problem = "pair_style cpmdc found no QM atom";
    if (problem.empty()) {
      std::vector<int> sorted(numbers);
      std::sort(sorted.begin(), sorted.end());
      if (session && sorted != composition) close_session();
      if (!session && open_session(problem) == 0) composition = sorted;
    }
    if (problem.empty()) {
      const double *h = domain->h;
      const double cell[9] = {h[0], 0.0, 0.0, h[5], h[1], 0.0, h[4], h[3], h[2]};
      double cell_bohr[9];
      for (int i = 0; i < 9; i++) cell_bohr[i] = cell[i] / length_per_bohr;
      const bool periodic = domain->xperiodic && domain->yperiodic && domain->zperiodic;
      std::vector<double> qm_positions(3 * qm_index.size());
      for (size_t i = 0; i < qm_index.size(); ++i)
        std::copy_n(&positions[3 * (size_t) qm_index[i]], 3, &qm_positions[3 * i]);
      // The set is installed on every step that has MM atoms, and cleared on
      // the first step that has none after one that had some, so charges
      // of deleted or retyped atoms never reach a later evaluation.
      if (!mm_index.empty() || external_active) {
        std::vector<double> mm_positions(3 * mm_index.size());
        std::vector<double> mm_charges(mm_index.size());
        std::vector<double> widths(mm_index.size(), external_width_ang);
        for (size_t i = 0; i < mm_index.size(); ++i) {
          const int k = mm_index[i];
          for (int d = 0; d < 3; ++d)
            mm_positions[3 * i + d] = positions[3 * (size_t) k + d] * kBohrToAngstrom;
          mm_charges[i] = charges[k];
        }
        if (cpmdc_session_set_external_charges(session, mm_index.size(), mm_positions.data(),
                                               mm_charges.data(), widths.data()) != 0) {
          const char *text = cpmdc_last_error();
          problem = text && *text ? text : "setting external charges failed";
        } else {
          external_active = !mm_index.empty();
        }
      }
      std::vector<double> qm_forces(3 * qm_index.size());
      if (problem.empty() &&
          cpmdc_session_calculate_au(session, (int) qm_index.size(), qm_positions.data(),
                                     numbers.data(), periodic ? cell_bohr : nullptr, &head[1],
                                     qm_forces.data(), &head[2]) == 0) {
        for (size_t i = 0; i < qm_index.size(); ++i)
          std::copy_n(&qm_forces[3 * i], 3, &forces[3 * (size_t) qm_index[i]]);
        if (!mm_index.empty()) {
          std::vector<double> mm_forces(3 * mm_index.size());
          if (cpmdc_session_last_external_charge_forces(session, mm_forces.data(),
                                                        mm_forces.size()) != 0) {
            const char *text = cpmdc_last_error();
            problem = text && *text ? text : "external charge forces are unavailable";
          } else {
            for (size_t i = 0; i < mm_index.size(); ++i)
              std::copy_n(&mm_forces[3 * i], 3, &forces[3 * (size_t) mm_index[i]]);
          }
        }
        if (problem.empty()) head[0] = 1.0;
      } else if (problem.empty()) {
        const char *text = cpmdc_last_error();
        problem = text && *text ? text : "evaluation failed";
      }
    }
    std::strncpy(why, problem.c_str(), sizeof(why) - 1);
  }
  MPI_Bcast(head, 11, MPI_DOUBLE, 0, world);
  if (head[0] == 0.0) {
    MPI_Bcast(why, (int) sizeof(why), MPI_CHAR, 0, world);
    error->all(FLERR, "pair_style cpmdc: {}", why);
  }
  if (comm->me != 0) tags.resize(natoms);
  MPI_Bcast(tags.data(), natoms, MPI_LONG_LONG, 0, world);
  MPI_Bcast(forces.data(), 3 * natoms, MPI_DOUBLE, 0, world);

  const double force_unit = energy_per_hartree / length_per_bohr;
  bad = 0;
  for (int i = 0; i < nlocal; i++) {
    const long long id = static_cast<long long>(tag[i]);
    const std::vector<long long>::const_iterator found =
        std::lower_bound(tags.begin(), tags.end(), id);
    if (found == tags.end() || *found != id) bad = 1;
  }
  MPI_Allreduce(&bad, &any_bad, 1, MPI_INT, MPI_MAX, world);
  if (any_bad)
    error->all(FLERR, "pair_style cpmdc lost an atom ID while distributing forces");
  for (int i = 0; i < nlocal; i++) {
    const long long id = static_cast<long long>(tag[i]);
    const std::vector<long long>::const_iterator found =
        std::lower_bound(tags.begin(), tags.end(), id);
    const int k = (int) std::distance(tags.cbegin(), found);
    for (int d = 0; d < 3; d++) f[i][d] += forces[3 * k + d] * force_unit;
  }
  if (comm->me == 0) {
    const double *w = &head[2];
    if (eflag_global) eng_vdwl += head[1] * energy_per_hartree;
    if (vflag_global) {
      virial[0] += w[0] * energy_per_hartree;
      virial[1] += w[4] * energy_per_hartree;
      virial[2] += w[8] * energy_per_hartree;
      virial[3] += 0.5 * (w[1] + w[3]) * energy_per_hartree;
      virial[4] += 0.5 * (w[2] + w[6]) * energy_per_hartree;
      virial[5] += 0.5 * (w[5] + w[7]) * energy_per_hartree;
    }
  }
}
