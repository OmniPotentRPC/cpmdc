/* pair_style cpmdc: cpmdc as a LAMMPS pair style, on the public C ABI. */
#ifdef PAIR_CLASS
// clang-format off
PairStyle(cpmdc,PairCPMDC);
// clang-format on
#else

#ifndef LMP_PAIR_CPMDC_H
#define LMP_PAIR_CPMDC_H

#include "pair.h"

#include <string>
#include <vector>

struct CPMDCSession;

namespace LAMMPS_NS {

class PairCPMDC : public Pair {
public:
  PairCPMDC(class LAMMPS *);
  ~PairCPMDC() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

private:
  std::string params_path;
  std::vector<unsigned char> params;
  std::vector<int> type_z;
  std::vector<int> composition;
  CPMDCSession *session;
  double length_per_bohr, energy_per_hartree;
  int adopted;

  std::vector<double> positions, forces;
  std::vector<int> numbers, owner;

  void allocate();
  void close_session();
  int open_session(std::string &why);
};

}    // namespace LAMMPS_NS

#endif
#endif
