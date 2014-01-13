#ifndef _MPOGEN_H
#define _MPOGEN_H

#include "include.h"
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"
#include <string.h>
#include <vector>
#include <fstream>

namespace btas { typedef SpinQuantum Quantum; };

using namespace btas;
using namespace mpsxx;

using std::string;
using std::vector;

class MPOGen {
protected:
  string input;
public:
  MPOGen(const char* m_input): input(m_input) {
    cout << "Read input file: " << input << endl;
  }
  virtual const MPO<Quantum> generate(int M = 0) = 0;
  virtual ~MPOGen() {}
};

enum class Spin: int {Up = 1, Down = -1};
void physical(Qshapes<Quantum>&, Dshapes&);
void zero(MPO<Quantum>&);
MPO<Quantum> create_op(const ColumnVector&, Spin);
MPO<Quantum> anni_op(const ColumnVector&, Spin);
Matrix rdm(const MPS<Quantum>&, Spin); // <a_i\sigma^\dagger a_j\sigma>
Matrix kappa(const MPS<Quantum>& A, bool symm = true); // <a_i\up a_j\down>

class MPOGen_Hubbard_BCS: public MPOGen {
private:
  // functions
  void read_matrix(Matrix&, string, std::ifstream&);
public:
  MPOGen_Hubbard_BCS(const char* m_input);
  virtual const MPO<Quantum> generate(int M = 0) final;
  ~MPOGen_Hubbard_BCS() {}
private:
  // data
  int nsite, ntei; // number of sites, number of 2body terms
  double U;
  Matrix H0, D0;
  vector<Matrix> f, g, d;
};

class MPOGen_Hubbard_UBCS: public MPOGen {
private:
  // functions
  void read_matrix(Matrix&, string, std::ifstream&);
public:
  MPOGen_Hubbard_UBCS(const char* m_input);
  virtual const MPO<Quantum> generate(int M = 0) final;
  ~MPOGen_Hubbard_UBCS() {}
private:
  // data
  int nsite, ntei; // number of sites, number of 2body terms
  double U;
  Matrix H0a, H0b, D0;
  vector<Matrix> fa, fb, ga, gb, da, db;
};

#endif
