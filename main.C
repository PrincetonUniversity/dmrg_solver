#include <iostream>
#include <iomanip>
#include <fstream>

#include "mpogen.h"
#include "mpsite.h"
#include "driver.h"
#include "dmrg.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <cstdlib>

using std::cout;
using std::endl;
using std::string;
using std::istringstream;
using std::ifstream;

using namespace prototype;
using namespace mpsxx;

double random_gen() { return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; }

void init(MpStorages& sites, MPS<Quantum>& A, int M) {
  int L = sites.size();
  for (int i = 0; i < L-1; ++i) {
    sites[i].wfnc = std::move(A[i]);
    A[i].clear();
  }

  auto qlast = A[L-1].qshape();
  qlast[2] = {A[L-1].q()};  
  Quantum qt = -(A[L-1].qshape()[2][0]);
  sites[L-1].wfnc.resize(qt, qlast, A[L-1].dshape());

  for (auto it = A[L-1].begin(); it != A[L-1].end(); ++it) {
    sites[L-1].wfnc.insert(it -> first, *(it -> second));
  }
  A[L-1].clear();
  //
  // canonicalize & renormalize
  //
  Qshapes<Quantum> qz(1, Quantum::zero());
  Dshapes dz(qz.size(), 1);  
  auto qshape = make_array( qz, qz, qz);
  auto dshape = make_array( dz, dz, dz);

  sites[L-1].ropr.resize(Quantum::zero(), qshape, dshape);
  sites[L-1].ropr = 1.0;
  for(int i = L-1; i > 0; --i) {
    QSDnormalize(sites[i].wfnc);
    prototype::Canonicalize(0, sites[i].wfnc, sites[i].rmps, M);
    QSDcopy(sites[i-1].wfnc, sites[i-1].lmps);
    prototype::ComputeGuess(0, sites[i].rmps, sites[i].wfnc, sites[i-1].lmps, sites[i-1].wfnc);
    sites[i-1].ropr.clear();
    Renormalize (0, sites[i].mpo, sites[i].ropr, sites[i].rmps, sites[i].rmps, sites[i-1].ropr);
  }
  btas::QSDnormalize(sites[0].wfnc);
  sites[0].lopr.resize(Quantum::zero(), qshape, dshape);
  sites[0].lopr = 1.0;
}

int main(int argc, char* argv[])
{
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(16);
  assert(argc > 2);
  int M = atoi(argv[2]);
  srand(time(NULL));
  
  string prefix(argv[1]);
  MPOGen_Hubbard_BCS hgen((prefix + "/DMRG.in").c_str());
  cout << "\tConstructing MPOs" << endl; 
  MPO<Quantum> H = hgen.generate(); // matrix product operator
  int nsite = H.size();
  //
  // define working space
  //
  MpStorages sites(nsite);
  for (int i = 0; i < nsite; ++i) {
    sites[i].mpo = H[i];
    H[i].clear();
  }

  cout << "\tIntializing MPS" << endl;
  Qshapes<Quantum> qp = {Quantum(1), Quantum(-1), Quantum(0)};
  Dshapes dp = {1, 1, 2};
  physical(qp, dp);
  MPS<Quantum> A = create(nsite, Quantum::zero(), qp, dp, M, random_gen);
  init(sites, A, M);

  cout << "\tCalling DMRG program ( two-site algorithm) " << endl;

  double energy = 0.0;

  energy = dmrg(sites, TWOSITE, M);
  cout.precision(16);  
  cout << "\tGround state energy (two-site) = " << setw(20) << fixed << energy << endl << endl;

  cout << "\tCalling DMRG program ( one-site algorithm) " << endl;

  energy = dmrg(sites, ONESITE, M);
  cout.precision(16);  
  cout << "\tGround state energy (one-site) = " << setw(20) << fixed << energy << endl << endl;
  
  // for program reading
  cout << "!DMRG_Energy " << setw(18) << fixed << energy << endl;
  // 
  // convert final wavefunction into MPS form
  //
  for (int i = 0; i < nsite-1; ++i) {
    A[i] = sites[i].lmps;
  }
  // last site is special 
  auto qshape = sites[nsite-1].wfnc.qshape();
  qshape[2] = {sites[nsite-1].wfnc.q()};
  A[nsite-1].resize(Quantum::zero(), qshape, sites[nsite-1].wfnc.dshape(), false);
  for (auto it = sites[nsite-1].wfnc.begin(); it != sites[nsite-1].wfnc.end(); ++it) {
    A[nsite-1].insert(it->first, *(it->second));
  }

  ofstream frdm;
  frdm.open(prefix+"/1RDM.A");
  frdm.precision(16);  
  frdm << "!RDM.A" << setw(4) << nsite << endl << setw(20) << fixed << hgen.rdm(A, Spin::Up);
  frdm.close();
  frdm.open(prefix+"/1RDM.B");
  frdm << "!RDM.B" << setw(4) << nsite << endl << setw(20) << fixed << hgen.rdm(A, Spin::Down);
  frdm.close();
  frdm.open(prefix+"/KAPPA");
  frdm << "!KAPPA" << setw(4) << nsite << endl << setw(20) << fixed << hgen.kappa(A);
  frdm.close();
  return 0;
}
