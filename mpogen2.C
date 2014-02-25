#include "mpogen.h"
#include <assert.h>
#include <sstream>      // std::istringstream
#include <iostream>
#include <omp.h>
#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::getline;
using boost::is_any_of;
using boost::token_compress_on;
using boost::trim;
using std::array;

MPOGen_Hubbard_Explicit::MPOGen_Hubbard_Explicit(const char* m_input): MPOGen(m_input) {
  // open input file
  std::ifstream in(input.c_str());
  string line;
  getline(in, line);
  trim(line);
  vector<string> tok;
  boost::split(tok, line, is_any_of(", \t"), token_compress_on);
  nsite = atoi(tok[2].c_str());
  getline(in, line);
  getline(in, line);

  getline(in, line);
  trim(line);
  boost::split(tok, line, is_any_of(", \t"), token_compress_on);
  if (boost::iequals(tok[0], "IUHF")) {
    restricted = false;
    H0a.ReSize(nsite);
    H0a = 0;
    H0b.ReSize(nsite);
    H0b = 0;  
    D0.ReSize(nsite, nsite);
    D0 = 0;
    vccdd_aa.ReSize(nsite * nsite);
    vccdd_aa = 0;
    vccdd_bb.ReSize(nsite * nsite);
    vccdd_bb = 0;
    vccdd_ab.ReSize(nsite * nsite);
    vccdd_ab = 0;
    vcccd_a.ReSize(nsite*(nsite-1)/2, nsite*nsite);
    vcccd_a = 0;
    vcccd_b.ReSize(nsite*(nsite-1)/2, nsite*nsite);
    vcccd_b = 0;
    vcccc.ReSize(nsite*(nsite-1)/2, nsite*(nsite-1)/2);
    vcccc = 0;
    getline(in, line);
    trim(line);
  } else {
    restricted = true;
    H0a.ReSize(nsite);
    H0a = 0;
    D0.ReSize(nsite, nsite);
    D0 = 0;
    vccdd_aa.ReSize(nsite * nsite);
    vccdd_aa = 0;
    vcccd_a.ReSize(nsite*(nsite-1)/2, nsite*nsite);
    vcccd_a = 0;
    vcccc.ReSize(nsite*(nsite-1)/2, nsite*(nsite-1)/2);
    vcccc = 0;
  }

  boost::split(tok, line, is_any_of(", \t"), token_compress_on);
  assert(boost::iequals(tok[0], "&END"));

  cout << "Start reading integrals" << endl;

  getline(in, line);
  trim(line);  
  int section = 0;
  while(line.size() != 0) {
    boost::split(tok, line, is_any_of(" \t"), token_compress_on);
    assert(tok.size() == 5);
    double value = atof(tok[0].c_str());
    int i = atoi(tok[1].c_str());
    int j = atoi(tok[2].c_str());
    int k = atoi(tok[3].c_str());
    int l = atoi(tok[4].c_str());
    if (i==0 && j==0 && k==0 && l==0) {
      section += 1;
    } else if (restricted) {
      if (section == 0) { // ccdd
        vccdd_aa((i-1)*nsite+j, (k-1)*nsite+l) = value;
      } else if (section == 1) { // cccd
        assert(i>j);
        vcccd_a((i-1)*(i-2)/2+j, (k-1)*nsite+l) = value;
      } else if (section == 2) { // cccc
        assert(i>j && k>l && (i-1)*(i-2)/2+j >= (k-1)*(k-2)/2+l);
        vcccc((i-1)*(i-2)/2+j, (k-1)*(k-2)/2+l) = value;
        vcccc((k-1)*(k-2)/2+l, (i-1)*(i-2)/2+j) = value;        
      } else if (section == 3) { // cd
        assert(k==0 && l==0);
        H0a(i, j) = value;
      } else if (section == 4) { // cc
        assert(k ==0 && l ==0);
        assert(i >= j);
        D0(i, j) = value;
        D0(j, i) = value;
      } else {
        abort();
      }
    } else {
      if (section == 0) { // ccdd_aa
        vccdd_aa((i-1)*nsite+j, (k-1)*nsite+l) = value;        
      } else if (section == 1) { // ccdd_bb
        vccdd_bb((i-1)*nsite+j, (k-1)*nsite+l) = value;        
      } else if (section == 2) { // ccdd_ab
        vccdd_ab((i-1)*nsite+j, (k-1)*nsite+l) = value;        
      } else if (section == 3) { // cccd_a
        assert(i>j);
        vcccd_a((i-1)*(i-2)/2+j, (k-1)*nsite+l) = value;
      } else if (section == 4) { // cccd_b
        assert(i>j);
        vcccd_b((i-1)*(i-2)/2+j, (k-1)*nsite+l) = value;
      } else if (section == 5) { // cccc
        assert(i>j && k>l);
        vcccc((i-1)*(i-2)/2+j, (k-1)*(k-2)/2+l) = value;
      } else if (section == 6) { // cd_a
        assert(k==0 && l==0);
        H0a(i, j) = value;
      } else if (section == 7) { // cd_b
        assert(k==0 && l==0);
        H0b(i, j) = value;
      } else if (section == 8) { // cc
        assert(k==0 && l==0);
        D0(i, j) = value;
      } else {
        abort();
      }
    }
    getline(in, line);
    trim(line); 
  }
  in.close();
}

const MPO<Quantum> MPOGen_Hubbard_Explicit::generate(int M) {
  MPO<Quantum> H(nsite);
  zero(H);
  if (H0a.NormFrobenius() > 1e-10) {
    for (int i = 0; i < nsite; ++i) {
      ColumnVector temp(nsite);
      temp = 0.;
      temp(i+1) = 1.;
      axpy(1., create_op(temp, Spin::Up) * anni_op(H0a.Row(i+1).t(), Spin::Up), H);
      if (restricted) {
        axpy(1., create_op(temp, Spin::Down) * anni_op(H0a.Row(i+1).t(), Spin::Down), H);
      } else {
        axpy(1., create_op(temp, Spin::Down) * anni_op(H0b.Row(i+1).t(), Spin::Down), H);
      }
      compress(H, MPS_DIRECTION::Left, 0);
      compress(H, MPS_DIRECTION::Right, 0);
    }
  }
  // then Delta0 part
  if (D0.NormFrobenius() > 1e-10) {
    for (int i = 0; i < nsite; ++i) {
      ColumnVector temp(nsite);
      temp = 0.;
      temp(i+1) = 1.;
      axpy(1., create_op(temp, Spin::Up) * create_op(D0.Row(i+1).t(), Spin::Down), H);
      axpy(1., anni_op(temp, Spin::Down) * anni_op(D0.Column(i+1), Spin::Up), H);
      compress(H, MPS_DIRECTION::Left, 0);
      compress(H, MPS_DIRECTION::Right, 0);
    }
  }

  // build components
  map<array<int, 2>, MPO<Quantum>> caca, cacb, cbcb, dada, dbda, dbdb, cbda, cadb;
  ColumnVector temp_i(nsite);
  ColumnVector temp_j(nsite);  
  for (int i = 0; i < nsite; ++i) {
    temp_i = 0.;
    temp_i(i+1) = 1.;
    for (int j = 0; j < nsite; ++j) {
      array<int, 2> idx = {i, j};
      temp_j = 0.;
      temp_j(j+1) = 1.;
      if (i != j) {
        caca.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, create_op(temp_i, Spin::Up) * create_op(temp_j, Spin::Up)));
        cbcb.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, create_op(temp_i, Spin::Down) * create_op(temp_j, Spin::Down)));
        dada.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, anni_op(temp_i, Spin::Up) * anni_op(temp_j, Spin::Up)));
        dbdb.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, anni_op(temp_i, Spin::Down) * anni_op(temp_j, Spin::Down)));
        compress(caca[idx], MPS_DIRECTION::Left, 0);
        compress(caca[idx], MPS_DIRECTION::Right, 0);
        compress(cbcb[idx], MPS_DIRECTION::Left, 0);
        compress(cbcb[idx], MPS_DIRECTION::Right, 0);
        compress(dada[idx], MPS_DIRECTION::Left, 0);
        compress(dada[idx], MPS_DIRECTION::Right, 0);
        compress(dbdb[idx], MPS_DIRECTION::Left, 0);
        compress(dbdb[idx], MPS_DIRECTION::Right, 0);
      }
      cacb.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, create_op(temp_i, Spin::Up) * create_op(temp_j, Spin::Down)));
      dbda.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, anni_op(temp_i, Spin::Down) * anni_op(temp_j, Spin::Up)));
      cbda.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, create_op(temp_i, Spin::Down) * anni_op(temp_j, Spin::Up)));
      cadb.insert(std::pair<array<int, 2>, MPO<Quantum>>(idx, create_op(temp_i, Spin::Up) * anni_op(temp_j, Spin::Down)));
      compress(cacb[idx], MPS_DIRECTION::Left, 0);
      compress(cacb[idx], MPS_DIRECTION::Right, 0);
      compress(dbda[idx], MPS_DIRECTION::Left, 0);
      compress(dbda[idx], MPS_DIRECTION::Right, 0);
      compress(cbda[idx], MPS_DIRECTION::Left, 0);
      compress(cbda[idx], MPS_DIRECTION::Right, 0);
      compress(cadb[idx], MPS_DIRECTION::Left, 0);
      compress(cadb[idx], MPS_DIRECTION::Right, 0);
    }
  }

  // now interacting part
  for (int i = 0; i < nsite; ++i) {
    for (int j = 0; j < nsite; ++j) {
      array<int, 2> idx1 = {i, j};
      array<int, 2> ridx1 = {j, i};      
      for (int k = 0; k < nsite; ++k) {
        for (int l = 0; l < nsite; ++l) {
          array<int, 2> idx2 = {k, l};
          array<int, 2> ridx2 = {l, k};          
          if (restricted) {
            if (fabs(vccdd_aa(i*nsite+l+1, j*nsite+k+1)) > 1e-10) {
              if (i != j && k != l) {
                axpy(0.5*vccdd_aa(i*nsite+l+1, j*nsite+k+1), caca[idx1] * dada[idx2], H);
                axpy(0.5*vccdd_aa(i*nsite+l+1, j*nsite+k+1), cbcb[idx1] * dbdb[idx2], H);
              }
              axpy(vccdd_aa(i*nsite+l+1, j*nsite+k+1), cacb[idx1] * dbda[idx2], H);
            }
            if (i > j && fabs(vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1)) > 1e-10) {
              axpy(vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1), caca[idx1] * cbda[idx2], H);
              axpy(-vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1), cbcb[idx1] * cadb[idx2], H);
              axpy(vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1), cadb[ridx2] * dada[ridx1], H);
              axpy(-vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1), cbda[ridx2] * dbdb[ridx1], H);
            }
          } else {
            if (i !=j && k != l && fabs(vccdd_aa(i*nsite+l+1, j*nsite+k+1)) > 1e-10) {
              axpy(0.5*vccdd_aa(i*nsite+l+1, j*nsite+k+1), caca[idx1] * dada[idx2], H);
            }
            if (i != j && k != l && fabs(vccdd_bb(i*nsite+l+1, j*nsite+k+1)) > 1e-10) {
              axpy(0.5*vccdd_bb(i*nsite+l+1, j*nsite+k+1), cbcb[idx1] * dbdb[idx2], H);
            }
            if (fabs(vccdd_ab(i*nsite+l+1, j*nsite+k+1)) > 1e-10) {
              axpy(vccdd_ab(i*nsite+l+1, j*nsite+k+1), cacb[idx1] * dbda[idx2], H);
            }
            if (i > j && fabs(vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1)) > 1e-10) {
              axpy(vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1), caca[idx1] * cbda[idx2], H);
              axpy(vcccd_a(i*(i-1)/2+j+1, k*nsite+l+1), cadb[ridx2] * dada[ridx1], H);
            }
            if (i > j && fabs(vcccd_b(i*(i-1)/2+j+1, k*nsite+l+1)) > 1e-10) {
              axpy(vcccd_b(i*(i-1)/2+j+1, k*nsite+l+1), cbcb[idx1] * cadb[idx2], H);
              axpy(vcccd_b(i*(i-1)/2+j+1, k*nsite+l+1), cbda[ridx2] * dbdb[ridx1], H);
            }
          }
          if (i > j && l > k && fabs(vcccc(i*(i-1)/2+j+1, l*(l-1)/2+k+1)) > 1e-10) {
            axpy(vcccc(i*(i-1)/2+j+1, l*(l-1)/2+k+1), caca[idx1] * cbcb[idx2], H);
            axpy(vcccc(i*(i-1)/2+j+1, l*(l-1)/2+k+1), dbdb[ridx2] * dada[ridx1], H);
          }
        }
        compress(H, MPS_DIRECTION::Left, 0);
        compress(H, MPS_DIRECTION::Right, 0);
      }
    }
  }
  return std::move(H);
}
