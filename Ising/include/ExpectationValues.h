// Struct file
#pragma once

struct ExpectationValues {
  double expecE, expecM, expecMabs, expecChi, expecChiAbs, expecCv;
  arma::vec Evec;
  int accepted;
};
