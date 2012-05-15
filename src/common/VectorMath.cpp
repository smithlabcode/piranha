
#include <iostream>
#include <assert.h>
#include <cmath>
#include "VectorMath.hpp"


using std::vector;
using std::cout;
using std::endl;

double
vectorDotProduct(const vector<double> a, const vector<double> b) {
  double res = 0;
  //assert(a.size() == b.size());
  for (size_t i=0; i<a.size(); i++) {
    res += a[i] * b[i]; 
  }
  return res;
}

vector<double>
vectorScalarProduct(const double s, const vector<double> v) {
  vector<double> res;
  for (size_t i=0; i<v.size(); i++) res.push_back(v[i] * s);
  return res;
}


vector<double>
vectorSubtraction(const vector<double> a, const vector<double> b) {
  vector<double> res;
  //assert(a.size() == b.size());
  for (size_t i=0; i<a.size(); i++) res.push_back(a[i] - b[i]);
  return res;
}


bool 
vectorUnderThreshold(const vector<double> v, double t) {
  for (size_t i=0; i<v.size(); i++) {
    if (v[i] > t) return false;
  }
  return true;
}


bool 
vectorUnderThresholdAbs(const vector<double> v, double t) {
  for (size_t i=0; i<v.size(); i++) {
    if (fabs(v[i]) > t) return false;
  }
  return true;
}

