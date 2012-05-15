#ifndef VMATH_HPP
#define VMATH_HPP

#include <vector>

double vectorDotProduct(const std::vector<double> a, 
                        const std::vector<double> b);
std::vector<double> vectorScalarProduct(const double s, 
                                        const std::vector<double> v);
std::vector<double> vectorSubtraction(const std::vector<double> a, 
                                      const std::vector<double> b);
bool vectorUnderThreshold(const std::vector<double> v, double t);
bool vectorUnderThresholdAbs(const std::vector<double> v, double t);


#endif

