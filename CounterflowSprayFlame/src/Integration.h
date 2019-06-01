#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <cmath>

double dtEff(const double dt, const double Beta);

double explicitDelta(const double phi, const double dtEff, const double Alpha, const double Beta);

double delta(const double phi, const double dt, const double Alpha, const double Beta);

#endif