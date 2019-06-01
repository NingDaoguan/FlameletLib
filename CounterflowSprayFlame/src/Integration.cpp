#include "Integration.h"


double dtEff(const double dt, const double Beta)
{
    return std::abs(Beta*dt) > 1e-10 ? (1 - std::exp(- Beta*dt))/Beta : dt;
}

double explicitDelta(const double phi, const double dtEff, const double Alpha, const double Beta)
{
    return (Alpha - Beta*phi)*dtEff;
}

double delta(const double phi, const double dt, const double Alpha, const double Beta)
{
    return explicitDelta(phi, dtEff(dt, Beta), Alpha, Beta);
}