#include "Lagrangian.h"
#include <iomanip>

namespace Cantera
{

Lagrangian::Lagrangian( const doublereal diameterInjection,
                        const doublereal mdotInjection,
                        const doublereal TInjection ) :
    diameterInjection_(diameterInjection),
    mdotInjection_(mdotInjection),
    TInjection_(TInjection),
    dt_(2e-5),
    delta_(1.0e-10),
    flow_(0),
    loopcnt_(0)
{

}


void Lagrangian::setupFlowField(const vector_fp& solution)
{
    this->resize();
    z_ = flow_->grid();
    rho_ = flow_->density();
    mu_ = flow_->viscosity();
    cp_ = flow_->cp();

    std::vector<size_t> fuelIndex(fuelName_.size());
    for (size_t ik=0;ik<fuelName_.size();ik++) {
        fuelIndex[ik] = flow_->componentIndex(fuelName_[ik]);
    }

    int cnt = 0;
    for (size_t i=0; i<solution.size(); i++) {
        if (i%(flow_->nsp()+c_offset_Y) == 0) {
            u_.push_back(solution[i]);
            cnt++;
        }
        else if (i%(flow_->nsp()+c_offset_Y) == 2) {
            T_.push_back(solution[i]);
            std::cout << std::setw(10) << std::setprecision(5) << T_[cnt-1];
            if (cnt%6==0) std::cout << std::endl;
        }
        else {
            for (size_t ik=0; ik<fuelName_.size(); ik++) {
                if (i%(flow_->nsp()+c_offset_Y) == fuelIndex[ik]) {
                    Y_[ik].push_back(solution[i]);
                }
            }
        }
    }
}


doublereal Lagrangian::linearInterpolate(const vector_fp& field, const doublereal z) const
{
    const size_t zSize = z_.size();
    if (field.size() != zSize) std::cerr << "LINEAR INTERPOLATE ERROR: UNEQUAL SIZE" << std::endl;
    size_t leftIndex=0;
    size_t rightIndex=0;
    doublereal leftLength=0.0;
    doublereal rightLength=0.0;
    doublereal sumLength=1.0;

    // check in case z is eactly the first/last point
    if (z == z_[0]) return field[0];
    else if (z == z_[zSize-1]) return field[zSize-1];
    else
    {
        for (int iz=0; iz<zSize; iz++)
        {
            if (z >= z_[iz]) continue;
            else
            {
                leftIndex = iz-1;
                rightIndex = iz;
                leftLength = z - z_[leftIndex];
                rightLength = z_[rightIndex] - z;
                sumLength = leftLength + rightLength;
                break;
            }
        }
        return field[leftIndex]*rightLength/sumLength + field[rightIndex]*leftLength/sumLength;
    }
}


void Lagrangian::calchmTrans(int ip)
{
    const doublereal oldP = position_[ip];
    const doublereal oldV = velocity_[ip];
    const doublereal oldD = diameter_[ip];
    const doublereal oldT = temperature_[ip];

    doublereal& newD = diameter_[ip];
    doublereal& newT = temperature_[ip];
    doublereal& newhTrans = hTrans_[ip];
    doublereal& newmTrans = mTrans_[ip];

    // gas-phase
    doublereal rhoGas = linearInterpolate(rho_,oldP);
    rhoGas = (rhoGas > delta_ ? rhoGas : 1.1);
    doublereal uGas = linearInterpolate(u_, oldP);
    doublereal muGas = linearInterpolate(mu_, oldP);
    muGas = (muGas > delta_ ? cpGas : 1.7e-5);
    doublereal cpGas = linearInterpolate(cp_, oldP);
    cpGas = (cpGas > delta_ ? cpGas : 1006.0);
    doublereal TGas = std::max( linearInterpolate(T_,oldP), 273.15 );
    vector_fp YGas(fuelName_.size(), 0.0);
    for (size_t i=0; i<fuelName_.size(); i++) {
        YGas[i] = linearInterpolate(Y_[i], oldP);   
    }

    // Non-equilibrium Langmuir-Knudsen evaporation law
    doublereal relativeRe = rhoGas * std::abs(uGas - oldV) * oldD / muGas;
    relativeRe = (relativeRe > delta_ ? relativeRe : delta_);

    const doublereal Nu = 2.0+0.552*std::pow(relativeRe, 0.5)*std::pow(Pr,0.333333); // Nusselt number
    const doublereal Sh = 2.0+0.552*std::pow(relativeRe, 0.5)*std::pow(Sc,0.333333); // Sherwood number


}


bool Lagrangian::solve()
{
    for (int num=1; num<20000; num++) {
        this->inject();

        for (int ip=0; ip<num; ip++) {
            const doublereal oldP = position_[ip];
            const doublereal oldV = velocity_[ip];
            const doublereal oldD = diameter_[ip];
            const doublereal oldT = temperature_[ip];
            doublereal& newP = position_[ip];
            doublereal& newV = velocity_[ip];
            doublereal& newD = diameter_[ip];
            doublereal& newT = temperature_[ip];

            // if particles leak out of the domain
            if (oldP > z_[z_.size()-1] || oldP < z_[0]) {
                continue;
            }
            else {
                this->calchmTrans(ip);
            }
        }

    }

    return false;
}


void Lagrangian::write() const
{
    std::ofstream fout(outfile_);
    fout<< "# Particle Tracking" << std::endl;
    fout<< "# Time (s), Position (m), Temperature (K), Diameter (m)" << std::endl;
    for (size_t it=0; it<position_.size(); it++)
    {
        fout<< it*dt_ << ","
            << position_[it] << ","
            << temperature_[it] << ","
            << diameter_[it] << std::endl;
    }
}






}