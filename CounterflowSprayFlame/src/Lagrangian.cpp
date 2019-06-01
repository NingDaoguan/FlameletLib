#include "Lagrangian.h"
#include "Integration.h"
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
    small(1.0e-10),
    p0_(101325.0),
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
    else {
        for (int iz=0; iz<zSize; iz++){
            if (z >= z_[iz]) continue;
            else {
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


void Lagrangian::calcTrans(int ip)
{
    const doublereal oldP = position_[ip];
    const doublereal oldV = velocity_[ip];
    const doublereal oldD = diameter_[ip];
    const doublereal oldT = temperature_[ip];
    const doublereal md = Pi*oldD*oldD*oldD/6.0 * rhod(oldT);

    // gas-phase
    doublereal rhoGas = linearInterpolate(rho_,oldP) > small ? rhoGas : 1.1;
    doublereal uGas = linearInterpolate(u_, oldP);
    doublereal muGas = linearInterpolate(mu_, oldP) > small ? muGas : 1.7e-5;
    doublereal cpGas = linearInterpolate(cp_, oldP) > small ? cpGas : 1006.0;
    doublereal TGas = std::max( linearInterpolate(T_,oldP), 273.15 );
    vector_fp YGas(fuelName_.size(), 0.0);
    for (size_t i=0; i<fuelName_.size(); i++) {
        YGas[i] = linearInterpolate(Y_[i], oldP);   
    }
    const doublereal Pr = 0.7;
    const doublereal Sc = 0.7;
    const doublereal D = muGas/(rhoGas*Sc); // diffusion coefficient
    const doublereal W_C = 28.97e-3; // molecular weight of air 

    // Non-equilibrium Langmuir-Knudsen evaporation law
    doublereal relativeRe = rhoGas * std::abs(uGas - oldV) * oldD / muGas;
    relativeRe = (relativeRe > small ? relativeRe : small);
    // Ranz-Marshall
    const doublereal Nu = 2.0+0.552*std::pow(relativeRe, 0.5)*std::pow(Pr,0.333333);
    const doublereal Sh = 2.0+0.552*std::pow(relativeRe, 0.5)*std::pow(Sc,0.333333);
    // time constant
    const doublereal taud = rhod(oldT)*(oldD + small)*(oldD + small)/(18.0*muGas);
    // iterative method to compute evaporation rate
    doublereal beta = 0.0;
    doublereal mdot = 0.0;
    for (int j=0; j<50; j++) {
        doublereal mdotDash = mdot;
        beta = 1.5*Pr*taud * (mdot/md);
        const doublereal Xs = (101325.0/p0_) * std::exp( (latentHeat(oldT)/(8.314/W_V) )
                            * (1/T_B - 1/oldT) );
        const doublereal Lk = muGas * std::sqrt(2.0*Pi*oldT*8.314/W_V) / (Sc*p0_);
        const doublereal Xsneq = Xs - 2.0*Lk*beta/oldD;
        doublereal Ysneq = Xsneq / ( Xsneq + (1.0 - Xsneq)*W_C/W_V );
        Ysneq = (Ysneq < 0.99999 ? Ysneq : 0.99999);
        doublereal B = (Ysneq - YGas[0]) / (1.0 - Ysneq); // single-component
        B = (B > 0 ? B : 0);
        mdot = Sh/(3.0*Sc) * (md/taud) * std::log(1.0 + B);
        if (std::abs(mdot - mdotDash)/std::max(small, mdotDash) < 1e-4) break;
    }
    beta = beta > small ? beta : small;
    doublereal mtc = (( beta / (std::exp(beta) - 1.0) ) / taud) * (cpGas / cpd(oldT));

    // update mass, diameter and the mass transfer
    doublereal mdNew = md - mdot*dt_;
    diameter_[ip] = std::pow(mdNew/rhod(oldT)*6.0/Pi, 1.0/3.0);
    mTrans_[0][ip] = mdot; // single-component

    // perform integration
    const doublereal bcp = Nu/(3.0*Pr)*mtc;
    const doublereal acp = bcp*TGas;
    const doublereal ancp = -latentHeat(oldT)*mdot/(md*cpd(oldT));

    const doublereal deltaT = delta(oldT, dt_, acp + ancp, bcp);
    const doublereal deltaTncp = ancp*dt_;
    const doublereal deltaTcp = deltaT - deltaTncp;

    // Calculate the new temperature and the enthalpy transfer
    doublereal Tnew = oldT + deltaT;
    Tnew = std::min(Tnew, T_B);
    temperature_[ip] = Tnew;
    hTrans_[ip] = md*cpd(oldT)*deltaTcp / dt_;

    if (diameter_[ip] > small) {
        // calculate force
        const doublereal dm = 0.5*(oldD + diameter_[ip]);
        const doublereal mm = 0.5*(md + mdNew);
        doublereal Red = rhoGas * std::abs(uGas - oldV) * dm / muGas;
        Red = (Red > small ? Red : small);
        doublereal Cd = 0.0;
        // Wen-Yu drag model--------------------------------------------------+
        // F = 0.5*Cd*rhoGas*(uGas-uDrop)*|uGas-uDrop|*Ap                     |
        // { Cd = 24.0/Red * (1.0+0.15*pow(Red,0.687)) # Re<1000              |
        // { Cd = 0.44 # Re>=1000                                             |
        // End----------------------------------------------------------------+
        if (Red<1000) {
            Cd = 24.0/Red;
            Cd *= (1.0+0.15*std::pow(Red,0.687));
        }
        else {
            Cd = 0.44;
        }
        // Second-order Runge-Kutta method
        doublereal F = 0.5*Cd*rhoGas*(uGas-oldV)*std::abs(uGas-oldV)*Pi*dm*dm/4.0;
        velocity_[ip] += F*dt_/mm;
        position_[ip] = 0.5*(oldV + velocity_[ip])*dt_ + oldP;
        uTrans_[ip] = -F;
    }
    else {
        diameter_[ip] = 0.0;
        velocity_[ip] = 0.0;
        uTrans_[ip] = 0.0;
    }
}


bool Lagrangian::solve()
{
    for (int num=1; num<10000; num++) {
        this->inject();
        if (diameter_[0] > small) {
            tp_.push_back(position_[0]);
            td_.push_back(diameter_[0]);
            tt_.push_back(temperature_[0]);
        }
        else break;
        for (int ip=0; ip<num; ip++) {
            // particle out of the domain
            if (position_[ip] > z_[z_.size()-1] || position_[ip] < z_[0]) {
                continue;
            }
            else {
                this->calcTrans(ip);
            }
        }
    }



    return false;
}


void Lagrangian::write() const
{
    std::ofstream fout(outfile_);
    fout<< "# Time (s), Position (m), Temperature (K), Diameter (m)" << std::endl;
    for (size_t it=0; it<tp_.size(); it++)
    {
        fout<< it*dt_ << ","
            << tp_[it] << ","
            << tt_[it] << ","
            << td_[it] << std::endl;
    }
}






}