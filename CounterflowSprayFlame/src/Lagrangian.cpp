#include "Lagrangian.h"
#include "Integration.h"
#include <iomanip>

namespace Cantera
{

Lagrangian::Lagrangian( const doublereal diameterInjection,
                        const doublereal mdotInjection,
                        const doublereal TInjection,
                        const doublereal p0 ) :
    diameterInjection_(diameterInjection),
    mdotInjection_(mdotInjection),
    TInjection_(TInjection),
    dt_(2e-5),
    small(1.0e-10),
    p0_(p0),
    flow_(0),
    loopcnt_(0)
{
    const doublereal d0 = diameterInjection_;
    const doublereal V =  Pi*d0*d0*d0/6.0;
    const doublereal dM = rhod(TInjection_) * V;
    const doublereal nPerSec_ = mdotInjection_ / dM;
    np_ = nPerSec_*dt_;
    std::cout << "#\tn particles per parcel\t:\t" << np_ << std::endl;
}


bool Lagrangian::setupFlowField(const vector_fp& solution)
{
    this->clear();
    z_ = flow_->grid();
    rho_ = flow_->density();
    mu_ = flow_->viscosity();
    cp_ = flow_->cp();
    std::cout << "rho: " << rho_[0] <<std::endl;

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
    std::cout << std::endl;

    this->evalRsd();
    return (rsd_ < 0.0001);
}


void Lagrangian::clear() {
    // clear flow field
    z_.clear();
    rho_.clear();
    mu_.clear();
    u_.clear();
    T_.clear();
    Y_.resize(fuelName_.size());
    for (size_t i=0; i<Y_.size(); i++) {
        Y_[i].clear();
    }
    // clear parcels
    position_.clear();
    velocity_.clear();
    diameter_.clear();
    temperature_.clear();
    hTrans_.clear();
    mTrans_.resize(fuelName_.size());
    for (size_t i=0; i<mTrans_.size(); i++) {
        mTrans_[i].clear();
    }

    // clear tracking data
    tp_.clear();
    tt_.clear();
    td_.clear();

    // clear transfer fields
    htf_.clear();
    mtf_.resize(fuelName_.size());
    for (size_t i=0; i<mtf_.size(); i++) {
        mtf_[i].clear();
    }
    tmtf_.clear();
}


void Lagrangian::relax() {
    vector_fp thtf(htf_);
    std::vector<std::vector<doublereal> > tmtf(mtf_);

    if (loopcnt_ == 0) {
        this->scale(htf_, 0.5*rlxf_);
        for (size_t i=0; i<fuelName_.size(); i++) {
            this->scale(mtf_[i], 1.3);
        }
        thtf = htf_;
        tmtf = mtf_;
    }
    else {
        this->relax(htf_, htfOld_);
        for (size_t i=0; i<fuelName_.size(); i++) {
            this->relax(mtf_[i], mtfOld_[i]);
        }
    }

    for (size_t i=0; i<fuelName_.size(); i++) {
        for (size_t j=0; j<tmtf_.size(); j++) {
            tmtf_[j] += mtf_[i][j];
        }
    }

    zOld_ = z_;
    htfOld_ = thtf;
    mtfOld_ = tmtf;
}


doublereal Lagrangian::linearInterpolate(const vector_fp& field, const doublereal z) const
{
    const size_t zSize = z_.size();
    if (field.size() != zSize) {
        std::cerr << "LINEAR INTERPOLATE ERROR: UNEQUAL SIZE" << std::endl;
    }
    size_t leftIndex=0;
    size_t rightIndex=0;
    doublereal leftLength=0.0;
    doublereal rightLength=0.0;
    doublereal sumLength=1.0;

    // check in case z is eactly the first/last point
    if (z == z_[0]) return field[0];
    else if (z == z_[zSize-1]) return field[zSize-1];
    else {
        for (size_t iz=0; iz<zSize; iz++){
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
        return field[leftIndex]*rightLength/sumLength
             + field[rightIndex]*leftLength/sumLength;
    }
}


doublereal Lagrangian::linearInterpolate(const vector_fp& field,
                                         const vector_fp& grid,
                                         const doublereal z) const
{
    const size_t zSize = grid.size();
    if (field.size()!=zSize) {
        std::cerr << "LINEAR INTERPOLATE ERROR: UNEQUAL SIZE" << std::endl;
    }
    size_t leftIndex=0;
    size_t rightIndex=0;
    doublereal leftLength=0.0;
    doublereal rightLength=0.0;
    doublereal sumLength=1.0;

    // check in case z is eactly the first/last point
    if (z == grid[0]) return field[0];
    else if (z == grid[zSize-1]) return field[zSize-1];
    else
    {
        for (size_t iz=0;iz<zSize;iz++)
        {
            if (z >= grid[iz]) continue;
            else
            {
                leftIndex = iz-1;
                rightIndex = iz;
                leftLength = z - grid[leftIndex];
                rightLength = grid[rightIndex] - z;
                sumLength = leftLength + rightLength;
                break;
            }
        }
        return field[leftIndex]*rightLength/sumLength
             + field[rightIndex]*leftLength/sumLength;
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
    doublereal rhoGas = linearInterpolate(rho_,oldP);
    rhoGas = rhoGas > small ? rhoGas : 1.1;
    doublereal uGas = linearInterpolate(u_, oldP);
    doublereal muGas = linearInterpolate(mu_, oldP);
    muGas = muGas > small ? muGas : 1.7e-5;
    doublereal cpGas = linearInterpolate(cp_, oldP);
    cpGas = cpGas > small ? cpGas : 1006.0;
    doublereal TGas = std::max( linearInterpolate(T_,oldP), 200.0 );
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
        const doublereal Xs = (101325.0/p0_)
                            * std::exp( (latentHeat(oldT)/(8.314/W_V) )
                            * (1.0/T_B - 1.0/oldT) );
        const doublereal Lk = muGas * std::sqrt(2.0*Pi*oldT*8.314/W_V) / (Sc*p0_);
        const doublereal Xsneq = Xs - 2.0*Lk*beta/oldD;
        doublereal Ysneq = Xsneq / ( Xsneq + (1.0 - Xsneq)*W_C/W_V );
        Ysneq = (Ysneq < 0.99999 ? Ysneq : 0.99999);
        doublereal B = (Ysneq - YGas[0]) / (1.0 - Ysneq); // single-component
        B = (B > -1.0 ? B : -0.99999);
        mdot = Sh/(3.0*Sc) * (md/taud) * std::log(1.0 + B);
        if (std::abs(mdot - mdotDash)/std::max(small, mdotDash) < 1.0e-4) break;
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
    }
    else {
        position_[ip] = -1.0;
        diameter_[ip] = 0.0;
        velocity_[ip] = 0.0;
    }
}


void Lagrangian::solve()
{
    for (int num=1; num<10000; num++) {
        this->inject();
        if (diameter_[0] > small) {
            this->track();
        }
        else break;
        for (int ip=0; ip<num; ip++) {
            if (position_[ip] > z_[z_.size()-1] || position_[ip] < z_[0]) {
                continue;
            }
            else {
                this->calcTrans(ip);
            }
        }
    }

    this->evalTrans();

    this->relax();

    loopcnt_++;
}


void Lagrangian::evalTrans()
{
    // resize to grid size
    htf_.resize(z_.size(), 0.0);
    for (size_t i=0; i<mtf_.size(); i++) {
        mtf_[i].resize(z_.size(), 0.0);
    }
    tmtf_.resize(z_.size(), 0.0);

    doublereal leftz;
    doublereal rightz;
    doublereal leftLength;
    doublereal rightLength;
    doublereal dz;
    doublereal sumH = 0.0;
    vector_fp sumM(fuelName_.size(), 0.0);

    // left-most point
    rightz = z_[1];
    // find
    sumH = 0.0;
    for (size_t is=0; is<sumM.size(); is++) sumM[is] = 0;
    for (size_t ip=0; ip<position_.size(); ip++) {
        dz = rightz - z_[0];
        if (position_[ip]>z_[0] && position_[ip]<=rightz) {
            rightLength = rightz - position_[ip];
            sumH += (rightLength/dz) * hTrans_[ip];
            for (size_t i=0;i<fuelName_.size();i++) {
                sumM[i] += (rightLength/dz) * mTrans_[i][ip];
            }
        }
    }
    dz = rightz - z_[0];
    htf_[0] = sumH/(0.5*dz);
    for (size_t i=0; i<mtf_.size(); i++) {
        mtf_[i][0] = sumM[i]/(0.5*dz);
    }

    // interior points
    for (size_t iz=1; iz<z_.size()-1; iz++) {
        leftz = z_[iz-1];
        rightz = z_[iz+1];
        // find
        sumH = 0.0;
        for (size_t is=0; is<sumM.size(); is++) sumM[is] = 0;
        for (size_t ip=0; ip<position_.size(); ip++) {
            dz = z_[iz] - leftz;
            if (position_[ip]>leftz && position_[ip]<=z_[iz]) {
                leftLength = position_[ip] - leftz;
                sumH += (leftLength/dz) * hTrans_[ip];
                for (size_t i=0; i<fuelName_.size(); i++) {
                    sumM[i] += (leftLength/dz) * mTrans_[i][ip];
                }
            }

            dz = rightz - z_[iz];
            if (position_[ip]>z_[iz] && position_[ip]<=rightz) {
                rightLength = rightz - position_[ip];
                sumH += (rightLength/dz) * hTrans_[ip];
                for (size_t i=0;i<fuelName_.size();i++) {
                    sumM[i] += (rightLength/dz) * mTrans_[i][ip];
                }
            }
        }
        dz = rightz - leftz;
        htf_[iz] = sumH/(csArea(iz)*0.5*dz);
        for (size_t i=0; i<mtf_.size(); i++) {
            mtf_[i][iz] = sumM[i]/(csArea(iz)*0.5*dz);
        }
    }

    // right-most point
    leftz = z_[z_.size()-2];
    // find
    sumH = 0.0;
    for (size_t is=0; is<sumM.size(); is++) sumM[is] = 0;
    for (size_t ip=0; ip<position_.size(); ip++) {
        dz = z_[z_.size()-1] - leftz;
        if (position_[ip]>leftz && position_[ip]<=z_[z_.size()-1]) {
            leftLength = position_[ip] - leftz;
            sumH += (leftLength/dz) * hTrans_[ip];
            for (size_t i=0;i<fuelName_.size();i++) {
                sumM[i] += (leftLength/dz) * mTrans_[i][ip];
            }
        }
    }
    dz = z_[z_.size()-1] - leftz;
    htf_[z_.size()-1] = sumH/(0.5*dz);
    for (size_t i=0; i<mtf_.size(); i++) {
        mtf_[i][z_.size()-1] = sumM[i]/(0.5*dz);
    }

    // multiplied by np_
    for (auto it=htf_.begin(); it!=htf_.end(); it++) {
        *it *= np_;
    }
    for (auto it1=mtf_.begin(); it1!=mtf_.end(); it1++) {
        for (auto it2=(*it1).begin(); it2!=(*it1).end(); it2++) {
            *it2 *= np_;
        }
    }
}


void Lagrangian::evalRsd()
{
    if (loopcnt_==0) {
        rsd_ = 1.0;
    }
    else {
        doublereal sum1 = 0.0;
        doublereal sum2 = 0.0;
        for (size_t j=0; j<T_.size(); j++) {
            sum1 += ( T_[j] - linearInterpolate(TOld_, zOld_, z_[j]) )
                  * ( T_[j] - linearInterpolate(TOld_, zOld_, z_[j]) );
            sum2 += linearInterpolate(TOld_, zOld_, z_[j])
                  * linearInterpolate(TOld_, zOld_, z_[j]);
        }
        sum1 /= z_.size()*fuelName_.size();
        sum2 /= z_.size()*fuelName_.size();
        rsd_ = std::sqrt(sum1) / std::sqrt(sum2 + small);
    }
    std::cout << "T RMS Residual\t:\t" << rsd_ << std::endl;

    TOld_ = T_;
}


doublereal Lagrangian::csArea(size_t iz) const {
    doublereal rho1;
    doublereal u1;
    doublereal rho2 = rho_[iz] > small ? rho_[iz] : 1.1;
    doublereal u2 = u_[iz];
    u2 = (std::abs(u2) > small ? u2 : small);

    if (u2 >=0 ) {
        rho1 = rho_[0] > small ? rho_[0] : 1.1;
        u1 = u_[0];
    }
    else {
        rho1 = rho_[z_.size()-1] > small ? rho_[z_.size()-1] : 1.1;
        u1 = u_[z_.size()-1];
    }
    return (rho1*u1) / (rho2*u2);
}


void Lagrangian::write() const
{
    std::ofstream fout1(outfile_);
    fout1 << "# t (s), z (m), T (K), D (m)" << std::endl;
    for (size_t it=0; it<tp_.size(); it++) {
        fout1 << it*dt_ << ","
              << tp_[it] << ","
              << tt_[it] << ","
              << td_[it] << std::endl;
    }

    std::ofstream fout2("trans" + outfile_);
    fout2 << "# z (m), hTrans (J/m^3/s), mTrans (kg/m^3/s)" << std::endl;
    for (size_t iz=0; iz<z_.size(); iz++) {
        fout2 << z_[iz] << ","
              << htf_[iz] << ","
              << mtf_[0][iz] << std::endl;
    }
}




}