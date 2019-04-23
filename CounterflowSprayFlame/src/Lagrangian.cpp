#include "Lagrangian.h"

namespace Cantera
{
Lagrangian::Lagrangian
(
    std::ofstream& LagrangianOutfile,
    size_t ignition,
    doublereal diameterInjection,
    doublereal mdotInjection
)
:
    ignition_(ignition),
    diameterInjection_(diameterInjection),
    mdotInjection_(mdotInjection),
    LagrangianOutfile_(LagrangianOutfile),
    dt_(2e-5),
    rhoDrop_(7.386e2),
    parcelNumber_(0),
    delta_(1.0e-12),
    loopCnt_(0),
    ok_(false)
{
    volumeDrop_ =  (4.0/3.0) * Pi * diameterInjection_ * diameterInjection_ * diameterInjection_/8.0;
    massDrop_ = rhoDrop_ * volumeDrop_;
    nPerSec_ = mdotInjection_ / massDrop_;
    nPerParcel_ = nPerSec_*dt_;
    std::cout << "\nNumber of particles per parcel: #\t" << nPerParcel_ << std::endl;
}


doublereal Lagrangian::linearInterpolate(const vector_fp& field, const doublereal z) const
{
    const size_t zSize = zField_.size();
    if (field.size()!=zSize) std::cerr << "LINEAR INTERPOLATE ERROR: UNEQUAL SIZE" << std::endl;
    size_t leftIndex=0;
    size_t rightIndex=0;
    doublereal leftLength=0.0;
    doublereal rightLength=0.0;
    doublereal sumLength=1.0;

    // check in case z is eactly the first/last point
    if (z == zField_[0]) return field[0];
    else if (z == zField_[zSize-1]) return field[zSize-1];
    else
    {
        for (int iz=0;iz<zSize;iz++)
        {
            if (z >= zField_[iz]) continue;
            else
            {
                leftIndex = iz-1;
                rightIndex = iz;
                leftLength = z - zField_[leftIndex];
                rightLength = zField_[rightIndex] - z;
                sumLength = leftLength + rightLength;
                break;
            }
        }
        return field[leftIndex]*rightLength/sumLength + field[rightIndex]*leftLength/sumLength;
    }
}


doublereal Lagrangian::linearInterpolate(const vector_fp& field, const vector_fp grid, const doublereal z) const
{
    const size_t zSize = grid.size();
    if (field.size()!=zSize) std::cerr << "LINEAR INTERPOLATE ERROR: UNEQUAL SIZE" << std::endl;
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
        for (int iz=0;iz<zSize;iz++)
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
        return field[leftIndex]*rightLength/sumLength + field[rightIndex]*leftLength/sumLength;
    }
}


void Lagrangian::evaluateTransferRateField()
{
    const size_t zSize = zField_.size();
    doublereal dz;
    doublereal sumLiquid;
    doublereal rhoGas;
    doublereal leftLength;
    doublereal rightLength;
    // left-most point
    rhoGas = rhoField_[0];
    if (rhoGas<delta_) rhoGas=0.8;
    sumLiquid=0.0;
    for (int iParcel=0;iParcel<parcelNumber_;iParcel++)
    {
        if (positions_[iParcel]>zField_[0] && positions_[iParcel] <= zField_[1])
        {
            rightLength = zField_[1] - positions_[iParcel];
            dz = zField_[1] - zField_[0];
            sumLiquid += (rightLength/dz) * (4.0*Pi/3.0) * particleDiameter_[iParcel]
                                                            * particleDiameter_[iParcel]
                                                            * particleDiameter_[iParcel]/8.0;
        }
    }
    liquidGasMassRatio_.push_back( (sumLiquid*nPerParcel_) / (0.5*dz) );
    heatTransferRateField_.push_back( 0 );
    totalMassTransferField_.push_back( 0 );
    for (int i=0;i<fuelNum_;i++)
        { massTransferRateField_[i].push_back( 0 ); }
    // interior points
    for (int iz=1;iz<zSize-1;iz++)
    {
        doublereal leftz = zField_[iz-1];
        doublereal rightz = zField_[iz+1];
        rhoGas = rhoField_[iz];
        if (rhoGas<delta_) rhoGas=0.8;
        // find
        sumLiquid = 0.0;
        doublereal sumQg = 0.0;
        doublereal sumEvap = 0.0;
        vector_fp sumJg(fuelNum_,0.0);

        for (int iParcel=0;iParcel<parcelNumber_;iParcel++)
        {
            if (positions_[iParcel]>leftz && positions_[iParcel]<=zField_[iz])
            {
                leftLength = positions_[iParcel] - leftz;
                dz = zField_[iz] - leftz;
                sumLiquid += (leftLength/dz) * (4.0*Pi/3.0) * particleDiameter_[iParcel]
                                                                * particleDiameter_[iParcel]
                                                                * particleDiameter_[iParcel]/8.0;
                sumQg += (leftLength/dz) * particleQg_[iParcel];
                for (int i=0;i<fuelNum_;i++)
                    { sumJg[i] += (leftLength/dz) * particleJg_[i][iParcel]; }
            }

            if (positions_[iParcel]>zField_[iz] && positions_[iParcel]<=rightz)
            {
                rightLength = rightz - positions_[iParcel];
                dz = rightz - zField_[iz];
                sumLiquid += (rightLength/dz) * (4.0*Pi/3.0) * particleDiameter_[iParcel]
                                                                * particleDiameter_[iParcel]
                                                                * particleDiameter_[iParcel]/8.0;
                sumQg += (rightLength/dz) * particleQg_[iParcel];
                for (int i=0;i<fuelNum_;i++)
                    { sumJg[i] += (rightLength/dz) * particleJg_[i][iParcel]; }
            }
        }
        dz = rightz - leftz;
        liquidGasMassRatio_.push_back( (sumLiquid*nPerParcel_) / (0.5*dz) );
        heatTransferRateField_.push_back( (sumQg*nPerParcel_) / (0.5*dz) );
        for (int i=0;i<fuelNum_;i++)
        {
            massTransferRateField_[i].push_back( (sumJg[i]*nPerParcel_) / (0.5*dz) );
            sumEvap += (sumJg[i]*nPerParcel_) / (0.5*dz);
        }
        totalMassTransferField_.push_back(sumEvap);
    }
    liquidGasMassRatio_.push_back( 0 );
    heatTransferRateField_.push_back( 0 );
    totalMassTransferField_.push_back( 0 );
    for (int i=0;i<fuelNum_;i++)
        { massTransferRateField_[i].push_back( 0 ); }

    doublereal residual = 0.0;
    // doublereal residualNorm = 0.0;
    if (loopCnt_!=0)
    {
        residual = this->evaluateResidual();
        // if (loopCnt_==1) residual0_ = residual;
        // residualNorm = residual/(residual0_+delta_);
        std::cout << "RMS residual:\t" << residual << std::endl;
    }
    
    if (loopCnt_ == 0 || residual > 0.01)
    {
        oldHeatTransferField_ = heatTransferRateField_;
        oldMassTransferField_ = massTransferRateField_;
        oldGrid_ = zField_;
    }
    else ok_ = true;

    if (loopCnt_ > 30) ok_ = true;
    loopCnt_++;
}


void Lagrangian::heatAndMassTransfer(const int iParcel)
{
    const doublereal zDrop = positions_[iParcel];
    const doublereal oldTemperature = particleTemperature_[iParcel];
    const doublereal oldDiameter = particleDiameter_[iParcel];
    const doublereal uDrop = particleVelocity_[iParcel];

    doublereal& newDiameter = particleDiameter_[iParcel];
    doublereal& newTemperature = particleTemperature_[iParcel];
    doublereal& newQg = particleQg_[iParcel];

    // gas phase information
    doublereal rhoGas = this->linearInterpolate(rhoField_,zDrop);
    if (rhoGas<delta_) rhoGas=0.8;
    doublereal uGas = this->linearInterpolate(uField_, zDrop);
    doublereal muGas = this->linearInterpolate(muField_, zDrop);
    if (muGas<delta_) muGas=1e-5;
    doublereal k = this->linearInterpolate(conField_, zDrop); // W/m*K
    k = (k > delta_ ? k : 0.0266);
    doublereal cpGas = this->linearInterpolate(cpField_, zDrop);
    cpGas = (cpGas > delta_ ? cpGas : 1006.0);
    doublereal TGas = std::max( this->linearInterpolate(TField_,zDrop), 273.15 );
    vector_fp YGas(fuelNum_,0);
    for (size_t i=0;i<fuelNum_;i++)
    {
        YGas[i] = this->linearInterpolate(YField_[i],zDrop);   
    }

    // ********************************************************************* //
    // thermophysical properties of Jet-A surrogate hard-coded
    // Non-equilibrium Langmuir-Knudsen evaporation model hard-coded
    doublereal relativeRe = rhoGas * std::abs(uGas - uDrop) * oldDiameter / muGas;
    relativeRe = (relativeRe > delta_ ? relativeRe : delta_);
    const doublereal mass = rhoDrop_ * (4.0/3.0)*Pi*oldDiameter*oldDiameter*oldDiameter/8.0;
    const doublereal W_C = 28.97e-3; // molecular weight of air 
    doublereal latentHeat = ( getLatentHeat(oldTemperature) > 1.0 ? getLatentHeat(oldTemperature) : 1.0 );
    const doublereal cpDrop = 3.84 * oldTemperature + 1056.88; // c_p droplet
    const doublereal Pr = 0.7;
    const doublereal Sc = 0.7;
    const doublereal D = muGas/(rhoGas*Sc); // diffusion coefficient
    doublereal beta = 0.0;
    doublereal JgTotal = 0.0;
    doublereal f2 = 0.0; // heat transfer correction fractor
    for (size_t i=0;i<fuelNum_;i++)
    {
        JgTotal += particleJg_[i][iParcel]; // total kg/s
    }
    beta = (rhoDrop_*Pr/(8*muGas)) * (4*JgTotal/(Pi*oldDiameter*rhoDrop_));
    if(beta<delta_) beta += delta_;
    f2 = beta / (std::exp(beta)-1);
    doublereal taud = rhoDrop_*oldDiameter*oldDiameter/(18.0*muGas); // response time
    vector_fp Xseq(fuelNum_,0); // equilibrium mole fractions
    vector_fp Xsneq(fuelNum_,0); // non-equilibrium mole fractions
    vector_fp Yseq(fuelNum_,0); // equilibrium mass fractions
    vector_fp Ysneq(fuelNum_,0); // non-equilibrium mass fractions
    vector_fp BH(4,0.0); // Spalding transfer numbers
    doublereal L_K; // Knudsen layer thickness
    for (size_t i=0;i<fuelNum_;i++)
    {
        Xseq[i] = (101325.0/p0_) * std::exp( (latentHeat/(8.314/W_V[i])) * (1/T_B[i] - 1/oldTemperature) );
        Yseq[i] = Xseq[i] / ( Xseq[i] + (1.0 - Xseq[i])*W_C/W_V[i] );
        L_K = (muGas*std::pow(2*Pi*oldTemperature*8.314/W_V[i],0.5)) / (Sc*p0_);
        Xsneq[i] = Xseq[i] - (2*L_K/oldDiameter)*beta;
        Ysneq[i] = Xsneq[i] / ( Xsneq[i] + (1.0 - Xsneq[i])*W_C/W_V[i] );
        Ysneq[i] = (Ysneq[i] < 0.99999 ? Ysneq[i] : 0.99999);
        BH[i] = (Ysneq[i] - YGas[i]) / (1.0 - Ysneq[i]);
        BH[i] = (BH[i] > 0 ? BH[i] : 0);
    }

    const doublereal Nu = 2.0+0.552*std::pow(relativeRe, 0.5)*std::pow(Pr,0.333333); // Nusselt number
    const doublereal Sh = 2.0+0.552*std::pow(relativeRe, 0.5)*std::pow(Sc,0.333333); // Sherwood number
    // ********************************************************************* //

    doublereal Qg = 0.0; // K/s
    doublereal sumMass = 0.0;
    vector_fp Jg(fuelNum_,0.0); // kg/s
    vector_fp massComp(fuelNum_,0.0);
    for (int i=0;i<fuelNum_;i++)
    {
        Jg[i] = Sh/(3*Sc) * (component_[i][iParcel]*mass/taud) * std::log(1.0+BH[i]); // kg/s
        if ( ( component_[i][iParcel]*mass - Jg[i]*dt_ ) <= 0.0 )
        {
            Jg[i] = component_[i][iParcel]*mass / dt_;
            massComp[i]=0.0;
        }
        else
        {
            massComp[i] = component_[i][iParcel]*mass - Jg[i]*dt_;
        }
        sumMass += massComp[i];
    }

    if (sumMass>delta_*delta_*delta_)
    {
        for (int i=0;i<fuelNum_;i++) { component_[i][iParcel] = massComp[i]/sumMass;}
        newDiameter = std::pow( (8.0*sumMass*3.0/(4.0*rhoDrop_*Pi)) , 1.0/3.0);
    }
    else newDiameter=0.0;

    // Qg (K/s)   receive from gas
    Qg = Nu/(3*Pr) * (cpGas/cpDrop) * (f2/taud) * (TGas - oldTemperature);
    if (newDiameter!=0.0)
    {
        doublereal sumJgLat = 0.0;
        for (int i=0;i<fuelNum_;i++)
            { sumJgLat += Jg[i]*latentHeat; }
        newTemperature = oldTemperature + dt_*(-sumJgLat)/(mass*cpDrop) + dt_*Qg;
    }
    else {newTemperature = oldTemperature;}

    // newQg = Qg*mass*cpDrop;
    newQg = cpDrop*mass*(newTemperature - oldTemperature) / dt_;
    for (int i=0;i<fuelNum_;i++)
    {
        particleJg_[i][iParcel] = Jg[i];
    }
}


doublereal Lagrangian::force(const int iParcel) const
{
    const doublereal zDrop = positions_[iParcel];
    const doublereal d0 = particleDiameter_[iParcel];
    const doublereal uDrop = particleVelocity_[iParcel];

    doublereal rhoGas = this->linearInterpolate(rhoField_,zDrop);
    if (rhoGas<delta_) rhoGas=0.8;
    doublereal uGas = this->linearInterpolate(uField_,zDrop);
    doublereal muGas = this->linearInterpolate(muField_,zDrop);
    if (muGas<delta_) muGas=1e-5;
    doublereal Red = rhoGas * std::abs(uGas - uDrop) * d0 / muGas;
    if (Red<delta_) return 0.0;

/*
    const doublereal mass = rhoDrop_ * (4.0/3.0)*Pi*d0*d0*d0/8.0;
    const doublereal taud = rhoDrop_*d0*d0/(18.0*muGas); // response time
    const doublereal a = 0.09+0.077*std::exp(-0.40*Red);
    const doublereal b = 0.40+0.770*std::exp(-0.04*Red);
    doublereal JgTotal = 0.0;
    for (size_t i=0;i<fuelNum_;i++)
    {
        JgTotal += particleJg_[i][iParcel]; // total kg/s
    }
    doublereal ub = JgTotal/(Pi*rhoGas*d0*d0);
    doublereal Reb = rhoGas * ub * d0 / muGas;
    doublereal f1 = (1+0.0545*Red+0.1*std::pow(Red,0.5)*(1-0.03*Red)) / (1+a*std::pow(Reb,b));
    return mass*(f1/taud)*(uGas-uDrop);
*/
    doublereal Cd = 0.0;

    // Wen-Yu drag model--------------------------------------------------+
    // F = 0.5*Cd*rhoGas*(uGas-uDrop)*|uGas-uDrop|*Ap                     |
    // { Cd = 24.0/Red * (1.0+0.15*pow(Red,0.687)) # Re<1000              |
    // { Cd = 0.44 # Re>=1000                                             |
    // End----------------------------------------------------------------+
    if (Red<1000)
    {
        Cd = 24.0/Red;
        Cd *= (1.0+0.15*std::pow(Red,0.687));
    }
    else
    {
        Cd = 0.44;
    }
    doublereal F = 0.5*Cd*rhoGas*(uGas-uDrop)*std::abs(uGas-uDrop)*Pi*d0*d0/4.0;
    return F;

}


void Lagrangian::solve()
{
    const size_t zSize = zField_.size();
    this->resize();
    doublereal mass = massDrop_;

    // outer loop: time marching    
    bool ok = false;
    while(!ok)
    {
        // inject a parcel every time step @ left most side with uInjection_
        particleVelocity_.push_back(uInjection_);
        particleDiameter_.push_back(diameterInjection_);
        particleTemperature_.push_back(TInjection_);
        particleQg_.push_back(0.0);
        for (int i=0;i<fuelNum_;i++)
            { particleJg_[i].push_back(0.0); }
        if (fuelNum_==1)
            component_[0].push_back(1.0);
        // Jet-A surrogate hard-coded
        if (fuelNum_==4)
            { component_[0].push_back(0.322); component_[1].push_back(0.446);
                component_[2].push_back(0.185); component_[3].push_back(0.047); }
        positions_.push_back(zField_[0]);
        parcelNumber_++;

        // inner loop: update each particle
        for (int iParcel=0;iParcel < parcelNumber_;iParcel++)
        {
            // store old values
            doublereal oldPosition = positions_[iParcel];
            doublereal oldTemperature = particleTemperature_[iParcel];
            doublereal oldDiameter = particleDiameter_[iParcel];
            doublereal oldVelocity = particleVelocity_[iParcel];

            // if particles leak out of the domain
            if (oldPosition > zField_[zSize-1] || oldPosition < zField_[0])
            {
                continue;
            }
            else
            {
                if (oldDiameter < delta_)
                {
                    particleQg_[iParcel] = 0.0;
                    for (int i=0;i<fuelNum_;i++)
                        { particleJg_[i][iParcel] = 0.0; }
                }
                if (oldDiameter > 0.0)
                {
                    this->heatAndMassTransfer(iParcel);
                    mass = rhoDrop_ * (4.0/3.0) * Pi * oldDiameter * oldDiameter * oldDiameter / 8.0;
                    // Second-order Runge-Kutta method
                    particleVelocity_[iParcel] += this->force(iParcel)*dt_/mass; // newVelocity
                    doublereal newVelocity = particleVelocity_[iParcel]; 
                    positions_[iParcel] = 0.5*(oldVelocity+newVelocity)*dt_ + oldPosition;
                }
            }
        }

        if (particleDiameter_[0] > delta_)
        {
            trackTemperature_.push_back(particleTemperature_[0]);
            trackDiameter_.push_back(particleDiameter_[0]);
            trackPosition_.push_back(positions_[0]);
            for (int i=0;i<fuelNum_;i++)
                { trackComponent_[i].push_back(component_[i][0]); }
        }
        if (parcelNumber_>20000) ok = true;
    }

    this->evaluateTransferRateField();
}


void Lagrangian::write()
{
    LagrangianOutfile_ << "#Transfer Rate Field#" << std::endl;
    if(fuelNum_==1)
    {
        for (int iz=0;iz<zField_.size();iz++)
        {
            LagrangianOutfile_ << zField_[iz] << "," 
                                << heatTransferRateField_[iz] << "," 
                                << massTransferRateField_[0][iz] << "," 
                                << totalMassTransferField_[iz]
                                << std::endl;
        }
    }
    else if (fuelNum_==4)
    {
        for (int iz=0;iz<zField_.size();iz++)
        {
            LagrangianOutfile_ << zField_[iz] << "," 
                                << heatTransferRateField_[iz] << "," 
                                << massTransferRateField_[0][iz] << "," 
                                << massTransferRateField_[1][iz] << "," 
                                << massTransferRateField_[2][iz] << "," 
                                << massTransferRateField_[3][iz] << "," 
                                << totalMassTransferField_[iz]
                                << std::endl;
        }
    }
    LagrangianOutfile_ << "#################################" << std::endl;
    LagrangianOutfile_ << "########Particle Tracking########" << std::endl;
    if(fuelNum_==1)
    {    
        for (int it=0;it<trackPosition_.size();it++)
        {
            LagrangianOutfile_ << it*dt_ << ","
                                << trackPosition_[it] << ","
                                << trackTemperature_[it]<< ","
                                << trackDiameter_[it] << ","
                                << trackComponent_[0][it] << std::endl;
        }
    }
    else if (fuelNum_==4)
    {
        for (int it=0;it<trackPosition_.size();it++)
        {
            LagrangianOutfile_ << it*dt_ << ","
                                << trackPosition_[it] << ","
                                << trackTemperature_[it]<< ","
                                << trackDiameter_[it] << ","
                                << trackComponent_[0][it] << "," 
                                << trackComponent_[1][it] << ","
                                << trackComponent_[2][it] << "," 
                                << trackComponent_[3][it] << std::endl;
        }
    }

    LagrangianOutfile_ << "#################################" << std::endl;
}


}
