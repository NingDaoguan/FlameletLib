//! @file Lagrangian.h

#ifndef CT_LAGRANGIAN_H
#define CT_LAGRANGIAN_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include "cantera/base/global.h"

namespace Cantera
{

class Lagrangian
{
public:
    // constructor
    Lagrangian( std::ofstream& LagrangianOutfile,
                size_t ignition=0,
                doublereal diameterInjection=100e-6,
                doublereal mdotInjection=0.01);

    Lagrangian() = default;
    // setters
    void setFuel(const std::vector<std::string>& fuelName)
    {
        fuelName_ = fuelName;
        fuelNum_ = fuelName.size();
        W_V.resize(fuelNum_);
        T_B.resize(fuelNum_);
        massTransferRateField_.resize(fuelNum_);
        oldMassTransferField_.resize(fuelNum_);
        component_.resize(fuelNum_);
        particleJg_.resize(fuelNum_);
        trackComponent_.resize(fuelNum_);

        // Liquid thermophysical properties
        if (fuelName_[0] == "C2H5OH")
        {
            W_V[0] = 46.0e-3;
            T_B[0] = 352.0;
            lat_a3 = -4.64447747e-03;
            lat_a2 =  3.69964540e+00;
            lat_a1 = -1.28942888e+03;
            lat_a0 =  1.05313435e+06;
            cp_a1 = 4.037;
            cp_a0 = 1124.9;
        }
        else if(fuelNum_ == 4)
        {
            W_V[0] = 170.0e-3;
            T_B[0] = 490.0;
            W_V[1] = 226.0e-3;
            W_V[2] = 138.0e-3;
            W_V[3] = 92.00e-3;
            T_B[1] = 513.0;
            T_B[2] = 460.0;
            T_B[3] = 384.0;
            lat_a3 = -6.72216465e-03;
            lat_a2 =  8.01138467e+00;
            lat_a1 = -3.62211999e+03;
            lat_a0 =  9.05036560e+05;
            cp_a1 = 3.84;
            cp_a0 = 1056.88;
        }
        else
        {
            W_V[0] = 170.0e-3;
            T_B[0] = 490.0;
            lat_a3 = -6.72216465e-03;
            lat_a2 =  8.01138467e+00;
            lat_a1 = -3.62211999e+03;
            lat_a0 =  9.05036560e+05;
            cp_a1 = 3.84;
            cp_a0 = 1056.88;
        }
    }

    void setInjectionVelocity(const doublereal uInjection)
    {
        uInjection_ = uInjection;
    }

    void setInjectionTemperature(const doublereal TInjection)
    {
        TInjection_ = TInjection;
    }

    void setRelaxationFractor(const doublereal relaxationFactor)
    {
        relaxationFactor_ = relaxationFactor;
    }

    void setZField(const vector_fp& zField)
    {
        zField_ = zField;
    }

    void setRhoField(const vector_fp& rhoField)
    {
        rhoField_ = rhoField;
    }

    void setVelocityField(const vector_fp& uField)
    {
        uField_ = uField;
        setInjectionVelocity(uField_[0]);
    }

    void setTemperatureField(const vector_fp& TField)
    {
        TField_ = TField;
    }

    void setMassFractionField(const std::vector<std::vector<double> >& YField)
    {
        YField_.resize(fuelNum_);

        for (int i=0;i<fuelNum_;i++)
            YField_[i] = YField[i];
    }

    void setMuField(const vector_fp& muField)
    {
        muField_ = muField;
    }

    void setConField(const vector_fp& conField)
    {
        conField_ = conField;
    }

    void setCpField(const vector_fp& cpField)
    {
        cpField_ = cpField;
    }

    void setPressure(const doublereal p)
    {
        p0_ = p;
    }


    // getters
    doublereal getHeatTransferRate(doublereal z) const
    {
        if (loopCnt_==0) return 0.0;
        else return linearInterpolate(heatTransferRateField_, z);
    }

    doublereal getTotalMassTransferRate(doublereal z) const
    {
        if (loopCnt_==0) return 0.0;
        else return linearInterpolate(totalMassTransferField_, z);
    }

    doublereal getMassTransferRate(size_t i, doublereal z) const
    {
        if (loopCnt_==0) return 0.0;
        else return linearInterpolate(massTransferRateField_[i], z);
    }

    doublereal getMassRatio(size_t iz) const
    {
        return liquidGasMassRatio_[iz];
    }

    bool checkOK() const
    {
        return ok_;
    }


    void solve();

    void write();

private:

    // private member functions
    void resize()
    {
        parcelNumber_ = 0;
        particleVelocity_.resize(0);
        positions_.resize(0);
        particleTemperature_.resize(0);
        particleDiameter_.resize(0);
        particleQg_.resize(0);

        trackPosition_.resize(0);
        trackDiameter_.resize(0);
        trackTemperature_.resize(0);

        liquidGasMassRatio_.resize(0);
        heatTransferRateField_.resize(0);
        totalMassTransferField_.resize(0);
        for(size_t i=0;i<fuelNum_;i++)
        {
            massTransferRateField_[i].resize(0);
            component_[i].resize(0);
            particleJg_[i].resize(0);
            trackComponent_[i].resize(0);
        }
    }

    void evaluateTransferRateField();

    doublereal evaluateResidual() const
    {
        const size_t zSize = zField_.size();
        doublereal sum = 0.0;
        doublereal sum2 = 0.0;
        doublereal dz = 0.0;
        for (int iz=0;iz<zSize-1;iz++)
        {
            dz = zField_[iz+1] - zField_[iz];
            for (int i=0;i<fuelNum_;i++)
            {
                sum += ( massTransferRateField_[i][iz] - linearInterpolate(oldMassTransferField_[i], oldGrid_, zField_[iz]) )
                     * ( massTransferRateField_[i][iz] - linearInterpolate(oldMassTransferField_[i], oldGrid_, zField_[iz]) )
                     * dz * dz;
                sum2 += massTransferRateField_[i][iz] * dz;
            }
        }
        sum /= zSize;
        sum2 /= zSize;
        return std::pow(sum,0.5) / (sum2+delta_);
    }

    void heatAndMassTransfer(const int iParcel);

    doublereal getLatentHeat(const doublereal T)
    {
        return lat_a0 + lat_a1*T + lat_a2*T*T + lat_a3*T*T*T;
    }

    doublereal getCpL(const doublereal T)
    {
        return cp_a0 + cp_a1*T;
    }

    doublereal linearInterpolate(const vector_fp& field, const doublereal z) const;

    doublereal linearInterpolate(const vector_fp& field, const vector_fp grid, const doublereal z) const;

    doublereal force(const int iParcel) const;


    // private member data
    doublereal p0_; // pressure
    doublereal TInjection_; // injection temperature
    doublereal mdotInjection_; // kg/m^2/s
    doublereal diameterInjection_; // m
    doublereal dt_;
    doublereal rhoDrop_; // kg/m^3
    doublereal uInjection_; // m/s
    doublereal volumeDrop_; // m^3
    doublereal massDrop_; // kg
    doublereal nPerSec_;
    doublereal nPerParcel_;
    doublereal relaxationFactor_;
    doublereal delta_;
    doublereal residual0_;
    size_t loopCnt_;
    size_t ignition_;
    bool ok_;


    // fuel
    std::vector<std::string> fuelName_;
    size_t fuelNum_;
    vector_fp W_V;
    vector_fp T_B;
    doublereal lat_a0;
    doublereal lat_a1;
    doublereal lat_a2;
    doublereal lat_a3;
    doublereal cp_a0;
    doublereal cp_a1;

    // parcel
    size_t parcelNumber_;
    vector_fp positions_; // position vector of each parcel
    vector_fp particleVelocity_;
    vector_fp particleTemperature_;
    vector_fp particleDiameter_;
    vector_fp particleQg_;
    std::vector<std::vector<double> > particleJg_;
    std::vector<std::vector<double> > component_;


    // grid
    vector_fp zField_;
    vector_fp rhoField_;
    vector_fp uField_;
    vector_fp TField_;
    std::vector<std::vector<double> > YField_;
    vector_fp muField_;
    vector_fp conField_;
    vector_fp cpField_;

    vector_fp liquidGasMassRatio_; // liquid phase
    vector_fp heatTransferRateField_; // J/m^s/s
    vector_fp oldHeatTransferField_;
    vector_fp totalMassTransferField_;
    vector_fp oldTotalMassTransferField_;
    std::vector<std::vector<double> > massTransferRateField_; // kg/m^3/s
    std::vector<std::vector<double> > oldMassTransferField_;
    vector_fp oldGrid_;


    // time
    vector_fp trackPosition_;
    vector_fp trackDiameter_;
    vector_fp trackTemperature_;
    std::vector<std::vector<double> > trackComponent_;

    std::ofstream& LagrangianOutfile_;
};

}

#endif
