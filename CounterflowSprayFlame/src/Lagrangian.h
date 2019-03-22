//! @file Lagrangian.h

#ifndef CT_LAGRANGIAN_H
#define CT_LAGRANGIAN_H

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


    // setters
    void setFuel(const std::vector<std::string> fuelName)
    {
        fuelName_ = fuelName;
        fuelNum_ = fuelName.size();
        massTransferRateField_.resize(fuelNum_);
        oldMassTransferField_.resize(fuelNum_);
        component_.resize(fuelNum_);
        particleJg_.resize(fuelNum_);
        trackComponent_.resize(fuelNum_);
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

    void setPressure(const doublereal p)
    {
        p0_ = p;
    }


    // getters
    doublereal getHeatTransferRate(doublereal z) const
    {
        if (loopCnt_==0) return 0.0;
        else if (loopCnt_==1) return 0.1*relaxationFactor_ * linearInterpolate(heatTransferRateField_, z);
        else return linearInterpolate(heatTransferRateField_, z);
    }

    doublereal getTotalMassTransferRate(doublereal z) const
    {
        if (loopCnt_==0) return 0.0;
        else if (loopCnt_==1) return relaxationFactor_ * linearInterpolate(totalMassTransferField_, z);
        else return linearInterpolate(totalMassTransferField_, z);
    }

    doublereal getMassTransferRate(size_t i, doublereal z) const
    {
        if (loopCnt_==0) return 0.0;
        else if (loopCnt_==1) return relaxationFactor_ * linearInterpolate(massTransferRateField_[i], z);
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
    std::vector<std::string> fuelName_;
    size_t fuelNum_;


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

    vector_fp liquidGasMassRatio_; // liquid phase
    vector_fp heatTransferRateField_; // J/m^s/s
    vector_fp oldHeatTransferField_;
    vector_fp totalMassTransferField_;
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