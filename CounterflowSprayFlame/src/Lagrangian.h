//! @file Lagrangian.h

#ifndef CT_LAGRANGIAN_H
#define CT_LAGRANGIAN_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include "cantera/base/global.h"
#include "StFlowSpr.h"

namespace Cantera
{
class StFlow;

class Lagrangian
{
public:
    Lagrangian( const doublereal diameterInjection = 50e-6,
                const doublereal mdotInjection = 0.0,
                const doublereal TInjection = 300.0 );
    Lagrangian() = default;

    void flow(StFlow& flow) {
        flow_ = &flow;
    }

    void setupFlowField(const vector_fp& solution);

    int loopcnt() const {
        return loopcnt_;
    }

    void setFuel(const std::vector<std::string>& fuelName) {
        fuelName_ = fuelName;
    }

    bool solve();

    void write() const;

    void outfile(const std::string& outfile) {
        outfile_ = outfile;
    }



private:

    // member functions
    void resize() {
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

    }


    void inject() {
        position_.push_back(z_[0]);
        velocity_.push_back(u_[0]);
        diameter_.push_back(diameterInjection_);
        temperature_.push_back(TInjection_);

        hTrans_.push_back(0.0);
        for (size_t i=0; i<mTrans_.size(); i++) {
            mTrans_[i].push_back(0.0);
        }

    }

    void calchmTrans();

    doublereal linearInterpolate(const vector_fp& field, const doublereal z) const;


    // initial data
    doublereal diameterInjection_;
    doublereal mdotInjection_;
    doublereal TInjection_;
    doublereal dt_;
    doublereal delta_;
    std::string outfile_;
    // fuel
    std::vector<std::string> fuelName_;


    // access to gas-phase
    StFlow* flow_;

    // gas-phase fields
    vector_fp z_;
    vector_fp rho_;
    vector_fp mu_;
    vector_fp cp_;
    vector_fp u_;
    vector_fp T_;
    std::vector<std::vector<doublereal> > Y_;


    // parcels
    vector_fp position_;
    vector_fp velocity_;
    vector_fp diameter_;
    vector_fp temperature_;

    vector_fp hTrans_;
    std::vector<std::vector<doublereal> > mTrans_;


    // loop
    int loopcnt_;

};

}

#endif
