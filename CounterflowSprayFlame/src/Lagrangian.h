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

        W_V = 46.0e-3;
        T_B = 352.0;
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
        uTrans_.clear();
        hTrans_.clear();
        mTrans_.resize(fuelName_.size());
        for (size_t i=0; i<mTrans_.size(); i++) {
            mTrans_[i].clear();
        }

        // clear tracking data
        tp_.clear();
        tt_.clear();
        td_.clear();
    }


    void inject() {
        position_.push_back(z_[0]);
        velocity_.push_back(u_[0]);
        diameter_.push_back(diameterInjection_);
        temperature_.push_back(TInjection_);

        uTrans_.push_back(0.0);
        hTrans_.push_back(0.0);
        for (size_t i=0; i<mTrans_.size(); i++) {
            mTrans_[i].push_back(0.0);
        }

    }

    void calcTrans(int ip);

    doublereal linearInterpolate(const vector_fp& field, const doublereal z) const;

    // ethanol
    doublereal rhod(doublereal T) const {
        return 70.1308387/std::pow(0.26395, 1 + std::pow(1 - T/516.25, 0.2367));
    }
    doublereal latentHeat(doublereal T) const {
        doublereal Tr = T/516.25;
        return 958345.091059064*std::pow(1 - Tr, ((0.0*Tr + 0.0)*Tr + 0.75362)*Tr + -0.4134);
    }
    doublereal cpd(doublereal T) const {
        return ( (((0.0*T + 0.0)*T + 5.20523562482363e-05)*T+ 0.00714146172046278)*T
                 + -1.21990926653498 )*T + 2052.57331394213;
    }


    // initial data
    doublereal diameterInjection_;
    doublereal mdotInjection_;
    doublereal TInjection_;
    doublereal dt_;
    doublereal small;
    doublereal p0_;
    std::string outfile_;
    // fuel
    std::vector<std::string> fuelName_;
    doublereal W_V;
    doublereal T_B;


    // access to gas-phase
    StFlow* flow_;

    // gas-phase
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

    vector_fp uTrans_;
    vector_fp hTrans_;
    std::vector<std::vector<doublereal> > mTrans_;


    // tracking
    vector_fp tp_;
    vector_fp td_;
    vector_fp tt_;

    // loop
    int loopcnt_;

};

}

#endif
