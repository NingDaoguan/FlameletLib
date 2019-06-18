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
                const doublereal TInjection = 300.0,
                const doublereal p0 = 101325.0 );
    Lagrangian() = default;

    void flow(StFlow& flow) {
        flow_ = &flow;
    }

    bool setupFlowField(const vector_fp& solution);

    void setRelaxation(const doublereal& rlxf) {
        rlxf_ = rlxf;
    }

    int loopcnt() const {
        return loopcnt_;
    }

    void setFuel(const std::vector<std::string>& fuelName) {
        fuelName_ = fuelName;

        W_V = 46.0e-3;
        T_B = 352.0;
    }

    std::vector<std::string> fuelName() const {
        return fuelName_;
    }

    void solve();

    void write() const;

    void outfile(const std::string& outfile) {
        outfile_ = outfile;
    }


    // access
    doublereal htr(const doublereal z) const {
        return  linearInterpolate(htf_, z);
    }

    doublereal mtr(const size_t i, const doublereal z) const {
        return std::max(linearInterpolate(mtf_[i], z), small);
    }

    doublereal tmtr(const doublereal z) const {
        return  std::max(linearInterpolate(tmtf_, z), small);
    }



private:

    // member functions
    void clear();

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

    void track() {
        tp_.push_back(position_[0]);
        td_.push_back(diameter_[0]);
        tt_.push_back(temperature_[0]);
    }

    void calcTrans(int ip);

    void evalTrans();

    void evalRsd();

    void relax();

    doublereal dNdz(int j) const {
        int jloc = (u_[j] > 0.0 ? j : j + 1);
        return (N_[jloc] - N_[jloc-1]) / (z_[jloc] - z_[jloc-1]);
    }
    // doublereal drdadz(int j) const {
    //     int jloc = (u_[j] > 0.0 ? j : j + 1);
    //     return (rhodAlpha_[jloc] - rhodAlpha_[jloc-1]) / (z_[jloc] - z_[jloc-1]);
    // }
    doublereal dudz(int j) const {
        int jloc = (u_[j] > 0.0 ? j : j + 1);
        return (u_[jloc] - u_[jloc-1]) / (z_[jloc] - z_[jloc-1]);
    }
    // void relax(vector_fp& field1, const vector_fp& field0) {
    //     for (size_t i=0; i<field1.size(); i++) {
    //         field1[i] = rlxf_*field1[i]
    //                     + (1 - rlxf_)*linearInterpolate(field0, zOld_, z_[i]);
    //     }
    // }

    void scale(vector_fp& field, const doublereal fc) {
        for (auto it=field.begin(); it!=field.end(); it++) {
            *it *= fc;
        }
    }

    // doublereal csArea(size_t iz) const;

    // doublereal sumevap(size_t iz) const;

    doublereal linearInterpolate(const vector_fp& field,
                                 const doublereal z) const;

    doublereal linearInterpolate(const vector_fp& field,
                                 const vector_fp& grid,
                                 const doublereal z) const;


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
    doublereal Vd_;
    doublereal dt_;
    doublereal small;
    doublereal p0_;
    doublereal rlxf_;
    std::string outfile_;
    doublereal np_;
    doublereal nPerSec_;
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
    vector_fp V_;
    vector_fp T_;
    std::vector<std::vector<doublereal> > Y_;


    // parcels
    vector_fp position_;
    vector_fp velocity_;
    vector_fp diameter_;
    vector_fp temperature_;
    vector_fp hTrans_;
    std::vector<std::vector<doublereal> > mTrans_;


    // tracking
    vector_fp tp_;
    vector_fp td_;
    vector_fp tt_;


    vector_fp N_;
    // vector_fp rhodAlpha_;

    // transfer fields
    vector_fp htf_;
    std::vector<std::vector<doublereal> > mtf_;
    vector_fp tmtf_;

    // old fields
    vector_fp zOld_;
    vector_fp TOld_;
    // vector_fp htfOld_;
    // std::vector<std::vector<doublereal> > mtfOld_;


    // loop
    doublereal rsd_;
    int loopcnt_;

};

}

#endif
