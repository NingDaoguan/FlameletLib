#include <fstream>

#include "cantera/thermo.h"
#include "Inlet1DFree.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "Sim1DFree.h"
#include "Sim1DFree.cpp"
#include "StFlowFree.h"
#include "StFlowFree.cpp"
#include "OneDimFree.h"
#include "OneDimFree.cpp"


using namespace Cantera;
using fmt::print;

void propagatingPremixedFlame()
{
    doublereal phi=1.0;
    doublereal C_atoms;
    doublereal H_atoms;
    std::string kinFileName;
    std::string kinPhase;
    std::vector<std::string> fuelName;
    std::vector<double> fuelX;
    std::ifstream infile("input.txt");
    doublereal inputValue;
    std::string name;
    size_t fuel=0;
    doublereal temp = 300.0; // K
    doublereal pressure = 1.0*OneAtm; //atm
    doublereal uin = 0.4; //m/sec
    while (infile >> name)
    {
        infile >> inputValue;
        if (name == "fuel")
            { fuel = size_t(inputValue); continue; }
        if (name == "temperature")
            { temp = inputValue; continue; }
        if (name == "phi")
            { phi = inputValue; continue; }
    }
    #include "FuelInfo.h"
    IdealGasMix gas(kinFileName, kinPhase);

    size_t nsp = gas.nSpecies();
    vector_fp x(nsp, 0.0);

    doublereal ax = C_atoms + H_atoms / 4.0;
    doublereal fa_stoic = 1.0 / (4.76 * ax);

    for (size_t i=0;i<fuelName.size();i++)
    {
        x[gas.speciesIndex(fuelName[i])] = fuelX[i];
    }
    if (fuel!=0)
    {
        x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic;
        x[gas.speciesIndex("N2")] = 0.78 / phi / fa_stoic;
        x[gas.speciesIndex("AR")] = 0.01 / phi / fa_stoic;
    }
    else
    {
        x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic;
        x[gas.speciesIndex("N2")] = 0.79 / phi / fa_stoic;
    }
    gas.setState_TPX(temp, pressure, x.data());
    doublereal rho_in = gas.density();

    vector_fp yin(nsp);
    gas.getMassFractions(&yin[0]);

    gas.equilibrate("HP");
    vector_fp yout(nsp);
    gas.getMassFractions(&yout[0]);
    doublereal rho_out = gas.density();
    doublereal Tad = gas.temperature();
    print("phi = {}, Tad = {}\n", phi, Tad);

    //=============  build each domain ========================


    //-------- step 1: create the flow -------------

    StFlow flow(&gas);
    flow.setFreeFlow();

    // create an initial grid
    int nz = 6;
    doublereal lz = 0.1;
    vector_fp z(nz);
    doublereal dz = lz/((doublereal)(nz-1));
    for (int iz = 0; iz < nz; iz++) {
        z[iz] = ((doublereal)iz)*dz;
    }

    flow.setupGrid(nz, &z[0]);

    // specify the objects to use to compute kinetic rates and
    // transport properties

    Transport* tr = newTransportMgr("Mix",&gas);
    flow.setTransport(*tr);

    flow.setKinetics(gas);
    flow.setPressure(pressure);

    //------- step 2: create the inlet  -----------------------

    Inlet1D inlet;

    inlet.setMoleFractions(x.data());
    doublereal mdot=uin*rho_in;
    inlet.setMdot(mdot);
    inlet.setTemperature(temp);


    //------- step 3: create the outlet  ---------------------

    Outlet1D outlet;

    //=================== create the container and insert the domains =====

    std::vector<Domain1D*> domains { &inlet, &flow, &outlet };
    Sim1D flame(domains);

    //----------- Supply initial guess----------------------

    vector_fp locs{0.0, 0.3, 0.7, 1.0};
    vector_fp value;

    double uout = inlet.mdot()/rho_out;
    value = {uin, uin, uout, uout};
    flame.setInitialGuess("u",locs,value);
    value = {temp, temp, Tad, Tad};
    flame.setInitialGuess("T",locs,value);

    for (size_t i=0; i<nsp; i++) {
        value = {yin[i], yin[i], yout[i], yout[i]};
        flame.setInitialGuess(gas.speciesName(i),locs,value);
    }

    inlet.setMoleFractions(x.data());
    inlet.setMdot(mdot);
    inlet.setTemperature(temp);

    // flame.showSolution();

    int flowdomain = 1;
    double ratio = 10.0;
    double slope = 0.2;
    double curve = 0.2;

    flame.setRefineCriteria(flowdomain,ratio,slope,curve);

    int loglevel=1;

    // Solve freely propagating flame

    // Linearly interpolate to find location where this temperature would
    // exist. The temperature at this location will then be fixed for
    // remainder of calculation.
    flame.setFixedTemperature(0.5 * (temp + Tad));
    flow.solveEnergyEqn();
    bool refine_grid = true;

    flame.solve(loglevel,refine_grid);
    double flameSpeed_mix = flame.value(flowdomain,flow.componentIndex("u"),0);
    print("Flame speed with mixture-averaged transport: {} m/s\n",
          flameSpeed_mix);

    vector_fp zvec,Tvec,COvec,CO2vec,Uvec;
    std::vector<std::vector<double> > SpeciesVec(nsp);

    print("\n{:9s}\t{:8s}\t{:5s}\t{:7s}\n",
          "z (m)", "T (K)", "U (m/s)", "Y(CO)");
    for (size_t n = 0; n < flow.nPoints(); n++) {
        Tvec.push_back(flame.value(flowdomain,flow.componentIndex("T"),n));
        COvec.push_back(flame.value(flowdomain,flow.componentIndex("CO"),n));
        CO2vec.push_back(flame.value(flowdomain,flow.componentIndex("CO2"),n));
        Uvec.push_back(flame.value(flowdomain,flow.componentIndex("u"),n));
        zvec.push_back(flow.grid(n));
        for (size_t k=0; k<nsp; k++)
        {
            SpeciesVec[k].push_back(flame.value(flowdomain,k+5,n));
        }
    }

    print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
    print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

    std::ofstream outfile("output.csv");
    size_t nCol = nsp+6;
    for (size_t n=0; n<nCol; n++)
    {
        if (n==0) outfile << "z";
        else outfile << flow.componentName(n-1);
        if (n!=nCol-1) outfile << ",";
    }
    outfile << std::endl;
    for (size_t n=0; n<flow.nPoints(); n++)
    {
        for (size_t k=0; k<nCol; k++)
        {
            if (k==0) outfile << flow.grid(n);
            else outfile << flame.value(flowdomain,k-1,n);
            if (k!=nCol-1) outfile << ",";
        }
        outfile << std::endl;
    }
}

int main()
{
    try {
        propagatingPremixedFlame();
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
}
