#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>

#include "cantera/thermo.h"
#include "Inlet1DSpr.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "Lagrangian.cpp"
#include "Sim1DSpr.cpp"
#include "StFlowSpr.cpp"
#include "OneDimSpr.cpp"

using namespace Cantera;

void counterflowDiffusionFlame(doublereal mdotF, doublereal mdotO, doublereal len)
{
    size_t tstepsSize = 6; // take tsteps[i] time steps
    const int tsteps[tstepsSize] = {10,20,10,20,20,20};
    doublereal p0 = OneAtm;
    std::string kinFileName;
    std::string kinPhase;
    std::vector<std::string> fuelName;
    std::vector<double> fuelX;
    doublereal C_atoms;
    doublereal H_atoms;
    std::string compLeft;
    std::string compRight;
    std::string compInit;
    // input from file
    std::ifstream infile("input.txt");
    if (!infile) std::cerr << "input.txt NOT FOUND!" << std::endl;
    std::string name; // name column
    doublereal value; // value column
    size_t fuel = 0; // 0 for Kerosene test
    doublereal mdotInjection; // kg/m^2/s
    doublereal phi = 1.0;
    doublereal diameterInjection = 100e-6; // m
    doublereal uInjection = 0.2; // m/s
    doublereal TInjection = 300; // K
    doublereal LagrangianRelaxationFactor = 1.0;
    doublereal TinLeft = 300; // fuel side
    doublereal TinRight = 300; // oxidizer side
    doublereal Teq;
    doublereal mdotLeft = mdotF; // kg/m^2/s
    doublereal mdotRight = mdotO; // kg/m^2/s
    size_t points = 51;
    doublereal width = len;
    doublereal a0 = 50;
    size_t ignition = 0; // 1 with ignition
    size_t rightType = 0; // 0 for air, 1 for product and 2 for fuel @ right side
    while (infile >> name)
    {
        infile >> value;
        if (name == "fuel")
            { fuel = size_t(value); continue; }
        if (name == "points")
            { points = size_t(value); continue; }
        if (name == "phi")
            { phi = value; continue; }
        if (name == "diameterInjection")
            { diameterInjection = value; continue; }
        if (name == "uInjection")
            { uInjection = value; continue; }
        if (name == "TInjection")
            { TInjection = value; continue; }
        if (name == "LagrangianRelaxationFactor")
            { LagrangianRelaxationFactor = value; continue; }
        if (name == "TinLeft")
            { TinLeft = value; continue; }
        if (name == "TinRight")
            { TinRight = value; continue; }
        if (name == "a0")
            { a0 = value; continue; }
        if (name=="ignition")
            { ignition = size_t(value); continue; }
        if (name=="rightType")
            { rightType = size_t(value); continue; }
    }
    doublereal minGrid = width/500;
    #include "FuelInfo.h"
    std::cout << "#Input Parameters" << std::endl;
    std::cout << "#\tDomain      :\t"
                << points << " "
                << width << std::endl;
    std::cout << "#\tLagrangian  :\t"
                << phi << " "
                << diameterInjection << " "
                << TInjection << std::endl;
    std::cout << "#\tTemperature :\t"
                << TinLeft << std::endl;
    std::cout << "#\tmdot        :\t"
                << mdotLeft << " "
                << mdotRight << std::endl;
    std::cout << std::endl;
    std::cout << "=====================================================================" << std::endl;



    // create gas, Lagrangian particle cloud and flow field
    IdealGasMix gas(kinFileName, kinPhase);
    size_t nsp = gas.nSpecies();
    std::cout << "\nNumber of species:\t" << nsp << std::endl;

    // calculate x according to phi
    vector_fp x(nsp, 0.0);
    vector_fp y(nsp, 0.0);
    doublereal mixtureFraction = 0.0;
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
    gas.setState_TPX(480, p0, x.data());
    gas.getMassFractions(&y[0]);
    for (size_t i=0;i<fuelName.size();i++)
    {
        mixtureFraction += y[gas.speciesIndex(fuelName[i])];
    }
    doublereal yI[nsp];
    doublereal* yIn = &yI[0];
    gas.getMassFractions(yIn);

    std::cout << "Mixture fraction:\t" << mixtureFraction << std::endl;
    mdotInjection = mdotLeft*mixtureFraction / (1 - mixtureFraction);
    std::cout << "mdotInjection:\t" << mdotInjection << std::endl;
    gas.equilibrate("HP");
    doublereal Tad = gas.temperature();


    std::stringstream ss1;
    ss1 << "LAGmf-" << mdotF << "mo-" << mdotO << "L-" << len << "_raw.txt";
    std::ofstream LagrangianOutfile(ss1.str());
    std::ofstream& LagrangianOutfileRef = LagrangianOutfile;
    Lagrangian particleCloud(LagrangianOutfileRef, ignition, diameterInjection, mdotInjection);
    particleCloud.setFuel(fuelName);
    particleCloud.setInjectionVelocity(uInjection);
    particleCloud.setInjectionTemperature(TInjection);
    particleCloud.setRelaxationFractor(LagrangianRelaxationFactor);
    particleCloud.setPressure(p0);
    Lagrangian& particleCloudRef = particleCloud;
    StFlow flow(particleCloudRef, &gas, nsp, points);
    flow.setAxisymmetricFlow();
    flow.setFuelName(fuelName);


    // flow configuration
    doublereal zpoints[points];
    doublereal* z = &zpoints[0];
    doublereal dz = width / ((doublereal)(points-1));
    std::cout << "\nGrid points:" << std::endl;
    for (int iz=0;iz<points;iz++)
    {
        z[iz] = ((doublereal)iz)*dz;
        std::cout << z[iz] << " ";
        if (iz == points-1)
            std::cout << std::endl;
    }
    flow.setupGrid(points,z);
    Transport* tr = newTransportMgr("Mix",&gas);
    flow.setTransport(*tr);
    flow.setKinetics(gas);
    Inlet1D left;
    left.setMoleFractions(compLeft);
    left.setMdot(mdotLeft);
    left.setTemperature(TinLeft);
    Inlet1D right;

    for (size_t i=0;i<fuelName.size();i++)
    {
        x[gas.speciesIndex(fuelName[i])] = fuelX[i];
    }
    if (fuel!=0)
    {
        x[gas.speciesIndex("O2")] = 0.21 / (phi) / fa_stoic;
        x[gas.speciesIndex("N2")] = 0.78 / (phi) / fa_stoic;
        x[gas.speciesIndex("AR")] = 0.01 / (phi) / fa_stoic;
    }
    else
    {
        x[gas.speciesIndex("O2")] = 0.21 / (phi) / fa_stoic;
        x[gas.speciesIndex("N2")] = 0.79 / (phi) / fa_stoic;
    }
    gas.setState_TPX(480, p0, x.data());
    vector_fp xR(nsp, 0.0);
    gas.equilibrate("HP");
    gas.getMoleFractions(&xR[0]);
    Teq = gas.temperature();
    if (rightType==0)
    {    
        right.setMoleFractions(compRight);
        right.setMdot(mdotRight);
        right.setTemperature(TinRight);
    }
    else if (rightType==1)
    {
        TinRight = Teq;
        compRight = "";
        for (size_t k=0; k<nsp; k++)
        {
            compRight.append(gas.speciesName(k));
            compRight.append(":");
            compRight.append(to_string(xR[k]));
            if (k!=nsp-1) compRight.append(", ");
        }
        right.setMoleFractions(compRight);
        right.setMdot(mdotRight);
        right.setTemperature(TinRight);
    }
    else
    {
        TinRight = 480.0;
        compRight = "";
        vector_fp xF(nsp,0.0);
        for (size_t i=0;i<fuelName.size();i++)
        {
            xF[gas.speciesIndex(fuelName[i])] = fuelX[i];
        }
        for (size_t k=0; k<nsp; k++)
        {
            compRight.append(gas.speciesName(k));
            compRight.append(":");
            compRight.append(to_string(xF[k]));
            if (k!=nsp-1) compRight.append(", ");
        }
        right.setMoleFractions(compRight);
        right.setMdot(mdotRight);
        right.setTemperature(TinRight);
    }
    std::vector<Domain1D*> domains;
    domains.push_back(&left);
    domains.push_back(&flow);
    domains.push_back(&right);
    Sim1D flame(domains);
    std::cout << "\nBoundary conditions @ left side: " << std::endl;
    left.showSolution(0);
    std::cout << "Boundary conditions @ right side: " << std::endl;
    right.showSolution(0);


    // init
    if (ignition == 0) Tad = TinLeft;
    vector_fp locs(7);
    vector_fp values(7);
    locs[0]=0;
    locs[1]=0.3;
    locs[2]=0.4;
    locs[3]=0.5;
    locs[4]=0.6;
    locs[5]=0.7;
    locs[6]=1;
    values[0]=TinLeft;
    values[1]=0.9*TinLeft+0.1*Tad;
    values[2]=0.7*Tad+0.3*TinLeft;
    values[3]=0.8*Tad+0.1*TinLeft+0.1*TinRight;
    values[4]=0.7*Tad+0.3*TinRight;
    values[5]=0.9*TinRight+0.1*Tad;
    values[6]=TinRight;
    flame.setInitialGuess("T",locs,values);
    values[0]=0.5*a0*(1*width);
    values[1]=0.5*a0*(0.2*width);
    values[2]=0.5*a0*(0.1*width);
    values[3]=0.5*a0*(0.0*width);
    values[4]=0.5*a0*(-0.1*width);
    values[5]=0.5*a0*(-0.2*width);
    values[6]=0.5*a0*(-1*width);
    flame.setInitialGuess("u",locs,values);
    for (int ik=0; ik<nsp;ik++)
    {
       values[0]=left.massFraction(ik);
       values[1]=0.8*left.massFraction(ik) + 0.2*yIn[ik];
       values[2]=0.5*( left.massFraction(ik) + yIn[ik] );
       values[3]=yIn[ik];
       values[4]=0.5*( right.massFraction(ik) + yIn[ik] );
       values[5]=0.8*right.massFraction(ik) + 0.2*yIn[ik];
       values[6]=right.massFraction(ik);
       flame.setInitialGuess(gas.speciesName(ik), locs, values);
    }


    // solve
    size_t domFlow = 1;
    flame.setMaxGridPoints(domFlow,2000);
    flame.setGridMin(domFlow, minGrid);
    flow.setPressure(p0);
    flow.solveEnergyEqn();
    flame.setRefineCriteria(domFlow,10.0,0.8,0.8,-0.1);
    // flame.setRefineCriteria(domFlow,10.0,1,1,-0.1);
    flame.setTimeStep(1.0e-5,tstepsSize,&tsteps[0]);
    flame.solve(1,true);


    // output
    std::stringstream ss2;
    ss2 << "mf-" << mdotF << "mo-" << mdotO << "L-" << len << "_raw.txt";
    std::ofstream outfile(ss2.str());
    int np = flow.nPoints();
    Array2D solution(np,nsp+7,0.0);
    vector_fp grid;
    grid.resize(np);
    for (int iz=0;iz<np;iz++)
    {
        grid[iz] = flow.grid(iz);
    }
    for (int n=0; n<np; n++)
    {
        solution(n,0) = grid[n];
        for (int i=1;i<5;i++)
        {
            solution(n,i) = flame.value(domFlow,i-1,n);
        }
        for (int i=5; i<nsp+6;i++)
        {
            if (flame.value(domFlow, i-1,n) > 0.0) {solution(n,i) = flame.value(domFlow, i-1,n);}
            else if (flame.value(domFlow, i-1,n) < 0.0)
            {
                solution(n,i) = 0.0;
                flame.setValue(domFlow, i-1,n,0.0);
            }
        }
        solution(n,nsp+6) = particleCloud.getMassRatio(n);
    }

    for (int j=0;j<solution.nColumns();j++)
    {
        if (j==0) outfile << "z";
        else if (j<solution.nColumns()-1) outfile << flow.componentName(j-1);
        else outfile << "alpha";
        if (j!=solution.nColumns()-1) outfile << ",";
    }
    outfile << std::endl;
    
    for (int i=0;i<solution.nRows();i++)
    {
        for (int j=0;j<solution.nColumns();j++)
        {
            outfile << solution(i,j);
            if (j!=solution.nColumns()-1) outfile << ",";
        }
        outfile << std::endl;
    }
}

int main(int argc, char *argv[])
{
    doublereal mdotF = 0.2;
    doublereal mdotO = 0.2;
    doublereal len = 0.1;

    if (argc == 4)
    {
        mdotF = atof(argv[1]);
        mdotO = atof(argv[2]);
        len = atof(argv[3]);
    }
    else
        return -1;

    try {
        counterflowDiffusionFlame(mdotF, mdotO, len);
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
}
