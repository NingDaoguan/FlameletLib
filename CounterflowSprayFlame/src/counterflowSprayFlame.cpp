#include "cantera/thermo.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "Inlet1DSpr.h"
#include "Sim1DSpr.h"
#include "StFlowSpr.h"
#include "Lagrangian.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>

using namespace Cantera;

void counterflowSprayFlame(doublereal mdotL, doublereal mdotR, doublereal len)
{
    // data
    doublereal p0(101325.0);
    std::string kinFileName("Ethanol_31.cti");
    doublereal C_atoms(2);
    doublereal H_atoms(6);
    doublereal O_atoms(1);
    std::vector<std::string> fuelName;
    std::vector<double> fuelX;
    std::string compLeft;
    std::string compRight;
    
    doublereal mdotInjection; // kg/m^2/s
    doublereal phi(1.0); // equivalence ratio
    doublereal diameterInjection(50e-06); // m
    doublereal uInjection(0.2); // m/s
    doublereal TInjection(300.0); // K
    doublereal LagrangianRelaxationFactor(1.0);
    doublereal TinLeft(300.0); // fuel side
    doublereal TinRight(300.0); // oxidizer side
    doublereal a0(40);
    int points(41);
    bool ignition(false); // 1 with ignition
    const doublereal mdotLeft = mdotL; // kg/m^2/s
    const doublereal mdotRight = mdotR; // kg/m^2/s
    const doublereal width = len; // m


    // read from input.txt
    std::ifstream infile("input.txt");
    if (!infile) throw std::runtime_error("input.txt NOT FOUND!");
    std::string name; // name column
    std::string value; // value column
    while (infile >> name)
    {
        infile >> value;
        if (name == "kinFileName")
            { kinFileName = value; continue; }
        else if (name == "C_atoms")
            { C_atoms = std::stod(value); continue; }
        else if (name == "H_atoms")
            { H_atoms = std::stod(value); continue; }
        else if (name == "O_atoms")
            { O_atoms = std::stod(value); continue; }
        else if (name == "fuelName") {
            if (value != "(") std::cerr << "fuelName ( fuel1 fuel2 )";
            while (infile >> value) {
                if (value == ")") break;
                else fuelName.push_back(value);
            }
        }
        else if (name == "fuelX") {
            if (value != "(") std::cerr << "fuelX ( fuelX1 fuelX2 )";
            while (infile >> value) {
                if (value == ")") break;
                else fuelX.push_back(std::stod(value));
            }
        }
        else if (name == "p0")
            { p0 = std::stod(value); continue; }
        else if (name == "compLeft")
            { compLeft = value; continue; }
        else if (name == "compRight")
            { compRight = value; continue; }
        else if (name == "phi")
            { phi = std::stod(value); continue; }
        else if ( name == "diameterInjection")
            { diameterInjection = std::stod(value); continue; }
        else if (name == "TInjection")
            { TInjection = std::stod(value); continue; }
        else if (name == "LagrangianRelaxationFactor")
            { LagrangianRelaxationFactor = std::stod(value); continue; }
        else if (name == "TinLeft")
            { TinLeft = std::stod(value); continue; }
        else if (name == "TinRight")
            { TinRight = std::stod(value); continue; }
        else if (name == "a0")
            { a0 = std::stod(value); continue; }
        else if (name == "points")
            { points = std::stoi(value); continue; }
        else if (name=="ignition") {
            ignition = (value == "true") ? true : false;
        }
        else continue;
    }
    doublereal minGrid = diameterInjection;
    std::cout << "#Input Parameters" << std::endl;
    std::cout << "#\tDomain      :\t"
                << points << "\t"
                << len << std::endl;
    std::cout << "#\tLagrangian  :\t"
                << phi << "\t"
                << diameterInjection << "\t"
                << TInjection << std::endl;
    std::cout << "#\tTemperature :\t"
                << TinLeft << "\t"
                << TinRight << std::endl;
    std::cout << "#\tmdot        :\t"
                << mdotLeft << "\t"
                << mdotRight << std::endl;


    // create gas
    IdealGasMix gas(kinFileName, "gas");
    size_t nsp = gas.nSpecies();
    std::cout << "#\tNumber of species\t:\t" << nsp << std::endl;

    // set equivalence ratio
    vector_fp x(nsp, 0.0);
    vector_fp y(nsp, 0.0);
    doublereal mixtureFraction(0.0);
    doublereal ax = C_atoms + H_atoms / 4.0 - O_atoms / 2.0;
    doublereal fa_stoic = 1.0 / (4.76 * ax);
    for (size_t i=0;i<fuelName.size();i++)
    {
        x[gas.speciesIndex(fuelName[i])] = fuelX[i];
    }
    x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic;
    x[gas.speciesIndex("N2")] = 0.78 / phi / fa_stoic;
    x[gas.speciesIndex("AR")] = 0.01 / phi / fa_stoic;
    gas.setState_TPX(TInjection, p0, x.data());
    gas.getMassFractions(&y[0]);
    for (size_t i=0;i<fuelName.size();i++)
    {
        mixtureFraction += y[gas.speciesIndex(fuelName[i])];
    }
    std::cout << "#\tMixture fraction\t:\t" << mixtureFraction << std::endl;
    mdotInjection = mdotLeft*mixtureFraction / (1 - mixtureFraction);
    std::cout << "#\tmdotInjection\t\t:\t" << mdotInjection << std::endl;
    gas.equilibrate("HP");
    doublereal Tad = gas.temperature();
    doublereal yI[nsp];
    doublereal* yIn = &yI[0];
    gas.getMassFractions(yIn);


    // cloud
    Lagrangian cloud(diameterInjection, mdotInjection, TInjection, p0);
    cloud.setFuel(fuelName);
    std::stringstream ss1;
    ss1 << "LAGml-" << mdotL << "mr-" << mdotR << "L-" << len << "_raw.csv";
    cloud.outfile(ss1.str());

    // setup flow
    StFlow flow(&gas, nsp, points);
    flow.setAxisymmetricFlow();
    flow.setupCloud(cloud);
    cloud.flow(flow);
    doublereal zpoints[points];
    doublereal* z = &zpoints[0];
    doublereal dz = len / (points-1);
    for (int iz=0; iz<points; iz++) {
        z[iz] = iz*dz;
    }
    flow.setupGrid(points,z);
    Transport* tr = newTransportMgr("Mix", &gas);
    flow.setTransport(*tr);
    flow.setKinetics(gas);
    Inlet1D left;
    left.setMoleFractions(compLeft);
    left.setMdot(mdotLeft);
    left.setTemperature(TinLeft);
    Inlet1D right;
    right.setMoleFractions(compRight);
    right.setMdot(mdotRight);
    right.setTemperature(TinRight);
    std::vector<Domain1D*> domains;
    domains.push_back(&left);
    domains.push_back(&flow);
    domains.push_back(&right);
    Sim1D flame(domains);
    flame.setupCloud(cloud);
    std::cout << "Boundary conditions @ left side: " << std::endl;
    left.showSolution(0);
    std::cout << "Boundary conditions @ right side: " << std::endl;
    right.showSolution(0);


    // init
    if (ignition == false) Tad = TinLeft;
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
    values[1]=0.8*TinLeft+0.2*Tad;
    values[2]=0.6*Tad+0.4*TinLeft;
    values[3]=0.8*Tad+0.1*TinLeft+0.1*TinRight;
    values[4]=0.6*Tad+0.4*TinRight;
    values[5]=0.8*TinRight+0.2*Tad;
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
    flame.setMaxGridPoints(domFlow,1000);
    flame.setGridMin(domFlow, minGrid);
    flow.setPressure(p0);
    flow.solveEnergyEqn();
    flame.setRefineCriteria(domFlow,10.0,0.6,0.6,-0.1);
    flame.solve(1,true);


    // output
    std::stringstream ss2;
    ss2 << "ml-" << mdotL << "mr-" << mdotR << "L-" << len << "_raw.csv";
    std::ofstream outfile(ss2.str());
    int np = flow.nPoints();
    Array2D solution(np,nsp+6,0.0);
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
    }

    for (int j=0;j<solution.nColumns();j++)
    {
        if (j==0) outfile << "z";
        else outfile << flow.componentName(j-1);

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
    doublereal mdotL = 0.2;
    doublereal mdotR = 0.2;
    doublereal len = 0.02;

    if (argc == 4) {
        mdotL = atof(argv[1]);
        mdotR = atof(argv[2]);
        len = atof(argv[3]);
    }
    else {
        std::cerr << "Usage: counterflowSprayFlame mdotL mdotR len" << std::endl;
        return -1;
    }

    try {
        counterflowSprayFlame(mdotL, mdotR, len);
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
}
