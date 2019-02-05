#include <iostream>
#include <ostream>
#include <fstream>
#include <string>

#include "cantera/thermo.h"
#include "Inlet1DNew.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "Lagrangian.h"
#include "Lagrangian.cpp"
#include "Sim1DNew.h"
#include "Sim1DNew.cpp"
#include "StFlowNew.h"
#include "StFlowNew.cpp"
#include "OneDimNew.h"
#include "OneDimNew.cpp"

using namespace Cantera;

void counterflowDiffusionFlame()
{
	size_t tstepsSize = 6; // take tsteps[i] time steps
	const int tsteps[tstepsSize] = {20,50,100,200,100,500};
	doublereal p0=OneAtm;
	std::string kinFileName;
	std::string kinPhase;
	std::string fuelName;
	std::string compLeft;
	std::string compRight;
	std::string compInit;
	// input from file
	std::ifstream infile("input.txt");
	std::string name; // name column
	doublereal value; // value column
	size_t fuel=0; // 0 for Kerosene and 1 for nDodecane
	doublereal mdotInjection = 0.01; // kg/m^2/s
	doublereal diameterInjection = 100e-6; // m
	doublereal uInjection = 0.2; // m/s
	doublereal TInjection = 300; // K
	doublereal LagrangianRelaxationFactor = 1.0;
	doublereal TinLeft = 300; // fuel
	doublereal TinRight = 1500; // oxidizer
	doublereal mdotLeft = 0.2; // kg/m^2/s
	doublereal mdotRight = 0.2; // kg/m^2/s
	size_t points = 41;
	doublereal width = 0.02;
	doublereal a0 = 20;
	size_t ignition = 0; // 1 with ignition
	while (infile >> name)
	{
		infile >> value;
		if (name == "fuel")
			{ fuel = size_t(value); continue; }
		if (name == "points")
			{ points = size_t(value); continue; }
		if (name == "width")
			{ width = value; continue; }
		if (name == "mdotInjection")
			{ mdotInjection = value; continue; }
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
		if (name == "mdotLeft")
			{ mdotLeft = value; continue; }
		if (name == "mdotRight")
			{ mdotRight = value; continue; }
		if (name == "a0")
			{ a0 = value; continue; }
		if (name == "ignition")
			{ ignition = size_t(value); continue; }
	}
	#include "FuelInfo.h"
	std::cout << "#Input Parameters" << std::endl;
	std::cout << "#\tDomain      :\t"
				<< points << " "
				<< width << std::endl;
	std::cout << "#\tLagrangian  :\t"
				<< mdotInjection << " "
				<< diameterInjection << " "
				<< uInjection  << " "
				<< TInjection << std::endl;
	std::cout << "#\tTemperature :\t"
				<< TinLeft << " "
				<< TinRight << std::endl;
	std::cout << "#\tmdot        :\t"
				<< mdotLeft << " "
				<< mdotRight << std::endl;
	std::cout << std::endl;
	std::cout << "=============================================" << std::endl;


	// create gas, Lagrangian particle cloud and flow field
	IdealGasMix gas(kinFileName, kinPhase);
	size_t nsp = gas.nSpecies();
	std::cout << "\nNumber of species: #\t" << nsp << std::endl;
	std::ofstream LagrangianOutfile("LagrangianOutput.csv");
	std::ofstream& LagrangianOutfileRef = LagrangianOutfile;
	Lagrangian particleCloud(LagrangianOutfileRef, ignition, diameterInjection, mdotInjection);
	particleCloud.setInjectionVelocity(uInjection);
	particleCloud.setInjectionTemperature(TInjection);
	particleCloud.setRelaxationFractor(LagrangianRelaxationFactor);
	Lagrangian& particleCloudRef = particleCloud;
	StFlow flow(particleCloudRef, &gas, nsp, points);
	flow.setAxisymmetricFlow();
	flow.setFuelName(fuelName);


	// flow configuration
	doublereal zpoints[points];
	doublereal* z = &zpoints[0];
	doublereal dz = width / ((doublereal)(points-1));
	std::cout << "\nGrid points: #" << std::endl;
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
	right.setMoleFractions(compRight);
	right.setMdot(mdotRight);
	right.setTemperature(TinRight);
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
	gas.setState_TPX(0.5*(TinLeft+TinRight),p0,compInit);
	try {
		gas.equilibrate("HP");
	} catch (CanteraError& err) {
		std::cout << err.what() << std::endl;
	}
	doublereal Tad = TinLeft;
	if (ignition != 0) Tad = 0.9*gas.temperature();
	else Tad = 0.6*TinLeft+0.4*TinRight;
	doublereal yI[nsp];
	doublereal* yIn = &yI[0];
	gas.getMassFractions(yIn);
	vector_fp locs(5);
	vector_fp values(5);
	locs[0]=0;
	locs[1]=0.3;
	locs[2]=0.4;
	locs[3]=0.5;
	locs[4]=1.0;
	values[0]=TinLeft;
	values[1]=0.9*TinLeft+0.1*Tad;
	values[2]=0.8*Tad+0.1*TinLeft+0.1*TinRight;
	values[3]=0.9*TinRight+0.1*Tad;
	values[4]=TinRight;
	flame.setInitialGuess("T",locs,values);
	values[0]=a0*(0.5*width);
	values[1]=a0*(0.05*width);
	values[2]=a0*(-0.05*width);
	values[3]=a0*(-0.1*width);
	values[4]=a0*(-0.8*width);
	flame.setInitialGuess("u",locs,values);
	for (int ik=0; ik<nsp;ik++)
	{
       values[0]=left.massFraction(ik);
       values[1]=0.5*( left.massFraction(ik) + yIn[ik] );
       values[2]=yIn[ik];
       values[3]=0.5*( right.massFraction(ik) + yIn[ik] );
       values[4]=right.massFraction(ik);
       flame.setInitialGuess(gas.speciesName(ik), locs, values);
	}


	// solve
	size_t domFlow = 1;
	flame.setMaxGridPoints(domFlow,1000);
	flow.setPressure(p0);
	flow.solveEnergyEqn();
	flame.setRefineCriteria(domFlow,10,1,1,-1);
	flame.setTimeStep(1.0e-5,tstepsSize,&tsteps[0]);
	flame.solve(1,true);


	// output
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
	std::ofstream outfile("output.csv");
	for (int j=0;j<solution.nColumns();j++)
	{
		if (j==0) outfile << "z";
		else if (j<solution.nColumns()-1) outfile << flow.componentName(j-1);
		else outfile << "LiquidMass (kg/m^3)";
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



int main()
{
	try {
		counterflowDiffusionFlame();
		return 0;
	} catch (CanteraError& err) {
		std::cout << err.what() << std::endl;
		return 1;
	}
}
