#include "cantera/thermo.h"
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>

using namespace Cantera;

void counterflowDiffusionFlame()
{
	std::cout << " 1D counterflow diffusion flame for CH4 " << std::endl;

	doublereal p0 = OneAtm;
	doublereal TinLeft = 300.0; // fuel
	doublereal TinRight = 300.0; // oxidizer
	doublereal mdotLeft = 0.24; // kg/m^2/s
	doublereal mdotRight = 0.72; // kg/m^2/s

	std::string compLeft = "CH4:0.999, AR:0.001";
	std::string compRight = "O2:0.21, N2:0.79";

	IdealGasMix gas("gri30.cti");
	size_t nsp = gas.nSpecies();
	
	size_t points = 21;
	doublereal width = 0.02;


	// create the flow
	StFlow flow(&gas,nsp,points);
	flow.setAxisymmetricFlow();


	// setup grid
	doublereal zpoints[points];
	doublereal* z = &zpoints[0];
	doublereal dz = width / ((doublereal)(points-1));

	for (int iz=0;iz<points;iz++)
	{
		z[iz] = ((doublereal)iz)*dz;
		std::cout << z[iz] << " ";
		if (iz == points-1)
			std::cout << std::endl;
	}

	flow.setupGrid(points,z);


	// transport model
	Transport* tr = newTransportMgr("Mix",&gas);
	flow.setTransport(*tr);
	flow.setKinetics(gas);


	// inlets
	Inlet1D left;
	left.setMoleFractions(compLeft);
	left.setMdot(mdotLeft);
	left.setTemperature(TinLeft);

	Inlet1D right;
	right.setMoleFractions(compRight);
	right.setMdot(mdotRight);
	right.setTemperature(TinRight);


	// domains
	std::vector<Domain1D*> domains;
	domains.push_back(&left);
	domains.push_back(&flow);
	domains.push_back(&right);

	Sim1D flame(domains);
	std::cout << "Boundary conditions @ left side: " << std::endl;
	left.showSolution(0);
	std::cout << "Boundary conditions @ right side: " << std::endl;
	right.showSolution(0);


	// init
	gas.setState_TPX(0.5*(TinLeft+TinRight),p0,"CH4:1.0, O2:2.0, N2:7.52");
	try {
		gas.equilibrate("HP");
	} catch (CanteraError& err) {
		std::cout << err.what() << std::endl;
	}
	doublereal Tad = gas.temperature();
	doublereal yI[nsp];
	doublereal* yIn = &yI[0];
	gas.getMassFractions(yIn);

	vector_fp locs;
	vector_fp values;
	locs.resize(3);
	values.resize(3);
	locs[0]=0;
	locs[1]=0;
	locs[2]=1.0;
	values[0]=TinLeft;
	values[1]=Tad;
	values[2]=TinRight;
	flame.setInitialGuess("T",locs,values);

	for (int ik=0; ik<nsp;ik++)
	{
       values[0]=left.massFraction(ik);
       values[1]=yIn[ik];
       values[2]=right.massFraction(ik);
       flame.setInitialGuess(gas.speciesName(ik), locs, values);
	}


	// solve
	size_t domFlow = 1;
	flame.setMaxGridPoints(domFlow,200);
	flow.setPressure(p0);
	flow.solveEnergyEqn();
	flame.setRefineCriteria(domFlow,200.0,0.1,0.2,0.02);
	flame.solve(1,true);


	// post-process
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
 
		for (int i=5; i<nsp+5;i++)
	   	{
	  		if (flame.value(domFlow, i-1,n) > 0.0) {solution(n,i) = flame.value(domFlow, i-1,n);}
	    	else if (flame.value(domFlow, i-1,n) < 0.0)
	    	{
	        	solution(n,i) = 0.0;
	        	flame.setValue(domFlow, i-1,n,0.0);
	    	}
		}		
	}

	std::ofstream outfile("output.txt");
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