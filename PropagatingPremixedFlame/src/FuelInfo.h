if (fuel == 0)
{
	std::cout << std::endl;
	std::cout << "=============================================================" << std::endl;
	std::cout << "      1D Freely-propagating Premixed Flame for KERO_BFE      " << std::endl;
	std::cout << "=============================================================" << std::endl;
	std::cout << "=============================================================" << std::endl;
	std::cout << "||                            |                            ||" << std::endl;
	std::cout << "||                            |                            ||" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "||                            |                            ||" << std::endl;
	std::cout << "||                            |                            ||" << std::endl;
	std::cout << "=============================================================" << std::endl;
	std::cout << "=============================================================" << std::endl;
	std::cout << std::endl;


	// input
	kinFileName = "2S_KEROSENE_BFER.cti";
	C_atoms = 10;
	H_atoms = 20;
	kinPhase = "gas";
	fuelName.push_back("KERO_LUCHE");
	fuelX.push_back(1.0);
}

if (fuel == 1)
{
	std::cout << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "         1D Freely-propagating Premixed Flame for n-Dodecane         " << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << std::endl;


	// input
	kinFileName = "nDodecane_CRECK.cti";
	C_atoms = 12;
	H_atoms = 26;
	kinPhase = "gas";
	fuelName.push_back("NC12H26");
	fuelX.push_back(1.0);
}

if (fuel == 2)
{
	std::cout << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "          1D Freely-propagating Premixed Flame for Kerosene          " << std::endl;
	std::cout << "                      Multi-component Surrogate                      " << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "||                                |                                ||" << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << "=====================================================================" << std::endl;
	std::cout << std::endl;

	// input
	kinFileName = "KEROSENE_CRECK231.cti";
	C_atoms = 12*0.3+16*0.36+10*0.246+7*0.094;
	H_atoms = 26*0.3+34*0.36+18*0.246+8*0.094;
	kinPhase = "gas";
	fuelName.push_back("NC12H26"); // component-0 n-Dodecane
	fuelName.push_back("IC16H34"); // component-1 Isocetane
	fuelName.push_back("DECALIN"); // component-2 Trans-decalin
	fuelName.push_back("C7H8"); // component-3 Toluene
	fuelX.push_back(0.300);
	fuelX.push_back(0.360);
	fuelX.push_back(0.246);
	fuelX.push_back(0.094);
}

