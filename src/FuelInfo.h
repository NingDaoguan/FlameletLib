if (fuel == 0)
{
	std::cout << std::endl;
	std::cout << "=============================================" << std::endl;
	std::cout << " 1D Counterflow Diffusion Flame for Kerosene " << std::endl;
	std::cout << "                   T-E-S-T                   " << std::endl;
	std::cout << "=============================================" << std::endl;
	std::cout << std::endl;


	// input
	kinFileName = "2S_KEROSENE_BFER.cti";
	kinPhase = "gas";
	fuelName = "KERO_LUCHE";
	// compLeft = "KERO_LUCHE:1.0";
	compLeft =  "O2:0.21, N2:0.79";
	compRight = "O2:0.21, N2:0.79";
	compInit = "KERO_LUCHE:1.0, O2:10.0, N2:37.6";
}

if (fuel == 1)
{
	std::cout << std::endl;
	std::cout << "=============================================" << std::endl;
	std::cout << "1D Counterflow Diffusion Flame for n-Dodecane" << std::endl;
	std::cout << "=============================================" << std::endl;
	std::cout << std::endl;


	// input
	kinFileName = "nDodecane_CRECK.cti";
	kinPhase = "gas";
	fuelName = "NC12H26";
	// compLeft = "NC12H26:1.0";
	compLeft =  "O2:0.21, N2:0.78, AR:0.01";
	compRight = "O2:0.21, N2:0.78, AR:0.01";
	compInit = "NC12H26:1.0, O2:18.5, N2:69.56";
}
