if (fuel == 0)
{
    std::cout << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << " 1D Counterflow Diffusion Flame for Kerosene " << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "||                    |                    ||" << std::endl;
    std::cout << "||                    |                    ||" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "||                    |                    ||" << std::endl;
    std::cout << "||                    |                    ||" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << std::endl;


    // input
    kinFileName = "2S_KEROSENE_BFER.cti";
    kinPhase = "gas";
    C_atoms = 10;
    H_atoms = 20;
    O_atoms = 0;
    Tvapor = 480;
    fuelName.push_back("KERO_LUCHE");
    compLeft =  "O2:0.21, N2:0.79";
    compRight = "O2:0.21, N2:0.79";
    fuelX.push_back(1.0);
}

else if (fuel == 1)
{
    std::cout << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "            1D Counterflow Diffusion Flame for n-Dodecane            " << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << std::endl;


    // input
    kinFileName = "nDodecane_CRECK.cti";
    kinPhase = "gas";
    C_atoms = 12;
    H_atoms = 26;
    O_atoms = 0;
    Tvapor = 480;
    fuelName.push_back("NC12H26");
    compLeft =  "O2:0.21, N2:0.78, AR:0.01";
    compRight = "O2:0.21, N2:0.78, AR:0.01";
    fuelX.push_back(1.0);
}

else if (fuel == 2)
{
    std::cout << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "             1D Counterflow Diffusion Flame for Kerosene             " << std::endl;
    std::cout << "                      Multi-component Surrogate                      " << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << std::endl;

    // input
    kinFileName = "KEROSENE_CRECK231.cti";
    kinPhase = "gas";
    C_atoms = 12*0.3+16*0.36+10*0.246+7*0.094;
    H_atoms = 26*0.3+34*0.36+18*0.246+8*0.094;
    O_atoms = 0;
    Tvapor = 480;
    fuelName.push_back("NC12H26"); // component-0 n-Dodecane
    fuelName.push_back("IC16H34"); // component-1 Isocetane
    fuelName.push_back("DECALIN"); // component-2 Trans-decalin
    fuelName.push_back("C7H8"); // component-3 Toluene
    compLeft  = "O2:0.21, N2:0.78, AR:0.01";
    compRight = "O2:0.21, N2:0.78, AR:0.01";
    fuelX.push_back(0.300);
    fuelX.push_back(0.360);
    fuelX.push_back(0.246);
    fuelX.push_back(0.094);
}

else if (fuel == 3)
{
    std::cout << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "              1D Counterflow Diffusion Flame for C2H5OH              " << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>|<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=<=" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "||                                |                                ||" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << std::endl;

    // input
    kinFileName = "Ethanol_31.cti";
    kinPhase = "gas";
    C_atoms = 2;
    H_atoms = 6;
    O_atoms = 1;
    Tvapor = 350;
    fuelName.push_back("C2H5OH");
    compLeft  = "O2:0.21, N2:0.78, AR:0.01";
    compRight = "O2:0.21, N2:0.78, AR:0.01";
    fuelX.push_back(1.0);
}