

#include "sp_band.h"

Sp_band::Sp_band() {

	this ->sp_band_mass = sp_band_mass_const * m_e_au;
   	this ->unit_cell_volume = unit_volume_au;
    	this ->n_sp_density = number_of_sp_electrons / unit_cell_volume;
    	this ->fermi_energy = fermi_energy_const *eV_to_Hartree;

}

Sp_band::Sp_band(int number_disc_points) {

	this ->number_disc_points = number_disc_points;
}

Sp_band::~Sp_band() {
}

double Sp_band::Get_sp_BandMass() {

    	return sp_band_mass;
}

int Sp_band::GetNumberDiscPoints() {

    	return number_disc_points;
}

double Sp_band::Get_n_sp_Density() {

    	return n_sp_density;
}

double Sp_band::GetUnitVolume() {

    	return unit_cell_volume;
}

double Sp_band::GetFermiEnergy() {

    	return fermi_energy;
}


//double Sp_band::get_sp_DOS(double energy){
// 
//    
//    return 2.0 * pow((2.0 * sp_band_mass), 1.5) / (4.0 * pow(PI, 2.0) * pow(hbar, 3.0)) * std::sqrt(energy);
// 
//}
//

double Sp_band::Get_sp_DOS(double energy) {

     	return (std::sqrt(2.0 * (energy - 0.0)) * std::pow(Get_sp_BandMass(), 1.5)) / (PI * PI * std::pow(hbar_au, 3.0));	 
}

double Sp_band::GetDispersionRelation(double wavevector) {

    	return hbar_au * hbar_au * wavevector * wavevector / (2.0 *Get_sp_BandMass());
}

double Sp_band::Get_sp_DOSWavevector(double wavevector) {

    	return Get_sp_DOS(GetDispersionRelation(wavevector));
}

double Sp_band::GetWavevector(double energy) {

    	return std::sqrt(2.0 * Get_sp_BandMass() * GetDispersionRelation(energy) * GetDispersionRelation(energy) / hbar_au);
}

void Sp_band::Read() {

    	std::string line;
    	int size = 12000;
    	double E, Dos;
    	std::vector<double> energy;
    	std::vector<double> dos;
	//double E_dos[size];
    	double s_dos[size];
    	double p_dos[size];
    	double sp_dos[size];
    	double d_dos[size];
    	std::fstream my_file("/home/ndione/NetBeansProjects/dielectric_function/example/dielectric/PDOS_S01_A0001.OUT");

    	while (my_file >> E >> Dos) {

        	if (!my_file.good()) break;

        		energy.push_back(E);
        		dos.push_back(Dos);

    	}
	
    	my_file.close();

	//    std::ofstream file("/home/ndione/NetBeansProjects/dielectric_function/example/dielectric/DOS_seb.OUT");
    	std::ofstream file("/home/ndione/Desktop/Root_finding_1D/DOS/DOS_states_Hartree_atom.OUT");
    	for (int i = 0; i < size; i++) {

        	s_dos[i] = 2.0 * dos[i + 0.0 * 12001];

        	p_dos[i] = 2.0 * (dos[i + 1.0 * 12001] + dos[i + 2.0 * 12001] + dos[i + 3.0 * 12001]);

        	sp_dos[i] = s_dos[i] + p_dos[i];

        	d_dos[i] = 2.0 * (dos [i + 4.0 * 12001] + dos [i + 5.0 * 12001] + dos [i + 6.0 * 12001] + dos [i + 7.0 * 12001] + dos [i + 8.0 * 12001]);

                std::cout<<"E= " << energy[i] << "\t \t" << "s_DOS= " << s_dos[i] << "\t" << "p_DOS= " << p_dos[i]
                        << "\t" << "sp_DOS= " << s_dos[i] + p_dos[i] << "\t"
                        << "d_DOS= " << d_dos[i] << std::endl;
                
                file<<std::scientific<<std::setprecision(15)
                        << energy[i]///eV_to_Hartree //+ 2.179872170950000e-18 
                        << " \t "
        	//                << s_dos[i]/Hartree_to_eV/unit_volume_SI
       	 	//                << " \t "
        	//                <<p_dos[i]
        	//                << " \t "
		//                        <<sp_dos[i]*eV_to_Hartree
		//                        <<" \t "
		//                        <<d_dos[i]*eV_to_Hartree 
		//                        <<" \t "
                        <<(sp_dos[i] + d_dos[i])//*eV_to_Hartree 
                        << std::endl; 
    	}


    	file.close();

}

void Sp_band::Readfile_s_dos(std::vector<double>* energy, std::vector<double>* s_dos, const char* filename) {

    	std::fstream my_file(filename);
    	double value1;
    	double value2;

    	for (int i = 0; i < 12000; i++) {

        	my_file >> value1 >> value2;
        	
		if (!my_file.good()) break;

        		energy->push_back(value1);
        		s_dos->push_back(value2);

    	}
    	my_file.close();
}

void Sp_band::Readfile_p_dos(std::vector<double>* energy, std::vector<double>* p_dos, const char* filename) {


    	std::fstream my_file(filename);
    	double value1;
    	double value2;

    	while (true) {

        	my_file >> value1 >> value2;
        	if (!my_file.good()) break;

        		energy->push_back(value1);
        		p_dos->push_back(value2);

    	}
    	my_file.close();
}

void::Sp_band::Readfile_sp_dos(std::vector<double>* energy, std::vector<double>* sp_dos, const char* filename) {

    exit(1);
}
