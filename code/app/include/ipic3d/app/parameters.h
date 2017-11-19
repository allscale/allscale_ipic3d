#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <regex>

#include "ipic3d/app/vector.h"

namespace ipic3d {

    /**
    * Our data types for double and int arrays
    */
    using dvector = std::vector<double>;
    using ivector = std::vector<int>;

	/**
	* An enumeration of use cases.
	*/
    enum class UseCase { Dipole, Test };

	struct Parameters {

		// light speed
		double c = 1.0;

		// time step
		double dt;

		// number of time cycles
		int ncycles;

		// simulation box length per direction
		Vector3<double> L;

		// object center per direction
		Vector3<double> objectCenter;
		
		// planet radius, assuming a cubic box
		double planetRadius;

		// number of cells per direction of problem domain
		Vector3<int> ncells;

		// grid spacing per direction
		Vector3<double> dspace;

		// number of species
		int ns;

		// number of particles per cell per direction
		ivector npcelx;
		ivector npcely;
		ivector npcelz;

		// charge to mass ratio array for different species
		dvector qom;

		// charge to mass ratio array for different species
		dvector rhoInit;
		// thermal velocity per direction
		dvector uth;
		dvector vth;
		dvector wth;
		// Drift velocity per direction
		dvector u0;
		dvector v0;
		dvector w0;


		// Case type
		UseCase useCase;

		// Output writing method
		std::string wmethod;

		// Simulation name
		std::string SimName;

		// Poisson correction flag
		std::string PoissonCorrection;

		// SaveDirName
		std::string SaveDirName;


		// GEM Challenge parameters
		// current sheet thickness
		double delta;
		// Amplitude of the field
		Vector3<double> B0;
		Vector3<double> B1;


		// Output for field
		int FieldOutputCycle;
		std::string  FieldOutputTag;
		std::string  MomentsOutputTag;

		// Output for particles
		int ParticlesOutputCycle;
		std::string ParticlesOutputTag;


		// TODO: the boundary conditions can be used in the alternative field solvers that are planned for the third year of the project
		// Boundary condition selection for BCFace for the electric field components
		ivector bcEx, bcEy, bcEz;
		// Boundary condition selection for BCFace for the magnetic field components
		ivector bcBx, bcBy, bcBz;

		// Boundary condition on particles 0 = exit 1 = perfect mirror 2 = riemission
		// Boundary Condition Particles: FaceRight and FaceLeft
		Vector3<int> bcPfaceRight;
		Vector3<int> bcPfaceLeft;


		// Field Boundary Condition 0 = Dirichlet Boundary Condition: specifies the valueto take pn the boundary of the domain 1 = Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain 2 = Periodic Condition
		// Boundary Condition Electrostatic Potential: PhiFaceRight and PhiFaceLeft
		Vector3<int> bcPHIfaceRight;
		Vector3<int> bcPHIfaceLeft;

		// Boundary Condition EM Field: EMFaceRight and EMFaceLeft
		Vector3<int> bcEMfaceRight;
		Vector3<int> bcEMfaceLeft;


		/**
		* Auxiliary function to split the line
		*/ 
		std::vector<std::string> split(const std::string & s, std::string delim = "\\s+") {
			std::vector<std::string> elems;

			std::regex rgx(delim);

			std::sregex_token_iterator iter(s.begin(), s.end(), rgx, -1);
			std::sregex_token_iterator end;
			for (; iter != end; ++iter) {
				elems.push_back(*iter);
			}

			return elems;
		}

		Parameters(std::string inputfile) {
			// open the config file and read all the inputs one by one
			std::ifstream in(inputfile);
			if(!in) {
				assert_fail() << "File not found: " << inputfile;
				exit(EXIT_FAILURE);
			}

			// read the input file line by line and parse parameters	
			std::string str;
			while( !std::getline(in, str).eof() ) {
				
				// ignore comments
				static const std::string comment = "#";
				str = str.substr( 0, str.find(comment) );

				// ignore those lines that do not follow the pattern
				static const std::string eqsign = "=";
				if ( str.find( eqsign ) == std::string::npos ) {
					continue;
				}

				if ( str.find("dt") != std::string::npos ) {
					dt = std::stod( split(str).back() );
					continue;
				}

				if ( str.find("ncycles") != std::string::npos ) {
					ncycles = std::stoi( split(str).back() );
					continue;
				}

				if ( str.find("Lx") != std::string::npos ) {
					L.x = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("Ly") != std::string::npos ) {
					L.y = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("Lz") != std::string::npos ) {
					L.z = std::stod( split(str).back() );
					continue;
				}

				if ( str.find("x_center") != std::string::npos ) {
					objectCenter.x = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("y_center") != std::string::npos ) {
					objectCenter.y = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("z_center") != std::string::npos ) {
					objectCenter.z = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("L_square") != std::string::npos ) {
					planetRadius = std::stod( split(str).back() );
					continue;
				}

				if ( str.find("delta") != std::string::npos ) {
					delta = std::stod( split(str).back() );
					continue;
				}

				if ( str.find("nxc") != std::string::npos ) {
					ncells.x = std::stoi( split(str).back() );
					continue;
				}
				if ( str.find("nyc") != std::string::npos ) {
					ncells.y = std::stoi( split(str).back() );
					continue;
				}
				if ( str.find("nzc") != std::string::npos ) {
					ncells.z = std::stoi( split(str).back() );
					continue;
				}

				if ( str.find("B0x") != std::string::npos ) {
					B0.x = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("B0y") != std::string::npos ) {
					B0.y = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("B0z") != std::string::npos ) {
					B0.z = std::stod( split(str).back() );
					continue;
				}

				if ( str.find("B1x") != std::string::npos ) {
					B1.x = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("B1y") != std::string::npos ) {
					B1.y = std::stod( split(str).back() );
					continue;
				}
				if ( str.find("B1z") != std::string::npos ) {
					B1.z = std::stod( split(str).back() );
					continue;
				}

				if ( str.find("Case") != std::string::npos ) {
					if ( split(str).back().compare("Dipole") == 0 )
						useCase = UseCase::Dipole;
					else 
						useCase = UseCase::Test;
					continue;
				}
				if ( str.find("SaveDirName") != std::string::npos ) {
					SaveDirName = split(str).back();
					continue;
				}

				if ( str.find("PoissonCorrection") != std::string::npos ) {
					PoissonCorrection = split(str).back();
					continue;
				}

				if ( str.find("WriteMethod") != std::string::npos ) {
					wmethod = split(str).back();
					continue;
				}

				if ( str.find("SimulationName") != std::string::npos ) {
					SimName = split(str).back();
					continue;
				}

				if ( str.find("ns = ") != std::string::npos ) {
					ns = std::stoi( split(str).back() );
					continue;
				}

				if ( str.find("npcelx") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						npcelx.push_back( std::stoi(*it) );
					continue;
				}
				if ( str.find("npcely") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						npcely.push_back( std::stoi(*it) );
					continue;
				}
				if ( str.find("npcelz") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						npcelz.push_back( std::stoi(*it) );
					continue;
				}

				if ( str.find("qom") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						qom.push_back( std::stod(*it) );
					continue;
				}

				if ( str.find("rhoINIT") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						rhoInit.push_back( std::stod(*it) );
					continue;
				}

				if ( str.find("uth") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						uth.push_back( std::stod(*it) );
					continue;
				}
				if ( str.find("vth") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						vth.push_back( std::stod(*it) );
					continue;
				}
				if ( str.find("wth") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						wth.push_back( std::stod(*it) );
					continue;
				}

				if ( str.find("u0") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						u0.push_back( std::stod(*it) );
					continue;
				}
				if ( str.find("v0") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						v0.push_back( std::stod(*it) );
					continue;
				}
				if ( str.find("w0") != std::string::npos ) {
					std::vector<std::string> vec = split(str);
					for (auto it = vec.cbegin() + 2; it != vec.cend(); ++it)
						w0.push_back( std::stod(*it) );
					continue;
				}

				if ( str.find("FieldOutputCycle") != std::string::npos ) {
					FieldOutputCycle = std::stoi( split(str).back() );
					continue;
				}
				if ( str.find("FieldOutputTag") != std::string::npos ) {
					FieldOutputTag = split(str).back();
					continue;
				}
				if ( str.find("MomentsOutputTag") != std::string::npos ) {
					MomentsOutputTag = split(str).back();
					continue;
				}

				if ( str.find("ParticlesOutputCycle") != std::string::npos ) {
					ParticlesOutputCycle = std::stoi( split(str).back() );
					continue;
				}
				if ( str.find("ParticlesOutputTag") != std::string::npos ) {
					ParticlesOutputTag = split(str).back();
					continue;
				}

			}
			
			dspace.x = L.x / ncells.x;
			dspace.y = L.y / ncells.y;
			dspace.z = L.z / ncells.z;

			in.close();
		}

	};

} // end namespace ipic3d
