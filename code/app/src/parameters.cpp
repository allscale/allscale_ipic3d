#include "ipic3d/parameters.h"

#include <array>
#include <cmath>

#include "ipic3d/ConfigFile.h"

namespace ipic3d {

	Parameters Parameters::read(std::string inputfile) {
		Parameters res;
		res.readInternal(inputfile);
		return res;
	}

	void Parameters::readInternal(std::string inputfile) {
		// Loading the input file
		ConfigFile config(inputfile);

			// the following variables are ALWAYS taken from inputfile
		dt = config.read < double >("dt");
		ncycles = config.read < int >("ncycles");
		th = config.read < double >("th",1.0);

		Smooth = config.read < double >("Smooth",1.0);
		SmoothNiter = config.read < int >("SmoothNiter",6);

		SaveDirName = config.read < std::string > ("SaveDirName","data");
		RestartDirName = config.read < std::string > ("RestartDirName","data");
		ns = config.read < int >("ns");
		nstestpart = config.read < int >("nsTestPart", 0);
		NpMaxNpRatio = config.read < double >("NpMaxNpRatio",1.5);
		// mode parameters for second order in time
		PushWithBatTime = config.read < double >("PushWithBatTime",0);
		PushWithEatTime = config.read < double >("PushWithEatTime",1);
		ImplSusceptTime = config.read < double >("ImplSusceptTime",0);
		ImplSusceptMode = 0;
		// GEM Challenge
		B0x = config.read <double>("B0x",0.0);
		B0y = config.read <double>("B0y",0.0);
		B0z = config.read <double>("B0z",0.0);

		// Earth parameters
		B1x = 0.0;
		B1y = 0.0;
		B1z = 0.0;
		B1x = config.read <double>("B1x",0.0);
		B1y = config.read <double>("B1y",0.0);
		B1z = config.read <double>("B1z",0.0);

		delta = config.read < double >("delta",0.5);

		Case              = config.read<std::string>("Case");
		wmethod           = config.read<std::string>("WriteMethod");
		SimName           = config.read<std::string>("SimulationName");
		PoissonCorrection = config.read<std::string>("PoissonCorrection");
		PoissonCorrectionCycle = config.read<int>("PoissonCorrectionCycle",10);

		using array_double = std::array<double,12>;

		rhoINIT = new double[12];
		array_double rhoINIT0 = config.read < array_double > ("rhoINIT");
		for(int i = 0; i < 12; i++) {
			rhoINIT[i] = rhoINIT0[i];
		}

		rhoINJECT = new double[ns];
		array_double rhoINJECT0 = config.read<array_double>( "rhoINJECT" );
		for(int i = 0; i < ns; i++) {
			rhoINJECT[i] = rhoINJECT0[i];
		}

		// take the tolerance of the solvers
		CGtol = config.read < double >("CGtol",1e-3);
		GMREStol = config.read < double >("GMREStol",1e-3);
		NiterMover = config.read < int >("NiterMover",3);
		// take the injection of the particless
		Vinj = config.read < double >("Vinj",0.0);

		// take the output cycles
		FieldOutputCycle = config.read < int >("FieldOutputCycle",100);
		ParticlesOutputCycle = config.read < int >("ParticlesOutputCycle",0);
		FieldOutputTag     =   config.read <std::string>("FieldOutputTag","");
		ParticlesOutputTag =   config.read <std::string>("ParticlesOutputTag","");
		MomentsOutputTag   =   config.read <std::string>("MomentsOutputTag","");
		TestParticlesOutputCycle = config.read < int >("TestPartOutputCycle",0);
		testPartFlushCycle = config.read < int >("TestParticlesOutputCycle",10);
		RestartOutputCycle = config.read < int >("RestartOutputCycle",5000);
		DiagnosticsOutputCycle = config.read < int >("DiagnosticsOutputCycle", FieldOutputCycle);
		CallFinalize = config.read < bool >("CallFinalize", true);

		restart_status = 0;
		last_cycle = -1;
		c = config.read < double >("c",1.0);

		Lx = config.read < double >("Lx",10.0);
		Ly = config.read < double >("Ly",10.0);
		Lz = config.read < double >("Lz",10.0);
		nxc = config.read < int >("nxc",64);
		nyc = config.read < int >("nyc",64);
		nzc = config.read < int >("nzc",64);
		XLEN = config.read < int >("XLEN",1);
		YLEN = config.read < int >("YLEN",1);
		ZLEN = config.read < int >("ZLEN",1);
		PERIODICX = config.read < bool >("PERIODICX",true);
		PERIODICY = config.read < bool >("PERIODICY",true);
		PERIODICZ = config.read < bool >("PERIODICZ",true);

		PERIODICX_P = config.read < bool >("PERIODICX_P",PERIODICX);
		PERIODICY_P = config.read < bool >("PERIODICY_P",PERIODICY);
		PERIODICZ_P = config.read < bool >("PERIODICZ_P",PERIODICZ);

		x_center = config.read < double >("x_center",5.0);
		y_center = config.read < double >("y_center",5.0);
		z_center = config.read < double >("z_center",5.0);
		L_square = config.read < double >("L_square",5.0);


		uth = new double[ns];
		vth = new double[ns];
		wth = new double[ns];
		u0 = new double[ns];
		v0 = new double[ns];
		w0 = new double[ns];
		array_double uth0 = config.read < array_double > ("uth");
		array_double vth0 = config.read < array_double > ("vth");
		array_double wth0 = config.read < array_double > ("wth");
		array_double u00 = config.read < array_double > ("u0");
		array_double v00 = config.read < array_double > ("v0");
		array_double w00 = config.read < array_double > ("w0");

		for(int i = 0; i<ns; i++) {
			uth[i] = uth0[i];
			vth[i] = vth0[i];
			wth[i] = wth0[i];
			u0[i] = u00[i];
			v0[i] = v00[i];
			w0[i] = w00[i];
		}


		if (nstestpart > 0) {
			pitch_angle = new double[nstestpart];
			energy      = new double[nstestpart];
			array_double pitch_angle0 = config.read < array_double > ("pitch_angle");
			array_double energy0 	  = config.read < array_double > ("energy");

			for(int i=0; i<nstestpart; i++) {
				pitch_angle[i] = pitch_angle0[i];
				energy[i] = energy0[i];
			}

		}

		npcelx = new int[ns+nstestpart];
		npcely = new int[ns+nstestpart];
		npcelz = new int[ns+nstestpart];
		qom = new double[ns+nstestpart];

		using array_int = std::array<int,12>;

		array_int npcelx0 = config.read < array_int > ("npcelx");
		array_int npcely0 = config.read < array_int > ("npcely");
		array_int npcelz0 = config.read < array_int > ("npcelz");
		array_double qom0 = config.read < array_double > ("qom");

		for(int i=0; i<ns+nstestpart; ++i) {
			npcelx[i] = npcelx0[i];
			npcely[i] = npcely0[i];
			npcelz[i] = npcelz0[i];
			qom[i]	= qom0[i];
		}

		// PHI Electrostatic Potential
		bcPHIfaceXright = config.read < int >("bcPHIfaceXright",1);
		bcPHIfaceXleft  = config.read < int >("bcPHIfaceXleft",1);
		bcPHIfaceYright = config.read < int >("bcPHIfaceYright",1);
		bcPHIfaceYleft  = config.read < int >("bcPHIfaceYleft",1);
		bcPHIfaceZright = config.read < int >("bcPHIfaceZright",1);
		bcPHIfaceZleft  = config.read < int >("bcPHIfaceZleft",1);

		// EM field boundary condition
		bcEMfaceXright = config.read < int >("bcEMfaceXright");
		bcEMfaceXleft  = config.read < int >("bcEMfaceXleft");
		bcEMfaceYright = config.read < int >("bcEMfaceYright");
		bcEMfaceYleft  = config.read < int >("bcEMfaceYleft");
		bcEMfaceZright = config.read < int >("bcEMfaceZright");
		bcEMfaceZleft  = config.read < int >("bcEMfaceZleft");

		/*  ---------------------------------------------------------- */
		/*  Electric and Magnetic field boundary conditions for BCface */
		/*  ---------------------------------------------------------- */
		// if bcEM* is 0: perfect conductor, if bcEM* is not 0: perfect mirror
		// perfect conductor: normal = free, perpendicular = 0
		// perfect mirror   : normal = 0,    perpendicular = free
		/*  ---------------------------------------------------------- */

		/* X component in faces Xright, Xleft, Yright, Yleft, Zright and Zleft (0, 1, 2, 3, 4, 5) */
		bcEx[0] = bcEMfaceXright == 0 ? 2 : 1;   bcBx[0] = bcEMfaceXright == 0 ? 1 : 2;
		bcEx[1] = bcEMfaceXleft  == 0 ? 2 : 1;   bcBx[1] = bcEMfaceXleft  == 0 ? 1 : 2;
		bcEx[2] = bcEMfaceYright == 0 ? 1 : 2;   bcBx[2] = bcEMfaceYright == 0 ? 2 : 1;
		bcEx[3] = bcEMfaceYleft  == 0 ? 1 : 2;   bcBx[3] = bcEMfaceYleft  == 0 ? 2 : 1;
		bcEx[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBx[4] = bcEMfaceZright == 0 ? 2 : 1;
		bcEx[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBx[5] = bcEMfaceZleft  == 0 ? 2 : 1;
		/* Y component */
		bcEy[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBy[0] = bcEMfaceXright == 0 ? 2 : 1;
		bcEy[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBy[1] = bcEMfaceXleft  == 0 ? 2 : 1;
		bcEy[2] = bcEMfaceYright == 0 ? 2 : 1;   bcBy[2] = bcEMfaceYright == 0 ? 1 : 2;
		bcEy[3] = bcEMfaceYleft  == 0 ? 2 : 1;   bcBy[3] = bcEMfaceYleft  == 0 ? 1 : 2;
		bcEy[4] = bcEMfaceZright == 0 ? 1 : 2;   bcBy[4] = bcEMfaceZright == 0 ? 2 : 1;
		bcEy[5] = bcEMfaceZleft  == 0 ? 1 : 2;   bcBy[5] = bcEMfaceZleft  == 0 ? 2 : 1;
		/* Z component */
		bcEz[0] = bcEMfaceXright == 0 ? 1 : 2;   bcBz[0] = bcEMfaceXright == 0 ? 2 : 1;
		bcEz[1] = bcEMfaceXleft  == 0 ? 1 : 2;   bcBz[1] = bcEMfaceXleft  == 0 ? 2 : 1;
		bcEz[2] = bcEMfaceYright == 0 ? 1 : 1;   bcBz[2] = bcEMfaceYright == 0 ? 2 : 1;
		bcEz[3] = bcEMfaceYleft  == 0 ? 1 : 1;   bcBz[3] = bcEMfaceYleft  == 0 ? 2 : 1;
		bcEz[4] = bcEMfaceZright == 0 ? 2 : 1;   bcBz[4] = bcEMfaceZright == 0 ? 1 : 2;
		bcEz[5] = bcEMfaceZleft  == 0 ? 2 : 1;   bcBz[5] = bcEMfaceZleft  == 0 ? 1 : 2;

		// Particles Boundary condition
		bcPfaceXright = config.read < int >("bcPfaceXright",1);
		bcPfaceXleft  = config.read < int >("bcPfaceXleft",1);
		bcPfaceYright = config.read < int >("bcPfaceYright",1);
		bcPfaceYleft  = config.read < int >("bcPfaceYleft",1);
		bcPfaceZright = config.read < int >("bcPfaceZright",1);
		bcPfaceZleft  = config.read < int >("bcPfaceZleft",1);

		// initialize some of the derived parameters
		/*! fourpi = 4 greek pi */
		fourpi = 16.0 * atan(1.0);
		/*! dx = space step - X direction */
		dx = Lx / (double) nxc;
		/*! dy = space step - Y direction */
		dy = Ly / (double) nyc;
		/*! dz = space step - Z direction */
		dz = Lz / (double) nzc;
		/*! npcel = number of particles per cell */
		npcel = new int[ns+nstestpart];
		for (int i = 0; i < (ns+nstestpart); i++) {
			npcel[i] = npcelx[i] * npcely[i] * npcelz[i];
		}
	}



} // end namespace ipic3d
