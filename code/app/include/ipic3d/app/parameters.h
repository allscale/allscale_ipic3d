#pragma once

#include <string>

namespace ipic3d {

	/**
	 * An enumeration of use cases.
	 */
	enum class UseCase {
		ParticleWave,
		Dipole,
		Test
	};


	struct Parameters {

		/*! light speed */
		double c;

		/*! 4 pi */
		double fourpi;

		/*! time step */
		double dt;

		//
		// parameters used to support second order accuracy in time
		//
		/*! decentering parameter */
		double th; // second-order for th=1/2, stable for 1/2 <= th <= 1
		/*! time of magnetic field used in particle push (0=initial, 1=final) */
		double PushWithBatTime; // 0=initial (default), 1=final
		/*! time of electric field used in particle push */
		double PushWithEatTime; // 0=initial, 1=final (default)
		/*! means of estimating time-advanced implicit susceptibility */
		int ImplSusceptMode; // 0 - "initial" (default), "explPredict", "implPredict"
		/*! time of implicit susceptibility used in field advance */
		double ImplSusceptTime; // 0=initial (default), 1=final
		//
		/*! Smoothing value */
		double Smooth;
		int SmoothNiter;
		/*! number of time cycles */
		int ncycles;
		/*! physical space dimensions */
		int dim;
		/*! simulation box length - X direction */
		double Lx;
		/*! simulation box length - Y direction */
		double Ly;
		/*! simulation box length - Z direction */
		double Lz;
		/*! object center - X direction */
		double x_center;
		/*! object center - Y direction */
		double y_center;
		/*! object center - Z direction */
		double z_center;
		/*! object size - assuming a cubic box */
		double L_square;
		// number of cells per direction of problem domain
		int nxc;
		int nyc;
		int nzc;
		/*! grid spacing - X direction */
		double dx;
		/*! grid spacing - Y direction */
		double dy;
		/*! grid spacing - Z direction */
		double dz;
		/*! number of MPI subdomains in each direction */
		int XLEN;
		int YLEN;
		int ZLEN;
		/*! periodicity in each direction */
		bool PERIODICX;
		bool PERIODICY;
		bool PERIODICZ;
		/*! Particle periodicity in each direction */
		bool PERIODICX_P;
		bool PERIODICY_P;
		bool PERIODICZ_P;

		/*! number of species */
		int ns;
		/*! number of test particle species */
		int nstestpart;
		/*! number of particles per cell */
		int *npcel;
		/*! number of particles per cell - X direction */
		int *npcelx;
		/*! number of particles per cell - Y direction */
		int *npcely;
		/*! number of particles per cell - Z direction */
		int *npcelz;
		// either make these of longid type or do not declare them.
		//int *np; /*! number of particles array for different species */
		//int *npMax; /*! maximum number of particles array for different species */
		/*! max number of particles */
		double NpMaxNpRatio;
		/*! charge to mass ratio array for different species */
		double *qom;
		/*! charge to mass ratio array for different species */
		double *rhoINIT;
		/*! density of injection */
		double *rhoINJECT;
		/*! thermal velocity - Direction X */
		double *uth;
		/*! thermal velocity - Direction Y */
		double *vth;
		/*! thermal velocity - Direction Z */
		double *wth;
		/*! Drift velocity - Direction X */
		double *u0;
		/*! Drift velocity - Direction Y */
		double *v0;
		/*! Drift velocity - Direction Z */
		double *w0;

		/*! Pitch Angle for Test Particles */
		double *pitch_angle;
		/*! Energy for Test Particles */
		double *energy;


		/*! Case type */
		UseCase useCase;
		/*! Output writing method */
		std::string wmethod;
		/*! Simulation name */
		std::string SimName;
		/*! Poisson correction flag */
		std::string PoissonCorrection;
		int PoissonCorrectionCycle;
		/*! TrackParticleID */
		//bool *TrackParticleID;
		/*! SaveDirName */
		std::string SaveDirName;
		/*! RestartDirName */
		std::string RestartDirName;
		/*! restart_status 0 --> no restart; 1--> restart, create new; 2--> restart, append; */
		int restart_status;
		/*! last cycle */
		int last_cycle;

		/*! Boundary condition selection for BCFace for the electric field components */
		int bcEx[6], bcEy[6], bcEz[6];
		/*! Boundary condition selection for BCFace for the magnetic field components */
		int bcBx[6], bcBy[6], bcBz[6];

		/*! Boundary condition on particles 0 = exit 1 = perfect mirror 2 = riemission */
		/*! Boundary Condition Particles: FaceXright */
		int bcPfaceXright;
		/*! Boundary Condition Particles: FaceXleft */
		int bcPfaceXleft;
		/*! Boundary Condition Particles: FaceYright */
		int bcPfaceYright;
		/*! Boundary Condition Particles: FaceYleft */
		int bcPfaceYleft;
		/*! Boundary Condition Particles: FaceYright */
		int bcPfaceZright;
		/*! Boundary Condition Particles: FaceYleft */
		int bcPfaceZleft;


		/*! Field Boundary Condition 0 = Dirichlet Boundary Condition: specifies the valueto take pn the boundary of the domain 1 = Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain 2 = Periodic Condition */
		/*! Boundary Condition Electrostatic Potential: FaceXright */
		int bcPHIfaceXright;
		/*! Boundary Condition Electrostatic Potential:FaceXleft */
		int bcPHIfaceXleft;
		/*! Boundary Condition Electrostatic Potential:FaceYright */
		int bcPHIfaceYright;
		/*! Boundary Condition Electrostatic Potential:FaceYleft */
		int bcPHIfaceYleft;
		/*! Boundary Condition Electrostatic Potential:FaceZright */
		int bcPHIfaceZright;
		/*! Boundary Condition Electrostatic Potential:FaceZleft */
		int bcPHIfaceZleft;

		/*! Boundary Condition EM Field: FaceXright */
		int bcEMfaceXright;
		/*! Boundary Condition EM Field: FaceXleft */
		int bcEMfaceXleft;
		/*! Boundary Condition EM Field: FaceYright */
		int bcEMfaceYright;
		/*! Boundary Condition EM Field: FaceYleft */
		int bcEMfaceYleft;
		/*! Boundary Condition EM Field: FaceZright */
		int bcEMfaceZright;
		/*! Boundary Condition EM Field: FaceZleft */
		int bcEMfaceZleft;


		/*! GEM Challenge parameters */
		/*! current sheet thickness */
		double delta;
		/* Amplitude of the field */
		double B0x;
		double B0y;
		double B0z;
		double B1x;
		double B1y;
		double B1z;


		/*! boolean value for verbose results */
		//bool verbose;
		/*! RESTART */
		bool RESTART1;

		/*! velocity of the injection from the wall */
		double Vinj;

		/*! CG solver stopping criterium tolerance */
		double CGtol;
		/*! GMRES solver stopping criterium tolerance */
		double GMREStol;
		/*! mover predictor correcto iteration */
		int NiterMover;

		/*! Output for field */
		int FieldOutputCycle;
		std::string  FieldOutputTag;
		std::string  MomentsOutputTag;
		/*! Output for particles */
		int ParticlesOutputCycle;
		std::string ParticlesOutputTag;
		/*! Output for test particles */
		int TestParticlesOutputCycle;
		/*! test particles are flushed to disk every testPartFlushCycle  */
		int testPartFlushCycle;
		/*! restart cycle */
		int RestartOutputCycle;
		/*! Output for diagnostics */
		int DiagnosticsOutputCycle;
		/*! Call Finalize() at end of program execution (true by default) */
		bool CallFinalize;

		static Parameters read(std::string inputfile);

	private:

		void readInternal(std::string inputfile);

	};

} // end namespace ipic3d
