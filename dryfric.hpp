#ifndef __DRY_H__
#define __DRY_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

#include "MersenneTwister.h"
#include "Helper.hpp"
#include "Vector.hpp"

// If defined, the trajectories are logged into a file
#define OUTPUT_TRAJECTORIES

#endif

/* Data structure for transfering simulation parameters from the main program */
class DRYOptions {
public:
	/* Number of realizations */
	int N;
	/* Diffusion constant */
	double D;
	/* Total time of the simulation run (excluding equilibration time) */
	double tf;
	/* Length of the equilibration interval */
	double teq;
	/* Time step of the simulation algorithm */
	double dt;
	/* Persistence time */
	double tau;
	/* Activity */
	double f0;
	/* Dry friction */
	double delta;
	/* Time interval for output */
	double tfTrajectory;

	/* This prefix is prepended before the filenames of all output file */
	std::string outputPrefix;
	
	DRYOptions() {
		/* Set some basic default values */
		N = 20;
		D = 1.0;
		dt = 0.001;
		tf = 2000;
		teq = 1000.0;
		tau = 1.0;
		f0 = 0.01;
		tfTrajectory = 1.0;
		outputPrefix = "../data/test";
	}
};

class DRY {
	
public:

	/* Number of realizations */
	int N; 
	/* Diffusionconstant */
	double D;
	/* Current simulation time */
	double t;
	/* Total time of the simulation run (excluding equilibration time) */
	double tf;
	/* Length of the equilibration interval */
	double teq;
	/* Time step of the simulation algorithm */
	double dt;
	/* Time-scaling */
	double tau;
	/* Activity */
	double f0;
	/* Dry friction */
	double delta;
	/* Time interval for output */
	double tfTrajectory;
	/* Variance of the white noise in Langevin equation */
	double varNoise;
	/* Variance of the white noise in activity equation */
	double activityNoise;
	/* Current positions of the particles */
	Vector<double> *x;
	/* Initial positions of the particles */
	Vector<double> *xinit;
	/* Current velocites of the particles */
	Vector<double> *v;
	/* Current ativities of the particles */
	Vector<double> *n;
	/* Probability distribution */
	Vector<double> *rho;
	int bin; // Number of bins in the interval [0, 1[
	int length_bin; // Interval [-l, l] for the probability distribution p(v)
	/* MersenneTwister Random Number Generator */
	MTRand randy;
	/* Output files*/
	FILE *miscMeasurementFile;
	#ifdef OUTPUT_TRAJECTORIES
	FILE *trajectoryFile;
	FILE *velocityFile;
	FILE *activityFile;
	#endif
	FILE *densityFile;
	double tTrajectory;
	std::string filePrefix;

	/* Constructor that initializes memory and sets parameters given by DRYOptions object */
	DRY(const DRYOptions &opt) {
		N = opt.N;
		x = new Vector<double>(N);
		v = new Vector<double>(N);
		n = new Vector<double>(N);
		xinit = new Vector<double>(N);
		bin = 10;
		length_bin = 101;

		rho = new Vector<double>(2* bin * length_bin);
	  
		dt = opt.dt;
		D = opt.D;
		tau = opt.tau;
		tf = opt.tf;
		delta = opt.delta;
		f0 = opt.f0;
		teq = opt.teq;
		
		filePrefix = opt.outputPrefix;
		tfTrajectory = opt.tfTrajectory;

	}
	
	~DRY() {
	}
	
	void init() {

		for (int i = 0; i < N; ++i) {
			(*x)[i] = 0.0;
			(*v)[i] = 0.0;
			(*n)[i] = 0.0;
		}

		for (int i = 0; i < 2 * bin * length_bin; ++i)
			(*rho)[i] = 0.0;
		
		/* Setup parameters for later use */
		varNoise = sqrt(2 * D * dt);
		activityNoise = sqrt(2 * dt / tau);
		t = 0.0;

	}

	double sign(double a) {
		if (a > 0) return 1.0;
		else if (a < 0) return -1.0;
		else return 0.0;
	}
	
	void step() {
		int i;

		for (i = 0; i < N; ++i) {
			(*x)[i] += dt * (*v)[i];
			(*v)[i] += dt * (-delta * sign((*v)[i])  + (*n)[i] * f0) + varNoise * randy.randNorm();
			(*n)[i] += dt * (-(*n)[i] / tau) + activityNoise * randy.randNorm();
			if (fabs((*v)[i]) < length_bin) (*rho)[round(((*v)[i] + length_bin)*bin)] += 1;
		}

		t += dt;

	}

	/* Do an equilibration run without measurements for the time given by equilibrationTime */
	void equilibrate(double equilibrationTime) {
		t = 0.0;
		while (t < equilibrationTime) {
			step();
		}
		(*x) = (*xinit);
	}
	
	/* Open all necessary files for logging the measured results */
	void openOutputFiles()  {
		std::string filename = filePrefix + "-density.csv";
		densityFile = fopen(filename.c_str(), "w");
		if (densityFile == NULL) {
			std::cerr << "ERROR: Unable to open activityFile file " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
		filename = filePrefix + "-misc_measurements.csv";
		miscMeasurementFile = fopen(filename.c_str(), "w");
		if (miscMeasurementFile == NULL) {
			std::cerr << "ERROR: Unable to open miscMeasurementFile file " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
		#ifdef OUTPUT_TRAJECTORIES
		filename = filePrefix + "-trajectories.csv";
		trajectoryFile = fopen(filename.c_str(), "w");
		if (trajectoryFile == NULL) {
			std::cerr << "ERROR: Unable to open trajectory file " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
		filename = filePrefix + "-velocity.csv";
		velocityFile = fopen(filename.c_str(), "w");
		if (velocityFile == NULL) {
			std::cerr << "ERROR: Unable to open velocityFile file " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
		filename = filePrefix + "-activity.csv";
		activityFile = fopen(filename.c_str(), "w");
		if (activityFile == NULL) {
			std::cerr << "ERROR: Unable to open activityFile file " << filename << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}
	
	/* Close all files */
	void closeOutputFiles() {
		fclose(miscMeasurementFile);
		fclose(densityFile);
		#ifdef OUTPUT_TRAJECTORIES
		fclose(trajectoryFile);
		fclose(velocityFile);
		fclose(activityFile);
		#endif
	}

	/* Run simulation for length of finalTime */
	void run(double finalTime) {
		tf = finalTime;
		t = 0.0;
		tTrajectory = 0.0;
		double tMSD = tfTrajectory;
		double time = 0.0;
		double MeanSquaredDisplacement;
		int i;
		/* Main simulation loop */
		while (t < tf) {
			step();
	  
			tTrajectory += dt;
			time += dt;
			#ifdef OUTPUT_TRAJECTORIES
			if (tTrajectory >= tfTrajectory) {
				tTrajectory = 0.0;
				fprintf(trajectoryFile, "%f, ", t);
				fprintf(velocityFile, "%f, ", t);
				fprintf(activityFile, "%f, ", t);
				for (i = 0; i < N-1; ++i) {
					fprintf(trajectoryFile, "%f, ", (*x)[i]);
					fprintf(velocityFile, "%f, ", (*v)[i]);
					fprintf(activityFile, "%f, ", (*n)[i]);
				}
				fprintf(trajectoryFile, "%f\n", (*x)[N-1]);
				fprintf(velocityFile, "%f\n", (*v)[N-1]);
				fprintf(activityFile, "%f\n", (*n)[N-1]);
			}
			#endif


			if (time >= tMSD) {
				time = 0.0;
				tMSD *= 1.1;
				fprintf(miscMeasurementFile, "%e, ", t);
				MeanSquaredDisplacement = 0.0;
				for (i = 0; i < N; ++i) MeanSquaredDisplacement += ((*x)[i] - (*xinit)[i]) * ((*x)[i] - (*xinit)[i]);
				fprintf(miscMeasurementFile, "%e\n", MeanSquaredDisplacement / N);
			}

		}

		double sum = 0.0;
		for (i = 0; i < 2 * bin * length_bin; i++)
			sum += (*rho)[i];
		for (i = 0; i < 2 * bin * length_bin; i++) {
			fprintf(densityFile, "%f, ", (i - bin*length_bin) / double(bin));
			fprintf(densityFile, "%e\n", (*rho)[i] * bin / sum);
		}
	}

};

