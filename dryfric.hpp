#ifndef DRY_H
#define DRY_H

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "MersenneTwister.h"

/*
 * Compile-time switches
 */

// Log particle trajectories
//#define OUTPUT_TRAJECTORIES


/*
 * Simulation parameters
 */
struct DRYOptions {
    // Number of particles
    int N = 20;

    // Diffusion constant
    double D = 1.0;

    // Total simulation time (excluding equilibration)
    double tf = 2000.0;

    // Equilibration time
    double teq = 1000.0;

    // Integration time step
    double dt = 0.001;

    // Time scale of the activity
    double tau = 1.0;

    // Activity strength
    double v0 = 0.01;

    // Dry friction coefficient
    double delta = 0.0;

    // Time interval between trajectory outputs
    double tfTrajectory = 1.0;

    // Output filename prefix
    std::string outputPrefix = "../data/test";
};


class DRY {

public:

    explicit DRY(const DRYOptions& opt)
        : N(opt.N),
        D(opt.D),
        t(0.0),
        tf(opt.tf),
        teq(opt.teq),
        dt(opt.dt),
        tau(opt.tau),
        v0(opt.v0),
        delta(opt.delta),
        tfTrajectory(opt.tfTrajectory),
        tTrajectory(0.0),
        varNoise(0.0),
        activityNoise(0.0),
        x(N, 0.0),
        y(N, 0.0),
        xinit(N, 0.0),
        yinit(N, 0.0),
        vx(N, 0.0),
        vy(N, 0.0),
        nx(N, 0.0),
        ny(N, 0.0),
        rho(2 * bin * length_bin, 0.0),
        randy(),
        filePrefix(opt.outputPrefix),
        trajectoryFile(nullptr),
        miscMeasurementFile(nullptr),
        densityFile(nullptr)
    {
    }


    ~DRY()
    {
        closeOutputFiles();
    }


    /*
     * Initialize the simulation.
     */
    void init()
    {
        std::fill(x.begin(), x.end(), 0.0);
        std::fill(y.begin(), y.end(), 0.0);

        std::fill(vx.begin(), vx.end(), 0.0);
        std::fill(vy.begin(), vy.end(), 0.0);

        std::fill(nx.begin(), nx.end(), 0.0);
        std::fill(ny.begin(), ny.end(), 0.0);

        std::fill(rho.begin(), rho.end(), 0.0);

        // Store the initial positions for MSD calculation.
        xinit = x;
        yinit = y;

        // Noise amplitudes.
        varNoise = std::sqrt(2.0 * D * dt);
        activityNoise = std::sqrt(2.0 * dt / tau);

        t = 0.0;
        tTrajectory = 0.0;
    }


    /*
     * Advance the simulation by one time step.
     */
    void step()
    {
        for (int i = 0; i < N; ++i) {

            /*
             * Position update.
             */
            x[i] += dt * vx[i];
            y[i] += dt * vy[i];


            /*
             * Current speed.
             *
             * Avoid division by zero in the dry-friction term.
             */
            const double v2 = vx[i] * vx[i] + vy[i] * vy[i];
            const double v = std::sqrt(std::max(v2, 1e-16));


            /*
             * Velocity update.
             *
             * Dry friction:
             *
             *     -delta * v / |v|
             *
             * Active force:
             *
             *     v0 * n
             */
            vx[i] += dt * (
                -delta * vx[i] / v
                + v0 * nx[i]
            );

            vy[i] += dt * (
                -delta * vy[i] / v
                + v0 * ny[i]
            );


            /*
             * Thermal/Langevin noise.
             */
            vx[i] += varNoise * randy.randNorm();
            vy[i] += varNoise * randy.randNorm();


            /*
             * Ornstein-Uhlenbeck activity dynamics.
             */
            nx[i] += dt * (-nx[i] / tau)
                   + activityNoise * randy.randNorm();

            ny[i] += dt * (-ny[i] / tau)
                   + activityNoise * randy.randNorm();


            /*
             * Velocity histogram.
             */
            if (v < length_bin) {
                const int index =
                    static_cast<int>(std::round(v * bin));

                if (index >= 0 &&
                    index < static_cast<int>(rho.size())) {
                    ++rho[index];
                }
            }
        }

        t += dt;
    }


    /*
     * Run equilibration without measurements.
     */
    void equilibrate(double equilibrationTime)
    {
        t = 0.0;

        while (t < equilibrationTime) {
            step();
        }

        /*
         * The equilibrated state becomes the reference
         * state for the subsequent MSD measurement.
         */
        xinit = x;
        yinit = y;

        /*
         * Reset simulation time after equilibration.
         */
        t = 0.0;
    }


    /*
     * Open output files.
     */
    void openOutputFiles()
    {
        closeOutputFiles();

        /*
         * Miscellaneous measurements.
         */
        {
            const std::string filename =
                filePrefix + "-misc_measurements.csv";

            miscMeasurementFile =
                std::fopen(filename.c_str(), "w");

            if (miscMeasurementFile == nullptr) {
                std::cerr
                    << "ERROR: Unable to open "
                    << filename
                    << '\n';

                std::exit(EXIT_FAILURE);
            }
        }


        /*
         * Trajectories.
         */
#ifdef OUTPUT_TRAJECTORIES

        {
            const std::string filename =
                filePrefix + "-trajectories.csv";

            trajectoryFile =
                std::fopen(filename.c_str(), "w");

            if (trajectoryFile == nullptr) {
                std::cerr
                    << "ERROR: Unable to open "
                    << filename
                    << '\n';

                std::exit(EXIT_FAILURE);
            }
        }


        /*
         * Velocity-density histogram.
         */
        {
            const std::string filename =
                filePrefix + "-density.csv";

            densityFile =
                std::fopen(filename.c_str(), "w");

            if (densityFile == nullptr) {
                std::cerr
                    << "ERROR: Unable to open "
                    << filename
                    << '\n';

                std::exit(EXIT_FAILURE);
            }
        }

#endif
    }


    /*
     * Close all output files.
     */
    void closeOutputFiles()
    {
        if (miscMeasurementFile != nullptr) {
            std::fclose(miscMeasurementFile);
            miscMeasurementFile = nullptr;
        }

#ifdef OUTPUT_TRAJECTORIES

        if (trajectoryFile != nullptr) {
            std::fclose(trajectoryFile);
            trajectoryFile = nullptr;
        }

        if (densityFile != nullptr) {
            std::fclose(densityFile);
            densityFile = nullptr;
        }

#endif
    }


    /*
     * Run the simulation.
     */
    void run(double finalTime)
    {
        tf = finalTime;

        t = 0.0;
        tTrajectory = 0.0;

        /*
         * First MSD measurement after tfTrajectory.
         * Measurements then become progressively less frequent.
         */
        double tMSD = tfTrajectory;
        double timeSinceMSD = 0.0;

        while (t < tf) {

            step();

            /*
             * Trajectory output timer.
             */
            tTrajectory += dt;
            timeSinceMSD += dt;


#ifdef OUTPUT_TRAJECTORIES

            if (tTrajectory >= tfTrajectory) {

                tTrajectory -= tfTrajectory;

                std::fprintf(
                    trajectoryFile,
                    "%f, ",
                    t
                );

                for (int i = 0; i < N; ++i) {

                    std::fprintf(
                        trajectoryFile,
                        "%f, %f",
                        x[i],
                        y[i]
                    );

                    if (i < N - 1) {
                        std::fprintf(
                            trajectoryFile,
                            ", "
                        );
                    }
                }

                std::fprintf(
                    trajectoryFile,
                    "\n"
                );
            }

#endif


            /*
             * MSD measurement.
             */
            if (timeSinceMSD >= tMSD) {

                timeSinceMSD = 0.0;
                tMSD *= 1.1;

                double meanSquaredDisplacement = 0.0;

                for (int i = 0; i < N; ++i) {

                    const double dx =
                        x[i] - xinit[i];

                    const double dy =
                        y[i] - yinit[i];

                    meanSquaredDisplacement +=
                        dx * dx + dy * dy;
                }

                meanSquaredDisplacement /= N;

                std::fprintf(
                    miscMeasurementFile,
                    "%f, %e\n",
                    t,
                    meanSquaredDisplacement
                );
            }
        }


        /*
         * Write velocity-density histogram.
         */
#ifdef OUTPUT_TRAJECTORIES

        double sum = 0.0;

        for (double value : rho) {
            sum += value;
        }

        if (sum > 0.0) {

            for (std::size_t i = 0; i < rho.size(); ++i) {

                const double velocity =
                    static_cast<double>(i) / bin;

                const double probabilityDensity =
                    rho[i] * bin / sum;

                std::fprintf(
                    densityFile,
                    "%f, %e\n",
                    velocity,
                    probabilityDensity
                );
            }
        }

#endif
    }

private:

    /*
     * Simulation parameters.
     */
    int N;

    double D;
    double t = 0.0;
    double tf;
    double teq;
    double dt;
    double tau;
    double v0;
    double delta;

    /*
     * Trajectory output.
     */
    double tfTrajectory;
    double tTrajectory = 0.0;

    /*
     * Noise amplitudes.
     */
    double varNoise = 0.0;
    double activityNoise = 0.0;

    /*
     * Particle state.
     */
    std::vector<double> x;
    std::vector<double> y;

    std::vector<double> xinit;
    std::vector<double> yinit;

    std::vector<double> vx;
    std::vector<double> vy;

    std::vector<double> nx;
    std::vector<double> ny;

    /*
     * Velocity histogram.
     */
    static constexpr int bin = 10;
    static constexpr int length_bin = 20;

    std::vector<double> rho;

    /*
     * Random-number generator.
     */
    MTRand randy;

    /*
     * Output.
     */
    std::string filePrefix;

    FILE* trajectoryFile = nullptr;
    FILE* miscMeasurementFile = nullptr;
    FILE* densityFile = nullptr;
};

#endif // DRY_H