#include <boost/program_options.hpp>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "dryfric.hpp"

namespace
{
    constexpr int ERROR_IN_COMMAND_LINE = 1;
    constexpr int SUCCESS = 0;
    constexpr int ERROR_UNHANDLED_EXCEPTION = 2;
}


int main(int argc, char** argv)
{
    DRYOptions opt;

    try
    {
        namespace po = boost::program_options;

        /*
         * Define command-line/config-file options.
         */
        po::options_description desc("Simulation options");

        desc.add_options()
            ("help,h",
             "Print help message")

            ("totaltime,t",
             po::value<double>(&opt.tf),
             "Total simulation time")

            ("particles,N",
             po::value<int>(&opt.N),
             "Number of particles")

            ("prefix,o",
             po::value<std::string>(&opt.outputPrefix),
             "Path prefix for output files")

            ("DiffusionConstant,D",
             po::value<double>(&opt.D),
             "Diffusion constant")

            ("timestep",
             po::value<double>(&opt.dt),
             "Time step for Langevin-Euler scheme")

            ("activity,v0",
             po::value<double>(&opt.v0),
             "Activity amplitude")

            ("delta",
             po::value<double>(&opt.delta),
             "Dry-friction strength")

            ("tau",
             po::value<double>(&opt.tau),
             "Activity time scale")

            ("tfTrajectory",
             po::value<double>(&opt.tfTrajectory),
             "Trajectory output interval")

            ("eqtime",
             po::value<double>(&opt.teq),
             "Equilibration time before sampling");


        /*
         * Parse command-line arguments and setup.cfg.
         */
        po::variables_map vm;

        try
        {
            po::store(
                po::parse_command_line(argc, argv, desc),
                vm
            );

            std::ifstream configFile("setup.cfg");

            if (configFile)
            {
                po::store(
                    po::parse_config_file(
                        configFile,
                        desc,
                        true
                    ),
                    vm
                );
            }

            /*
             * --help
             */
            if (vm.count("help"))
            {
                std::cout
                    << "Dry-friction particle simulation\n\n"
                    << desc
                    << '\n';

                return SUCCESS;
            }

            /*
             * Apply parsed values to DRYOptions.
             */
            po::notify(vm);
        }
        catch (const po::error& e)
        {
            std::cerr
                << "ERROR: "
                << e.what()
                << "\n\n"
                << desc
                << '\n';

            return ERROR_IN_COMMAND_LINE;
        }


        /*
         * Print simulation parameters.
         */
        std::cout
            << "\nSimulation parameters:\n"
            << "----------------------\n"
            << "Total time:          " << opt.tf << '\n'
            << "Particle number:     " << opt.N << '\n'
            << "Output prefix:       " << opt.outputPrefix << '\n'
            << "Time step:           " << opt.dt << '\n'
            << "Diffusion constant:  " << opt.D << '\n'
            << "Activity v0:         " << opt.v0 << '\n'
            << "Dry friction delta:  " << opt.delta << '\n'
            << "Activity tau:        " << opt.tau << '\n'
            << "Equilibration time:  " << opt.teq << '\n'
            << "Trajectory interval: " << opt.tfTrajectory << '\n'
            << '\n';


        /*
         * Create and initialize the simulation.
         */
        DRY sys(opt);

        sys.init();


        /*
         * Equilibration.
         */
        std::cout << "Equilibrating...\n";

        sys.equilibrate(opt.teq);

        std::cout << "Equilibration complete.\n";


        /*
         * Open output files only after equilibration.
         */
        sys.openOutputFiles();


        /*
         * Run the actual simulation and measure runtime.
         */
        const auto start =
            std::chrono::high_resolution_clock::now();

        sys.run(opt.tf);

        const auto finish =
            std::chrono::high_resolution_clock::now();


        /*
         * Close output files.
         */
        sys.closeOutputFiles();


        /*
         * Report runtime.
         */
        const std::chrono::duration<double> elapsed =
            finish - start;

        std::cout
            << "Simulation complete.\n"
            << "Elapsed time: "
            << elapsed.count()
            << " s\n";
    }
    catch (const std::exception& e)
    {
        std::cerr
            << "Unhandled exception reached the top of main:\n"
            << e.what()
            << "\n";

        return ERROR_UNHANDLED_EXCEPTION;
    }

    return SUCCESS;
}