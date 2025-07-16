#include <boost/program_options.hpp> 
#include <boost/program_options/parsers.hpp>
#include <boost/filesystem/fstream.hpp>
#include <cstdlib>
#include <iostream> 
#include <string> 

#include "timemeasure.h"
#include "dryfric.hpp"

namespace 
{ 
	const size_t ERROR_IN_COMMAND_LINE = 1; 
	const size_t SUCCESS = 0; 
	const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
	
} 

int main(int argc, char** argv) 
{ 
	DRYOptions opt;
	try 
	{ 
		/** 
		 * Define and parse the program options 
		 */ 
		
		namespace po = boost::program_options; 
		namespace fs = boost::filesystem; 

		po::options_description desc("Options"); 
		desc.add_options() 
		("help", "Print help messages") 
		("totaltime,t", po::value<double>(&opt.tf), "Total time of simulation.")
		("particles,N", po::value<int>(&opt.N), "Number of realizations.")
		("prefix,o", po::value<std::string>(&opt.outputPrefix), "Path prefix for output files.")
		("DiffusionConstant,D", po::value<double>(&opt.D), "Diffusion constant.")
		("timestep", po::value<double>(&opt.dt), "Timestep for Langevin-Euler-Scheme.")
		("activity,f0", po::value<double>(&opt.f0), "Activity amplitude.")
		("delta", po::value<double>(&opt.delta), "Strength of Coulomb friction.")
		("tau", po::value<double>(&opt.tau), "Persistence time.")
		("tfTrajectory", po::value<double>(&opt.tfTrajectory), "Output frequency.")
		("eqtime", po::value<double>(&opt.teq), "Equilibration time before sampling starts.")
		;
		
		po::variables_map vm; 
		try 
		{ 
			po::store(po::parse_command_line(argc, argv, desc), vm); // can throw 
			
			fs::ifstream in("setup.cfg");
			po::store(po::parse_config_file(in, desc, true), vm);
			
			/** --help option 
			 */ 
			
			if ( vm.count("help")  ) 
			{ 
				std::cout << "Basic Command Line Parameter App" << std::endl 
				<< desc << std::endl; 
				return SUCCESS; 
			} 
			
			po::notify(vm); // throws on error, so do after help in case 
			// there are any problems 
		} 
		catch(po::error& e) 
		{ 
			std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
			std::cerr << desc << std::endl; 
			return ERROR_IN_COMMAND_LINE; 
		} 
		
		std::cout << "Simulation parameters:" << std::endl;
		std::cout << "Time: " << opt.tf << std::endl;
		std::cout << "Number of realizations: " << opt.N << std::endl;
		std::cout << "Path prefix: " << opt.outputPrefix << std::endl;
		std::cout << "dt: " << opt.dt << std::endl;
		std::cout << "D: " << opt.D << std::endl;
		std::cout << "f0: " << opt.f0 << std::endl;
		std::cout << "delta: " << opt.delta << std::endl;
		std::cout << "tau: " << opt.tau << std::endl;
		std::cout << "teq: " << opt.teq << std::endl;
		std::cout << "tf: " << opt.tf << std::endl;
		std::cout << "tfTrajectory: " << opt.tfTrajectory << std::endl;
		
		DRY sys(opt);
		
		sys.init();
		sys.equilibrate(opt.teq);
		sys.openOutputFiles();
		auto start = std::chrono::high_resolution_clock::now();
		sys.run(opt.tf);
		auto finish = std::chrono::high_resolution_clock::now();
		sys.closeOutputFiles();

		std::chrono::duration<double> elapsed = finish - start;

		std::cout << "Elapsed time: " << elapsed.count() << " s\n";
		
	} 
	catch(std::exception& e) 
	{ 
		std::cerr << "Unhandled Exception reached the top of main: " 
		<< e.what() << ", application will now exit" << std::endl; 
		return ERROR_UNHANDLED_EXCEPTION; 
		
	} 
	
	return EXIT_SUCCESS;
	
}
