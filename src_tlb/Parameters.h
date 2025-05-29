#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <time.h>
#include <math.h>
#include <string>

#define OPT_CUT_TYPE_BENDERS 1

class Parameters
{

public:
	
	static void SetBeta(double b){ beta = b; }
	static void SetTimeInterval(int t){ time_interval = t; }
	static void SetTypeOfOptimalityCuts(int t){opt_cut_type = t;}
	static void SetTmax(int t){ Tmax=t; }
	static void SetScenarioCount(int e){ nb_scenarios=e;}
	static void SetRealScenarioCount(int e){ nb_real_scenarios=e;}
	static void SetSolver(int i) {mcfp_solver=i;}
	static void SetInstanceFormat(std::string str){ instance_format = (str == "trips") ? 'T' : 'R';}
	
	static int GetSolver() {return mcfp_solver;}
	static char * GetSolverName() {
		switch (mcfp_solver) {
			case 1:
				return (char*)"MCFP_SOLVER_UNIPI_SIMPLEX";
			case 2:
				return (char*)"MCFP_SOLVER_UNIPI_RELAXIV";
			case 3:
				return (char*)"MCFP_SOLVER_LEMON_NETSIMPLEX";
			case 4:
				return (char*)"MCFP_SOLVER_LEMON_CSCALING";
			default:
				return (char*)"Unknown Solver";
		}
	}
	static int GetModel(){ return model; }
	static int GetTmax(){ return Tmax; }
	static double GetBeta(){ return beta; }
	static int GetTimeInterval(){ return time_interval; }
	static int GetTypeOfOptimalityCuts(){return opt_cut_type;}
	static char* GetInstanceFileName(){return instance_file.size() == 0?NULL:(char *)instance_file.c_str();}
	static char GetInstanceFormat(){return instance_format;}
	static int GetScenarioCount(){ return nb_scenarios; }
	static int GetRealScenarioCount(){ return nb_real_scenarios; }
	static char* GetReFileName(){ return (char *)re_file_name.c_str(); }
	static char* GetOutputFileName(){ return (char *)output_file_name.c_str(); }
	static std::string GetJsonFileName(){ return json_file_name; }
	static int GetDelta(){ return delta; }
	//UnbiasedEstimator
	static int UnbiasedEstimator(){return calculate_unbiased_estimator;}
	static char* GetScenariosFileName(){ return (char*)scenarios_file_name.c_str(); } 
	static char* GetTargetsFileName(){ return (char*)targets_file_name.c_str(); } 
	//Capacities for the DP
	static int CalculateDP(){ return calculate_DP; }
	//AvgScenario
	static char* GetAvgScenarioFileName(){ return (char*)avg_scenario_file_name.c_str(); }
	static int CalculateAvgScenario(){ return AvgScenario; }
	//OptStatCap Problem
	static double GetBudget(){ return Budget; }

	static double GetQtot(){ return Qtot; }
	
	static double ParseNumber(const std::string& str);
	
	//Function to read all input parameters
	void Read(int arg, char ** argv);
	
private:

	static double beta;
	static int time_interval;
	static int opt_cut_type;
	static std::string instance_file;
	static char instance_format;
	static int Tmax;
	static int nb_scenarios;
	static int nb_real_scenarios;
	static int mcfp_solver;
	static int model;
	static std::string re_file_name;
	static std::string output_file_name;
	static std::string json_file_name;
	static int delta;
	
	//UnbiasedEstimator
	static int calculate_unbiased_estimator;
	static std::string scenarios_file_name;
	static std::string targets_file_name;
	//DP
	static int calculate_DP;
	static double Qtot;
	//AvgScenario
	static int AvgScenario;
	static std::string avg_scenario_file_name;
	//OptStatCap Problem
	static double Budget;

};


#endif
