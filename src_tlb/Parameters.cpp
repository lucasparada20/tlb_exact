#include "Parameters.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <string.h>
#include <vector>
#include <limits>
#include<algorithm>

int Parameters::Tmax = 0;
double Parameters::beta = 0;
int Parameters::time_interval = 0;
int Parameters::nb_scenarios = 0;
int Parameters::nb_real_scenarios = 0;
int Parameters::opt_cut_type = OPT_CUT_TYPE_BENDERS;
std::string Parameters::instance_file;
char Parameters::instance_format;
int Parameters::mcfp_solver = 0;
int Parameters::model = 1;
std::string Parameters::re_file_name;
std::string Parameters::output_file_name;
std::string Parameters::json_file_name = "";
int Parameters::delta = 0;
//Unbiased estimator
int Parameters::calculate_unbiased_estimator;
std::string Parameters::targets_file_name;
std::string Parameters::scenarios_file_name;
//DP
int Parameters::calculate_DP;
double Parameters::Qtot = 1.0;
//AvgScenario
int Parameters::AvgScenario = 0;
std::string Parameters::avg_scenario_file_name;
//OptStatCap Problem
double Parameters::Budget = 0; 

double Parameters::ParseNumber(const std::string& str) {
    size_t pos = str.find('/');
    if (pos != std::string::npos) {
        // It's a fraction
        double numerator = std::stod(str.substr(0, pos));
        double denominator = std::stod(str.substr(pos + 1));
        if (denominator == 0) {
			printf("Division by zero in fraction. Exiting"); exit(1);
        }
        return numerator / denominator;
    } else {
        // It's a standard floating-point number
        return std::stod(str);
    }
}

void Parameters::Read(int arg, char ** argv)
{
	printf("Reading parameters\n");
	for(int i=0;i<arg;i++)
	{
		char * first = strtok (argv[i]," ;=");
		char * second = strtok (NULL, " ;=");
		printf ("Parameter:%s value:%s\n",first,second);
		if(second == NULL) continue;
		
		if(strcmp(first, "instance_file") == 0)
		{
			instance_file = std::string(second);
		}
		else if(strcmp(first, "mcfp_solver") == 0)
		{
			sscanf(second, "%d", &mcfp_solver);
			printf("Solver:%s\n",GetSolverName());
		}
		else if(strcmp(first, "instance_format") == 0)
		{
			std::string instance_format_input = std::string(second);
			SetInstanceFormat(instance_format_input);
		}
		else if(strcmp(first, "model")==0)
		{
			model = std::stoi(second);
		}
		else if(strcmp(first, "re_file")==0)
		{
			re_file_name = std::string(second);
		}
		else if(strcmp(first, "output_file")==0)
		{
			output_file_name = std::string(second);
		}
		else if(strcmp(first,"unbiased_estimator")==0)
		{
			calculate_unbiased_estimator = std::stoi(second);
		}
		else if(strcmp(first,"scenarios_file_name")==0)
		{
			scenarios_file_name = std::string(second);
		}
		else if(strcmp(first,"targets_file_name")==0)
		{
			targets_file_name = std::string(second);
		}
		else if(strcmp(first,"calculate_DP")==0)
		{
			calculate_DP = std::stoi(second);
		}
		else if(strcmp(first,"Qtot")==0)
		{
			Qtot = ParseNumber(second);
		}
		else if(strcmp(first,"AvgScenario")==0)
		{
			AvgScenario = std::stoi(second);
		}
		else if(strcmp(first,"avg_scenario_file_name")==0)
		{
			avg_scenario_file_name = std::string(second);
		}
		else if(strcmp(first,"budget")==0)
		{
			Budget = std::stod(second);
		}
		else if(strcmp(first,"json")==0)
		{
			json_file_name = std::string(second);
		}
		else if(strcmp(first,"delta")==0)
		{
			delta = std::stoi(second);
		}

	}
}