#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <csignal>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string.h> //strcm

#include "Parameters.h"
#include "ProblemDefinition.h"
#include "RecourseCalculation.h"
#include "NodeTLB.h"
#include "LoadTLB.h"
#include "Solution.h"
#include <omp.h> // omp_get_wtime();

int main(int arg, char ** argv)
{
	printf("Main tlb test ...\n");
	Parameters param;
	param.Read(arg,argv);
	
	Prob pr; Scenarios scs; Scenarios scs_400; LoadTLB LoadObject; 
	LoadObject.Load(&pr,&scs,Parameters::GetInstanceFileName());
	
	std::vector<int> targets;
	printf("Targets filename:%s\n",
			Parameters::GetTargetsFileName());
	LoadObject.LoadTargets(targets,Parameters::GetTargetsFileName());

	LoadObject.LoadScenarios(&pr,&scs_400,Parameters::GetScenariosFileName()); 
	
	Sol sol(&pr);
	 
	sol.SetScenarioCount(scs_400.GetScenarioCount());
	sol.SetODtripCount(scs_400.GetTripsCount());
	printf("Scenario_count:%d Trip_count:%d\n",scs_400.GetScenarioCount(),scs_400.GetTripsCount());

	RecourseCalculation r(&scs_400,&pr);	
	
	//sol.SetMatrixTargetLevels(scs_400.GetScenarioCount(),pr.GetNodeCount()); //Required to store WS targets
	//double WS_in400 = r.CalculateWS(&sol);
	//printf("WS:%.1lf\n",WS_in400);
	
	double start_time = omp_get_wtime();
	//r.show = true;
	
	pr.SetTgtLvlVectors(scs_400.GetScenarioCount()); // Reserve memory to store the node target levels from Schuijbroek (2017)
	r.calculate_target_levels = true;
	double UnbiasedEstimator = r.Calculate( targets );
	double end_time = omp_get_wtime();
	printf("Unbiased Estimator value:%.2lf time_taken:%.1lf\n",UnbiasedEstimator,end_time-start_time);
	
	std::vector<std::vector<double>> pick_tgt_lvls = pr.GetPickTgtLvls();
	std::vector<std::vector<double>> del_tgt_lvls = pr.GetDelTgtLvls();
	
	// Pick target level bins
	int total_non_neg_pick = 0;
	int count_leq_085_pick = 0;
	int count_085_090_pick = 0;
	int count_090_095_pick = 0;
	int count_095_097_pick = 0;
	int count_097_099_pick = 0;
	int count_099_100_pick = 0;

	for (std::vector<double> & row : pick_tgt_lvls) {
		for (double val : row) {
			if (val >= 0.0) {
				total_non_neg_pick++;
				if (val <= 0.85)
					count_leq_085_pick++;
				else if (0.85 <= val && val <= 0.90)
					count_085_090_pick++;
				else if (0.90 <= val && val <= 0.95)
					count_090_095_pick++;
				else if (0.95 <= val && val <= 0.97)
					count_095_097_pick++;
				else if (0.97 <= val && val <= 0.99)
					count_097_099_pick++;
				else
					count_099_100_pick++;
			}
		}
	}

	// Delivery target level bins
	int total_non_neg_del = 0;
	int count_leq_085_del = 0;
	int count_085_090_del = 0;
	int count_090_095_del = 0;
	int count_095_097_del = 0;
	int count_097_099_del = 0;
	int count_099_100_del = 0;

	for (std::vector<double> & row : del_tgt_lvls) {
		for (double val : row) {
			if (val >= 0.0) {
				total_non_neg_del++;
				if (val <= 0.85)
					count_leq_085_del++;
				else if (0.85 <= val && val <= 0.90)
					count_085_090_del++;
				else if (0.90 <= val && val <= 0.95)
					count_090_095_del++;
				else if (0.95 <= val && val <= 0.97)
					count_095_097_del++;
				else if (0.97 <= val && val <= 0.99)
					count_097_099_del++;
				else
					count_099_100_del++;
			}
		}
	}

	// Print results
	printf("\n=== Pick target level breakdown (%% of total non-negatives) ===\n");
	printf("<= 0.85             : %.1lf%%\n", 100.0 * count_leq_085_pick / total_non_neg_pick);
	printf("[0.85 – 0.90]       : %.1lf%%\n", 100.0 * count_085_090_pick / total_non_neg_pick);
	printf("[0.90 – 0.95]       : %.1lf%%\n", 100.0 * count_090_095_pick / total_non_neg_pick);
	printf("[0.95 – 0.97]       : %.1lf%%\n", 100.0 * count_095_097_pick / total_non_neg_pick);
	printf("[0.97 – 0.99]       : %.1lf%%\n", 100.0 * count_097_099_pick / total_non_neg_pick);
	printf("[0.99 – 1.00]       : %.1lf%%\n", 100.0 * count_099_100_pick / total_non_neg_pick);

	printf("\n=== Delivery target level breakdown (%% of total non-negatives) ===\n");
	printf("<= 0.85             : %.1lf%%\n", 100.0 * count_leq_085_del / total_non_neg_del);
	printf("[0.85 – 0.90]       : %.1lf%%\n", 100.0 * count_085_090_del / total_non_neg_del);
	printf("[0.90 – 0.95]       : %.1lf%%\n", 100.0 * count_090_095_del / total_non_neg_del);
	printf("[0.95 – 0.97]       : %.1lf%%\n", 100.0 * count_095_097_del / total_non_neg_del);
	printf("[0.97 – 0.99]       : %.1lf%%\n", 100.0 * count_097_099_del / total_non_neg_del);
	printf("[0.99 – 1.00]       : %.1lf%%\n", 100.0 * count_099_100_del / total_non_neg_del);



	return 0;
}

