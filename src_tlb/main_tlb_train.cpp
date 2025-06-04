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
#include "NodeTLB.h"
#include "LoadTLB.h"
#include "RecourseCalculation.h"
#include "DP.h"
#include "Solution.h"
#include "ExactTLB.h"

#include <omp.h> // omp_get_wtime();

int main(int arg, char ** argv)
{
	Parameters param;
	param.Read(arg,argv);
	
	Prob pr; Scenarios scs; LoadTLB LoadObject; 
	LoadObject.Load(&pr,&scs,Parameters::GetInstanceFileName());

	Sol sol(&pr);
	sol.SetScenarioCount(scs.GetScenarioCount());
	sol.SetODtripCount(scs.GetTripsCount());
	printf("Scenario_count:%d Trip_count:%d\n",scs.GetScenarioCount(),scs.GetTripsCount());
	
	sol.SetMatrixTargetLevels(scs.GetScenarioCount(),pr.GetNodeCount());
	sol.SetMatrixTargetLbLevels(scs.GetScenarioCount(),pr.GetNodeCount());
	sol.SetMatrixTargetUbLevels(scs.GetScenarioCount(),pr.GetNodeCount());

	printf("In main ...\n");	
	RecourseCalculation r(&scs,&pr);	
	pr.SetUpperBound(r.CalculateWS(&sol));
	printf("Setting upper bound as:%.1lf and unsatisfiedRequests:%.1lf\n",pr.GetUpperBound(),sol.GetWsUnsatisfiedReq());

	if(Parameters::GetModel()==2 || Parameters::GetModel()==6)
	{
		double start_time = omp_get_wtime();
		r.CalculateGlobalTargetLbUb();
		double end_time = omp_get_wtime();
		double global_bounds_time = end_time-start_time;
		
		r.CalculateTargetBounds();
		sol.SetIndividualBoundsTime(r.lb_i_ub_i_time);
		sol.SetTotalBoundsTime(global_bounds_time + r.lb_i_ub_i_time);
		printf("Time for TargetBounds:%.1lf TotalBounds:%.1lf\n",r.lb_i_ub_i_time,r.lb_i_ub_i_time+global_bounds_time);
	}
	if(Parameters::GetModel()==3 || Parameters::GetModel()==7)
	{	
		double start_time = omp_get_wtime();
		r.CalculateMinMaxBounds();
		double end_time = omp_get_wtime();
		sol.SetTotalBoundsTime(end_time-start_time);
		printf("MinMaxBounds time:%.1lf\n",end_time-start_time);
	}
	
	ExactTlb ex;
	ex.max_time = 1800;
	if(Parameters::GetModel()==11) ex.BandCheck = true;
	else ex.BandCheck = false;
	ex.Solve(&sol,&r);	
	
	//UnbiasedEstimator computations
	Scenarios other_scs;
	LoadObject.LoadScenarios(&pr,&other_scs,Parameters::GetScenariosFileName()); //N=400 scenarios

	RecourseCalculation other_r(&other_scs,&pr);
	
	sol.SetMatrixTargetLevels(other_scs.GetScenarioCount(),pr.GetNodeCount());
	sol.SetAvgTrips400(other_r.total_trip_count/(double)other_scs.GetScenarioCount());
	
	double WS_in400 = other_r.CalculateWS(&sol);
	sol.SetWS400( WS_in400 );
	printf("in400 AvgTripCount:%.1lf WS:%.1lf\n",other_r.total_trip_count/(double)other_scs.GetScenarioCount(),WS_in400); //exit(1);
	
	double start_time = omp_get_wtime();
	std::vector<int> WS_targets;
	//printf("WS_targets:\n");
	for(int i=0;i<pr.GetNodeCount();i++)
	{
		int target=0;
		for(int e=0;e<scs.GetScenarioCount();e++)
		{
			//printf("(%d,%d):%d ",e,i,sol.GetTarget(e,i));
			target += sol.GetTarget(e,i);
		}
			
		target /= scs.GetScenarioCount(); //One of many ways of averaging the scenario targets.
		target = target <= pr.GetNode(i)->stationcapacity ? target : pr.GetNode(i)->stationcapacity; 
		WS_targets.push_back( target );
		sol.SetWSTarget( target );
		//printf("t:%d ",target);
	}
	//printf("\n");
	
	double UnbiasedEstimator_WS = other_r.Calculate( WS_targets );
	double end_time = omp_get_wtime();
	printf("Scenarios for UnbiasedEstimator_WS:%d Value:%.2lf time_taken:%.1lf\n",other_scs.GetScenarioCount(),UnbiasedEstimator_WS,end_time-start_time);
	sol.SetUnbiasedWS( UnbiasedEstimator_WS );

	start_time = omp_get_wtime(); 
	std::vector<int> hp_vec; std::vector<int> hm_vec;
	if( Parameters::GetModel() == 5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
	{
		hp_vec = sol.GetHpVec(); hm_vec = sol.GetHmVec();
		for(int i=0;i<hp_vec.size();i++)
			printf("hp_vec[%d]:%d ",i,hp_vec[i]);
		printf("\n");
		for(int i=0;i<hm_vec.size();i++)
			printf("hm_vec[%d]:%d ",i,hm_vec[i]);
		printf("\n");			
	}
	double UnbiasedEstimator_Exact = other_r.Calculate(sol.GetExactTargets(), hm_vec, hp_vec);
	end_time = omp_get_wtime();
	printf("Scenarios for UnbiasedEstimator_Exact:%d Value:%.2lf time_taken:%.1lf\n",other_scs.GetScenarioCount(),UnbiasedEstimator_Exact,end_time-start_time);
	sol.SetUnbiasedExact( UnbiasedEstimator_Exact );
	
	clock_t start_time_dp = clock();
	DP dp;
	std::vector<int> capacities;
	for(int i=0; i<pr.GetNodeCount();i++)
		capacities.push_back( pr.GetNode(i)->stationcapacity );
		
	std::vector<int> y;
	double cost = dp.GetCostDP(capacities,pr.GetQtot(),y);
	clock_t end_time_dp = clock();
	double time_taken = (double)(end_time_dp - start_time_dp) / CLOCKS_PER_SEC;
	printf("DP cost:%.2lf TimeTaken:%.1lf\n",cost,time_taken);
	
	start_time = omp_get_wtime(); 
	double UnbiasedEstimator_DP = other_r.Calculate(y);
	end_time = omp_get_wtime();
	printf("Scenarios for UnbiasedEstimator_DP:%.2lf time_taken:%.1lf\n",
			UnbiasedEstimator_DP,end_time-start_time);
	sol.SetUnbiasedDP( UnbiasedEstimator_DP );			

	//Trips Lb	
	std::vector<double> lbs = other_r.CalculateLb(sol.GetExactTargets());
	sol.SetUnbiasedLb(lbs[0]);
	sol.SetUnbiasedRejectedPicks(lbs[1]);
	sol.SetUnbiasedRejectedDrops(lbs[2]);
	
	printf("Unbiased Lb:%.1lf Ub:%.1lf Rejected Picks:%.1lf Drops:%.1lf\n",
		sol.GetUnbiasedLb(),
		sol.GetUnbiasedExact(),
		sol.GetUnbiasedRejectedPicks(),
		sol.GetUnbiasedRejectedDrops());
	
	//the names need to be something like: n%d_e%d.txt
	sol.WriteReFile(Parameters::GetReFileName());
	sol.WriteExactTargetsToFile(Parameters::GetOutputFileName());
	
	return 0;
}

