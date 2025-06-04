
#ifndef _SOLUTION_H
#define _SOLUTION_H

#include <stdio.h>
#include <stdlib.h>
#include "ProblemDefinition.h"
#include "Parameters.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip> //std::setprecision()
#include <numeric>      // std::accumulate
#include <algorithm> //std::sort

class Sol
{
	public:
		Sol(Prob * prob): _prob(prob), scenario_count(-1), ODtrip_count(-1), lb_i_ub_i_time(-1.0), total_bounds_time(-1.0), nb_benders_cuts(-1),re(-1),cplex_nodes(-1),exact_ub(-1.0),exact_lb(-1.0),exact_time(-1.0),UnbiasedEstimator_DP(-1.0), UnbiasedEstimator_Exact(-1.0), UnbiasedEstimator_WS(-1.0), Avg_WS(-1.0), Avg_Exact(-1.0), Avg_DP(-1.0), WS_in400(-1.0) {
			for(int i=0;i<_prob->GetNodeCount();i++)
				nodes.push_back(*_prob->GetNode(i));
		}
		Sol(const Sol * other) {
			matrix_target_lb_levels = other->matrix_target_lb_levels;
			matrix_target_ub_levels = other->matrix_target_ub_levels;
			nodes = other->nodes;
		}
		Prob * GetProblem(){ return _prob;}
		Node * GetNode(int i){ return &nodes[i]; }
		int GetNodeCount(){ return (int)nodes.size(); }
		
		void SetODtripCount(int v){ODtrip_count = v;}
		void SetScenarioCount(int v){scenario_count=v;}
		
		void SetMatrixTargetLevels(int scenarios, int nb_stations) 
		{ 
			matrix_target_levels.resize(scenarios,std::vector<int>());
			for(int e=0;e<scenarios;e++)
				matrix_target_levels[e].resize(nb_stations,-1);
			scenario_objective.resize(scenarios,0.0);
		}
		void SetMatrixTargetLbLevels(int scenarios, int nb_stations) 
		{ 
			matrix_target_lb_levels.resize(scenarios,std::vector<int>());
			for(int e=0;e<scenarios;e++)
				matrix_target_lb_levels[e].resize(nb_stations,-1);
			scenario_bike_lb.resize(scenarios,0.0);
		}
		void SetMatrixTargetUbLevels(int scenarios, int nb_stations) 
		{ 
			matrix_target_ub_levels.resize(scenarios,std::vector<int>());
			for(int e=0;e<scenarios;e++)
				matrix_target_ub_levels[e].resize(nb_stations,9999);
			scenario_bike_ub.resize(scenarios,0.0);
		}
		void SetTarget(int i, int target){ nodes[i].target = target;}
		void SetTarget(int sce, int i, int target)
		{
			if(sce>(int)matrix_target_levels.size() || i>(int)matrix_target_levels[sce].size() ) 
			{
				printf("sce:%d i:%d tgt:%d\n",sce,i,target); 
				printf("In Solution found Phill Collins (1989). Exiting ...\n "); exit(1);
			}
			matrix_target_levels[sce][i]=target;
		}
		void SetTargetLb(int sce, int i, int target_lb)
		{
			if(sce>(int)matrix_target_lb_levels.size() || i>(int)matrix_target_lb_levels[sce].size() )
			{
				printf("In Solution found Phill Collins (1989). Exiting ...\n "); exit(1);
			}				
			matrix_target_lb_levels[sce][i]=target_lb;
		}
		void SetTargetUb(int sce, int i, int target_ub)
		{
			if(sce>(int)matrix_target_ub_levels.size() || i>(int)matrix_target_ub_levels[sce].size() )
			{
				printf("In Solution found Phill Collins (1989). Exiting ...\n "); exit(1);
			}				
			matrix_target_ub_levels[sce][i]=target_ub;
		}
		int GetTarget(int i){ return nodes[i].target;}
		int GetTarget(int sce, int i)
		{
			if (sce >= int(matrix_target_levels.size()) || i >= (int)(matrix_target_levels[sce].size()))
			{
				printf("sce:%d i:%d\n",sce,i); 
				printf("In Solution found Phill Collins (1989). Exiting ...\n "); exit(1);
			}				
			return matrix_target_levels[sce][i];
		}
		int GetTargetLb(int sce, int i)
		{
			if(sce>(int)matrix_target_lb_levels.size() || i>(int)matrix_target_lb_levels[sce].size() ) 
			{
				printf("In Solution found Phill Collins (1989). Exiting ...\n "); exit(1);
			}
			return matrix_target_lb_levels[sce][i];
		}
		int GetTargetUb(int sce, int i)
		{
			if(sce>(int)matrix_target_ub_levels.size() || i>(int)matrix_target_ub_levels[sce].size() ) 
			{
				printf("In Solution found Phill Collins (1989). Exiting ...\n "); exit(1);
			}
			return matrix_target_ub_levels[sce][i];
		}
		void SetObjective(int scenario, double value){ scenario_objective[scenario]=value;}
		double GetObjective(int scenario){ return scenario_objective[scenario];}
		void SetBikeLb(int scenario, double value){ scenario_bike_lb[scenario]=value;}
		double GetBikeLb(int scenario){ return scenario_bike_lb[scenario];}
		void SetBikeUb(int scenario, double value){ scenario_bike_ub[scenario]=value;}
		double GetBikeUb(int scenario){ return scenario_bike_ub[scenario];}
		
		void SetStationLbVec(int i){station_bike_lb.resize(i,-1);}
		void SetStationUbVec(int i){station_bike_ub.resize(i,-1);}
		
		void SetStationLb(int station, int value){ station_bike_lb[station]=value;}
		int GetStationLb(int station){ return station_bike_lb[station];}
		void SetStationUb(int station, int value){ station_bike_ub[station]=value;}
		int GetStationUb(int station){ return station_bike_ub[station];}		
		
		void SetExactTarget(int t){exact_targets.push_back(t);}
		void SetHp(int hp){ hp_vec.push_back(hp); }
		void SetHm(int hm){ hm_vec.push_back(hm); }
		std::vector<int> & GetHpVec(){ return hp_vec; }
		std::vector<int> & GetHmVec(){ return hm_vec; }
		
		int GetExactTarget(int i){return exact_targets[i];}
		void SetWSTarget(int t){WS_targets.push_back(t);}
		int GetWSTarget(int i){return WS_targets[i];}
		std::vector<int> & GetExactTargets(){ return exact_targets; }
		std::vector<int> & GetWSTargets(){ return WS_targets; }
		
		void WriteTargetToFile(const char* filename) 
		{
			std::ofstream outputFile(filename); 
			if (outputFile.is_open()) 
			{
				int nb=0;
				for (const std::vector<int>& row : matrix_target_levels)
				{
					outputFile << std::fixed << std::setprecision(2) << scenario_objective[nb] <<'\n'; 
					for (int value : row) {
						outputFile << value << ' '; 
					}
					outputFile << '\n';
					nb++;
				}
				outputFile.close();  // Close the file
				std::cout << "Data has been written to " << filename << std::endl;
			} else {
					std::cerr << "Unable to open the file: " << filename << std::endl;
				}
		}
		void WriteExactTargetsToFile(const char* filename) 
		{
			std::ofstream outputFile(filename); 
			if (outputFile.is_open()) 
			{
				for (int i=0;i<exact_targets.size();i++)
					outputFile << exact_targets[i] <<'\n'; 
				outputFile.close();  // Close the file
				std::cout << "Exact targets have been written to " << filename << std::endl;
			} else {
					std::cerr << "Unable to open the file: " << filename << std::endl;
				}
		}
		void SetIndividualBoundsTime(double t){lb_i_ub_i_time=t;}
		void SetTotalBoundsTime(double t){total_bounds_time=t;}
		void SetRe(int _re){re=_re;}
		int GetRe(){return re;}
		void SetCplexNodes(int n){cplex_nodes=n;}
		void SetNbBendersCuts(int _nb_benders_cuts){nb_benders_cuts=_nb_benders_cuts;}
		void SetExactUb(double _exact_ub){exact_ub = _exact_ub;}
		double GetExactUb(){return exact_ub;}
		void SetExactLb(double _exact_lb){exact_lb = _exact_lb;}
		double GetExactLb(){return exact_lb;}
		
		void SetExactTime(double _exact_time){exact_time = _exact_time;}
		
		void SetUnsatisfiedRequests(const std::vector<int> & vec){
			unsastisfied_requests_fromExact = vec; 
			avg_exact_unsastifiedReq = std::accumulate(unsastisfied_requests_fromExact.begin(), unsastisfied_requests_fromExact.end(), 0.0) / (double)unsastisfied_requests_fromExact.size();
		}
		double GetExactUnsatisfiedReq(){return avg_exact_unsastifiedReq;}
		void SetUnsatisfiedRequestsFromWs(const std::vector<int> & vec){
			unsastisfied_requests_fromWs = vec;
			avg_Ws_unsastifiedReq = std::accumulate(unsastisfied_requests_fromWs.begin(),unsastisfied_requests_fromWs.end(), 0.0) / unsastisfied_requests_fromWs.size();
		}
		double GetWsUnsatisfiedReq(){return avg_Ws_unsastifiedReq;}
		
		void WriteReFile(const char * filename)
		{
			int requiredBikes=0; double avg_fill_rate = 0.0; 
			std::vector<double> fill_rates;
			for(int i=0;i<exact_targets.size();i++) 
			{
				requiredBikes += exact_targets[i]; 
				if (_prob->GetNode(i)->stationcapacity > 0) {
							fill_rates.push_back((double)exact_targets[i] / (double)_prob->GetNode(i)->stationcapacity);
						}				
			}
			std::sort(fill_rates.begin(), fill_rates.end());
			// Calculate the indices for the different quantiles
			int idx25 = (int)fill_rates.size() * 0.25;
			int idx50 = (int)fill_rates.size() * 0.50;
			int idx75 = (int)fill_rates.size() * 0.75;			

			// Count how many exact_targets fall into each quantile
			int count25 = 0, count50 = 0, count75 = 0, count100 = 0;
			for(int i = 0; i < exact_targets.size(); i++) {
				if (_prob->GetNode(i)->stationcapacity > 0) {
					double fill_rate = (double)exact_targets[i] / (double)_prob->GetNode(i)->stationcapacity;
					if (fill_rate <= fill_rates[idx25]) count25++;
					else if (fill_rate <= fill_rates[idx50]) count50++;
					else if (fill_rate <= fill_rates[idx75]) count75++;
					else count100++;
				}
			}			
			double cntr = count25 + count50 + count75 + count100;
			double count25p100 = (double)count25/cntr; double count50p100 = (double)count50/cntr; double count75p100 = (double)count75/cntr; double count100p100 = (double)count100/cntr;
			
			double avg_real_unsastifiedReq = std::accumulate(unsastisfied_requests_realTrips.begin(), unsastisfied_requests_realTrips.end(), 0.0) / unsastisfied_requests_realTrips.size();
			
			int nb_increased = std::accumulate( hp_vec.begin(), hp_vec.end(), 0.0 ); int nb_decreased = std::accumulate( hm_vec.begin(), hm_vec.end(), 0.0 );
			
			std::ofstream outputFile(filename);
			if(outputFile.is_open())
			{
				outputFile << Parameters::GetSolver() << ";"<< _prob->GetNodeCount() << ";" << _prob->GetQtot() << ";" << scenario_count << ";" << ODtrip_count << ";" << std::fixed << std::setprecision(2) << lb_i_ub_i_time << ";" << total_bounds_time << ";" << _prob->GetUpperBound() << ";" << avg_Ws_unsastifiedReq << ";" ;
				
				outputFile << re << ";";
				
				outputFile << std::fixed << std::setprecision(2) << exact_ub << ";" << exact_lb << ";" << exact_time << ";" << requiredBikes << ";" << avg_exact_unsastifiedReq << ";" << std::fixed << std::setprecision(2) << count25p100 << ";" << count50p100 << ";" << count75p100 << ";" << count100p100 << ";" ; 
				
				outputFile << nb_benders_cuts << ";" << cplex_nodes << ";";
				
				//if(Parameters::UnbiasedEstimator()) //Main.cpp will always compute the unbiased estimator now ...
				{
					outputFile << std::fixed << std::setprecision(2) << Avg_trips_400 << ";" << WS_in400 << ";" << UnbiasedEstimator_WS << ";" << UnbiasedEstimator_Exact << ";" << UnbiasedEstimator_DP << ";" << UnbiasedEstimator_LB << ";" << UnbiasedRejectedPicks << ";" << UnbiasedRejectedDrops << ";";
				}
				if(Parameters::CalculateAvgScenario())
				{
					outputFile << std::fixed << std::setprecision(2) << Avg_WS << ";" << Avg_Exact << ";" << Avg_DP << ";";
				}				
				if(Parameters::GetModel() == 5 || Parameters::GetModel() == 6 || Parameters::GetModel() == 7)
				{
					outputFile << nb_increased << ";" << nb_decreased << ";";
				}
				outputFile << "\n";
				
				std::cout << "The Re file has been written to " << filename << std::endl;
				outputFile.close();  // Close the file
			} else { 
				std::cerr << "Unable to open the file " << filename << std::endl;
			}
			
		}
		
		//Setting the W,Z vars from the exact class
		void SetW(int w){w_vars.push_back(w);}
		int GetW(int i){return w_vars[i];}
		void SetZ(int z){z_vars.push_back(z);}
		int GetZ(int i){return z_vars[i];}
		
		void SetWOmegaVars(int e){w_omega_vars.resize(e);}
		void SetZOmegaVars(int e){z_omega_vars.resize(e);}
		void SetWOmega(int e, int w){w_omega_vars[e].push_back(w);}
		int GetWOmega(int e, int i){return w_omega_vars[e][i];}
		
		void SetZOmega(int e, int z){z_omega_vars[e].push_back(z);}
		int GetZOmega(int e, int i){return z_omega_vars[e][i];}
		
		//For the real trips
		void SetNbRealTrips(int i){ nb_real_trips = i;}
		void SetNbRealDays(int i){nb_real_days = i;}
		void SetRecourseCost(double d){recourse_cost = d;}
		void SetRealUnsatisfiedRequests(std::vector<int> vec){
			unsastisfied_requests_realTrips = vec;
			avg_real_unsastifiedReq = std::accumulate(unsastisfied_requests_realTrips.begin(), unsastisfied_requests_realTrips.end(), 0.0) / unsastisfied_requests_realTrips.size();
		}
		double GetRealUnsatisfiedReq(){return avg_real_unsastifiedReq;}
		
		//UnbiasedEstimator
		void SetAvgTrips400(double d){ Avg_trips_400 = d; }
		void SetWS400(double d){ WS_in400=d; }
		void SetUnbiasedWS(double d){ UnbiasedEstimator_WS=d;}
		void SetUnbiasedExact(double d){ UnbiasedEstimator_Exact=d ;}
		double GetUnbiasedExact(){ return UnbiasedEstimator_Exact;}
		void SetUnbiasedLb(double d){ UnbiasedEstimator_LB=d; }
		double GetUnbiasedLb(){ return UnbiasedEstimator_LB; }
		void SetUnbiasedDP(double d){ UnbiasedEstimator_DP=d ;}
		
		void SetUnbiasedRejectedPicks(double d){ UnbiasedRejectedPicks = d;}
		void SetUnbiasedRejectedDrops(double d){ UnbiasedRejectedDrops = d;}
		double GetUnbiasedRejectedPicks(){ return UnbiasedRejectedPicks; }
		double GetUnbiasedRejectedDrops(){ return UnbiasedRejectedDrops; }
		
		//AvgScenario
		void SetAvgWS(double d){ Avg_WS=d;}
		void SetAvgExact(double d){ Avg_Exact=d ;}
		void SetAvgDP(double d){ Avg_DP=d ;}		

		//From the Lb,Ub methods
		std::vector<std::vector<int>> matrix_target_levels;
		std::vector<std::vector<int>> matrix_target_lb_levels;
		std::vector<std::vector<int>> matrix_target_ub_levels;
		
	private:
		Prob * _prob;
		std::vector<Node> nodes;
		
		int scenario_count; int ODtrip_count;
		
		//For z,w vars from the exact class
		std::vector<int> z_vars;
		std::vector<int> w_vars;
		std::vector<std::vector<int>> z_omega_vars;
		std::vector<std::vector<int>> w_omega_vars;
		
		std::vector<double> scenario_objective;
		std::vector<double> scenario_bike_lb;
		std::vector<double> scenario_bike_ub;
		std::vector<int> station_bike_lb;
		std::vector<int> station_bike_ub;
		
		//Targets
		std::vector<int> WS_targets;
		std::vector<int> exact_targets;
		std::vector<std::vector<int>> scenario_targets;
		
		std::vector<int> hp_vec;
		std::vector<int> hm_vec;
		
		double avg_exact_unsastifiedReq = -1.0; double avg_Ws_unsastifiedReq=-1.0; double avg_real_unsastifiedReq =-1.0;
		
		//For the real trips
		int nb_real_trips=-1;
		int nb_real_days=-1;
		double recourse_cost=-1.0;
		std::vector<int> unsastisfied_requests_fromWs;
		std::vector<int> unsastisfied_requests_fromExact;
		std::vector<int> unsastisfied_requests_realTrips;
		
		double lb_i_ub_i_time;
		double total_bounds_time;
		
		//From the exact class
		int nb_benders_cuts;
		int re;
		int cplex_nodes;
		double exact_ub;
		double exact_lb;
		double exact_time;

		//UnbiasedEstimators
		double WS_in400;
		double UnbiasedEstimator_WS;
		double UnbiasedEstimator_Exact;
		double UnbiasedEstimator_DP;
		double UnbiasedEstimator_LB;
		double Avg_trips_400;
		double UnbiasedRejectedPicks;
		double UnbiasedRejectedDrops;
		//AvgScenario
		double Avg_WS;
		double Avg_Exact;
		double Avg_DP;

};

#endif
