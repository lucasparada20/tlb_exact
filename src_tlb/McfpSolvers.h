#ifndef MCFP_SOLVERS
#define MCFP_SOLVERS

#include "ScenarioGraph.h"
#include "Network.h"

//Headers for LEMON
#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>


using namespace lemon;

#define WS 1
#define GlobalTargetLbUb 2
#define TargetBounds 3
#define Recourse 4
#define OptStatCapRecourse 5
#define WSnoCap 6
#define SingleSourceRecourse 7
#define AllPositive 8
#define xZero 10

class McfpSolvers {
public:
    double Solve(ScenarioGraph * UserGraph, bool show, bool storeDuals, int graph_type, int bound_type = -1, int station_no = -1, const std::vector<int>& targets = std::vector<int>(), const std::vector<int>& hm = std::vector<int>(), const std::vector<int>& hp = std::vector<int>());
	
	void MakeLemonGraph(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & lowerCapacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, int graph_type, int bound_type = -1, int station_no = 1, const std::vector<int>& targets = std::vector<int>(), const std::vector<int>& hm = std::vector<int>(), const std::vector<int>& hp = std::vector<int>());
	
	void DebugInfGraph(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & lowerCapacities, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, int graph_type, ListDigraph::NodeMap<int> & nodePotentials);
	
	int GetDualCost(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, ListDigraph::NodeMap<int> & nodePotentials, const std::vector<int> & targets);
	
	int GetDualCostSingleSource(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, ListDigraph::NodeMap<int> & nodePotentials, const std::vector<int> & targets, int graph_type);	
	
	void PrintGraph(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, char * graphName);
	
	void StoreInSol(Sol * sol, int scenario);
	int GetFlow(int i){return flows[i];}
	double GetObjectiveValue(){return objective_value;}
	int GetNbSatisfiedTrips(){return numberOfSatisfiedTrips;}
	
	std::vector<int> y0;
	
	bool calculate_target_levels = false;
	void CalculateTargetLevels(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & costs,NetworkSimplex<ListDigraph> & NetSimplex, int scenario);
	
private:
	Prob * prob;
	std::vector<int> flows;
	double objective_value;
	int numberOfSatisfiedTrips;
};

#endif