#ifndef SCENARIO_GRAPH
#define SCENARIO_GRAPH

#include <cstdint>

#include <vector>
//#include "mathfunc.h"
#include "NodeTLB.h"
#include "ProblemDefinition.h"
#include "Parameters.h"
#include "Network.h"
#include "Solution.h"
#include "Scenarios.h"

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Simple hash combine function
        return h1 ^ h2;
    }
};

class MCF_arc{
public:
	int32_t from; int32_t to; int32_t cap; 
	//int idx; 
	int16_t from_cust_no; int16_t to_cust_no; 
	int16_t cust_no = -1; //Needed for Model() == 5
	int8_t type; // type 1 =  source outgoing; type 2 = holding; type 3 = trip;
	float dual = 0.0; // Duals need to be >= 0
	float lowerDual = 0.0; // Duals need to be >= 0
	void SetDual(double v){dual=v;};
	void Show(){
		printf("a_MCF from:%d to:%d cap:%d type:%d cust_no:%d dual:%.1lf\n",from,to,cap,type,cust_no,dual);
		//printf("a_MCF from:%d to:%d from_cust_no:%d to_cust_no:%d cap:%d type:%c idx:%d dual:%.1lf\n",from,to,from_cust_no,to_cust_no,cap,type,idx,dual);
	}
};

class ScenarioGraph
{
public:
	ScenarioGraph(Prob * _prob) : prob(_prob) {}
	ScenarioGraph(Prob * _prob, Scenario * _scenario_object) : prob(_prob), scenario_object(_scenario_object) {
		MakeMCFGraph(scenario_object);
	}
    // Copy constructor
    ScenarioGraph(const ScenarioGraph& other) : prob(other.prob), scenario_object(other.scenario_object) {
        MakeMCFGraph(scenario_object); // Recreate the MCF graph for the copied object
    }
	
	void MakeMCFGraph(Scenario * scenario_object);
	
	void PrintGraph(char * filename);
	void PrintGraph(char * filename, std::vector<int> & deficits, int nb_graph_nodes);
	void PrintFlowGraph(char * filename);
	
	int GetNodeCount(){return nb_nodes;}
	int GetArcCount(){return (int)arcs.size();}
	int GetNetworkArcCount(){return nb_network_arcs;}
	int GetODtripCount(){return nb_trip_arcs;}
	
	MCF_arc * GetArc(int i){return &arcs[i];}
	
	void SetDual(int i, double v){arcs[i].SetDual(v);}
	void SetDualObjectiveValue(double v) { dual_objective_value = v; }
	void SetObjectiveValue(double v) { objective_value = v; }
	double GetObjectiveValue() { return objective_value; }
	double GetDualObjectiveValue() { return dual_objective_value; }
	
	void StoreInSol(Sol * sol);
	
	void ShortestPathToSink(); //computes the shortest path from each node to the sink and store in dist.
	void ShortestPathFromSource(); //Dido, but from the source to each node and stores in dist.
	
	std::vector<int> deficits;
	std::vector<int> flows;
	int GetFlow(int i){return flows[i];}
	std::vector<int> pi;
	std::vector<int> dist;
	
	std::vector<std::vector<int>> OutArcs;
	std::vector<std::vector<int>> InArcs;
	int GetOutArcCount(int i){ return (int)OutArcs[i].size(); }
	int GetInArcCount(int i){ return (int)InArcs[i].size(); }
	int GetOutArc(int i, int k){ return OutArcs[i][k]; }
	int GetInArc(int i, int k){ return InArcs[i][k]; }
	
	//for AllPositive graph
	std::vector<int> lemonDeficits;
	//-----------------------
	
	Prob * GetProblem(){return prob;}
	
	int scenario_no;
	void ResetDuals(){ for(MCF_arc & a: arcs) {a.dual=0.0; a.lowerDual=0.0;}}
	void AddArc(MCF_arc & a){ arcs.push_back(a); }
	void SetNodeCount( int c ){ nb_nodes = c; }
	void GetDeficits( std::vector<int> & d ){ d = deficits; }
	int GetNbEventTimes(int station){ return (int)times[station].size(); }
	
private:
	
	Prob * prob;
	Scenario * scenario_object;
	std::vector<MCF_arc> arcs;
	std::vector<uint16_t> t;
	int nb_nodes;
	int nb_network_arcs;
	int nb_trip_arcs;
	double objective_value; 
	double dual_objective_value;
	std::vector<std::vector<int>> times;
};

#endif