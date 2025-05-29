/*
 * Copyright Jean-Francois Cote 2012
 *
 * The code may be used for academic, non-commercial purposes only.
 *
 * Please contact me at cotejean@iro.umontreal.ca for questions
 *
 * If you have improvements, please contact me!
*/


#ifndef NETWORKH
#define NETWORKH

#include <map>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

struct network_arc_t
{
	int id;
	int from;
	int to;
	char marked;
	double value;
	int cost;
	char can_visit;
	int cap;
	int lower_capacity;
	network_arc_t() : id(-1),from(-1),to(-1),marked(0),value(0),can_visit(1),cost(0),cap(-1), lower_capacity(-1) {}

	void Show(){printf("NetworkArc from:%d to:%d val:%.1lf lowerCap:%d upperCap:%d\n",from,to,value,lower_capacity,cap);}
	void init()
	{
		marked = 0;
	}
} ;
typedef network_arc_t* network_arc_ptr;

struct network_node_t
{
	int id;
	std::vector<network_arc_ptr> arcs;
	char visited;
	char marked;
	int out_arcs_count;
	int demand;
	network_node_t():id(0),arcs(0),visited(0),marked(0),out_arcs_count(0),demand(0){}
	
	void Show(){printf("Network node id:%d dmd:%d arcs:%d out_arcs:%d\n",id,demand,(int)arcs.size(),out_arcs_count);};
	
	int GetArcCount(){return (int)arcs.size();}
	network_arc_ptr GetArc(int i){return arcs[i];}

	void init()
	{
		visited = 0;
		marked = 0;
		out_arcs_count = 0;
	}
};
typedef network_node_t* network_node_ptr;

struct network_bfs_report_t
{
	int nb_marked_nodes;
	int nb_marked_arcs;
	int nb_positive_arcs;
	int nb_marked_positive_arcs;
	int nbnodes;
	int nbarcs;

	network_bfs_report_t():nb_marked_nodes(0),nb_marked_arcs(0),nb_positive_arcs(0),nb_marked_positive_arcs(0),
							nbnodes(0),nbarcs(0){}
	void init()
	{
		nbnodes = 0;
		nbarcs = 0;
		nb_marked_nodes = 0;
		nb_marked_arcs = 0;
		nb_marked_positive_arcs = 0;
		nb_positive_arcs = 0;
	}
	void print()
	{
		printf("nbnodes:%d nbarcs:%d nb_marked_nodes:%d nb_marked_arcs:%d nb_positive_arcs:%d nb_marked_positive_arcs:%d\n",
				nbnodes, nbarcs,nb_marked_nodes,nb_marked_arcs,nb_positive_arcs,nb_marked_positive_arcs);
	}
};
typedef network_bfs_report_t* network_bfs_report_ptr;

class Network
{
	public:
	Network(int nbnodes);
	Network(std::vector<int>demands);
	
	~Network();

	//add an Arc
	void AddArc(int from, int to,int id);
	void AddArc(int from, int to);

	int ArcCount(){return (int)arcs.size();}

	//set the basic information of an arc
	void SetArc(int id, char can_visit, double value);
	void SetArc(int id, char can_visit, double value, double cost);
	void SetArc(int id, char can_visit, int cap, double cost, int lower_capacity);
	void SetArc(int id, char can_visit);
	void CloseOutGoingArcs(int node);
	network_arc_ptr GetArc(int id);
	network_node_ptr GetNode(int from);
	void CalculateOutArcs();
	void PrintGraphViz(char * filename);

	//run a Breadth-First Search from an origin node
	void BFS(int from);
	void BFS(int from, network_bfs_report_ptr report);

	bool IsVisited(int node);

	int GetNodeCount(){return (int)nodes.size();}
	int GetArcCount(){return (int)arcs.size();}
	
	void SetFlowGraph(){is_flow_graph=true;}
	void SetInputGraph(){is_input_graph=true;}
	
	
	private:
		std::vector<network_node_t> nodes;
		std::vector<network_arc_ptr> arcs;
		std::map<int, network_arc_ptr> map_arcs;
		
		bool is_flow_graph = false;
		bool is_input_graph = false;
		char * graph_name;
};




#endif
