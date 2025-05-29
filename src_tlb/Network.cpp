
#include "Network.h"


Network::Network(int nbnodes): nodes(nbnodes)
{

}
Network::Network(std::vector<int>demands): nodes((int)demands.size())
{
	for(int i=0;i<nodes.size();i++)
	{
		nodes[i].demand = demands[i];
		nodes[i].id = i;
	}
		
}

Network::~Network()
{
	for(size_t i = 0;i<arcs.size();i++)
		delete arcs[i];
	arcs.clear();
	nodes.clear();
	map_arcs.clear();
}

//add an Arc
void Network::AddArc(int from, int to)
{
	AddArc(from, to, (int)arcs.size());
}
void Network::AddArc(int from, int to, int id)
{
	//printf("networkArc from:%d to:%d id:%d\n", from,to,id);
	network_arc_ptr arc = new network_arc_t();
	arc->id = (int)arcs.size();
	arc->from = from;
	arc->to = to;
	arc->id = id;
	arcs.push_back(arc);
	nodes[from].arcs.push_back(arc);
	map_arcs[id] = arc;
	
}

//set the basic information of an arc
void Network::SetArc(int id, char can_visit, double value)
{
	network_arc_ptr arc = map_arcs[id];
	arc->can_visit = can_visit;
	arc->value = value;
	arc->id = id;
}
void Network::SetArc(int id, char can_visit, double cap, double cost)
{
	network_arc_ptr arc = map_arcs[id];
	arc->can_visit = can_visit;
	arc->cap = cap;
	arc->cost = cost;
	arc->id = id;
	//printf("setNetworkArc from:%d to:%d id:%d cap:%d cost:%d\n",arc->from,arc->to,id,arc->cap,arc->cost);
}

void Network::SetArc(int id, char can_visit, int cap, double cost, int lower_capacity)
{
	network_arc_ptr arc = map_arcs[id];
	arc->can_visit = can_visit;
	arc->cap = cap;
	arc->cost = cost;
	arc->id = id;
	arc->lower_capacity = lower_capacity;
	printf("setNetworkArc from:%d to:%d id:%d cap:%d lower_capacity:%d cost:%d\n",arc->from,arc->to,id,arc->cap,arc->lower_capacity,arc->cost);
}

//set the basic information of an arc
void Network::SetArc(int id, char can_visit)
{
	network_arc_ptr arc = map_arcs[id];
	arc->can_visit = can_visit;
}

void Network::CloseOutGoingArcs(int from)
{
	network_node_ptr node = &nodes[from];
	for(int i=0;i<node->GetArcCount();i++)
		node->GetArc(i)->can_visit = 0;
}

network_arc_ptr Network::GetArc(int id)
{
	return map_arcs[id];
}

network_node_ptr Network::GetNode(int from)
{
	return &nodes[from];
}

//run a Breadth-First Search from an origin node and returns the number of marked arcs
void Network::BFS(int from)
{
	network_bfs_report_t report;
	BFS(from, &report);
}
void Network::BFS(int from, network_bfs_report_ptr report)
{
	for(size_t i = 0;i<arcs.size();i++)
		arcs[i]->init();
	for(size_t i = 0;i<nodes.size();i++)
		nodes[i].init();

	report->init();
	report->nbnodes = (int)nodes.size();
	report->nbarcs = (int)arcs.size();
	std::vector<network_node_ptr> to_visit;
	to_visit.push_back( &nodes[from] );

	int nb_visited = 0;
	while(nb_visited < (int)to_visit.size())
	{
		network_node_ptr node = to_visit[nb_visited];
		node->visited = 1;
		node->marked = 1;

		for(int i=0;i<node->GetArcCount();i++)
		{
			network_arc_ptr arc = node->GetArc(i);
			if(arc->can_visit == 0) continue;

			arc->marked = 1;

			if(nodes[ arc->to ].marked == 0)
			{
				nodes[ arc->to ].marked = 1;
				to_visit.push_back( &nodes[ arc->to ]);
			}
		}
		nb_visited++;
	}
	//printf("nbvisited:%d\n", nb_visited);
	for(size_t i = 0;i<arcs.size();i++)
	{
		report->nb_marked_arcs += arcs[i]->marked;
		if(arcs[i]->value > 0.0001) report->nb_positive_arcs++;
		if(arcs[i]->value > 0.0001 && arcs[i]->can_visit == 1)
			report->nb_marked_positive_arcs+= arcs[i]->marked;
	}

	for(size_t i = 0;i<nodes.size();i++)
		report->nb_marked_nodes += nodes[i].marked;
}

bool Network::IsVisited(int node)
{
	return nodes[node].visited == 1;
}

void Network::CalculateOutArcs()
{
	for(size_t i = 0;i<nodes.size();i++)
		nodes[i].out_arcs_count = 0;
	for(size_t i = 0;i<nodes.size();i++)
	{
		network_node_ptr node = &nodes[i];
		for(int j=0;j<node->GetArcCount();j++)
		{
			network_arc_ptr arc = node->GetArc(j);
			nodes[arc->to].out_arcs_count++;
		}
	}
}

void Network::PrintGraphViz(char * filename)
{
	CalculateOutArcs();
	FILE * f = fopen(filename, "w");
	if(f == NULL) return;

	fprintf(f,"digraph G {\n");
 	fprintf(f,"size = \"1000, 1000\";\n");
 	fprintf(f,"overlap = scale;\n");
 	fprintf(f,"splines = true;\n");
 	for(size_t i = 0;i<nodes.size();i++)
	{
		if(nodes[i].GetArcCount() >= 1 || nodes[i].out_arcs_count >= 1)
		{
			//nodes[i].Show();
			//fprintf(f,"%d [pos = \"%d,%d!\"]\n", (int)i,0,(int)i);
			// Use the 'label' attribute to add the demand as a label
            //fprintf(f, "%d [pos = \"%d,%d!\", label = \"%d\"]\n", (int)i, 0, (int)i, nodes[i].demand);
			fprintf(f, "%d [pos = \"%d,%d!\", label = \"%d (%d)\"]\n", (int)i, 0, (int)i,(int)i, nodes[i].demand);
			//printf("%d [pos = \"%d,%d!\", label = \"%d (Dmd:%d)\"]\n", (int)i, 0, (int)i,(int)i, nodes[i].demand);
		}
	}
	for(size_t i = 0;i<nodes.size();i++)
	{
		network_node_ptr node = &nodes[i];
		for(int j=0;j<node->GetArcCount();j++)
		{
			network_arc_ptr arc = node->GetArc(j);
			//For the flow graph
			if(is_flow_graph)
				fprintf(f,"%d -> %d[label=\"%.2lf x %d\"];\n", arc->from,arc->to,arc->value, (int) arc->cost);
			//For the input graph
			if(is_input_graph)
			{
				//OG
				//fprintf(f,"%d -> %d[label=\"%.0lf[0,%.0lf]\"];\n",arc->from,arc->to,arc->cost,arc->value);
				//arc->Show();
				fprintf(f,"%d -> %d[label=\"%d[%d,%d]\"];\n",arc->from,arc->to,arc->cost,arc->lower_capacity,arc->cap);
				//fprintf(f,"%d -> %d[label=\"%d[0,%d]\"];\n",arc->from,arc->to,arc->cost,arc->cap);
				//printf("%d -> %d[label=\"%d[0,%d]\"];\n",arc->from,arc->to,(int)arc->cost,arc->cap);
			}
				
				
		}
	}

	fprintf(f,"}\n");
	fclose(f);
}
