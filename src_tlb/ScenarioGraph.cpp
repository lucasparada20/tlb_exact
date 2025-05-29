#include "ScenarioGraph.h"
#include <unordered_map>
#include <functional> // for std::hash
#include <deque>


void ScenarioGraph::MakeMCFGraph(Scenario * scenario_object) {
    
	scenario_no = scenario_object->GetScenarioNo();
	scenario_object = scenario_object;
	nb_network_arcs = 0;
	nb_trip_arcs = 0;
	nb_nodes = 2;
    arcs.clear();
	
	//printf("MakeGraph scenario:%d\n",scenario_no);
	
	//Source
	deficits.push_back(-1*prob->GetQtot()); 
	//Sink
	deficits.push_back(prob->GetQtot());
	// Create nodes at t = -1 
   	//NETWORK ARCS
	for(int i=0; i<prob->GetNodeCount(); i++)
	{
		Node * station = prob->GetNode(i);
		
		MCF_arc start_depot_arc; start_depot_arc.from = 1; start_depot_arc.to = i+3; start_depot_arc.cap = station->stationcapacity; start_depot_arc.type= 1;
		//start_depot_arc.from_cust_no = 10000; start_depot_arc.to_cust_no = i;
		arcs.push_back(start_depot_arc);
		//start_depot_arc.Show();
		nb_nodes++;
	}
	
	//std::vector<std::vector<int>> times(prob->GetNodeCount(), std::vector<int>(0));
	times.resize(prob->GetNodeCount(), std::vector<int>(0));
	scenario_object->CalculateEventTimes(times);
	//Uncomment to print the event times ...
	/*printf("Scenario:%d\n",scenario_no);
	for (int i = 0; i < prob->GetNodeCount(); i++) 
	{
		printf("sta:%d times:",i);
		for (int k = 0; k < times[i].size(); k++) 
			printf("%d ",times[i][k]);
		printf("\n");
	}*/
	nb_nodes++;
	//printf("nb_nodes:%d before creating ALL nodes\n",nb_nodes);
	//Source is nb_nodes = 1 and sink is nb_nodes = 2
	std::map<std::pair<int, int>, int> node_idx;
	std::map<int,std::pair<int,int>> reverse_node_idx;
	for (int i = 0; i < prob->GetNodeCount(); i++) {
		node_idx[{i,-1}] = i+3;
		reverse_node_idx[i+3] = {i,-1};
		for (int k = 0; k < times[i].size(); k++) {
			node_idx[{i, times[i][k]}] = nb_nodes;
			reverse_node_idx[nb_nodes] = {i, times[i][k]};
			nb_nodes++;
		}
		node_idx[{i,prob->GetTmax()+1}] = 2;
		reverse_node_idx[2] = {i,prob->GetTmax()+1};
	}
	nb_nodes--;
	//printf("ScenarioGraph: sce:%d nb_nodes:%d after creating ALL nodes\n",scenario_no,nb_nodes);
	// Print the contents of node_idx
	/*for (const auto& entry : node_idx) {
		std::cout << "Node Index for (" << entry.first.first << ", " << entry.first.second << "): " << entry.second << std::endl;
	}*/

	for (auto it = node_idx.begin(); std::next(it) != node_idx.end(); ++it) {
		const auto& currentPair = it->first;
		const auto& nextPair = std::next(it)->first;

		int currentNodeIdx = it->second;
		int nextNodeIdx = std::next(it)->second;

		if (currentNodeIdx == 2) continue;

		// Create an arc between current and next nodes
		MCF_arc arc;
		arc.from = currentNodeIdx;
		arc.to = nextNodeIdx;
		arc.cap = prob->GetNode(currentPair.first)->stationcapacity; 
		arc.type = 2; 
		//arc.from_cust_no = currentPair.first; 
		//arc.to_cust_no = nextPair.first; 
		//printf("%d %d\n",currentPair.first,nextPair.first);
		arc.cust_no = currentPair.first;
		//arc.Show(); 
		//getchar();
		arcs.push_back(arc);
	}
	MCF_arc arc; arc.from = 1; arc.to = 2; arc.cap = prob->GetQtot(); arc.type=1; 
	//arc.from_cust_no = 10000; arc.to_cust_no = 10000;
	arcs.push_back(arc);

    nb_network_arcs = arcs.size();
	

	// Important: Stations in the instances start from 1 (station 0 is the depot in the .txt file)
	// Need to store (and later search) for Node+1 ...
	// Process ODTrips
	nb_trip_arcs = 0;
	for (int i = 0; i < scenario_object->GetODTripCount(); i++) {
		Trip* trip = scenario_object->GetODTrip(i);

		// Convert trip->start_no and trip->end_no once outside the loop
		int start_no = trip->start_no - 1;
		int end_no = trip->end_no - 1;

		auto from_node_iter = node_idx.find({start_no, trip->start_t});
		auto to_node_iter = node_idx.find({end_no, trip->end_t});

		if ((Parameters::GetInstanceFormat() == 'T' && (from_node_iter == node_idx.end() || to_node_iter == node_idx.end()))) {
			// Keys NOT found, handle the error
			
			std::cerr << "start_no:" << start_no << " end_no:" << end_no << " Error: Keys not found in ODtrip_nodeMap\n";
			trip->Show();
			exit(1);
		}

		MCF_arc arc;
		arc.from = from_node_iter->second;
		arc.to = to_node_iter->second;
		arc.cap = 1;
		arc.type = 3;
		
		arc.from_cust_no = reverse_node_idx[ from_node_iter->second ].first;
		arc.to_cust_no = reverse_node_idx[ to_node_iter->second ].first;
		
		arcs.push_back(arc);
		//arc.Show();
		nb_trip_arcs++;
	}

	//for(int i=0;i<arcs.size();i++)
		//arcs[i].idx = i; 		

	arcs.resize( arcs.size() );
	
	//Some output ...
	//printf("MakeGraph sce:%d nodes:%d network arcs:%d Trip arcs:%d Total Arcs in MCF_graph:%d\n",e,(int)nodes.size(),network_arcs,(int)arcs.size()-network_arcs,(int)arcs.size());
	//getchar();
	//adjust nb_nodes
	
	for(int i=0;i<arcs.size();i++)
		if((arcs[i].from<0 ||arcs[i].from>nb_nodes) || (arcs[i].to<0 ||arcs[i].to>nb_nodes))
			{ printf("Wrong 3d arc to MCF arc... Exiting\n"); arcs[i].Show(); exit(1); }
		
	//------ For ParetoOpt ------
	OutArcs.resize(nb_nodes,std::vector<int>(0));
	InArcs.resize(nb_nodes,std::vector<int>(0));
	for(int i=0;i<arcs.size();i++)
	{
		InArcs[ arcs[i].to-1 ].push_back(i);
		OutArcs[ arcs[i].from-1 ].push_back(i);
	}
	
	/*for(int i=0;i<nb_nodes;i++)
	{
		printf("Node %d InArcs: ",i);
		for(int k=0;k<InArcs[i].size();k++)
			printf("%d ",InArcs[i][k]);
		printf("\n");
		printf("Node %d OutArcs: ",i);
		for(int k=0;k<OutArcs[i].size();k++)
			printf("%d ",OutArcs[i][k]);
		printf("\n");

	}*/
		
	//------ ------------- ------		
	
	//ShortestPathToSink();
	//ShortestPathFromSource();

}

void ScenarioGraph::StoreInSol(Sol * sol)
{
	for(int i=0; i<prob->GetNodeCount();i++)
		sol->SetTarget(scenario_object->GetScenarioNo(), i, flows[i]);
	sol->SetObjective(scenario_object->GetScenarioNo(), objective_value);	
}

void ScenarioGraph::PrintGraph(char * filename)
{	
	//Need to test for Lemon solvers ...
	//Recall that node 0 needs to have index one in Unipi Solvers. So need to -=1 to all arc indices.
	std::vector<int> demands(nb_nodes,0);
	for(int i=0;i<deficits.size();i++)
		demands[i] = -1*deficits[i];
	
	Network * n = new Network(demands);
	n->SetInputGraph();
	
	for(size_t i=0;i<arcs.size();i++)
		n->AddArc(arcs[i].from-1,arcs[i].to-1,i);
	for(size_t i=0;i<arcs.size();i++)
		if(arcs[i].type == 3)
			n->SetArc(i,1,arcs[i].cap,-1);
		else n->SetArc(i,1,arcs[i].cap,0);
	n->PrintGraphViz(filename);
	delete n;
}

void ScenarioGraph::PrintGraph(char * filename, std::vector<int> & deficits, int nb_graph_nodes)
{	
	for(int i=0;i<deficits.size();i++)
		printf(" %d:%d",i,deficits[i]);
	printf("\n");
	//printf("PrintGraph Arcs:%d Nodes:%d\n",(int)arcs.size(),nb_graph_nodes);
	//Need to test for Lemon solvers ...
	Network * n = new Network(deficits);
	n->SetInputGraph();
	
	//If calling from Lemon, nodes already starting from 0,1,...
	for(size_t i=0;i<arcs.size();i++)
	{
		//arcs[i].Show();
		n->AddArc(arcs[i].from,arcs[i].to,i);
	}
	for(int i=0;i<arcs.size();i++)
	{
		if(arcs[i].type == 1 && arcs[i].to > 1 && arcs[i].to <= prob->GetNodeCount() + 2)
		{
			n->SetArc(i, 1, arcs[i].cap, 0, arcs[i].cap);
			arcs[i].Show(); 
			printf("tgt[%d]:%d\n",arcs[i].to, deficits[ arcs[i].to ]);
			getchar();
		}
			
		else
		{
			n->SetArc(i, 1, arcs[i].cap, arcs[i].type ==3 ? -1 : 0, 0 );
			arcs[i].Show(); getchar();
		}
			
	}

	n->PrintGraphViz(filename);
	delete n;
}

void ScenarioGraph::ShortestPathFromSource() {
    clock_t start_time = clock();

    // The source node is 0
    int source = 0;

    // Initialize distances and parents
    dist.resize(nb_nodes + 1, 999999); 
    std::vector<std::vector<int>> parents(nb_nodes + 1, std::vector<int>(999999)); // To keep track of shortest paths

    // Adjacency list for the original graph
    std::vector<std::vector<std::pair<int, int>>> graph(nb_nodes + 1);
    for (const auto& arc : arcs) {
        int u = arc.from - 1; // 0-based indices to match Lemon graphs
        int v = arc.to - 1;   // 0-based indices to match Lemon graphs
        int cost = (arc.type == 3) ? -1 : 0; // Negative weight for type 3 arcs
        graph[u].emplace_back(v, cost); // Add the arc as-is (not reversed)
    }

    // Start BFS from the source node (node 0)
    std::deque<int> queue;
    dist[source] = 0;
    queue.push_back(source);

    while (!queue.empty()) {
        int u = queue.front();
        queue.pop_front();

        for (auto& neighbor : graph[u]) {
            int v = neighbor.first;
            int weight = neighbor.second;

            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                parents[v].clear(); // Clear old parents since we found a shorter path
                parents[v].push_back(u); // Record the parent
                if (weight == -1) {
                    queue.push_front(v);  // Negative weight: higher priority
                } else {
                    queue.push_back(v);   // Zero weight: normal priority
                }
            } else if (dist[u] + weight == dist[v]) {
                // If a path with equal distance is found, add this node to parents
                parents[v].push_back(u);
            }
        }
    }

    // Verifying triangle inequality for all reachable nodes with c_ij <= 0:
    // Triangle inequality: D_j >= D_i
    // Violated triangle inequality: D_j < D_i -> exit
    for (auto& arc : arcs) {
        if (arc.to == 0) continue; // Skip depot arcs
        int u = arc.from;  // remember to adjust for zero-based index
        int v = arc.to;    // remember to adjust for zero-based index
        int cost = (arc.type == 3) ? -1 : 0;

        if (dist[u - 1] != 999999 && dist[v - 1] != 999999) {
            if (dist[v - 1] > dist[u - 1] + cost) { // From needs to be ALWAYS BIGGER
                printf("Inconsistency found: D[%d] = %d > D[%d] = %d + %d\n",
                       v - 1, dist[v - 1], u - 1, dist[u - 1], cost);
                arc.Show();
                for (int i = 0; i <= nb_nodes; i++)
                    printf("D%d:%d\n", i, dist[i]);
                PrintGraph((char*)"SSgraph.dot");
                // getchar();
                exit(1);
            }
        }
    }

    clock_t end_time = clock();
    double elapsedSeconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("ShortestPathFromSource time:%.1lf scenario:%d nodes:%d arcs:%d trips:%d\n", elapsedSeconds, scenario_no, nb_nodes, (int)arcs.size(), nb_trip_arcs);
    for (int i = 0; i <= nb_nodes; i++)
        printf("D%d:%d\n", i, dist[i]);
}

void ScenarioGraph::ShortestPathToSink() {
    clock_t start_time = clock();

    // The sink node is 1
    int sink = 1;

    // Initialize distances and parents
    dist.resize(nb_nodes + 1, 999999);
    std::vector<std::vector<int>> parents(nb_nodes + 1);

    // Adjacency list for the reversed graph
    std::vector<std::vector<std::pair<int, int>>> reversedGraph(nb_nodes + 1);
    for (const auto& arc : arcs) {
        int u = arc.from - 1; // 0-based indices to match Lemon graphs
        int v = arc.to - 1;   // 0-based indices to match Lemon graphs
        int cost = (arc.type == 3) ? -1 : 0;
        reversedGraph[v].emplace_back(u, cost); // Reverse the arc direction
    }

    // Start BFS from the sink node (node 1)
    std::deque<int> queue;
    dist[sink] = 0;
    queue.push_back(sink);

    while (!queue.empty()) {
        int u = queue.front();
        queue.pop_front();

        for (auto& neighbor : reversedGraph[u]) {
            int v = neighbor.first;
            int weight = neighbor.second;

            if (dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                parents[v].clear(); // Clear old parents since we found a shorter path
                parents[v].push_back(u); // Record the parent
                if (weight == -1) {
                    queue.push_front(v);  // Negative weight: higher priority
                } else {
                    queue.push_back(v);   // Zero weight: normal priority
                }
            } else if (dist[u] + weight == dist[v]) {
                // If a path with equal distance is found, add this node to parents
                parents[v].push_back(u);
            }
        }
    }

    // Verifying triangle inequality for all reachable nodes with c_ij <= 0:
    for (auto& arc : arcs) {
        if (arc.from == 1) continue; // Skip depot arcs
        int u = arc.from;  // Adjust for zero-based index
        int v = arc.to;    // Adjust for zero-based index
        int cost = (arc.type == 3) ? -1 : 0;

        if (dist[u - 1] != 999999 && dist[v - 1] != 999999) {
            if (dist[v - 1] < dist[u - 1] + cost) {
                printf("Inconsistency found: D[%d] = %d > D[%d] = %d + %d\n",
                       v - 1, dist[v - 1], u - 1, dist[u - 1], cost);
                arc.Show();
                for (int i = 0; i <= nb_nodes; i++)
                    printf("D%d:%d\n", i, dist[i]);
                PrintGraph((char*)"SSgraph.dot");
                exit(1);
            }
        }
    }

    clock_t end_time = clock();
    double elapsedSeconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("ShortestPathToSink time:%.1lf scenario:%d nodes:%d arcs:%d trips:%d\n", elapsedSeconds, scenario_no, nb_nodes, (int)arcs.size(), nb_trip_arcs);
}
