#include "McfpSolvers.h"
#include <iostream>
#include <ctime>
#include <omp.h>

double McfpSolvers::Solve(ScenarioGraph * UserGraph, bool show, bool storeDuals, int graph_type, int bound_type, int station_no, const std::vector<int> & targets, const std::vector<int>& hm, const std::vector<int>& hp) {
    
	//printf("In McfpSolvers sizeOf hm:%d hp:%d graph_type:%d\n",(int)hm.size(),(int)hp.size(),graph_type);
	prob = UserGraph->GetProblem();
	flows.clear();
	flows.resize(prob->GetNodeCount(),-1);
	objective_value = 0.0;	
	
	double startTime = omp_get_wtime();
	
    // Update the Lemon graph
    ListDigraph LemonGraph;
	ListDigraph::ArcMap<int> capacities(LemonGraph);
    ListDigraph::ArcMap<int> costs(LemonGraph);
	ListDigraph::ArcMap<int> lemon_flows(LemonGraph);
    ListDigraph::NodeMap<int> supplies(LemonGraph);
	ListDigraph::NodeMap<int> nodePotentials(LemonGraph);
	ListDigraph::ArcMap<int> lowerCapacities(LemonGraph);
	
	MakeLemonGraph(UserGraph,LemonGraph,capacities,lowerCapacities,costs,supplies,graph_type,bound_type,station_no,targets,hm,hp);
	
	switch (Parameters::GetSolver()) {
		case 3: {
			NetworkSimplex<ListDigraph> NetSimplex(LemonGraph);
			NetSimplex.upperMap(capacities);
			NetSimplex.costMap(costs);
			NetSimplex.supplyMap(supplies);
			
			if(graph_type == SingleSourceRecourse || graph_type == xZero) // LOOOOOOOOOOOOOOOOLLLLL!
				NetSimplex.lowerMap(lowerCapacities);
			
			NetworkSimplex<ListDigraph>::ProblemType status = NetSimplex.run();

			if (status == NetworkSimplex<ListDigraph>::OPTIMAL) {
				objective_value = -1*NetSimplex.totalCost();
				UserGraph->SetObjectiveValue(objective_value);

				for (int i = 0; i < prob->GetNodeCount(); i++) {
					ListDigraph::Arc lemonArc = LemonGraph.arcFromId(i);
					int flow = NetSimplex.flow(lemonArc);
					flows[i] = flow;

				}
				if (storeDuals) {
					NetSimplex.flowMap(lemon_flows);
					NetSimplex.potentialMap(nodePotentials);
				}
				
				numberOfSatisfiedTrips = 0;
				for(int i=UserGraph->GetNetworkArcCount();i<UserGraph->GetArcCount();i++)
				{
					ListDigraph::Arc lemonArc = LemonGraph.arcFromId(i);
					if(graph_type != xZero)
						numberOfSatisfiedTrips += NetSimplex.flow(lemonArc);
					else if( graph_type == xZero && NetSimplex.flow(lemonArc) > 0 )
						numberOfSatisfiedTrips++;
				}
				
			} else if (status == NetworkSimplex<ListDigraph>::UNBOUNDED) {
				std::cout << "MCF problem unbounded." << std::endl;
				UserGraph->PrintGraph((char*)"UnboundedGraph.dot");
				exit(1);
			} else if (status == NetworkSimplex<ListDigraph>::INFEASIBLE) {
				std::cout << "MCF problem infeasible." << std::endl;
				for(ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node)
					nodePotentials[node] = 9999;
				DebugInfGraph(UserGraph,LemonGraph,lowerCapacities,capacities,costs,supplies,graph_type,nodePotentials);
				exit(1);
			} else {
				std::cout << "Error in the MCF Lemon solver." << std::endl;
				exit(1);
			}
			if (calculate_target_levels)
				CalculateTargetLevels(UserGraph,LemonGraph,costs,NetSimplex,UserGraph->scenario_no);
		}
		break;

		case 4: {
			CostScaling<ListDigraph> CScaling(LemonGraph);
			CScaling.upperMap(capacities);
			CScaling.costMap(costs);
			CScaling.supplyMap(supplies);
			CostScaling<ListDigraph>::ProblemType status = CScaling.run();

			if (status == CostScaling<ListDigraph>::OPTIMAL) {
				objective_value = -1 * CScaling.totalCost();
				UserGraph->SetObjectiveValue(objective_value);
				
				for (int i = 0; i < prob->GetNodeCount(); i++) {
					ListDigraph::Arc lemonArc = LemonGraph.arcFromId(i);
					int flow = CScaling.flow(lemonArc);
					flows[i] = flow;						
				}

				if (storeDuals) {
					CScaling.flowMap(lemon_flows);
					CScaling.potentialMap(nodePotentials);
				}
			} else if (status == CostScaling<ListDigraph>::UNBOUNDED) {
				std::cout << "MCF problem unbounded." << std::endl;
				UserGraph->PrintGraph((char*)"UnboundedGraph.dot");
				exit(1);
			} else if (status == CostScaling<ListDigraph>::INFEASIBLE) {
				std::cout << "MCF problem infeasible." << std::endl;
				//UserGraph->PrintGraph((char*)"InfGraph.dot");
				for(ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node)
					nodePotentials[node] = 9999;
				//DebugInfGraph(UserGraph,LemonGraph,lowerCapacities,capacities,costs,supplies,graph_type,nodePotentials);
				exit(1);
			} else {
				std::cout << "Error in the MCF Lemon solver." << std::endl;
				exit(1);
			}
		}
		break;

		default: {
			std::cout << "Invalid solver specified." << std::endl;
			exit(1);
		}
		break;
	}

	if (show) 
	{
		double lastTime = omp_get_wtime();
		double elapsedSeconds = lastTime - startTime;
		std::cout << "Solver:" << Parameters::GetSolverName() << " GraphType:" << graph_type << " Scenario:" << UserGraph->scenario_no << " TotArcs:" << countArcs(LemonGraph) << " Odtrips:" << UserGraph->GetODtripCount() << " Nodes:" << countNodes(LemonGraph) << " Dfcit at source:" << UserGraph->deficits[0] << " NbSatisfiedTrips:" << numberOfSatisfiedTrips << " Optimal Objective Function value = " << std::fixed << std::setprecision(1) << objective_value << " Elapsed[s]:" << elapsedSeconds << std::endl;
	}	
	
    if(storeDuals)
	{
		int dual_obj = 0;
		if(graph_type == SingleSourceRecourse || graph_type == xZero)
			dual_obj = GetDualCostSingleSource(UserGraph, LemonGraph, capacities, costs, supplies, nodePotentials, targets, graph_type);
		else dual_obj = GetDualCost(UserGraph, LemonGraph, capacities, costs, supplies, nodePotentials, targets);
		//printf("Scenario:%d GraphType:%d Primal:%d Dual:%d\n",UserGraph->scenario_no, graph_type, (int)objective_value, dual_obj);
		
		/*int mcfpCost = 0;
		if(graph_type == xZero)
		{
			ListDigraph NewGraph;
			ListDigraph::ArcMap<int> capacities(NewGraph);
			ListDigraph::ArcMap<int> costs(NewGraph);
			ListDigraph::ArcMap<int> lemon_flows(NewGraph);
			ListDigraph::NodeMap<int> supplies(NewGraph);
			ListDigraph::NodeMap<int> nodePotentials(NewGraph);
			ListDigraph::ArcMap<int> lowerCapacities(NewGraph);
			
			MakeLemonGraph(UserGraph,NewGraph,capacities,lowerCapacities,costs,supplies,SingleSourceRecourse,bound_type,station_no,targets,hm,hp);
			NetworkSimplex<ListDigraph> NetSimplex(NewGraph);
			NetSimplex.upperMap(capacities);
			NetSimplex.costMap(costs);
			NetSimplex.supplyMap(supplies);
			
			NetSimplex.lowerMap(lowerCapacities);
			NetworkSimplex<ListDigraph>::ProblemType status = NetSimplex.run();
			
			int RyBar = NetSimplex.totalCost();
			int x0 = prob->GetQtot();
			mcfpCost = -1*objective_value - x0 * RyBar;
			printf("Mcfp e:%d Obj:%d R(yBar):%d x0:%d mcfpCost:%d dual:%d\n",UserGraph->scenario_no,(int)objective_value,RyBar,x0,mcfpCost,dual_obj);	
		}*/
		
		if ( std::abs(dual_obj - objective_value) > 0.1 )
		{
			printf("graph_type:%s\n",graph_type == SingleSourceRecourse ? "SingleSource" : "xZero");
			for(int i=0;i<prob->GetNodeCount();i++)
				printf("t%d:%d ",i,targets[i]);
			printf("\n");
			for(int i=0;i<prob->GetNodeCount();i++)
				printf("y0_%d:%d ",i,y0[i]);
			printf("\n");
			DebugInfGraph(UserGraph, LemonGraph, lowerCapacities, capacities, costs, supplies, graph_type, nodePotentials);
			printf("Wrong primal:%.1lf - dual:%.1lf values, exiting ...\n",objective_value, (double)dual_obj);
			exit(1);
		}
		
		UserGraph->SetDualObjectiveValue( (double)dual_obj );
		
		/*int nbNewDualVars = 0; int nbNewPotentials = 0;
		if(Parameters::GetModel() == 1 || Parameters::GetModel() == 3)
		{
			//printf("Magnanti86 Scenario:%d\n",UserGraph->scenario_no);
			UserGraph->ResetDuals();
			ListDigraph::NodeMap<int> newNodePotentials(LemonGraph);
			
			//Trying ... Violates optimality conditions!
			newNodePotentials[ LemonGraph.nodeFromId(0) ] = 0;
			newNodePotentials[ LemonGraph.nodeFromId(1) ] = -1*objective_value;
			
			//Nope! objective is wrong: c*(pi[Src] - pi[Snk]), c>1.
			//if(objective_value != 0)
			//{
				//newNodePotentials[ LemonGraph.nodeFromId(1) ] = -1*objective_value; //Think I need to change the source with sink ... to reflect the -1 constraint!!
				//newNodePotentials[ LemonGraph.nodeFromId(0) ] = 0;
			//} else {
				//newNodePotentials[ LemonGraph.nodeFromId(1) ] = nodePotentials[ LemonGraph.nodeFromId(1) ];
				//newNodePotentials[ LemonGraph.nodeFromId(0) ] = nodePotentials[ LemonGraph.nodeFromId(0) ];
			//}
			
			//for(int i=0;i<UserGraph->GetNodeCount();i++)
				//printf("nId:%d D:%d\n",i,UserGraph->dist[i]);
			
			ListDigraph::NodeMap<int> delta(LemonGraph);
			for(ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node)
				delta[node] = newNodePotentials[ LemonGraph.nodeFromId(1) ] - UserGraph->dist[LemonGraph.id(node)];
				//delta[node] = std::min( dual_obj - UserGraph->dist[LemonGraph.id(node)], newNodePotentials[ LemonGraph.nodeFromId(1) ] - UserGraph->dist[LemonGraph.id(node)]  );
				
			//Checking ...
			for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
			{
				ListDigraph::Node tail = LemonGraph.source(arc); ListDigraph::Node head = LemonGraph.target(arc);
				if(delta[head] > delta[tail] + costs[arc])
				{
					printf("delta[%d]:%d > delta[%d]:%d + cost:%d\n",LemonGraph.id(head),delta[head],LemonGraph.id(tail),delta[tail],costs[arc]) ; exit(1);
				}
						
			}				
			
			printf("Scenario:%d ObjectiveValue:%.1lf\n",UserGraph->scenario_no,objective_value);
			printf("NodeSource: oldPi:%d newPi:%d\n",nodePotentials[LemonGraph.nodeFromId(0)], newNodePotentials[LemonGraph.nodeFromId(0)]);
			printf("NodeSink: oldPi:%d newPi:%d\n",nodePotentials[LemonGraph.nodeFromId(1)], newNodePotentials[LemonGraph.nodeFromId(1)]);
			for(ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node)
			{
				if(LemonGraph.id( node ) == 0 || LemonGraph.id( node ) == 1) continue;
				newNodePotentials[node] = std::min( delta[node], nodePotentials[node]  );
				printf("Node:%d oldPi:%d newPi:%d delta:%d D:%d\n",LemonGraph.id( node ),nodePotentials[node],newNodePotentials[node],delta[node],UserGraph->dist[LemonGraph.id(node)]);
				//printf("Node:%d oldPi:%d delta:%d D:%d\n",LemonGraph.id( node ),nodePotentials[node],delta[node],UserGraph->dist[LemonGraph.id(node)]);
				if(newNodePotentials[node] != nodePotentials[node])
				{
					nbNewPotentials++;
					//printf("Node:%d oldPi:%d newPi:%d delta:%d D:%d\n",LemonGraph.id( node ),nodePotentials[node],newNodePotentials[node],delta[node],UserGraph->dist[LemonGraph.id(node)]);
				}
					
			} //getchar();
			
			int new_dual_obj = GetDualCostSingleSource(UserGraph, LemonGraph, capacities, costs, supplies, newNodePotentials, targets, graph_type);
			
			ListDigraph::ArcMap<int> reducedCosts(LemonGraph);
			ListDigraph::ArcMap<int> newReducedCosts(LemonGraph);
			
			//Non-depot arcs
			for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
			{
				ListDigraph::Node tail = LemonGraph.source(arc);
				ListDigraph::Node head = LemonGraph.target(arc);
				reducedCosts[arc] = costs[arc] + nodePotentials[tail] - nodePotentials[head];
				newReducedCosts[arc] = costs[arc] + newNodePotentials[tail] - newNodePotentials[head];
				
				if(LemonGraph.id(tail) == 0 && LemonGraph.id(head) != 1 && LemonGraph.id(tail) < prob->GetNodeCount()+2) continue; //Skip source to station arcs
				
				if(newReducedCosts[arc] < 0) //otherwise, dual will be zero and was already reset with ResetDuals()
				{
					int arcId = LemonGraph.id(arc);
					int dual = std::max( 0, -1*reducedCosts[arc]); //a positive number
					int newDual = std::max( 0, -1*newReducedCosts[arc]);
					if(newDual > dual)
					{
						printf("Critical Comparison:\n");
						printf("NonDepotArc Id:%d from:%d to:%d piFrom:%d piTo:%d rC:%d cap:%d dual:%d\n",arcId,LemonGraph.id(tail),LemonGraph.id(head),nodePotentials[tail],nodePotentials[head],reducedCosts[arc],capacities[arc],dual);
						printf("NewArc Id:%d from:%d to:%d piFrom:%d piTo:%d rC:%d cap:%d dual:%d\n",arcId,LemonGraph.id(tail),LemonGraph.id(head),newNodePotentials[tail],newNodePotentials[head],newReducedCosts[arc],capacities[arc],newDual);
						exit(1);
					}
					UserGraph->SetDual(arcId,newDual);
					if (dual < newDual) 
					{
						printf("Sce:%d newDual:%d oldDual:%d from:%d to:%d cap:%d cost:%d\n",UserGraph->scenario_no,newDual,dual,LemonGraph.id(tail),LemonGraph.id(head),capacities[arc],costs[arc]);
						nbNewDualVars++;
					}
					
				}
			}
			
			//The source - station arcs
			for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
			{
				ListDigraph::Node tail = LemonGraph.source(arc);
				ListDigraph::Node head = LemonGraph.target(arc);
				reducedCosts[arc] = costs[arc] + nodePotentials[tail] - nodePotentials[head];
				newReducedCosts[arc] = costs[arc] + newNodePotentials[tail] - newNodePotentials[head];
				
				if (!(LemonGraph.id(tail) == 0 && LemonGraph.id(head) != 1 && LemonGraph.id(head) < prob->GetNodeCount() + 2)) continue;
				
				int dual = 0; int newDual = 0;
				//if reducedCost[arc] < 0 ---> UpperCap constraint is active and LowerCap is not ---> Increases the objective
				//if reducedCost[arc] > 0 ---> LoweCap constraint is active and UpperCap is not ---> Reduced the objective
				dual = -1*reducedCosts[arc];
				newDual = -1*newReducedCosts[arc];
				
				int arcId = LemonGraph.id(arc);
				
				if(newDual > dual)
				{
					printf("DepotArc Id:%d from:%d to:%d piFrom:%d piTo:%d rC:%d cap:%d dual:%d\n",arcId,LemonGraph.id(tail),LemonGraph.id(head),nodePotentials[tail],nodePotentials[head],reducedCosts[arc],capacities[arc],dual);
					printf("NewArc Id:%d from:%d to:%d piFrom:%d piTo:%d rC:%d cap:%d dual:%d\n",arcId,LemonGraph.id(tail),LemonGraph.id(head),newNodePotentials[tail],newNodePotentials[head],newReducedCosts[arc],capacities[arc],newDual);
					exit(1);
				}
				if (dual < newDual) 
				{
					printf("Sce:%d newDual:%d oldDual:%d from:%d to:%d cap:%d cost:%d\n",UserGraph->scenario_no,newDual,dual,LemonGraph.id(tail),LemonGraph.id(head),capacities[arc],costs[arc]);
					nbNewDualVars++;
				}
					
				UserGraph->SetDual(arcId,newDual);
			}
			
			if ( std::abs(new_dual_obj - objective_value) > 1.0 )
			{
				DebugInfGraph(UserGraph, LemonGraph, capacities, costs, supplies, graph_type);
				printf("Wrong primal:%.1lf - newDual:%.1lf - oldDual:%.1lf values, exiting ...\n",objective_value, (double)new_dual_obj, (double) dual_obj);
				printf("OLD -> pi(Src):%d pi(Snk):%d\n",nodePotentials[LemonGraph.nodeFromId(0)],nodePotentials[LemonGraph.nodeFromId(1)]);
				printf("NEW -> pi(Src):%d pi(Snk):%d\n",newNodePotentials[LemonGraph.nodeFromId(0)],newNodePotentials[LemonGraph.nodeFromId(1)]);
				exit(1);
			}
			if(nbNewDualVars) printf("Magnanti86 Sce:%d nbNewDualVars:%d\n",UserGraph->scenario_no,nbNewDualVars);
			if(nbNewDualVars) printf("Magnanti86 Sce:%d nbNewPotentials:%d\n",UserGraph->scenario_no,nbNewPotentials);
		}*/  //End Magnating 
		
	}
	return objective_value;
}


void McfpSolvers::MakeLemonGraph(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & lowerCapacities,ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, int graph_type, int bound_type, int station_no, const std::vector<int> & targets, const std::vector<int>& hm, const std::vector<int>& hp)
{
	Prob * prob = UserGraph->GetProblem();
	int x0 = prob->GetQtot();
	//printf("In MakeLemonGraph graphType:%d sizeOf hm:%d hp:%d\n",graph_type,(int)hm.size(),(int)hp.size());
	
	std::vector<int> NewCap(prob->GetNodeCount(),0);
	
	//Nodes
	ListDigraph::Node lemonSource = LemonGraph.addNode();
	ListDigraph::Node lemonSink = LemonGraph.addNode();
	if(graph_type == SingleSourceRecourse) 
	{
		supplies[lemonSink] = -1*prob->GetQtot(); 
		supplies[lemonSource] = prob->GetQtot();
	}		
	if(graph_type == xZero)
	{
		supplies[lemonSink] = -( prob->GetQtot()+(prob->GetQtot()*x0) ); 
		supplies[lemonSource] = prob->GetQtot()+(prob->GetQtot()*x0); //piSource is positive ALWAYS.
		//printf("Sce:%d sourceDmd:%d sinkDmd:%d\n",UserGraph->scenario_no,supplies[lemonSource],supplies[lemonSink]);

	}
	if (graph_type == WS || graph_type == GlobalTargetLbUb || graph_type == TargetBounds || graph_type == WSnoCap)
	{
		supplies[lemonSource] = prob->GetQtot();
		supplies[lemonSink] = -1 * prob->GetQtot();
	}
	
	for (int i = 0; i < UserGraph->GetNodeCount()-2; i++)
		ListDigraph::Node lemonNode = LemonGraph.addNode();
	if (graph_type == WS || graph_type == GlobalTargetLbUb || graph_type == TargetBounds || graph_type == SingleSourceRecourse || graph_type == xZero)
	{
		for (int i = 2; i < UserGraph->GetNodeCount(); i++) {
			ListDigraph::Node lemonNode = LemonGraph.nodeFromId(i);
			supplies[lemonNode] = 0;
		}
		
	}
	
	int sum_deficits = 0;
	if (graph_type == Recourse || graph_type == OptStatCapRecourse)
	{
		for (int i = 0; i < targets.size(); i++)
			sum_deficits += targets[i];
		
		if(Parameters::GetModel() > 3 && sum_deficits > prob->GetQtot())
		{
			printf("Sum_deficits:%d Qtot:%d. The graph will be unfeasible\n",sum_deficits,prob->GetQtot()); 
			for(int i=0;i<targets.size();i++)
				printf("t%d:%d ",i,targets[i]);
			printf("\n");
			exit(1);
		}
		
		supplies[lemonSource] = prob->GetQtot() - sum_deficits;
		supplies[lemonSink] = -prob->GetQtot();
		
		
		for (int i = 0; i < targets.size(); i++) {
			ListDigraph::Node lemonNode = LemonGraph.nodeFromId(i + 2);
				supplies[lemonNode] = targets[i];
			//printf("tgtNode%d tgt:%d hp:%d hm:%d NewCap:%d\n",i,targets[i],hp[i],hm[i],targets[i] + hp[i] - hm[i]);
		}
	}		
		
	//Network Arcs
	{
		//Source Outgoing arcs
		for (int i = 0; i < prob->GetNodeCount(); i++) 
		{
			MCF_arc* arc = UserGraph->GetArc(i);
			ListDigraph::Arc lemonArc = LemonGraph.addArc(LemonGraph.nodeFromId(arc->from - 1), LemonGraph.nodeFromId(arc->to - 1));
			
			if (graph_type == Recourse || graph_type == OptStatCapRecourse) {
				capacities[lemonArc] = 0;
	} else if (graph_type == WSnoCap || (graph_type == TargetBounds || graph_type == GlobalTargetLbUb) && (Parameters::GetModel() == 6 || Parameters::GetModel() == 7)) {
				capacities[lemonArc] = 9999;
			} else {
				capacities[lemonArc] = prob->GetNode(i)->stationcapacity;
				if(prob->GetNode(i)->stationcapacity < prob->GetNode(i)->target)
				{
					printf("MCFP Something wrong ... Exiting.\n"); exit(1);
				}	
			}
			
			if (graph_type == TargetBounds && arc->to == station_no + 3) {
				costs[lemonArc] = (bound_type == 1) ? 1 : ((bound_type == 2) ? -1 : 1); //Lb=1; Ub=2
			}
			else
				costs[lemonArc] = 0;
			
			if(bound_type == 1 && arc->to == station_no + 3 && costs[lemonArc] != 1 || bound_type == 2 && arc->to == station_no + 3 && costs[lemonArc] != -1)
			{
				printf("Bound:%s Station:%d LemonArc: %d -> %d, Capacity: %d, Cost: %d\n",(bound_type==1)?"Lb":"Ub",station_no,LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),capacities[lemonArc],costs[lemonArc]); exit(1);
			}
			if(graph_type == SingleSourceRecourse)
			{
				capacities[lemonArc] = targets[i];
				lowerCapacities[lemonArc] = targets[i];
			}
			if(graph_type == xZero)
			{
				capacities[lemonArc] = y0[i] + targets[i]*x0;
				lowerCapacities[lemonArc] = y0[i] + targets[i]*x0;
				//printf("Sce: %d SrcArc: %d -> %d, lowerCapacity: %d, Capacity: %d, Cost: %d\n", UserGraph->scenario_no, LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),lowerCapacities[lemonArc],capacities[lemonArc],costs[lemonArc]);
			}
			
			//if(graph_type == 7)
				//printf("SrcStationArc: %d -> %d, Capacity: %d, Cost: %d\n", LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),capacities[lemonArc],costs[lemonArc]);
			
			//if(graph_type== TargetBounds && arc->to == station_no + 3)
			//{	
				//arc->Show();
				
				//printf("Bound:%s Station:%d LemonArc: %d -> %d, Capacity: %d, Cost: %d\n",(bound_type==1)?"Lb":"Ub",//station_no,LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),capacities[lemonArc],costs[lemonArc]); 
				//getchar();
			//}
		}
		//printf("Non-depot arcs:\n");
		for (int i = prob->GetNodeCount(); i < UserGraph->GetNetworkArcCount()-1; i++)
		{
			MCF_arc* arc = UserGraph->GetArc(i);
			//arc->Show();
			ListDigraph::Arc lemonArc = LemonGraph.addArc(LemonGraph.nodeFromId(arc->from - 1),LemonGraph.nodeFromId(arc->to - 1));
			
			costs[lemonArc] = 0;
			if( graph_type == OptStatCapRecourse )
			{
				int node_no = arc->cust_no;
				if(node_no < 0 || node_no > prob->GetNodeCount() || node_no > hp.size() || node_no > hm.size() )
				{
					arc->Show(); std::cerr << "Wrong arc in McfpSolvers. node_no:" << node_no << " hm.size():" << hm.size() << " hp.size():" << hp.size() << " Exiting." << std::endl; exit(1);
				}
				//printf("node_no:%d tgt:%d hm:%d hp:%d NewCap:%d\n",node_no,targets[node_no],hm[node_no],hp[node_no],NewCap[node_no]);
				capacities[lemonArc] = arc->cap + hp[ node_no ] - hm[ node_no ];
				if(arc->cap + hp[ node_no ] - hm[ node_no ] < targets[node_no])
				{
					printf("Wrong holding arc cap. Node:%d LemonId:%d OldCap:%d h+:%d h-:%d NewCap:%d tgt:%d\n",node_no,node_no+2,arc->cap,hp[node_no],hm[node_no],arc->cap + hp[ node_no ] - hm[ node_no ],targets[node_no]); 
					exit(1);
				}
			} else 
				if (graph_type == WSnoCap || 
				((graph_type == TargetBounds || graph_type == GlobalTargetLbUb) && 
				 (Parameters::GetModel() == 6 || Parameters::GetModel() == 7))) 
				 {
					capacities[lemonArc] = 9999;
				 } else {
					capacities[lemonArc] = arc->cap;
				}
			if(graph_type == xZero)
				capacities[lemonArc] = arc->cap + arc->cap*x0;
			//arc->Show();
			//printf("HoldingLemonArc: %d -> %d, Capacity: %d, Cost: %d\n", LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),capacities[lemonArc],costs[lemonArc]);	
		}
		
		//The Source to Sink arc.
		MCF_arc* arc = UserGraph->GetArc(UserGraph->GetNetworkArcCount()-1);
		//arc->Show(); getchar();
		if (arc->from != 1 || arc->to != 2) 
		{
			arc->Show();
			printf("Wrong arc in Lemon. Exiting\n"); exit(1);
		}

		ListDigraph::Arc lemonArc = LemonGraph.addArc(LemonGraph.nodeFromId(arc->from - 1), LemonGraph.nodeFromId(arc->to - 1));
		capacities[lemonArc] = prob->GetQtot();
		if(graph_type == xZero)
			capacities[lemonArc] = prob->GetQtot() + prob->GetQtot()*x0;
		
		if (graph_type == GlobalTargetLbUb)
		{
			costs[lemonArc] = (bound_type == 1) ? -1 : (bound_type == 2) ? 1 : 0; //Lb=1; Ub=2
			//printf("Bound:%d LemonArc: %d -> %d, Capacity: %d, Cost: %d\n", bound_type,LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),capacities[lemonArc],costs[lemonArc]); getchar();
		} else {
			costs[lemonArc] = 0;
		}
		//printf("Source2Sink LemonArc: %d -> %d, Capacity: %d, Cost: %d\n", LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),capacities[lemonArc],costs[lemonArc]); getchar();
	}
	//Trip Arcs
	for(int i=UserGraph->GetNetworkArcCount();i<UserGraph->GetArcCount();i++)
	{
		MCF_arc * arc = UserGraph->GetArc(i);
		if(arc->type != 3) 
		{
			arc->Show(); printf("Not a trip arc ... \n"); exit(1);
		}	
		ListDigraph::Arc lemonArc = LemonGraph.addArc(LemonGraph.nodeFromId(arc->from - 1), LemonGraph.nodeFromId(arc->to - 1));
		
		if(graph_type == xZero)
			capacities[lemonArc] = 1 + 1*x0;
		else capacities[lemonArc] = 1;
		
		if(graph_type == GlobalTargetLbUb || graph_type == TargetBounds) 
			costs[lemonArc] = -1000;
		else costs[lemonArc] = -1;
		//printf("TripLemonArc: %d -> %d, Capacity: %d, Cost: %d\n", LemonGraph.id(LemonGraph.source(lemonArc)), LemonGraph.id(LemonGraph.target(lemonArc)),capacities[lemonArc],costs[lemonArc]);
		
	}	
}

void McfpSolvers::StoreInSol(Sol * sol, int scenario)
{	
	printf("Storing scenario:%d in sol ...\n",scenario);
	
	for(int i=0; i<prob->GetNodeCount();i++)
		sol->SetTarget(scenario, i, flows[i]);
	sol->SetObjective(scenario, objective_value);	
}

void McfpSolvers::DebugInfGraph(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & lowerCapacities, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, int graph_type, ListDigraph::NodeMap<int> & nodePotentials)
{
	ScenarioGraph SceGraph(prob);
	std::vector<int> deficits(countNodes(LemonGraph));
	printf("In DbgGraph: sce:%d graph_type:%d\n",UserGraph->scenario_no,graph_type);
	for (ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc) 
	{			
		int rC = costs[arc] + nodePotentials[LemonGraph.source(arc)] - nodePotentials[LemonGraph.source(arc)];
		
		//printf("Arc: %d -> %d, lowerCap: %d, Capacity: %d, Cost:, %d rC: %d\n", LemonGraph.id(LemonGraph.source(arc)), LemonGraph.id(LemonGraph.target(arc)), lowerCapacities[arc], capacities[arc], costs[arc], rC);
		
		MCF_arc a;
		a.from = LemonGraph.id(LemonGraph.source(arc));
		a.to = LemonGraph.id(LemonGraph.target(arc));
		a.cap = capacities[arc];
		if (costs[arc] == -1 || costs[arc] == -1000) {
			a.type = 3;
		} else if (a.from == 0) {
			a.type = 1;
		} else {
			a.type = 2;
		}
		SceGraph.AddArc( a );
	}
	int sum_deficits = 0;
	// Print node supplies
	for (ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node) 
	{
		deficits[ LemonGraph.id(node) ] = supplies[node];
		//if(LemonGraph.id(node)>prob->GetNodeCount()+2 && graph_type != AllPositive) continue;
		
		sum_deficits += supplies[node];
	
		printf("Node %d, Supply: %d, pi: %d\n", LemonGraph.id(node), supplies[node], nodePotentials[node]);
	}
	printf("SumDficts Needs to be 0! sum_deficits:%d graphType:%d model:%d\n",sum_deficits,graph_type,Parameters::GetModel());
	
	printf("Printing the graph ...\n");
	SceGraph.PrintGraph((char*)"InfGraph.dot", deficits, countNodes(LemonGraph)); }

int McfpSolvers::GetDualCost(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, ListDigraph::NodeMap<int> & nodePotentials, const std::vector<int> & targets)
{
	int dual_obj = 0; int dual_obj_nodes = 0; int dual_obj_arcs = 0;
	UserGraph->ResetDuals();
	ListDigraph::ArcMap<int> reducedCosts(LemonGraph);
	
	for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
	{
		ListDigraph::Node tail = LemonGraph.source(arc);
		ListDigraph::Node head = LemonGraph.target(arc);
		reducedCosts[arc] = costs[arc] + nodePotentials[tail] - nodePotentials[head];
		
		if(reducedCosts[arc] < 0) //otherwise, dual will be zero and was already reset with ResetDuals()
		{
			int dual = std::max( 0, -1*reducedCosts[arc]);
			int arcId = LemonGraph.id(arc);
			UserGraph->SetDual(arcId,dual);
			dual_obj_arcs += dual * capacities[arc];
			//printf("NonZeroDualArc Id:%d from:%d to:%d piFrom:%d piTo:%d rC:%d cap:%d dual:%d\n",arcId,LemonGraph.id(tail),LemonGraph.id(head),nodePotentials[tail],nodePotentials[head],reducedCosts[arc],capacities[arc],(int)(UserGraph->GetArc(arcId)->dual));
		}
	}
	
	Prob * prob = UserGraph->GetProblem();
	
	UserGraph->pi.resize(prob->GetNodeCount() + 2, 0.0);
	for (ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node) 
	{
		int nodeId = LemonGraph.id(node);
		if (nodeId > prob->GetNodeCount() + 1) continue;

		//printf("Node:%d Supply:%d Pi:%d\n", nodeId, supplies[node], nodePotentials[node]);
		UserGraph->pi[nodeId] = nodePotentials[node];
		//printf("Node:%d Pi:%d\n", nodeId, UserGraph->pi[nodeId]);
	}
	
	ListDigraph::Node lemonSource = LemonGraph.nodeFromId(0); ListDigraph::Node lemonSink = LemonGraph.nodeFromId(1);
	//printf("SourcePi:%d SinkPi:%d\n",nodePotentials[lemonSource],nodePotentials[lemonSink]);
	
	for(int i=0;i<prob->GetNodeCount();i++)
			dual_obj_nodes += (UserGraph->pi[i+2] - UserGraph->pi[0]) * targets[i];
	
	dual_obj_nodes += prob->GetQtot() * ( UserGraph->pi[0]  - UserGraph->pi[1] ); 
	
	dual_obj = dual_obj_arcs + dual_obj_nodes;
	//printf("In Solver Sce:%d Primal:%.1lf Dual:%d\n",UserGraph->scenario_no, objective_value, dual_obj );
	
	return dual_obj;
}

int McfpSolvers::GetDualCostSingleSource(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, ListDigraph::NodeMap<int> & nodePotentials, const std::vector<int> & targets, int graph_type)
{
	int dual_obj = 0; int dual_obj_nodes = 0; int dual_obj_arcs = 0;
	Prob * prob = UserGraph->GetProblem();
	UserGraph->ResetDuals();
	ListDigraph::ArcMap<int> reducedCosts(LemonGraph);
	
	for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
	{
		ListDigraph::Node tail = LemonGraph.source(arc);
		ListDigraph::Node head = LemonGraph.target(arc);
		reducedCosts[arc] = costs[arc] + nodePotentials[tail] - nodePotentials[head];
		//To match Lemon where we have: c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
		
		if(LemonGraph.id(tail) == 0 && LemonGraph.id(head) != 1 && LemonGraph.id(tail) < prob->GetNodeCount()+2) continue; //Skip source to station arcs
		
		if(reducedCosts[arc] < 0) //otherwise, dual will be zero and was already reset with ResetDuals()
		{
			int dual = std::max( 0, -1*reducedCosts[arc]); //a positive number
			int arcId = LemonGraph.id(arc);
			UserGraph->SetDual(arcId,dual);
			dual_obj_arcs += dual * capacities[arc];
			//printf("NonDepotArc Sce:%d Id:%d from:%d to:%d piFrom:%d piTo:%d rC:%d cap:%d dual:%d\n",UserGraph->scenario_no,arcId,LemonGraph.id(tail),LemonGraph.id(head),nodePotentials[tail],nodePotentials[head],reducedCosts[arc],capacities[arc],(int)(UserGraph->GetArc(arcId)->dual));
		}
	}
	
	//The source - station arcs
	for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
	{
		ListDigraph::Node tail = LemonGraph.source(arc);
		ListDigraph::Node head = LemonGraph.target(arc);
		reducedCosts[arc] = costs[arc] + nodePotentials[tail] - nodePotentials[head];
		//To match Lemon where we have: c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
		
		if (!(LemonGraph.id(tail) == 0 && LemonGraph.id(head) != 1 && LemonGraph.id(head) < prob->GetNodeCount() + 2)) continue;
		
		int dual = 0;
		//if reducedCost[arc] < 0 ---> UpperCap constraint is active and LowerCap is not ---> Increases the objective
		//if reducedCost[arc] > 0 ---> LoweCap constraint is active and UpperCap is not ---> Reduced the objective
		dual = -1*reducedCosts[arc];
		int arcId = LemonGraph.id(arc);
		dual_obj_arcs += dual * capacities[arc];
		
		UserGraph->SetDual(arcId,dual);
		
		//if(dual!=0) 
			//printf("DepotArc Sce:%d Id:%d from:%d to:%d piFrom:%d piTo:%d rC:%d tgt:%d cap:%d dual:%d\n",UserGraph->scenario_no,arcId,LemonGraph.id(tail),LemonGraph.id(head),nodePotentials[tail],nodePotentials[head],reducedCosts[arc],targets[LemonGraph.id(head)-2],capacities[arc],(int)(UserGraph->GetArc(arcId)->dual));
	}
	
	graph_type == SingleSourceRecourse ? UserGraph->pi.resize(prob->GetNodeCount() + 2, 0.0) : UserGraph->pi.resize(countNodes(LemonGraph) + 2, 0.0);
	int trip_nodes = 0;
	for (ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node) 
	{
		int nodeId = LemonGraph.id(node);
		if (nodeId > prob->GetNodeCount() + 1 && (graph_type == SingleSourceRecourse || graph_type == xZero)) continue;

		//printf("Sce:%d Node:%d Supply:%d Pi:%d\n", UserGraph->scenario_no, nodeId, supplies[node], nodePotentials[node]);
		UserGraph->pi[nodeId] = nodePotentials[node];
		//printf("Node:%d Pi:%d\n", nodeId, UserGraph->pi[nodeId]);
		if(graph_type == AllPositive && nodeId >= 2)
		{
			trip_nodes += nodePotentials[node] * supplies[node];
			dual_obj_nodes += nodePotentials[node] * supplies[node];
		}
	}
	dual_obj_nodes += prob->GetQtot() * ( UserGraph->pi[0] - UserGraph->pi[1] ) ; 
	dual_obj = dual_obj_arcs + dual_obj_nodes;
	
	ListDigraph::Node lemonSource = LemonGraph.nodeFromId(0); ListDigraph::Node lemonSink = LemonGraph.nodeFromId(1);
	//printf("SourcePi:%d SinkPi:%d\n",nodePotentials[lemonSource],nodePotentials[lemonSink]);
	
	//printf("Targets: ");
	//for(int i=0;i<prob->GetNodeCount();i++)
		//printf("t%d:%d ",i,targets[i]);
	//printf("\n");
	//if(graph_type == SingleSourceRecourse)
		//printf("In DualCost Sce:%d Primal:%.1lf Dual:%d pi(Src):%d pi(Snk):%d dual_obj_nodes:%d dual_obj_arcs:%d\n",UserGraph->scenario_no, objective_value, dual_obj, UserGraph->pi[0], UserGraph->pi[1], dual_obj_nodes, dual_obj_arcs );
	//if(graph_type == AllPositive)
		//printf("In DualCost Sce:%d Primal:%.1lf Dual:%d pi(Src):%d pi(Snk):%d dual_obj_nodes:%d trip_nodes:%d dual_obj_arcs:%d\n",UserGraph->scenario_no, objective_value, dual_obj, UserGraph->pi[0], UserGraph->pi[1], dual_obj_nodes, trip_nodes, dual_obj_arcs );
	//if(graph_type == xZero)
		//printf("In DualCost Sce:%d Primal:%.1lf Dual:%d pi(Src):%d pi(Snk):%d dual_obj_nodes:%d dual_obj_arcs:%d\n",UserGraph->scenario_no, objective_value, dual_obj, UserGraph->pi[0], UserGraph->pi[1], dual_obj_nodes, dual_obj_arcs );
	
	return dual_obj;	
}

void McfpSolvers::PrintGraph(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & capacities, ListDigraph::ArcMap<int> & costs, ListDigraph::NodeMap<int> & supplies, char * graphName)
{
	printf("Mcfp Print Graph %s nodes:%d arcs:%d\n",graphName,countNodes(LemonGraph),countArcs(LemonGraph));
	std::vector<int> allSupplies(countNodes(LemonGraph),-1);
	for(ListDigraph::NodeIt node(LemonGraph); node != INVALID; ++node)
		allSupplies[ LemonGraph.id(node) ] = supplies[node];
	
	Network * n = new Network(allSupplies);
	int cntr = 0;
	for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
	{
		n->AddArc(LemonGraph.id( LemonGraph.source(arc) ), LemonGraph.id( LemonGraph.target(arc) ), cntr); cntr++;
	}
	cntr = 0;
	for(ListDigraph::ArcIt arc(LemonGraph); arc != INVALID; ++arc)
	{
		n->SetArc(cntr,1,capacities[arc],costs[arc]); cntr++;
	}
	n->PrintGraphViz(graphName);
	delete n;	
}

void McfpSolvers::CalculateTargetLevels(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & costs, NetworkSimplex<ListDigraph> & NetSimplex, int scenario)
{
	Prob * pr = UserGraph->GetProblem();	
	std::vector<int> pick_num(pr->GetNodeCount(),0);
	std::vector<int> pick_denom(pr->GetNodeCount(),0);
	std::vector<int> del_num(pr->GetNodeCount(),0);
	std::vector<int> del_denom(pr->GetNodeCount(),0);

	// USE THIS LOOP!
	//Trip Arcs
	for(int i=UserGraph->GetNetworkArcCount();i<UserGraph->GetArcCount();i++)
	{
		ListDigraph::Arc arc = LemonGraph.arcFromId(i);
		if(costs[arc] == -1) // Only trip arcs have a cost of -1
		{
			MCF_arc * user_arc = UserGraph->GetArc(i);
			int flow = NetSimplex.flow(arc);
			if(flow > 0)
			{
				pick_num[ user_arc->from_cust_no ]++;
				del_num[ user_arc->to_cust_no ]++;
			} 
			pick_denom[ user_arc->from_cust_no ]++;
			del_denom[ user_arc->to_cust_no ]++;
		}
	}
	
	std::vector<double> pick_tgt_level(pr->GetNodeCount(),-1);
	std::vector<double> del_tgt_level(pr->GetNodeCount(),-1);
	for(int i=0;i<pr->GetNodeCount();i++)
	{
		if(pick_denom[i] > 0)
			pick_tgt_level[i] = pick_num[i]/(double)pick_denom[i];
		
		if(del_denom[i] > 0)
			del_tgt_level[i] = del_num[i]/(double)del_denom[i];
	}
	
	printf("Storing tgt levels of scenario:%d\n",scenario);
	pr->StorePickTgtLvls(scenario,pick_tgt_level);
	pr->StoreDelTgtLvls(scenario,del_tgt_level);
}


//OG
/*void McfpSolvers::CalculateTargetLevels(ScenarioGraph * UserGraph, ListDigraph & LemonGraph, ListDigraph::ArcMap<int> & costs, NetworkSimplex<ListDigraph> & NetSimplex, int scenario)
{
	std::vector<int> pick_num(countNodes(LemonGraph),0);
	std::vector<int> pick_denom(countNodes(LemonGraph),0);
	std::vector<int> del_num(countNodes(LemonGraph),0);
	std::vector<int> del_denom(countNodes(LemonGraph),0);

	// USE THIS LOOP!
	//Trip Arcs
	for(int i=UserGraph->GetNetworkArcCount();i<UserGraph->GetArcCount();i++)
	{
		ListDigraph::Arc arc = LemonGraph.arcFromId(i);
		if(costs[arc] == -1) // Only trip arcs have a cost of -1
		{
			int flow = NetSimplex.flow(arc);
			if(flow > 0)
			{
				pick_num[ LemonGraph.id(LemonGraph.source(arc)) ]++;
				del_num[ LemonGraph.id(LemonGraph.target(arc)) ]++;
			} 
			
			pick_denom[ LemonGraph.id(LemonGraph.source(arc)) ]++;
			del_denom[ LemonGraph.id(LemonGraph.target(arc)) ]++;
		}
		
	}
	
	std::vector<double> pick_tgt_level(countNodes(LemonGraph),-1);
	std::vector<double> del_tgt_level(countNodes(LemonGraph),-1);
	
	for(int i=0;i<countNodes(LemonGraph);i++)
	{
		if(pick_denom[i] > 0)
			pick_tgt_level[i] = pick_num[i]/(double)pick_denom[i];
		
		if(del_denom[i] > 0)
			del_tgt_level[i] = del_num[i]/(double)del_denom[i];
	}
	
	Prob * prob = UserGraph->GetProblem();
	
	printf("Storing tgt levels of scenario:%d\n",scenario);
	prob->StorePickTgtLvls(scenario,pick_tgt_level);
	prob->StoreDelTgtLvls(scenario,del_tgt_level);
}*/