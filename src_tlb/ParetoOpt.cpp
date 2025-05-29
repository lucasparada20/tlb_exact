#include "ParetoOpt.h"

void ParetoOpt::Init(IloEnv env, ScenarioGraph * _UserGraph, Prob * _prob, std::vector<int> & _y0, std::vector<int> & _yBar)
{
	//printf("ParetoOpt Init sce:%d\n",_UserGraph->scenario_no);
	UserGraph = _UserGraph;
	prob = _prob;
	yBar = _yBar;
	y0 = _y0;
	
	model = IloModel(env);
	x0 = IloNumVar(env,0,IloInfinity,ILOFLOAT);
	x0.setName("x0");
	
	flow = IloNumVarArray(env,UserGraph->GetArcCount(),0,IloInfinity,ILOFLOAT);
	for(int i=0; i<UserGraph->GetArcCount(); i++) 
	{
		MCF_arc * a = UserGraph->GetArc(i);
		
		char name[40];
		if(a->type == 1)
			sprintf(name, "f%d", UserGraph->GetArc(i)->to-1);
		else if(a->type == 2)
			sprintf(name, "h%d_%d", UserGraph->GetArc(i)->from-1, UserGraph->GetArc(i)->to-1);
		else if(a->type == 3)
			sprintf(name, "t%d_%d", UserGraph->GetArc(i)->from-1, UserGraph->GetArc(i)->to-1);
		flow[i].setName(name);
		//printf("Sce:%d Setting name for flow[%d]: %s\n", UserGraph->scenario_no, i, name); 
		//UserGraph->GetArc(i)->Show();
	}

	
	flow_constrs = IloRangeArray(env);
	balance_constrs = IloRangeArray(env);
	
	for(int i=0;i<UserGraph->GetNodeCount();i++)
	{
		IloExpr expr(env);
		for(int k=0;k<UserGraph->GetInArcCount(i);k++)
			expr += flow[ UserGraph->GetInArc(i,k) ];
		for(int k=0;k<UserGraph->GetOutArcCount(i);k++)
			expr -= flow[ UserGraph->GetOutArc(i,k) ];
		
		if(i==0)
			expr += prob->GetQtot()+prob->GetQtot()*x0;
		else if(i==1)
			expr -= (prob->GetQtot()+prob->GetQtot()*x0);
		
		balance_constrs.add( expr == 0 );
	}
	for(int i=0;i<prob->GetNodeCount();i++)
		flow_constrs.add( flow[i] - yBar[i]*x0 == y0[i] );
	
	for (int i = prob->GetNodeCount(); i < UserGraph->GetArcCount(); i++)
	{
		MCF_arc * a = UserGraph->GetArc(i);		
		flow_constrs.add( flow[i] - a->cap*x0 <= a->cap );
	}
	
	model.add(balance_constrs);
	model.add(flow_constrs);
}

void ParetoOpt::Solve(IloEnv env)
{
	IloExpr obj1(env);
	for(int i=0;i<UserGraph->GetArcCount();i++)
		if( UserGraph->GetArc(i)->type == 3 )
			obj1 -= flow[i];
	obj1 -= RyBar*x0;
	obj_func = IloMinimize(env, obj1);
	model.add(obj_func);
	obj1.end();	
	
	cplex = IloCplex(model);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
	cplex.setParam(IloCplex::Param::TimeLimit, 999999);
	cplex.setOut(env.getNullStream());	
	
	std::string fileName;
	//fileName = "ParetoOptSce" + std::to_string(UserGraph->scenario_no) + ".lp";
	//cplex.exportModel(fileName.c_str());

	clock_t start_time_clock = clock();
	re = cplex.solve();
	cplex_status = (int)cplex.getCplexStatus();
	sol_value = re?cplex.getObjValue():0; //This the Ub from cplex
	cplex_nb_nodes = (int)cplex.getNnodes();
	clock_t end_time_clock = clock();

	double time_taken = (double) ( end_time_clock - start_time_clock )/CLOCKS_PER_SEC;
	
	printf("ParetoOpt e:%d re:%d Obj:%.2lf R(yBar):%.2lf x0:%.2lf Arcs:%d Nodes:%d Trips:%d Time:%.2lf\n",UserGraph->scenario_no,re,sol_value,RyBar,(double)cplex.getValue(x0),UserGraph->GetArcCount(),UserGraph->GetNodeCount(),UserGraph->GetODtripCount(),time_taken);
	
	//Store dual vars ...
    if (re && cplex_status == 1) {
        
		IloNumArray duals_balance(env);
        IloNumArray duals_flow(env);

        cplex.getDuals(duals_balance, balance_constrs);
        cplex.getDuals(duals_flow, flow_constrs);
		
		UserGraph->ResetDuals();
		
		UserGraph->pi.resize(2, 0.0);
		UserGraph->pi[0] = duals_balance[0];
		UserGraph->pi[1] = duals_balance[1];
		printf("Sce:%d piSource:%d piSink:%d\n",UserGraph->scenario_no,UserGraph->pi[0],UserGraph->pi[1]);
		for(int i=0;i<UserGraph->GetArcCount();i++)
		{
			MCF_arc * a = UserGraph->GetArc(i);
			UserGraph->SetDual(i,-duals_flow[i]);
			printf("%s flow:%.1lf dual:%.1lf\n",flow[i].getName(),cplex.getValue(flow[i]),-duals_flow[i]);
		}
		
        duals_balance.end();
        duals_flow.end();
    } else if (cplex_status != 1)
	{
		printf("Wrong ParetoOpt LP, exiting ...\n"); 
		std::string fileName;
		fileName = "GraphSce" + std::to_string(UserGraph->scenario_no) + ".dot";
		UserGraph->PrintGraph((char*)fileName.c_str(), yBar, (int)balance_constrs.getSize());
		exit(1);
	}		
	
}

void ParetoOpt::SolveMCF(IloEnv env)
{
	IloExpr obj1(env);
	for(int i=0;i<UserGraph->GetArcCount();i++)
		if( UserGraph->GetArc(i)->type == 3 )
			obj1 -= flow[i];
	obj_func = IloMinimize(env, obj1);
	model.add(obj_func);
	obj1.end();	
	
	x0.setUB(prob->GetQtot());
	x0.setLB(prob->GetQtot());
	
	cplex = IloCplex(model);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
	cplex.setParam(IloCplex::Param::TimeLimit, 999999);
	cplex.setOut(env.getNullStream());	

	std::string fileName;
	//fileName = "ParetoOptSce" + std::to_string(UserGraph->scenario_no) + ".lp";
	//cplex.exportModel(fileName.c_str());

	clock_t start_time_clock = clock();
	

	cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Network);
	cplex.setParam(IloCplex::Param::NodeAlgorithm, IloCplex::Network); 	
 	IloInt rootAlg = cplex.getParam(IloCplex::Param::RootAlgorithm);
    //std::cout << "Root node algorithm parameter: " << rootAlg << std::endl; getchar();
	
	re = cplex.solve();
	cplex_status = (int)cplex.getCplexStatus();
	sol_value = re?cplex.getObjValue():0; //This the Ub from cplex
	double mcfpCost = sol_value - RyBar*cplex.getValue(x0);
	cplex_nb_nodes = (int)cplex.getNnodes();
	clock_t end_time_clock = clock();

	double time_taken = (double) ( end_time_clock - start_time_clock )/CLOCKS_PER_SEC;
	
	//Store dual vars ...
	int dual_obj = 0; int dual_obj_arcs = 0; int dual_obj_nodes = 0;
    if (re && cplex_status == 1) {
        
		IloNumArray duals_balance(env);
        IloNumArray duals_flow(env);

        cplex.getDuals(duals_balance, balance_constrs);
        cplex.getDuals(duals_flow, flow_constrs);
		
		/*for(int i=0;i<prob->GetNodeCount();i++)
			printf("t%d:%d ",i,yBar[i]);
		printf("\n");
		for(int i=0;i<prob->GetNodeCount();i++)
			printf("y0_%d:%d ",i,y0[i]);
		printf("\n");*/
		
		UserGraph->ResetDuals();
		
		UserGraph->pi.resize(2, 0.0);
		UserGraph->pi[0] = duals_balance[0];
		UserGraph->pi[1] = duals_balance[1];
		//printf("Sce:%d piSource:%d piSink:%d\n",UserGraph->scenario_no,UserGraph->pi[0],UserGraph->pi[1]);
		dual_obj_nodes = prob->GetQtot() * ( UserGraph->pi[0] - UserGraph->pi[1] );
		
		//Building the dual objective ..
		
		std::vector<int> cap(UserGraph->GetArcCount(),0);
		for(int i=0;i<prob->GetNodeCount();i++)
			cap[i] = y0[i] + yBar[i]*prob->GetQtot(); 
		for (int i = prob->GetNodeCount(); i < UserGraph->GetArcCount(); i++)
		{
			MCF_arc * a = UserGraph->GetArc(i);
			cap[i] = a->cap*prob->GetQtot() + a->cap;
		}
		
		for(int i=0;i<UserGraph->GetArcCount();i++)
		{
			MCF_arc * a = UserGraph->GetArc(i);
			UserGraph->SetDual(i,-duals_flow[i]);
			//if(duals_flow[i] != 0)
				//printf("%s flow:%.1lf dual:%.1lf cap:%d\n",flow[i].getName(),cplex.getValue(flow[i]),-duals_flow[i],cap[i]);
			dual_obj_arcs += duals_flow[i] * cap[i];
		}

        duals_balance.end();
        duals_flow.end();
		
		dual_obj += dual_obj_arcs + dual_obj_nodes;
		if(std::abs(dual_obj - sol_value) > 0.01 )
		{
			printf("Wrong dual in sce:%d dual:%d dual_obj_nodes:%d dual_obj_arcs:%d mcfpCost:%.2lf\n", UserGraph->scenario_no,dual_obj,dual_obj_nodes,dual_obj_arcs,sol_value); exit(1);
		}
			
		//getchar();
		
    } else if (cplex_status != 1)
	{
		printf("Wrong ParetoOpt LP, exiting ...\n"); 
		std::string fileName;
		fileName = "GraphSce" + std::to_string(UserGraph->scenario_no) + ".dot";
		UserGraph->PrintGraph((char*)fileName.c_str(), yBar, (int)balance_constrs.getSize());
		exit(1);
	}
	
	//printf("ParetoOpt e:%d re:%d Obj:%.2lf R(yBar):%.2lf x0:%.2lf McfpCost:%.2lf dual:%d dualNodes:%d dualArcs:%d Arcs:%d Nodes:%d Trips:%d Time:%.2lf\n",UserGraph->scenario_no,re,sol_value,RyBar,(double)cplex.getValue(x0),mcfpCost,dual_obj,dual_obj_nodes,dual_obj_arcs,UserGraph->GetArcCount(),UserGraph->GetNodeCount(),UserGraph->GetODtripCount(),time_taken);
	
	//RyBar doesn't need to be the same as mcfpCost ...
	//if(std::abs(RyBar - mcfpCost) > 0.01)
	//{
		//printf("ParetoOpt sce:%d Wrong RyBar:%.2lf McfpCost:%.2lf\n",UserGraph->scenario_no,RyBar,mcfpCost); //exit(1);
	//}
}