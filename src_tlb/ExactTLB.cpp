#include "ExactTLB.h"
#include "DP.h"
#include <iostream>
#include <cstring>
#include <numeric> // std::accumulate
#include <omp.h>

void ExactTlb::Solve(Sol * _sol, RecourseCalculation * _r)
{
	prob = _sol->GetProblem();
	sol =_sol;
	scs = _r->GetScenarios();
	r = _r;
	
	double total_startTime = omp_get_wtime();	

	nb_vars = (Parameters::GetModel() == 2 || Parameters::GetModel() == 6) ? scs->GetScenarioCount() :
				  (Parameters::GetModel() == 3 || Parameters::GetModel() == 7) ? prob->GetNodeCount() : 0;
	
	IloEnv env;
	InitMaster(env);
	if(BandCheck) SolveBranchAndCheck(env);
	else SolveProblem(env);
	
	double total_lastTime = omp_get_wtime();
	double elapsedSeconds = total_lastTime - total_startTime;
	
	//Store relevant data in the solution object
	sol->SetRe(cplex_status == 1 ? 1 : 0);
	sol->SetNbBendersCuts(nb_benders_cuts);
	sol->SetCplexNodes(cplex_nb_nodes);
	sol->SetExactUb(best_upper_bound);
	sol->SetExactLb(sol_value);
	sol->SetExactTime(elapsedSeconds);

	if (Parameters::GetModel() == 2 || Parameters::GetModel() == 3 || Parameters::GetModel() == 6 || Parameters::GetModel() == 7) {
		sol->SetWOmegaVars(nb_vars);
		sol->SetZOmegaVars(nb_vars);

		for (int k = 0; k < nb_vars; k++) {
			sol->SetW((int)(cplex.getValue(w[k]) + 0.1));
			sol->SetZ((int)(cplex.getValue(z[k]) + 0.1));

			if (Parameters::GetModel() == 2) {
				for (int i = 0; i < prob->GetNodeCount(); i++) {
					sol->SetWOmega(k, (int)(cplex.getValue(w_array[i][k]) + 0.1));
					sol->SetZOmega(k, (int)(cplex.getValue(z_array[i][k]) + 0.1));
				}
			}
		}
	}

	//station #, lb_i1, ub_i1, lb_i2, ub_i2,..., z_i1, w_i1, ..., z_i4, w_i4
	// Open a file for writing
	/*FILE *file = fopen("bounds.txt", "w");
	if (file != NULL) {
		fprintf(file, "i=0,...,%d:\n", prob->GetNodeCount() - 1);
		fprintf(file, "e=0,...,%d:\n", scs->GetScenarioCount() - 1);
		fprintf(file, "#;OptTgt;");

		for (int e = 0; e < scs->GetScenarioCount(); e++) {
			fprintf(file, "lb_i%d;ub_i%d;", e, e);
		}

		for (int e = 0; e < scs->GetScenarioCount(); e++) {
			fprintf(file, "z_i%d;w_i%d;", e, e);
		}

		fprintf(file, "\n");

		for (int i = 0; i < prob->GetNodeCount(); i++) {
			fprintf(file, "%d;%d;", i,(int)(cplex.getValue(y[i])+0.1));

			for (int e = 0; e < Parameters::GetScenarioCount(); e++) {
				fprintf(file, "%d;%d;", sol->GetTargetLb(e, i), sol->GetTargetUb(e, i));
			}

			for (int e = 0; e < Parameters::GetScenarioCount(); e++) {
				fprintf(file, "%d;%d;", (int)(cplex.getValue(z_array[i][e]) + 0.1), (int)(cplex.getValue(w_array[i][e]) + 0.1));
			}

			fprintf(file, "\n");
		}

		// Close the file when done writing
		fclose(file);
	} else {
		printf("Unable to create the file.\n");
	}
	exit(1);*/	
	
	Clear();
	env.end();	
}

void ExactTlb::InitMaster(IloEnv env)
{
	printf("ExactClass InitMaster");
	model = IloModel(env);
	int Qtot=0;
	if(prob->GetUpperBound() < 0.0 )
	{
		printf("Wrong UB in exact class:%.1lf\n. Exiting ...\n",prob->GetUpperBound()); exit(1);
	}
	theta = IloNumVar(env,0,prob->GetUpperBound()+0.1,ILOFLOAT); theta.setName("t");
	y = IloNumVarArray(env,prob->GetNodeCount(),0,IloInfinity,ILOINT);
	for(int i=0; i<prob->GetNodeCount(); i++)
	{
		Node * n = prob->GetNode(i);
		//n->Show();
		char name[40];
		sprintf(name,"y%d",n->no);
		y[i].setName(name);
		
		//if(n->stationcapacity < n->max_scenario_out_trips)
			y[i].setUB(n->stationcapacity);
		//else y[i].setUB(n->max_scenario_out_trips);
	
		Qtot += n->stationcapacity;
	}
	printf("Model:%d nbCust:%d SystemAvailaleBikes:%d SumStationCap:%d\n",Parameters::GetModel(),prob->GetNodeCount(),prob->GetQtot(),Qtot);
	
	//Constraints for all 3 Models
	{
		IloExpr expr(env);
		for(int i=0;i<prob->GetNodeCount();i++)
			expr += y[i];
		if(Parameters::GetDelta() == 0) model.add(expr<=prob->GetQtot());
		else {
			model.add(expr <= prob->GetQtot());
			model.add(expr >= prob->GetQtot() - Parameters::GetDelta());
			
		}
		expr.end();
		
		//IloExpr expr1(env);
		//for(int i=0;i<prob->GetNodeCount();i++)
			//expr1 += y[i];
		//model.add(expr1>=1);
		//expr1.end();		
	} 
	
	/*IloExpr obj1(env);
	obj1 += theta;
	obj_func = IloMaximize(env, obj1);
	model.add(obj_func);
	obj1.end();*/
	
	z = IloNumVarArray(env,nb_vars,0,IloInfinity,ILOFLOAT);
	w = IloNumVarArray(env,nb_vars,0,IloInfinity,ILOFLOAT);	
	
	for(int k=0;k<nb_vars;k++)
	{
		char name_z[40]; char name_w[40];
		sprintf(name_z,"z%d",k); sprintf(name_w,"w%d",k);
		//printf("Naming z%d w%d\n",k,k);
		z[k].setName(name_z); w[k].setName(name_w);		
	}

	if(Parameters::GetModel() == 2 || Parameters::GetModel() == 6)
	{
		z_array = IloArray<IloNumVarArray>(env,prob->GetNodeCount());
		w_array = IloArray<IloNumVarArray>(env,prob->GetNodeCount());	
		for(int i=0; i<prob->GetNodeCount(); i++)
		{	
			z_array[i] = IloNumVarArray(env,scs->GetScenarioCount(),0,IloInfinity,ILOFLOAT);
			w_array[i] = IloNumVarArray(env,scs->GetScenarioCount(),0,IloInfinity,ILOFLOAT);
			
			for(int e=0; e<scs->GetScenarioCount(); e++)
			{	
				char name_z_array[40]; char name_w_array[40];
				sprintf(name_z_array,"z%d_%d",i+1,e); sprintf(name_w_array,"w%d_%d",i+1,e);
				//printf("Naming z%d_%d w%d_%d\n",i+1,e,i+1,e);
				z_array[i][e].setName(name_z_array); w_array[i][e].setName(name_w_array);
			}
				
		}		
	}
	
	//Constraints for Models 2
	if (Parameters::GetModel() == 2 || Parameters::GetModel() == 6) 
	{
		//In model 2: Eq. (19)-(20)
		IloExpr expr(env);
		for (int i = 0; i < prob->GetNodeCount(); i++)
			expr += y[i];
		for (int k = 0; k < scs->GetScenarioCount(); k++) 
		{			
			model.add(expr + z[k] >= scs->GetBikeLb(k));
			model.add(expr - w[k] <= scs->GetBikeUb(k));
		}
		expr.end();
		
		//Eq. (21)
		IloExpr expr1(env);
		for(int e=0;e<scs->GetScenarioCount();e++)
			expr1 += (1.0/(double)scs->GetScenarioCount())*(z[e]+w[e]);
		model.add(theta - (prob->GetUpperBound()+0.01) + expr1 <= 0);
		expr1.end();
		
		//Eq. (22)-(23)
		for(int i=0;i<prob->GetNodeCount();i++)
		{
			IloExpr expr2(env);
			expr2 += y[i];
			for(int e=0;e<scs->GetScenarioCount();e++)
			{
				model.add( expr2 + z_array[i][e] >= scs->GetTargetLb(e,i) );
				model.add( expr2 - w_array[i][e] <= scs->GetTargetUb(e,i) );
			}
			expr2.end();
		}
		
		//Eq. (24)
		IloExpr expr3(env);
		for(int i=0;i<prob->GetNodeCount();i++)
			for(int e=0;e<scs->GetScenarioCount();e++)
				expr3 += (1.0/(double)scs->GetScenarioCount())*(z_array[i][e]+w_array[i][e]);
		model.add(theta - (prob->GetUpperBound()+0.01) + expr3 <= 0);
		expr3.end();		
	}
	
	//Constraints for Models 3
	if(Parameters::GetModel()==3 || Parameters::GetModel() == 7)
	{
		//Eq. (28)-(29)
		for(int i=0;i<prob->GetNodeCount();i++)
		{	
			model.add(y[i] + z[i] >= prob->GetStationLb(i));
			model.add(y[i] - w[i] <= prob->GetStationUb(i));
		}		
		//Eq. (30)
		IloExpr expr(env);
		expr += theta - (prob->GetUpperBound()+0.01);
		for(int i=0;i<prob->GetNodeCount();i++)
			expr += ( z[i] + w[i] );
		model.add( expr <= 0 );
		expr.end();
	}

	//Constraints for Models 5
	if(Parameters::GetModel()==5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
	{
		//theta.setUB( IloInfinity ); //WS with no cap is calculated!
		//hm = IloNumVarArray(env,prob->GetNodeCount(),0,IloInfinity,ILOINT);
		hm = IloNumVarArray(env,prob->GetNodeCount(),0,0,ILOINT);
		hp = IloNumVarArray(env,prob->GetNodeCount(),0,IloInfinity,ILOINT);
		IloExpr exprSumH(env);
		for(int i=0;i<prob->GetNodeCount();i++)
		{
			Node * n = prob->GetNode(i);

			char name_hp[40]; char name_hm[40];
			sprintf(name_hp,"hp%d",i); sprintf(name_hm,"hm%d",i);
			hm[i].setName(name_hm); hp[i].setName(name_hp);	
			
			exprSumH += hp[i] + hm[i];
			
			y[i].setUB( IloInfinity );
			
			IloExpr expr(env);
			expr += y[i] - hp[i] + hm[i];
			model.add( expr <= n->stationcapacity );
			model.add( 0 <= expr );
			expr.end();
		}
		printf("Budget:%.2lf\n",Parameters::GetBudget());
		if(Parameters::GetBudget() < 0.01)
		{
			printf("Budget is zero. Forgot to give budget?\n"); exit(1);
		}
		model.add( exprSumH <= std::floor( Parameters::GetBudget() * prob->GetCapTot() ) );
		exprSumH.end();
	}
	
	/*if(Parameters::GetModel()==1 || Parameters::GetModel()==3)
	{
		std::vector<int> WStargets(prob->GetNodeCount(),0);
		//Store the average targets
		for(int i=0;i<prob->GetNodeCount();i++)
		{
			int target=0;
			for(int e=0;e<scs->GetScenarioCount();e++)
				target += sol->GetTarget(e,i);
			target = std::ceil( target / (double)scs->GetScenarioCount() );	
			WStargets[i] = target <= prob->GetNode(i)->stationcapacity ? target : prob->GetNode(i)->stationcapacity;
		}
		double WStargetsCost = r->Calculate(WStargets);
		// Vector to store pairs of (event time, index)
		std::vector<std::pair<int, int>> event_times(prob->GetNodeCount());

		for(int i = 0; i < prob->GetNodeCount(); i++) {
			int sum_event_time = 0;
			for(int e = 0; e < scs->GetScenarioCount(); e++) {
				ScenarioGraph * UserGraph = r->GetUserGraph(e);
				sum_event_time += UserGraph->GetNbEventTimes(i);
			}
			event_times[i] = {sum_event_time, i}; 
		}
		// Sort in descending order by sum_event_time (the first element of the pair)
		std::sort(event_times.begin(), event_times.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
			return a.first > b.first;  // Compare by sum_event_time in descending order
		});
		
		std::cout << "Sorted by nb event times:\n";
		for(int i = 0; i < prob->GetNodeCount(); i++) {
			std::cout << "Index: " << event_times[i].second 
					  << ", TotalEvents: " << event_times[i].first << ", AvgEvents: " << std::ceil(event_times[i].first/(double)scs->GetScenarioCount())
					  << ", Cap: " << prob->GetNode(i)->stationcapacity
					  << ", WStarget: " << WStargets[event_times[i].second] << std::endl;
		}
		
		
		for (int k = 0; k < std::min(5, prob->GetNodeCount()); k++) {
			int no = event_times[k].second;
			std::vector<int> v1 = WStargets;
			std::vector<int> v2 = WStargets;
			
			if (WStargets[no] > 0 && WStargets[no] < prob->GetNode(no)->stationcapacity) {
				 v1[no] = std::min(WStargets[no] + 1, prob->GetNode(no)->stationcapacity); 
				v2[no] = std::max(WStargets[no] - 1, 0);  
			}

			for (int cut = 0; cut < 2; ++cut) {
				std::vector<int> v = (cut == 0) ? v1 : v2;  // Select v1 for the first iteration, v2 for the second
				std::vector<double> duals(prob->GetNodeCount() + 1, 0.0);
				std::vector<double> capacity_duals;
				std::vector<int> hp_int_vec;
				std::vector<int> hm_int_vec;
				
				double RecCost = r->Calculate(v, duals, capacity_duals, hm_int_vec, hp_int_vec);
				if (RecCost > 0.1 ) {
					IloExpr expr(env);
					double rhs = 0.0;
					
					for (int i = 0; i < prob->GetNodeCount(); i++) {
						expr -= duals[i] * y[i];
						rhs += duals[i] * v[i];
					}
					
					expr -= duals[prob->GetNodeCount()];
					expr += theta;
					rhs += duals[prob->GetNodeCount()];
					
					model.add(expr <= 0);
					
					const char* cut_type = (cut == 0) ? "v1" : "v2";
					//printf("%s targets: ",cut_type);
					//for(int i=0;i<v.size();i++)
						//printf("t%d:%d ",i,v[i]);
					//printf("\n");
					printf("BendersCut for %s with varying station %d at the root with rhs: %.2lf\n", cut_type, no, rhs);
					std::cout << expr << " <= 0" << std::endl;
					
					expr.end();
				}
			}
		}
	}*/

	IloExpr obj1(env);
	obj1 += theta;	
	obj_func = IloMaximize(env, obj1);
	model.add(obj_func);
	obj1.end();
		
	cplex = IloCplex(model);
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
	cplex.setParam(IloCplex::Param::TimeLimit, max_time);
	
	//Model 1 is unbounded unless Optimality Cuts are added
	if(!BandCheck)
	{
		int model = Parameters::GetModel();
		if (model == 5 || model == 6 || model == 7) {
			lazy_call = new (env) ExactTlbLazyCallBack(env, prob, r, y, hm, hp, theta);
		} else if (model <= 3) {
			lazy_call = new (env) ExactTlbLazyCallBack(env, prob, r, y, theta);
		} 
		cplex.use(lazy_call);
	}
		
}

void ExactTlb::SolveBranchAndCheck(IloEnv env)
{
	best_sol = 0.0; best_gap = 100.0;
	int itr_cntr = 0; clock_t start_time_clock = clock();
	cplex.setOut(env.getNullStream());
	std::vector<int> targets;
	std::vector<int> hp_int_vec; std::vector<int> hm_int_vec;
	do
	{
		re = cplex.solve();
		sol_value = re?cplex.getObjValue():0;
		cplex_status = (int)cplex.getCplexStatus();
		clock_t end_time_clock = clock();
		double elapsedSeconds = (double)(end_time_clock - start_time_clock) / CLOCKS_PER_SEC;
		
		IloNumArray values(env);
		cplex.getValues(values,y);
		targets.clear();
		for(int i=0; i<prob->GetNodeCount(); i++)
		{
			targets.push_back((int)(values[i]+0.1));
			sol->SetTarget(i,(int)(values[i]+0.1));
			
			if(Parameters::GetModel()==5)
			{
				hm_int_vec.push_back( (int)cplex.getValue( hm[i] )+0.1 ); hp_int_vec.push_back( (int)cplex.getValue( hp[i] )+0.1 );
			}
		}		
		itr_cntr++;
		double RecCost = r->Calculate(targets,hm_int_vec,hp_int_vec); //computes and stores the unsatisfied requests
		double gap = 100.0*((double)cplex.getValue(theta) - RecCost) /  RecCost;
		//cplex.exportModel("B&Check.lp");		
		
		if(RecCost > best_sol)
		{
			best_sol = RecCost;
			best_gap = gap;
			best_solution.clear();
			printf("new_sol cost:%.1lf\n",RecCost);
			for(int i = 0; i < sol->GetNodeCount(); i++)
				best_solution.push_back( targets[i] );
		}			
		if(itr_cntr % 10 == 0)
			printf("B&Check iter:%d BestRec:%.3lf Obj:%.3lf Gap%%:%.3lf ElapsedT:%.1lf\n",itr_cntr,best_sol,sol_value,best_gap,elapsedSeconds);	

	if(Parameters::GetModel()==5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
	{
		IloNumArray valueshp(env); IloNumArray valueshm(env);
		cplex.getValues(valueshp,hp); cplex.getValues(valueshm,hm);
		for(int i=0;i<hp.getSize();i++)
		{
			hp_int_vec.push_back( (int)valueshp[i] );
			hm_int_vec.push_back( (int)valueshm[i] );
		}
		valueshm.end(); valueshp.end();
	}			
		
	}while(AddOptCut(env,targets,hm_int_vec,hp_int_vec));

	double recCost=0.0;
	clock_t end_time_clock = clock();
	double time_taken = (double)(end_time_clock - start_time_clock) / CLOCKS_PER_SEC;
	
	std::vector<int> exact_start_targets;
	for(int i=0; i<prob->GetNodeCount(); i++)
	{
		sol->SetExactTarget((int)(cplex.getValue(y[i])+0.1));
		exact_start_targets.push_back((int)cplex.getValue(y[i]+0.1));
		if(Parameters::GetModel()==5)
		{
			hm_int_vec.push_back( (int)cplex.getValue( hm[i] )+0.1 ); hp_int_vec.push_back( (int)cplex.getValue( hp[i] )+0.1 );
		}		
	}	
	recCost = r->Calculate(exact_start_targets,hm_int_vec,hp_int_vec); //computes and stores the unsatisfied requests
	sol->SetUnsatisfiedRequests(r->GetUnsatisfiedReqVec());

	printf("McfpSolver:%s B&Check Model:%d re:%d status:%d CplexSol:%.3lf Theta:%.3lf WS:%.2lf RecCost:%.3lf unsatisfiedReq:%.1lf BendersCuts:%d time:%.3lf\n",Parameters::GetSolverName(), Parameters::GetModel(),(int)re, cplex_status, sol_value, (double)cplex.getValue(theta),prob->GetUpperBound(), recCost, sol->GetExactUnsatisfiedReq(), nb_benders_cuts, time_taken);
	printf("Target levels:\n");
	for(int i=0; i<sol->GetNodeCount(); i++)
		printf("y[%d]:%d ",i,(int)cplex.getValue(y[i]));
	printf("\n");
	
}

void ExactTlb::SolveProblem(IloEnv env)
{
    int numVariables = cplex.getNcols();
    int numConstraints = cplex.getNrows();
    std::cout << "Model: " << Parameters::GetModel() << " Number of Variables before Cplex prepro: " << numVariables << std::endl;
    std::cout << "Model: " << Parameters::GetModel() << " Number of Constraints before Cplex prepro: " << numConstraints << std::endl;

	try
	{	
		//if(Parameters::GetModel() == 1 || Parameters::GetModel() == 11)
			//SetMIPstartDP(env);		

		std::string fileName;
		int model = Parameters::GetModel();
		fileName = "sbrpod_model" + std::to_string(model) + ".lp";
		cplex.exportModel(fileName.c_str());
		
		clock_t start_time_clock = clock();
		re = cplex.solve();
		cplex_status = (int)cplex.getCplexStatus();
		sol_value = re?cplex.getObjValue():0; //This the Lb from cplex
		best_upper_bound = re?cplex.getBestObjValue():prob->GetUpperBound(); //This the Ub from cplex
		cplex_nb_nodes = (int)cplex.getNnodes();
		clock_t end_time_clock = clock();
		
		time_taken = (double) ( end_time_clock - start_time_clock )/CLOCKS_PER_SEC;
		
		if(lazy_call)
		{
			nb_benders_cuts = lazy_call->nb_benders_opt_cuts;
			
			if (cplex_status == 11) // Callbacks timed out	
			{
				// Retrieve best-known upper bound
				exact_lb = lazy_call->best_sol;
				printf("CallBack T/O best_sol:%.2lf nb_benders_cuts:%d\n", exact_lb, nb_benders_cuts);
			}			
		} else {
			printf("No benders cuts because no lazy_call was used ...\n");
		}
		
		printf("Target levels (double) at end of Cplex:\n");
		for(int i=0; i<prob->GetNodeCount(); i++)
		{
			printf("y[%d]:%.2lf ",i,(double)cplex.getValue(y[i]));
		}printf("\n");
		
		double recCost=0.0;
		std::vector<int> hm_int_vec; std::vector<int> hp_int_vec;
		std::vector<int> exact_start_targets;
		for(int i=0; i<prob->GetNodeCount(); i++)
		{
			sol->SetExactTarget((int)(cplex.getValue(y[i])+0.1));
			exact_start_targets.push_back((int)cplex.getValue(y[i]+0.1));
			if(Parameters::GetModel()==5|| Parameters::GetModel()==6 || Parameters::GetModel()==7)
			{
				hm_int_vec.push_back( (int)cplex.getValue( hm[i] )+0.1 ); hp_int_vec.push_back( (int)cplex.getValue( hp[i] )+0.1 );
				
				sol->SetHp( (int)cplex.getValue( hp[i] )+0.1 ); sol->SetHm( (int)cplex.getValue( hm[i] )+0.1 );
			}
		}
		printf("Finalizing exact model procedure, computing the recourse of the best solution found ...\n"); //getchar();
		if(Parameters::GetModel()==5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
			recCost = r->Calculate(exact_start_targets,hm_int_vec,hp_int_vec); //computes and stores the unsatisfied requests
		else recCost = r->Calculate(exact_start_targets); //computes and stores the unsatisfied requests
		sol->SetUnsatisfiedRequests(r->GetUnsatisfiedReqVec());
		
		int total_bikes = std::accumulate(exact_start_targets.begin(),exact_start_targets.end(),0);
		
		if( Parameters::GetModel() == 5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
			printf("McfpSolver:%s GraphType%s Stations:%d Scenarios:%d Model:%d re:%d status:%d UpperBound:%.3lf CplexSol:%.3lf Theta:%.3lf WS:%.2lf RecCost:%.3lf unsatisfiedReq:%.1lf Budget:%.2lf%%CapTot nbnodes:%d BendersCuts:%d time:%.3lf\n",Parameters::GetSolverName(),r->graph_type == 4 ? "MultiSourceRecourse" : "SingleSourceRecourse", prob->GetNodeCount(),r->scs->GetScenarioCount(),Parameters::GetModel(),(int)re, cplex_status, best_upper_bound,sol_value, (double)cplex.getValue(theta),prob->GetUpperBound(), recCost, sol->GetExactUnsatisfiedReq(), Parameters::GetBudget(), cplex_nb_nodes, nb_benders_cuts, time_taken);
		else
			printf("McfpSolver:%s GraphType:%s Stations:%d Scenarios:%d Model:%d re:%d status:%d UpperBound:%.3lf CplexSol:%.3lf Theta:%.3lf WS:%.2lf RecCost:%.3lf unsatisfiedReq:%.1lf nbnodes:%d BendersCuts:%d time:%.3lf\n",Parameters::GetSolverName(),r->graph_type == 4 ? "MultiSourceRecourse" : "SingleSourceRecourse", prob->GetNodeCount(),r->scs->GetScenarioCount(),Parameters::GetModel(),(int)re, cplex_status, best_upper_bound,sol_value, (double)cplex.getValue(theta),prob->GetUpperBound(), recCost, sol->GetExactUnsatisfiedReq(), cplex_nb_nodes, nb_benders_cuts, time_taken);
		
		printf("CapTot:%d Qtot:%d TotalBikesUsed:%d Target levels:\n",prob->GetCapTot(),prob->GetQtot(),total_bikes);
		for(int i=0; i<sol->GetNodeCount(); i++)
			printf("y[%d]:%d ",i,exact_start_targets[i]);
		printf("\n");
		if(Parameters::GetModel()==5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
		{
			int increased = 0; int decreased = 0;
			for(int i=0; i<sol->GetNodeCount(); i++)
				if((int)cplex.getValue(hm[i]) > 0 || (int)cplex.getValue(hp[i]) > 0)
				{
					printf("hm[%d]:%d hp[%d]:%d ",i,(int)cplex.getValue(hm[i]),i,(int)cplex.getValue(hp[i]));
					increased += (int)cplex.getValue(hp[i]); decreased += (int)cplex.getValue(hm[i]);
					//prob->ModifyCap( i, (int)cplex.getValue(hp[i]) - (int)cplex.getValue(hm[i]) );
				}
			printf("\nTotalDocks Increased:%d Decreased:%d\n",increased,decreased);
		}
		//exit(1);

	//End try
	}catch (const IloException& ex) {
        // Handle IloException
        std::cerr << "IloException caught: " << ex.getMessage() << std::endl;
    } catch (const std::exception& ex) {
        // Handle other standard C++ exceptions
        std::cerr << "Standard C++ exception caught: " << ex.what() << std::endl;
    } catch (...) {
        // Handle any other unhandled exceptions
        std::cerr << "Unknown exception caught." << std::endl;
    }
	
}

void ExactTlb::Clear()
{
    if (lazy_call != NULL)
    {
		cplex.remove(lazy_call);
		delete lazy_call;
    }
	
	
}

void ExactTlb::SetMIPstartDP(IloEnv env)
{
	std::vector<int> capacities;
	for(int i=0;i<prob->GetNodeCount();i++)
		capacities.push_back(prob->GetNode(i)->stationcapacity);
	
	int sum_capacities = std::accumulate(capacities.begin(),capacities.end(),0);
	
	clock_t start_time_dp = clock();
	DP dp;
	std::vector<int> start_targets;
	double cost = dp.GetCostDP(capacities,prob->GetQtot(),start_targets);
	clock_t end_time_dp = clock();
	double time_taken = (double)(end_time_dp - start_time_dp) / CLOCKS_PER_SEC;
	printf("DP cost:%.2lf Qtot:%d sum h_i:%d TimeTaken:%.1lf NewTargets:\n",cost,prob->GetQtot(),sum_capacities,time_taken);

	//x_i
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for(int i=0;i<prob->GetNodeCount();i++)
	{
		startVar.add(y[ i ]); 
		startVal.add(start_targets[ i ]);	
		if(start_targets[i] > prob->GetNode(i)->stationcapacity)
		{
			//printf("i:%d start_target:%d cap:%d\n",i,start_targets[i],prob->GetNode(i)->stationcapacity); 
			
			printf("Wrong start target\n"); exit(1);
		}
	}

	//Always repair!
	cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartRepair);

	startVal.end();
	startVar.end();

	std::vector<int> hp_int_vec; std::vector<int> hm_int_vec;
	if(Parameters::GetModel()==5)
	{
		IloNumArray valueshp(env); IloNumArray valueshm(env);
		cplex.getValues(valueshp,hp); cplex.getValues(valueshm,hm);
		for(int i=0;i<hp.getSize();i++)
		{
			hp_int_vec.push_back( (int)valueshp[i] );
			hm_int_vec.push_back( (int)valueshm[i] );
		}
		valueshm.end(); valueshp.end();
	}		

	bool addedOpt = AddOptCut(env,start_targets,hm_int_vec,hm_int_vec);
	printf("AddOptCut in Mip start from DP:%d\n",addedOpt); //getchar();
}

void ExactTlb::SetMIPstart(IloEnv env, Sol & s)
{
  printf("Setting Mip start ...\n");
  std::vector<int> start_targets(prob->GetNodeCount(),0);
  for(int i=0;i<prob->GetNodeCount();i++)
  {
	int target=0;
	for(int e=0;e<scs->GetScenarioCount();e++)
		target+=s.GetTarget(e,i);
	start_targets[i] = (int)(std::round(target/scs->GetScenarioCount()));
  }
  
  //x_i
  IloNumVarArray startVar(env);
  IloNumArray startVal(env);
  for(int i=0;i<prob->GetNodeCount();i++)
  {
	startVar.add(y[ i ]); 
	startVal.add(start_targets[ i ]);	
  }

  printf("Setting total rec:%.1lf\n",prob->GetUpperBound()+0.01);
  startVar.add(theta);
  startVal.add(prob->GetUpperBound()+0.01);
  //cplex.addMIPStart(startVar, startVal,IloCplex::MIPStartSolveFixed);
  
  //Always repair!
  cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartRepair);
   
  startVal.end();
  startVar.end();

}

void ExactTlb::UpdateMIPStart(IloEnv env)
{
	//printf("Updating Mip start ... \n");
	//x_i
	IloNumVarArray startVar(env);
	IloNumArray startVal(env);
	for(int i=0;i<prob->GetNodeCount();i++)
	{
		startVar.add(y[ i ]); 
		startVal.add((cplex.getValue(y[ i ])+0.1));	
	}

	printf("Updating total rec:%.1lf\n",(double)cplex.getObjValue()+0.01);
	startVar.add(theta);
	startVal.add(cplex.getObjValue()+0.01);
	
    //Always repair!	
	cplex.addMIPStart(startVar, startVal, IloCplex::MIPStartRepair);
	//cplex.addMIPStart(startVar, startVal,IloCplex::MIPStartSolveFixed);
	//cplex.addMIPStart(startVar,startVal);

	startVal.end();
	startVar.end();	
}

bool ExactTlb::AddOptCut(IloEnv env, std::vector<int> & targets, std::vector<int> & hm, std::vector<int> & hp)
{
	printf("In add OptCut ...\n");
	
	bool AddedCut = false;
		
	//Separate the Benders cut
	std::vector<double> duals(prob->GetNodeCount()+1,0.0);
	std::vector<double> capacity_duals;
	
	double RecCost = r->Calculate(targets,duals,capacity_duals,hm,hp);
	IloExpr expr(env);
	double rhs = 0.0; double lhs = 0.0;
	for(int i=0;i<sol->GetNodeCount();i++)
	{
		expr -= duals[i] * y[ i ];
		rhs += duals[i] * targets[i];
	}
	//printf("Coefficients:\n");
	//for(int i=0;i<sol->GetNodeCount()+1;i++)
		//printf("%.1lf ",duals[i]);
	//printf("\n");
	
	expr -= duals[prob->GetNodeCount()];	
	expr += theta;
		
	// lhs = cplex.getValue(theta);
	lhs = r->Calculate(targets,hm,hp);
	rhs += duals[(int)targets.size()];
	
	if(lhs > rhs + 0.0001)
	{
		model.add(expr<=0);
		std::cout << "lhs: " << lhs << " rhs:" << rhs << " BCut:" << expr << "<= 0" << std::endl;
		//getchar();
		nb_benders_cuts++;
		AddedCut = true;
		if(nb_benders_cuts%250 ==0)
			printf("nb benders so far:%d\n",nb_benders_cuts);
	}
	expr.end();
	//getchar();
	return AddedCut;
}