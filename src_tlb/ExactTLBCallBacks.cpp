#include "ExactTLBCallBacks.h"

ExactTlbLazyCallBack::ExactTlbLazyCallBack(IloEnv env, Prob * _prob, RecourseCalculation * _r, IloNumVarArray _y, IloNumVarArray _hm, IloNumVarArray _hp, IloNumVar _theta) : IloCplex::LazyConstraintCallbackI(env), prob(_prob), r(_r), y(_y), hm(_hm), hp(_hp), theta(_theta)
{
	best_sol = 0.0;
	nb_benders_opt_cuts = 0;
}

ExactTlbLazyCallBack::ExactTlbLazyCallBack(IloEnv env, Prob * _prob, RecourseCalculation * _r, IloNumVarArray _y, IloNumVar _theta) : IloCplex::LazyConstraintCallbackI(env), prob(_prob), r(_r), y(_y), theta(_theta)
{
	best_sol = 0.0;
	nb_benders_opt_cuts = 0;
}

void ExactTlbLazyCallBack::main()
{
	/*printf("ExactSbrpodLazyCallBack obj:%.3lf\n", (double)getObjValue());
	printf("vars:\n");
	for (int i = 0; i < y.getSize(); i++) {
        double value = getValue(y[i]);
		if(value>0.01)
			printf("y[%i]:%.1lf ",i,value);
    } printf("\n");
	printf("theta:%.2lf\n",(double)getValue(theta));*/
	
	//Parameters::GetModel()==5
	/*printf("size of hp:%d hm:%d\n",(int)hm.getSize(),(int)hp.getSize());
	for(int i=0;i<hp.getSize();i++)
	{
		printf("hm[%d]:%d hp[%d]:%d ",i,(int)getValue(hm[i]),i,(int)getValue(hp[i]));
	}printf("\n");*/
	//getchar();	
	
	IloNumArray values(getEnv());
	getValues(values,y);
	std::vector<int> targets;
	for(int i = 0 ; i < y.getSize();i++)
	{
		targets.push_back((int)(values[i]+0.1));
	}
	values.end();
	
	//Separate the Benders cut
	std::vector<double> duals(prob->GetNodeCount()+1,0.0);
	//printf("nodes:%d y_vars:%d targets:%d\n",prob->GetNodeCount(),(int)y.getSize(),(int)targets.size());
	std::vector<double> capacity_duals;
	std::vector<int> hp_int_vec; std::vector<int> hm_int_vec;
	if(Parameters::GetModel()==5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
	{
		IloNumArray valueshp(getEnv()); IloNumArray valueshm(getEnv());
		getValues(valueshp,hp); getValues(valueshm,hm); //getValues returns, populates doubles not integers!
		for(int i=0;i<hp.getSize();i++)
		{
			hp_int_vec.push_back( (int)(valueshp[i]+0.5) );
			hm_int_vec.push_back( (int)(valueshm[i]+0.5) );
		}
		valueshm.end(); valueshp.end();
		if(hp_int_vec.size() != prob->GetNodeCount() || hm_int_vec.size() != prob->GetNodeCount())
		{
			std::cerr << "hm_int_vec.size():" << hm_int_vec.size() << "  hp_int_vec.size():" << hp_int_vec.size() << std::endl; exit(1);
		}
		/*for(int i=0;i<prob->GetNodeCount();i++)
			if( hm_int_vec[i]>0 || hp_int_vec[i]>0 )
			printf("hm%d:%d hp%d:%d ",i,hm_int_vec[i],i,hp_int_vec[i]);
		printf("\n");*/
	}	
	
	if(nb_benders_opt_cuts==0)
	{
		r->y0.resize(prob->GetNodeCount(),0);
		for(int i = 0 ; i < prob->GetNodeCount();i++)
			r->y0[i] = targets[i];
	}
	double RecCost = r->Calculate(targets,duals,capacity_duals,hm_int_vec,hp_int_vec);
	//double RecCost = r->CalculateParetoOpt(targets,duals);

	IloExpr expr(getEnv());
	double rhs = 0.0; double lhs = 0.0; 
	for(int i=0;i<prob->GetNodeCount();i++)
	{
		expr -= duals[i] * y[ i ];
		rhs += duals[i] * targets[i];
	}	
	expr -= duals[prob->GetNodeCount()];	
	expr += theta;

	lhs = (double)getValue(theta);
	rhs += duals[prob->GetNodeCount()];
	
	if(Parameters::GetModel()==5 || Parameters::GetModel()==6 || Parameters::GetModel()==7)
	{
		for(int i=0;i<prob->GetNodeCount();i++)
		{
			expr += ( - hp[i] + hm[i] ) * capacity_duals[i];
			rhs += ( hp_int_vec[i] - hm_int_vec[i] ) * capacity_duals[i];
		}
	}
	//for(int i = 0; i < prob->GetNodeCount(); i++)
		//printf("t%d:%d ",i,targets[i]);
	//printf("\n");	
	//for(int i=0;i<=prob->GetNodeCount();i++)
		//printf("duals[%d]:%.1lf ",i,duals[i]);
	//printf("\n");

	//printf("Callback RecCost:%.3lf lhs:%.3lf should be <= than rhs:%.3lf\n",RecCost,lhs,rhs);
	//To DBG ...
	//double goodRecCost = r->Calculate(targets);
	//if( std::abs(goodRecCost-rhs) > 0.01 )
	//{
		//printf("Callback Wrong goodRecCost:%.3lf rhs:%.3lf\n",goodRecCost,rhs); exit(1);
	//}
	if(lhs > rhs + 0.0001)
	{
		add(expr<=0);
		// Output with three-point precision
		//std::cout << "lhs: " << std::fixed << std::setprecision(3) << lhs 
        //<< " rhs: " << std::fixed << std::setprecision(3) << rhs 
        //<< " BCut nb:" << nb_benders_opt_cuts << " : " << std::fixed << std::setprecision(1) << expr << " <= 0" << std::endl;
		nb_benders_opt_cuts++;
		if(nb_benders_opt_cuts%250 ==0)
			printf("nb benders so far:%d\n",nb_benders_opt_cuts);
		//getchar();
	}
	expr.end();
	
	if(rhs > best_sol)
	{
		r->y0.clear(); r->y0.resize(prob->GetNodeCount(),0);
		//best_sol = RecCost;
		best_sol = rhs;
		best_solution.clear();
		//printf("new_sol RecCost:%.3lf\n",RecCost);
		printf("new_sol RecCost:%.3lf\n",rhs);
		for(int i = 0; i < prob->GetNodeCount(); i++)
		{
			best_solution.push_back( targets[i] );
			r->y0[i] = targets[i];
			//printf("t%d:%d ",i,targets[i]);
		}
		//printf("\n");
			
	}		
}

