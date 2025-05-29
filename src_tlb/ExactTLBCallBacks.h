#ifndef TLB_CALLBACKS_H
#define TLB_CALLBACKS_H

#include <algorithm>
#include <vector>
#include <ilcplex/ilocplex.h>
#include "Scenarios.h"
#include "RecourseCalculation.h"


class ExactTlbLazyCallBack : public IloCplex::LazyConstraintCallbackI
{
	public:
		ExactTlbLazyCallBack(IloEnv env, Prob * prob, RecourseCalculation * r, IloNumVarArray y, IloNumVar theta);
		ExactTlbLazyCallBack(IloEnv env, Prob * prob, RecourseCalculation * r, IloNumVarArray y, IloNumVarArray hm, IloNumVarArray hp, IloNumVar theta);
		
		//~ExactTlbLazyCallBack();

		IloCplex::CallbackI *duplicateCallback() const
		{
        	return new (getEnv()) ExactTlbLazyCallBack(*this);
        }

        void main();
		double best_sol;
		int nb_benders_opt_cuts;
		std::vector<int> best_solution;
		//IloConstraintArray added_constraints_array;

    private:
		Prob * prob;
		RecourseCalculation * r;
		IloNumVarArray y;
		IloNumVar theta;
		IloNumVarArray hm;
		IloNumVarArray hp;
				
};

#endif