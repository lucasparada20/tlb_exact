#include "NodeTLB.h"
#include "Parameters.h"
#include "Scenarios.h"
#include "ProblemDefinition.h"
#include "Trip.h"
#include <map>

class LoadTLB
{
	public:
		void Load(Prob * pr, Scenarios * scs, char * filename);
		void LoadRealTrips(Prob * pr, Scenarios * scs, char * filename);
		
		void LoadScenarios(Prob * pr, Scenarios * scs, char * filename);
		void LoadTargets(std::vector<int> & targets, char * filename);
		
	private:
		std::map<int,int> bss_id_map;
		std::vector<int> bss_id; 
	
};
