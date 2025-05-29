#ifndef PROBLEM_DEFINITION
#define PROBLEM_DEFINITION

#include "NodeTLB.h"

class Prob
{
public:
	Prob() : Qtot(-1), ScenarioCount(-1), ODtripCount(-1), Tmax(-1) {}
	
	Node * GetNode(int i){ return &nodes[i];}
	void ModifyCap( int i, int new_cap ){ nodes[i].ModifyCap(new_cap); }
	int GetNodeCount(){ return (int)nodes.size(); }
	void AddNode(Node n){nodes.push_back(n);}
	void ShowNodes(){for(Node & n : nodes) n.Show();}
	
	void SetStationLb(int i, int v){ nodes[i].SetLb(v);}
	void SetStationUb(int i, int v){ nodes[i].SetUb(v);}
	int GetStationLb(int i){ return nodes[i].GetLb();}
	int GetStationUb(int i){ return nodes[i].GetUb();}
	
	void SetQtot(int i){ Qtot = i; }
	int GetQtot(){ return Qtot; }
	
	void SetCapTot(int i){ hTot = i; }
	int GetCapTot(){ return hTot; }
	
	void SetUpperBound(double ub){ UpperBound = ub; }
	double GetUpperBound(){ return UpperBound; }
	
	void SetTmax(int t){Tmax = t;}
	int GetTmax(){return Tmax;}
	
	void SetTgtLvlVectors(int e)
	{
		pick_tgt_lvls.resize(e,std::vector<double>());
		del_tgt_lvls.resize(e,std::vector<double>());
	}
	void StorePickTgtLvls(int e, std::vector<double> & pick_tgt_level)
	{ 
		if(e>=pick_tgt_lvls.size())
		{
			printf("in Prob e:%d pickTgtSize:%d. Phil Collins (1989).\n",e,(int)pick_tgt_lvls.size()); exit(1);
		}
		pick_tgt_lvls[e] = std::move(pick_tgt_level); 
	}
	void StoreDelTgtLvls(int e, std::vector<double> & del_tgt_level)
	{ 
		if(e>=del_tgt_lvls.size())
		{
			printf("in Prob e:%d delTgtSize:%d. Phil Collins (1989).\n",e,(int)del_tgt_lvls.size()); exit(1);
		}
		del_tgt_lvls[e] = std::move(del_tgt_level); 
	}
	
	std::vector<std::vector<double>> & GetPickTgtLvls() {return pick_tgt_lvls;}
	std::vector<std::vector<double>> & GetDelTgtLvls() {return del_tgt_lvls;}
	
private:
	
	int Tmax;
	std::vector<Node> nodes;
	int Qtot;
	int hTot;
	int ScenarioCount;
	int ODtripCount;
	double UpperBound;
	
	std::vector<std::vector<double>> pick_tgt_lvls;
	std::vector<std::vector<double>> del_tgt_lvls;	
};

#endif