
#ifndef NODE_SBRPOD
#define NODE_SBRPOD

#include "constants.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Trip.h"
#include <cstring>
#include "Parameters.h"
#include <unordered_map>
#include <unordered_set>

class Node
{
	public:
		Node() : id(-1), bss_id(-1), stationcapacity(-1), occupancy(-1), no(-1), x(0), y(0), target(-1), min_station_bikes(-1), max_station_bikes(1000) {}
	
		int id;			//from 0 to n-1
		int bss_id;
		int stationcapacity;   // Station capacity
		int occupancy;   //Initial occupational level of the station
		int no;			//personnal identifier, put what you want
		double x;
		double y;
		int target;

		int min_station_bikes; int max_station_bikes;		

		void SetLb(int v){ min_station_bikes=v;}
		void SetUb(int v){ max_station_bikes=v;}
		int GetLb(){ return min_station_bikes;}
		int GetUb(){ return max_station_bikes;}
		
		void ModifyCap(int new_cap){ stationcapacity = new_cap; }
		void ModifyTgt(int new_tgt){ target = new_tgt; }
		void Show()
		{
			printf("Node:%d cap:%d tgt:%d\n",id,stationcapacity,target);
			if(target>stationcapacity)
				printf("In Node: Wrong tgt > cap. Exiting");
			//printf("Node:%d cap:%d out_OD:%d in_OD:%d realOTrips:%d realITrips:%d\n", id, stationcapacity, (int)out_ODtrips.size(),(int)in_ODtrips.size(),(int)real_out_ODtrips.size(),(int)real_in_ODtrips.size());
		}

};


#endif
