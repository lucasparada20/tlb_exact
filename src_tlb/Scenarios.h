#ifndef SCENARIOS
#define SCENARIOS

#include <algorithm>
#include "Parameters.h"
#include "Trip.h"

class Scenario
{
public:
	Scenario(int nb_stations, int e): scenario_int(e), artificial(false), station_count(nb_stations), lb(std::vector<int>(nb_stations,-1)), ub(std::vector<int>(nb_stations,1000)), bike_lb(-1), bike_ub(-1) {}
	Scenario(){}
	
	void Init(int nb_stations, int e){ station_count = nb_stations; scenario_int = e;
	artificial = false; lb.resize(nb_stations,-1); ub.resize(nb_stations,-1);
	bike_lb = -1; bike_ub =-1;
	}
	
	bool AreTripsArtificial() {return artificial;}
	int GetScenarioNo(){return scenario_int;}
	
	//Trips
	int GetODTripCount(){ return (int)trips.size(); }
	Trip * GetODTrip(int i){ return &trips[i];}
	std::vector<Trip> * GetODTrips(){return &trips;}
	void AddODtrip(Trip trip){trips.push_back(trip);}
	
	//Getters and Setters for LB and UB
	void SetTargetLb(int station, int v){lb[station] = v;}
	void SetTargetUb(int station, int v){ub[station] = v;}
	int GetTargetLb(int station){ return lb[station];}
	int GetTargetUb(int station){ return ub[station];}
	
	void SetBikeLb(int v){bike_lb=v;}
	void SetBikeUb(int v){bike_ub=v;}
	int GetBikeLb(){ return bike_lb;}
	int GetBikeUb(){ return bike_ub;}

	void CalculateEventTimes(std::vector<std::vector<int>> & events)
	{
		std::vector<Trip> trips_start = trips;
		std::vector<Trip> trips_end = trips;
		
		//printf("SizeOf events:%d trips_start:%d trips_end:%d\n",(int)events.size(),(int)trips_start.size(),(int)trips_end.size());
		
		std::sort(trips_start.begin(), trips_start.end(), [](const Trip& t1, const Trip& t2) {
			return t1.start_t < t2.start_t;
		});
		
		std::sort(trips_end.begin(), trips_end.end(), [](const Trip& t1, const Trip& t2) {
			return t1.end_t < t2.end_t;
		});

		std::vector<std::vector<int>> events_start;
		events_start.resize((int)events.size());
		
		for (Trip & t : trips_start) 
		{
			if(t.start_no <= 0 || t.start_no-1 > (int)events_start.size()){
				printf("The following trip has time out of events.size(). Most likely situation: You forgot to do N++ in the first line of the instance (int the TLB, the depot is 0 and the first station is 1 ...)\n");
				t.Show();
				printf("t.start_no-1:%d events_start.size():%d\n",
						t.start_no-1,(int)events_start.size());
				printf("SizeOf events:%d trips_start:%d trips_end:%d scenario:%d\n",
						(int)events.size(),(int)trips_start.size(),(int)trips_end.size(),scenario_int);
				exit(1);
			}
			if (events_start[t.start_no-1].empty() || events_start[t.start_no-1].back() < t.start_t)
				events_start[t.start_no-1].push_back(t.start_t);
		}
		
		std::vector<std::vector<int>> events_end;
		events_end.resize((int)events.size());
		for (const Trip & t : trips_end) {
			if (events_end[t.end_no-1].empty() || events_end[t.end_no-1].back() < t.end_t)
				events_end[t.end_no-1].push_back(t.end_t);
		}

			// Merge and remove duplicates
			for (int i = 0; i < station_count; ++i) {
				std::merge(events_start[i].begin(), events_start[i].end(),
						   events_end[i].begin(), events_end[i].end(),
						   std::back_inserter(events[i]));

			// Remove duplicates from events[i]
			auto last = std::unique(events[i].begin(), events[i].end());
			events[i].erase(last, events[i].end());
		}
	}
	
private:
	int scenario_int;
	bool artificial;
	std::vector<Trip> trips;
	int station_count;
	int bike_lb; int bike_ub;
	std::vector<int> lb; std::vector<int> ub;
	
};

class Scenarios {
public:	
    Scenarios() : nbTrips(-1), artificial(true) {}

    void SetTripsCount(int v) { nbTrips = v; }
    int GetTripsCount() { return nbTrips; }	

    void AddScenario(Scenario &sce) {trips_vec.push_back(sce);}
    void SetScenarios(int e) { trips_vec.resize(e); }
    int GetScenarioCount() { return (int)trips_vec.size(); }
	
    Scenario *GetScenarioObject(int e) { return &trips_vec[e]; }

    bool GetTypeOfTrips() { return artificial; }
    void SetTypeOfTrips(bool b) { artificial = b; }

    void SetTargetLb(int scenario, int station, int v) { trips_vec[scenario].SetTargetLb(station, v); }
    void SetTargetUb(int scenario, int station, int v) { trips_vec[scenario].SetTargetUb(station, v); }
    int GetTargetLb(int scenario, int station) { return trips_vec[scenario].GetTargetLb(station); }
    int GetTargetUb(int scenario, int station) { return trips_vec[scenario].GetTargetUb(station); }

    void SetBikeLb(int scenario, int v) { trips_vec[scenario].SetBikeLb(v); }
    void SetBikeUb(int scenario, int v) { trips_vec[scenario].SetBikeUb(v); }
    int GetBikeLb(int scenario) { return trips_vec[scenario].GetBikeLb(); }
    int GetBikeUb(int scenario) { return trips_vec[scenario].GetBikeUb(); }

private:
    bool artificial;
    std::vector<Scenario> trips_vec;
    int nbTrips;
};


#endif
