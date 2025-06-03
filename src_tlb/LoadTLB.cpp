
#include "LoadTLB.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <float.h>
#include <iterator>
#include <sstream>
#include <ctime>
#include <iomanip>      // std::setprecision
#include "json.hpp"

void LoadTLB::LoadTargets(std::vector<int> & targets, char * filename)
{
	printf("In the Load Targets:\n");
	FILE * ff = fopen(filename,"r");
 	if(!ff)
	{
		printf("Error in the input filename:%s\n", filename);
		exit(1);
	}
	char line[100];
	while (fgets(line, 100, ff)) 
	{
		int t;
		if(sscanf(line,"%d",&t)==1)
		{
			targets.push_back(t);
		} else{
			printf("Error parsing line:%s Load Targets\n", line);
		}		
	}
}

void LoadTLB::LoadScenarios(Prob * pr, Scenarios * scs, char * filename)
{
	printf("In the Load Scenarios:\n");
	FILE * ff = fopen(filename,"r");
 	if(!ff)
	{
		printf("Error in the input filename:%s\n", filename);
		exit(1);
	}

	//std::string to to retrieve the integer next to the last 'e' character
    std::string strFilename(filename);
    // Find the last position of the 'e' character
    size_t ePos = strFilename.rfind('e');

    // Check if 'e' is found
	int nb_scenarios = 0;
    if (Parameters::CalculateAvgScenario()==0 && ePos != std::string::npos) {
        // Extract the substring after the last 'e'
        std::string substring = strFilename.substr(ePos + 1);
        
        //std::cout << "Extracted substring: " << substring << std::endl;

        try {
            int intValue = std::stoi(substring);

            // Set the scenario count using Parameters::SetScenarioCount
            nb_scenarios = intValue;
	
            std::cout << "Parsed integer of nb_scenario:" << nb_scenarios << std::endl;
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Failed to parse nb_scenario integer." << std::endl;
        }
    } else {
        std::cerr << "Error: 'e' character not found." << std::endl;
    } 
	if(Parameters::CalculateAvgScenario()==1) nb_scenarios = 1;

	printf("Instance format to read:%c\n",Parameters::GetInstanceFormat());
	std::map<int,std::vector<Trip>> ODtrip_map;
	int idx=0;
	char line[100];
	while (fgets(line, 100, ff)) 
	{
		char type; int start_no, end_no, start_t, end_t, sce;
		if (sscanf(line, "%d %d %d %d %d", &start_no, &end_no, &start_t, &end_t, &sce) == 5 && Parameters::GetInstanceFormat() == 'T') 
		{
			// Check if the values fit into 16 bits
			if (start_no < std::numeric_limits<int16_t>::min() || start_no > std::numeric_limits<int16_t>::max() ||
				end_no < std::numeric_limits<int16_t>::min() || end_no > std::numeric_limits<int16_t>::max() ||
				start_t < std::numeric_limits<int16_t>::min() || start_t > std::numeric_limits<int16_t>::max() ||
				end_t < std::numeric_limits<int16_t>::min() || end_t > std::numeric_limits<int16_t>::max() ||
				sce < std::numeric_limits<int16_t>::min() || sce > std::numeric_limits<int16_t>::max()) 
			{
				// Values are out of range for int16_t
				printf("Error: Integers out of range for int16_t.\nExiting ...");
				exit(1);
			}

			Trip OD;
			OD.start_no = start_no;
			OD.end_no = end_no;
			
			if(OD.start_no == -1 || OD.end_no == -1)
			{
				printf("Could not match bss_id of:%s in Load and instance format:%c\nExiting ...",line,Parameters::GetInstanceFormat()); exit(1);
			}
			OD.start_t = (int32_t)start_t; OD.end_t = (int32_t)end_t;
			if(OD.start_t == end_t)
			{
				continue; //OD.Show(); throw MyException("Trip with same start, end time");
			}
					
			// Push the modified Trip struct to the vector
			OD.idx = idx; 
			idx++;
			//OD.Show();
			ODtrip_map[sce].push_back(OD);
		} 
		else if (sscanf(line, "%d %d %d %d %d", &start_no, &end_no, &start_t, &end_t, &sce) == 5 && Parameters::GetInstanceFormat() == 'R') 
		{
			// Check if the values fit into 16 bits
			if (start_no < std::numeric_limits<int16_t>::min() || start_no > std::numeric_limits<int16_t>::max() ||
				end_no < std::numeric_limits<int16_t>::min() || end_no > std::numeric_limits<int16_t>::max() ||
				start_t < std::numeric_limits<int16_t>::min() || start_t > std::numeric_limits<int16_t>::max() ||
				end_t < std::numeric_limits<int16_t>::min() || end_t > std::numeric_limits<int16_t>::max() ||
				sce < std::numeric_limits<int16_t>::min() || sce > std::numeric_limits<int16_t>::max()) 
			{
				// Values are out of range for int16_t
				printf("Error: Integers out of range for int16_t.\nExiting ...");
				exit(1);
			}

			auto it_start = bss_id_map.find(start_no);
			auto it_end = bss_id_map.find(end_no);
			
			if(it_start == bss_id_map.end() || it_end == bss_id_map.end() || start_t == end_t )
			{
				printf("Could not match bss_id or equal start, end times in %s. In Load and instance format:%c\nExiting ...",line,Parameters::GetInstanceFormat()); exit(1);
			}
			Trip pickup; Trip delivery;
			pickup.start_no = (int16_t)it_start->second;
			delivery.end_no = (int16_t)it_end->second;
			pickup.start_t = (int16_t)start_t; delivery.end_t = (int16_t)end_t;

			pickup.idx = idx;
			idx++;
			delivery.idx = idx;
			idx++;
			//printf("%s\n",line);
			//pickup.Show(); delivery.Show();
			ODtrip_map[sce].push_back(pickup);
			ODtrip_map[sce].push_back(delivery);
		}		
		else if(sscanf(line, "%c %d %d %d", &type,&start_no,&start_t,&sce)==4)
		{
			// Check if the values fit into 16 bits
			if (start_no < std::numeric_limits<int16_t>::min() || start_no > std::numeric_limits<int16_t>::max() ||
				start_t < std::numeric_limits<int16_t>::min() || start_t > std::numeric_limits<int16_t>::max() ||
				sce < std::numeric_limits<int16_t>::min() || sce > std::numeric_limits<int16_t>::max()) 
			{
				// Values are out of range for int16_t
				printf("Error: Integers out of range for int16_t.\nExiting ...");
				exit(1);
			}

			Trip OD;
			auto it_start = bss_id_map.find(start_no);
			if(type=='p')
			{
				OD.start_no = (int16_t)start_no;
				OD.start_t = (int16_t)start_t;
			}else{//type=='q'
				OD.end_no = (int16_t)start_no;
				OD.end_t = (int16_t)start_t;				
			}
			OD.idx = idx;

			idx++;
			
			ODtrip_map[sce].push_back(OD);
			
		}
		else{
			printf("Error parsing line:%s in instance format:%c\n", line,Parameters::GetInstanceFormat());
		}
	}	
	
	//Create scenario data structures
	scs->SetScenarios(nb_scenarios);	
	scs->SetTypeOfTrips(false);
	
	printf("Nb scenarios:%d %s:\n",scs->GetScenarioCount(),(Parameters::GetInstanceFormat() == 'T')? "trips" : "requests");
	//printf("Size of Bss_id_map:%d\n",(int)bss_id_map.size());
	int nb_total_trips = 0; int stations = pr->GetNodeCount(); //bss_id_map.size();
	for(int e=0;e<nb_scenarios;e++)
	{
		Scenario * sce = scs->GetScenarioObject(e);
		//sce->Init(stations-1,e);
		sce->Init(stations,e);
		//printf("%d ",nb);
		for (int k = 0; k < ODtrip_map[e].size(); k++)
			sce->AddODtrip(ODtrip_map[e][k]);
		
		printf("sce:%d trips:%d\n", e, sce->GetODTripCount());
		nb_total_trips += sce->GetODTripCount();	
	}
	scs->SetTripsCount( nb_total_trips );

	
}

void LoadTLB::Load(Prob * pr, Scenarios * scs, char * filename)
{
	printf("In the Load instance:\n");
	FILE * ff = fopen(filename,"r");
 	if(!ff)
	{
		printf("Error in the input filename:%s\n", filename);
		exit(1);
	}

	//std::string to to retrieve the integer next to the last 'e' character
    std::string strFilename(filename);
    // Find the last position of the 'e' character
    size_t ePos = strFilename.rfind('e');

    // Check if 'e' is found
	int nb_scenarios = 0;
    if (ePos != std::string::npos) {
        // Extract the substring after the last 'e'
        std::string substring = strFilename.substr(ePos + 1);
        
        //std::cout << "Extracted substring: " << substring << std::endl;

        try {
            int intValue = std::stoi(substring);

            // Set the scenario count using Parameters::SetScenarioCount
            nb_scenarios = intValue;
	
            std::cout << "Parsed integer of nb_scenario:" << nb_scenarios << std::endl;
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Failed to parse nb_scenario integer." << std::endl;
        }
    } else {
        std::cerr << "Error: 'e' character not found." << std::endl;
    }
	
	int stations, Q, Qtot = 0;
	char line[100];
	fgets(line,100,ff); //Read N
	sscanf(line, "%d\n", &stations);
	fgets(line,100,ff); //Read Q
	sscanf(line, "%d\n", &Qtot);
	printf("From file Qtot:%d N:%d\n",Qtot,stations);
	pr->SetQtot(Qtot);
	
	char big_line[10000];
	fgets(big_line,10000,ff); //Read stations cap
	//std::cout << big_line;
	std::string line_str(big_line);
	std::istringstream is_cap( line_str );
	std::vector<int> station_cap( ( std::istream_iterator<int>( is_cap ) ), ( std::istream_iterator<int>() ) );	
	
	long last_pos = ftell(ff);  // Save position before first read
	fgets(big_line, 500, ff);
	char type; int start_no, end_no, start_t, end_t, sce;
	while (sscanf(big_line, "%d %d %d %d %d", &start_no, &end_no, &start_t, &end_t, &sce) != 5)
	{
		last_pos = ftell(ff);  // Update position before reading again
		if (!fgets(big_line, 500, ff)) break; // Protect against EOF
	}

	// Go back to the start of the valid line
	fseek(ff, last_pos, SEEK_SET);
	
	
	/*std::map<std::string, int> map_str_id; int cntr = 0;
	std::vector<int> initial_bikes(stations-1,0);
	fgets(big_line, 500, ff); //the depot id is `-1' -> skip	
	while(cntr < stations-1)
	{
		fgets(big_line, 500, ff);
		big_line[strcspn(big_line, "\n")] = '\0';  // Remove trailing instance file newline. Finds the position of \n and replaces it with \0
		map_str_id[std::string(big_line)] = cntr;
		//printf("cntr:%d id:%s", cntr, big_line);
        cntr++;		
	}*/
	/*if(!Parameters::GetJsonFileName().empty())
	{
		std::ifstream file_json(Parameters::GetJsonFileName());
		if (!file_json)
		{
			printf("Error in the input filename: %s\n", Parameters::GetJsonFileName().c_str());
			exit(1);
		}

		// Parse the JSON data
		nlohmann::json jsonData;
		file_json >> jsonData;
		nlohmann::json json_stations = jsonData["data"]["stations"];
		printf("Loading initial capacities from:%s NbStations:%d\n", Parameters::GetJsonFileName().c_str(), (int)json_stations.size());
		int found = 0;
		int counter = 0;
		for (const auto& json_station : json_stations)
		{
			std::string station_id = json_station["station_id"];
			//printf("station_id:%s\n",station_id.c_str());
			counter++;
			if (counter % 50 == 0)
			{
				printf("Looped through %d stations from the json file ...\n", counter);
			}

			int bikes = 0;
			auto it = map_str_id.find(station_id);
			if (it != map_str_id.end())
			{
				//std::cout << "iterator first:" << it->first << " second:" << it->second << std::endl;
				if (json_station.contains("num_bikes_available"))
				{
					bikes = json_station["num_bikes_available"];
					initial_bikes[it->second] = bikes;
					found++;
				}
			}
		}
		printf("Found %d/%d stations from the json station_status file!\n", found, stations);		
		if(Parameters::GetDelta() > 0)
		{
			Qtot = 0;
			for(int i=0;i<stations;i++)
				Qtot += initial_bikes[i];
			pr->SetQtot(Qtot + Parameters::GetDelta());
			printf("Setting Qtot as the sum of initial bikes +- Delta:%d New Qtot:%d\n",Parameters::GetDelta(),pr->GetQtot());
		}
	}*/ 	
	//exit(1);
	
	//Make intervals
	int nb_intervals = (int)960 + 1; //Max time in the generator was modified from 840 to 960 to include 2 more hours!
	std::vector<int> intervals;
	for(int i=0; i<nb_intervals; i++)
		intervals.push_back(i);
	//printf("nb_intervals:%d\n",nb_intervals);
	//printf("nb_intervals:%d\nintervals:",nb_intervals);
	//for(int i=0; i<nb_intervals; i++) printf("%d ",intervals[i]);
	//printf("\n");
	Parameters::SetTmax(nb_intervals);
	pr->SetTmax(nb_intervals);
	
	//Load ODtrips or Requests
	printf("Instance format to read:%c\n",Parameters::GetInstanceFormat());
	std::map<int,std::vector<Trip>> ODtrip_map;
	int idx=0;
	while (fgets(line, 100, ff)) 
	{
		//std::cout << line;
		//char type; int start_no, end_no, start_t, end_t, sce;
		//Remember to remove
		line[strcspn(line, "\n")] = 0;
		if (sscanf(line, "%d %d %d %d %d", &start_no, &end_no, &start_t, &end_t, &sce) == 5 && Parameters::GetInstanceFormat() == 'T') 
		{
			
			// Check if the values fit into 16 bits
			//if (start_no < std::numeric_limits<int16_t>::min() || start_no > std::numeric_limits<int16_t>::max() ||
			//	end_no < std::numeric_limits<int16_t>::min() || end_no > std::numeric_limits<int16_t>::max() ||
			//	start_t < std::numeric_limits<int16_t>::min() || start_t > std::numeric_limits<int16_t>::max() ||
			//	end_t < std::numeric_limits<int16_t>::min() || end_t > std::numeric_limits<int16_t>::max() ||
			//	sce < std::numeric_limits<int16_t>::min() || sce > std::numeric_limits<int16_t>::max()) 
			//{
				// Values are out of range for int16_t
			//	printf("Error: Integers out of range for int16_t.\nExiting\n ...");
			//	printf("%d %d %d %d %d\n", start_no, end_no, start_t, end_t, sce);
			//	exit(1);
			//}
			
			if( start_no < 0 || start_no > stations || end_no < 0 || end_no > stations)
			{
				printf("Wrong start_no:%d end_no:%d stations:%d, skipping ...\n",start_no,end_no,stations);
				continue;
				//exit(1);
			}

			Trip OD;
			OD.start_no = start_no;
			OD.end_no = end_no;
			
			OD.start_t = (int16_t)start_t; OD.end_t = (int16_t)end_t;
			if(OD.start_t == end_t)
			{
				continue; //OD.Show(); throw MyException("Trip with same start, end time");
			}
					
			// Push the modified Trip struct to the vector
			OD.idx = idx; 
			idx++;
			//OD.Show();
			ODtrip_map[sce].push_back(OD);
		} 
		else{
			printf("Error parsing line:%s in instance format:%c\n", line,Parameters::GetInstanceFormat());
		}
	}	
	printf("Trips read:%d\n",idx);
	for(int i = 0; i< stations-1; i++)
	{
		Node n;
		n.id = i;
		n.no = i+1;
		n.stationcapacity = station_cap[i+1];
		//n.occupancy = initial_bikes[i] > station_cap[i+1] ? station_cap[i+1] : initial_bikes[i];
		//n.Show();		
		pr->AddNode(n);
		
	}
	int sum_station_cap = 0;
	for (int i = 0; i < pr->GetNodeCount(); i++)
		sum_station_cap += pr->GetNode(i)->stationcapacity;
	pr->SetCapTot( sum_station_cap );
	printf("Setting CapTot:%d required for budget models ...\n",pr->GetCapTot());
	printf("If you are calculating these models, then budget is set to:%d\n",(int)std::floor( Parameters::GetBudget() * pr->GetCapTot() ));
	/*if(sum_station_cap < pr->GetQtot())
	{
		printf("SumStationCap:%d < Qtot:%d -> Setting Qtot:%d\n",sum_station_cap,pr->GetQtot(),sum_station_cap);
		Qtot = std::min(sum_station_cap, pr->GetQtot());
		pr->SetQtot( Qtot );
	}*/
	if (Qtot < 1)
	{
		printf("Phil Collins (1989) : Wrong Qtot:%d in Load",Qtot); exit(1);
	}

	if(Parameters::GetQtot() > 0.0 && Parameters::GetQtot() < 1.0)
	{
		pr->SetQtot( (double)Qtot * Parameters::GetQtot() );
		printf("Reducing Qtot:%d to %.1lf%% which is:%d\n", Qtot, Parameters::GetQtot(), pr->GetQtot());
		//getchar();
	}

	//Create scenario data structures
	scs->SetScenarios(nb_scenarios);	
	scs->SetTypeOfTrips(false);
	
	printf("Nb scenarios:%d %s:\n",scs->GetScenarioCount(),(Parameters::GetInstanceFormat() == 'T')? "trips" : "requests");
	
	int nb_total_trips = 0;
	for(int e=0;e<nb_scenarios;e++)
	{
		Scenario * sce = scs->GetScenarioObject(e);
		sce->Init(stations-1,e);
		//printf("%d ",nb);
		for (int k = 0; k < ODtrip_map[e].size(); k++)
			sce->AddODtrip(ODtrip_map[e][k]);
		
		printf("sce:%d trips:%d\n", e, sce->GetODTripCount());
		nb_total_trips += sce->GetODTripCount();	
	}
	scs->SetTripsCount( nb_total_trips );
		
	printf("Nb node objects created:%d\n",pr->GetNodeCount());
	
	fclose(ff);
	
}
