#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <utility>
#include <set>
#include <cstdlib>
#include <algorithm>

#include <nlohmann/json.hpp>
#include <ortools/linear_solver/linear_solver.h>

#include "instance.h"

using json = nlohmann::json;
using namespace std;

Instance::Instance(const string& filename) : input_filename(filename) {
  cout << "Loading instance... \n";
  ifstream input_file("instances/" + input_filename);
  if(input_file.is_open()) {
    json json_file;
    input_file >> json_file;

    buildingCosts.resize(3);
    productionCosts.resize(3);
    routingCosts.resize(2);
    capacities.resize(2);

    buildingCosts[0] = float(json_file["parameters"]["buildingCosts"]["productionCenter"]);
    buildingCosts[1] = float(json_file["parameters"]["buildingCosts"]["automationPenalty"]);
    buildingCosts[2] = float(json_file["parameters"]["buildingCosts"]["distributionCenter"]);
    
    productionCosts[0] = float(json_file["parameters"]["productionCosts"]["productionCenter"]);
    productionCosts[1] = float(json_file["parameters"]["productionCosts"]["automationBonus"]);
    productionCosts[2] = float(json_file["parameters"]["productionCosts"]["distributionCenter"]);

    routingCosts[0] = float(json_file["parameters"]["routingCosts"]["primary"]);
    routingCosts[1] = float(json_file["parameters"]["routingCosts"]["secondary"]);

    capacityCost = float(json_file["parameters"]["capacityCost"]);

    capacities[0] = float(json_file["parameters"]["capacities"]["productionCenter"]);
    capacities[1] = float(json_file["parameters"]["capacities"]["automationBonus"]);

    clients.reserve(json_file["clients"].size());
    sites.reserve(json_file["sites"].size());

    for(auto& client : json_file["clients"]) {
      clients.push_back({client["demand"], {client["coordinates"][0], client["coordinates"][1]}});
    }

    for(auto& site : json_file["sites"]) {
      sites.push_back({site["coordinates"][0], site["coordinates"][1]});
    }

    siteSiteDistances.reserve(sites.size());
    siteClientDistances.reserve(sites.size());

    for(auto& distances : json_file["siteSiteDistances"]) {
      vector<float> row;
      row.reserve(sites.size());
      for(auto& distance : distances) {
        row.push_back(distance);
      }
      siteSiteDistances.push_back(row);
    }

    for(auto& distances : json_file["siteClientDistances"]) {
      vector<float> row;
      row.reserve(clients.size());
      for(auto& distance : distances) {
        row.push_back(distance);
      }
      siteClientDistances.push_back(row);
    }
    /**
    for(auto& cost : buildingCosts)
      cout << cost << " ";
    cout << endl;
    for(auto& cost : productionCosts)
      cout << cost << " ";
    cout << endl;
    for(auto& cost : routingCosts)
      cout << cost << " ";
    cout << endl << capacityCost << endl;

    cout << "Clients:" << clients.size() << endl;
    for(auto& client : clients)
      cout << client.first << " : " << client.second.first << ", " << client.second.second << " " << endl;
    
    cout << "Sites:" << endl;
    for(auto& site : sites)
      cout << site.first << ", " << site.second << endl;

    for(auto& row : siteSiteDistances) {
      for(auto& d : row)
        cout << d << " ";
      cout << endl;
    }

    for(auto& row : siteClientDistances) {
      for(auto& d : row)
        cout << d << " ";
      cout << endl;
    }

    **/

    input_file.close();
  }
  solution = Solution(sites.size(), clients.size());
  visited_sites = vector<bool>(sites.size(), false);
  cout << "Instance " + input_filename + " successfully loaded!\n";

  cout << "Center building costs:\n";
  cout << "c_P^b = " << (c_Pb = buildingCosts[0]) << endl;
  cout << "c_A^b = " << (c_Ab = buildingCosts[1]) << endl;
  cout << "c_D^b = " << (c_Db = buildingCosts[2]) << endl;

  cout << "Unit production costs:\n";
  cout << "c_P^p = " << (c_Pp = productionCosts[0]) << endl;
  cout << "c_A^p = " << (c_Ap = productionCosts[1]) << endl;
  cout << "c_D^p = " << (c_Dp = productionCosts[2]) << endl;

  cout << "Unit routing costs:\n";
  cout << "c_1^r = " << (c_1r = routingCosts[0]) << endl;
  cout << "c_2^r = " << (c_2r = routingCosts[1]) << endl;

  cout << "Center capacities:\n";
  cout << "u_P = " << (u_P = capacities[0]) << endl;
  cout << "u_A = " << (u_A = capacities[1]) << endl;
  cout << "c^u = " << (cu = capacityCost) << endl;
}

void Instance::solve() {
  cout << "Solving instance... \n";
  time_t begin, end;

  begin = clock();

  const int n = sites.size();
  const int m = clients.size();
  // write heuristic solver here
  
  //step 0: compute the demand
  int demand = 0;
  for(auto& client : clients)
    demand += client.first;

  //step 1: compute w_ij = d_j*\Delta_siteclient(i,j)
  vector<vector<float>> w(n, vector<float>(m));
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < m; ++j)
      w[i][j] = clients[j].first*siteClientDistances[i][j];

  //step 2: compute the order of {w_ij}_{j client} \forall possible site i
  vector<vector<int>> order_siteclient(n, vector<int>(m));
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      order_siteclient[i][j] = j;
    }
    auto compare_client = [*this, i] (float x, float y) {
      return this->siteClientDistances[i][x] < this->siteClientDistances[i][y];
    };
    sort(order_siteclient[i].begin(), order_siteclient[i].end(), compare_client);
  }

  //step 3: assign clients and distribution sites
  int n_visited_clients = 0;
  vector<bool> visited_clients(m, false);
  vector<bool> visited_sites(n, false);
  float alpha = 0.1;
  float beta = 0.3;
  float _alpha = alpha*(float)(u_P + u_A);
  float _beta = beta*(float)(u_P + u_A);

  while(n_visited_clients <= m) {
    if(demand > alpha) {
      
    } else {

    }
  }

  end = clock();

  cout << "Instance successfully solved in " << (end - begin)/(float)CLOCKS_PER_SEC << " seconds.\n";
}

void Instance::save() {
  time_t rawtime;
  time(&rawtime);
  struct tm * local_time = localtime(&rawtime);

  json json_solution = json({});

  json_solution["productionCenters"] = json::array();
  json_solution["distributionCenters"] = json::array();
  json_solution["clients"] = json::array();

  for(int i = 0; i < solution.P.size(); ++i) {
    if(solution.P[i])
      json_solution["productionCenters"].push_back({{"id", i+1}, {"automation", solution.a[i] ? 1 : 0}});
  }

  for(int i = 0; i < solution.D.size(); ++i) {
    if(solution.D[i])
      json_solution["distributionCenters"].push_back({{"id", i+1}, {"parent", solution.p[i]+1}});
  }

  for(int i = 0; i < solution.s.size(); ++i) {
    json_solution["clients"].push_back({{"id", i+1}, {"parent", solution.s[i]+1}});
  }

  string solution_filename = "solutions/solution-" + to_string(local_time->tm_hour) + ':' + to_string(local_time->tm_min) + ':' + to_string(local_time->tm_sec) + '-'  + input_filename;
  ofstream output_file(solution_filename);

  if(output_file.is_open()) {
    output_file << json_solution << std::endl;
    output_file.close();
  }
  cout << "Solution successfully saved in " << solution_filename << "\n\n";
}

Solution::Solution(int sites_number, int clients_number) {
  P = vector<bool>(sites_number, false);
  D = vector<bool>(sites_number, false);
  a = vector<bool>(sites_number, false);
  p = vector<int>(sites_number, 0);
  s = vector<int>(clients_number, 0);
}