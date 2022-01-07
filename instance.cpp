#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <utility>
#include <set>
#include <cstdlib>

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

//First, we attribute every client to a distribution center in D
map<int, vector<int> > clients_distribution(vector<bool> D, vector<pair<int, pair<float, float> > > clients, vector<vector<float> > scd){
    int k = 0;
    for (int i = 0; i < D.size(); i++){
        k += D[i];
    }

    int m = clients.size();
    int ratio = m/k;
    if (m%k != 0){
        ratio ++;
    }
    map<int, vector<int >> cd;
    vector<int> cl;
    for (int i = 0; i < clients.size(); i++){
        cl.push_back(i);
    }
    for (int i = 0; i < D.size(); i++){
        if (D[i]){
            for (int _ = 0; _ < ratio; _ ++){
                if (cl.size()!=0){
                    int min = 0;
                    for (int j = 0; j < cl.size(); j++){
                        if (scd[i][cl[j]] < scd[i][min]){
                            min = j;
                        }
                    }
                    cd[i].push_back(cl[min]);
                    cl.erase(cl.begin() + min);
                }
            }
        }
    }
    return cd;
}

// Returns the vector s
vector<int> supply_client(vector<bool> D, vector<pair<int, pair<float, float> > > clients, vector<vector<float> > scd){
    map<int, vector<int> > cd = clients_distribution(D, clients, scd);
    vector<int> sc;
    for (int i = 0; i < clients.size(); i++){
        for (auto it = cd.begin(); it != cd.end(); it++){
            for (int j = 0; j<it->second.size(); j++){
                if (i == it->second[j]){
                    sc.push_back(it->first);
                    break;
                }
            }
        }
    }
    return sc;
}

// Second, we attribute each distribution center to a production center
map<int, vector<int> > distribution_production (vector<bool> P, vector<bool> D, vector<vector<float> > ssd){
    int k = 0;
    int l = 0;
    for (int i = 0; i < P.size(); i++){
        k += D[i];
        l += P[i];
    }

    int ratio = k/l;
    if (k%l != 0){
        ratio ++;
    }
    map<int, vector<int > > dp;
    vector<int> distribution;
    for (int i = 0; i < D.size(); i++){
        if (D[i]){
            distribution.push_back(i);
        }
    }

    for (int i = 0; i < P.size(); i++){
        if(P[i]){
            for (int _ = 0; _ < ratio; _ ++){
                if (distribution.size() !=0){
                    int min = 0;
                    for (int j = 0; j < distribution.size(); j++){
                        if (ssd[i][distribution[j]] < ssd[i][distribution[min]]){
                            min = j;
                        }
                    }
                    dp[i].push_back(distribution[min]);
                    distribution.erase(distribution.begin() + min);
                }
            }
        }
    }
    return dp;

}

// returns the vector p
vector<int> distribution_supply(vector<bool> P, vector<bool> D, vector<vector<float> > ssd){
    map<int, vector<int> > dp = distribution_production(P, D, ssd);
    vector<int> ds;
    for (int i = 0; i < P.size(); i++){
        if (P[i]){
            ds.push_back(i);
        }
        if (D[i]){
            for (auto it = dp.begin(); it != dp.end(); it++){
                for (int j = 0; j<it->second.size(); j++){
                    if (i == it->second[j]){
                        ds.push_back(it->first);
                        break;
                    }
                }
            }
        }
        if (!P[i] && !D[i]){
            ds.push_back(-1);
        }
    }
    return ds;
}


// We define the demand for every production center
map<int, int> demand_production (map<int, vector<int> > pd, map<int, vector<int> > dc, vector<pair<int, pair<float, float> > > clients){
    map<int, int> m;
    for (auto i = pd.begin(); i != pd.end(); i++){
        int di = 0;
        for (int j = 0; j < (i->second).size(); j++){
            for (int k = 0; k < dc[i->second[j]].size(); k++){
                di += clients[dc[i->second[j]][k]].first;
            }
        }
        m[i->first] = di;
    }
    return m;
}

// We define the demand left (if there is one) not satsified by a production center
map<int, int> demand_left (map<int, int> dp, int capacity_max){
    map<int, int> dl;
    for (auto it = dp.begin(); it != dp.end(); it++){
        dl[it->first] = max(0, it->second - capacity_max);
    }
    return dl;
}

// We return the vector of the automated sites
vector<bool> automated(map<int, int> dl, vector<bool> P){
    vector<bool> a;
    for (int i = 0; i < P.size(); i++){
        bool ai = 0;
        if (P[i] && dl[i]>0){
            ai = 1;
        }
        a.push_back(ai);
    }
    return a;
}


// Given two numbers k and l, we generate two vectors P and D of k production centers and l distribution centers
vector<bool> production_sites(int k, int s){
    vector<bool> P;
    while (accumulate(P.begin(), P.end(), 0) != k){
        P.clear();
        for (int j = 0; j < s; j++){
            P.push_back(rand()%2);
        }
    }
    return P;
}

vector<bool> distribution_sites(int l, vector<bool> P){
    vector<bool> D;
    while (accumulate(D.begin(), D.end(), 0) != l){
        D.clear();
        for (int j = 0; j < P.size(); j++){
            if (P[j]){
                D.push_back(0);
            }
            else {
                D.push_back(rand()%2);
            }
        }
    }
    return D;
}

void Instance::solve() {
  cout << "Solving instance... \n";
  time_t begin, end;

  begin = clock();

  const int n = sites.size();
  const int m = clients.size();
  // write heuristic solver here
    /**
    for(int i = 0; i < sites.size(); ++i) {
      float seed1 = rand()/(float)INT_MAX;
      float seed2 = rand()/(float)INT_MAX;
      if(seed1 < 0.3) {
        if(seed2 > 0.7)
          solution.P[i] = true;
        else
          solution.D[i] = true;
      }
    }
    **/
    //vector<pair<int, pair<float, float>>> c = {{10, {0, 0}}, {34, {3, 2}}, {12, {-1, 3}}};
    //vector<vector<float>> distances = {vector<float> {1, 3, 2}, vector<float> {3, 4, 2}, vector<float> {6, 3, 8}, vector<float> {2, 4, 3}, vector<float> {5, 2, 3}};

    //vector<vector<float>> distances_sites = {vector<float> {0, 3, 2, 1, 4}, vector<float> {3, 0, 2, 7, 3}, vector<float> {6, 3, 0, 6, 2}, vector<float> {2, 4, 3, 0, 7}, vector<float> {5, 2, 3, 5, 0}};
    //vector<bool> P = {0, 1, 0, 0, 1};
    solution.P = production_sites((int)(0.5*sites.size()), sites.size());
    solution.D = distribution_sites(max(1, (int)(0.3*sites.size())), solution.P);
    map<int, vector<int>> mm = clients_distribution(solution.D, clients, siteClientDistances);
    //show(mm);
    cout << "here i am\n";
    map<int, vector<int>> mp = distribution_production(solution.P, solution.D, siteSiteDistances);
    //show(mp);

    map<int, int> dm = demand_production(mp, mm, clients);
    //show(dm);

    map<int, int> dl = demand_left(dm, cu);
    //show(dl);

    solution.a = automated(dl, solution.P);
    //show(a);

    solution.s = supply_client(solution.D, clients, siteClientDistances);
    //show(s);

    solution.p = distribution_supply(solution.P, solution.D, siteSiteDistances);
    //show(ds);

    //vector<bool> Pi = production_sites(4, 12);
    //show(Pi);

    //vector<bool> Di = distribution_sites(7, Pi);
    //show(Di);

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