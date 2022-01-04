#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <utility>
#include <set>
#include <cstdlib>

#include <nlohmann/json.hpp>

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
}

float Instance::dist(int id_client, int id_site) {
    float x = clients[id_client].second.first-sites[id_site].first;
    float y = clients[id_client].second.second-sites[id_site].second;
    float d = x*x + y*y;
    return sqrt(d) + 100/(float)clients[id_client].first;
}

float dist2(pair<float, float>& x, pair<float, float>& y) {
    float a = x.first - y.first;
    float b = x.second - y.second;
    float d = a*a + b*b;
    return sqrt(d);
}

int Instance::plus_proche_site(pair<float, float>& bary) {
    int id = 0;
    for (int i = 1; i < sites.size(); i++){
      if(visited_sites[i]) continue;
      if (dist2(sites[id], bary) > dist2(sites[i], bary)) {
          id = i;
      }
    }
    return id;
}


void Instance::compute_cost() {
  current_cost = 0;
  for(int i = 0; i < sites.size(); ++i) {
    current_cost += solution.P[i]*(capacities[0] + solution.a[i]*capacities[1]);
  }
  for(int i = 0; i < clients.size(); ++i) {
    current_cost += clients[i].first*(productionCosts[0]-solution.a[solution.s[i]]*productionCosts[1]);
    current_cost += clients[i].first*(routingCosts[1]*dist(i, solution.s[i]));
  }
}

void Instance::solve() {
  cout << "Solving instance... \n";
  time_t begin, end;

  begin = clock();
  /**
  solution.P[0] = true;
  solution.D[1] = true;

  solution.a[0] = true;
  solution.p[1] = 0;

  solution.s[0] = 1;
  solution.s[1] = 1;
  solution.s[2] = 0;
  **/

  Solution test(sites.size(), clients.size());

  current_clients.reserve(clients.size());

  for(int i = 0; i < clients.size(); ++i)
    current_clients.push_back(i);
  int number_sites = 0;
  while(!current_clients.empty() && number_sites < sites.size()) {
    pair<float, float> barycenter = barycentre();

    int nearest_site = plus_proche_site(barycenter);

    tri_liste_clients(nearest_site);

    affecte_client_site(nearest_site);
    ++number_sites;
  }



  end = clock();

  cout << "Instance successfully solved in " << (end - begin)/(float)CLOCKS_PER_SEC << " seconds.\n";
}

pair<float, float> Instance::barycentre() {
    int n = current_clients.size();
    pair<float, float> bar(0, 0);
    for(auto it = current_clients.begin(); it != current_clients.end(); ++it) {
      bar.first += clients[*it].second.first;
      bar.second += clients[*it].second.second;
    }
    bar.first /= n;
    bar.second /= n;
    return bar;
}

void Instance::tri_liste_clients(int site_id){
    int n = current_clients.size();
    vector<float> liste_distances;
    for(auto it = current_clients.begin(); it != current_clients.end(); ++it) {
        liste_distances.push_back(dist(*it, site_id));
    }
    for (int i=0; i<n; i++){
        for (int j=0; j < n-i-1; j++){
            if (liste_distances[j+1] < liste_distances[i]){
                swap(current_clients[j], current_clients[j+1]);
            }
        }
    }
}

void Instance::affecte_client_site(int id_site) {
    solution.P[id_site] = true;
    visited_sites[id_site] = true;
    int demande = 0;
    int id = current_clients.back();
    while (demande + clients[id].first <= (capacities[0] + solution.a[id_site]*capacities[1])){
      solution.s[id] = id_site;
      demande += clients[id].first;
      current_clients.pop_back();
      if(current_clients.empty()) break;
      id = current_clients.back();
    }
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
      json_solution["productionCenters"].push_back({{"id", i+1}, {"automation", int(solution.a[i])}});
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
  cout << "Solution successfully saved in " << solution_filename << '\n';
}

Solution::Solution(int sites_number, int clients_number) {
  P = vector<bool>(sites_number, false);
  D = vector<bool>(sites_number, false);
  a = vector<bool>(sites_number, true);
  p = vector<int>(sites_number, -1);
  s = vector<int>(clients_number, -1);
}

