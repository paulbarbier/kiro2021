#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <utility>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <queue>

#include <nlohmann/json.hpp>
#include <ortools/linear_solver/linear_solver.h>

#include "instance.h"

using json = nlohmann::json;
using namespace std;
using namespace operations_research;

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

  cout << "clients: " << clients.size() << " sites: " << sites.size() << endl;
}

float dist2(pair<float, float>& x, pair<float, float>& y) {
    float a = x.first - y.first;
    float b = x.second - y.second;
    float d = a*a + b*b;
    return d;
}

void Instance::solve() {
  cout << "Solving instance... \n";
  time_t begin, end;

  begin = clock();

  const int n = sites.size();
  const int m = clients.size();
  // write heuristic solver here
  
  //step 0: compute the demand
  float demand = 0;
  for(auto& client : clients)
    demand += client.first;

  int min_production_center;
  cout << "minimum production centers needed: " << (min_production_center = (int)(demand/(u_P + u_A)) + 2) << endl;


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
    auto compare_client = [*this, i] (int x, int y) {
      return this->siteClientDistances[i][x] < this->siteClientDistances[i][y];
    };
    sort(order_siteclient[i].begin(), order_siteclient[i].end(), compare_client);
  }

  //step 3: assign clients and distribution sites
  int max_siteclients = 500;
  int max_sitesites = 40;

  //vector<bool> visited_clients(m, false);
  //vector<bool> visited_sites(n, false);
  vector<vector<int>> X_result(n);
  vector<float> tracked_obj(n, 0);

  float alpha = 0.95;
  float beta = 1.0;
  float _alpha = alpha*(float)(u_P + u_A);
  float _beta = beta*(float)(u_P + u_A);

  int n_visited_clients = 0;
  int n_visited_distribution = 0;

  while(n_visited_clients < m) {
    if(demand > (u_P + u_A)) {
      vector<float> obj(n, numeric_limits<float>::infinity());
      for(int i = 0; i < n; ++i) {
        if(solution.D[i] || solution.P[i]) continue;
        
        unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));
        if (!solver) {
          LOG(WARNING) << "SCIP solver unavailable.";
          return;
        }
        //cout << "SCIP solver successfully loaded!\n";
        vector<const MPVariable*> X(m, nullptr);

        MPConstraint* constraint = solver->MakeRowConstraint(_alpha*min( 1.0, 70.0/(n_visited_clients+1) ), _beta, "");
        MPObjective* const objective = solver->MutableObjective();

        for(int j = 0; j < max_siteclients; ++j) {
          if(solution.s[order_siteclient[i][j]] >= 0) continue;
          X[order_siteclient[i][j]] = solver->MakeBoolVar("");
          constraint->SetCoefficient(X[order_siteclient[i][j]], clients[order_siteclient[i][j]].first);
          objective->SetCoefficient(X[order_siteclient[i][j]], w[i][order_siteclient[i][j]]);
        }

        objective->SetMinimization();
        const MPSolver::ResultStatus result_status = solver->Solve();
        obj[i] = objective->Value();
        //cout << "obj:" << (obj[i] = objective->Value()) << "\n";
        for(int k = 0; k < m; ++k) {
          if(X[k] != nullptr && int(X[k]->solution_value()) == 1) {
            X_result[i].push_back(k);
          }
        }

        if (result_status != MPSolver::OPTIMAL) {
          LOG(FATAL) << "The problem does not have an optimal solution.";
        }
        //cout << i << " site processed.\n";
        solver.reset();
      }
      int i_0 = min_element(obj.begin(), obj.end()) - obj.begin();
      if(solution.D[i_0]) continue;
      solution.D[i_0] = true;
      ++n_visited_distribution;
      for(auto& j : X_result[i_0]) {
        if(solution.s[j] >= 0) continue;
        solution.s[j] = i_0;
        tracked_obj[i_0] += clients[j].first;
        ++n_visited_clients;
        demand -= clients[j].first;
      }

      if(tracked_obj[i_0] > 0.2*(u_P+u_A)) {
        solution.D[i_0] = false;
        solution.P[i_0] = true;
        solution.a[i_0] = !(tracked_obj[i_0] <= u_P);
      }

      std::cout << m-n_visited_clients << " clients remaining.\n";

    } else { // we compute the barycenter of the remaining clients to be supplied
      pair<float, float> barycenter({0, 0});
      for(int i = 0; i < m; ++i) {
        if(solution.s[i] >= 0) continue;

        barycenter.first += clients[i].second.first;
        barycenter.second += clients[i].second.second;
      }
      barycenter.first /= (float)(m-n_visited_clients);
      barycenter.second /= (float)(m-n_visited_clients);

      int id = 0;
      for(; id < n; ++id) {
        if(!solution.D[id] && !solution.P[id]) {
          break;
        }
      }

      for (int i = id+1; i < n; i++){
        if(solution.D[i] || solution.P[i]) continue;
        if (dist2(sites[id], barycenter) > dist2(sites[i], barycenter)) {
          id = i;
        }
      }

      solution.P[id] = true;
      solution.a[id] = !(demand <= u_P);
      //++n_visited_distribution;
      for(int i = 0; i < m; ++i) {
        if(solution.s[i] >= 0) continue;
        solution.s[i] = id;
        //tracked_obj[id] += clients[i].first;
        ++n_visited_clients;
      }
    }
  }
  std::cout << "number of distribution sites: " << n_visited_distribution << endl;

  //step 5: assign production sites and bind them with distribution center
    //preliminary: ordering the sites by site distance
  vector<vector<int>> order_sitesite(n);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(i != j) {
        order_sitesite[i].push_back(j);
      }
    }
    auto compare_site = [*this, i] (int x, int y) {
      return this->siteSiteDistances[i][x] < this->siteSiteDistances[i][y];
    };
    sort(order_sitesite[i].begin(), order_sitesite[i].end(), compare_site);
  }

  // step 5: compute barycenter and assign production sites to distribution site
  unordered_set<int> distribution_centers;
  for(int i = 0; i < n; ++i) {
    if(solution.D[i]) distribution_centers.insert(i);
  }

  while(!distribution_centers.empty()) {
    int first = *distribution_centers.begin();
    distribution_centers.erase(first);

    vector<int> neighbours;
    neighbours.push_back(first);
    
    pair<float, float> barycenter(sites[first]);
    float capacity = tracked_obj[first];

    for(int j = 0; (j < n-1) && (capacity <= (u_P + u_A)); ++j) {
      if(distribution_centers.find(order_sitesite[first][j]) != distribution_centers.end()) {
        neighbours.push_back(order_sitesite[first][j]);
        distribution_centers.erase(order_sitesite[first][j]);

        barycenter.first += sites[order_sitesite[first][j]].first;
        barycenter.second += sites[order_sitesite[first][j]].second;
        capacity += tracked_obj[order_sitesite[first][j]];
      }
    }
    barycenter.first /= (float)(neighbours.size());
    barycenter.second /= (float)(neighbours.size());

    int production_id = 0;
    for(; production_id < n; ++production_id) {
      if(!solution.P[production_id] && !solution.D[production_id])
        break;
    }

    for(int k = production_id+1; k < n; ++k) {
      if(solution.P[k] || solution.D[k]) continue;
      if(dist2(sites[production_id], barycenter) > dist2(sites[k], barycenter)) {
        production_id = k;
      }
    }

    solution.P[production_id] = true;
    solution.a[production_id] = !(capacity <= u_P);

    for(auto& distribution : neighbours) {
      solution.p[distribution] = production_id;
    }
  }
  /**
  //vector<bool> visited_distribution(n, false);
  while(n_visited_distribution > 0) {
    int first = 0;
    for(; first < n; ++first) {
      if(solution.D[first])
        break;
    }
    //visited_distribution[first] = true;
    vector<int> neighbours;
    neighbours.push_back(first);
    pair<float, float> barycenter(sites[first]);

    float capacity = tracked_obj[first];

    for(int j = 0; (j < n-1) && (capacity <= u_P + u_A); ++j) {
      if(!visited_distribution[order_sitesite[first][j]] && !solution.P[order_sitesite[first][j]] && solution.D[order_sitesite[first][j]]) {
        visited_distribution[order_sitesite[first][j]] = true;
        neighbours.push_back(order_sitesite[first][j]);

        barycenter.first += sites[order_sitesite[first][j]].first;
        barycenter.second += sites[order_sitesite[first][j]].second;

        capacity += tracked_obj[order_sitesite[first][j]];
      }
    }
    barycenter.first /= (float)(neighbours.size());
    barycenter.second /= (float)(neighbours.size());

    int production_id = 0;
    for(; production_id < n; ++production_id) {
      if(!solution.P[production_id] && !solution.D[production_id])
        break;
    }

    for(int k = production_id+1; k < n; ++k){
      if(solution.P[k] || solution.D[k]) continue;
      if(dist2(sites[production_id], barycenter) > dist2(sites[k], barycenter)) {
        production_id = k;
      }
    }

    solution.P[production_id] = true;
    solution.a[production_id] = !(capacity <= u_P);

    for(auto& distribution : neighbours) {
      solution.p[distribution] = production_id;
      --n_visited_distribution;
    }
  }
  **/
  end = clock();

  std::cout << "Instance successfully solved in " << (end - begin)/(float)CLOCKS_PER_SEC << " seconds.\n";
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
      json_solution["productionCenters"].push_back({{"id", i+1}, {"automation", (solution.a[i] ? 1 : 0)}});
  }

  for(int i = 0; i < solution.D.size(); ++i) {
    if(solution.D[i])
      json_solution["distributionCenters"].push_back({{"id", i+1}, {"parent", (solution.p[i]+1)}});
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
  p = vector<int>(sites_number, -1);
  s = vector<int>(clients_number, -1);
}

// Compute the cost of a solution
float Instance::compute_cost(){
    float building_cost = 0;
    float capacity_cost = 0;
    for (int i = 0; i < solution.P.size(); i++){
        if (solution.P[i]){
            building_cost += c_Pb + solution.a[i]*c_Ab;
            int demand = 0;
            for (int j = 0; j < clients.size(); j++){
                if (solution.s[j] == i || solution.p[solution.s[j]] == i){
                    demand += clients[j].first;
                }
            }
            float t = cu*(demand - u_P -  solution.a[i]*u_A);
            if (t > 0){
                capacity_cost += t;
            }
        }
        if (solution.D[i]){
            building_cost += c_Db;
        }
    }

    float production_cost = 0;
    float routing_cost = 0;
    for (int i = 0; i < clients.size(); i++){
        if (solution.P[solution.s[i]]){
            production_cost += clients[i].first * (c_Pp - solution.a[solution.s[i]]*c_Ap);
            routing_cost += clients[i].first * c_2r * siteClientDistances[solution.s[i]][i];
        }
        if (solution.D[solution.s[i]]){
            production_cost += clients[i].first * (c_Pb - solution.a[solution.p[solution.s[i]]]*c_Ap + c_Dp);
            routing_cost += clients[i].first * (c_1r * siteSiteDistances[solution.p[solution.s[i]]][solution.s[i]] + c_2r * siteClientDistances[solution.s[i]][i]);
        }
    }

    return building_cost + capacity_cost + production_cost + routing_cost;

}

// Check une solution
bool Instance::is_valid(){
    for (int i = 0; i < solution.P.size(); i++){
        if (solution.P[i] + solution.D[i] == 2){
            return false;
        }
        if (solution.a[i] > solution.P[i]){
            return false;
        }
    }
    for (int i = 0; i < solution.p.size(); i++){
        if (solution.P[solution.p[i]] == 0){
            return false;
        }
    }

    for (int i = 0; i < solution.s.size(); i++){
        if (solution.P[solution.s[i]] + solution.D[solution.s[i]] != 1){
            return false;
        }
    }

    return true;
}