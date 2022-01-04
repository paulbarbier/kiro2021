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
  visited_sites = vector<bool>(sites.size(), false);
  cout << "Instance " + input_filename + " successfully loaded!\n";
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

  // write MILP solver here using or-tools

  //init solver
  unique_ptr<MPSolver> solver(MPSolver::CreateSolver("SCIP"));

  if (!solver) {
    LOG(WARNING) << "SCIP solver unavailable.";
    return;
  }
  cout << "SCIP solver successfully loaded!\n";

  const size_t n = sites.size();
  const size_t m = clients.size();
  const double infinity = solver->infinity();

  //init variables
  vector<const MPVariable*> P(n), D(n), a(n), t(n);
  vector<vector<const MPVariable*>> X(n, vector<const MPVariable*>(n));
  vector<vector<const MPVariable*>> C_P(n, vector<const MPVariable*>(m));
  vector<vector<const MPVariable*>> C_D(n, vector<const MPVariable*>(m));
  vector<vector<const MPVariable*>> I_P(n, vector<const MPVariable*>(m));
  vector<vector<const MPVariable*>> C_Pa(n, vector<const MPVariable*>(m));
  vector<vector<const MPVariable*>> I_Pa(n, vector<const MPVariable*>(m));
  
  for(auto& var : P)
    var = solver->MakeBoolVar("");
  for(auto& var : D)
    var = solver->MakeBoolVar("");
  for(auto& var : a)
    var = solver->MakeBoolVar("");
  for(auto& var : t)
    var = solver->MakeNumVar(0.0, infinity, "");

  for(auto& row : X)
    for(auto& var : row)
      var = solver->MakeBoolVar("");
  for(auto& row : C_P)
    for(auto& var : row)
      var = solver->MakeBoolVar("");
  for(auto& row : C_D)
    for(auto& var : row)
      var = solver->MakeBoolVar("");
  for(auto& row : I_P)
    for(auto& var : row)
      var = solver->MakeBoolVar("");
  for(auto& row : C_Pa)
    for(auto& var : row)
      var = solver->MakeBoolVar("");
  for(auto& row : I_Pa)
    for(auto& var : row)
      var = solver->MakeBoolVar("");

  // define the constraints



  //define the objective function
  MPObjective* const objective = solver->MutableObjective();
  /**
  for (int j = 0; j < data.num_vars; ++j) {
    objective->SetCoefficient(x[j], data.obj_coeffs[j]);
  }
  **/
  objective->SetMinimization();
  
  //solve the minimization problem
  const MPSolver::ResultStatus result_status = solver->Solve();

  
  //retrieve the solution and put the data in solution


  //end of MILP implementation
  //
  //

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
  P = vector<bool>(sites_number);
  D = vector<bool>(sites_number);
  a = vector<bool>(sites_number);
  p = vector<int>(sites_number);
  s = vector<int>(clients_number);
}

