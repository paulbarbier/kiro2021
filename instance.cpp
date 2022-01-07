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

    input_file.close();
  }
  solution = Solution(sites.size(), clients.size());
  cout << "Instance " + input_filename + " successfully loaded!\n\n";
  
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

  cout << "\n\n";
}

void Instance::solve() {
  cout << "Solving instance... \n";
  time_t begin, end;

  begin = clock();

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
  vector<vector<vector<const MPVariable*>>> v(n, vector<vector<const MPVariable*>>(m, vector<const MPVariable*>(n)));
  
  for(auto& var : P)
    var = solver->MakeBoolVar("");
  for(auto& var : D)
    var = solver->MakeBoolVar("");
  for(auto& var : a)
    var = solver->MakeBoolVar("");
  for(auto& var : t)
    var = solver->MakeNumVar(0.0, infinity, "");

  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      X[i][j] = solver->MakeBoolVar("");

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
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < m; ++j)
      for(int k = 0; k < n; ++k)
        v[i][j][k] = solver->MakeBoolVar("");

  // define the constraints
  /**
  // X_ii = 0
  for(int i = 0; i < n; ++i) {
    MPConstraint* constraint = solver->MakeRowConstraint(0.0, 0.0, "");
    constraint->SetCoefficient(X[i][i], 1);
  }

  // v_iji = 0
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      MPConstraint* constraint = solver->MakeRowConstraint(0.0, 0.0, "");
      constraint->SetCoefficient(v[i][j][i], 1);
    }
  }
  **/
  // P_i + D_i <= 1
  for(int i = 0; i < n; ++i) {
    MPConstraint* constraint = solver->MakeRowConstraint(-infinity, 1.0, "");
    constraint->SetCoefficient(P[i], 1);
    constraint->SetCoefficient(D[i], 1);
  }
  
  // \sum_i X_ij = D_j, i \neq j
  for(int j = 0; j < n; ++j) {
    MPConstraint* constraint = solver->MakeRowConstraint(0.0, 0.0, "");
    for(int i = 0; i < n; ++i) {
      if(i != j)
        constraint->SetCoefficient(X[i][j], 1);
    }
    constraint->SetCoefficient(D[j], -1);
  }

  // \sum_j C_P_ij + C_D_ij = 1
  for(int j = 0; j < m; ++j) {
    MPConstraint* constraint = solver->MakeRowConstraint(1.0, 1.0, "");
    for(int i = 0; i < n; ++i) {
      constraint->SetCoefficient(C_P[i][j], 1);
      constraint->SetCoefficient(C_D[i][j], 1);
    }
  }

  // a_i <= P_i
  for(int i = 0; i < n; ++i) {
    MPConstraint* constraint = solver->MakeRowConstraint(-infinity, 0.0, "");
    constraint->SetCoefficient(a[i], 1);
    constraint->SetCoefficient(P[i], -1);
  }
  
  // X_ij <= P_i, i \neq j
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(i != j) {
        MPConstraint* constraint = solver->MakeRowConstraint(-1.0, 0.0, "");
        constraint->SetCoefficient(X[i][j], 1);
        constraint->SetCoefficient(P[i], -1);
      }
    }
  }
  
  // C_P_ij <= P_i
  // C_D_ij <= D_i
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      MPConstraint* constraint1 = solver->MakeRowConstraint(-1.0, 0.0, "");
      constraint1->SetCoefficient(C_P[i][j], 1);
      constraint1->SetCoefficient(P[i], -1);

      MPConstraint* constraint2 = solver->MakeRowConstraint(-1.0, 0.0, "");
      constraint2->SetCoefficient(C_D[i][j], 1);
      constraint2->SetCoefficient(D[i], -1);
    }
  }
  /**
  // for testing
  for(int i = 0; i < n; ++i) {
    MPConstraint* c1 = solver->MakeRowConstraint(-infinity, 0.0, "");
    MPConstraint* c2 = solver->MakeRowConstraint(0.0, infinity, "");
    for(int j = 0; j < m; ++j) {
      c1->SetCoefficient(C_P[i][j], 1/(double)m);
      c2->SetCoefficient(C_P[i][j], 1);
    }
    for(int j = 0; j < n; ++j) {
      if(i != j) {
        c1->SetCoefficient(X[i][j], 1/(double)n);
        c2->SetCoefficient(X[i][j], 1);
      }
    }
    c1->SetCoefficient(P[i], -1);
    c2->SetCoefficient(P[i], -1);
  }

  
  for(int i = 0; i < n; ++i) {
    MPConstraint* c1 = solver->MakeRowConstraint(-infinity, 0.0, "");
    MPConstraint* c2 = solver->MakeRowConstraint(0.0, infinity, "");
    for(int j = 0; j < m; ++j) {
      c1->SetCoefficient(C_D[i][j], 1/(double)m);
      c2->SetCoefficient(C_D[i][j], 1);
    }
    c1->SetCoefficient(D[i], -1);
    c2->SetCoefficient(D[i], -1);
  }
  //
  **/
  
  // I_P_ij = \sum_k v_ijk, k \neq i
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      MPConstraint* constraint = solver->MakeRowConstraint(0.0, 0.0, "");
      for(int k = 0; k < n; ++k) {
        if(k != i)
          constraint->SetCoefficient(v[i][j][k], 1);
      }
      constraint->SetCoefficient(I_P[i][j], -1);
    }
  }
  
  /**
   * implemented constraints : v_ijk = X_ik*C_D_kj with :
   * (X_ik + C_D_kj - v_ijk) <= 1
   * v_ijk <= 0.5*(X_ik + C_D_kj)
   * i \neq k
   **/
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      for(int k = 0; k < n; ++k) {
        if(i != k) {
          MPConstraint* constraint1 = solver->MakeRowConstraint(-infinity, 1.0, "");
          MPConstraint* constraint2 = solver->MakeRowConstraint(0.0, infinity, "");

          constraint1->SetCoefficient(C_D[k][j], 1);
          constraint1->SetCoefficient(X[i][k], 1);
          constraint1->SetCoefficient(v[i][j][k], -1);

          constraint2->SetCoefficient(C_D[k][j], 0.5);
          constraint2->SetCoefficient(X[i][k], 0.5);
          constraint2->SetCoefficient(v[i][j][k], -1);
        }
      }
    }
  }

  // a_i + C_P_ij - C_Pa_ij <= 1
  // 0 <= 0.5*(a_i + C_P_ij) - C_Pa_ij
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      MPConstraint* constraint1 = solver->MakeRowConstraint(-infinity, 1.0, "");
      MPConstraint* constraint2 = solver->MakeRowConstraint(0.0, infinity, "");
      
      constraint1->SetCoefficient(a[i], 1);
      constraint1->SetCoefficient(C_P[i][j], 1);
      constraint1->SetCoefficient(C_Pa[i][j], -1);

      constraint2->SetCoefficient(a[i], 0.5);
      constraint2->SetCoefficient(C_P[i][j], 0.5);
      constraint2->SetCoefficient(C_Pa[i][j], -1);
    }
  }

  // a_i + I_P_ij - I_Pa_ij <= 1
  // 0 <= 0.5*(a_i + I_P_ij) - I_Pa_ij
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      MPConstraint* constraint1 = solver->MakeRowConstraint(-infinity, 1.0, "");
      MPConstraint* constraint2 = solver->MakeRowConstraint(0.0, infinity, "");
      
      constraint1->SetCoefficient(a[i], 1);
      constraint1->SetCoefficient(I_P[i][j], 1);
      constraint1->SetCoefficient(I_Pa[i][j], -1);
      
      /**
      for(int k = 0; k < n; ++k) {
        if(i != k)
          constraint1->SetCoefficient(v[i][j][k], 1);
      }
      **/

      constraint2->SetCoefficient(a[i], 0.5);
      constraint2->SetCoefficient(I_P[i][j], 0.5);
      constraint2->SetCoefficient(I_Pa[i][j], -1);
      
      /**
      for(int k = 0; k < n; ++k) {
        if(i != k)
          constraint2->SetCoefficient(v[i][j][k], 0.5);
      }
      **/
    }
  }
  for(int i = 0; i < n; ++i) {
    MPConstraint* constraint = solver->MakeRowConstraint(-infinity, cu*u_P, "");

    constraint->SetCoefficient(t[i], -1);
    constraint->SetCoefficient(a[i], -cu*u_A);
    for(int j = 0; j < m; ++j) {
      constraint->SetCoefficient(C_P[i][j], cu*clients[j].first);
      constraint->SetCoefficient(I_P[i][j], cu*clients[j].first);
      /**
      for(int k = 0; k < n; ++k) {
        if(i != k)
          constraint->SetCoefficient(v[i][j][k], cu*clients[j].first);
      }
      **/
    }
  }
  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < m; ++j) {
      MPConstraint* constraint = solver->MakeRowConstraint(0.0, 1.0, "");
      constraint->SetCoefficient(C_P[i][j], 1);
      constraint->SetCoefficient(I_P[i][j], 1);
      /**
      for(int l = 0; l < n; ++l) {
        if(i != l)
          constraint->SetCoefficient(v[i][j][l], 1);
      }
      **/
    }
  }
  
  
  //define the objective function
  MPObjective* const objective = solver->MutableObjective();
  
  // building cost
  for(int i = 0; i < n; ++i) {
    objective->SetCoefficient(P[i], c_Pb);
    objective->SetCoefficient(D[i], c_Db);
    objective->SetCoefficient(a[i], c_Ab);
  }
  
  // production cost
  for(int j = 0; j < m; ++j) {
    for(int i = 0; i < n; ++i) {
      // here we put also the coefficient in front of C_D in routing cost
      objective->SetCoefficient(C_D[i][j], clients[j].first*(c_Dp + siteClientDistances[i][j]*c_2r));
      objective->SetCoefficient(C_Pa[i][j], -clients[j].first*c_Ap);
      objective->SetCoefficient(I_Pa[i][j], -clients[j].first*c_Ap);
    }
  }

  //routing cost
  for(int j = 0; j < m; ++j) {
    for(int i = 0; i < n; ++i) {
      objective->SetCoefficient(C_P[i][j], clients[j].first*siteClientDistances[i][j]*c_2r);
      for(int k = 0; k < n; ++k) {
        /**
        for(int l = 0; l < n; ++l) {
          if(k != l)
            objective->SetCoefficient(v[k][j][l], clients[j].first*siteSiteDistances[k][i]*c_1r);
        }
        **/
       //objective->SetCoefficient(I_P[k][j], clients[j].first*siteSiteDistances[k][i]*c_1r);
        
        // here I put X_ki instead of I_P_kj
        if(k != i)
          objective->SetCoefficient(X[k][i], clients[j].first*siteSiteDistances[k][i]*c_1r);
      }
    }
  }
  
  // capacity cost
  for(int i = 0; i < n; ++i) {
    objective->SetCoefficient(t[i], 1);
  }

  objective->SetMinimization();
  //additional constant cost to get the right optimal cost value
  int constant_cost = 0;
  for(int j = 0; j < m; ++j)
    constant_cost += clients[j].first;
  constant_cost *= c_Pp;
  
  //solve the minimization problem
  const MPSolver::ResultStatus result_status = solver->Solve();
  if (result_status != MPSolver::OPTIMAL) {
    LOG(FATAL) << "The problem does not have an optimal solution.";
  }
  LOG(INFO) << "Optimal objective value = " << (objective->Value() + constant_cost)/1e4;
  /**
  cout << "C_P:\n";
  for(auto& row : C_P) {
    for(auto& var : row)
      cout << int(var->solution_value()) << " ";
    cout << endl;
  }
  cout << "\n\n";
  cout << "C_D:\n";
  for(auto& row : C_D) {
    for(auto& var : row)
      cout << int(var->solution_value()) << " ";
    cout << endl;
  }
  cout << "X:\n";
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j)
      cout << int(X[i][j]->solution_value()) << " ";
    cout << endl;
  }
  cout << "\nP D a\n";
  **/
  //retrieve the solution and put the data in solution
  for(int i = 0; i < n; ++i) {
    //cout << int(P[i]->solution_value()) << " " << int(D[i]->solution_value()) << " " << int(a[i]->solution_value()) << endl;
    solution.P[i] = (int(P[i]->solution_value()) == 1);
    solution.D[i] = (int(D[i]->solution_value()) == 1);
    solution.a[i] = (int(a[i]->solution_value()) == 1);

    if(solution.D[i] == false) {
      solution.p[i] = -1;
    } else {
      for(int j = 0; j < n; ++j) {
        if( j != i && int(X[j][i]->solution_value()) == 1) {
          solution.p[i] = j;
          break;
        }
      }
    }
  }
  for(int j = 0; j < m; ++j) {
    for(int i = 0; i < n; ++i) {
      if(int(C_P[i][j]->solution_value()) + int(C_D[i][j]->solution_value()) == 1) {
        solution.s[j] = i;
        break;
      }
    }
  }
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

