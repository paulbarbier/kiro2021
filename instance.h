#pragma once

#include <string>
#include <vector>
#include <utility>
#include "ortools/linear_solver/linear_solver.h"

using namespace std;
using namespace operations_research;

struct Solution{
    vector<bool> P;
    vector<bool> D;
    vector<bool> a;
    vector<int> p;
    vector<int> s;
    Solution(){};
    Solution(int sites_number, int clients_number);
};

class Instance {
  vector<float> buildingCosts;
  vector<float> productionCosts;
  vector<float> routingCosts;
  float capacityCost;
  vector<float> capacities;
  vector<pair<int, pair<float, float>>> clients;
  vector<pair<float, float>> sites;
  vector<vector<float>> siteSiteDistances;
  vector<vector<float>> siteClientDistances;
  
  vector<int> current_clients;
  vector<bool> visited_sites;

  Solution solution;

  string input_filename;
  public:
    pair<float, float> barycentre();
    Instance(const string& filename);
    void solve();
    void save();
};
