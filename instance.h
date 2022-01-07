#pragma once

#include <string>
#include <vector>
#include <utility>

using namespace std;

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

  float c_Pb, c_Ab, c_Db, c_Pp, c_Ap, c_Dp, c_1r, c_2r, u_P, u_A, cu;
  
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
