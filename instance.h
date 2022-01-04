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
  
  int current_cost;
  vector<int> current_clients;
  vector<bool> visited_sites;

  Solution solution;

  string input_filename;
  public:
    pair<float, float> barycentre();
    Instance(const string& filename);
    void tri_liste_clients(int site_id);
    void affecte_client_site(int id_site);
    void compute_cost();
    int plus_proche_site(pair<float, float>& bary);
    float dist(int id_client, int id_site);
    void solve();
    void save();
};
