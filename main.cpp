#include "instance.h"
#include <string>
#include <array>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
  srand(time(nullptr));
  array<string, 1> filenames = {"KIRO-large.json"};
  for(string filename : filenames) {
    Instance instance(filename);
    instance.solve();
    instance.save();
    cout << "validity: " << (instance.is_valid() ? "true" : "false") << endl;
    cout << "resulting cost: " << instance.compute_cost()/1e4 << endl;
  }
  return 0;
}
