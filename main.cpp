#include "instance.h"
#include <string>
#include <array>

using namespace std;

int main(int argc, char** argv) {
  srand(time(nullptr));
  array<string, 1> filenames = {"KIRO-large.json"};
  for(string filename : filenames) {
    Instance instance(filename);
    instance.solve();
    instance.save();
  }
  return 0;
}
