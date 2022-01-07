#include "instance.h"
#include <string>
#include <array>

using namespace std;

int main(int argc, char** argv) {
  srand(time(nullptr));
  array<string, 2> filenames = {"KIRO-tiny.json", "KIRO-small.json"};
  for(string filename : filenames) {
    Instance instance(filename);
    instance.solve();
    instance.save();
  }
  return 0;
}
