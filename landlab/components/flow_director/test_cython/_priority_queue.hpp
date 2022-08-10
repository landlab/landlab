// priority_queue
#include <queue>
// pair
#include <utility>
// function
#include <functional>
using namespace std;

using _priority_queue = std::priority_queue<std::pair<long,double>,std::vector<std::pair<long,double>>,std::function<bool(std::pair<long,double>,std::pair<long,double>)>>;