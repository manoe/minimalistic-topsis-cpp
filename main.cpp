#include "topsis.h"
#include <iostream>

using namespace arma;

int main() {
    TopsisEngine te(2,2);
    te.addAlternative({1,0.0},{2,3});
    te.addAlternative({2, 0.0},{3,1});
    for (auto i : te.getRanking()) {
        std::cout<<i.id<<" "<<i.rank<<std::endl;
    }
    return 0;
}
