#include "topsis.h"
#include <iostream>

using namespace arma;

int main() {
    TopsisEngine te(3,2);
    te.addAlternative({1,0.0},{0.125,0.125}); // a
    te.addAlternative({2, 0.0},{0.75,0.5});  // b
    te.addAlternative({3, 0.0},{0.5,0.75});  // c
    te.addBenefits({false, true});
    te.addWeights({6,8});
    for (auto [id, rank] : te.getRanking()) {
        std::cout<<id<<" "<<rank<<std::endl;
    }
    return 0;
}
