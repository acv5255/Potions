#include "potions.hpp"

PotionsRunType ModelInputs::runType() {
    return _runType;
}

ModelInputs ModelInputs::ReadInputs(const string& inputName) {
    throw NotImplemented();
}