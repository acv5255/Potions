#include "potions.hpp"

int main(int argc, char* argv[]) {
    // Handle command line input
    RunOptions runOptions = RunOptions::FromArguments(argc, argv);

    // Read the model inputs
    ModelInputs modelInputs = ModelInputs::ReadInputs(runOptions.outputName());

    // Get model output directory
    const path outputFilename = getOutputFilepath(runOptions.outputName());

    // Construct the data structures for solving the problem
    if (modelInputs.runType() == EQUILIBRIUM) {
        std::cout << "Beginning equilibrium solution\n";
        EquilibriumSolution solver = EquilibriumSolution::Initialize(modelInputs);
        auto primConcInit = solver.primConc();
        const auto[primConc, conc] = solver.solve(primConcInit);
        solver.updatePrimConc(primConc);
        solver.updateConcentration(conc);

        // Now, write the outputs
        solver.saveSolution(outputFilename);

        std::cout << "Successfully computed equilibrium solution\n";
        std::cout << "Wrote solution to: '" << outputFilename.string() << "'\n";
    }
    else if (modelInputs.runType() == KINETIC)
    {
        std::cout << "Beginning kinetic solution\n";

        throw NotImplemented();
    }
    else {
        std::cerr << "Error: expected EQUILIBRIUM or KINETIC run types\n";
        exit(-1);
    }

    // Plot the solution
    std::cout << "Successfully ran Potions simulation\n";

    return 0;
}