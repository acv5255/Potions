#include "potions.hpp"
using std::cout;

int main(int argc, char* argv[]) {
    // Handle command line input
    RunOptions runOptions = RunOptions::FromArguments(argc, argv);

    // Read the model inputs
    ModelInputs modelInputs = ModelInputs::ReadInputs(runOptions.outputName());
    std::cout << "Model folder name: " << runOptions.outputName() << "\n";

    // Get model output directory
    const string outputFilename = getOutputFilepath(runOptions.outputName());
    std::cout << "Writing results to: " << outputFilename << "\n";

    // Construct the data structures for solving the problem
    if (modelInputs.runType() == EQUILIBRIUM) {
        std::cout << "\nBeginning equilibrium solution\n";
        std::cout << "================================\n";
        const ChemicalState chemInit = modelInputs.initChemState();

        std::cout << "Total concentrations:\n";
        for (auto x: modelInputs.chem().primarySpecies()) {
            std::cout << x.first << ": " << x.second << "\n";
        }
        std::cout << std::endl;

        EquilibriumConstants eqConsts = modelInputs.equilibriumConstants();

        // Print out equilibrium parameters
        {
            cout.precision(3);
            cout << "Equilibrium parameters:\n";
            cout << "Stoichiometry matrix\n";
            for (int i = 0; i < eqConsts.stoichMat.n_rows; i++) {
                for (int j = 0; j < eqConsts.stoichMat.n_cols; j++) {
                    cout << eqConsts.stoichMat(i,j) << " ";
                }
                cout << "\n";
            }
            cout << "\n";
            cout << "Log10 of Equilibrium constants\n";
            for (int i = 0; i < eqConsts.eqConsts.size(); i++) {
                cout << eqConsts.eqConsts[i] << "\n";
            }
            cout << "\n";
        }

        TotalConstants totConsts = modelInputs.totalConstants();

        {
            // Print out total constant parameters
            cout << "Total species parameters\n";
            for (int i = 0; i < totConsts.tot_mat.n_rows; i++) {
                for (int j = 0; j < totConsts.tot_mat.n_cols; j++) {
                    cout << totConsts.tot_mat(i,j) << " ";
                }
                cout << "\n";
            }
            cout << "\n";
        }

        const ChemicalState chemFinal = SolveEquilibrium(
            chemInit.totalConcentration,
            eqConsts,
            totConsts
        );

        SaveEquilibriumResults(chemFinal, modelInputs.speciesNames(), outputFilename);
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