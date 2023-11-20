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
    cout << "Writing results to: " << outputFilename << "\n";

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
        cout << "Beginning kinetic solution\n";
        cout << "================================\n";

        // 1) Construct initial chemical state
        ChemicalState chem = modelInputs.initChemState();
        const map<string, double> surfaceAreaMap = modelInputs.surfaceAreas();
        const map<string, unsigned int> minMap = modelInputs.chem().mineralMap();
        vec surfaceArea = arma::zeros(surfaceAreaMap.size());

        cout << "Mineral surface areas: \n";
        for (auto x: surfaceAreaMap) {
            auto index = minMap.at(x.first);
            surfaceArea[index] = surfaceAreaMap.at(x.first);
            cout << x.first << ": " << surfaceArea[index];
        }
        cout << "\n";
        const TotalConstants totConsts = modelInputs.totalConstants();
        const EquilibriumConstants eqConsts = modelInputs.equilibriumConstants();
        const KineticConstants kinConsts = modelInputs.kineticConstants();

        chem = SolveEquilibrium(chem.totalConcentration, eqConsts, totConsts);

        // 2) Get model time steps
        const int startTime = 0.0;
        const int endTime = modelInputs.chem().endTime();
        const int numSteps = modelInputs.chem().numSteps();
        const double dt = (endTime - startTime) / (double)numSteps;
        const vec timeSteps = arma::linspace(dt, endTime, numSteps);

        // 3) Construct model outputs
        vector<pair<double, ChemicalState>> results(modelInputs.chem().numSteps() + 1);
        results[0] = {0.0, chem};

        // 4) Run the model
        for (int i = 0; i < numSteps; i++) {
            chem = SolveKineticEquilibrium(chem, surfaceArea, kinConsts, eqConsts, totConsts, dt);
            results[i+1] = {timeSteps[i], chem};
        }

        // 5) Save outputs
        SaveKineticResults(results, modelInputs.speciesNames(), outputFilename);

        // 6) Plot outputs
        std::cerr << "ERROR: cannot yet plot model outputs\n";

        // throw NotImplemented();
    }
    else {
        std::cerr << "Error: expected EQUILIBRIUM or KINETIC run types\n";
        exit(-1);
    }

    // Plot the solution
    std::cout << "Successfully ran Potions simulation\n";

    return 0;
}