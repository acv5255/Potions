#include <chrono>
#include "potions.hpp"
using std::cout;
using namespace std::chrono;

int main(int argc, char* argv[]) {
    // Handle command line input
    RunOptions runOptions = RunOptions::from_arguments(argc, argv);

    // Read the model inputs
    ModelInputs modelInputs = ModelInputs::read_inputs(runOptions.output_name());
    std::cout << "Model folder name: " << runOptions.output_name() << "\n";

    // Get model output directory
    const string outputFilename = get_output_file_path(runOptions.output_name());
    cout << "Writing results to: " << outputFilename << "\n";

    auto start = high_resolution_clock::now();

    // Construct the data structures for solving the problem
    if (modelInputs.run_type() == EQUILIBRIUM) {
        std::cout << "\nBeginning equilibrium solution\n";
        std::cout << "================================\n";
        const ChemicalState chemInit = modelInputs.initial_chem_state();

        std::cout << "Total concentrations:\n";
        for (auto x: modelInputs.chem().primary_species()) {
            std::cout << x.first << ": " << x.second << "\n";
        }
        std::cout << std::endl;

        EquilibriumConstants eqConsts = modelInputs.equilibrium_constants();

        // Print out equilibrium parameters
        {
            cout.precision(3);
            cout << "Equilibrium parameters:\n";
            cout << "Stoichiometry matrix\n";
            for (int i = 0; i < eqConsts.stoichiometry_matrix.n_rows; i++) {
                for (int j = 0; j < eqConsts.stoichiometry_matrix.n_cols; j++) {
                    cout << eqConsts.stoichiometry_matrix(i,j) << " ";
                }
                cout << "\n";
            }
            cout << "\n";
            cout << "Log10 of Equilibrium constants\n";
            for (int i = 0; i < eqConsts.log_eq_consts.size(); i++) {
                cout << eqConsts.log_eq_consts[i] << "\n";
            }
            cout << "\n";
        }

        TotalConstants totConsts = modelInputs.total_constants();

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

        const ChemicalState chemFinal = solve_equilibrium(
            chemInit.total_concentration,
            eqConsts,
            totConsts
        );

        auto stop = high_resolution_clock::now();
        auto duration = (long)duration_cast<milliseconds>(stop - start).count();
        double seconds = (double)duration / 1000.0;
        cout << "Time for equilibrium simulation: " << seconds << " seconds\n";

        save_equilibrium_results(chemFinal, modelInputs.species_names(), outputFilename);
    }
    else if (modelInputs.run_type() == KINETIC)
    {
        cout << "Beginning kinetic solution\n";
        cout << "================================\n";

        // 1) Construct initial chemical state
        ChemicalState chem = modelInputs.initial_chem_state();
        const map<string, double> surfaceAreaMap = modelInputs.surface_areas();
        const map<string, unsigned int> minMap = modelInputs.chem().mineral_map();
        vec surfaceArea = arma::zeros(surfaceAreaMap.size());

        cout << "Mineral surface areas: \n";
        for (auto x: surfaceAreaMap) {
            auto index = minMap.at(x.first);
            surfaceArea[index] = surfaceAreaMap.at(x.first);
            cout << x.first << ": " << surfaceArea[index];
        }
        cout << "\n";
        const TotalConstants totConsts = modelInputs.total_constants();
        const EquilibriumConstants eqConsts = modelInputs.equilibrium_constants();
        const KineticConstants kinConsts = modelInputs.kinetic_constants();

        chem = solve_equilibrium(chem.total_concentration, eqConsts, totConsts);

        // 2) Get model time steps
        const int startTime = 0.0;
        const int endTime = modelInputs.chem().end_time();
        const int numSteps = modelInputs.chem().num_steps();
        const double dt = (endTime - startTime) / (double)numSteps;
        const vec timeSteps = arma::linspace(dt, endTime, numSteps);

        // 3) Construct model outputs
        vector<pair<double, ChemicalState>> results(modelInputs.chem().num_steps() + 1);
        results[0] = {0.0, chem};

        // 4) Run the model
        for (int i = 0; i < numSteps; i++) {
            chem = solve_kinetic_equilibrium(chem, surfaceArea, kinConsts, eqConsts, totConsts, dt);
            results[i+1] = {timeSteps[i], chem};
        }
        auto stop = high_resolution_clock::now();
        auto duration = (long)duration_cast<milliseconds>(stop - start).count();
        double seconds = (double)duration / 1000.0;
        cout << "Time for kinetic simulation: " << seconds << " seconds\n";

        // 5) Save outputs
        save_kinetic_results(results, modelInputs.species_names(), outputFilename);

        // 6) Plot outputs
        // std::cerr << "ERROR: cannot yet plot model outputs\n";
        plot_results(results, modelInputs.species_names());

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