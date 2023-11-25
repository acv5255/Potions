# BioRT++

## Compiling the code 
To compile the code on any platform, open a 
command prompt in the root directory and change
into that directory. Then, run the following CMake commands:

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

## Running the model
To run the model, place folder containing the required `cdbs.yaml` and `chem.yaml` in the `input` directory. The format of the input folders are present in the 
documentation folder. Run the model using the following command:

On *nix
```
./potions <input_folder_name>
```
 
On Windows:
```
potions.exe <input_folder_name>
```

## Input Data Format
The input data format follows the standard YAML format and uses YAML maps and arrays, and the structure of these YAML data structures can be found anywhere online.

### Input Data: (chem.yaml)
The 'chem.yaml' file is used to specify the chemical universe, meaning all chemical species, equilibrium reactions, and kinetic reactions. It consists of multiple parts, and *certain* the header sections are required. The required sections are:
- **SimulationType**: (required string) either 'Kinetic' or 'Equiilibrium', specifying the type of simulation to run
- **EndTime** (optional float/int) For a kinetic simulation, this is the number (in seconds) to run the simulation
- **NumSteps**:  (optional integer) The number of points to run the simulation at
- **PrimarySpecies**: (required) a map with the key being a primary species and the value being the total concentration for this species
  - Note that these species must be present in the database
  - Note that H+ must be set to zero because the mass balance on this species is actually a charge balance.
- **SecondarySpecies**: (required) A list of the names of secondary species for equilibrium simulations to run
  - Note that these species must all have entries in the database
- **MineralSpecies**: (required) A map with the key being the name of a mineral (that must exist in the database), and the value being the surface area of the mineral in the system, in terms of m^2
  - Note that this section must have at least one entry for a kinetic simulation but can be empty for an equilibrium simulation, which does not use kinetic species.


### Chemical Database (cdbs.yaml)
The chemical database contains all of the constant emprical data that we can use to set up a simulation. All three sections in this file are required, and the structure follows standard YAML constructs:
- **PrimarySpecies**: An array of all the primary species that can be used in a simulation.
- **Secondary Species**: A map of secondary species entries with the key being the species name and the value having the following entries:
  - **stoichiometry**: A map describing the chemical reaction, with the key being a primary species and the value being the stoichiometric coefficient in the equilibrium reaction
  - **logK**: The base-10 logarithm of the equilibrium constant describing this chemical reaction
- **Mineral Species**: A map containing the information about the kinetic mineral reactions. The key is the mineral name, and the value has the following required values:
  - **stoichiometry**: A map describing the chemical reaction, with the key being a primary species and the value being the stoichiometric coefficient in the equilibrium reaction
  - **logK**: The base-10 logarithm of the equilibrium constant describing this chemical reaction
  - **rate**: The base-10 logarithm of the kinetic rate constant for this kinetic reaction
  - **molarVolume**: the volume of a single mole of this mineral with units cm^3/mol
  - **molarMass**: The mass of a single mole of this mineral with units g/mol