# BioRT++

## Compiling the code 
To compile the code on any platform, open a 
command prompt in the root directory and change
into that directory. Then, run the following CMake commands:

```
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