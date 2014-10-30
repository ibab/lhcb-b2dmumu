
# Decfiles

Decfiles are used to configure MC generation with Gauss.
The official LHCb decfiles (2012) can be retrieved with

    SetupProject Gauss v45r5 --build-env
    getpack Gen/DecFiles head

Each decfile has to receive a unique numeric code (see `general.md`).

When generating B_c in the initial state, the special generator `BcVegPy` has to be used.
When running Gauss, the usual Pythia options then have to be replaced with the BcVegPy options.

Copies of the decfiles used for this analysis are located under `mc/decfiles`.

