# Search for B→Dμμ at LHCb

Rare decays of the form B→Dμμ could provide an opportunity to check our understanding of the Standard Model and provide hints of New Physics.
None of the possible B→Dμμ decays have previously been discovered.
With the help of the LHCb detector at CERN, we might now be able to discover and study these decays for the first time.

## Setup

To repeat the analysis, install the dependencies with
```
pip install -r ./requirements.txt
```

Define the path to the datastore:
```
export DATASTORE=/path/to/data/
```

Run the analysis pipeline with
```
python pipeline.py
```

## Repository layout

The top level of the repository contains

 - The entry point to the analysis pipeline: `pipeline.py`
 - The `logbook` directory containing information on how the analysis was planned and performed
 - The `lhcb` directory containing customized LHCb software and configuration for fetching all analysed data from the grid
 - The `analysis` python package containing code used by the analysis pipeline


