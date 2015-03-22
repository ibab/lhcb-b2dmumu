# B→Dμμ Analysis

This is the repository for my master's project at TU Dortmund.

During the course of my project, I will be looking for the rare decay Bd → Dbar mu mu and similar decays with [LHCb](http://lhcb.web.cern.ch/lhcb/) data.

Should the decay be found, a measurement of its branching fraction might be possible.

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


