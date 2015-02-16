
# Analysis Outline

## Goal

Determination of

    - B0→D~0mumu
    - B0→D~0J/psi
    - (B+→D\*+mumu)
    - (B+→D\*+J/psi)
    - (B→Ds+mumu)
    - (B→Ds+J/psi)

branching fractions.

## For each decay

- Retrieve dataset through B2Xmumu stripping line
  - Calculate isolation variable
- Retrieve signal simulation
- Apply blinding by cut on B signal region
- Apply cuts to PID variables
  - Select for true Kaons/pions
- Apply multivariate classifier (AdaBoost)
  - Use decay topology/kinematics variables
  - Use isolation variable
  - Use PID variables?
  - Choose optimization figure (look into similar analyses)
  - Apply cut on discriminant according to optimization figure
- Perform fit of B mass
  - Include D mass (D lifetime) for better control of K pi background
  - Profile likelihood ratio to account for nuisance parameters
  - CLs method to test signal/background-only hypotheses
- Determine trigger efficiency through TISTOS
- Determine PID efficiency through MC resampling
- Determine selection efficiency of multivariate classification
- Determine LHCb acceptance efficency
- Find normalization channel (B0→K\*J/psi?)
  - Keep selection similar
  - Perform fit of B mass
  - Incorporate into result of analysis to reduce systematics
      - Use known branching fraction to calculate result
      - Or: only calculate ratio of branching fractions
- Calculate estimate of signal branching ratio from signal yield, efficiencies and normalization channel
- Calculate systematic uncertainties
  - Number of simulated signal candidates
  - Choice of binning scheme for PID MC resampling
  - Systematics related to 

