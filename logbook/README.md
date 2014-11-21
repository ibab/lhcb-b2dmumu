
# Chronological log

### 2014-10-31

- Started this chronological log.
- Familiarized myself with MC generation at LHCb using Gauss.
- Wrote a decfile for Bc+->Ds+MuMu.

### 2014-11-03

- Request for B->Dmumu Monte Carlo sent out.
 - They asked if we also want B+->D+mumu
- Sent a mail to LHCb Gauss Managers asking if my Bc->Dsmumu decfile is OK.
- Sketched out first steps of analysis:
 - Generate dataset from B->Xmumu stripping line.
 - Request Monte Carlo for interesting channels
 - Start training classifier with B->K\*mumu MC

### 2014-11-11

- Nearly finished configuration for pulling data (`data` subdirectory)
- Finished configuration for pulling B->D0mumu Monte Carlo (`mc` subdirectory)
- Submitted jobs for pulling B->D0mumu Monte Carlo dataset
- Something to keep in mind for later: uniform classifiers [https://github.com/anaderi/lhcb_trigger_ml](https://github.com/anaderi/lhcb_trigger_ml)

### 2014-11-21

- Started building an analysis pipeline
  - Reduce variables
  - Blind signal region
  - Perform preselection
  - Add mis-id mass hypotheses
  - Classify signal/background with AdaBoost
  - Create plots of all variables
  - Still missing: perform MLE on the resulting mass distribution
- Bd -> K\*mumu data might be used to evaluate the BDT, especially wrt to bias in mass variable
  - Implemented a second pipeline for running the Bd -> D0mumu classifier on Bd -> K\*mumu data
  - Clear mass peak appears when applying the classifier
- Started assembling a list of preselection criteria


