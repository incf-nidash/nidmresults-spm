
# NIDM-Results for SPM

Export mass-univariate neuroimaging results computed in SPM as NIDM-Results packs.

A *NIDM-Results pack* is a compressed file containing a NIDM-Results serialization and some or all of the referenced image data files in compliance with [NIDM-Results specification](http://nidm.nidash.org/specs/nidm-results.html).

##### Usage
In Matlab, open `SPM`
```
spm fmri
```
 1. Open the `Batch Editor` by clicking on the `Batch` button in the SPM12 `Menu` window
 1. Open the menu `SPM` > `Stats` > `Results Report`
 
 1. In the batch window:
  - Fill in information about the results you are interested in (in particular `SPM.mat` file, contrast number, threshold, etc.)
  - In `Export results using the Neuroimaging Data Model (NIDM)?`, selected `Yes`
  - Fill in information about your analysis (`Single-subject analysis`/`Group analysis`, `Modality`, `Reference space`, etc.)
  
##### Requirements
  - [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

##### Installation
Download the latest release (**12.575ac2c**)
```
curl -sLo -  https://github.com/incf-nidash/nidmresults-spm/archive/12.575ac2c.tar.gz | tar xzvf -
export exporter_path=`pwd`/nidmresults-spm-12.575ac2c/exporter
```
Copy the core SPM files modified by the exporter (replace by `<PATH_TO_SPM>` by the full path to your spm installation)
```
export spm_path="<PATH_TO_SPM>"
cp -p $spm_path/spm_cfg_results.m $spm_path/spm_cfg_results_ORIGINAL.m
cp -p $spm_path/spm_run_results.m $spm_path/spm_run_results_ORIGINAL.m
cp -p $exporter_path/spm_cfg_results.m $spm_path/spm_cfg_results.m
cp -p $exporter_path/spm_run_results.m $spm_path/spm_run_results.m
```
In Matlab, update your path (replace by `<PATH_TO_SPM>` by the output of `echo $spm_path` and `<PATH_TO_EXPORTER>` by the output of `echo $exporter_path`)
 ```
 addpath('<PATH_TO_SPM>')
 addpath('<PATH_TO_EXPORTER>')
 ```

 


##### Tests

Testing procedures for SPM NIDM-Results export.

Test data is available at https://github.com/incf-nidash/nidmresults-examples/, you will need a local copy of this repository stored with [git lfs](https://git-lfs.github.com/):
```
git clone https://github.com/incf-nidash/nidmresults-examples.git
git lfs install
```

To run the test battery, follow those steps:
 1. Run the NIDM-Results export on your machine (within Matlab)
```
nidm_export_all('LOCAL_PATH_TO_NIDMRES_EX/nidmresults-examples/', 'LOCAL_PATH_TO_NIDMRES_SPM/nidm-results_spm/spmexport')
``` 
 2. Push the updates ttl and provn file to GitHub:
```
    git add -u
    git commit -m "description of updated feature"
    git push origin master
```
