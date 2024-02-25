

## Synthetic Data Summary

<img width="991" alt="image" src="https://github.com/JungiinChoi/prostate-active-surveillance-hier/assets/86948245/1ed67658-90c4-46df-beb9-1e293b69ef10">

## Data Format and Description

   - DX (Diagnostic)

| DX              | Description | Note |
| :---------: | :----------------------: | :-----------------: |
| ID        |   Patient ID   | |
| gprimp           |   True gleason primary   |can be NA |
| gsecondp    |  True gleason secondary   | can be NA |
| rp |  RP result   | can be NA |
| vol        |   Volume   | Optional |
| age_diag           |   Age at diagnosis (years)   |  |
| dxyear    |  Year of diagnosis   |  |
| dxrpdays |  Days from diag date to RP   |  |

   - BX (Biopsy)
| BX              | Description | Note |
| :---------: | :----------------------: | :-----------------: |
| ID        |   Patient ID   | |
| bxpos           |   Biopsy positive   |Set to 0 if gprim or gsecond is NA |
| gprim    |  Gleason primary   |  |
| gsecond |  Gleason secondary   |  |
| corestaken        |   Biopsy cores number dissected   | Optional |
| corespos           |Biopsy cores number pos| Optional |
| dxBXdays    |  Days from diag date to BX   |  |


   - PSA
| PSA              | Description | Note |
| :---------: | :----------------------: | :-----------------: |
| ID        |   Patient ID   | |
| psa           |   PSA ng/ml   | |
| dxPSAdays    |  Days from diag date to PSA   |  |


   - MRI

| MRI              | Description | Note |
| :---------: | :----------------------: | :-----------------: |
| ID        |   Patient ID   | |
| mriprosvol           |   MRI Prostate volume(cc)   ||
| mripirads    |  MRI PI-RADS  |  |
| dxMRIdays |  Days from diag date to MRI   |  |



## Notes on Data Format
- List the patients with gleason grade known from RP should be listed at the top of the diagnostic (dx) dataset.
- Only data collected before radical prostatectomy(RP) are used in all datasets (dx, bx, MRI, and PSA). 

## How to Run the Model?
1) Customize the variable `$work_dir` in the `submit.sh` script to match your working directory.
2) Running the model involves the following steps:
   - Execute `submit.sh` to run our model without MRI information.
   - Execute `submit.sh MRI` to include MRI-targeted information into biopsies.
   - Execute `submit.sh MRI_ST` to run the model which makes use of actual PIRADS score rather than (PIRADS <=3 or PIRADS >3 )
3) Ensure that all data is stored in the `data` directory in the specified format.
4) The fitted parameters and estimated cancer states will be automatically saved in the `generated-files` directory.
5) After obtaining estimated cancer states, compute the Area Under the Curve (AUC) by running the `AUC.R` script.

## Workflow

<img width="1322" alt="image" src="https://github.com/JungiinChoi/prostate-active-surveillance-hier/assets/86948245/ec3930fa-c49b-4961-880a-8b9be2b0774a">

  
## Libraries required
Model estimation:
lme4
splines
gtools
bayesm
rjags
R2jags

Model assessment:
mixtools
lme4
splines
gtools
bayesm
coda
ROCR
pROC
scales
grDevices
