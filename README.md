Yates Coley
rycoley@gmail.com
December 27, 2017

Contents of this readme:

1. Abbreviations and definitions used in code and annotation
2. Workflow for model estimation
3. Workflow for model cross-validation
4. Workflow for model evaluation
5. Libraries required


1. Abbreviations and defintions
## Abbreviations
bin: binary, used to indicate ordinal variables (bx or pathological grade) that are dichotomized 
bx: biopsy
cv: cross validation
dob, dod: date of birth, death
dx: diagnosis 
dt: date
ea: each
fup: follow-up
GS: Gleason score
IOP: informative observation process. This describes possible informative nature of biopsy and surgery decisions. Bx and surg decisions were not related to PGG in our analysis, so these components were not included in final estimation but they can be seen in the data shaping code (and are commented out in model estimation scripts)
LR: low risk (vs. very low risk) for bx dx volume criteria
lat: laterality (i.e., uni- or bilateral cancer)
LTF: loss to follow-up
MPC: maximum percent cancer involvement for any positive cores on biopsy
ncs: number of cores sampled on biopsy
NPC: number of positive cores (found on biopsy)
PGG: pathological grade group, sometimes written PPGG (pathological prognostic grade group) or PGS (prognostic Gleason score). (Our terminology has evolved as we’ve worked on this project for several years.)
pt: patient
rc: biopsy grade reclassication 
rm: “remove”, an indicator variable used in data processing to mark records for removal (e.g., patient doesn’t meet criteria, needs to be removed)
RP: radical prostatectomy, sometimes referred to as RRP (radical retropubic prostatectomy) or surg (surgery)
subj: subject, a different unique patient identifier that uses patient ordering needed for JAGS run 
surg: surgery
time.int: time interval (for annual biopsy intervals)
tx: treatment
vol: volume

## Definitions
B: length of saved posterior chain
eta: PGG
n: number of individuals


2. Workflow for model estimation
## Data
User should start with four datasets (name given in my code):
1. patient-level dataset with demographic information (including date of birth, death) (“demo.data”)
2. biopsy-level dataset with one row per biopsy (diagnostic and follow-up biopsies included) with biopsy date and findings specified (“bx.data”)
3. PSA-level dataset with one row per PSA observations (including pre-diagnosis) with test date and result (“psa.data”)
4. treatment-level dataset with one row per treatment received. There may be multiple rows (treatments) per patient if more than one treatment were received. Treatment type and, if surgery, post-RP PGG should be noted. (“tx.data”)

## Model estimation
User can run model estimation scripts in order by executing “to-run.R”. The user should amend this script to correctly reference file directories, file names, and the date of the data pull. Model estimation was performed using JAGS and R on a SGE high performance computing cluster and submitted via batch script. Each task was assigned an id “SGE_TASK_ID”, which is used to set the starting seed of model estimation (and ensure reproducitbility.)

Each script called by to-run.R is described in turn below:

Script first runs data-load-check-and-shaping.R to check data values and shape 4 datasets into format for model estimation. Warning messages will be produced to describe data processing and any problems identified in data. Execution will only be halted for flaws that prevent model estimation from proceeding. This creates the file “IOP-data-shaping-work-spaced.Rdata”.

Next, data-prep-for-jags.R performs further data manipulation to translate the formatted datasets into individual variables, vectors, and matrices for JAGS estimation. Variable names are shortened for brevity and generally correspond to variables and parameter labels given in Coley et al. (2017) Biometrics description of this mode.  

Third, argument-prep-for-jags.R defines the arguments needed for the jags function (data, initial parameter values, parameters to track, length and number of chains, length of burn-in, and number to thin. 

Fourth, the JAGS model is written to a txt file by sourcing JAGS-prediction-model.R. 

After sourcing these scripts, to-run.R calls the jags() function to perform model estimation.



3. Workflow for model cross-validation
Leave-one-out cross-validation follows a similar workflow to primary model estimation (and same data sources are needed). Run cv-to-run.R to complete data shaping and formatting, argument definition for JAGS, JAGS model definition, and model estimation via jags(). R scripts particular to cross validation have prefix cv-. In cv-data-prep-for-jags.R the post-RP PGG for one patient (with true state known) is masked in order to repeat model estimation and generate a true state prediction for that individual. Renamed variables are propagated through cv-argument-prep-for-jags.R and cv-JAGS-prediction-model.R.



4. Workflow for model evaluation

For initial model assessment, run get-model-estimation-assessment.R. This script will run the following programs (in addition to data load and prep scripts from to-run.R):
calibration-plot-observables.R generates calibration and ROC plots for biopsy grade reclassification predictions
convergence-diagnostics.R generates trace plots, cumulative quantile plots, and posterior density plots for model parameter in order to evaluate model convergence
psa-lme-plot.R creates a scatterplot of individual-level intercept and slope estimates with credible intervals for the posterior bivariate normal distribution of mean intercept and slope for patients with PGG=1 and PGG>1.

For predictive accuracy assessment, run cv-get-predictive-accuracy.R. This script will calculate and save the AUC, MAE, and MSE and create calibration and ROC curve plots. 

To generate patient-level PGG and biopsy grade reclassification predictions, run get-patient-predictions.R,  which will run the following R scripts: patient-predictions-pgs.R and patient-predictions-biopsy-reclassification.R.

The following scripts provide additional model summaries and assessments presented in the paper and appendix:
examine-predictions-associations-eta.R generates a figure showing the population-level estimated probability of PGG>1 for a patient enrolled in AS given age and volume at diagnosis. This script also produces summary estimates of the PGG distribution (probability PGG=1, =2, =3, and =4) for the average age at diagnosis and low or high cancer volume at diagnosis. 

examine-predictions-associations.R generates plots examining the association between patient-level PGG predictions and PSA intercept and slope estimates (Figure 2-4 in PSA paper) and residuals from the mixed effects model. get-patient-predictions.R must be run first. 



5. Libraries required
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
