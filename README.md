# ReliableRadiomicsVS

Code and partial data for *"Bayesian feature selection for radiomics
using reliability metrics"*

### Data

* Case study data is available [here](/CaseStudyData).
  * Data processing is also done in [HN_Aerts_Data.R](MatlabCode/HN_Aerts_Data.R)
* Reliability information is available [here](/CaseStudyData/LungProjectData). 
* Simulation data is created by the files beginning with `sim` in [MatLabCode](MatlabCode/)

### Code 

* `Harness.m` and `Harness_unbalanced.m` are the files that run the 25 simulations and record their results in the balanced and unbalanced settings, respectively. 
* `InputCaseStudy.m` sets the hyperpriors, calls `mainprogGBM.m` to run the chains, stores and pools the results, computes the marginal probabilites for each variable, then computes the predictions for each observation in the test set.  
* `mainprogRR.m` initializes the class specific means and calls the variable selection code.
* `bvsRR.m` is the primary variable selection code, calculating the pieces of the marginal log likelihood and then calling `metroRR.m`, the metropolis hastings code, then storing the results. 
* `metroRR.m` is the Metropolis Hastings step, containing the code for the add/delete/swap of the latent variable and updating the means, alphas, and intercepts as needed
  * `mmlRRdelta.m` and `mmlRgamma.m` compute the marginal log likelihoods of the means and the latent variables, respectively. 
