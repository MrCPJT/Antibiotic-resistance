# Mathematically-Modelling-Microbial-Interactions
Final year project concerned with understanding the microbial interactions between _Staphylococcus aureus_ and _Pseudomonas aeruginosa_ in co-culture.

## Summary
- Performed an extensive literature review, implemented existing models and reproduced results to establish an initial understanding of the field
- Preprocessed microbial growth data before fitting to an existing mathematical model using the Hamiltonian Monte Carlo algorithm, validating fits using MCMC convergence diagnostics
- Investigated underlying dynamics using steady-state analysis and bifurcation theory, finding a transcritical bifurcation to underpin changes in interactions

## Project Details
> Please see PDF file for additional details.

### Experimental Data
- For this project I had the opportunity to work with novel _S. aureus_ and _P. aeruginosa_ mono and co-culture data. The data (.csv) is featured above and includes a mono and co-culture case for both _S. aureus_ and _P. aeruginosa_ with each case including 6 replicates

### Literature Review 
- Gentle introduction to surrounding literature: key definitions (microbes, microbial communities, etc.), introduction of _S. aureus_ and _P. aeruginosa_ (what they are, why they are important, motivations behind understanding interactions between the two, etc.), background information and a brief history on microbial interactions

### Reproduction of Literature Results
> The purpose of this step was to improve my understanding of existing models by implementing them myself. Not all of the points mentioned below made it into the final report but were interesting to investigate nonetheless.

- 11D generalized Lotka-Volterra (gLV) model [1,2,3]
  - Solved the system of ODEs using initial conditions and parameter estimates from existing literature [1]
  - Investigated the importance of timing when issuing faecal microbiota transplant (FMT) [3]
  - Analysed system dynamics (mostly identification and classification of steady states) [2,3]
  
 - Mechanistic and pairwise models [4]
   - Evaluated and compared pairwise and mechanistic equations for the following cases:
     
     - 2 Species interacting via a single consumable metabolite
     - 2 Species interacting via two consumable metabolites

### Preparing the Data
- Calculated averages (partly to eliminate null values) and performed normalisation improving the data quality

### Model Discussion
- Introduced the model of choice (gLV) and discussed relevant reparameterisations and parameter constraints
     
### Parameter Estimation
- Performed Bayesian inferencing (using RStudio and Stan) to estimate unknown model weights
  
  - MCMC sampling/estimates are enabled through RStan (powerful and robust - usually requiring C++)
  - Sampling requires a master `.R script` and a respective `.stan file` to call (see [5] for more information)

### Steady-State Analysis
- Analytically investigated the gLV equations, finding the nullclines and equilibria of our model
  
### Dynamical Systems Analysis
- Classified stability of equilibrium points and discussed what different behaviour cases translates to
- Touched on bifurcation theory and what qualitative changes we might expect our system to demonstrate before performing bifurcation analysis using MATCONT
 
 ## References
 - [1] https://doi.org/10.1371/journal.pcbi.1003388
 - [2] https://doi.org/10.1371/journal.pcbi.1006001
 - [3] https://doi.org/10.1103/PhysRevE.99.032403
 - [4] http://dx.doi.org/10.7554/eLife.25051
 - [5] https://mc-stan.org/rstan/reference/stan.html
 
 ## Useful Resources for R, RStan, Stan
 - https://shug3502.github.io/blog/DifferentialEqnsStan
 - https://betanalpha.github.io/assets/case_studies/markov_chain_monte_carlo.html#1_hello_monte_carlo_our_old_friend
 - https://alexanderetz.com/2015/07/25/understanding-bayes-updating-priors-via-the-likelihood/
 - https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html
 - General documentation + Stack exchange
