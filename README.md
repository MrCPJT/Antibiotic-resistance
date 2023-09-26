# Final year project
Final year project concerned with understanding the microbial interactions between _Staphylococcus aureus_ and _Pseudomonas aeruginosa_ in co-culture.
## Project Details

### Experimental Data
- For this project I had the opportunity to work with novel _S. aureus_ and _P. aeruginosa_ mono and co-culture data. The data (.csv) is featured above and features a mono and co-culture case for both _S. aureus_ and _P. aeruginosa_ with each case including 6 replicates to identify any outliers.

### Literature review 
- Gentle introduction to general literature: key definitions (microbes, microbial communities, etc.), introduction of _S. aureus_ and _P. aeruginosa_ (what they are, why they are important, motivations behind understanding interactions between the two, etc.), background information and a brief history on microbial interactions.
### Reproduction of literature results:
The purpose of this step was to improve my understanding of existing models by implementing them myself. Not all of the points mentioned below made it into the final report but were interesting to investigate nonetheless.

- 11D generalized Lotka-Volterra (gLV) model [1,2,3]
  - Solved the system of ODEs using initial conditions and parameter estimates from existing literature [1]
  - Extended the basic gLV equations to include a perturbation / sensitivity term [1]
  - Investigated the importance of timing when issuing feacal microbiota transplant (FMT) [3]
  - Considered varying interactivity parameters at each time step using a 1D random walk
  - Analysed system dynamics (mostly identification and classification of steady states) [2,3]
  
 - Mechanistic and pairwise models [4]
    - Evaluated and compared pairwise and mechanistic equations for the following cases:
      - 2 Species interacting via a single consumable metabolite (S1 -> C1 -> S2)*
      - 2 Species interacting via two consumable metabolites (More complex interactions)**

     - *In our first case, the dynamics are relatively straight forward and didn't warrant further investigation. One may additionally consider the case of **reusable** chemical mediators, but such was not discussed in the literature.
     - **In the second case, we considered different scenarios (mutualism, competition, predation), based on initial conditions and interactivity parameters, and investigated how such can affect the system's dynamics.
     
### Parameter Estimation
Given novel mono and co-culture data for two microbial species we estimated parameter values for existing models (gLV):
- Least-squares regression (MATLAB)
  - Make sensible assumptions and approximations to parameter values -> Iterate through pairs of estimate parameters to find minimum distance -> Update search
  - Requires a lot of manual exploration (updating of parameters, exploring parameter space for minimum, etc)

- Bayesian Inferencing (R, Stan)
  - Inferencing conducted using R and RStan package.
  - MCMC sampling/estimates are enabled through RStan (powerful and robust - usually requiring C++)
  - Sampling requires a master 'R script' and a respective 'Stan file' to call (details can be found in the repository)
 
 ### Systems Analysis
 Having generated some parameter estimates we can investigate the underlying dynamics behind S.Aureus and P.Aeruginosa interactions. We first analytically investigate the generalised Lotka-Volterra equations and discuss some of the expected behaviour for our 2D system. 
 
 ## References
 - [1] https://doi.org/10.1371/journal.pcbi.1003388
 - [2] https://doi.org/10.1371/journal.pcbi.1006001
 - [3] https://doi.org/10.1103/PhysRevE.99.032403
 - [4] http://dx.doi.org/10.7554/eLife.25051
 
 ## Useful resources for R, RStan, Stan
 - https://shug3502.github.io/blog/DifferentialEqnsStan
 - https://betanalpha.github.io/assets/case_studies/markov_chain_monte_carlo.html#1_hello_monte_carlo_our_old_friend
 - https://alexanderetz.com/2015/07/25/understanding-bayes-updating-priors-via-the-likelihood/
 - https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html
 - General documentation + Stack exchange
 
