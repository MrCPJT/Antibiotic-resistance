# Final year project
Final year project investigating the quantification of microbial interactions, particularly in quantifying dynamic growth rates and interactivity parameters. 

## Timeline

### Literature review
- Reading of relevant literature - Getting a feel for the field
### Reproduction of literature results:
- 11D generalized Lotka-Volterra (gLV) model [1,2,3]
  - Solved the system of ODEs using initial conditions and parameter estimates from existing literature [1]
  - Extended the basic gLV equations to include a perturbation / sensitivity term [1]
  - Investigated the importance of timing when issuing fecal microbiota transplant (FMT) [3]
  - Considered varying interactivity parameters at each time step using a 1D random walk
  - Analysed system dynamics (mostly identification and classifcation of steady states) [2,3]
  
 - Mechanistic and pairwise models [4]
    - Evaluated and compared pairwise and mechanistic equations for the following cases:
      - 2 Species interacting via a single consumable metabolite (S1 -> C1 -> S2)
      - 2 Species interacting via two consumable metabolites (More complex interactions

     - In our first case, the dynamics are relatively straight forward and didn't warrant further investigation. One may additionally consider the case of **reusable** chemical mediators, but such was not discussed in the literature.
     - In the second case, we considered different scenarios (mutualism, competition, predation), based on initial conditions and interactivity parameters, and investigated how such can affect the system's dynamics.
     
### Data fitting
- Given novel mono and co-culture data for two bacterial species; attempted to fit the data to existing models (gLV and pairwise) using a collection of different methods, including:
  - Ordinary least-squares regression (MATLAB)
  - Bayesian Inferencing (R, Stan)
 
 
 
 ## References
 - [1] https://doi.org/10.1371/journal.pcbi.1003388
 - [2] https://doi.org/10.1371/journal.pcbi.1006001
 - [3] https://doi.org/10.1103/PhysRevE.99.032403
 - [4] http://dx.doi.org/10.7554/eLife.25051
 
