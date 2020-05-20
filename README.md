The MATLAB® code in this reposirtory runs the SEPIA epidemiological
model and reproduces the scenarios shown in figure 2 and S1 of the
manuscript

["The geography of COVID-19 spread in Italy and implications for the
relaxation of confinement
measures"](https://www.medrxiv.org/content/10.1101/2020.04.30.20083568v2)

Launch `SEPIA_main.m` to run the model. The code reads the data and a
sample of the posterior distribution of parameters and runs the model
producing the two figures contained in the folder `Expected_Results`.

To produce different scenarios, the user can change the parameter
`betaInc` (line 99 of `SEPIA_main.m`). Figure 2 and S1 of the preprint
have been obtained using `betaInc=1`, `betaInc=1.2` and `betaInc=1.4`.

The code runs 2000 realizations of the model and takes less than a
minute in a standard laptop or desktop computer. To obtain a more
robust and smoother estimation of the probability distribution of the
output, the user can increase the number of posterior samples used
(parameter `NSample`, line 93 of `SEPIA_main.m`).

The code has been tested on MATLAB versions R2018a and R2019b.
