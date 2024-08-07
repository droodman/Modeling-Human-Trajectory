Code and data archive for working paper version of Roodman, "Modeling the Human Trajectory"

"GWP.xlsx" exhibits the data set construction and its "Data" tab is accessed directly by "Modeling GWP.do"

"Modeling GWP.do" generates all graphs, tables, and other estimates in the paper, except those generated by...

"asdfSuperNL.ado" is a helper program run after superexponential fits with asdf. It adds to the estimation set estimates of probability of no takeoff, median takeoff, etc.

"Multivariate simulator.do", which generates the multivariate simulations, with and without endogenous natural resources

"Plots of Feller solutions.do" plots the Feller and noncentral chi2 distributions (main text) and the corresponding diffusions (appendix)

"graphs and estimates.zip" is an archive of graphs, tables, and Stata estimation files generated by the do files, most of which do not appear in the text.

"AggregateBLSIPP.do" is a version of the Matlab program AggregateBLSIPP.m from Bloom et al. (2020) data+code archive, which adapts from that paper an estimate of \phi_A.

"Modeling GWP.do" requires the asdf package, which is in a primitive state but can be installed in Stata with
  "net install asdf, from(https://raw.github.com/droodman/asdf/v0.1.0) replace"

The package requires Stata 16 or later.
