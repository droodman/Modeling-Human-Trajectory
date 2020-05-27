Code and data archive for working paper version of Roodman, "Modeling the Human Trajectory"

"GWP.xlsx" exhibits the data set construction and its "Data" tab is accessed directly by "Modeling GWP.do"
"Modeling GWP.do" generates all graphs, tables, and other estimates in the paper, except those generated by...
"Multivariate simulator.do", which generates the multivariate simulations, with and without endogenous natural resources
"Results.zip" is an archive of graphs, tables, and Stata estimation files generated by the do files, most of which do not appear in the text.

"Modeling GWP.do" requires the asdf package, which is in a primitive state but can be installed in Stata with
  "net install asdf, from(https://raw.github.com/droodman/asdf/v0.1.0) replace"
Currently the package requires Stata 16 or later.

Also needed are Ben Jann's grstyle, estout, and colorpalette packages.