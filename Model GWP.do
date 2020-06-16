*** Estimation code for Roodman, "Modeling the Human Trajectory" as of May 22, 2020
*** Requires Ben Jann's estout and grstyle packages, as well as the asdf package distributed with this file
*** asdf version 0.1.0 can be installed with "net install asdf, from(https://raw.github.com/droodman/asdf/v0.1.0) replace"
***
*** The estimation code below is complex, but performing a basic estimate only takes a few lines:
***   PrepData if Year>=-10000, depvar(GWP) historicpop(PopMcEvedyJones) prehistoricpop(PopDeevey)  // make data set starting in 10,000 BCE
***   mat init = 0,0,-1,-1                                                                      // starting estimates for lna, b, nu, gamma
***   asdf GWP [aw=HYDEwt], model(bernoudiff) reflect(0) from(init)       // fit non-reflecting Bernoulli diffusion with HYDE-based weights


set scheme s1mono
set autotabgraphs on
graph set window fontface Cambria
grstyle init
grstyle linewidth plineplot vthin
set seed 902381476

scalar T = 5000       // end-year for plots, e.g., 5000 CE
scalar Tsamp = 10000  // number of evenly spaced observations that simulator should return for each path
scalar Tres = 10      // number of time steps per returned time point
scalar M = 10000      // number of paths to generate when plotting path distributions

cap program drop PrepData
program define PrepData

  global GWPvar GWP
  global Popvar Pop
  global GWPcapvar GWPcap
  global GDPcapvar FRAGDPcapMaddison2018

	import excel GWP, firstrow sheet(Data) cellrange(A2:BR160) clear

	keep Year PopMaddison* PopDeevey PopHYDE* PopUN* GWPcapMaddison* PopMcEvedyJones* *WEO* FRAGDPcapMaddison2018 PopKremer PopgrowthKremer GWPDeLong GWPHansonreconstructed
  ren GWPHansonreconstructed GWPHanson

	replace Year = 0 if Year==1
	tsset Year

	tempvar n tmp lnpop

	gen double `tmp' = cond(Year<=-10000, 1, cond(Year==0, .75, cond(Year==1700, .25, cond(Year==1900, .05, cond(Year>=2000, .01, .)))))  // HYDE 3.2.1 uncertainty bounds. E.g. .75 = +/-75%. Bounds are a linear spline w.r.t. year
	ipolate `tmp' Year, gen(HYDEsd)
	gen double HYDEwt = HYDEsd*HYDEsd
	replace HYDEwt = 1 / (1+2*HYDEwt)

	syntax [if], depvar(string) historicpop(varname) prehistoricpop(varname) [subsistence(real 0) gwpcap0(real 400) civstart(real -10000) gwpstart(real 1500) historystart(real -10000) region(string)]

	gen double Pop = cond(Year<`historystart', cond("`region'"=="",`prehistoricpop',.), cond(Year<`gwpstart', `historicpop'`region', PopMaddison2010`region'))
	gen long `n' = _n
	sum `n' if Pop<.
	replace Pop = Pop[r(max)] * PopUN`region' / PopUN`region'[r(max)] if _n > r(max)  // extrapolate forward from Maddison population series using UN

	gen double GWPcap = GWPcapMaddison2010`region' if Year>=`gwpstart'
	replace GWPcap = `gwpcap0' if Year<`civstart' & Pop < .
	foreach range in "1000 `gwpstart'" "0 1000" "`civstart' 0" {  // interpolate using Maddision (2010) values for 1 and 1000 AD, and subsistence
		local ti: word 1 of `range'
		local tf: word 2 of `range'
		sum `n' if Year>=`ti' & Year<=`tf', meanonly
		scalar lnP0 = ln(Pop[r(min)])
		scalar lnPf = ln(Pop[r(max)])
		scalar yi = cond(`ti'>`civstart', GWPcapMaddison2010`region'[r(min)], `gwpcap0')
		scalar yf = cond(`tf'<`gwpstart', GWPcapMaddison2010`region'[r(max)], GWPcap[r(max)])
		replace GWPcap = yi * (yf / yi) ^ ((ln(Pop) - lnP0)/(lnPf - lnP0)) if Year>=`ti' & Year<=`tf'
	}
	replace GWPcap = GWPcap - `subsistence'  // optionally take GWP/capita relative to a zero-investment level

	gen double GWP = GWPcap * Pop / 1000
	replace GWP = GWP[_n-1] * (1 + GWPgrowthWEO) if Year>=2000 & GWP==.
	replace GWPcap = GWP / Pop * 1000 if GWPcap==. & GWP<.
	replace $GDPcapvar = . if Pop==.
  replace $GDPcapvar = $GWPcapvar / 17647.0052262412 * 29031 if Year<=-5000  // take early obs of France GDP/capita from GWP, adjusting 1990 $ -> 2011 $ using Maddison project (2013, 2018) figs for France 1990
  gen n = _n
  sum n if $GDPcapvar < ., meanonly
  replace $GDPcapvar = $GDPcapvar[_n-1] * (1 + FRAGDPcapitagrowthWEO) if _n > r(max)
  
	char $GWPvar[unit] "$ billion"
	char $Popvar[unit] "Population (million)"
	char $GWPcapvar[unit] "GWP/capita ($)"
	char $GDPcapvar[unit] "France GDP/capita ($)"

	cap keep `if'
	cap keep if `depvar'<.

	gen strYear = "" + cond(Year, string(abs(Year)) + cond(Year<0, " BCE", ""), "1") + ""
  gen tdelta = cond(_n==1, 0, Year - Year[_n-1])  // putting a 0 in first entry prevents premature clipping of sample in static regressions: first entry has info needed to set up those regressions
	gen double ydot = (`depvar'[_n+1] / `depvar') ^ (1 / tdelta[_n+1]) - 1

	gen double ln`depvar' = ln(`depvar')
	sum Year, meanonly
	scalar ti = r(min)  // starting time
	scalar tf = r(max)  // final time
	sum `depvar' if Year==ti, meanonly
	scalar Yi = r(mean)
	scalar lnfloor = ln(1000)*ceil(ln(r(max))/ln(1000)-1)  // for clipping in graphs
	sum `depvar' if Year==tf, meanonly
	scalar Yf = r(max)
	scalar lnceiling = ln(1000)*ceil(ln(r(max))/ln(1000)+1)
end


* given logged variables, return appropriate argument for a graph xlabels() or ylabels() option
cap program drop myNiceLogLabels
program define myNiceLogLabels, rclass
	syntax varlist [if] [in]
	tempvar min max
	qui egen double `min' = rowmin(`varlist') `if' `in'
	qui egen double `max' = rowmax(`varlist') `if' `in'
	sum `min' `if' `in', meanonly
	local lo = r(min)
	sum `max' `if' `in', meanonly
	local hi = r(max) + ln(.8)  // hack to delete uppermost label if it is too close to top of data/graph
	local step = (`hi' - `lo') / ln(10)
	local step = cond(`step' < 8, 1, cond(`step' < 48, 3, 6)) * ln(10) // go one or three orders of magnitude at a time
	local lo = ceil(`lo'/`step'-.1) * `step'
	local hi = floor(`hi'/`step'+.1) * `step'
	forvalues i=`lo'(`step')`hi' {
		local text = string(10^(round(`i'/ln(10))), "%10.1gc")
		local labels `labels' `=`i''  "`text'"
	}
	return local labels `labels'
end


* After asdf estimation, make path and path distribution plots
* depvar() is short name of variable modeled--GWP, Pop, GWPcap, GDPcap
* modelname(), tsamplename(), region() are for labelling saved graph files too
* mediantakeoff() is Y in ln(Y - Year) when scaling time axis in logs
* nov prevents taking of distinct parameter values for each path
* xlabels() allows custom labelling of x axis tick marks
* noclose prevents closing each graph after it is generated and saved; closing can save memory

cap program drop Plot
program define Plot
	syntax, depvar(string) mediantakeoff(string) [nov modelname(string) tsamplename(string) region(string) XLABels(string asis) noclose]

	asdfSim, t0(`=ti') y0(`=Yi') tsamp(`=Tsamp') tres(`=Tres') m(`=M') nq(20) ns(99) `v' t(`=T-ti')  // simulate from initial obs
  mata lnYq = ln(Yq); lnYs = ln(Ys)
	getmata Yeari=t (Yqi*)=lnYq (Ys*)=lnYs, double force replace

	* plot sample paths
	qui foreach var of varlist Ys* {
		replace `var' = cond(`var'[_n-1]>=lnceiling               , ., lnceiling) if `var'>=lnceiling & `var'<.  // for display, censor extreme values 
		replace `var' = cond(`var'[_n-1]<=lnfloor | `var'[_n-1]==., ., lnfloor  ) if `var'<=lnfloor
	}
	myNiceLogLabels Ys* ln${`depvar'var}

	local ylabels `r(labels)'
	local name `modelname'Paths`depvar'`tsamplename'`v'`region'
	line Ys* Yeari, legend(off) xtitle(Year) ylabels(`ylabels', notick angle(hor)) xlabel(`xlabels') ytitle("") graphregion(margin(l=0)) subtitle(`:char ${`depvar'var}[unit]', pos(11) span size(medsmall)) name(`name', replace) || ///
		connected ln${`depvar'var} Year if ln${`depvar'var}<., lcolor(red) mcolor(red) lstyle(solid) lwidth(vthin) msize(tiny) msymbol(O) plotregion(lwidth(none) margin(zero)) graphregion(margin(0 7 0 0))
	graph save `name', replace
	graph export `name'.png, width(2680) height(1552) replace
  if "`close'"=="" graph close `name'

	if `mediantakeoff' < . {
    asdfSim, t0(`=tf') y0(`=Yf') tsamp(`=Tsamp') tres(`=Tres') m(`=M') nq(20) `v' t(`=T-tf')  // simulate from final obs
    mata lnYq = ln(Yq)
    getmata Yearf=t (Yqf*)=lnYq, double force replace
    cap drop lnYear
    cap gen double lnYear = ln(`mediantakeoff' - Year)
		cap drop lnYeari
		gen double lnYeari = ln(`mediantakeoff' - Yeari)
		cap drop lnYearf
		gen double lnYearf = ln(`mediantakeoff' - Yearf)
		myNiceLogLabels lnYeari lnYearf if Yearf<`mediantakeoff'-min(.99, (`mediantakeoff' - tf)/2)
		local xloglabels `r(labels)'
		local name `modelname'Paths`depvar'`tsamplename'Log`v'`region'
		line Ys* lnYeari, legend(off) || ///
			connected ln${`depvar'var} lnYear if ln${`depvar'var}<., lcolor(red) mcolor(red) lstyle(solid) lwidth(vthin) msize(tiny)	msymbol(O) ///
			xscale(reverse) xtitle(Years till `=round(`mediantakeoff')') ylabels(`ylabels', angle(hor) notick) ytitle("") subtitle(`:char ${`depvar'var}[unit]', pos(11) span size(medsmall)) ///
			xlabels(`xloglabels') ///
			plotregion(lwidth(none) margin(zero)) graphregion(margin(0 7 0 0)) ///
			name(`name', replace)
		graph save `name', replace
		graph export `name'.png, width(2680) height(1552) replace
    if "`close'"=="" graph close `name'
	}

	* plot evolution of distribution of simulations
	preserve
		if `mediantakeoff' < . {  // depict distribution starting from final observation too?
			cap drop n
			gen long n = _n
			foreach if in i f {
				sum n if Year`if' > `mediantakeoff' - 1, meanonly
				qui forvalues q=1/19 {
					// shift first obs clipped on right to right edge, interpolating to prettify graph
					replace Yq`if'`q' = (`mediantakeoff' - 1 - Year`if'[_n-1]) / (Year`if' - Year`if'[_n-1]) * (Yq`if'`q' - Yq`if'`q'[_n-1]) + Yq`if'`q'[_n-1] in `r(min)'
					replace Yq`if'`q' = lnceiling if Yq`if'`q'>=lnceiling  // for display, censor extreme values 
					replace Yq`if'`q' = lnfloor   if Yq`if'`q'<=lnfloor | (Yq`if'`q'==lnfloor & Yq`if'`q'==.)
				}
				replace Year`if' = `mediantakeoff' - 1 in `r(min)'
				replace lnYear`if' = ln(`mediantakeoff' - Year`if') in `r(min)'
				replace lnYear`if' = . if _n > `r(min)'
			}
			myNiceLogLabels Yqi* Yqf*
			local ylabels `r(labels)'

			local name `modelname'Dist`depvar'`tsamplename'Log`v'`region'
			twoway rarea Yqi19 Yqi1 lnYeari, fcolor(black*.10) lwidth(none) || rarea Yqi18 Yqi2 lnYeari, fcolor(black*.20) lwidth(none) || rarea Yqi17 Yqi3 lnYeari, fcolor(black*.30) lwidth(none) || rarea Yqi16 Yqi4 lnYeari, fcolor(black*.40) lwidth(none) || rarea Yqi15 Yqi5 lnYeari, fcolor(black*.50) lwidth(none) || rarea Yqi14 Yqi6 lnYeari, fcolor(black*.60) lwidth(none) || rarea Yqi13 Yqi7 lnYeari, fcolor(black*.70) lwidth(none) || rarea Yqi12 Yqi8 lnYeari, fcolor(black*.80) lwidth(none) || rarea Yqi11 Yqi9 lnYeari, fcolor(black*.90) lwidth(none) || /// 
				line Yqi10 lnYeari if Yqi10[_n-1]<lnceiling | _n==1, lcolor(black) || ///
				rarea Yqf19 Yqf1 lnYearf, fcolor(black*.10) lwidth(none) || rarea Yqf18 Yqf2 lnYearf, fcolor(black*.20) lwidth(none) || rarea Yqf17 Yqf3 lnYearf, fcolor(black*.30) lwidth(none) || rarea Yqf16 Yqf4 lnYearf, fcolor(black*.40) lwidth(none) || rarea Yqf15 Yqf5 lnYearf, fcolor(black*.50) lwidth(none) || rarea Yqf14 Yqf6 lnYearf, fcolor(black*.60) lwidth(none) || rarea Yqf13 Yqf7 lnYearf, fcolor(black*.70) lwidth(none) || rarea Yqf12 Yqf8 lnYearf, fcolor(black*.80) lwidth(none) || rarea Yqf11 Yqf9 lnYearf, fcolor(black*.90) lwidth(none) || /// 
				line Yqf10 lnYearf if Yqf10 & (Yqi10[_n-1]<lnceiling | _n==1), lcolor(black) || ///
				connected ln${`depvar'var} lnYear if Year>=ti, lcolor(red) mcolor(red) lstyle(solid) lwidth(vthin) msize(tiny) msymbol(O) ///
				legend(off) ///
				ylabel(`ylabels', notick angle(hor)) ytitle("") subtitle(`:char ${`depvar'var}[unit]', pos(11) span size(medsmall)) ///
				xtitle(Years till `=round(`mediantakeoff')') xscale(reverse) ///
				xlabels(`xloglabels') ///
				plotregion(lwidth(none) margin(zero)) graphregion(margin(0 7 0 0)) ///
				name(`name', replace) `graphopt'
			graph save `name', replace
			graph export `name'.png, width(2680) height(1552) replace
      if "`close'"=="" graph close `name'
		}
		else {
			qui forvalues q=1/19 {
				// shift first obs clipped on right to right edge, interpolating to prettify graph
				replace Yqi`q' = lnceiling if (Yqi`q' >= lnceiling & Yqi`q' < .) | (Yqi`q'[_n-1]==lnceiling & Yqi`q'==.)  // for display, censor extreme values 
				replace Yqi`q' = lnfloor   if  Yqi`q' <= lnfloor                 | (Yqi`q'[_n-1]==lnfloor   & Yqi`q'==.)
			}
			myNiceLogLabels Yqi* ln${`depvar'var} if Year>=ti
			local ylabels `r(labels)'
		}

		local name `modelname'Dist`depvar'`tsamplename'`v'`region'
		twoway rarea Yqi19 Yqi1 Yeari, fcolor(black*.10) lwidth(none) || rarea Yqi18 Yqi2 Yeari, fcolor(black*.20) lwidth(none) || rarea Yqi17 Yqi3 Yeari, fcolor(black*.30) lwidth(none) || rarea Yqi16 Yqi4 Yeari, fcolor(black*.40) lwidth(none) || rarea Yqi15 Yqi5 Yeari, fcolor(black*.50) lwidth(none) || rarea Yqi14 Yqi6 Yeari, fcolor(black*.60) lwidth(none) || rarea Yqi13 Yqi7 Yeari, fcolor(black*.70) lwidth(none) || rarea Yqi12 Yqi8 Yeari, fcolor(black*.80) lwidth(none) || rarea Yqi11 Yqi9 Yeari, fcolor(black*.90) lwidth(none) || /// 
			line Yqi10 Yeari if Yqi10 & Yqi10>lnfloor & Yqi10<lnceiling, lcolor(black) legend(off) || ///
			connected ln${`depvar'var} Year if Year>=ti, lcolor(red) mcolor(red) lstyle(solid) lwidth(vthin) msize(tiny) msymbol(O) ///
			ylabel(`ylabels', notick angle(hor)) xlabel(`xlabels', notick) ytitle("") subtitle(`:char ${`depvar'var}[unit]', pos(11) span size(medsmall)) ///
			xtitle(Year) ///
			plotregion(lwidth(none) margin(zero)) graphregion(margin(0 7 0 0)) ///
			name(`name', replace) `graphopt'
		graph save `name', replace
		graph export `name'.png, width(2680) height(1552) replace
    if "`close'"=="" graph close `name'
	restore
	
	drop Yq*
end

global region // blank, Eurasia, Africa, Americas, or Oceania
local region = substr("$region",1,6)

constraint 1 /nu + /gamma = 0  // defines CEV model within Bernoulli model
constraint 2                   // no constraint for general model

***
*** Main estimation code, for 4 outcomes-- GWP, Pop, GWP/cap (GWPcap), and frontier GDP/cap (GDPcap), 
***                           4 samples--starting 1 m BCE ("All") or 10,000 BCE ("12K"), decennial ("Dec") or annual data after 1950
***                           2 or 3 models (Bernoulli, NLS, and for regressions reported in text, time-consuming CEV)
*** For each non-NLS estimate, plots are generated of 99 sample paths ("Paths")
***                               and of distribution of 10,000 paths ("Dist") with time scale in logs ("Log") or not
***                                                                            with distinct parameter draws for each path or not ("nov")
***                               and of quantiles of data points in simulated paths
*** Produces a LOT of graphs in png and Stata gph format. File naming is <Model>[Paths/Dist]<Depvar><Sample>[Log][nov] and <Model>CDF<Depvar><Sample>.
*** In paper: BernouPathsGWPAllnov, BernouDistGWPAllnov, BernouPathsGWP12KDecLog, BernouDistGWP12KDecLog, BernouCDFGWP12KDec, and,
*** in combined figure, BernouDistPop12KDecLog, BernouDistGWPcap12KDecLog, BernouDistGDPcap12KDecLog, BernouCDFPop12KDec, BernouCDFGWPcap12KDec, BernouCDFGDPcap12KDec.
***
*** Also produces Stata estimation results (.ster) files with names <Model><Depvar><Sample>. Versions with "nl" suffixes add nonlinear derived estimates
*** (s, B, delta, sigma, median takeoff times and takeoff probabilities) 
***
*** And produces results tables "SDE fits [depvar] .rtf" like paper's Table 2; and paper's Table 3, "SDE fits 12KDecnl .rtf".
***

forvalues v=1/4 {
	local ests
	local depvar: word `v' of GWP Pop GWPcap GDPcap

	forvalues s=1/4 {
		local tsample: word `s' of 1 "Year<=1950 | inlist(Year,1960,1970,1980,1990,2000,2010,2019)" "Year>=-10000" "Year>=-10000 & (Year<=1950 | inlist(Year,1960,1970,1980,1990,2000,2010,2019))"
		local tsamplename: word `s' of All AllDec 12K 12KDec

		if "`depvar'"=="GDPcap" local xlabels 0 "1 CE" 1e3 "5000 CE" 2e3 "2000 CE" 3e3 "3000 CE" 4e3 "4000 CE" 5e3 "5000 CE"
		                   else local xlabels: word `=2-inlist(`s',1,2)' of `"-1e6 "1,000,000 BCE" -7.5e5 "750,000 BCE" -5e5 "500,000 BCE" -2.5e5 "250,000 BCE" 0 "1 CE" "' ///
		                                                                    `"-1e4 "10000 BCE" -5e3 "5000 BCE" 0 "1 CE" 5e3 "5000 CE" "'

		PrepData if `tsample', depvar(${`depvar'var}) historicpop(PopMcEvedyJones) prehistoricpop(PopDeevey) region($region)
		gen byte NLSsample = _n < _N

    * Kremer-style NLS
    asdf ${`depvar'var}, model(bernounls)
    scalar s     = [/s]
    scalar B     = [/B]
    scalar delta = [/delta]
    scalar lna   = 2 * ([/lnsig] + ln(B)) - 1.62e42fefa39efX-001 /*ln(2)*/
    scalar b     = -B * delta
    scalar nu    = (1 - (s+s) / exp(2*[/lnsig])) / B
    scalar gamma = -1 / B

    constraint 10 /s = s  // for simulation and plotting of NLS results, estimate lnsig by fitting diffusion model with all other parameters fixed at NLS estimates
    constraint 11 /B = B
    constraint 12 /delta = delta
    mat init = s, B, delta, [/lnsig]
    cap noi asdf ${`depvar'var} [aw=HYDEwt], model(bernoudiff2) reflect(`=nu>0') from(init) iter(25) constr(10 11 12)
    scalar lna    = 2*([/lnsig] + ln(B)) - 1.62e42fefa39efX-001 /*ln(2)*/
    scalar nu     = (1 - (s+s) / exp([/lnsig]+[/lnsig])) / B
    scalar mediantakeoff = tf - ln1m(b*Yf^(-B) / exp(lna) / invgammap(-nu, .5))/b
*    Plot, nov modelname(NLS) tsamplename(`tsamplename') depvar(`depvar') region(`region') xlabels(`xlabels') mediantakeoff(mediantakeoff)

		mat NLSinit = lna, b, nu, gamma
		forvalues m=`=1+(`v'!=1 & `s'!=4)'/2 {  // do the time-consuming CEV only for regressions reported in text
			local modelname : word `m' of CEV Bernou
      
			* diffusion model--try estimating from a few starting points, reflecting and non-reflecting variants
      asdf ${`depvar'var} [aw=HYDEwt], model(bernoudiff) reflect(0) from(NLSinit) iter(`=cond(`m'==1,1000,50)') constraint(`m')
      if !e(converged) {
      	mat init = 0,0,-1,-1
      	asdf ${`depvar'var} [aw=HYDEwt], model(bernoudiff) reflect(0) from(init) iter(`=cond(`m'==1,1000,50)') constraint(`m')
      }
      est sto absorb
      if !e(converged) {
        if [/nu] >= -1 mat init = e(b)
                  else mat init = 0, 0, -.5,  1
        asdf ${`depvar'var} [aw=HYDEwt], model(bernoudiff) reflect(1) from(init) iter(2000) constraint(`m')
        if !e(converged) est restore absorb
      }
			est store `modelname'`depvar'`tsamplename'`region'
			est save  `modelname'`depvar'`tsamplename'`region', replace
			scalar converged = e(converged)
			scalar reflect   = e(reflect)

			* compute nonlinear derived quantities including s, B, delta, sigma, no-take-off probability, median take-off year, with standard errors
      scalar s = exp([/lna]) * [/gamma] * ([/nu] + [/gamma])
			local nlcomcmd (lna  : [/lna]) (b:[/b]) (nu:[/nu]) (gamma:[/gamma]) ///
										 (s    : `=cond(e(rank)<e(k), "0", "exp([/lna]) * [/gamma] * ([/nu] + [/gamma]) / s")') ///  // under constrained model, s=0 by fiat and nlcom complains
										 (B    : -1 / [/gamma]) ///
										 (delta: [/b] * [/gamma]) ///
										 (sigma: sqrt(2) * exp([/lna]/2) * abs([/gamma]))
			if ([/nu]+[/gamma]) / [/b] < 0 local nlcomcmd `nlcomcmd' (Y_b  : (-exp([/lna]) * ([/nu]+[/gamma]) / [/b]) ^ [/gamma])   // zero-drift level

			foreach if in i f {  // initial, final
				if [/nu]<=0 & [/gamma]<0 {
					local notakeoffprob gammap(-[/nu], [/b] / exp([/lna]) * Y`if'^(1/[/gamma]))
					scalar notakeoffprob`if' = `notakeoffprob'
					if notakeoffprob`if' local nlcomcmd `nlcomcmd' (notakeoffprob`if': `notakeoffprob' / `=notakeoffprob`if'')  // trick to avoid crash stata.com/statalist/archive/2009-03/msg01244.html
													else local nlcomcmd `nlcomcmd' (notakeoffprob`if': 0                                     )

					local takeoffb0 Y`if'^(1/[/gamma]) / exp([/lna]) / invgammap(-[/nu], .5)
					if notakeoffprob`if'<.5 {
						local mediantakeoff`if' cond([/b], -ln1m([/b]*`takeoffb0')/[/b], `takeoffb0')
						scalar mediantakeoff`if' = `mediantakeoff`if''
						local nlcomcmd `nlcomcmd' (mediantakeoff`if': `mediantakeoff`if'' / mediantakeoff`if')
					}
				}
				else {
					scalar notakeoffprobi = 1
					scalar notakeoffprobf = 1
				}
			}
			nlcom `nlcomcmd', post

			local nlcomcmd
			foreach var in `:colnames e(b)' {
				     if inlist(substr("`var'",1,9), "s", "notakeoff") local nlcomcmd `nlcomcmd' (`var': _b[`var'] * `var'                           )  // finish trick to avoid crash
				else if        substr("`var'",1,13)=="mediantakeoff"  local nlcomcmd `nlcomcmd' (`var': _b[`var'] * `var' + t`=substr("`var'",14,1)')
				else                                                  local nlcomcmd `nlcomcmd' (`var': _b[`var']                                   )
			}
			nlcom `nlcomcmd', post

			cap scalar mediantakeoff = _b[mediantakeofff]
			if _rc scalar mediantakeoff = .

			estadd scalar converged = converged
			estadd scalar reflect   = reflect
			est sto   `modelname'`depvar'`tsamplename'`region'nl
			est resto `modelname'`depvar'`tsamplename'`region'

			* plot sample paths and distribution thereof for full-sample fit
      foreach nov in nov "" {  // with and without variability around parameter estimates
        Plot, `nov' modelname(`modelname') tsamplename(`tsamplename') depvar(`depvar') region(`region') xlabels(`xlabels') mediantakeoff(mediantakeoff)
      }
      
      * goodness-of-fit tests and graph
      cap drop p
      predict double p, cdf
      di "Kolmogorov-Smirnov test that p values from `est' model are uniformly distributed:"
      ksmirnov p = p
      estadd scalar ksmirnovp = r(p), replace: `modelname'`depvar'`tsamplename'`region'nl
      cap drop Lp
      gen double Lp = p[_n-1]
      reg p Lp [aw=HYDEwt]
      test Lp
      estadd scalar corrp = r(p), replace: `modelname'`depvar'`tsamplename'`region'nl
			est restore `modelname'`depvar'`tsamplename'`region'nl
			est save    `modelname'`depvar'`tsamplename'`region'nl, replace

      cap drop lnYear
      gen lnYear = ln(tf + 10 - Year)
      myNiceLogLabels lnYear if p<.
      local xloglabels `r(labels)'
      twoway connected p lnYear if p<. , lwidth(vthin) lcolor(gs8) msize(tiny) ||  ///
        scatter p lnYear if p<. & (Year<1950 | mod(Year,10)==0 | Year==2019), mlab(strYear) mlabgap(0) mlabsize(small) msym(none)  ///
        legend(off) xscale(noline reverse range(1.75 .)) xlab(none) ymtick(0(.01)1, notick grid glwidth(vthin) glcolor(gs15)) ///
        ylab(0 "0" .1 "0.1" .2 "0.2" .3 "0.3" .4 "0.4" .5 "0.5" .6 "0.6" .7 "0.7" .8 "0.8" .9 "0.9" 1 "1", notick angle(hor) format(%3.1f) grid glwidth(thin) glcolor(gs14)) ///
        plotregion(lwidth(none) margin(zero)) graphregion(margin(0 0 0 1)) ytitle(Quantile) xtitle(Year) xscale(range(3.7 8.1)) ///
        name(`modelname'CDF`depvar'`tsamplename', replace)
      graph save "`modelname'CDF`depvar'`tsamplename'", replace
      graph export `modelname'CDF`depvar'`tsamplename'.png, replace width(2680) height(1552)
		}
    
    * LR test of CEV vs Bernoulli diffusion, when CEV done
		cap lrtest Bernou`depvar'`tsamplename'`region' CEV`depvar'`tsamplename'`region', df(1)
		if !_rc {
    	estadd scalar chi2CEV  = r(chi2), replace: Bernou`depvar'`tsamplename'`region'nl
      estadd scalar chi2pCEV = r(p)   , replace: Bernou`depvar'`tsamplename'`region'nl
			est save Bernou`depvar'`tsamplename'`region'nl, replace
    }
		local ests `ests' Bernou`depvar'`tsamplename'`region'nl
	}
  
  * Estimation table for each depvar (only that for GWP in text, Table 2)
	esttab `ests' using "SDE fits `depvar' `region'.rtf", replace ///
		keep (lna b nu gamma s B delta sigma Y_b notakeoffprobi mediantakeoffi notakeoffprobf mediantakeofff) ///
		order(lna b nu gamma s B delta sigma Y_b notakeoffprobi mediantakeoffi notakeoffprobf mediantakeofff) ///
		scalars(reflect converged chi2CEV chi2pCEV ksmirnovp corrp N) noobs ///
		se nostar varlabels(delta \u0948? sigma \u0963? gamma \u0947? nu \u0957? lna "log a")
}

* Cross-var table for sample starting 10,000 BCE, decennial after 1950 (Table 3)
esttab BernouGWP12KDec`region'nl BernouPop12KDec`region'nl BernouGWPcap12KDec`region'nl BernouGDPcap12KDec`region'nl using "SDE fits 12KDecnl `region'.rtf", replace ///
	keep (lna b nu gamma s B delta sigma Y_b notakeoffprobi mediantakeoffi notakeoffprobf mediantakeofff) ///
	order(lna b nu gamma s B delta sigma Y_b notakeoffprobi mediantakeoffi notakeoffprobf mediantakeofff) ///
	scalars(reflect converged chi2CEV chi2pCEV ksmirnovp corrp N) noobs ///
	se nostar varlabels(delta \u0948? sigma \u0963? gamma \u0947? nu \u0957? lna "log a") ///
	mtitles(GWP Pop GWPcap GDPcap) nonumbers

* combined figure in paper for Pop, GWPcap, and GDPcap
graph combine BernouDistPop12KDecLog BernouDistGWPcap12KDecLog BernouDistGDPcap12KDecLog BernouCDFPop12KDec BernouCDFGWPcap12KDec BernouCDFGDPcap12KDec, ///
  cols(2) colfirst xsize(7) ysize(8.5) graphregion(margin(0 2 0 2))
graph export Combined.png, replace


***
*** Robustness test: Fit to De Long and Hanson series, as reported in footnote ~41
***

mat init = 0,0,-1,-1
local e 0
forvalues s=1/2 {
  local tsample: word `s' of 1 Year>=-10000
  local tsamplename: word `s' of All90Dec 12K90Dec
  foreach GWPvar in GWP GWPDeLong GWPHanson {
    PrepData if `tsample' & (Year<=1950 | inlist(Year,1960,1970,1980,1990,2000)), depvar(`GWPvar') historicpop(PopMcEvedyJones) prehistoricpop(PopDeevey)
    asdf `GWPvar' [aw=HYDEwt], model(bernoudiff) reflect(0) from(init)
    scalar s = exp([/lna]) * [/gamma] * ([/nu] + [/gamma])
    nlcom (lna:[/lna]) (b:[/b]) (nu:[/nu]) (gamma:[/gamma]) ///
                       (s    : exp([/lna]) * [/gamma] * ([/nu] + [/gamma]) / s)  ///
                       (B    : -1 / [/gamma]) ///
                       (delta: [/b] * [/gamma]) ///
                       (sigma: sqrt(2) * exp([/lna]/2) * abs([/gamma])), post
    eststo est`++e': nlcom (lna:_b[lna]) (b:_b[b]) (nu:_b[nu]) (gamma:_b[gamma]) (s:_b[s]*s) (B:_b[B]) (delta:_b[delta])  (sigma:_b[sigma]), post
  }
}
esttab est? using "vsDeLongHanson.rtf", se replace


***
*** Robustness test: reduce or increase GWP by HYDE-based factor, as reported in text
***

PrepData if Year>=-10000 & (Year<=1950 | inlist(Year,1960,1970,1980,1990,2000,2010,2019)), depvar(GWP) historicpop(PopMcEvedyJones) prehistoricpop(PopDeevey)
est use BernouGWP12KDec
estadd scalar B = -1/[/gamma]
estadd scalar mediantakeoff = Year[`=_N'] - ln1m([/b] * GWP[`=_N'] ^ (1/[/gamma]) / exp([/lna]) / invgammap(-[/nu], .5)) / [/b]
mat b = e(b)
cap drop _GWP
gen double _GWP = GWP * exp(-HYDEsd)
eststo lo: asdf _GWP [aw=HYDEwt], model(bernoudiff) reflect(0) from(b)
estadd scalar B = -1/[/gamma]
estadd scalar mediantakeoff = Year[`=_N'] - ln1m([/b] * _GWP[`=_N'] ^ (1/[/gamma]) / exp([/lna]) / invgammap(-[/nu], .5)) / [/b]
replace _GWP = GWP * exp(HYDEsd)
eststo hi: asdf _GWP [aw=HYDEwt], model(bernoudiff) reflect(0) from(b)
estadd scalar B = -1/[/gamma]
estadd scalar mediantakeoff = Year[`=_N'] - ln1m([/b] * _GWP[`=_N'] ^ (1/[/gamma]) / exp([/lna]) / invgammap(-[/nu], .5)) / [/b]
est table BernouGWP12KDec lo hi, stat(B mediantakeoff)


***
*** Compute quantile of each obs in predicted distribution from regression on previous ones
*** among the saved graphs, BernoudiffPredGWP12KDec used in paper
***

forvalues v=1/4 {
	local depvar: word `v' of GWP Pop GWPcap GDPcap

	forvalues s=4/4 {
		local tsample: word `s' of 1 "Year<=1950 | inlist(Year,1960,1970,1980,1990,2000,2010,2019)" "Year>=-10000" "Year>=-10000 & (Year<=1950 | inlist(Year,1960,1970,1980,1990,2000,2010,2019))"
		local tsamplename: word `s' of All AllDec 12K 12KDec

		PrepData if `tsample', depvar(${`depvar'var}) historicpop(PopMcEvedyJones) prehistoricpop(PopDeevey) region($region)

		* get quantiles of predictions of next observation conditional on previous; incorporates modelling uncertainty via simulation
		foreach est in nls diff diffstat {
      cap drop ptile_`est'
      gen double ptile_`est' = .
      mata st_view(ptile_`est'=., ., "ptile_`est'")
    }
    mata Y = st_data(., "${`depvar'var}"); tdelta = st_data(., "tdelta"); Ydot = (Y[|2\.|] :/ Y[|.\rows(Y)-1|]) :^ (1 :/ tdelta[|2\.|]) :- 1
		gen byte sample = _n < _N

    est use Bernou`depvar'`tsamplename'
    mat init = e(b)

    forvalues n=`=_N'(-1)6 {	
			* Kremer-style NLS
			asdf GWP if _n<`n', model(bernounls)
      mata b = st_matrix("e(b)"); symeigensystem(st_matrix("e(V)"), X=., L=.)
      mata bsim = b :+ rnormal(10000, 4, 0, 1) * (X :* sqrt(edittozero(L,1)))'
      mata ptile_nls[`n'] = mean(normal((Ydot[`n'-1] :- bsim[,1] :* Y[`n'-1] :^ bsim[,2] :- bsim[,3]) :/ (exp(bsim[,4]) :/ tdelta[`n'])))
    
			* ...and diffusion model, dynamic
			cap noi asdf ${`depvar'var} if _n<`n' [aw=HYDEwt], model(bernoudiff) iter(50) reflect(0) from(init)  // tight iter limit because if it takes many, probably unstable
			if !_rc & e(converged) {
        mat init = e(b)
        asdfSim, t0(`=Year[`n'-1]') y0(`=${`depvar'var}[`n'-1]') t(`=tdelta[`n']') tsamp(1) tres(10000) m(10000)
        mata ptile_diff[`n'] = mean(Yf :< Y[`n'])
      }

/*			* ...and diffusion model, static
      cap noi asdf ${`depvar'var} if _n<`n' [aw=HYDEwt*tdelta], model(bernoudiff) reflect(0) from(NLSinit) iter(100) modeltype(static)
      if  !_rc | e(converged)==0 {
        if !_rc {
          est sto best
          scalar bestll = e(ll)
        }
        else scalar bestll = .
        cap noi asdf ${`depvar'var} if _n<`n' [aw=HYDEwt*tdelta], model(bernoudiff) reflect(0) from(init) iter(100) modeltype(static)
        if _rc | e(ll) < bestll est restore best
      }
      asdfSim, t0(`=Year') y0(`=${`depvar'var}') t(`=Year[`n']-Year') tsamp(1) tres(10000) m(10000)  // simulate even if no convergence...
      mata ptile_diffstat[`n'] = mean(Yf() :< Y[`n'])*/
    }
    
    est use Bernou`depvar'`tsamplename'`region'nl
    cap drop lnYear
    gen double lnYear = ln(tf + 10 - Year)
    format ptile_* %3.2f
    foreach est in nls diff /*diffstat*/ {
      di "Kolmogorov-Smirnov test that p values from `est' model are uniformly distributed:"
      ksmirnov ptile_`est'=ptile_`est'
      twoway connected ptile_`est' lnYear if ptile_`est'<. , lwidth(vthin) lcolor(gs8) msize(tiny) || ///
        scatter ptile_`est' lnYear if ptile_`est'<. & (Year<1950 | mod(Year,10)==0 | Year==2019), mlab(strYear) mlabgap(0) mlabsize(medium) msym(none) ///
        legend(off) xscale(noline reverse range(1.75 .)) xlab(none) ymtick(0(.01)1, notick grid glwidth(vthin) glcolor(gs15)) ///
        ylab(0 "0" .1 "0.1" .2 "0.2" .3 "0.3" .4 "0.4" .5 "0.5" .6 "0.6" .7 "0.7" .8 "0.8" .9 "0.9" 1 "1", notick angle(hor) format(%3.1f) grid glwidth(thin) glcolor(gs14)) ///
        plotregion(lwidth(none) margin(0 0 1 0)) graphregion(margin(0 0 0 3)) ytitle("") xtitle(Year) ///
        name(Bernou`est'Pred`depvar'`tsamplename', replace)
      graph save "Bernou`est'Pred`depvar'`tsamplename'", replace
      graph export Bernou`est'Pred`depvar'`tsamplename'.png, replace width(2680) height(1552)
    }
  }
}



***
*** adapt graphs for blog post -- runs way faster if you shrink the graph window to nothing
***

set scheme s1mono
grstyle init blogscheme, replace
grstyle linewidth plineplot vthin
grstyle color heading    "68 68 68"
grstyle gsize heading    medium
grstyle color axis_title "68 68 68"
grstyle color tick       "68 68 68"
grstyle color tickline   "68 68 68"
grstyle color tick_label "68 68 68"
grstyle yesno title_span yes
grstyle yesno alt_yaxes yes
graph set window fontface Merriweather

graph use BernouPathsGWP12KDecnov, scheme(blogscheme)
gr_edit .style.editstyle margin(10 0 0 0) editcopy
gr_edit .subtitle.DragBy -.5948275285189022 108.5560239546959
// gr_edit .plotregion1.plot100.style.editstyle line(color("48 155 166")) editcopy
gr_edit .plotregion1.plot100.style.editstyle line(width(medium)) editcopy
// gr_edit .plotregion1.plot100.style.editstyle marker(fillcolor("48 155 166")) editcopy
// gr_edit .plotregion1.plot100.style.editstyle marker(linestyle(color("48 155 166"))) editcopy
gr_edit .yaxis1.style.editstyle linestyle(pattern(blank)) editcopy
gr_edit .yaxis1.style.editstyle majorstyle(tickstyle(show_ticks(yes))) editcopy
gr_edit .title.text.Arrpush "Sample paths from Bernoulli diffusion model for GWP"
forvalue p=21/99 {
  gr_edit .plotregion1.plot`p'.style.editstyle line(color(none)) editcopy
}
graph save BernouPathsGWP12KDecnovBlog, replace
graph export BernouPathsGWP12KDecnovBlog.png, replace width(1440)

graph use BernouDistGWP12KDecLog, scheme(blogscheme)
gr_edit .style.editstyle margin(5 0 0 0) editcopy
gr_edit .subtitle.DragBy -.5948275285189022 108.5560239546959
// gr_edit .plotregion1.plot21.style.editstyle line(color("48 155 166")) editcopy
gr_edit .plotregion1.plot21.style.editstyle line(width(medium)) editcopy
// gr_edit .plotregion1.plot21.style.editstyle marker(fillcolor("48 155 166")) editcopy
// gr_edit .plotregion1.plot21.style.editstyle marker(linestyle(color("48 155 166"))) editcopy
gr_edit .yaxis1.style.editstyle linestyle(pattern(blank)) editcopy
gr_edit .yaxis1.style.editstyle majorstyle(tickstyle(show_ticks(yes))) editcopy
gr_edit .xaxis1.title.style.editstyle margin(0 0 0 1) editcopy
gr_edit .title.text.Arrpush "Bernoulli diffusion model for GWP, incorporating"
gr_edit .title.text.Arrpush "modeled stochasticity and modeling uncertainty"
graph save BernouDistGWP12KDecLogBlog, replace
graph export BernouDistGWP12KDecLogBlog.png, replace width(1440)

graph use BernouDiffPredGWP12KDec, scheme(blogscheme)
gr_edit .style.editstyle margin(1 0 0 0) editcopy
gr_edit .plotregion1.plot1.style.editstyle line(color("48 155 166")) editcopy
gr_edit .plotregion1.plot1.style.editstyle line(width(medium)) editcopy
gr_edit .plotregion1.plot1.style.editstyle marker(linestyle(color("48 155 166"))) editcopy
gr_edit .plotregion1.plot1.style.editstyle marker(fillcolor("48 155 166")) editcopy
gr_edit .yaxis1.edit_tick 0 0 "0%", tickset(major)
gr_edit .yaxis1.edit_tick 1 0.1 "10%", tickset(major)
gr_edit .yaxis1.edit_tick 2 0.2 "20%", tickset(major)
gr_edit .yaxis1.edit_tick 3 0.3 "30%", tickset(major)
gr_edit .yaxis1.edit_tick 4 0.4 "40%", tickset(major)
gr_edit .yaxis1.edit_tick 5 0.5 "50%", tickset(major)
gr_edit .yaxis1.edit_tick 6 0.6 "60%", tickset(major)
gr_edit .yaxis1.edit_tick 7 0.7 "70%", tickset(major)
gr_edit .yaxis1.edit_tick 8 0.8 "80%", tickset(major)
gr_edit .yaxis1.edit_tick 9 0.9 "90%", tickset(major)
gr_edit .yaxis1.edit_tick 10 1 "100%", tickset(major)
gr_edit .title.style.editstyle margin(small) editcopy
gr_edit .title.text.Arrpush Percentile of GWP in distribution when model fit to previous data
graph export BernouDiffPredGWP12KDecBlog.png, replace width(1440)


***
*** Monte Carlo test of Bernoulli diffusion model and Bernoulli NLS
***

cap program drop sim
program define sim, rclass
	drop _all
	est restore fit
  asdfSim, t0(-10000) t(25000) y0(1.6) tsamp(250000) tres(1) m(1) nq(0) ns(1) nov
	mata keep = ceil(Years * (colmax(selectindex(Y:<100000)) / 12019)); keep[1] = 1  // draw from range < $100 trillion
	mata Y = Y[keep]; t = t[keep]
	getmata Y t, force replace double

	cap noi asdf Y t, model(bernounls)
	if _rc {
		return scalar convergednls = 0
		mat init = 0,0,-1,-1
	}
	else {
    mat init = 2*([/lnsig] + ln([/B])) - 1.62e42fefa39efX-001 , -[/B]*[/delta] , (1 - 2*[/s] / exp(2*[/lnsig])) / [/B] , -1/[/B]
		return scalar Bnls = [/B]
		return scalar convergednls = 1
	}

  cap noi asdf Y t, model(bernoudiff) reflect(`=init[1,3]>0') from(init) iter(100)
  if _rc return scalar convergeddiff = 0
  else {
    return scalar Bdiff = -1/[/gamma]
    return scalar convergeddiff = e(converged)
  }
  return scalar icdiff = e(ic)
end

PrepData if Year>=-10000 & (Year<=1950 | inlist(Year,1960,1970,1980,1990,2000,2010,2019)), depvar(GWP) historicpop(PopMcEvedyJones) prehistoricpop(PopDeevey)
mat init = 0,0,-1,-1
asdf GWP [aw=HYDEwt], model(bernoudiff) reflect(0) from(init)
est sto fit
mata Years = st_data(., "Year") :+ 10000
drop _all

simulate Bnls=r(Bnls) convergednls=r(convergednls) Bdiff=r(Bdiff) convergeddiff=r(convergeddiff), reps(5000) seed(2394857): sim
save "Diffusion vs NLS simulation results", replace

sum if convergeddiff  // summary stats for simualtions in which Bernoulli fit converged
est resto fit
scalar B = -1/[/gamma]
mean Bnls if convergednls
di "NLS bias = " _b[Bnls] - B
test _b[Bnls] = `=B'
mean Bdiff if convergeddiff
di "Bernoulli ML bias = " _b[Bdiff] - B
test _b[Bdiff] = `=B'

twoway kdensity Bnls if abs(Bnls-B)<=.4, lcolor(gray) || kdensity Bdiff if abs(Bdiff-B)<=.4, lcolor(gray) xline(`=B') legend(off) ///
  plotregion(lwidth(none) margin(zero)) graphregion(margin(3 3 0 0)) yscale(off range(0 9)) xtitle(Estimate of {it:B}) name(DiffvNLSMonteCarlo, replace) ///
  xlabels(.3 "0.3" .4 "0.4" .5 "0.5" .6 "0.6" .7 "0.7" .8 "0.8" .9 "0.9") ///
  text(5 .61 "Bernoulli diffusion/" "maximum likelihood", just(left) place(e)) ///
  text(2.5 .28 "Nonlinear least" "squares", just(left) place(e)) ///
  text(8.5 .555 "True value", just(left) place(e))
graph save DiffvNLSMonteCarlo, replace
graph export DiffvNLSMonteCarlo.png, width(2680) height(1552) replace


***
*** replicate a couple of Kremer tables
*** doesn't include D-W test because it's not formally meaningful with unevenly spaced observations
***

PrepData if Year<=1990, depvar(Pop) historicpop(PopKremer) prehistoricpop(PopKremer) gwpstart(.)
replace PopgrowthKremer = PopgrowthKremer * 100
replace PopKremer = PopKremer / 1000
eststo II_1: reg PopgrowthKremer PopKremer
eststo II_2: reg PopgrowthKremer PopKremer if Year>=-200
eststo II_3: reg PopgrowthKremer PopKremer               [aw=Year[_n+1]-Year]
eststo II_4: reg PopgrowthKremer PopKremer if Year>=-200 [aw=Year[_n+1]-Year]
esttab II_?, r2 obslast se nostar

eststo VI_1: nl (PopgrowthKremer = {delta=4.51e4} + {s=.493   } * PopKremer ^ {B=1.03}) [aw=Year[_n+1]-Year]
eststo VI_2: nl (PopgrowthKremer = {delta=6.25e4} + {s=.507   } * PopKremer ^ {B=1.22}) [aw=Year[_n+1]-Year] if Year<=1960
eststo VI_3: nl (PopgrowthKremer = {delta=-.038 } + {s=1.18e-9} * PopKremer ^ {B=1.43})
eststo VI_4: nl (PopgrowthKremer = {delta=-.036 } + {s=2.13e-6} * PopKremer ^ {B=.907})                      if Year<=1960
esttab VI_?, r2 obslast se nostar
