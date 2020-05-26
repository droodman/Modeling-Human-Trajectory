*** Simulation code for Roodman, "Modeling the Human Trajectory" as of May 22, 2020
*** Requires Ben Jann's grstyle and colorpalette package
*** Generates graphs based on simulations of multivariate system, with and without endogenous natural resources

set scheme s1color
grstyle init
colorpalette Set1, nograph
grstyle set color "`r(p5)'" "`r(p4)'" "`r(p2)'" "`r(p7)'" "`r(p3)'" "`r(p1)'"  // Use ColorBrewer's Set1 palette, but save red for Y and green for R and avoid yellow
graph set window fontface Cambria


* simulation engine
* B, s, delta are the objects defined in the paper
* ly0 is the starting value of y, in logs
* dt is the simulation time increment, and steps the number of steps
* Tres is the number of steps per time point in the return value
* Returns a matrix of path values, with steps/Tres+1 rows and length(s) columns

cap mata mata drop sim()
mata
mata set matastrict on
mata set mataoptimize on
mata set matalnum off
real matrix sim(real matrix B, real colvector s, real colvector delta, real colvector ly0, real scalar dt, real scalar steps, real scalar Tres) {
	real scalar i, j, repsteps; real matrix BT, lny; real rowvector sTdt, lnyi, deltaTdt
	
	repsteps = steps/Tres + 1  // reported steps

	sTdt = s'dt
	BT = B'
	(lny = J(repsteps, rows(s), 0))[1,] = lnyi = ly0'
	deltaTdt = delta'dt
	for (j=2; j<=repsteps; j++) {
		for (i=Tres; i; i--)
			lnyi = lnyi + ln1p(deltaTdt + sTdt :* exp(lnyi * BT))
		lny[j,] = lnyi
	}
	return(lny)
}
end


foreach r in 0 1 { // flag to endogenize natural resource

	mata alpha = .3 \ .3 \ .3 \ .1 \ 1                  // exponents in production for K, P, H, R, A
	mata phi = 0 \ -.1 \ .1 \ -.01 \ .4                 // elasticities w.r.t. technology
	mata s = .25 \ .2 \ .04 \ -.01*`r' \ .025           // reinvestment coefficients
	mata delta = -.03 \ -.02 \ -.02 \ .001*`r' \ -.001  // exogenous depreciation/appreciation
	mata (y0 = J(5,1,1))[(1\3\5)] = J(3,1,.10374)       // initial values

	mata k = length(phi); iota = J(k,1,1); (Phi = J(k,k,0))[,k] = phi; B = iota*alpha' - I(k) + Phi
	mata eigensystem(B, V=., L=.)
	mata "eigenvectors"
	mata V
	mata "eigenvalues"
	mata L

	mata dt = 2e-5
	mata steps = 1e8
	mata Tres = 1000  // steps between each returned state

	mata lny1 = sim(B, s, delta, ln(y0), dt, steps, Tres)
	mata lnY  = lny1 * alpha
	mata T = (dt*Tres) * (0::rows(lny1))
	mata theta = acos(abs((lny1 * Re(V[,1])) :/ sqrt(quadrowsum(lny1:*lny1,1))))  // angle between state vector and largest eigenvalue
	mata z = s' :* exp(lny1 * B') :+ delta'

	clear
	getmata t=T ly1_0=lnY (ly1_*)=lny1 (z*)=z theta, double force replace

	if !`r' {
		mata (y0 = J(5,1,1))[(1\3\5)] = J(3,1,.10364)
		mata lny2 = sim(B, s, delta, ln(y0), dt, steps, Tres); lnY2  = lny2 * alpha
		getmata ly2_0=lnY2 (ly2_*)=lny2, double force replace
	}

	* simulate explosion at finer time scale
	mata tr   = max((1::rows(z)) :* (z[,1]:<1))  // transition point to higher-time-res simulation
	mata tinf = max((1::rows(z)) :* (rownonmissing(lny1):>0))
	mata ly0 = lny1[tr,]'
	mata steps = (T[tinf]-T[tr]+1)/dt * 1.01
	mata lny1 = lny1[|.,.\tr-1,.|] \ sim(B, s, delta, ly0, dt, steps, Tres=1)  // tack on higher-time-res simulation near explosion
	mata T    = T   [|.,.\tr-1,.|] \ ((dt*Tres) * (0::steps) :+ T[tr])

	mata z = s' :* exp(lny1 * B') :+ delta'
	mata zdot = (z :- delta') :* z*B'
	mata lnY  = lny1 * alpha
	mata Z    = z * alpha
	mata Zdot = zdot * alpha

	getmata _t=T (ly1_*)=lny1 ly1_0=lnY (z*)=z z0=Z (zdot*)=zdot zdot0=Zdot theta, double force replace  // 0 index for output
	drop if ly1_1==. & ly1_2==. & ly1_3==. & ly1_4==. & ly1_5==.

	foreach var of varlist z? zdot? {
		cap drop l`var'
		gen double l`var' = ln(`var')
	}

	forvalues i=0/5 {
		cap noi regress lz`i' ly1_`i' if lz`i'>5
	}

	if `r' {
		replace ly1_0 = -6 if (ly1_0==. & ly1_4[_n-1] < 0) | ly1_4 < -6 // censor plunging output for graphing
		replace ly1_4 = -6 if (ly1_4==. & ly1_4[_n-1] < 0) | ly1_4 < -6 // censor plunging resources for graphing
	}

  line ly1_1 _t if ly1_1 < 5      , lwidth(medthick) || ///
  line ly1_2 _t if ly1_2 < 5      , lwidth(medthick) || ///
  line ly1_3 _t if ly1_3 < 5      , lwidth(medthick) || ///
  line ly1_5 _t if ly1_5 < 5      , lwidth(medthick) || ///
  line ly1_4 _t if ly1_4 < 5 & `r', lwidth(medthick) || ///
  line ly1_0 _t if ly1_0 < 5      , lwidth(medthick)    ///
    `=cond(`r', "", "|| line ly2_1-ly2_3 ly2_5 ly2_0 t, lstyle(p1 p2 p3 p4 p6) lwidth(medthick..)")' ///
     plotregion(margin(zero) lwidth(none)) graphregion(margin(0 5 0 3)) ytitle("Level (in natural logarithms)") ///
     ylab(-5 "{&minus}5" 0 "0" 5 "5", angle(hor)) ///
     xtitle(Years) ///
     legend(order(1 2 3 4 `=cond(`r',"5","")' 6) rows(1) region(lwidth(none)) bmargin(zero) margin(zero) span symxsize(*.2) keygap(*.2) size(medsmall) lab(1 Physical capital ({it:K})) lab(2 Population ({it:P})) lab(3 Human capital ({it:H})) lab(4 Technology ({it:A})) lab(5 Resources ({it:R})) lab(6 Output ({it:Y}))) ///
     name(lny_v_t_r`r', replace)
	graph export lny_v_t_r`r'.png, width(2680) height(1552) replace

  * version of graph for blog
  graph set window fontface Merriweather
	line ly1_1 _t if ly1_1 < 5      , lwidth(medthick) || ///
  line ly1_2 _t if ly1_2 < 5      , lwidth(medthick) || ///
  line ly1_3 _t if ly1_3 < 5      , lwidth(medthick) || ///
  line ly1_5 _t if ly1_5 < 5      , lwidth(medthick) || ///
  line ly1_4 _t if ly1_4 < 5 & `r', lwidth(medthick) || ///
  line ly1_0 _t if ly1_0 < 5      , lwidth(medthick)    ///
     plotregion(margin(zero) lwidth(none)) graphregion(margin(5 0 0 2)) ytitle("") ///
     ylab(`=ln(.01)' "0.01" `=ln(.1)' "0.1" 0 "1"  `=ln(10)' "10"  `=ln(100)' "100", angle(hor) labgap(small) labcolor("68 68 68")) ///
     xlab(0(500)2000, labcolor("68 68 68")) xtitle(Years, margin(tiny) color("68 68 68")) ///
     yscale(alt noline lcolor("68 68 68")) xscale(lcolor("68 68 68")) ///
     legend(order(1 2 3 4 `=cond(`r',"5","")' 6) rows(1) color("68 68 68") region(lwidth(none)) bmargin(zero) margin(zero) span symxsize(*.2) keygap(*.2) size(medsmall) lab(1 Capital) lab(2 Labor) lab(3 Human capital) lab(4 Technology) lab(5 Resources) lab(6 GWP)) ///
     name(lny_v_t_r`r'Blog, replace)
	graph export lny_v_t_r`r'Blog.png, width(2680) height(1552) replace
  graph set window fontface Cambria

  cap drop n
  gen long n = _n
  sum n if z1[_n-1]<=0 & z1>0 & z1<., meanonly  // cut out early history, start when technology growth becomes permanently positive
  if r(N) {
    keep if _n >= r(max)
    replace n = _n
  }

  forvalues i=0/5 {
    sum n if lz`i'<25 & ly1_`i'<., meanonly
    local x`i' = ly1_`i'[r(max)]
    local y`i' = lz`i'[r(max)]
  }
  twoway connected lz1 ly1_1 if lz1>-10 & lz1<25      , msize(tiny) lwidth(vthin) || ///
         connected lz2 ly1_2 if lz2>-10 & lz2<25      , msize(tiny) lwidth(vthin) || ///
         connected lz3 ly1_3 if lz3>-10 & lz3<25      , msize(tiny) lwidth(vthin) || ///
         connected lz5 ly1_5 if lz5>-10 & lz5<25      , msize(tiny) lwidth(vthin) || ///
         connected lz4 ly1_4 if lz5>-10 & lz5<25 & `r', msize(tiny) lwidth(vthin) || ///
         connected lz0 ly1_0 if lz0>-10 & lz0<25      , msize(tiny) lwidth(vthin)    ///
         , name(lnz_v_lny_r`r', replace) ///
         plotregion(margin(zero) lwidth(none)) graphregion(margin(0 3 0 3))  legend(off) ///
         ylab(-10 "{&minus}10" 0 "0" 10 "10" 20 "20", angle(hor)) ///
         xtitle("Level (in natural logarithms)") ytitle("Growth rate (in natural logarithms)") ///
         text(`y1' `x1' "{it:K}" `y2' `x2' "{it:P}" `y3' `x3' "{it:H}" `y5' `x5' "{it:A}" `y0' `x0' "{it:Y}" `=cond(`r',`"`y4' `x4' "{it:R}""',"")', place(n) size(medsmall) margin(tiny))
  graph export lnz_v_lny_r`r'.png, width(2680) height(1552) replace

  forvalues i=0/5 {
    sum n if lz`i'<. & lzdot`i'<-5, meanonly
    local x`i' = lz`i'[r(max)]
    local y`i' = lzdot`i'[r(max)]
  }
  twoway connected lzdot1 lz1 if lzdot1>-20 & lz1>-10 & lzdot1<-5      , msize(tiny) lwidth(vthin) || ///
         connected lzdot2 lz2 if lzdot2>-20 & lz2>-10 & lzdot2<-5      , msize(tiny) lwidth(vthin) || ///
         connected lzdot3 lz3 if lzdot3>-20 & lz3>-10 & lzdot3<-5      , msize(tiny) lwidth(vthin) || ///
         connected lzdot5 lz5 if lzdot5>-20 & lz5>-10 & lzdot5<-5      , msize(tiny) lwidth(vthin) || ///
         connected lzdot5 lz5 if lzdot5>-20 & lz5>-10 & lzdot5<-5 & `r', msize(tiny) lwidth(vthin) || ///
         connected lzdot0 lz0 if lzdot0>-20 & lz0>-10 & lzdot0<-5      , msize(tiny) lwidth(vthin)    ///
         , name(lnzdot_v_lnz_r`r', replace) ///
         plotregion(margin(zero) lwidth(none)) graphregion(margin(0 3 0 3)) ///
         ylab(-15  "{&minus}20" -10  "{&minus}10" -5  "{&minus}5", angle(hor)) ///
         xlab(-10  "{&minus}10" -5  "{&minus}5") legend(off) ///
         xtitle("Growth rate (in natural logarithms)") ytitle("Time derivative of growth rate (in natural logarithms)") ///
         text(`y1' `x1' "{it:K}" `y2' `x2' "{it:P}" `y3' `x3' "{it:H}" `y5' `x5' "{it:A}" `y0' `x0' "{it:Y}" `=cond(`r',`"`y4' `x4' "{it:R}""',"")', place(n) size(medsmall) margin(tiny)) 
  graph export lnzdot_v_lnz_r`r'.png, width(2680) height(1552) replace
}
