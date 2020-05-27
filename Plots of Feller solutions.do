*** Code for Roodman, "Modeling the Human Trajectory" as of May 22, 2020.
*** Plot the Feller and noncentral chi2 distributions (in main text) and the corresponding diffusions (appendix)
*** The code for the latter comes in two versions, one relying on R because Stata was crashing when trying to 
*** use the viridis color palette with so many gradations.
***
*** In Stata, requires Ben Jann's grstyle, colrspace, and palettes packages and E.F. Haghish's rcall
*** In R, requires libraries ggplot2, viridis, gridExtra, and latex2exp

cap mata mata drop fstar()
cap mata mata drop f()

mata
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

real matrix fstar(real colvector X, real colvector X0, real rowvector t, real scalar a, real scalar b, real scalar nu, real scalar reflect) {
	real matrix retval; real colvector x, lambda; real scalar i, j, _x, _lambda; real rowvector C
	C = (b? b :/ expm1(b*t) : 1:/t) / a
	x = X * C
	lambda = X0 * (C :* exp(b:*t))
	retval = J(rows(X), cols(t), .)
	for (i=rows(X);i;i--)
		for (j=cols(t);j;j--) {
			_x      =      x[min((i,rows(     x))), min((j,cols(     x)))]
			_lambda = lambda[min((i,rows(lambda))), min((j,cols(lambda)))]
			retval[i,j] = lnPDFFeller(nu, ln(_lambda), ln(_x), reflect)
		}
	return (exp(retval) :* C)
}

// crude implementation that works even when diffusion is negative
real colvector f(real colvector lambda, real colvector x, real scalar nu, real scalar reflect) {
	real colvector m; complex colvector m1, m2
  m = 0..1e3; m1 = reflect? m : C(m:-nu); m2 = reflect? C(m:+nu) : m
  return (quadrowsum(Re(exp(ln(lambda)*m1 :+ ln(x)*m2 :- lngamma(m1:+1) :- lngamma(m2:+1)))) :/ exp(lambda) :/ exp(x))
}
end


*** sample curves for distribution (figure in text)
set scheme s1color
grstyle init
grstyle set linewidth .5mm
grstyle set color hue, n(`=(3.0 - -3.0)/.5+1') dir(-1)
* global n05color `r(p`=3/.5')' // snag color for sticky distribution nu = -.5

mata lambda=1;x=rangen(0.0001,exp(-5/3),100) \ exp((-500::1000)/300)
drop _all
mata p = J(rows(x), (3.0 - -3.0)/.5+1, .)
forvalues r=0/1 {
	mata c = 0
	mata for (nu=-3; nu<=3; nu=nu+.5) p[,++c] = f(lambda, x, nu, `r')
	getmata x (pn3_`r' pn25_`r' pn2_`r' pn15_`r' pn1_`r' pn05_`r' p0_`r' p05_`r' p1_`r' p15_`r' p2_`r' p25_`r' p3_`r') = p, force replace
}

/*scalar mu = 1
scalar X0 = 1
scalar t = 1
gen double rho = 2 * mu * (1 - normal((X+X0+2*mu*t)/sqrt(t))) / normalden(X+X0+2*mu*t, sqrt(t))
gen double pn05_s = pn05_0 + rho * (pn05_1 - pn05_0)*/

forvalues r=0/1 {
	qui foreach var of varlist p* {
		replace `var' = . if `var' < -.1 | `var' > .4
	}
}

line p*_0 x, lwidth(thin..) || line p*_1 x, lpat(shortdash..) lwidth(thin..) pstyle(p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13) /*|| line p*_s x, lcolor("$n05color") lpat(longdash..) lwidth(thin..)*/ ///
  || if x<11, ///
	plotregion(lwidth(none) margin(zero)) graphregion(margin(0 0 0 3)) ///
	ylabel(-.1  "{fontface Cambria:{&minus}0.1}" 0  "{fontface Cambria:0}" .1 "{fontface Cambria:0.1}" .2 "{fontface Cambria:0.2}" .3 "{fontface Cambria:0.3}" .4 "{fontface Cambria:0.4}", labsize(medium) angle(hor)) ///
	xlabel(0 "{fontface Cambria:0}" 5 "{fontface Cambria:5}" 10 "{fontface Cambria:10}", labsize(medium)) ///
	xtitle("{it:{fontface Cambria:x}}") ///
	legend(off) xscale(lcolor(none)) yscale(lcolor(none)) ///
	text(0 1.5 "{fontface Cambria:{it:{&nu}} = {&minus}3}") text(0.15 6.5 "{fontface Cambria:{it:{&nu}} = +3}") ///
	name(distsamples, replace)
graph export dist_samples.png, width(2680) height(1552) replace



*** plot diffusions in Stata

set scheme s1color
grstyle init
grstyle set linewidth .1mm

colorpalette hcl, n(180) viridis nograph
local ccolors `"`r(p)'"'

global Xstep = 25 * 4/exp(1)  // assure X includes 4 exactly, max graphed
mata a=1; X0=1; X=exp((200::-500)/$Xstep); t=(1..51)/50:-.019999  // X ordered from high to low to facilitate computation of medians below
local g 0
foreach b in -1 1 {
	forvalues nu=-1.5/.5 {
		forvalues r=`=`nu'>0'/`=`nu'>-1' {
			local hasmode = `r'==0 | `nu'>=0
			cap drop _all
			mata C = (`b'? `b' :/ expm1(`b'*t) : 1:/t) / a
			mata lambda = X0 * exp(`b'*t) :* C
			mata mean = 1 :+ `nu' :+ lambda
			mata if (!`r') mean = mean :* gammap(-(`nu'), lambda) :+ lambda :* gammaden(-(`nu'),1,0,lambda);;
			mata p = (t \ mean:/C)'#J(length(X),1,1), J(length(t),1,X), vec(fstar(X, X0, t, a, `b', `nu', `r'))
			getmata (t mean X p) = p, force replace
			gsort t -X  // data already sorted this way; tell Stata

			gen A = X*p/$Xstep  // area of bar in Riemann integral
			bysort t: gen runsumA = sum(A)  // running sum, from high to low X
			gen byte medial = runsumA[_n-1] < .5 & (runsumA >= .5 | _n==_N)
			gen median = cond(_n==_N, 0, (.5 - runsumA[_n-1]) / A * (X - X[_n-1]) + X[_n-1]) if medial // interpolate median estimate, for smoothing

			if `hasmode' {
				by t: egen peak = max(p)
				by t: gen byte modal = p==peak
				// quadratic interpolation of mode using 3 points, for smoothing, unless highest point is lowest X value
				by t: gen mode = cond(_n==_N, X, .5 * (X[_n-1]^2*p + X^2*p[_n+1] + X[_n+1]^2*p[_n-1] - X[_n-1]^2*p[_n+1] - X^2*p[_n-1] - X[_n+1]^2*p) / (X[_n-1]*p + X*p[_n+1] + X[_n+1]*p[_n-1] - X[_n-1]*p[_n+1] - X*p[_n-1] - X[_n+1]*p)) if modal
			}
			else gen byte modal = 0

			if `g'==0 {  // use same color cuts in all graphs
				local ccuts 0(.00625)1
				sum p if p, meanonly
				forvalues i=1/19 {
					local ccuts	`ccuts' `=r(max)^(`i'/20)'
				}
			}

			local centerleft = `nu'==-.5 & !`r'
			twoway contour p X t if X<=4.001, ccuts(`ccuts') legend(off) clegend(off) interp(none) /*ccolors(`ccolors')*/  ///
				xtitle(`=cond(`centerleft', "{fontface Cambria:{it:t}}", `"" ""')', margin(zero) size(medlarge))  ///
				ytitle(`=cond(`r'==1-(`nu'==.5), `"" ""', `"{fontface Cambria:{it:{&nu}} = `=cond(`nu'<0,"{&minus}","+")+string(abs(`nu'),"%3.1f")' ({it:c} = `=cond(`nu'<-1,"{&minus}","+")+string(abs(`nu'+1),"%3.1f")')}"')', size(medlarge) margin(zero)) ///
				xscale(noline) yscale(noline) ztitle("")  ///
				ylabel(1 "{fontface Cambria:1}" 2 "{fontface Cambria:2}" 3 "{fontface Cambria:3}" 4 "{fontface Cambria:4}", angle(hor) `=cond(`centerleft',"","labcolor(white)")') ///
				xlabel(.2 "{fontface Cambria:0.2}" .4 "{fontface Cambria:0.4}" .6 "{fontface Cambria:0.6}" .8 "{fontface Cambria:0.8}" 1 "{fontface Cambria:1}", `=cond(`centerleft', "", "labcolor(white)")') ///
				plotregion(lwidth(none) margin(zero)) graphregion(margin(0 3 0 0)) ||  ///
				line mean   t if X==X[1] &  mean<4.001 & t>.01, lcolor(gray) ||  ///
				line median t if medial & median<4.001 & t>.01, lcolor(gray)  ///
				`=cond(`hasmode', "|| line mode t if modal & t>.01, lcolor(gray) ","")' ///
				name(contour`++g', replace)  ///
				title(`=cond(`nu'==-1.5+`r', `"{fontface Cambria:`=cond(`r', "Noncentral {it:{&chi}{superscript:2}}", "Feller")' diffusion}"', `"" ""')', margin(zero) size(medlarge)) ///
				nodraw
		}
	}
}
graph combine contour1 contour2 contour3 contour4, holes(2 5) cols(2) xsize(7.5) ysize(8.67) imargin(tiny) name(diffusions_bn1, replace)  // figure for -1
graph export diffusions_bn1.png, width(720) replace
graph combine contour5 contour6 contour7 contour8, holes(2 5) cols(2) xsize(7.5) ysize(8.67) imargin(tiny) name(diffusions_bp1, replace)  // figure for +1
graph export diffusions_bp1.png, width(720) replace

*** plot diffusions using R for graphics

rcall clear
rcall: library(ggplot2); library(viridis); library(gridExtra); library(latex2exp)
global Xstep 50
mata a=1; t=(1..201)/200:-1/200*.99; X0=1; X=(2000::1)/$Xstep  // X ordered from high to low to facilitate computation of medians below
local g 0
foreach b in -1 1 {
	forvalues nu=-1.5/.5 {
		forvalues r=`=`nu'>0'/`=`nu'>-1' {
			local hasmode = `r'==0 | `nu'>=0
			cap drop _all
			mata C = (`b'? `b' :/ expm1(`b'*t) : 1:/t) / a
			mata lambda = X0 * exp(`b'*t) :* C
			mata mean = 1 :+ `nu' :+ lambda
			mata if (!`r') mean = mean :* gammap(-(`nu'), lambda) :+ lambda :* gammaden(-(`nu'),1,0,lambda);;
			mata p = (t \ mean:/C)'#J(length(X),1,1), J(length(t),1,X), vec(fstar(X, X0, t, a, `b', `nu', `r'))
			getmata (t mean X p) = p, force replace
			gsort t -X  // data already sorted this way; tell Stata

			gen double A = p/$Xstep  // area of bar in Riemann integral
			bysort t: gen double runsumA = sum(A)  // running sum, from high to low X
			gen byte medial = runsumA[_n-1] < .5 & (runsumA >= .5 | _n==_N)
			gen double _median = cond(_n==_N, 0, (.5 - runsumA[_n-1]) / A * (X - X[_n-1]) + X[_n-1]) if medial  // interpolate median estimate, for smoothing
			by t: egen double median = min(_median)
			keep if X<=4
			drop _median

			if `hasmode' {
				by t: egen double peak = max(p)
				by t: gen byte modal = p==peak
				// quadratic interpolation of mode using 3 points, for smoothing, unless highest point is lowest X value
				by t: gen double _mode = cond(_n==_N, X, .5 * (X[_n-1]^2*p + X^2*p[_n+1] + X[_n+1]^2*p[_n-1] - X[_n-1]^2*p[_n+1] - X^2*p[_n-1] - X[_n+1]^2*p) / (X[_n-1]*p + X*p[_n+1] + X[_n+1]*p[_n-1] - X[_n-1]*p[_n+1] - X*p[_n-1] - X[_n+1]*p)) if modal
				by t: egen double mode = min(_mode)
				drop _mode
			}
			else gen byte modal = 0

			if `g++'==0 {
				sum p if p, meanonly
				scalar maxp = r(max)
			}
			replace p = p^(1/9/ln(maxp)) if p > 1  // trim extremes of p for coloring -- for p > 1 give log p same range width as 0 - 1/9.

			rcall: mydata <- subset(st.data(), X==X[1])  // one entry per time value, for graphing mean, median, mode
			rcall: g`g' <- ggplot(st.data(), aes(x = t, y = X)) + ///
			                      geom_raster(interpolate=TRUE, aes(fill = p)) + ///
							              scale_fill_viridis(option="viridis") + ///
							              scale_x_continuous(breaks=seq(.2,1,.2), expand=c(0, 0), name="t") + ///
							              scale_y_continuous(breaks=seq(1,4,1)  , expand=c(0, 0)) + ///
									          theme(legend.position = "none", ///
														      plot.title = element_text(hjust = 0.5), ///
										              axis.title.x = element_text(face="italic"), ///
											            axis.title.y = element_text(face="italic")) + ///
										        coord_cartesian(ylim=c(0,4)) + ///
										        geom_line(data=mydata, aes(x=t, y=mean  ), color="gray75", size=.25) + ///
										        geom_line(data=mydata, aes(x=t, y=median), color="gray75", size=.25) + ///
														labs(title = ifelse(`nu'==.5, " ", TeX("`=cond(`nu'==-.5, "Reflecting boundary", "Absorbing boundary")'"))) ///
			                      `=cond(`hasmode', `"+ geom_line(data=mydata, aes(x=t, y=mode), color="gray75", size=.25)"', "")'; ///
			       if (`nu'==-.5 & !`r')  g`g' <- g`g' + theme(plot.title=element_text(color="white")) ///
			         else                    g`g' <- g`g' + theme(axis.text.x=element_text(color="white"), axis.title.x=element_text(color="white"), axis.text.y=element_text(color="white")); ///
			       if (`r'==1-(`nu'==.5)) g`g' <- g`g' + ylab(" ") ///
			         else                    g`g' <- g`g' + ylab(TeX("$\\nu = `=string(`nu',"%3.1f")'$"))
			rcall: ggsave("diffusions_g`g'.png", g`g', width=7.5, height=4.5)
		}       
	}
}
rcall: ggsave("diffusions_bn1R.png", grid.arrange(g1, g2, g3, g4, layout_matrix=rbind(c(1,NA),c(2,3),c(NA,4))), width=7.5, height=8.48)
rcall: ggsave("diffusions_bp1R.png", grid.arrange(g5, g6, g7, g8, layout_matrix=rbind(c(1,NA),c(2,3),c(NA,4))), width=7.5, height=8.48)


----
* Cross-sections as t->\infty
mata

delta = -.0003903
s = .0003668
B = .4267431
sigma2 = .0042699^2
lna    = ln(sigma2 * B*B) - 1.62e42fefa39efX-001 /*ln(2)*/
a = exp(lna)
b      = -B * delta
nu     = (1 - (s+s) / sigma2) / B
gamma  = -1 / B

t = 1000000
Y0 = 1
Y = exp((-10000::10000)/10000)

ttilde = (1 - exp(-b*t)) / b
X  = Y:^-gamma
X0 = Y0^-gamma
Z = exp(-b*t) * X / a
Z0 = X0 / a

p = exp(ln(abs(B)) :+ ln(Y):*(-B-1) :- ln(ttilde) :+ lnPDFFeller(nu, ln(Z0/ttilde), ln(Z/ttilde), 0))

end
getmata p Y, force replace double
line p Y