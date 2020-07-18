*
* After "superexponential" diffusion ML fitting, compute nonlinear derived quantities including s, B, delta, sigma, no-take-off probability, median take-off year
* Uses Stata's nlcom command and *replaces* current estimation results
*

cap program drop asdfSuperNL
program define asdfSuperNL
  tempname s notakeoffprobi notakeoffprobf mediantakeoffi mediantakeofff ti tf Yi Yf
  
  scalar `ti' = e(ti)
  scalar `tf' = e(tf)
  scalar `Yi' = e(depvari)
  scalar `Yf' = e(depvarf)

  scalar `s' = exp([/lna]) * [/gamma] * ([/nu] + [/gamma])
  local nlcomcmd (lna  : [/lna]) (b:[/b]) (nu:[/nu]) (gamma:[/gamma]) ///
                 (s    : `=cond(e(rank)<e(k), "0", "exp([/lna]) * [/gamma] * ([/nu] + [/gamma]) / `s'")') ///  // under constrained model, s=0 by fiat and nlcom complains
                 (B    : -1 / [/gamma]) ///
                 (delta: [/b] * [/gamma]) ///
                 (sigma: sqrt(2) * exp([/lna]/2) * abs([/gamma]))
  if ([/nu]+[/gamma]) / [/b] < 0 local nlcomcmd `nlcomcmd' (Y_b  : (-exp([/lna]) * ([/nu]+[/gamma]) / [/b]) ^ [/gamma])   // zero-drift level
	if [/gamma] local nlcomcmd `nlcomcmd' (phi_A : [/gamma] - 1 / [/gamma])

  foreach if in i f {  // initial, final
    if [/nu] <= 0 & [/gamma] < 0 {
     local notakeoffprobeq gammap(-[/nu], [/b] / exp([/lna]) * `Y`if''^(1/[/gamma]))
      scalar `notakeoffprob`if'' = `notakeoffprobeq'
      if `notakeoffprob`if'' local nlcomcmd `nlcomcmd' (notakeoffprob`if': `notakeoffprobeq' / `notakeoffprob`if'')  // trick to avoid crash stata.com/statalist/archive/2009-03/msg01244.html
                        else local nlcomcmd `nlcomcmd' (notakeoffprob`if': 0                                      )

      local takeoffb0 `Y`if''^(1/[/gamma]) / exp([/lna]) / invgammap(-[/nu], .5)
      if `notakeoffprob`if''<.5 {
        local mediantakeoffeq cond([/b], -ln1m([/b]*`takeoffb0')/[/b], `takeoffb0')
        scalar `mediantakeoff`if'' = `mediantakeoffeq'
        local nlcomcmd `nlcomcmd' (mediantakeoff`if': `mediantakeoffeq' / `mediantakeoff`if'')
      }
    }
    else {
      scalar `notakeoffprobi' = 1
      scalar `notakeoffprobf' = 1
    }
  }
  qui nlcom `nlcomcmd', post

  local nlcomcmd
  foreach var in `:colnames e(b)' {
         if inlist(substr("`var'",1,9), "s", "notakeoff") local nlcomcmd `nlcomcmd' (`var': _b[`var'] * ``var''                             )  // finish trick to avoid crash
    else if        substr("`var'",1,13)=="mediantakeoff"  local nlcomcmd `nlcomcmd' (`var': _b[`var'] * ``var'' + `t`=substr("`var'",14,1)'')
    else                                                  local nlcomcmd `nlcomcmd' (`var': _b[`var']                                       )
  }
  qui nlcom `nlcomcmd', post
  nlcom
end
