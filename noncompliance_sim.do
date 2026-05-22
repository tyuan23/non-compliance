/******************************************************************************
  Simulation study: Comparison of RPSFTM and IPCW methods for treatment
                    non-compliance in survival analysis

  Version : 17  
******************************************************************************/

clear all
set more off


/*============================================================================
 SET PARAMETERS
  --------------------------------------------------------------------------
============================================================================*/

global N_OBS        ""     
global N_REPS       ""     
global SEED         ""    
global LAMBDA_NS    ""     
global LAMBDA_SEL   ""     
global NC_INTERCEPT ""      
global NEWV_COEF    ""      
global TRT_COEF     ""      
global NODE_COEF    ""     
global TUMOR_COEF   ""     
global AFT_COEF     ""     
global GAMMA_S      ""     


/*============================================================================
                    PART A. NO SELECTION BIAS SCENARIO
============================================================================*/

/*----------------------------------------------------------------------------
  Program 1 of 6: RPSFTM (Cox and Weibull), no selection bias
----------------------------------------------------------------------------*/

capture program drop simstudy_rpsftm
program define simstudy_rpsftm, rclass

    syntax [, obs(int $N_OBS)             ///
              trtandcoef(real $TRT_COEF) ///
              nodecoef(real $NODE_COEF)           ///
              tumorcoef(real $TUMOR_COEF)          ///
              aftcoef(real $AFT_COEF)           ///
              lambda(real $LAMBDA_NS)       ///
              gammas(real $GAMMA_S)             ///
            ]

    clear
    pr drop _all

    /*--- 1. Generate baseline covariates and survival times ---*/

    rbinary age toxicity new_v tumor_new node_new,                         ///
        means(0.33, 0.48, 0.45, 0.58, 0.28)                                ///
        corr(1    , 0.4 , 0.68, 0   , 0    \                               ///
             0.4  , 1   , 0.68, 0   , 0    \                               ///
             0.68 , 0.68, 1   , 0   , 0    \                               ///
             0    , 0   , 0   , 1   , 0.5  \                               ///
             0    , 0   , 0   , 0.5 , 1)                                   ///
        n($N_OBS)

    set obs `obs'
    gen trtrand = rbinomial(1, 0.5)

    survsim timeDFS_y, distribution(weibull) lambda(`lambda') gammas(`gammas') ///
        cov(trtrand `trtandcoef' node_new `nodecoef' tumor_new `tumorcoef')

    /*--- Administrative censoring at 6 years ---*/

    gen id      = _n
    gen admin   = 6
    gen admin_w = admin * 52
    gen timeDFS_w = timeDFS_y * 52
    replace timeDFS_w = round(timeDFS_w)
    replace timeDFS_w = 1 if timeDFS_w == 0
    replace admin_w   = round(admin_w)

    gen dead = 1 if timeDFS_w <  admin_w
    replace dead = 0 if timeDFS_w >= admin_w
    replace timeDFS_w = admin_w if timeDFS_w > admin_w

    /*--- Non-compliance indicator (CRT arm only) ---*/

    gen p     = invlogit($NC_INTERCEPT + 0.6 * age + 0.7 * toxicity) if trtrand == 1
    gen noncp = rbinomial(1, p)                            if trtrand == 1
    replace noncp = 0 if noncp == .
    drop p

    /*--- Time of first non-compliance (in weeks) ---*/

    gen noncp_time1 = rbinomial(1, 0.24) if noncp == 1
    gen noncp_time2 = rbinomial(1, 0.34) if noncp == 1 & noncp_time1 == 0
    gen noncp_time3 = rbinomial(1, 0.6)  if noncp == 1 & noncp_time1 == 0 & noncp_time2 == 0
    gen noncp_time  = 1  if noncp_time1 == 1
    replace noncp_time = 4  if noncp_time2 == 1
    replace noncp_time = 11 if noncp_time3 == 1 & noncp == 1
    replace noncp_time = 7  if noncp_time == . & noncp == 1
    replace noncp_time = .  if noncp_time > timeDFS_w & trt == 1 & noncp == 1
    replace noncp      = 0  if noncp_time == . & trt == 1

    /*--- 2. New (shorter) survival times for non-compliers ---*/

    gen xoOSgainobs = timeDFS_w - noncp_time if noncp == 1
    replace xoOSgainobs = round(xoOSgainobs * `aftcoef')
    gen timeDFS_w2 = cond(noncp == 1, xoOSgainobs + noncp_time, timeDFS_w)

    gen died = 0
    replace died = 1 if timeDFS_w2 < admin_w & dead == 1
    replace timeDFS_w2 = admin_w if timeDFS_w2 > admin_w
    replace xoOSgainobs = . if noncp_time >  admin_w
    replace noncp       = . if noncp_time >= admin_w
    replace noncp_time  = . if noncp_time >= admin_w

    /*--- 3. True dRMST, rRMST and MSD from Weibull ---*/

    scalar lambda = `lambda'
    scalar gamma  = `gammas'
    scalar beta1  = `trtandcoef'

    gen timeDFS_y1 = timeDFS_y
    replace timeDFS_y1 = 6 if timeDFS_y1 > 6
    stset timeDFS_y1, failure(dead = 1)
    sort _t
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    gen H0      = lambda * _t^(gamma)
    gen H1      = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen au_S0_gp = lag * S0_true
    gen au_S1_gp = lag * S1_true

    foreach i in 3 {
        egen le_S0_gp`i' = sum(au_S0_gp) if _t <= `i'
        egen le_S1_gp`i' = sum(au_S1_gp) if _t <= `i'
        summarize le_S0_gp`i'
        return scalar le0_gp`i' = r(mean)
        summarize le_S1_gp`i'
        return scalar le1_gp`i' = r(mean)
    }

    * Difference in survival at 3 years
    drop H0 H1 S1_true S0_true
    gen H0 = lambda * 3^(gamma)
    gen H1 = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen S_diff_true = S1_true - S0_true
    sum S_diff_true
    return scalar S_diff_true = r(mean)

    /*=== RPSFTM (Cox-based test, default) ===*/

    replace noncp      = 0          if noncp == .
    replace noncp_time = timeDFS_w2 if noncp == 0
    replace noncp_time = timeDFS_w2 if trt   == 0

    * Hazard ratio
    stset timeDFS_w2, failure(died == 1)
    strbee trtrand, xo1(noncp_time noncp) adj(node_new tumor_new) ///
        endstudy(admin_w) gen(c1)
    strbee, hr

    return scalar rpsftm_hr        = r(HR_adj)
    return scalar rpsftm_hr_LB     = r(HR_adj_low)
    return scalar rpsftm_hr_UB     = r(HR_adj_upp)
    return scalar rpsftm_aft       = exp(r(psi))
    return scalar rpsftm_aft_LB    = exp(r(psi_low))
    return scalar rpsftm_aft_UB    = exp(r(psi_upp))
    return scalar itt_statistics   = r(Z_ITT)

    gen cov_rpsftm_hr = 0
    replace cov_rpsftm_hr = 1 if exp(`trtandcoef') > return(rpsftm_hr_LB) & exp(`trtandcoef') < return(rpsftm_hr_UB)
    return scalar cov_rpsftm_hr = cov_rpsftm_hr[1]

    gen cov_rpsftm_aft = 0
    replace cov_rpsftm_aft = 1 if `aftcoef' > return(rpsftm_aft_LB) & `aftcoef' < return(rpsftm_aft_UB)
    return scalar cov_rpsftm_aft = cov_rpsftm_aft[1]

    gen power_rpsftm_hr = 0
    replace power_rpsftm_hr = 1 if 1 > return(rpsftm_hr_LB) & 1 < return(rpsftm_hr_UB)
    return scalar power_rpsftm_hr = power_rpsftm_hr[1]

    gen power_rpsftm_aft = 0
    replace power_rpsftm_aft = 1 if 1 > return(rpsftm_aft_LB) & 1 < return(rpsftm_aft_UB)
    return scalar power_rpsftm_aft = power_rpsftm_aft[1]

    /*--- RMST and dSurvival, WITH recensoring ---*/

    gen newtime = .
    replace newtime = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft) if noncp == 1
    replace newtime = timeDFS_w2 if newtime == .
    replace newtime = admin_w    if trt == 1 & newtime > admin_w

    gen dd = dc1  if trt == 1
    replace dd = died if trt == 0
    gen newtime1 = newtime / 52
    stset newtime1, failure(dd == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr

    foreach i in 3 {
        predict S0_rpsftm1, basesurv
        gen S1_rpsftm1     = S0_rpsftm1^exp(_b[trt])
        gen au_S0_rpsftm1  = lag * S0_rpsftm1
        gen au_S1_rpsftm1  = lag * S1_rpsftm1
        egen le_S0_rpsftm1`i' = sum(au_S0_rpsftm1) if _t <= `i'
        egen le_S1_rpsftm1`i' = sum(au_S1_rpsftm1) if _t <= `i'
        summarize le_S0_rpsftm1`i'
        return scalar le0_rpsftm1`i' = r(mean)
        summarize le_S1_rpsftm1`i'
        return scalar le1_rpsftm1`i' = r(mean)

        sort _t
        gen t1 = `i' - _t
        replace t1 = 999999 if _t > `i'
        egen min_t1 = min(t1)
        sum S0_rpsftm1 if t1 == min_t1
        return scalar S0_rpsftm1 = r(mean)
        sum S1_rpsftm1 if t1 == min_t1
        return scalar S1_rpsftm1 = r(mean)
    }

    /*--- RMST and dSurvival, WITHOUT recensoring ---*/

    gen newtime_wr = .
    replace newtime_wr = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft) if noncp == 1
    replace newtime_wr = timeDFS_w2 if newtime_wr == .
    gen newtime1_wr = newtime_wr / 52
    stset newtime1_wr, failure(died == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr
    return scalar rpsftm_hr_nore    = exp(_b[trtrand])
    return scalar rpsftm_hr_nore_SE = _b[trtrand] / return(itt_statistics)
    return scalar rpsftm_hr_nore_LB = exp(_b[trtrand] - 1.96 * return(rpsftm_hr_nore_SE))
    return scalar rpsftm_hr_nore_UB = exp(_b[trtrand] + 1.96 * return(rpsftm_hr_nore_SE))

    gen cov_rpsftm_hr_nore = 0
    replace cov_rpsftm_hr_nore = 1 if exp(`trtandcoef') > return(rpsftm_hr_nore_LB) & exp(`trtandcoef') < return(rpsftm_hr_nore_UB)
    return scalar cov_rpsftm_hr_nore = cov_rpsftm_hr_nore[1]

    gen power_rpsftm_hr_nore = 0
    replace power_rpsftm_hr_nore = 1 if 1 > return(rpsftm_hr_nore_LB) & 1 < return(rpsftm_hr_nore_UB)
    return scalar power_rpsftm_hr_nore = power_rpsftm_hr_nore[1]

    foreach i in 3 {
        predict S0_rpsftm1_wr, basesurv
        gen S1_rpsftm1_wr     = S0_rpsftm1_wr^exp(_b[trt])
        gen au_S0_rpsftm1_wr  = lag * S0_rpsftm1_wr
        gen au_S1_rpsftm1_wr  = lag * S1_rpsftm1_wr
        egen le_S0_rpsftm1`i'_wr = sum(au_S0_rpsftm1_wr) if _t <= `i'
        egen le_S1_rpsftm1`i'_wr = sum(au_S1_rpsftm1_wr) if _t <= `i'
        summarize le_S0_rpsftm1`i'_wr
        return scalar le0_rpsftm1`i'_wr = r(mean)
        summarize le_S1_rpsftm1`i'_wr
        return scalar le1_rpsftm1`i'_wr = r(mean)

        sort _t
        gen t1_wr = `i' - _t
        replace t1_wr = 999999 if _t > `i'
        egen min_t1_wr = min(t1_wr)
        sum S0_rpsftm1_wr if t1_wr == min_t1_wr
        return scalar S0_rpsftm1_wr = r(mean)
        sum S1_rpsftm1_wr if t1_wr == min_t1_wr
        return scalar S1_rpsftm1_wr = r(mean)
    }

    /*=== RPSFTM (Weibull-based test) ===*/

    drop cc1 dd dc1 c1
    stset timeDFS_w2, failure(died == 1)
    strbee trtrand, xo1(noncp_time noncp) adj(node_new tumor_new) ///
        endstudy(admin_w) gen(c1) test(weibull)
    strbee, hr

    return scalar rpsftm_hr_w        = r(HR_adj)
    return scalar rpsftm_hr_LB_w     = r(HR_adj_low)
    return scalar rpsftm_hr_UB_w     = r(HR_adj_upp)
    return scalar rpsftm_aft_w       = exp(r(psi))
    return scalar rpsftm_aft_LB_w    = exp(r(psi_low))
    return scalar rpsftm_aft_UB_w    = exp(r(psi_upp))
    return scalar itt_statistics_w   = r(Z_ITT)

    gen cov_rpsftm_hr_w = 0
    replace cov_rpsftm_hr_w = 1 if exp(`trtandcoef') > return(rpsftm_hr_LB_w) & exp(`trtandcoef') < return(rpsftm_hr_UB_w)
    return scalar cov_rpsftm_hr_w = cov_rpsftm_hr_w[1]

    gen cov_rpsftm_aft_w = 0
    replace cov_rpsftm_aft_w = 1 if `aftcoef' > return(rpsftm_aft_LB_w) & `aftcoef' < return(rpsftm_aft_UB_w)
    return scalar cov_rpsftm_aft_w = cov_rpsftm_aft_w[1]

    gen power_rpsftm_hr_w = 0
    replace power_rpsftm_hr_w = 1 if 1 > return(rpsftm_hr_LB_w) & 1 < return(rpsftm_hr_UB_w)
    return scalar power_rpsftm_hr_w = power_rpsftm_hr_w[1]

    gen power_rpsftm_aft_w = 0
    replace power_rpsftm_aft_w = 1 if 1 > return(rpsftm_aft_LB_w) & 1 < return(rpsftm_aft_UB_w)
    return scalar power_rpsftm_aft_w = power_rpsftm_aft_w[1]

    /*--- RMST and dSurvival, Weibull, WITH recensoring ---*/

    gen newtime_w = .
    replace newtime_w = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft_w) if noncp == 1
    replace newtime_w = timeDFS_w2 if newtime_w == .
    replace newtime_w = admin_w    if trt == 1 & newtime_w > admin_w

    gen dd = dc1  if trt == 1
    replace dd = died if trt == 0
    gen newtime1_w = newtime_w / 52
    stset newtime1_w, failure(dd == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr

    foreach i in 3 {
        predict S0_rpsftm1_w, basesurv
        gen S1_rpsftm1_w    = S0_rpsftm1_w^exp(_b[trt])
        gen au_S0_rpsftm1_w = lag * S0_rpsftm1_w
        gen au_S1_rpsftm1_w = lag * S1_rpsftm1_w
        egen le_S0_rpsftm1`i'_w = sum(au_S0_rpsftm1_w) if _t <= `i'
        egen le_S1_rpsftm1`i'_w = sum(au_S1_rpsftm1_w) if _t <= `i'
        summarize le_S0_rpsftm1`i'_w
        return scalar le0_rpsftm1`i'_w = r(mean)
        summarize le_S1_rpsftm1`i'_w
        return scalar le1_rpsftm1`i'_w = r(mean)

        sort _t
        gen t1_w = `i' - _t
        replace t1_w = 999999 if _t > `i'
        egen min_t1_w = min(t1_w)
        sum S0_rpsftm1_w if t1_w == min_t1_w
        return scalar S0_rpsftm1_w = r(mean)
        sum S1_rpsftm1_w if t1_w == min_t1_w
        return scalar S1_rpsftm1_w = r(mean)
    }

    /*--- RMST and dSurvival, Weibull, WITHOUT recensoring ---*/

    gen newtime_wr_w = .
    replace newtime_wr_w = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft_w) if noncp == 1
    replace newtime_wr_w = timeDFS_w2 if newtime_wr_w == .
    gen newtime1_wr_w = newtime_wr_w / 52
    stset newtime1_wr_w, failure(died == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr
    return scalar rpsftm_hr_nore_wr_w    = exp(_b[trtrand])
    return scalar rpsftm_hr_nore_SE_wr_w = _b[trtrand] / return(itt_statistics_w)
    return scalar rpsftm_hr_nore_LB_wr_w = exp(_b[trtrand] - 1.96 * return(rpsftm_hr_nore_SE_wr_w))
    return scalar rpsftm_hr_nore_UB_wr_w = exp(_b[trtrand] + 1.96 * return(rpsftm_hr_nore_SE_wr_w))

    gen cov_rpsftm_hr_nore_wr_w = 0
    replace cov_rpsftm_hr_nore_wr_w = 1 if exp(`trtandcoef') > return(rpsftm_hr_nore_LB_wr_w) & exp(`trtandcoef') < return(rpsftm_hr_nore_UB_wr_w)
    return scalar cov_rpsftm_hr_nore_wr_w = cov_rpsftm_hr_nore_wr_w[1]

    gen power_rpsftm_hr_nore_wr_w = 0
    replace power_rpsftm_hr_nore_wr_w = 1 if 1 > return(rpsftm_hr_nore_LB_wr_w) & 1 < return(rpsftm_hr_nore_UB_wr_w)
    return scalar power_rpsftm_hr_nore_wr_w = power_rpsftm_hr_nore_wr_w[1]

    foreach i in 3 {
        predict S0_rpsftm1_wr_w, basesurv
        gen S1_rpsftm1_wr_w    = S0_rpsftm1_wr_w^exp(_b[trt])
        gen au_S0_rpsftm1_wr_w = lag * S0_rpsftm1_wr_w
        gen au_S1_rpsftm1_wr_w = lag * S1_rpsftm1_wr_w
        egen le_S0_rpsftm1`i'_wr_wr_w = sum(au_S0_rpsftm1_wr_w) if _t <= `i'
        egen le_S1_rpsftm1`i'_wr_w    = sum(au_S1_rpsftm1_wr_w) if _t <= `i'
        summarize le_S0_rpsftm1`i'_wr_w
        return scalar le0_rpsftm1`i'_wr_w = r(mean)
        summarize le_S1_rpsftm1`i'_wr_w
        return scalar le1_rpsftm1`i'_wr_w = r(mean)

        sort _t
        gen t1_wr_w = `i' - _t
        replace t1_wr_w = 999999 if _t > `i'
        egen min_t1_wr_w = min(t1_wr_w)
        sum S0_rpsftm1_wr_w if t1_wr_w == min_t1_wr_w
        return scalar S0_rpsftm1_wr_w = r(mean)
        sum S1_rpsftm1_w   if t1_wr_w == min_t1_wr_w
        return scalar S1_rpsftm1_wr_w = r(mean)
    }

end

simulate                                                                         ///
    le0_gp3             = r(le0_gp3)              le1_gp3             = r(le1_gp3)              ///
    S_diff_true         = r(S_diff_true)                                                        ///
    rpsftm_hr           = r(rpsftm_hr)            rpsftm_aft          = r(rpsftm_aft)           ///
    cov_rpsftm_hr       = r(cov_rpsftm_hr)        cov_rpsftm_aft      = r(cov_rpsftm_aft)       ///
    power_rpsftm_hr     = r(power_rpsftm_hr)      power_rpsftm_aft    = r(power_rpsftm_aft)     ///
    le0_rpsftm13        = r(le0_rpsftm13)         le1_rpsftm13        = r(le1_rpsftm13)         ///
    S0_rpsftm1          = r(S0_rpsftm1)           S1_rpsftm1          = r(S1_rpsftm1)           ///
    rpsftm_hr_nore      = r(rpsftm_hr_nore)                                                     ///
    cov_rpsftm_hr_nore  = r(cov_rpsftm_hr_nore)   power_rpsftm_hr_nore= r(power_rpsftm_hr_nore) ///
    le0_rpsftm13_wr     = r(le0_rpsftm13_wr)      le1_rpsftm13_wr     = r(le1_rpsftm13_wr)      ///
    S0_rpsftm1_wr       = r(S0_rpsftm1_wr)        S1_rpsftm1_wr       = r(S1_rpsftm1_wr)        ///
    rpsftm_hr_w         = r(rpsftm_hr_w)          rpsftm_aft_w        = r(rpsftm_aft_w)         ///
    cov_rpsftm_hr_w     = r(cov_rpsftm_hr_w)      cov_rpsftm_aft_w    = r(cov_rpsftm_aft_w)     ///
    power_rpsftm_hr_w   = r(power_rpsftm_hr_w)    power_rpsftm_aft_w  = r(power_rpsftm_aft_w)   ///
    le0_rpsftm13_w      = r(le0_rpsftm13_w)       le1_rpsftm13_w      = r(le1_rpsftm13_w)       ///
    S0_rpsftm1_w        = r(S0_rpsftm1_w)         S1_rpsftm1_w        = r(S1_rpsftm1_w)         ///
    rpsftm_hr_nore_wr_w = r(rpsftm_hr_nore_wr_w)                                                ///
    cov_rpsftm_hr_nore_wr_w   = r(cov_rpsftm_hr_nore_wr_w)                                      ///
    power_rpsftm_hr_nore_wr_w = r(power_rpsftm_hr_nore_wr_w)                                    ///
    le0_rpsftm13_wr_w   = r(le0_rpsftm13_wr_w)    le1_rpsftm13_wr_w   = r(le1_rpsftm13_wr_w)    ///
    S0_rpsftm1_wr_w     = r(S0_rpsftm1_wr_w)      S1_rpsftm1_wr_w     = r(S1_rpsftm1_wr_w),     ///
    reps($N_REPS) saving(simstudy_rpsftm, replace) seed($SEED):                     ///
    simstudy_rpsftm


/*----------------------------------------------------------------------------
  Program 2 of 6: IPCW with logistic weights, no selection bias
----------------------------------------------------------------------------*/

capture program drop simstudy_ipcw
program define simstudy_ipcw, rclass

    syntax [, obs(int $N_OBS)             ///
              trtandcoef(real $TRT_COEF) ///
              nodecoef(real $NODE_COEF)           ///
              tumorcoef(real $TUMOR_COEF)          ///
              aftcoef(real $AFT_COEF)           ///
              lambda(real $LAMBDA_NS)       ///
              gammas(real $GAMMA_S)             ///
            ]

    clear
    pr drop _all

    /*--- 1. Generate baseline covariates and survival times ---*/

    rbinary age toxicity new_v tumor_new node_new,                         ///
        means(0.33, 0.48, 0.45, 0.58, 0.28)                                ///
        corr(1    , 0.4 , 0.68, 0   , 0    \                               ///
             0.4  , 1   , 0.68, 0   , 0    \                               ///
             0.68 , 0.68, 1   , 0   , 0    \                               ///
             0    , 0   , 0   , 1   , 0.5  \                               ///
             0    , 0   , 0   , 0.5 , 1)                                   ///
        n($N_OBS)

    set obs `obs'
    gen trtrand = rbinomial(1, 0.5)

    survsim timeDFS_y, distribution(weibull) lambda(`lambda') gammas(`gammas') ///
        cov(trtrand `trtandcoef' node_new `nodecoef' tumor_new `tumorcoef')

    /*--- Administrative censoring ---*/

    gen id      = _n
    gen admin   = 6
    gen admin_w = admin * 52
    gen timeDFS_w = timeDFS_y * 52
    replace timeDFS_w = round(timeDFS_w)
    replace timeDFS_w = 1 if timeDFS_w == 0
    replace admin_w   = round(admin_w)
    gen dead = 1 if timeDFS_w <  admin_w
    replace dead = 0 if timeDFS_w >= admin_w
    replace timeDFS_w = admin_w if timeDFS_w > admin_w

    /*--- Simulate non-compliance ---*/

    gen p     = invlogit($NC_INTERCEPT + 0.6 * age + 0.7 * toxicity) if trtrand == 1
    gen noncp = rbinomial(1, p)                            if trtrand == 1
    replace noncp = 0 if noncp == .
    drop p

    /*--- Time of first non-compliance (in weeks) ---*/

    gen noncp_time1 = rbinomial(1, 0.24) if noncp == 1
    gen noncp_time2 = rbinomial(1, 0.34) if noncp == 1 & noncp_time1 == 0
    gen noncp_time3 = rbinomial(1, 0.6)  if noncp == 1 & noncp_time1 == 0 & noncp_time2 == 0
    gen noncp_time  = 1  if noncp_time1 == 1
    replace noncp_time = 4  if noncp_time2 == 1
    replace noncp_time = 11 if noncp_time3 == 1 & noncp == 1
    replace noncp_time = 7  if noncp_time == . & noncp == 1
    replace noncp_time = .  if noncp_time > timeDFS_w & trt == 1 & noncp == 1
    replace noncp      = 0  if noncp_time == . & trt == 1

    /*--- 2. New survival times for non-compliers ---*/

    gen xoOSgainobs = timeDFS_w - noncp_time if noncp == 1
    replace xoOSgainobs = round(xoOSgainobs * `aftcoef')
    gen timeDFS_w2 = cond(noncp == 1, xoOSgainobs + noncp_time, timeDFS_w)

    gen died = 0
    replace died = 1 if timeDFS_w2 < admin_w & dead == 1
    replace timeDFS_w2 = admin_w if timeDFS_w2 > admin_w
    replace xoOSgainobs = . if noncp_time >  admin_w
    replace noncp       = . if noncp_time >= admin_w
    replace noncp_time  = . if noncp_time >= admin_w

    /*--- 3. True dRMST, rRMST and MSD ---*/

    scalar lambda = `lambda'
    scalar gamma  = `gammas'
    scalar beta1  = `trtandcoef'

    gen timeDFS_y1 = timeDFS_y
    replace timeDFS_y1 = 6 if timeDFS_y1 > 6
    stset timeDFS_y1, failure(dead = 1)
    sort _t
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    gen H0      = lambda * _t^(gamma)
    gen H1      = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen au_S0_gp = lag * S0_true
    gen au_S1_gp = lag * S1_true

    foreach i in 3 {
        egen le_S0_gp`i' = sum(au_S0_gp) if _t <= `i'
        egen le_S1_gp`i' = sum(au_S1_gp) if _t <= `i'
        summarize le_S0_gp`i'
        return scalar le0_gp`i' = r(mean)
        summarize le_S1_gp`i'
        return scalar le1_gp`i' = r(mean)
    }

    drop H0 H1 S1_true S0_true
    gen H0 = lambda * 3^(gamma)
    gen H1 = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen S_diff_true = S1_true - S0_true
    sum S_diff_true
    return scalar S_diff_true = r(mean)

    /*=== IPCW (logistic) ===*/

    /*--- Build the time-varying panel (15 scheduled visits per patient,
          trimmed to each patient's follow-up time) ---*/

    expand 15
    sort id
    by id: gen visit = _n
    gen cumweek = 1   if visit == 1
    replace cumweek = 4   if visit == 2
    replace cumweek = 7   if visit == 3
    replace cumweek = 11  if visit == 4
    replace cumweek = 15  if visit == 5
    replace cumweek = 19  if visit == 6
    replace cumweek = 35  if visit == 7
    replace cumweek = 51  if visit == 8
    replace cumweek = 67  if visit == 9
    replace cumweek = 91  if visit == 10
    replace cumweek = 115 if visit == 11
    replace cumweek = 139 if visit == 12
    replace cumweek = 163 if visit == 13
    replace cumweek = 215 if visit == 14
    replace cumweek = 267 if visit == 15

    drop if cumweek > timeDFS_w2
    sort id
    by id: gen firstobs = 0
    by id: replace firstobs = 1 if _n == 1
    by id: gen finalobs = 0
    by id: replace finalobs = 1 if _n == _N
    by id: drop if cumweek > noncp_time
    by id: replace finalobs = 0
    by id: replace finalobs = 1 if _n == _N

    gen censOS = 0
    replace censOS = 1 if finalobs == 1 & died == 0 & visit == 15 & timeDFS_w2 == 312
    sum censOS
    replace died = . if censOS == 1

    /*--- Time-dependent compliance and death indicators ---*/

    gen compliance_tdo = 0
    replace compliance_tdo = 1 if trtrand == 1 & noncp_time == cumweek & noncp == 1
    by id: replace died = . if (compliance_tdo == 1 | censOS == 1)

    gen death_tdo = .
    sort id
    by id: replace death_tdo = 0 if finalobs == 0
    by id: replace death_tdo = 1 if died == 1 & finalobs == 1
    by id: replace death_tdo = . if compliance_tdo == 1 & finalobs == 1 & died == 1

    /*--- Entry cumulative time ---*/

    gen cumweek_entry = .
    sort id
    by id: replace cumweek_entry = 0              if _n == 1
    sort id
    by id: replace cumweek_entry = cumweek[_n - 1] if _n != 1

    /*--- Weights ---*/

    * Denominator: probability of remaining uncensored given baseline + time-varying covariates
    xi: logistic compliance_tdo age i.toxicity if trtrand == 1
    predict p1 if e(sample)
    gen pa = p1 * compliance_tdo + (1 - p1) * (1 - compliance_tdo)
    replace pa = 1 if pa == .
    sort id cumweek
    by id: replace pa = pa * pa[_n - 1] if _n != 1
    rename pa p_denom

    * Numerator: probability given baseline covariates only
    xi: logistic compliance_tdo age if trtrand == 1
    predict p2 if e(sample)
    gen pa2 = p2 * compliance_tdo + (1 - p2) * (1 - compliance_tdo)
    replace pa2 = 1 if pa2 == .
    sort id cumweek
    by id: replace pa2 = pa2 * pa2[_n - 1] if _n != 1
    rename pa2 p_num

    gen weight  = 1     / p_denom
    gen sweight = p_num / p_denom

    * Set weight = 1 for RT (control) arm
    replace sweight = 1 if trt == 0
    replace weight  = 1 if trt == 0

    /*--- Stabilised weight distribution ---*/

    summarize sweight if noncp == 1, detail
    return scalar sw_mean1   = r(mean)
    return scalar sw_median1 = r(p50)
    return scalar sw_se1     = r(sd)
    return scalar sw_min1    = r(min)
    return scalar sw_max1    = r(max)

    summarize sweight if noncp != 1, detail
    return scalar sw_mean   = r(mean)
    return scalar sw_median = r(p50)
    return scalar sw_se     = r(sd)
    return scalar sw_min    = r(min)
    return scalar sw_max    = r(max)

    summarize sweight if trt == 1, detail
    return scalar sw_mean2   = r(mean)
    return scalar sw_median2 = r(p50)
    return scalar sw_se2     = r(sd)
    return scalar sw_min2    = r(min)
    return scalar sw_max2    = r(max)

    summarize sweight if trt == 0, detail
    return scalar sw_mean3   = r(mean)
    return scalar sw_median3 = r(p50)
    return scalar sw_se3     = r(sd)
    return scalar sw_min3    = r(min)
    return scalar sw_max3    = r(max)

    /*--- Weighted Cox model ---*/

    stset cumweek [pweight = sweight], failure(death_tdo == 1) ///
        entry(cumweek_entry) exit(cumweek)
    stcox trtrand node_new tumor_new, vce(robust)

    return scalar ipcw_cox_hr    = exp(_b[trtrand])
    return scalar ipcw_cox_hr_SE = exp(_b[trtrand]) * _se[trtrand]
    return scalar ipcw_cox_hr_LB = exp(_b[trtrand] - 1.96 * _se[trtrand])
    return scalar ipcw_cox_hr_UB = exp(_b[trtrand] + 1.96 * _se[trtrand])

    gen cov_ipcw_cox_hr = 0
    replace cov_ipcw_cox_hr = 1 if exp(`trtandcoef') > return(ipcw_cox_hr_LB) & exp(`trtandcoef') < return(ipcw_cox_hr_UB)
    return scalar cov_ipcw_cox_hr = cov_ipcw_cox_hr[1]

    gen power_ipcw_cox_hr = 0
    replace power_ipcw_cox_hr = 1 if 1 > return(ipcw_cox_hr_LB) & 1 < return(ipcw_cox_hr_UB)
    return scalar power_ipcw_cox_hr = power_ipcw_cox_hr[1]

    /*--- RMST and dSurvival from Weibull fit + IPCW HR ---*/

    collapse (max) trtrand timeDFS_w2 node_new tumor_new new_v xoOSgainobs ///
        noncp_time died, by(id)
    by id: replace noncp_time = 0 if noncp_time == .

    preserve
    stset timeDFS_w2, failure(died) id(id)
    xi: streg node_new tumor_new if trtrand == 0, distribution(weibull) nohr

    gen cons = _b[_cons]
    gen lnp  = _b[/ln_p]

    expand 2
    sort id
    gen time_w = 0
    replace time_w = _n * (312 / _N) if _n > 1
    gen time_y = time_w / 52

    sort time_w
    gen survf = 1
    replace survf = exp(-exp(cons) * time_w^exp(lnp)) if _n > 1
    gen hazf = 0
    replace hazf = 1 - survf / survf[_n - 1] if _n > 1
    gen hazfexp = hazf * return(ipcw_cox_hr)
    gen survfexp = 1
    replace survfexp = survfexp[_n - 1] * (1 - hazfexp) if _n > 1

    gen S0_ipcw = survf
    gen S1_ipcw = survfexp

    gen lag = time_y - time_y[_n - 1] if _n > 1
    replace lag = time_y if _n == 1
    gen au_S0_ipcw = lag * survf
    gen au_S1_ipcw = lag * survfexp

    foreach i in 3 {
        egen le_S0_ipcw`i' = sum(au_S0_ipcw) if _t <= `i'
        egen le_S1_ipcw`i' = sum(au_S1_ipcw) if _t <= `i'
        summarize le_S0_ipcw`i'
        return scalar le0_ipcw`i' = r(mean)
        summarize le_S1_ipcw`i'
        return scalar le1_ipcw`i' = r(mean)

        sort _t
        gen t1 = `i' - _t
        replace t1 = 999999 if _t > `i'
        egen min_t1 = min(t1)
        sum S0_ipcw if t1 == min_t1
        return scalar S0_ipcw1 = r(mean)
        sum S1_ipcw if t1 == min_t1
        return scalar S1_ipcw1 = r(mean)
    }

end

simulate                                                              ///
    sw_mean1   = r(sw_mean1)   sw_median1 = r(sw_median1)              ///
    sw_se1     = r(sw_se1)     sw_min1    = r(sw_min1)   sw_max1 = r(sw_max1) ///
    sw_mean    = r(sw_mean)    sw_median  = r(sw_median)                ///
    sw_se      = r(sw_se)      sw_min     = r(sw_min)    sw_max  = r(sw_max)  ///
    sw_mean2   = r(sw_mean2)   sw_median2 = r(sw_median2)               ///
    sw_se2     = r(sw_se2)     sw_min2    = r(sw_min2)   sw_max2 = r(sw_max2) ///
    sw_mean3   = r(sw_mean3)   sw_median3 = r(sw_median3)               ///
    sw_se3     = r(sw_se3)     sw_min3    = r(sw_min3)   sw_max3 = r(sw_max3) ///
    le0_gp3    = r(le0_gp3)    le1_gp3    = r(le1_gp3)                  ///
    ipcw_cox_hr= r(ipcw_cox_hr)                                          ///
    cov_ipcw_cox_hr   = r(cov_ipcw_cox_hr)                               ///
    power_ipcw_cox_hr = r(power_ipcw_cox_hr)                             ///
    le0_ipcw3  = r(le0_ipcw3)  le1_ipcw3  = r(le1_ipcw3)                ///
    S0_ipcw1   = r(S0_ipcw1)   S1_ipcw1   = r(S1_ipcw1),                ///
    reps($N_REPS) saving(simstudy_ipcw, replace) seed($SEED): ///
    simstudy_ipcw


/*----------------------------------------------------------------------------
  Program 3 of 6: IPCW with Cox weights, no selection bias
----------------------------------------------------------------------------*/

capture program drop simstudyipcw_cox
program define simstudyipcw_cox, rclass

    syntax [, obs(int $N_OBS)             ///
              trtandcoef(real $TRT_COEF) ///
              nodecoef(real $NODE_COEF)           ///
              tumorcoef(real $TUMOR_COEF)          ///
              aftcoef(real $AFT_COEF)           ///
              lambda(real $LAMBDA_NS)       ///
              gammas(real $GAMMA_S)             ///
            ]

    clear
    pr drop _all

    /*--- 1. Generate baseline covariates and survival times ---*/

    rbinary age toxicity new_v tumor_new node_new,                         ///
        means(0.33, 0.48, 0.45, 0.58, 0.28)                                ///
        corr(1    , 0.4 , 0.68, 0   , 0    \                               ///
             0.4  , 1   , 0.68, 0   , 0    \                               ///
             0.68 , 0.68, 1   , 0   , 0    \                               ///
             0    , 0   , 0   , 1   , 0.5  \                               ///
             0    , 0   , 0   , 0.5 , 1)                                   ///
        n($N_OBS)

    set obs `obs'
    gen trtrand = rbinomial(1, 0.5)

    survsim timeDFS_y, distribution(weibull) lambda(`lambda') gammas(`gammas') ///
        cov(trtrand `trtandcoef' node_new `nodecoef' tumor_new `tumorcoef')

    gen id      = _n
    gen admin   = 6
    gen admin_w = admin * 52
    gen timeDFS_w = timeDFS_y * 52
    replace timeDFS_w = round(timeDFS_w)
    replace timeDFS_w = 1 if timeDFS_w == 0
    replace admin_w   = round(admin_w)
    gen dead = 1 if timeDFS_w <  admin_w
    replace dead = 0 if timeDFS_w >= admin_w
    replace timeDFS_w = admin_w if timeDFS_w > admin_w

    /*--- Non-compliance ---*/

    gen p     = invlogit($NC_INTERCEPT + 0.6 * age + 0.7 * toxicity) if trtrand == 1
    gen noncp = rbinomial(1, p)                            if trtrand == 1
    replace noncp = 0 if noncp == .
    drop p

    /*--- Time of first non-compliance (in weeks) ---*/

    gen noncp_time1 = rbinomial(1, 0.24) if noncp == 1
    gen noncp_time2 = rbinomial(1, 0.34) if noncp == 1 & noncp_time1 == 0
    gen noncp_time3 = rbinomial(1, 0.6)  if noncp == 1 & noncp_time1 == 0 & noncp_time2 == 0
    gen noncp_time  = 1  if noncp_time1 == 1
    replace noncp_time = 4  if noncp_time2 == 1
    replace noncp_time = 11 if noncp_time3 == 1 & noncp == 1
    replace noncp_time = 7  if noncp_time == . & noncp == 1
    replace noncp_time = .  if noncp_time > timeDFS_w & trt == 1 & noncp == 1
    replace noncp      = 0  if noncp_time == . & trt == 1

    /*--- 2. New survival times for non-compliers ---*/

    gen xoOSgainobs = timeDFS_w - noncp_time if noncp == 1
    replace xoOSgainobs = round(xoOSgainobs * `aftcoef')
    gen timeDFS_w2 = cond(noncp == 1, xoOSgainobs + noncp_time, timeDFS_w)

    gen died = 0
    replace died = 1 if timeDFS_w2 < admin_w & dead == 1
    replace timeDFS_w2 = admin_w if timeDFS_w2 > admin_w
    replace xoOSgainobs = . if noncp_time >  admin_w
    replace noncp       = . if noncp_time >= admin_w
    replace noncp_time  = . if noncp_time >= admin_w

    /*--- 3. True dRMST, rRMST and MSD ---*/

    scalar lambda = `lambda'
    scalar gamma  = `gammas'
    scalar beta1  = `trtandcoef'

    gen timeDFS_y1 = timeDFS_y
    replace timeDFS_y1 = 6 if timeDFS_y1 > 6
    stset timeDFS_y1, failure(dead = 1)
    sort _t
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    gen H0      = lambda * _t^(gamma)
    gen H1      = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen au_S0_gp = lag * S0_true
    gen au_S1_gp = lag * S1_true

    foreach i in 3 {
        egen le_S0_gp`i' = sum(au_S0_gp) if _t <= `i'
        egen le_S1_gp`i' = sum(au_S1_gp) if _t <= `i'
        summarize le_S0_gp`i'
        return scalar le0_gp`i' = r(mean)
        summarize le_S1_gp`i'
        return scalar le1_gp`i' = r(mean)
    }

    drop H0 H1 S1_true S0_true
    gen H0 = lambda * 3^(gamma)
    gen H1 = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen S_diff_true = S1_true - S0_true
    sum S_diff_true
    return scalar S_diff_true = r(mean)

    /*=== IPCW (Cox-based weights) ===*/

    /*--- Build the time-varying panel (15 scheduled visits per patient,
          trimmed to each patient's follow-up time) ---*/

    expand 15
    sort id
    by id: gen visit = _n
    gen cumweek = 1   if visit == 1
    replace cumweek = 4   if visit == 2
    replace cumweek = 7   if visit == 3
    replace cumweek = 11  if visit == 4
    replace cumweek = 15  if visit == 5
    replace cumweek = 19  if visit == 6
    replace cumweek = 35  if visit == 7
    replace cumweek = 51  if visit == 8
    replace cumweek = 67  if visit == 9
    replace cumweek = 91  if visit == 10
    replace cumweek = 115 if visit == 11
    replace cumweek = 139 if visit == 12
    replace cumweek = 163 if visit == 13
    replace cumweek = 215 if visit == 14
    replace cumweek = 267 if visit == 15

    drop if cumweek > timeDFS_w2
    sort id
    by id: gen firstobs = 0
    by id: replace firstobs = 1 if _n == 1
    by id: gen finalobs = 0
    by id: replace finalobs = 1 if _n == _N
    by id: drop if cumweek > noncp_time
    by id: replace finalobs = 0
    by id: replace finalobs = 1 if _n == _N

    gen censOS = 0
    replace censOS = 1 if finalobs == 1 & died == 0 & visit == 15 & timeDFS_w2 == 312
    sum censOS
    replace died = . if censOS == 1

    gen compliance_tdo = 0
    replace compliance_tdo = 1 if trtrand == 1 & noncp_time == cumweek & noncp == 1
    by id: replace died = . if (compliance_tdo == 1 | censOS == 1)

    gen death_tdo = .
    sort id
    by id: replace death_tdo = 0 if finalobs == 0
    by id: replace death_tdo = 1 if died == 1 & finalobs == 1
    by id: replace death_tdo = . if compliance_tdo == 1 & finalobs == 1 & died == 1

    gen cumweek_entry = .
    sort id
    by id: replace cumweek_entry = 0              if _n == 1
    sort id
    by id: replace cumweek_entry = cumweek[_n - 1] if _n != 1

    /*--- Denominator (Cz): Cox model with baseline + time-varying covariates ---*/

    stset cumweek if trt == 1, failure(compliance_tdo == 1) enter(time cumweek_entry)
    stcox age i.toxicity if trt == 1
    predict xb_cz, xb
    gen exp_xb_cz = exp(xb_cz)
    predict basech_cz, basechazard
    by id: gen surv_cz = exp(-basech_cz[_n] * exp_xb_cz[_n])
    stset cumweek if trt == 1, failure(compliance_tdo == 1) enter(time cumweek_entry)

    /*--- Numerator (C0): Cox model with baseline covariates only ---*/

    stcox age if trt == 1
    predict xb_c0, xb
    gen exp_xb_c0 = exp(xb_c0)
    predict basech_c0, basechazard
    by id: gen surv_c0 = exp(-basech_c0[_n] * exp_xb_c0[_n])

    gen sweight = surv_c0 / surv_cz
    replace sweight = 1 if trt == 0

    /*--- Stabilised weight distribution ---*/

    summarize sweight if noncp == 1, detail
    return scalar sw_mean1   = r(mean)
    return scalar sw_median1 = r(p50)
    return scalar sw_se1     = r(sd)
    return scalar sw_min1    = r(min)
    return scalar sw_max1    = r(max)

    summarize sweight if noncp != 1, detail
    return scalar sw_mean   = r(mean)
    return scalar sw_median = r(p50)
    return scalar sw_se     = r(sd)
    return scalar sw_min    = r(min)
    return scalar sw_max    = r(max)

    summarize sweight if trt == 1, detail
    return scalar sw_mean2   = r(mean)
    return scalar sw_median2 = r(p50)
    return scalar sw_se2     = r(sd)
    return scalar sw_min2    = r(min)
    return scalar sw_max2    = r(max)

    summarize sweight if trt == 0, detail
    return scalar sw_mean3   = r(mean)
    return scalar sw_median3 = r(p50)
    return scalar sw_se3     = r(sd)
    return scalar sw_min3    = r(min)
    return scalar sw_max3    = r(max)

    /*--- Weighted Cox model ---*/

    stset cumweek [pweight = sweight], failure(death_tdo == 1) ///
        entry(cumweek_entry) exit(cumweek)
    stcox trtrand node_new tumor_new, vce(robust)

    return scalar ipcw_cox_hr    = exp(_b[trtrand])
    return scalar ipcw_cox_hr_SE = exp(_b[trtrand]) * _se[trtrand]
    return scalar ipcw_cox_hr_LB = exp(_b[trtrand] - 1.96 * _se[trtrand])
    return scalar ipcw_cox_hr_UB = exp(_b[trtrand] + 1.96 * _se[trtrand])

    gen cov_ipcw_cox_hr = 0
    replace cov_ipcw_cox_hr = 1 if exp(`trtandcoef') > return(ipcw_cox_hr_LB) & exp(`trtandcoef') < return(ipcw_cox_hr_UB)
    return scalar cov_ipcw_cox_hr = cov_ipcw_cox_hr[1]

    gen power_ipcw_cox_hr = 0
    replace power_ipcw_cox_hr = 1 if 1 > return(ipcw_cox_hr_LB) & 1 < return(ipcw_cox_hr_UB)
    return scalar power_ipcw_cox_hr = power_ipcw_cox_hr[1]

    /*--- RMST and dSurvival ---*/

    collapse (max) trtrand timeDFS_w2 node_new tumor_new new_v xoOSgainobs ///
        noncp_time died, by(id)
    by id: replace noncp_time = 0 if noncp_time == .

    preserve
    stset timeDFS_w2, failure(died) id(id)
    xi: streg node_new tumor_new if trtrand == 0, distribution(weibull) nohr

    gen cons = _b[_cons]
    gen lnp  = _b[/ln_p]

    expand 2
    sort id
    gen time_w = 0
    replace time_w = _n * (312 / _N) if _n > 1
    gen time_y = time_w / 52

    sort time_w
    gen survf = 1
    replace survf = exp(-exp(cons) * time_w^exp(lnp)) if _n > 1
    gen hazf = 0
    replace hazf = 1 - survf / survf[_n - 1] if _n > 1
    gen hazfexp = hazf * return(ipcw_cox_hr)
    gen survfexp = 1
    replace survfexp = survfexp[_n - 1] * (1 - hazfexp) if _n > 1

    gen S0_ipcw = survf
    gen S1_ipcw = survfexp

    gen lag = time_y - time_y[_n - 1] if _n > 1
    replace lag = time_y if _n == 1
    gen au_S0_ipcw = lag * survf
    gen au_S1_ipcw = lag * survfexp

    foreach i in 3 {
        egen le_S0_ipcw`i' = sum(au_S0_ipcw) if _t <= `i'
        egen le_S1_ipcw`i' = sum(au_S1_ipcw) if _t <= `i'
        summarize le_S0_ipcw`i'
        return scalar le0_ipcw`i' = r(mean)
        summarize le_S1_ipcw`i'
        return scalar le1_ipcw`i' = r(mean)

        sort _t
        gen t1 = `i' - _t
        replace t1 = 999999 if _t > `i'
        egen min_t1 = min(t1)
        sum S0_ipcw if t1 == min_t1
        return scalar S0_ipcw1 = r(mean)
        sum S1_ipcw if t1 == min_t1
        return scalar S1_ipcw1 = r(mean)
    }

end

simulate                                                              ///
    sw_mean1   = r(sw_mean1)   sw_median1 = r(sw_median1)              ///
    sw_se1     = r(sw_se1)     sw_min1    = r(sw_min1)   sw_max1 = r(sw_max1) ///
    sw_mean    = r(sw_mean)    sw_median  = r(sw_median)                ///
    sw_se      = r(sw_se)      sw_min     = r(sw_min)    sw_max  = r(sw_max)  ///
    sw_mean2   = r(sw_mean2)   sw_median2 = r(sw_median2)               ///
    sw_se2     = r(sw_se2)     sw_min2    = r(sw_min2)   sw_max2 = r(sw_max2) ///
    sw_mean3   = r(sw_mean3)   sw_median3 = r(sw_median3)               ///
    sw_se3     = r(sw_se3)     sw_min3    = r(sw_min3)   sw_max3 = r(sw_max3) ///
    le0_gp3    = r(le0_gp3)    le1_gp3    = r(le1_gp3)                  ///
    ipcw_cox_hr= r(ipcw_cox_hr)                                          ///
    cov_ipcw_cox_hr   = r(cov_ipcw_cox_hr)                               ///
    power_ipcw_cox_hr = r(power_ipcw_cox_hr)                             ///
    le0_ipcw3  = r(le0_ipcw3)  le1_ipcw3  = r(le1_ipcw3)                ///
    S0_ipcw1   = r(S0_ipcw1)   S1_ipcw1   = r(S1_ipcw1),                ///
    reps($N_REPS) saving(simstudyipcw_cox, replace) seed($SEED): ///
    simstudyipcw_cox


/*============================================================================
                     PART B. SELECTION BIAS SCENARIO
============================================================================*/

/*----------------------------------------------------------------------------
  Program 4 of 6: RPSFTM (Cox and Weibull), selection bias
----------------------------------------------------------------------------*/

capture program drop simstudy_rpsftm
program define simstudy_rpsftm, rclass

    syntax [, obs(int $N_OBS)             ///
              trtandcoef(real $TRT_COEF) ///
              nodecoef(real $NODE_COEF)           ///
              tumorcoef(real $TUMOR_COEF)          ///
              aftcoef(real $AFT_COEF)           ///
              newvcoeff(real $NEWV_COEF)    ///
              lambda(real $LAMBDA_SEL)      ///
              gammas(real $GAMMA_S)             ///
            ]

    clear
    pr drop _all

    /*--- 1. Generate baseline covariates and survival times ---*/

    rbinary age toxicity new_v tumor_new node_new,                         ///
        means(0.33, 0.48, 0.45, 0.58, 0.28)                                ///
        corr(1    , 0.4 , 0.68, 0   , 0    \                               ///
             0.4  , 1   , 0.68, 0   , 0    \                               ///
             0.68 , 0.68, 1   , 0   , 0    \                               ///
             0    , 0   , 0   , 1   , 0.5  \                               ///
             0    , 0   , 0   , 0.5 , 1)                                   ///
        n($N_OBS)

    set obs `obs'
    gen trtrand = rbinomial(1, 0.5)

    * Selection bias is induced by including new_v in the survival model
    survsim timeDFS_y, distribution(weibull) lambda(`lambda') gammas(`gammas') ///
        cov(trtrand `trtandcoef' node_new `nodecoef' tumor_new `tumorcoef'    ///
            new_v `newvcoeff')

    /*--- Administrative censoring ---*/

    gen id      = _n
    gen admin   = 6
    gen admin_w = admin * 52
    gen timeDFS_w = timeDFS_y * 52
    replace timeDFS_w = round(timeDFS_w)
    replace timeDFS_w = 1 if timeDFS_w == 0
    replace admin_w   = round(admin_w)
    gen dead = 1 if timeDFS_w <  admin_w
    replace dead = 0 if timeDFS_w >= admin_w
    replace timeDFS_w = admin_w if timeDFS_w > admin_w

    /*--- Non-compliance ---*/

    gen p     = invlogit($NC_INTERCEPT + 0.6 * age + 0.7 * toxicity) if trtrand == 1
    gen noncp = rbinomial(1, p)                            if trtrand == 1
    replace noncp = 0 if noncp == .
    drop p

    /*--- Time of first non-compliance (in weeks) ---*/

    gen noncp_time1 = rbinomial(1, 0.24) if noncp == 1
    gen noncp_time2 = rbinomial(1, 0.34) if noncp == 1 & noncp_time1 == 0
    gen noncp_time3 = rbinomial(1, 0.6)  if noncp == 1 & noncp_time1 == 0 & noncp_time2 == 0
    gen noncp_time  = 1  if noncp_time1 == 1
    replace noncp_time = 4  if noncp_time2 == 1
    replace noncp_time = 11 if noncp_time3 == 1 & noncp == 1
    replace noncp_time = 7  if noncp_time == . & noncp == 1
    replace noncp_time = .  if noncp_time > timeDFS_w & trt == 1 & noncp == 1
    replace noncp      = 0  if noncp_time == . & trt == 1

    /*--- 2. New survival times for non-compliers ---*/

    gen xoOSgainobs = timeDFS_w - noncp_time if noncp == 1
    replace xoOSgainobs = round(xoOSgainobs * `aftcoef')
    gen timeDFS_w2 = cond(noncp == 1, xoOSgainobs + noncp_time, timeDFS_w)

    gen died = 0
    replace died = 1 if timeDFS_w2 < admin_w & dead == 1
    replace timeDFS_w2 = admin_w if timeDFS_w2 > admin_w
    replace xoOSgainobs = . if noncp_time >  admin_w
    replace noncp       = . if noncp_time >= admin_w
    replace noncp_time  = . if noncp_time >= admin_w

    /*--- 3. True dRMST, rRMST and MSD ---*/

    scalar lambda = `lambda'
    scalar gamma  = `gammas'
    scalar beta1  = `trtandcoef'

    gen timeDFS_y1 = timeDFS_y
    replace timeDFS_y1 = 6 if timeDFS_y1 > 6
    stset timeDFS_y1, failure(dead = 1)
    sort _t
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    gen H0      = lambda * _t^(gamma)
    gen H1      = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen au_S0_gp = lag * S0_true
    gen au_S1_gp = lag * S1_true

    foreach i in 3 {
        egen le_S0_gp`i' = sum(au_S0_gp) if _t <= `i'
        egen le_S1_gp`i' = sum(au_S1_gp) if _t <= `i'
        summarize le_S0_gp`i'
        return scalar le0_gp`i' = r(mean)
        summarize le_S1_gp`i'
        return scalar le1_gp`i' = r(mean)
    }

    drop H0 H1 S1_true S0_true
    gen H0 = lambda * 3^(gamma)
    gen H1 = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen S_diff_true = S1_true - S0_true
    sum S_diff_true
    return scalar S_diff_true = r(mean)

    /*=== RPSFTM (Cox-based test) ===*/

    replace noncp      = 0          if noncp == .
    replace noncp_time = timeDFS_w2 if noncp == 0
    replace noncp_time = timeDFS_w2 if trt   == 0

    stset timeDFS_w2, failure(died == 1)
    strbee trtrand, xo1(noncp_time noncp) adj(node_new tumor_new) ///
        endstudy(admin_w) gen(c1)
    strbee, hr

    return scalar rpsftm_hr        = r(HR_adj)
    return scalar rpsftm_hr_LB     = r(HR_adj_low)
    return scalar rpsftm_hr_UB     = r(HR_adj_upp)
    return scalar rpsftm_aft       = exp(r(psi))
    return scalar rpsftm_aft_LB    = exp(r(psi_low))
    return scalar rpsftm_aft_UB    = exp(r(psi_upp))
    return scalar itt_statistics   = r(Z_ITT)

    gen cov_rpsftm_hr = 0
    replace cov_rpsftm_hr = 1 if exp(`trtandcoef') > return(rpsftm_hr_LB) & exp(`trtandcoef') < return(rpsftm_hr_UB)
    return scalar cov_rpsftm_hr = cov_rpsftm_hr[1]

    gen cov_rpsftm_aft = 0
    replace cov_rpsftm_aft = 1 if `aftcoef' > return(rpsftm_aft_LB) & `aftcoef' < return(rpsftm_aft_UB)
    return scalar cov_rpsftm_aft = cov_rpsftm_aft[1]

    gen power_rpsftm_hr = 0
    replace power_rpsftm_hr = 1 if 1 > return(rpsftm_hr_LB) & 1 < return(rpsftm_hr_UB)
    return scalar power_rpsftm_hr = power_rpsftm_hr[1]

    gen power_rpsftm_aft = 0
    replace power_rpsftm_aft = 1 if 1 > return(rpsftm_aft_LB) & 1 < return(rpsftm_aft_UB)
    return scalar power_rpsftm_aft = power_rpsftm_aft[1]

    /*--- RMST and dSurvival, WITH recensoring ---*/

    gen newtime = .
    replace newtime = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft) if noncp == 1
    replace newtime = timeDFS_w2 if newtime == .
    replace newtime = admin_w    if trt == 1 & newtime > admin_w

    gen dd = dc1  if trt == 1
    replace dd = died if trt == 0
    gen newtime1 = newtime / 52
    stset newtime1, failure(dd == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr

    foreach i in 3 {
        predict S0_rpsftm1, basesurv
        gen S1_rpsftm1    = S0_rpsftm1^return(rpsftm_hr)
        gen au_S0_rpsftm1 = lag * S0_rpsftm1
        gen au_S1_rpsftm1 = lag * S1_rpsftm1
        egen le_S0_rpsftm1`i' = sum(au_S0_rpsftm1) if _t <= `i'
        egen le_S1_rpsftm1`i' = sum(au_S1_rpsftm1) if _t <= `i'
        summarize le_S0_rpsftm1`i'
        return scalar le0_rpsftm1`i' = r(mean)
        summarize le_S1_rpsftm1`i'
        return scalar le1_rpsftm1`i' = r(mean)

        sort _t
        gen t1 = `i' - _t
        replace t1 = 999999 if _t > `i'
        egen min_t1 = min(t1)
        sum S0_rpsftm1 if t1 == min_t1
        return scalar S0_rpsftm1 = r(mean)
        sum S1_rpsftm1 if t1 == min_t1
        return scalar S1_rpsftm1 = r(mean)
    }

    /*--- RMST and dSurvival, WITHOUT recensoring ---*/

    gen newtime_wr = .
    replace newtime_wr = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft) if noncp == 1
    replace newtime_wr = timeDFS_w2 if newtime_wr == .
    gen newtime1_wr = newtime_wr / 52
    stset newtime1_wr, failure(died == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr
    return scalar rpsftm_hr_nore    = exp(_b[trtrand])
    return scalar rpsftm_hr_nore_SE = _b[trtrand] / return(itt_statistics)
    return scalar rpsftm_hr_nore_LB = exp(_b[trtrand] - 1.96 * return(rpsftm_hr_nore_SE))
    return scalar rpsftm_hr_nore_UB = exp(_b[trtrand] + 1.96 * return(rpsftm_hr_nore_SE))

    gen cov_rpsftm_hr_nore = 0
    replace cov_rpsftm_hr_nore = 1 if exp(`trtandcoef') > return(rpsftm_hr_nore_LB) & exp(`trtandcoef') < return(rpsftm_hr_nore_UB)
    return scalar cov_rpsftm_hr_nore = cov_rpsftm_hr_nore[1]

    gen power_rpsftm_hr_nore = 0
    replace power_rpsftm_hr_nore = 1 if 1 > return(rpsftm_hr_nore_LB) & 1 < return(rpsftm_hr_nore_UB)
    return scalar power_rpsftm_hr_nore = power_rpsftm_hr_nore[1]

    stcox trt node_new tumor_new, nohr
    foreach i in 3 {
        predict S0_rpsftm1_wr, basesurv
        gen S1_rpsftm1_wr    = S0_rpsftm1_wr^return(rpsftm_hr_nore)
        gen au_S0_rpsftm1_wr = lag * S0_rpsftm1_wr
        gen au_S1_rpsftm1_wr = lag * S1_rpsftm1_wr
        egen le_S0_rpsftm1`i'_wr = sum(au_S0_rpsftm1_wr) if _t <= `i'
        egen le_S1_rpsftm1`i'_wr = sum(au_S1_rpsftm1_wr) if _t <= `i'
        summarize le_S0_rpsftm1`i'_wr
        return scalar le0_rpsftm1`i'_wr = r(mean)
        summarize le_S1_rpsftm1`i'_wr
        return scalar le1_rpsftm1`i'_wr = r(mean)

        sort _t
        gen t1_wr = `i' - _t
        replace t1_wr = 999999 if _t > `i'
        egen min_t1_wr = min(t1_wr)
        sum S0_rpsftm1_wr if t1_wr == min_t1_wr
        return scalar S0_rpsftm1_wr = r(mean)
        sum S1_rpsftm1_wr if t1_wr == min_t1_wr
        return scalar S1_rpsftm1_wr = r(mean)
    }

    /*=== RPSFTM (Weibull-based test) ===*/

    drop cc1 dd dc1 c1
    stset timeDFS_w2, failure(died == 1)
    strbee trtrand, xo1(noncp_time noncp) adj(node_new tumor_new) ///
        endstudy(admin_w) gen(c1) test(weibull)
    strbee, hr

    return scalar rpsftm_hr_w        = r(HR_adj)
    return scalar rpsftm_hr_LB_w     = r(HR_adj_low)
    return scalar rpsftm_hr_UB_w     = r(HR_adj_upp)
    return scalar rpsftm_aft_w       = exp(r(psi))
    return scalar rpsftm_aft_LB_w    = exp(r(psi_low))
    return scalar rpsftm_aft_UB_w    = exp(r(psi_upp))
    return scalar itt_statistics_w   = r(Z_ITT)

    gen cov_rpsftm_hr_w = 0
    replace cov_rpsftm_hr_w = 1 if exp(`trtandcoef') > return(rpsftm_hr_LB_w) & exp(`trtandcoef') < return(rpsftm_hr_UB_w)
    return scalar cov_rpsftm_hr_w = cov_rpsftm_hr_w[1]

    gen cov_rpsftm_aft_w = 0
    replace cov_rpsftm_aft_w = 1 if `aftcoef' > return(rpsftm_aft_LB_w) & `aftcoef' < return(rpsftm_aft_UB_w)
    return scalar cov_rpsftm_aft_w = cov_rpsftm_aft_w[1]

    gen power_rpsftm_hr_w = 0
    replace power_rpsftm_hr_w = 1 if 1 > return(rpsftm_hr_LB_w) & 1 < return(rpsftm_hr_UB_w)
    return scalar power_rpsftm_hr_w = power_rpsftm_hr_w[1]

    gen power_rpsftm_aft_w = 0
    replace power_rpsftm_aft_w = 1 if 1 > return(rpsftm_aft_LB_w) & 1 < return(rpsftm_aft_UB_w)
    return scalar power_rpsftm_aft_w = power_rpsftm_aft_w[1]

    /*--- RMST and dSurvival, Weibull, WITH recensoring ---*/

    gen newtime_w = .
    replace newtime_w = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft_w) if noncp == 1
    replace newtime_w = timeDFS_w2 if newtime_w == .
    replace newtime_w = admin_w    if trt == 1 & newtime_w > admin_w

    gen dd = dc1  if trt == 1
    replace dd = died if trt == 0
    gen newtime1_w = newtime_w / 52
    stset newtime1_w, failure(dd == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr

    foreach i in 3 {
        predict S0_rpsftm1_w, basesurv
        gen S1_rpsftm1_w    = S0_rpsftm1_w^return(rpsftm_hr_w)
        gen au_S0_rpsftm1_w = lag * S0_rpsftm1_w
        gen au_S1_rpsftm1_w = lag * S1_rpsftm1_w
        egen le_S0_rpsftm1`i'_w = sum(au_S0_rpsftm1_w) if _t <= `i'
        egen le_S1_rpsftm1`i'_w = sum(au_S1_rpsftm1_w) if _t <= `i'
        summarize le_S0_rpsftm1`i'_w
        return scalar le0_rpsftm1`i'_w = r(mean)
        summarize le_S1_rpsftm1`i'_w
        return scalar le1_rpsftm1`i'_w = r(mean)

        sort _t
        gen t1_w = `i' - _t
        replace t1_w = 999999 if _t > `i'
        egen min_t1_w = min(t1_w)
        sum S0_rpsftm1_w if t1_w == min_t1_w
        return scalar S0_rpsftm1_w = r(mean)
        sum S1_rpsftm1_w if t1_w == min_t1_w
        return scalar S1_rpsftm1_w = r(mean)
    }

    /*--- RMST and dSurvival, Weibull, WITHOUT recensoring ---*/

    gen newtime_wr_w = .
    replace newtime_wr_w = noncp_time + (timeDFS_w2 - noncp_time) / return(rpsftm_aft_w) if noncp == 1
    replace newtime_wr_w = timeDFS_w2 if newtime_wr_w == .
    gen newtime1_wr_w = newtime_wr_w / 52
    stset newtime1_wr_w, failure(died == 1) scale(1) id(id)
    sort _t
    drop _t_lag lag
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    stcox trt node_new tumor_new, nohr
    return scalar rpsftm_hr_nore_wr_w    = exp(_b[trtrand])
    return scalar rpsftm_hr_nore_SE_wr_w = _b[trtrand] / return(itt_statistics_w)
    return scalar rpsftm_hr_nore_LB_wr_w = exp(_b[trtrand] - 1.96 * return(rpsftm_hr_nore_SE_wr_w))
    return scalar rpsftm_hr_nore_UB_wr_w = exp(_b[trtrand] + 1.96 * return(rpsftm_hr_nore_SE_wr_w))

    gen cov_rpsftm_hr_nore_wr_w = 0
    replace cov_rpsftm_hr_nore_wr_w = 1 if exp(`trtandcoef') > return(rpsftm_hr_nore_LB_wr_w) & exp(`trtandcoef') < return(rpsftm_hr_nore_UB_wr_w)
    return scalar cov_rpsftm_hr_nore_wr_w = cov_rpsftm_hr_nore_wr_w[1]

    gen power_rpsftm_hr_nore_wr_w = 0
    replace power_rpsftm_hr_nore_wr_w = 1 if 1 > return(rpsftm_hr_nore_LB_wr_w) & 1 < return(rpsftm_hr_nore_UB_wr_w)
    return scalar power_rpsftm_hr_nore_wr_w = power_rpsftm_hr_nore_wr_w[1]

    stcox trt node_new tumor_new, nohr
    foreach i in 3 {
        predict S0_rpsftm1_wr_w, basesurv
        gen S1_rpsftm1_wr_w    = S0_rpsftm1_wr_w^return(rpsftm_hr_nore_wr_w)
        gen au_S0_rpsftm1_wr_w = lag * S0_rpsftm1_wr_w
        gen au_S1_rpsftm1_wr_w = lag * S1_rpsftm1_wr_w
        egen le_S0_rpsftm1`i'_wr_wr_w = sum(au_S0_rpsftm1_wr_w) if _t <= `i'
        egen le_S1_rpsftm1`i'_wr_w    = sum(au_S1_rpsftm1_wr_w) if _t <= `i'
        summarize le_S0_rpsftm1`i'_wr_w
        return scalar le0_rpsftm1`i'_wr_w = r(mean)
        summarize le_S1_rpsftm1`i'_wr_w
        return scalar le1_rpsftm1`i'_wr_w = r(mean)

        sort _t
        gen t1_wr_w = `i' - _t
        replace t1_wr_w = 999999 if _t > `i'
        egen min_t1_wr_w = min(t1_wr_w)
        sum S0_rpsftm1_wr_w if t1_wr_w == min_t1_wr_w
        return scalar S0_rpsftm1_wr_w = r(mean)
        sum S1_rpsftm1_w   if t1_wr_w == min_t1_wr_w
        return scalar S1_rpsftm1_wr_w = r(mean)
    }

end

simulate                                                                         ///
    le0_gp3             = r(le0_gp3)              le1_gp3             = r(le1_gp3)              ///
    S_diff_true         = r(S_diff_true)                                                        ///
    rpsftm_hr           = r(rpsftm_hr)            rpsftm_aft          = r(rpsftm_aft)           ///
    cov_rpsftm_hr       = r(cov_rpsftm_hr)        cov_rpsftm_aft      = r(cov_rpsftm_aft)       ///
    power_rpsftm_hr     = r(power_rpsftm_hr)      power_rpsftm_aft    = r(power_rpsftm_aft)     ///
    le0_rpsftm13        = r(le0_rpsftm13)         le1_rpsftm13        = r(le1_rpsftm13)         ///
    S0_rpsftm1          = r(S0_rpsftm1)           S1_rpsftm1          = r(S1_rpsftm1)           ///
    rpsftm_hr_nore      = r(rpsftm_hr_nore)                                                     ///
    cov_rpsftm_hr_nore  = r(cov_rpsftm_hr_nore)   power_rpsftm_hr_nore= r(power_rpsftm_hr_nore) ///
    le0_rpsftm13_wr     = r(le0_rpsftm13_wr)      le1_rpsftm13_wr     = r(le1_rpsftm13_wr)      ///
    S0_rpsftm1_wr       = r(S0_rpsftm1_wr)        S1_rpsftm1_wr       = r(S1_rpsftm1_wr)        ///
    rpsftm_hr_w         = r(rpsftm_hr_w)          rpsftm_aft_w        = r(rpsftm_aft_w)         ///
    cov_rpsftm_hr_w     = r(cov_rpsftm_hr_w)      cov_rpsftm_aft_w    = r(cov_rpsftm_aft_w)     ///
    power_rpsftm_hr_w   = r(power_rpsftm_hr_w)    power_rpsftm_aft_w  = r(power_rpsftm_aft_w)   ///
    le0_rpsftm13_w      = r(le0_rpsftm13_w)       le1_rpsftm13_w      = r(le1_rpsftm13_w)       ///
    S0_rpsftm1_w        = r(S0_rpsftm1_w)         S1_rpsftm1_w        = r(S1_rpsftm1_w)         ///
    rpsftm_hr_nore_wr_w = r(rpsftm_hr_nore_wr_w)                                                ///
    cov_rpsftm_hr_nore_wr_w   = r(cov_rpsftm_hr_nore_wr_w)                                      ///
    power_rpsftm_hr_nore_wr_w = r(power_rpsftm_hr_nore_wr_w)                                    ///
    le0_rpsftm13_wr_w   = r(le0_rpsftm13_wr_w)    le1_rpsftm13_wr_w   = r(le1_rpsftm13_wr_w)    ///
    S0_rpsftm1_wr_w     = r(S0_rpsftm1_wr_w)      S1_rpsftm1_wr_w     = r(S1_rpsftm1_wr_w),     ///
    reps($N_REPS) saving(simstudy_rpsftm, replace) seed($SEED):                       ///
    simstudy_rpsftm


/*----------------------------------------------------------------------------
  Program 5 of 6: IPCW with logistic weights, selection bias
----------------------------------------------------------------------------*/

capture program drop simstudy_ipcw
program define simstudy_ipcw, rclass

    syntax [, obs(int $N_OBS)             ///
              trtandcoef(real $TRT_COEF) ///
              nodecoef(real $NODE_COEF)           ///
              tumorcoef(real $TUMOR_COEF)          ///
              aftcoef(real $AFT_COEF)           ///
              newvcoeff(real $NEWV_COEF)    ///
              lambda(real $LAMBDA_SEL)      ///
              gammas(real $GAMMA_S)             ///
            ]

    clear
    pr drop _all

    /*--- 1. Generate baseline covariates and survival times ---*/

    rbinary age toxicity new_v tumor_new node_new,                         ///
        means(0.33, 0.48, 0.45, 0.58, 0.28)                                ///
        corr(1    , 0.4 , 0.68, 0   , 0    \                               ///
             0.4  , 1   , 0.68, 0   , 0    \                               ///
             0.68 , 0.68, 1   , 0   , 0    \                               ///
             0    , 0   , 0   , 1   , 0.5  \                               ///
             0    , 0   , 0   , 0.5 , 1)                                   ///
        n($N_OBS)

    set obs `obs'
    gen trtrand = rbinomial(1, 0.5)

    survsim timeDFS_y, distribution(weibull) lambda(`lambda') gammas(`gammas') ///
        cov(trtrand `trtandcoef' node_new `nodecoef' tumor_new `tumorcoef'    ///
            new_v `newvcoeff')

    gen id      = _n
    gen admin   = 6
    gen admin_w = admin * 52
    gen timeDFS_w = timeDFS_y * 52
    replace timeDFS_w = round(timeDFS_w)
    replace timeDFS_w = 1 if timeDFS_w == 0
    replace admin_w   = round(admin_w)
    gen dead = 1 if timeDFS_w <  admin_w
    replace dead = 0 if timeDFS_w >= admin_w
    replace timeDFS_w = admin_w if timeDFS_w > admin_w

    /*--- Non-compliance ---*/

    gen p     = invlogit($NC_INTERCEPT + 0.6 * age + 0.7 * toxicity) if trtrand == 1
    gen noncp = rbinomial(1, p)                            if trtrand == 1
    replace noncp = 0 if noncp == .
    drop p

    /*--- Time of first non-compliance (in weeks) ---*/

    gen noncp_time1 = rbinomial(1, 0.24) if noncp == 1
    gen noncp_time2 = rbinomial(1, 0.34) if noncp == 1 & noncp_time1 == 0
    gen noncp_time3 = rbinomial(1, 0.6)  if noncp == 1 & noncp_time1 == 0 & noncp_time2 == 0
    gen noncp_time  = 1  if noncp_time1 == 1
    replace noncp_time = 4  if noncp_time2 == 1
    replace noncp_time = 11 if noncp_time3 == 1 & noncp == 1
    replace noncp_time = 7  if noncp_time == . & noncp == 1
    replace noncp_time = .  if noncp_time > timeDFS_w & trt == 1 & noncp == 1
    replace noncp      = 0  if noncp_time == . & trt == 1

    /*--- 2. New survival times for non-compliers ---*/

    gen xoOSgainobs = timeDFS_w - noncp_time if noncp == 1
    replace xoOSgainobs = round(xoOSgainobs * `aftcoef')
    gen timeDFS_w2 = cond(noncp == 1, xoOSgainobs + noncp_time, timeDFS_w)

    gen died = 0
    replace died = 1 if timeDFS_w2 < admin_w & dead == 1
    replace timeDFS_w2 = admin_w if timeDFS_w2 > admin_w
    replace xoOSgainobs = . if noncp_time >  admin_w
    replace noncp       = . if noncp_time >= admin_w
    replace noncp_time  = . if noncp_time >= admin_w

    /*--- 3. True dRMST, rRMST and MSD ---*/

    scalar lambda = `lambda'
    scalar gamma  = `gammas'
    scalar beta1  = `trtandcoef'

    gen timeDFS_y1 = timeDFS_y
    replace timeDFS_y1 = 6 if timeDFS_y1 > 6
    stset timeDFS_y1, failure(dead = 1)
    sort _t
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    gen H0      = lambda * _t^(gamma)
    gen H1      = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen au_S0_gp = lag * S0_true
    gen au_S1_gp = lag * S1_true

    foreach i in 3 {
        egen le_S0_gp`i' = sum(au_S0_gp) if _t <= `i'
        egen le_S1_gp`i' = sum(au_S1_gp) if _t <= `i'
        summarize le_S0_gp`i'
        return scalar le0_gp`i' = r(mean)
        summarize le_S1_gp`i'
        return scalar le1_gp`i' = r(mean)
    }

    drop H0 H1 S1_true S0_true
    gen H0 = lambda * 3^(gamma)
    gen H1 = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen S_diff_true = S1_true - S0_true
    sum S_diff_true
    return scalar S_diff_true = r(mean)

    /*=== IPCW (logistic) ===*/

    /*--- Build the time-varying panel (15 scheduled visits per patient,
          trimmed to each patient's follow-up time) ---*/

    expand 15
    sort id
    by id: gen visit = _n
    gen cumweek = 1   if visit == 1
    replace cumweek = 4   if visit == 2
    replace cumweek = 7   if visit == 3
    replace cumweek = 11  if visit == 4
    replace cumweek = 15  if visit == 5
    replace cumweek = 19  if visit == 6
    replace cumweek = 35  if visit == 7
    replace cumweek = 51  if visit == 8
    replace cumweek = 67  if visit == 9
    replace cumweek = 91  if visit == 10
    replace cumweek = 115 if visit == 11
    replace cumweek = 139 if visit == 12
    replace cumweek = 163 if visit == 13
    replace cumweek = 215 if visit == 14
    replace cumweek = 267 if visit == 15

    drop if cumweek > timeDFS_w2
    sort id
    by id: gen firstobs = 0
    by id: replace firstobs = 1 if _n == 1
    by id: gen finalobs = 0
    by id: replace finalobs = 1 if _n == _N
    by id: drop if cumweek > noncp_time
    by id: replace finalobs = 0
    by id: replace finalobs = 1 if _n == _N

    gen censOS = 0
    replace censOS = 1 if finalobs == 1 & died == 0 & visit == 15 & timeDFS_w2 == 312
    sum censOS
    replace died = . if censOS == 1

    gen compliance_tdo = 0
    replace compliance_tdo = 1 if trtrand == 1 & noncp_time == cumweek & noncp == 1
    by id: replace died = . if (compliance_tdo == 1 | censOS == 1)

    gen death_tdo = .
    sort id
    by id: replace death_tdo = 0 if finalobs == 0
    by id: replace death_tdo = 1 if died == 1 & finalobs == 1
    by id: replace death_tdo = . if compliance_tdo == 1 & finalobs == 1 & died == 1

    gen cumweek_entry = .
    sort id
    by id: replace cumweek_entry = 0              if _n == 1
    sort id
    by id: replace cumweek_entry = cumweek[_n - 1] if _n != 1

    /*--- Logistic weights ---*/

    xi: logistic compliance_tdo age i.toxicity if trtrand == 1
    predict p1 if e(sample)
    gen pa = p1 * compliance_tdo + (1 - p1) * (1 - compliance_tdo)
    replace pa = 1 if pa == .
    sort id cumweek
    by id: replace pa = pa * pa[_n - 1] if _n != 1
    rename pa p_denom

    xi: logistic compliance_tdo age if trtrand == 1
    predict p2 if e(sample)
    gen pa2 = p2 * compliance_tdo + (1 - p2) * (1 - compliance_tdo)
    replace pa2 = 1 if pa2 == .
    sort id cumweek
    by id: replace pa2 = pa2 * pa2[_n - 1] if _n != 1
    rename pa2 p_num

    gen weight  = 1     / p_denom
    gen sweight = p_num / p_denom
    replace sweight = 1 if trt == 0
    replace weight  = 1 if trt == 0

    /*--- Stabilised weight distribution ---*/

    summarize sweight if noncp == 1, detail
    return scalar sw_mean1   = r(mean)
    return scalar sw_median1 = r(p50)
    return scalar sw_se1     = r(sd)
    return scalar sw_min1    = r(min)
    return scalar sw_max1    = r(max)

    summarize sweight if noncp != 1, detail
    return scalar sw_mean   = r(mean)
    return scalar sw_median = r(p50)
    return scalar sw_se     = r(sd)
    return scalar sw_min    = r(min)
    return scalar sw_max    = r(max)

    summarize sweight if trt == 1, detail
    return scalar sw_mean2   = r(mean)
    return scalar sw_median2 = r(p50)
    return scalar sw_se2     = r(sd)
    return scalar sw_min2    = r(min)
    return scalar sw_max2    = r(max)

    summarize sweight if trt == 0, detail
    return scalar sw_mean3   = r(mean)
    return scalar sw_median3 = r(p50)
    return scalar sw_se3     = r(sd)
    return scalar sw_min3    = r(min)
    return scalar sw_max3    = r(max)

    /*--- Weighted Cox model ---*/

    stset cumweek [pweight = sweight], failure(death_tdo == 1) ///
        entry(cumweek_entry) exit(cumweek)
    stcox trtrand node_new tumor_new, vce(robust)

    return scalar ipcw_cox_hr    = exp(_b[trtrand])
    return scalar ipcw_cox_hr_SE = exp(_b[trtrand]) * _se[trtrand]
    return scalar ipcw_cox_hr_LB = exp(_b[trtrand] - 1.96 * _se[trtrand])
    return scalar ipcw_cox_hr_UB = exp(_b[trtrand] + 1.96 * _se[trtrand])

    gen cov_ipcw_cox_hr = 0
    replace cov_ipcw_cox_hr = 1 if exp(`trtandcoef') > return(ipcw_cox_hr_LB) & exp(`trtandcoef') < return(ipcw_cox_hr_UB)
    return scalar cov_ipcw_cox_hr = cov_ipcw_cox_hr[1]

    gen power_ipcw_cox_hr = 0
    replace power_ipcw_cox_hr = 1 if 1 > return(ipcw_cox_hr_LB) & 1 < return(ipcw_cox_hr_UB)
    return scalar power_ipcw_cox_hr = power_ipcw_cox_hr[1]

    /*--- RMST and dSurvival ---*/

    collapse (max) trtrand timeDFS_w2 node_new tumor_new new_v xoOSgainobs ///
        noncp_time died, by(id)
    by id: replace noncp_time = 0 if noncp_time == .

    preserve
    stset timeDFS_w2, failure(died) id(id)
    xi: streg node_new tumor_new if trtrand == 0, distribution(weibull) nohr

    gen cons = _b[_cons]
    gen lnp  = _b[/ln_p]

    expand 2
    sort id
    gen time_w = 0
    replace time_w = _n * (312 / _N) if _n > 1
    gen time_y = time_w / 52

    sort time_w
    gen survf = 1
    replace survf = exp(-exp(cons) * time_w^exp(lnp)) if _n > 1
    gen hazf = 0
    replace hazf = 1 - survf / survf[_n - 1] if _n > 1
    gen hazfexp = hazf * return(ipcw_cox_hr)
    gen survfexp = 1
    replace survfexp = survfexp[_n - 1] * (1 - hazfexp) if _n > 1

    gen S0_ipcw = survf
    gen S1_ipcw = survfexp

    gen lag = time_y - time_y[_n - 1] if _n > 1
    replace lag = time_y if _n == 1
    gen au_S0_ipcw = lag * survf
    gen au_S1_ipcw = lag * survfexp

    foreach i in 3 {
        egen le_S0_ipcw`i' = sum(au_S0_ipcw) if _t <= `i'
        egen le_S1_ipcw`i' = sum(au_S1_ipcw) if _t <= `i'
        summarize le_S0_ipcw`i'
        return scalar le0_ipcw`i' = r(mean)
        summarize le_S1_ipcw`i'
        return scalar le1_ipcw`i' = r(mean)

        sort _t
        gen t1 = `i' - _t
        replace t1 = 999999 if _t > `i'
        egen min_t1 = min(t1)
        sum S0_ipcw if t1 == min_t1
        return scalar S0_ipcw1 = r(mean)
        sum S1_ipcw if t1 == min_t1
        return scalar S1_ipcw1 = r(mean)
    }

end

* NOTE: in the original code this simulate command included references to
* le03/le13/S0/S1 that are not returned by the program. They have been removed
* so the call runs. If you intended to return those scalars, add the matching
* `return scalar` lines inside the program.
simulate                                                              ///
    sw_mean1   = r(sw_mean1)   sw_median1 = r(sw_median1)              ///
    sw_se1     = r(sw_se1)     sw_min1    = r(sw_min1)   sw_max1 = r(sw_max1) ///
    sw_mean    = r(sw_mean)    sw_median  = r(sw_median)                ///
    sw_se      = r(sw_se)      sw_min     = r(sw_min)    sw_max  = r(sw_max)  ///
    sw_mean2   = r(sw_mean2)   sw_median2 = r(sw_median2)               ///
    sw_se2     = r(sw_se2)     sw_min2    = r(sw_min2)   sw_max2 = r(sw_max2) ///
    sw_mean3   = r(sw_mean3)   sw_median3 = r(sw_median3)               ///
    sw_se3     = r(sw_se3)     sw_min3    = r(sw_min3)   sw_max3 = r(sw_max3) ///
    le0_gp3    = r(le0_gp3)    le1_gp3    = r(le1_gp3)                  ///
    ipcw_cox_hr= r(ipcw_cox_hr)                                          ///
    cov_ipcw_cox_hr   = r(cov_ipcw_cox_hr)                               ///
    power_ipcw_cox_hr = r(power_ipcw_cox_hr)                             ///
    le0_ipcw3  = r(le0_ipcw3)  le1_ipcw3  = r(le1_ipcw3)                ///
    S0_ipcw1   = r(S0_ipcw1)   S1_ipcw1   = r(S1_ipcw1),                ///
    reps($N_REPS) saving(simstudy_ipcw, replace) seed($SEED): ///
    simstudy_ipcw


/*----------------------------------------------------------------------------
  Program 6 of 6: IPCW with Cox weights, selection bias
----------------------------------------------------------------------------*/

capture program drop simstudyipcw_cox
program define simstudyipcw_cox, rclass

    syntax [, obs(int $N_OBS)             ///
              trtandcoef(real $TRT_COEF) ///
              nodecoef(real $NODE_COEF)           ///
              tumorcoef(real $TUMOR_COEF)          ///
              aftcoef(real $AFT_COEF)           ///
              newvcoeff(real $NEWV_COEF)    ///
              lambda(real $LAMBDA_SEL)      ///
              gammas(real $GAMMA_S)             ///
            ]

    clear
    pr drop _all

    /*--- 1. Generate baseline covariates and survival times ---*/

    rbinary age toxicity new_v tumor_new node_new,                         ///
        means(0.33, 0.48, 0.45, 0.58, 0.28)                                ///
        corr(1    , 0.4 , 0.68, 0   , 0    \                               ///
             0.4  , 1   , 0.68, 0   , 0    \                               ///
             0.68 , 0.68, 1   , 0   , 0    \                               ///
             0    , 0   , 0   , 1   , 0.5  \                               ///
             0    , 0   , 0   , 0.5 , 1)                                   ///
        n($N_OBS)

    set obs `obs'
    gen trtrand = rbinomial(1, 0.5)

    survsim timeDFS_y, distribution(weibull) lambda(`lambda') gammas(`gammas') ///
        cov(trtrand `trtandcoef' node_new `nodecoef' tumor_new `tumorcoef'    ///
            new_v `newvcoeff')

    gen id      = _n
    gen admin   = 6
    gen admin_w = admin * 52
    gen timeDFS_w = timeDFS_y * 52
    replace timeDFS_w = round(timeDFS_w)
    replace timeDFS_w = 1 if timeDFS_w == 0
    replace admin_w   = round(admin_w)
    gen dead = 1 if timeDFS_w <  admin_w
    replace dead = 0 if timeDFS_w >= admin_w
    replace timeDFS_w = admin_w if timeDFS_w > admin_w

    /*--- Non-compliance ---*/

    gen p     = invlogit($NC_INTERCEPT + 0.6 * age + 0.7 * toxicity) if trtrand == 1
    gen noncp = rbinomial(1, p)                            if trtrand == 1
    replace noncp = 0 if noncp == .
    drop p

    /*--- Time of first non-compliance (in weeks) ---*/

    gen noncp_time1 = rbinomial(1, 0.24) if noncp == 1
    gen noncp_time2 = rbinomial(1, 0.34) if noncp == 1 & noncp_time1 == 0
    gen noncp_time3 = rbinomial(1, 0.6)  if noncp == 1 & noncp_time1 == 0 & noncp_time2 == 0
    gen noncp_time  = 1  if noncp_time1 == 1
    replace noncp_time = 4  if noncp_time2 == 1
    replace noncp_time = 11 if noncp_time3 == 1 & noncp == 1
    replace noncp_time = 7  if noncp_time == . & noncp == 1
    replace noncp_time = .  if noncp_time > timeDFS_w & trt == 1 & noncp == 1
    replace noncp      = 0  if noncp_time == . & trt == 1

    /*--- 2. New survival times for non-compliers ---*/

    gen xoOSgainobs = timeDFS_w - noncp_time if noncp == 1
    replace xoOSgainobs = round(xoOSgainobs * `aftcoef')
    gen timeDFS_w2 = cond(noncp == 1, xoOSgainobs + noncp_time, timeDFS_w)

    gen died = 0
    replace died = 1 if timeDFS_w2 < admin_w & dead == 1
    replace timeDFS_w2 = admin_w if timeDFS_w2 > admin_w
    replace xoOSgainobs = . if noncp_time >  admin_w
    replace noncp       = . if noncp_time >= admin_w
    replace noncp_time  = . if noncp_time >= admin_w

    /*--- 3. True dRMST, rRMST and MSD ---*/

    scalar lambda = `lambda'
    scalar gamma  = `gammas'
    scalar beta1  = `trtandcoef'

    gen timeDFS_y1 = timeDFS_y
    replace timeDFS_y1 = 6 if timeDFS_y1 > 6
    stset timeDFS_y1, failure(dead = 1)
    sort _t
    gen _t_lag = _t[_n - 1]
    gen lag    = _t - _t_lag
    replace lag = _t if lag == .

    gen H0      = lambda * _t^(gamma)
    gen H1      = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen au_S0_gp = lag * S0_true
    gen au_S1_gp = lag * S1_true

    foreach i in 3 {
        egen le_S0_gp`i' = sum(au_S0_gp) if _t <= `i'
        egen le_S1_gp`i' = sum(au_S1_gp) if _t <= `i'
        summarize le_S0_gp`i'
        return scalar le0_gp`i' = r(mean)
        summarize le_S1_gp`i'
        return scalar le1_gp`i' = r(mean)
    }

    drop H0 H1 S1_true S0_true
    gen H0 = lambda * 3^(gamma)
    gen H1 = H0 * exp(beta1)
    gen S0_true = exp(-H0)
    gen S1_true = exp(-H1)
    gen S_diff_true = S1_true - S0_true
    sum S_diff_true
    return scalar S_diff_true = r(mean)

    /*=== IPCW (Cox-based weights) ===*/

    /*--- Build the time-varying panel (15 scheduled visits per patient,
          trimmed to each patient's follow-up time) ---*/

    expand 15
    sort id
    by id: gen visit = _n
    gen cumweek = 1   if visit == 1
    replace cumweek = 4   if visit == 2
    replace cumweek = 7   if visit == 3
    replace cumweek = 11  if visit == 4
    replace cumweek = 15  if visit == 5
    replace cumweek = 19  if visit == 6
    replace cumweek = 35  if visit == 7
    replace cumweek = 51  if visit == 8
    replace cumweek = 67  if visit == 9
    replace cumweek = 91  if visit == 10
    replace cumweek = 115 if visit == 11
    replace cumweek = 139 if visit == 12
    replace cumweek = 163 if visit == 13
    replace cumweek = 215 if visit == 14
    replace cumweek = 267 if visit == 15

    drop if cumweek > timeDFS_w2
    sort id
    by id: gen firstobs = 0
    by id: replace firstobs = 1 if _n == 1
    by id: gen finalobs = 0
    by id: replace finalobs = 1 if _n == _N
    by id: drop if cumweek > noncp_time
    by id: replace finalobs = 0
    by id: replace finalobs = 1 if _n == _N

    gen censOS = 0
    replace censOS = 1 if finalobs == 1 & died == 0 & visit == 15 & timeDFS_w2 == 312
    sum censOS
    replace died = . if censOS == 1

    gen compliance_tdo = 0
    replace compliance_tdo = 1 if trtrand == 1 & noncp_time == cumweek & noncp == 1
    by id: replace died = . if (compliance_tdo == 1 | censOS == 1)

    gen death_tdo = .
    sort id
    by id: replace death_tdo = 0 if finalobs == 0
    by id: replace death_tdo = 1 if died == 1 & finalobs == 1
    by id: replace death_tdo = . if compliance_tdo == 1 & finalobs == 1 & died == 1

    gen cumweek_entry = .
    sort id
    by id: replace cumweek_entry = 0              if _n == 1
    sort id
    by id: replace cumweek_entry = cumweek[_n - 1] if _n != 1

    /*--- Denominator (Cz) ---*/

    stset cumweek if trt == 1, failure(compliance_tdo == 1) enter(time cumweek_entry)
    stcox age i.toxicity if trt == 1
    predict xb_cz, xb
    gen exp_xb_cz = exp(xb_cz)
    predict basech_cz, basechazard
    by id: gen surv_cz = exp(-basech_cz[_n] * exp_xb_cz[_n])
    stset cumweek if trt == 1, failure(compliance_tdo == 1) enter(time cumweek_entry)

    /*--- Numerator (C0) ---*/

    stcox age if trt == 1
    predict xb_c0, xb
    gen exp_xb_c0 = exp(xb_c0)
    predict basech_c0, basechazard
    by id: gen surv_c0 = exp(-basech_c0[_n] * exp_xb_c0[_n])

    gen sweight = surv_c0 / surv_cz
    replace sweight = 1 if trt == 0

    /*--- Stabilised weight distribution ---*/

    summarize sweight if noncp == 1, detail
    return scalar sw_mean1   = r(mean)
    return scalar sw_median1 = r(p50)
    return scalar sw_se1     = r(sd)
    return scalar sw_min1    = r(min)
    return scalar sw_max1    = r(max)

    summarize sweight if noncp != 1, detail
    return scalar sw_mean   = r(mean)
    return scalar sw_median = r(p50)
    return scalar sw_se     = r(sd)
    return scalar sw_min    = r(min)
    return scalar sw_max    = r(max)

    summarize sweight if trt == 1, detail
    return scalar sw_mean2   = r(mean)
    return scalar sw_median2 = r(p50)
    return scalar sw_se2     = r(sd)
    return scalar sw_min2    = r(min)
    return scalar sw_max2    = r(max)

    summarize sweight if trt == 0, detail
    return scalar sw_mean3   = r(mean)
    return scalar sw_median3 = r(p50)
    return scalar sw_se3     = r(sd)
    return scalar sw_min3    = r(min)
    return scalar sw_max3    = r(max)

    /*--- Weighted Cox model ---*/

    stset cumweek [pweight = sweight], failure(death_tdo == 1) ///
        entry(cumweek_entry) exit(cumweek)
    stcox trtrand node_new tumor_new, vce(robust)

    return scalar ipcw_cox_hr    = exp(_b[trtrand])
    return scalar ipcw_cox_hr_SE = exp(_b[trtrand]) * _se[trtrand]
    return scalar ipcw_cox_hr_LB = exp(_b[trtrand] - 1.96 * _se[trtrand])
    return scalar ipcw_cox_hr_UB = exp(_b[trtrand] + 1.96 * _se[trtrand])

    gen cov_ipcw_cox_hr = 0
    replace cov_ipcw_cox_hr = 1 if exp(`trtandcoef') > return(ipcw_cox_hr_LB) & exp(`trtandcoef') < return(ipcw_cox_hr_UB)
    return scalar cov_ipcw_cox_hr = cov_ipcw_cox_hr[1]

    gen power_ipcw_cox_hr = 0
    replace power_ipcw_cox_hr = 1 if 1 > return(ipcw_cox_hr_LB) & 1 < return(ipcw_cox_hr_UB)
    return scalar power_ipcw_cox_hr = power_ipcw_cox_hr[1]

    /*--- RMST and dSurvival ---*/

    collapse (max) trtrand timeDFS_w2 node_new tumor_new new_v xoOSgainobs ///
        noncp_time died, by(id)
    by id: replace noncp_time = 0 if noncp_time == .

    preserve
    stset timeDFS_w2, failure(died) id(id)
    xi: streg node_new tumor_new if trtrand == 0, distribution(weibull) nohr

    gen cons = _b[_cons]
    gen lnp  = _b[/ln_p]

    expand 2
    sort id
    gen time_w = 0
    replace time_w = _n * (312 / _N) if _n > 1
    gen time_y = time_w / 52

    sort time_w
    gen survf = 1
    replace survf = exp(-exp(cons) * time_w^exp(lnp)) if _n > 1
    gen hazf = 0
    replace hazf = 1 - survf / survf[_n - 1] if _n > 1
    gen hazfexp = hazf * return(ipcw_cox_hr)
    gen survfexp = 1
    replace survfexp = survfexp[_n - 1] * (1 - hazfexp) if _n > 1

    gen S0_ipcw = survf
    gen S1_ipcw = survfexp

    gen lag = time_y - time_y[_n - 1] if _n > 1
    replace lag = time_y if _n == 1
    gen au_S0_ipcw = lag * survf
    gen au_S1_ipcw = lag * survfexp

    foreach i in 3 {
        egen le_S0_ipcw`i' = sum(au_S0_ipcw) if _t <= `i'
        egen le_S1_ipcw`i' = sum(au_S1_ipcw) if _t <= `i'
        summarize le_S0_ipcw`i'
        return scalar le0_ipcw`i' = r(mean)
        summarize le_S1_ipcw`i'
        return scalar le1_ipcw`i' = r(mean)

        sort _t
        gen t1 = `i' - _t
        replace t1 = 999999 if _t > `i'
        egen min_t1 = min(t1)
        sum S0_ipcw if t1 == min_t1
        return scalar S0_ipcw1 = r(mean)
        sum S1_ipcw if t1 == min_t1
        return scalar S1_ipcw1 = r(mean)
    }

end

simulate                                                              ///
    sw_mean1   = r(sw_mean1)   sw_median1 = r(sw_median1)              ///
    sw_se1     = r(sw_se1)     sw_min1    = r(sw_min1)   sw_max1 = r(sw_max1) ///
    sw_mean    = r(sw_mean)    sw_median  = r(sw_median)                ///
    sw_se      = r(sw_se)      sw_min     = r(sw_min)    sw_max  = r(sw_max)  ///
    sw_mean2   = r(sw_mean2)   sw_median2 = r(sw_median2)               ///
    sw_se2     = r(sw_se2)     sw_min2    = r(sw_min2)   sw_max2 = r(sw_max2) ///
    sw_mean3   = r(sw_mean3)   sw_median3 = r(sw_median3)               ///
    sw_se3     = r(sw_se3)     sw_min3    = r(sw_min3)   sw_max3 = r(sw_max3) ///
    le0_gp3    = r(le0_gp3)    le1_gp3    = r(le1_gp3)                  ///
    ipcw_cox_hr= r(ipcw_cox_hr)                                          ///
    cov_ipcw_cox_hr   = r(cov_ipcw_cox_hr)                               ///
    power_ipcw_cox_hr = r(power_ipcw_cox_hr)                             ///
    le0_ipcw3  = r(le0_ipcw3)  le1_ipcw3  = r(le1_ipcw3)                ///
    S0_ipcw1   = r(S0_ipcw1)   S1_ipcw1   = r(S1_ipcw1),                ///
    reps($N_REPS) saving(simstudyipcw_cox, replace) seed($SEED): ///
    simstudyipcw_cox

	
	
	
	
/*--- Carry-over effect -----------------------------------------------------*/
/* aftcoef = 1- (1-aftcoef)*0.9
 

/*--- Non-compliance intercept variants -------------------------------------*/
/* gen p = invlogit( 0.4 + 0.6*age + 0.7*toxicity) if trtrand == 1           */
/* changing the beta0



/* End of do-file ----------------------------------------------------------- */





