# Version 1.1 released by David Romps on August 25, 2022.
# 	(Added the bloodflow routine)
# Version 1.0 released by David Romps on May 12, 2022.
# 	(Initial release)
# 
# When using this code, please cite:
# 
# @article{20heatindex,
#   Title   = {Extending the Heat Index},
#   Author  = {Yi-Chuan Lu and David M. Romps},
#   Journal = {Journal of Applied Meteorology and Climatology},
#   Year    = {2022},
#   Volume  = {61},
#   Number  = {10},
#   Pages   = {1367--1383},
#   Year    = {2022},
# }
#
# This headindex function returns the Heat Index in Kelvin. The inputs are:
# - T, the temperature in Kelvin
# - rh, the relative humidity, which is a value from 0 to 1
# - verbose is an optional logical flag. If true, the function returns the physiological state.
# - debug is an optional logical flag. If true, code checks its answer.

common_variables <- list(
   Ttrip   = 273.16,    # K
   ptrip   = 611.65,    # Pa
   E0v     = 2.3740e6,  # J/kg
   E0s     = 0.3337e6,  # J/kg
   rgasa   = 287.04,    # J/kg/K 
   rgasv   = 461,       # J/kg/K 
   cva     = 719,       # J/kg/K
   cvv     = 1418,      # J/kg/K 
   cvl     = 4119,      # J/kg/K 
   cvs     = 1861,      # J/kg/K 
   cpa     = 719 + 287.04,  # J/kg/K
   cpv     = 1418 + 461,    # J/kg/K
   L       = 2.3740e6 + (1418-4119)*(310-273.16) + 461*310,  # J/kg
   pa0     = 1.6e3,  # Pa, reference pressure for Heat Index
   phi0    = 0.84,   # Boundary between regions I and II
   Mc      = 83.6,   # kg, fryar2018 mean adult mass
   H       = 1.69,   # m, fryar2018 mean adult height
   cpc     = 3492,   # J/kg/K, heat capacity of core at constant pressure, gagge1972 
   A       = 0.202*(83.6^0.425)*(1.69^0.725),  # m^2, equation 8.2 of parsons2014
   Cc      = 83.6*3492/(0.202*(83.6^0.425)*(1.69^0.725)),  # J / kg / K / m^2, fryar2018 mean adult mass, gagge1972 specific heat
   Q       = 180,  # W / m^2
   phisalt = 0.9,
   p       = 1.013e5,  # Pa
   r       = 124,  # Pa / K
   eta     = 1.43e-6,  # kg / J
   sigma   = 5.67e-8,  # W / m^2 / K^4
   Tcfatal = 315   # K
)

pvstar <- Vectorize(function(T) {
   if (T > with(common_variables,Ttrip)) {
      return(with(common_variables, ptrip * (T/Ttrip)^((cpv-cvl)/rgasv) *
         exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/T) ) ))
   } else {
      return(with(common_variables, ptrip * (T/Ttrip)^((cpv-cvs)/rgasv) *
         exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1/Ttrip - 1/T) ) ))
   }
})

clothed_variables <- list(
   dTcdt = 0,   # K / s
   Tc = 310,   # K
   Rs = 0.0387,   # m^2 K / W
   Zs = 52.1,   # m^2 Pa / W
   hc = 17.4,   # W / K / m^2
   hcbar = 11.6,   # W / K / m^2
   phirad = 0.85,
   phiradbar = 0.79,
   Za = 60.6/17.4,   # Pa m^2 / W
   Zabar = 60.6/11.6,   # Pa m^2 / W
   pc = with(common_variables,phisalt*pvstar(310)),  # Pa
   epsilon = 0.97,
   epsilonbar = 0.97
)

covering_variables <- list(
   dTcdt = 0,   # K / s
   Tc = 310,   # K
   Rs = 0.0387,   # m^2 K / W
   Zs = 52.1,   # m^2 K / W
   hc = 17.4,   # W / K / m^2
   phirad = 0.85,
   Za = 60.6/17.4,   # Pa m^2 / W
   pc = with(common_variables,phisalt*pvstar(310)),  # Pa
   psbar = with(common_variables,phisalt*pvstar(310)),  # Pa
   Tsbar = 310,   # K
   epsilon = 0.97,
   Rf = Inf,
   Zf = Inf
)

naked_variables <- list(
   dTcdt = 0,   # K / s
   Tc = 310,   # K
   pc = with(common_variables,phisalt*pvstar(310)),  # Pa
   hc = 12.3,   # W / K / m^2
   Za = 60.6/12.3,   # m^2 Pa / W
   phirad = 0.8,
   epsilon = 0.97,
   phi = 0,
   Rf = 0,
   Zf = 0
)

warming_variables <- list(
   Tc = 310,   # K
   pc = with(common_variables,phisalt*pvstar(310)),  # Pa
   hc = 12.3,   # W / K / m^2
   Za = 60.6/12.3,   # m^2 Pa / W
   phirad = 0.8,
   epsilon = 0.97,
   phi = 0,
   Rf = 0,
   Zf = 0,
   Rs = 0,
   Zs = 0,
   ps = with(common_variables,phisalt*pvstar(310)),  # Pa
   psbar = with(common_variables,phisalt*pvstar(310)),  # Pa
   pf = with(common_variables,phisalt*pvstar(310)),  # Pa
   Ts = 310,
   Tsbar = 310,
   Tf = 310
)


find_covering_state <- function(Ta,rha) {

   if (rha < 0 | rha > 1) {
      stop('The rha provided as input is outside of bounds')
   }

   # Get parameters
   params <- c(common_variables,covering_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa <- rha*pvstar(Ta)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   
   # Solve for ps
   ps <- (Za*pc + Zs*pa) / (Zs + Za)
   
   # Solve for Ts
   error <- function(Ts) {
      Ra <- 1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta) + hc)
      return((Tc-Ts)/Rs-(Ts-Ta)/Ra-(ps-pa)/Za)
   }
   Ts <- uniroot(error,c(100,10000))$root
   Ra <- 1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta) + hc)

   # Solve for phi
   phi <- 1 - Rs*(Q-Qv)/(Tc-Ts)

   Tf <- Ta
   pf <- pa

   state <- list(
      dTcdt=dTcdt,Q=Q,Qv=Qv,
      Rs=Rs,Rf=Rf,Ra=Ra,
      Zs=Zs,Zf=Zf,Za=Za,
      phi=phi,
      pc=pc,ps=ps,psbar=psbar,pf=pf,pa=pa,
      Tc=Tc,Ts=Ts,Tsbar=Tsbar,Tf=Tf,Ta=Ta )

   return(state)

}


find_clothed_state <- function(Ta,rha) {

   if (rha < 0 | rha > 1) {
      stop('The rha provided as input is outside of bounds')
   }

   pa <- rha*pvstar(Ta)

   # Get parameters
   params <- c(common_variables,clothed_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])
   
   phi <- phi0
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   
   # Solve for ps
   ps <- (Za*pc + Zs*pa) / (Zs + Za)
   
   # Solve for Ts
   Ra <- function(Ts) {
      return(1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta) + hc))
   }
   error <- function(Ts) {
      return((Tc-Ts)/Rs - (Ts-Ta)/Ra(Ts) - (ps-pa)/Za)
   }
   Ts <- uniroot(error,c(0,1e4),tol=1e-6)$root
   Ra <- Ra(Ts)
   
   # Solve for Tsbar
   Tsbar <- Tc + (Rs/phi)*(Qv - Q + (1-phi)*(Tc-Ts)/Rs)

   # Solve for Tf
   error <- function(Tf) {
      # Solve for Rabar, Rf, Zf, psbar and pf in terms of Tf
      if (Tf==Ta) {
        return((Tc-Tsbar)/Rs)
      } else {
         Rabar <- 1/( phiradbar*epsilonbar*sigma*(Tf^2+Ta^2)*(Tf+Ta) + hcbar )
         Rf <- Rabar*(Tsbar-Tf)/(Tf-Ta)
         psbar <- pc + Zs*(pa-pc)/(Zabar+r*Rf+Zs)
         pf <- (Zabar*pc+(r*Rf+Zs)*pa)/(Zabar+r*Rf+Zs)
         Zf <- r*Rf
         # Return the net cooling of the clothed skin
         return((Tc-Tsbar)/Rs - (Tsbar-Tf)/Rf - (psbar-pf)/Zf)
      }
   }
   if (sign(error(Ta)) != sign(error(1e4))) {
      Tf <- uniroot(error,c(Ta,1e4),tol=1e-6)$root
   } else {
      Tf <- uniroot(error,c(-100,Ta),tol=1e-6)$root
   }

   Rabar <- 1/( phiradbar*epsilonbar*sigma*(Tf^2+Ta^2)*(Tf+Ta) + hcbar )
   Rf <- Rabar*(Tsbar-Tf)/(Tf-Ta)
   psbar <- pc + Zs*(pa-pc)/(Zabar+r*Rf+Zs)
   pf <- (Zabar*pc+(r*Rf+Zs)*pa)/(Zabar+r*Rf+Zs)
   Zf <- r*Rf

   state <- list(
      dTcdt=dTcdt,Q=Q,Qv=Qv,
      Rs=Rs,Rf=Rf,Ra=Ra,
      Zs=Zs,Zf=Zf,Za=Za,
      phi=phi,
      pc=pc,ps=ps,psbar=psbar,pf=pf,pa=pa,
      Tc=Tc,Ts=Ts,Tsbar=Tsbar,Tf=Tf,Ta=Ta )

   return(state)

}


find_naked_state <- function(Ta,rha) {

   if (rha < 0 | rha > 1) {
      stop('The rha provided as input is outside of bounds')
   }

   # Get parameters
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa <- rha*pvstar(Ta)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))

   # Calculate Ts
   error <- Vectorize(function(Ts) {
      Ra <- 1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta)+hc)
      # Solve for Rs
      Rs <- (Tc-Ts)/(Q-Qv)
      Zs <- 6.0e8*Rs^5
      # Solve for ps
      ps <- min((Za*pc + Zs*pa)/(Za+Zs),phisalt*pvstar(Ts))
      # Calculate error
      error <- Q-Qv - (Ts-Ta)/Ra - (ps-pa)/Za
      return(error)
   })
   if (sign(error(0)) != sign(error(Tc))) {
      Ts <- uniroot(error,c(0,Tc),tol=1e-6)$root
   } else {
      Ts <- uniroot(error,c(Tc,1e4),tol=1e-6)$root
   }
   Ra <- 1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta)+hc)
   # Solve for Rs
   Rs <- (Tc-Ts)/(Q-Qv)
   Zs <- 6.0e8*Rs^5
   # Solve for ps
   ps <- min(phisalt*pvstar(Ts),(Za*pc + Zs*pa)/(Za+Zs))

   psbar <- ps
   Tsbar <- Ts
   pf <- ps
   Tf <- Ts

   state <- list(
      dTcdt=dTcdt,Q=Q,Qv=Qv,
      Rs=Rs,Rf=Rf,Ra=Ra,
      Zs=Zs,Zf=Zf,Za=Za,
      phi=phi,
      pc=pc,ps=ps,psbar=psbar,pf=pf,pa=pa,
      Tc=Tc,Ts=Ts,Tsbar=Tsbar,Tf=Tf,Ta=Ta )

   return(state)

}


find_warming_state <- function(Ta,rha) {

   if (rha < 0 | rha > 1) {
      stop('The rha provided as input is outside of bounds')
   }

   # Get parameters
   params <- c(common_variables,warming_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa <- rha*pvstar(Ta)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)

   dTcdt <- (Q-Qv-(Tc-Ta)/Ra-(pc-pa)/Za)/Cc

   state <- list(
      dTcdt=dTcdt,Q=Q,Qv=Qv,
      Rs=Rs,Rf=Rf,Ra=Ra,
      Zs=Zs,Zf=Zf,Za=Za,
      phi=phi,
      pc=pc,ps=ps,psbar=psbar,pf=pf,pa=pa,
      Tc=Tc,Ts=Ts,Tsbar=Tsbar,Tf=Tf,Ta=Ta )

   return(state)

}


find_covering_Ta <- function(phi,pa=pa0) {
   if (phi < 0 | phi > find_covering_state(0,0)$phi) {
      stop('The phi provided as input is outside of bounds')
   }
   error <- function(Ta) {
      rh <- min(pa/pvstar(Ta),1)
      rh[which(Ta==0)] <- 0
      return(find_covering_state(Ta,rh)$phi - phi)
   }
   Ta <- uniroot(error,c(0,300),tol=1e-6)$root
   return(Ta)
}


find_clothed_Rfinfty_Ta <- function(pa=pa0) {

   # Get parameters
   params <- c(common_variables,clothed_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa_requested <- pa
   phi <- phi0

   # Solve for Ta
   error <- function(Ta) {
      pa <- min(pvstar(Ta),pa_requested)
      # Solve for ps
      ps <- (Za*pc + Zs*pa) / (Zs + Za)
      # Solve for Ts
      Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(p*rgasv)*(pc-pa))
      Ts <- Tc - Rs*(Q-Qv)/(1-phi)
      Ra <- 1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta) + hc)
      return((Tc-Ts)/Rs - (Ts-Ta)/Ra - (ps-pa)/Za)
   }
   Ta <- uniroot(error,c(100,10000),tol=1e-6)$root
   return(Ta)

}


find_clothed_Rfzero_Ta <- function(pa=pa0) {

   # Get parameters
   params <- c(common_variables,clothed_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa_requested <- pa
   phi <- phi0

   # Solve for Ta
   error <- function(Ta) {
      pa <- min(pvstar(Ta),pa_requested)
      # Solve for ps and pf
      ps <- (Za   *pc + Zs*pa) / (Zs + Za   )
      pf <- (Zabar*pc + Zs*pa) / (Zs + Zabar)
      # Define Ra and Rabar
      Ra <- function(Ts) {
         return(1/(phirad   *epsilon   *sigma*(Ts^2+Ta^2)*(Ts+Ta) + hc   ))
      }
      Rabar <- function(Tf) {
         return(1/(phiradbar*epsilonbar*sigma*(Tf^2+Ta^2)*(Tf+Ta) + hcbar))
      }
      # Solve for Ts
      error2 <- function(Ts) {
         return((Tc-Ts)/Rs - (Ts-Ta)/Ra   (Ts) - (ps-pa)/Za   )
      }
      Ts <- uniroot(error2,c(100,10000),tol=1e-6)$root
      # Solve for Tf
      error2 <- function(Tf) {
         return((Tc-Tf)/Rs - (Tf-Ta)/Rabar(Tf) - (pf-pa)/Zabar)
      }
      Tf <- uniroot(error2,c(100,10000),tol=1e-6)$root
      Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(p*rgasv)*(pc-pa))
      return(Q-Qv-phi*(Tc-Tf)/Rs-(1-phi)*(Tc-Ts)/Rs)
   }
   Ta <- uniroot(error,c(100,10000),tol=1e-6)$root
   return(Ta)

}


find_clothed_Ta <- function(Rf,pa=pa0) {
   pa0 <- common_variables$pa0
   if (Rf < 0) {
      stop('The Rf provided as input is outside of bounds')
   }
   Ta_lower <- find_clothed_Rfinfty_Ta(pa) + 1e-4
   Ta_upper <- find_clothed_Rfzero_Ta(pa) - 1e-4
   if (Rf > find_clothed_state(Ta_lower,min(1,pa/pvstar(Ta_lower)))$Rf) {
      return(Ta_lower)
   }
   if (Rf < find_clothed_state(Ta_upper,min(1,pa/pvstar(Ta_upper)))$Rf) {
      return(Ta_upper)
   }
   error <- function(Ta) {
      rh <- min(pa/pvstar(Ta),1)
      return(find_clothed_state(Ta,rh)$Rf - Rf)
   }
   Ta <- uniroot(error,c(Ta_lower,Ta_upper),tol=1e-6)$root
   return(Ta)
}


find_naked_Rszero_Ta <- function(pa=pa0) {

   # Get parameters
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa_requested <- pa

   error <- function(Ta) {
      pa <- min(pa_requested,pvstar(Ta))
      Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(p*rgasv)*(pc-pa))
      Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)
      return(Q-Qv-(Tc-Ta)/Ra-(pc-pa)/Za)
   }
   Ta <- uniroot(error,c(300,400),tol=1e-12)$root

   return(Ta)

}


find_naked_Ta <- function(Rs,pa=pa0) {
   pa0 <- common_variables$pa0
   Ta_lower <- find_clothed_Rfzero_Ta(pa) #+ 1e-4
   Ta_upper <- find_naked_Rszero_Ta(pa) #- 1e-4
   error <- function(Ta) {
      rh <- min(pa/pvstar(Ta),1)
      return(find_naked_state(Ta,rh)$Rs - Rs)
   }
   if (sign(error(Ta_lower)) != sign(error(Ta_upper))) {
      Ta <- uniroot(error,c(Ta_lower,Ta_upper),tol=1e-6)$root
   } else if (sign(error(Ta_upper)) != sign(error(9e3))) {
      Ta <- uniroot(error,c(Ta_upper,9e3),tol=1e-6)$root
   } else {
      Ta <- uniroot(error,c(0,Ta_lower),tol=1e-6)$root
   }
   return(Ta)
}


find_warming_Ta <- function(dTcdt,pa=pa0) {
   pa0 <- common_variables$pa0
   error <- function(Ta) {
      rh <- min(pa/pvstar(Ta),1)
      return(find_warming_state(Ta,rh)$dTcdt - dTcdt)
   }
   Ta <- uniroot(error,c(100,1e4),tol=1e-6)$root
   return(Ta)
}


# Return the Ta at rh that gives the equilibrium Tc
warming_equilibrium_Ta <- function(Tc,rh) {

   # Get parameters
   Tc_save <- Tc
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])
   Tc <- Tc_save

   pc <- phisalt*pvstar(Tc)

   error <- function(Ta) {
      pa <- rh*pvstar(Ta)
      Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(p*rgasv)*(pc-pa))
      Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)
      return((Q-Qv-(Tc-Ta)/Ra-(pc-pa)/Za)/Cc)
   }
   Ta <- uniroot(error,c(100,10000),tol=1e-6)$root

   return(Ta)

}


# Return the equilibrium Tc at Ta and rh
warming_equilibrium_Tc <- Vectorize(function(Ta,rh) {

   # Get parameters
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa <- rh*pvstar(Ta)

   # Solve for Tc
   error <- function(Tc) {
      pc <- phisalt*pvstar(Tc)
      Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(p*rgasv)*(pc-pa))
      Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)
      return((Q-Qv-(Tc-Ta)/Ra-(pc-pa)/Za)/Cc)
   }
   Tc <- uniroot(error,c(100,10000),tol=1e-6)$root

   return(Tc)

})
   

# Return time to death
time_to_death <- function(Ta,rh) {

   # Get parameters
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   if (warming_equilibrium_Tc(Ta,rh) <= Tcfatal) { return(Inf) }

   pa <- rh*pvstar(Ta)

   dt <- 1
   t <- 0
   Tc <- 310
   while (Tc < Tcfatal) {

      pc <- phisalt*pvstar(Tc)
      Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(p*rgasv)*(pc-pa))
      Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)

      dTcdt <- (Q-Qv-(Tc-Ta)/Ra-(pc-pa)/Za)/Cc
      Tc <- Tc + dTcdt*dt
      t <- t + dt

   }

   return(t)

}

is_covering <- function(Ta,rha) {
   if (rha < 0 | rha > 1) {
      stop('The rha provided as input is outside of bounds')
   }
   pa <- rha*pvstar(Ta)
   return(Ta < find_clothed_Rfinfty_Ta(pa))
}


is_clothed <- function(Ta,rha) {
   if (rha < 0 | rha > 1) {
      stop('The rha provided as input is outside of bounds')
   }
   pa <- rha*pvstar(Ta)
   return( Ta >= find_clothed_Rfinfty_Ta(pa) &
           Ta <= find_clothed_Rfzero_Ta (pa) )
}


is_naked <- function(Ta,rha) {

   return(!is_warming(Ta,rha) & !is_covering(Ta,rha) & !is_clothed(Ta,rha))

}


is_warming <- function(Ta,rha) {

   # Get parameters
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   pa <- rha*pvstar(Ta)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)

   return(Q-Qv-(Tc-Ta)/Ra-(pc-pa)/Za > 0)

}


condition <- Vectorize(function(Ta,rha,verbose=FALSE) {

   # Get parameters
   params <- c(common_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   if (is_warming(Ta=Ta,rha=rha)) {
      if (verbose) { cat('Region VI, warming up (dTc/dt > 0)\n') }
      return(6)
   } else if (is_naked(Ta=Ta,rha=rha)) {
      state <- find_naked_state(Ta=Ta,rha=rha)
      if (state$ps < phisalt*pvstar(state$Ts)) {
         if (verbose) { cat('Region IV, naked (variable Rs, ps < phisalt*pvstar)\n') }
         return(4)
      } else {
         if (verbose) { cat('Region V, naked dripping sweat (variable Rs, ps = phisalt*pvstar)\n') }
         return(5)
      }
   } else if (is_clothed(Ta=Ta,rha=rha)) {
      state <- find_clothed_state(Ta=Ta,rha=rha)
      Ta <- find_clothed_Ta(Rf=state$Rf,pa=pa0)
      if (pvstar(Ta) < pa0) {
         if (verbose) { cat('Region II, clothed (variable Rf, pa = pvstar)\n') }
         return(2)
      } else {
         if (verbose) { cat('Region III, clothed (variable Rf, pa = pref)\n') }
         return(3)
      }
   } else if (is_covering(Ta=Ta,rha=rha)) {
      if (verbose) { cat('Region I, covering up (variable phi)\n') }
      return(1)
   } else {
      stop('Region not found')
   }
})


region_bounds <- function() {

   # Get parameters
   for (varname in c('Tcfatal','phi0','pa0')) {
      assign(varname,common_variables[[varname]])
   }

   # Warm edge of region I
   hi <- find_covering_Ta(phi0,pa0)

   # Warm edge of region III
   hi <- c(hi,find_clothed_Ta(1e-9,pa0))
   
   # Warm edge of region V
   hi <- c(hi,find_warming_Ta(0,pa0))

   # Within region VI, boundary between sick and dying
   error <- function(Ta) {
      return(Tcfatal-warming_equilibrium_Tc(Ta,pa0/pvstar(Ta)))
   }
   hi <- c(hi,uniroot(error,c(350,375),tol=1e-12)$root)
   
   return(hi)

}


validate_covering <- function(state) {

   # Get parameters
   params <- c(common_variables,covering_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   # Get state
   for (varname in c('Ta','pa','Ts','ps','phi')) {
      assign(varname,state[[varname]])
   }

   # Calculate equation residuals
   Ra <- 1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta)+hc)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   cat(
      'Covering 1: ',Q-Qv-(1-phi)*(Tc-Ts)/Rs,'\n',
      'Covering 2: ',(Tc-Ts)/Rs-(Ts-Ta)/Ra-(ps-pa)/Za,'\n',
      'Covering 3: ',(pc-ps)/Zs-(ps-pa)/Za,'\n',sep='')

}


validate_clothed <- function(state) {

   # Get parameters
   params <- c(common_variables,clothed_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   # Get state
   for (varname in c('ps','psbar','pf','pa',
                     'Ts','Tsbar','Tf','Ta',
                     'Rf')) {
      assign(varname,state[[varname]])
   }
   Zf <- r*Rf

   Ra    <- 1/(phirad   *epsilon   *sigma*(Ts^2+Ta^2)*(Ts+Ta) + hc   )
   Rabar <- 1/(phiradbar*epsilonbar*sigma*(Tf^2+Ta^2)*(Tf+Ta) + hcbar)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   cat(
      'Clothed 1: ',Q-Qv-phi*(Tc-Tsbar)/Rs-(1-phi)*(Tc-Ts)/Rs,'\n',
      'Clothed 2: ',(Tc-Ts)/Rs-(Ts-Ta)/Ra-(ps-pa)/Za,'\n',
      'Clothed 3: ',(Tc-Tsbar)/Rs-(Tsbar-Tf)/Rf-(psbar-pf)/Zf,'\n',
      'Clothed 4: ',(Tsbar-Tf)/Rf-(Tf-Ta)/Rabar,'\n',
      'Clothed 5: ',(pc-ps)/Zs-(ps-pa)/Za,'\n',
      'Clothed 6: ',(pc-psbar)/Zs-(psbar-pf)/Zf,'\n',
      'Clothed 7: ',(psbar-pf)/Zf-(pf-pa)/Zabar,'\n',sep='')

}


validate_naked <- function(state) {

   # Get parameters
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   # Get state
   for (varname in c('Ta','pa','Ts','ps','Rs')) {
      assign(varname,state[[varname]])
   }

   # Calculate equation residuals
   Zs <- 6.0e8*Rs^5
   Ra <- 1/(phirad*epsilon*sigma*(Ts^2+Ta^2)*(Ts+Ta)+hc)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   cat(
      'Naked 1 : ',Q-Qv-(Tc-Ts)/Rs,'\n',
      'Naked 2 : ',(Tc-Ts)/Rs-(Ts-Ta)/Ra-(ps-pa)/Za,'\n',sep='')
   if (ps < phisalt*pvstar(Ts)) {
      cat(
      'Naked 3a: ',(pc-ps)/Zs-(ps-pa)/Za,'\n',sep='')
   } else {
      cat(
      'Naked 3b: ',ps - phisalt*pvstar(Ts),'\n',sep='')
   }

}


validate_warming <- function(state) {

   # Get parameters
   params <- c(common_variables,warming_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   # Get state
   for (varname in c('Ta','pa','dTcdt')) {
      assign(varname,state[[varname]])
   }

   # Calculate equation residuals
   Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)
   Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(rgasv*p)*(pc-pa))
   cat(
      'Warming 1: ',Cc*dTcdt-Q+Qv+(Tc-Ta)/Ra+(pc-pa)/Za,'\n',sep='')

}


heatindex <- Vectorize(function(T,rh,verbose=FALSE,debug=FALSE) {

   if (rh < 0 | rh > 1) {
      stop('The rh provided as input is outside of bounds')
   }

   # Get parameters
   params <- c(common_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   region <- condition(Ta=T,rha=rh)
   if (region==6) {
      state <- find_warming_state(Ta=T,rha=rh)
      if (debug) { validate_warming(state) }
      if (verbose) {
         cat('Region VI, warming up (dTc/dt > 0)\n')
         cat('dTc/dt =',state$dTcdt*3600,'K/hour\n')
      }
      Ta <- find_warming_Ta(dTcdt=state$dTcdt,pa=pa0)
      if (debug) { cat('dTcdt error:',find_warming_state(Ta,rh=min(1,pa0/pvstar(Ta)))$dTcdt-state$dTcdt,'\n') }
   } else if (region %in% 4:5) {
      state <- find_naked_state(Ta=T,rha=rh)
      if (debug) { validate_naked(state) }
      if (state$ps < phisalt*pvstar(state$Ts)) {
         if (verbose) { cat('Region IV, naked (variable k, ps < phisalt*pvstar)\n') }
      } else {
         if (verbose) { cat('Region V, naked dripping sweat (variable k, ps = pvstar)\n') }
      }
      if (verbose) {
         bloodflow <- (1/state$Rs - 5.28) / (1e3*4184/A)   # m^3/s
         cat('Blood flow is',round(bloodflow*1e3*60,dig=2),'l/min\n')
      }
      Ta <- find_naked_Ta(Rs=state$Rs,pa=pa0)
      if (debug) { cat('Rs error:',find_naked_state(Ta,rh=min(1,pa0/pvstar(Ta)))$Rs-state$Rs,'\n') }
   } else if (region %in% 2:3) {
      state <- find_clothed_state(Ta=T,rha=rh)
      if (debug) { validate_clothed(state) }
      Ta <- find_clothed_Ta(Rf=state$Rf,pa=pa0)
      if (pvstar(Ta) < pa0) {
         if (verbose) { cat('Region II, clothed (variable kbar, pa = pvstar)\n') }
      } else {
         if (verbose) { cat('Region III, clothed (variable kbar, pa = pref)\n') }
      }
      if (debug) { cat('Rf error:',find_clothed_state(Ta,rh=min(1,pa0/pvstar(Ta)))$Rf-state$Rf,'\n') }
   } else if (region == 1) {
      state <- find_covering_state(Ta=T,rha=rh)
      if (debug) { validate_covering(state) }
      if (verbose) { cat('Region I, covering (variable phi)\n') }
      Ta <- find_covering_Ta(phi=state$phi,pa=pa0)
      if (debug) { cat('phi error:',find_covering_state(Ta,rh=min(1,pa0/pvstar(Ta)))$phi-state$phi,'\n') }
   }
   
   return(Ta)

})


heatindex_of_equilibrium_Tc <- Vectorize(function(Tc) {
   if (Tc < 310) {
      return(NA)
   } else {
      return(heatindex(warming_equilibrium_Ta(Tc,0),0))
   }
})
   

# Return time to death
time_to_death <- function(Ta,rha) {

   if (rha < 0 | rha > 1) {
      stop('The rha provided as input is outside of bounds')
   }

   # Get parameters
   params <- c(common_variables,naked_variables)
   for (i in 1:length(params)) assign(names(params)[i], params[[i]])

   if (warming_equilibrium_Tc(Ta,rha) <= Tcfatal) { return(Inf) }

   pa <- rha*pvstar(Ta)

   dt <- 1
   t <- 0
   Tc <- 310
   while (Tc < Tcfatal) {

      pc <- phisalt*pvstar(Tc)
      Qv <- eta*Q*(cpa*(Tc-Ta) + L*rgasa/(p*rgasv)*(pc-pa))
      Ra <- 1/(phirad*epsilon*sigma*(Tc^2+Ta^2)*(Tc+Ta)+hc)

      dTcdt <- (Q-Qv-(Tc-Ta)/Ra-(pc-pa)/Za)/Cc
      Tc <- Tc + dTcdt*dt
      t <- t + dt

   }

   return(t)

}


# Return the Ta corresponding to hi and rh
Ta_of_hi <- function(hi,rh) {
   error <- function(Ta) {
      return(hi-heatindex(Ta,rh))
   }
   return(uniroot(error,c(0,9e3))$root)
}


# Bloodflow (m^3/s)
bloodflow <- Vectorize(function(hi) {
   # bounds <- region_bounds()
   # Hard-coded bounds to speed up the calculation
   bounds <- c(238.48963454513125271,298.43929324762166289,
      344.65446922486057701,366.44386593525092621)
   rh <- min(1,with(common_variables,pa0)/pvstar(hi))
   if (hi < bounds[1]) {
      Rs <- find_covering_state(hi,rh)$Rs
   } else if (hi < bounds[2]) {
      Rs <- find_clothed_state(hi,rh)$Rs
   } else if (hi >= bounds[3]) {
      return(Inf)
   } else {
      Rs <- find_naked_state(hi,rh)$Rs
   }   
   bloodflow <- (1/Rs - 5.28) / (1e3*4184/with(common_variables,A))   # m^3/s
   return(bloodflow)
})
