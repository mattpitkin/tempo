c      $Id$
	subroutine bnrydds(torb,fctn)

c  Damour et Deruelle modele pour le chronometrage des temps d'arrive 
c  au barycentre du systeme solaire, a la premiere approximation
c  post-newtonienne de la Relativite Generale.

c  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch.
c  Pulsar proper time is then TP=T+TORB.
c  Units are such that c=G=1. Thus masses have units of seconds, with
c  one solar mass = 4.925490947 usec.

c  Also computes the binary orbit-related values of fctn: partial
c  derivatives of each arrival time residual with respect to the model
c  parameters.

c  Initial guesses for all parameters must be placed in common/orbit/ by the
c  calling program. 

c  Changes by Wex for Double Pulsar timing
c  - numerical inversion of the timing model
c  - retardation effect in Shapiro delay - Kopeikin & Schäfer 1999
c  - correction of Shapiro due to lensing - Klioner & Zschocke 2010
c  - bending/lens-rotational delay with retardation - Doroshenko & Kopeikin 1995, Rafikov & Lai 2006
c    -> latitudinal contribution (Rafikov & Lai 2006) not included. Depends on profile
c    details and is highly covariant with other parameters (see Kramer et al. 2021)
c  Note: approximations used are (still) more than sufficient, which has
c        been checked in analytical calculations and numerical simulations

c Last modified: 2021-Apr-08

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'acom.h'
	real*8 fctn(NPAP1),k,m2,theta
	parameter (twopi=6.28318530717958648d0,SUNMASS=4.925490947d-6)
	parameter (RAD=360.d0/twopi)
	include 'dp.h'
	include 'orbit.h'
      common /CRDS/ RCB(6),RBE(6),RCE(6),RSE(3),RCA(3),RSA(3),RBA(3),
     +             REA(3),PSID,EPSD,PC,PS,TSID,PRN(3,3),ATUT,UT1UT,ETAT,
     +             DTGR,TDIS,BCLT      
	include 'trnsfr.h'

	real*8 frb, tt0, tt, orbits

	frb = 1.d0/pb(1)
	an  = twopi/pb(1)
	k   = omdot/(RAD*365.25d0*86400.d0*an)
	m2  = am2*SUNMASS

	tt0 = (ct - t0(1))*86400.d0

	xp  = a1(1) + xdot*tt0
	ecc = e(1)  + edot*tt0
	er  = ecc*(1.d0 + dr)
	eth = ecc*(1.d0 + dth)

	sqr1me2 = DSQRT(1.d0 - ecc**2)
		
c  --> Inversion of timing model by iteration: begin of loop
	epsNum = 1.0d-12 ! numerical precision of inversion in seconds
	delta  = 0.0d0
 10	delta_old = delta        
	tt  = tt0 - delta
	orbits  = tt*frb - 0.5d0*(pbdot+xpbdot)*(tt*frb)**2
	norbits = orbits
	if(orbits.lt.0.d0) norbits = norbits - 1
	phase = TWOPI*(orbits - norbits)

c  Compute eccentric anomaly u by iterating Kepler's equation.
	u  = phase + ecc*DSIN(phase)*(1.d0 + ecc*DCOS(phase))
 100	su = DSIN(u)
	cu = DCOS(u)
	du = (phase - (u - ecc*su))/(1.d0 - ecc*cu)
	u  = u + du
	if(DABS(du).gt.1.d-14) goto 100

	onemecu = 1.d0 - ecc*cu
	cume    = cu - ecc

	cae  = (cu-ecc)/onemecu
	sae  = sqr1me2*su/onemecu
	ae1  = DATAN2(sae,cae)
	if(ae1.lt.0.d0) ae1 = ae1 + TWOPI  ! ae1 in [0,TWOPI)
	ae   = ae1 + norbits*TWOPI

	if (useannorb) then ! IHS 20210921 Note hardwiring of coti sign -- use of this option is not recommended unless you are prepared to recompile tempo to suit your system.
	  deltai0 = -rca(1)*dsin(pra) + rca(2)*dcos(pra)
          deltaj0 = (-rca(1))*dsin(pdec)*dcos(pra) +
     +              (-rca(2))*dsin(pdec)*dsin(pra) +
     +              ( rca(3))*dcos(pdec)
          if (usefixeddist) then  ! use fixed distance in kpc
            dist = 499.00478364D0/(twopi/fixeddist/1000/3600/360)
          else                    ! use parallax distance
            dist = 499.00478364D0/(twopi*px/1000/3600/360)
          endif
          Omkopeikin = twopi/4-(PAAscNode/360*twopi)
          omegax = -(1.d0/si)/dist*
     +           (deltai0*dcos(Omkopeikin)+deltaj0*dsin(Omkopeikin))
          coti = -dsqrt(1/si**2-1)  ! hack -- hard-wired for negative cot i
          xpterm = coti/dist *
     +           (deltai0*dsin(Omkopeikin)-deltaj0*dcos(Omkopeikin))
	  xp = xp*(1.d0+xpterm)
          omega=omz(1)/RAD + k*ae + omegax
        else
          omega=omz(1)/RAD + k*ae 
        endif
        sw=dsin(omega)
        cw=dcos(omega)

	omega = omz(1)/RAD + k*ae
	sw    = DSIN(omega)
	cw    = DCOS(omega)
	
	psi  = omega + ae1 ! angle w.r.t. ascending node
	spsi = DSIN(psi)
	cpsi = DCOS(psi)

c *** Roemer delay (DD)
	alpha = xp*sw
	beta  = xp*DSQRT(1.d0 - eth**2)*cw
	dRoe  = alpha*(cu - er) + beta*su

c *** Einstein delay (DD)
	dEin  = gamma * su

c *** Shapiro delay 
	sidds = 1.d0 - DEXP(-1.d0*shapmax) ! see Kramer et al. 2006, AnP
	brace = onemecu - sidds*(sw*cume + sqr1me2*cw*su)

	if(nshapho.eq.0)then

		dlogbr = DLOG(brace)
		dSha   = -2.d0*m2*dlogbr

	else ! NW: higher order corrections related to light propagation (for Double Pulsar only!)

 	   ratiompmc = 1.0714d0  ! mass ratio mp/mc (= xc/xp) for Double Pulsar 
       xc = xp * ratiompmc
       xR = xp + xc;         ! aR*sini/c [s]
       aR = xR/sidds         ! aR [s] 

c   --> Account for lensing contribution to propagation time 
c       - simplified version (cf. Zschocke & Klioner 2010, eq. (73))
       
	   epsLen = 2.d0*m2/aR

c  --> 1.5pN contribution to Shapiro, i.e. leading order velocity dependence 
c      (Kopeikin & Schäfer 1999, eq. (130))

	   epsVel = an*xp/sidds*ratiompmc*ecc*su 
     :          - an*xp*sidds*ratiompmc/sqr1me2
     :            * (sw*cume + sqr1me2*cw*su)
     :            * (ecc*cw + (cw*cume - sqr1me2*sw*su)/onemecu)
       
c  --> Shapiro delay incl. higher order corrections
	   braceho =  brace + (epsLen + epsVel) * shaphof
	   dlogbr  =  DLOG(braceho)
	   dShaho  = -2.d0*m2*dlogbr

c  --> Bending/lens-rotational delay 
c      - Doroshenko & Kopeikin 1995 approximation,  
c        including retardation to leading order (shift in c's position)
c      - assumes pulsar spin to be perpendicular to the orbital plane,
c        i.e. nearly perpendicular to the line-of-sight
       dfdt     = an*sqr1me2/onemecu**2
       dpsi     = dfdt*ratiompmc*dRoe
       cpsiRet  = cpsi - spsi*dpsi
       braceRet = brace + epsVel * shaphof
       dRotDefl = 2.d0*m2/(TWOPI*f0)/xR * cpsiRet/braceRet

c   --> Sum of all the signal propagation contributions
      
       dSha = dShaho + dRotDefl * shaphof
        
	endif

c *** Aberration (DD)
	dAbe = a0*(spsi + ecc*sw) + b0*(cpsi + ecc*cw)

c *** Sum of delays
	delta = dRoe + dEin + dSha + dAbe

	diff = DABS(delta - delta_old)
	if(diff.gt.epsNum) goto 10
c  --> Inversion of timing model by iteration: end of loop              

	torb = -delta

c  Now we need the partial derivatives. Use DD equations 62a - 62k,
c  plus DDS specific modifications
	csigma = xp*(-sw*su+sqr1me2*cw*cu)/onemecu
	ce     = su*csigma-xp*sw-ecc*xp*cw*su/sqr1me2
	cx     = sw*cume+sqr1me2*cw*su
	comega = xp*(cw*cume-sqr1me2*sw*su)
	cgamma = su
	cdth   = -ecc*ecc*xp*cw*su/sqr1me2
	cm2    = -2.d0*dlogbr
	if(nshapho.eq.0)then
		csi = 2.d0*m2*(sw*cume+sqr1me2*cw*su)/brace
	else
		csi = 2.d0*m2*(sw*cume+sqr1me2*cw*su)/braceho
	endif
	cshapmax = csi * DEXP(-1.d0*shapmax)
	if(nshapho.eq.0)then
		cshaphof = 0.d0
	else
		cshaphof = -2.d0*m2/braceho*(epsLen + epsVel) + dRotDefl
	endif

	fctn(9)=cx*f0
	fctn(10)=ce*f0
	fctn(11)=-csigma*an*f0*86400.d0
	fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
	fctn(13)=comega*f0
	fctn(14)=ae*fctn(13)/(an*RAD*365.25d0*86400.d0)
	fctn(15)=cgamma*f0
	fctn(18)=0.5d-6*tt0*fctn(12)
	fctn(20)=cshapmax*f0
 1	fctn(39)=cshaphof*f0
	fctn(22)=cm2*f0*SUNMASS
	fctn(23)=cdth*f0
	fctn(24)=tt0*fctn(9)
	fctn(25)=tt0*fctn(10)

	return
	end


