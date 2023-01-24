        PROGRAM TestBenford
C        
C "Blind" program for computing the
C * exact size and power *  of several (GOF) tests of Benford's law.
C N.B.: To be optimized for speed and efficiency
C USE IMSL ROUTINES FOR RANDOM NUMBER GENERATION + OTHER TASKS
C 
C Date: January 2023
C     
C ** FASTER ALGORITHM FOR EXACT P-VALUE COMPUTATION OF COMBINED TESTS **
C  * maxn = 10**4 * TO REDUCE MEMORY STORAGE IN SIMULATIONS ==> THIS
C    COULD BE CHANGED (maxn > 10**4) FOR EMPIRICAL ANALISIS OF BIG DATA
C    (idistr < 0) ==> set *maxsim2=1*
C                     set *maxsim1<10**6*
C  * COMPUTATION OF P-VALUES IS IN *SINGLE PRECISION* TO REDUCE MEMORY
C    STORAGE (sorted statistics and p-values are currently saved with 
C    more digits than required by single precision)      
C
C Introduce new dummy *itest*  ==>  if *itest=1* only
C     *individual tests* are performed in power computations:
C     skip the loop for computing exact p-values of combined
C     tests (but for simplicity retain subsequent computations)
C
C List of tests currently implemented (version - April 2022):
C     different arrays *powind()* and *powcomb()* are used  
C     for *individual* and *combined* tests 
C     correspondingly array *perc()* is replaced by arrays
C     *percind()* and *perccomb()*
C     ==> less elegant but more flexible code when new tests are introduced
C Individual tests:
C powind(1): X2      = 1st digit Chi-square test
C powind(2): X2_2d   = 2-digit Chi-square test
C powind(3): Q       = Hotelling test based on sum invariance
C powind(4): KS      = Kolmogorov-Smirnov test (6-dimensional vector output of DKSONE)
C powind(5): X2DIF   = Q - X2 (one tail)
C powind(6): X2DIF2T = |Q-X2| (two tails) 
C powind(7): SCORE1  = Score test based on 1 trigonometric component
C powind(8): SCORE2  = Score test based on 2 trigonometric components
C powind(9): SCOREBIC= Score test based on data-driven (BIC) choice of the
C                      number of trigonometric components
C powind(10): XLRT1  = LR test based on 1 trigonometric component
C powind(11): XLRT2  = LR test based on 2 trigonometric components
C powind(12): XLRTBIC= LR test based on data-driven (BIC) choice of the
C                      number of trigonometric components
C powind(13): VK     = Kuiper test (uses output from DKSONE)
C powind(14): ZZZZZ - TO DO
C powind(15): ZZZZZ - TO DO
C
C Combined tests: 
C powcomb(1): X2Q = Combined X2-Q
C powcomb(2): QKS = Combined Q-KS
C powcomb(3): X2QKS = Combined X2-Q-KS
C powcomb(4):X2DIFKS = Combined X2DIF-Q-KS (new 21/05/2021 - previous X2DIFKS: 
C                                          Combined X2DIF-KS ==> N.B.: variable 
C                                          and array names are *not* modified)
C powcomb(5):X2DIF2TKS  = Combined X2DIF2T-Q-KS - new 09/06/2021 
C powcomb(6):X2DIFQSCORE= Combined X2DIF-Q-SCOREBIC (new February 2022)
C powcomb(7):X2DIFQXLRT = Combined X2DIF-Q-XLRTBIC (new February 2022)
C powcomb(8):X2DIFQVK   = Combined X2DIF-Q-VK (new February 2022)
C powcomb(9):X2DIFKSFR  = ZZZZZ - TO DO
C powcomb(10):X2DIFVKFR = ZZZZZ - TO DO
C
C * iapprox = 1 *: approximate algorithm for computation of power of combined
C                  tests ==> faster under close alternatives
C * icomb* = 1: check the exact percentile of p-values of combined tests ==> if
C   exact percentile < testsize(1) the approximate algoritm for the corresponding
C   combined test is expected to exact
C
C      
C MAX # SIMULATIONS FOR NULL DF AND PERCENTILE ESTIMATION: 
C * maxsim1 = 10**6 * (for memory constraints - use real*4 vector)
C     with some memory savings, maxsim1 up 1.5*10**6 
C
C Test statistics are either simulated or read from file according to
C the value of *iread* ==> iread = 1: read form file; otherwise simulate
C
C MAX # SIMULATIONS FOR POWER AND SIZE ESTIMATION: 
C * maxsim2 = 10**5 * 
C COMPARES EXACT POWER-SIZE WITH ASYMPTOTIC POWER-SIZE 
C (single precision is used)
C
C nsim1 = ACTUAL # SIMULATIONS FOR NULL DF AND PERCENTILE ESTIMATION 
C         (read from file)
C nsim2 = ACTUAL # SIMULATIONS FOR POWER AND SIZE ESTIMATION 
C         (read from file)
C rcont = CONTAMINATION RATE (read from file) ==> rcont = 1 for *power analysis*
C         (including check of size: idistr=0); 
C         0 < rcont < 1 for outlier analysis under Benford null distribution
C nqprop = 6: # OF TEST SIZES to be read from file
C npar   = 6: # OF PARAMETERS OF ALTERNATIVE DISTRIBUTIONS to be read from file
C             N.B.: NOT ALL THE PARAMETERS NEED TO BE ACTUALLY USED FOR  
C             POWER COMPUTATION (e.g. location and/or scale) ==> see
C             options for the specific value of *idistr*
C             N.B.: FOR SOME DISTRIBUTIONS npar MUST BE AN *INTEGER*
C ntestind=15:  # OF INDIVIDUAL TESTS (see list above) - *new April 2022*
C ntestcomb=10: # OF COMBINED TESTS (see list above) - *new April 2022*
C
C
C SEVERAL ALTERNATIVE DISTRIBUTIONS FOR POWER COMPARISON 
C (according to the value of idistr):
C idistr.lt.0 ==> read data from file for empirical analysis
C idistr < 0: Read data from file (npartrue=1 and set nsim2=1, iapprox=0)
C             (data are jittered to avoid ties in DKSONE)
C
C distinguish between *idistr.lt.(-1)* and *idistr.eq.(-1)*
C if idistr.lt.(-1) ==> keep *nsim2* as in file .in1 to allow
C                       analysis of multiple (simulated) 
C                       data sets each of size n: compute
c                       power + save statistics and p-values 
C                       (must check the compatibility of
C                        maxim2 and maxsim1)                     
C if idistr.eq.(-1) ==> set nsim2=1 as in previous versions
C     N.B.: although less efficient, ensure compatibility with 
C           previous code
C
C if idistr.lt.0 ==> introduce new dummy *ijit*: 
C                    if *ijit.gt.0* jittering (as before)
C                    otherwise no jittering in the data
C                    in empirical analysis                                                              
C
C introduce new dummy *idig1*  ==>  if *idig1=1* 1st digit
C     is Benford in the selected distribution (as in BenfGBL)
C     ==> General *Manipulated Benford distribution*
C
C idistr = 0: Benford (size estimation) - No parameters: npartrue=1
C idistr = 1: Normal - npartrue=npar: npar/2 different values of the scale  
C             parameter (standard deviation) * 2 different values of the
C             location parameter
C idistr = 2: Lognormal - npartrue=npar-1 different values of the shape  
C             parameter (scale read from file and fixed)
C idistr = 3: (standard) Weibull - npartrue=npar different values of the the shape  
C             parameter (scale is fixed and set = 1: default of routine DRNWIB)
C idistr = 4: (standard) Gamma - npartrue=npar different values of the the shape  
C             parameter (scale is fixed and set = 1: default of routine DRNGAM)
C idistr = 5: Beta - npartrue=npar: npar/2 different values of shape parameter alpha  
C             * 2 different values of shape parameter beta
C idistr = 6: Uniform - npartrue=npar different values of the range 
C             of interval (0,a)
C idistr = 7: Stable - npartrue=npar-1 different values of the shape  
C             parameter alpha (skewness beta read from file and fixed)
C     choice of the simulation algorithm for Stable distribution
C idistr = 8: Generalized Benford - npartrue=npar different values of the 
C             GBL parameter
C idistr = 9: BenfGBL = First-Digit Benford + Generalized Benford - as above
C idistr =10: 2-component Generalized Benford Mixture - npartrue=npar: npar/2 
C             different values of the GBL parameter (with opposite sign 
C             in the two components) * 2 different values of the mixing 
C             proportion (refers to component 1)
C idistr =11: Dirac distribution for the first digit (mainly useful for 
C             outlier analysis) - npartrue = npar: fixed value of the replacing digit
C idistr =12: Truncated Benford distribution (Last-Digit(s) Dirac distribution)
C             - npartrue = npar: ==> the Benford significand is truncated 
C             at the *7th* decimal place (8 digits) and the *last* *npar=6* digits
C             are successively replaced with a fixed digit value 
C                                ==> the resulting significand is jittered
C                                    in the 10th decimal place because of the 
C                                    continuity requirement of routine DKSONE
C more jittering in large samples (n>200)
C more jittering also when n=200
C
C idistr=120: NEW Truncated Benford distribution (Last-Digit(s) Dirac distribution)
C             - npartrue = 1: ==> the Benford significand is truncated at the *6th*
C              decimal place (7 digits) and the *last* *npar=6* digits are replaced
C              with the fixed digit value 0 in a *given proportion of cases*
C              (these given proportions are read as input parameters)
C                                ==> the resulting significand is jittered
C                                    in the 10th decimal place because of the 
C                                    continuity requirement of routine DKSONE
C
C even more jittering for big n and idistr.eq.120: two options 
C (very minor effect - see trader Trader 7966819015)
C
C idistr =13: Cardioid distribution - npartrue=npar: npar/2 different values of
C             the concentration parameter 0<a<1 * 2 different values of the
C             location parameter b (to be multiplied by 2*pi)
C
C     *25/06/21* The Cardioid distribution is introduced   
C                + other minor modifications (e.g.: npartrue/2 ==> npar2)
C
C idistr =14: NNTS (Non Negative Trigonometric Sums: Fernandez-Duran, Biometrics, 2004)
C             npartrue=1 ==> read: *numNNTS = actual number of trigonometric sums*
C                                             must be <= npar - 1 = 5 
C                                             (npar = max numNNTS + intercept) 
C                            read: *parNNTS = complex parameter vector*
C                                             of max dimension npar, 
C                                             including intercept (always real)  
C
C     *15/07/2021* ALSO SAVE DATA WITH SIGNIFICANDS IN THE SAME FILE
C                  (for value distributions: see before instruction 1001)
C 
C     *16/07/21* The NNTS distribution is introduced  
C
C     *17/07/21* The 2-component Normal Mixture distribution is introduced  
C
C OLD VERSION idistr =15: 2-component Normal Mixture - npartrue=(2/3)*npar: npar/3 pairs 
C             of location parameters * 2 different values of the mixing 
C             proportion (refers to component 1)
C             scale is *fixed* and equal in both components ==> read from file
C *SEE PROGRAM VERSION OF 08/05/2022*
C
C NEW VERSION *OCTOBER 2022*
C idistr =15: 2-component heteroschedastic Normal Mixture - npartrue=2: 
c             2 location parameters (1 for each mixture component)
C             2 scale parameters (1 for each mixture component)
C             2 mixing proportions (refer to component 1)
C             ==> npartrue = (1 heteroschedastic mixture) * 2 mixing proportions
C
C 
C idistr =16: 3-component Normal Mixture - npartrue=1: 3 fixed location parameters 
C                                                      3 fixed scale parameters
C                                                      3 fixed mixing proportions (1/3)
C
C     *08/12/21* The 3-component Normal Mixture distribution is introduced 
C 
C     *February 2022* New combined tests using X2DIF are introduced
C
C     *22/03/2022* New option for saving test statistics + minor related changes
C                  * always save if idistr < 0 (empirical analysis)
C                  * only save when *isavetest=1* if idistr >= 0 (simulated data)
C     OLD OPTION: *isavetest = isavesig*
C     NEW OPTION: *isavetest* read from file
C
C     *April 2022* New (individual and combined) tests based on the significand
C                  fractional part
C
C     POSSIBLY NEED TO COMPUTE EXACT P-VALUES IN DOUBLE PRECISION 
C     TO AVOID INACCURACIES WHEN THE P-VALUE IS CLOSE TO 1 - TO DO
C
C MAX SAMPLE SIZE: * maxn = 10**4 * IN SiMULATIONS 
C n = ACTUAL SAMPLE SIZE (read from file)
C
c
c        use CHIDF_INT
c        use imsl_libraries
c     
c      USE CHIIN_INT
c
      IMPLICIT NONE
      real*8 dchiin,dchidf,drnunf,drnnof,tol,xdif
      real*8 pi
      integer*4 nca,nra,lda,ldginv
      integer*4 irank
      integer*4 nsim1,nsim2,n,idistr,iread,ncont,ngood,iapprox,idig1,
     + itest
      real*8 rcont
      integer*4 i,j,jj,jjj,jjjj,jsimul
      integer*4 isavesig,isavetest,jout,ijit
      integer*4 ij
      integer*4 maxsim1,maxsim2,maxn,nqprop,npar,ntest,ntestind,
     + ntestcomb
      integer*4 numNNTS
      integer*4 iseed,ioptMT,iseed0,nmiss
      integer*4 loopn,loopp
      integer*4 ncomps,ncompx,ncompmax
      real*8 xncompmax
c      
c FOR SIMULATION AND EMPIRICAL ANALYSIS OF MODERATE DATA SETS (standard simulation setting)
      parameter (maxsim1=10**6)
c      parameter (maxsim1=1.5*10**6)
      parameter (maxsim2=10**5)
      parameter (maxn=10**4)     
c FOR EMPIRICAL ANALYSIS OF BIG DATA SETS (idistr < 0)
!      parameter (maxsim1=10**3)
!      parameter (maxsim2=1)
!      parameter (maxn=10**6)
c FOR SIMULATION AND EMPIRICAL ANALYSIS OF VERY BIG (BENFORD) DATA SETS 
!      parameter (maxsim1=100)
!      parameter (maxsim2=10**5)
!      parameter (maxn=10**6)
c FOR SIMULATION OF TEST STATISTICS UNDER THE ALTERNATIVE (isavetest = 1)
!      parameter (maxsim1=10**3)
!      parameter (maxsim2=10**6)
!      parameter (maxn=10**4)
c
      parameter (ioptMT=9)
      parameter (nqprop=6)
      parameter (npar=6)
      parameter (ntestind=15)
      parameter (ntestcomb=10) 
      parameter (ntest=ntestind+ntestcomb)     
      parameter (loopn=10000)
      parameter (loopp=1000)
      parameter (tol=1.0d-10)
      parameter (pi=3.14159265358979d0)
c
      real*4 testsize(nqprop)
      real*4 qprop(nqprop),xlo(nqprop),xhi(nqprop)
      real*8 rn,rnb,rng,ssb,ssg,rn1,rn2,xx1,xx2,ss,ss1,ss2,fn
      integer*4 iss,mm
      integer*8 iss8
c defines vectors of parameters - npar = max # of parameters, including 
c location/scale; npartrue = actual # of parameters used in power computation
      real*8 par(npar),partrue,xloc,scal,tolstab(npar),shift
      real*8 scal1,scal2
      integer*4 istab
      integer*4 npartrue,npar2
c parameters and other quantities of NNTS distribution      
      complex*16 parNNTS(npar)
      complex*16 expNNTS,sumNNTS
      real*8 absparNNTS(npar),abssumNNTS,sumabsparNNTS,eNNTS,prodNNTS   
c defines vectors of simulated data, digits and significand
      real*8 simdata(maxn),s(maxn),sFR(maxn)
      integer*4 d1(maxn),d2(maxn),d12(maxn,2),d1b,d1g,d1dif,
     + ipartrue,ipval,ineg
c defines Benford probabilities, frequencies, significand and moments
      real*8 benf1p(9),benf12p(9,10),
     + benf1n(9),benf12n(9,10),benfsig,
     + ee,benfE,benfVAR(9,9),benfVARinv(9,9)
      real*8 dd,dd2
c defines test statistics (in single precision for quantile computation)
      real*8 X2,X2_2d,Q,KS(6),KS1,X2DIF,X2DIF2T,checkX2DIF,checkX2DIF2T,
     + XLRT1,XLRT2,XLRTBIC,SCORE1,SCORE2,SCOREBIC,VK,
     + KSFR(6),VKFR
      real*4 X2s,X2_2ds,Qs,KSs(6),X2DIFs,X2DIF2Ts,
     + XLRT1s,XLRT2s,XLRTBICs,SCORE1s,SCORE2s,SCOREBICs,VKs,
     + KSFRs(6),VKFRs
      real*4 X2ss,Qss,KSss,X2DIFss,X2DIF2Tss,XLRTBICss,SCOREBICss,VKss,
     + KSFRss,VKFRss
c not (yet) used
c      real*4 X2_2dss,XLRT1ss,XLRT2ss,SCORE1ss,SCORE2ss
c defines summary statistics
      real*8 corX2Q(2,2),X2Qmean(2),corX2DIF(2,2),X2DIFmean(2),
     + corX2DIF2T(2,2),X2DIF2Tmean(2),fn1,fn2
      real*8 avencompscoreb,avencompxlrtb,varncompscoreb,varncompxlrtb,
     + maxncompscoreb,maxncompxlrtb
      real*8 avencompscores(npar),avencompxlrts(npar),
     + varncompscores(npar),varncompxlrts(npar),
     + maxncompscores(npar),maxncompxlrts(npar)    
c defines (dummy) vector of simulated statistics for null DF + percentile
c estimation ==> simulated statistics are saved in files and then read 
c                in dummy vector *benfarray* for null percentile estimation; 
c                if iread=1 null percentiles are read from file
c                use *single precision* to allow a larger value of maxsim1
      real*4 benfarray(maxsim1)
      real*4 pval1(maxsim1),pval2(maxsim1),pval3(maxsim1)
      real*4 pvalmin,pvalX2,pvalQ,pvalKS,pvalX2DIF,pvalX2DIF2T
      real*4 pvalSCOREBIC,pvalXLRTBIC,pvalVK
      real*4 pvalKSFR,pvalVKFR
c not (yet) used
c      real*4 pvalX2_2d,pvalSCORE1,pvalSCORE2,pvalXLRT1,pvalXLRT2
      real*4 pvalX2Q,pvalQKS,pvalX2QKS,pvalX2DIFKS,pvalX2DIF2TKS
      real*4 pvalX2DIFQSCORE,pvalX2DIFQXLRT,pvalX2DIFQVK
      real*4 pvalX2DIFKSFR,pvalX2DIFVKFR
      real*4 pvalcom
      integer*4 iperm1,iperm2,iperm3
c use the same permutation vector *ipermstat* for all the test statistics 
c to save memory space
      integer*4 ipermstat(maxsim1)
      integer*4 istop,istopX2,istopQ,istopKS,istopX2DIF,istopX2DIF2T
      integer*4 istopSCOREBIC,istopXLRTBIC,istopVK
      integer*4 istopKSFR,istopVKFR
      integer*4 icomb(ntestcomb)
c defines arrays of dimension:
c   nqprop*ntestind  (individual tests) 
c   nqprop*ntestcomb (combined tests) 
c of estimated + asymptotic null 
c percentiles ==> use *single precision* for exact percentiles (the statistic 
c                 values are read from file in single precision)
c             ==> asymptotic percentiles are used for *individual tests* only
      real*4 percind(nqprop,ntestind),perccomb(nqprop,ntestcomb),
     + pp(nqprop)
      integer*4 iperc(nqprop),iperc2(nqprop)
      real*8 qprop2(nqprop),testsize2(nqprop),
     + df,qq,qqq,percas(nqprop,ntestind)
      real*8 pvalVKas,pvalVKFRas
c defines dummy array of null percentiles ==> to be used (in double precision)
c in power computation: arrays for individual + combined tests
      real*8 pind(nqprop,ntestind)
c defines vectors of power-size estimates (exact + asymptotic)
c N.B.: use array pow of dimension = npar*nqprop*ntest for all distributions
c       ==> the third dimension of *pow* corresponds to test (see list above)
      real*8 powind(npar,nqprop,ntestind),
     + powcomb(npar,nqprop,ntestcomb),
     + powas(npar,nqprop,ntestind)
      real*8 betaprime,tprob,tailprob,betap(npar)
c defines vectors for checking parameter/digit values      
      real*8 checkstab(npar),checkstabt(npar)
      real*8 checkmean(npar),checkvar(npar),checkmeant(npar),
     + checkvart(npar)
      real*8 dgamma,gg1,gg2
      real*8 checksig(9,npar),checkdig(9,npar),
     + checksigt(9,npar),checkdigt(9,npar)
      integer*4 ntrunc(npar),cumtrunc(npar)     
c      
      EXTERNAL BENFCDF
      EXTERNAL BENFCDFFR      
C
c Reads input parameters and check values
      open(unit=1,file='benf_in1.txt')
      read(1,*) nsim1,nsim2,n,idistr,iread,iapprox,rcont,isavesig,
     + isavetest,ijit,idig1,itest 
      close(1)
c      isavetest=isavesig
      if(nsim1.lt.1.or.nsim1.gt.maxsim1) then
        write(*,*) "Error in *nsim1*"
        go to 9999
      end if
      if(nsim2.lt.1.or.nsim2.gt.maxsim2) then 
        write(*,*) "Error in *nsim2*"
        go to 9999
      end if
      if(n.lt.1.or.n.gt.maxn) then 
        write(*,*) "Error in *n*"
        go to 9999
      end if
      if(rcont.lt.0.d0.or.rcont.gt.1.d0) then
        write(*,*) "Error in *rcont*"
        go to 9999
      end if        
      open(unit=2,file='benf_in2.txt')
      read(2,*) (testsize(j), j=1,nqprop)
      close(2)
      open(unit=2,file='benf_in2.txt')
      read(2,*) (testsize2(j), j=1,nqprop)
      close(2)
      fn1=dfloat(nsim1)
      fn2=dfloat(nsim2)
      fn=dfloat(n)
c Define maximum number of components for BIC selection according to 
c results of Inglot and Ledwina (1996) ==> the correct order of magnitude 
c is guaranteed; the constant is chosen in order to have 
c ncompmax=10 when n=100 and 8 <= ncompmax <= 13 when 50 <= n <= 200
      xncompmax=2.d0*(1.5d0+tol)+1.d0
      xncompmax=1.d0/xncompmax
      xncompmax=dfloat(n)**(xncompmax)
      ncompmax=floor(xncompmax*4.d0-2.d0)      
c Define quantile proportions = 1 - test sizes in vector *qprop* 
c of size *nqprop*
c Define the order statistics corresponding to quantile proprtions in  
c vector *iperc* ==> to be used for percentile computation of *nsim1*
c                    simulated statistics after sorting
      do j=1,nqprop
        qprop(j)=1.0-testsize(j)
        qprop2(j)=1.d0-testsize2(j)  
        iperc(j)=floor(qprop(j)*float(nsim1+1))
        if(iperc(j).lt.1) iperc(j)=1
        iperc2(j)=floor(testsize(j)*float(nsim1+1))
        if(iperc2(j).lt.1) iperc2(j)=1
      end do 
c Some checks are printed to screen
      write(*,*) nsim1,nsim2,n,idistr,iread,iapprox,rcont,isavesig,
     + ijit,idig1,itest 
      write(*,*) (testsize(j), j=1,nqprop)
      write(*,*) (qprop(j), j=1,nqprop)
      write(*,*) (iperc(j), j=1,nqprop)
      write(*,*) (iperc2(j), j=1,nqprop)      
c Compute asymptotic percentiles for *individual* tests (see list above)
c N.B.: For the Kolmogorov-Smirnov test do not compute percentiles but use 
c the asymptotic p-value from DKSONE
      do j=1,nqprop
        do jj=1,ntestind
          percas(j,jj)=0.d0
        end do
      end do
      do j=1,nqprop
        qq=qprop2(j)
        df=8.d0
        qqq=dchiin(qq,df)
        percas(j,1)=qqq
        df=89.d0
        qqq=dchiin(qq,df)
        percas(j,2)=qqq
        df=9.d0
        qqq=dchiin(qq,df)
        percas(j,3)=qqq
        percas(j,4)=0.d0
        df=1.d0
        qqq=dchiin(qq,df)
        percas(j,5)=qqq
        percas(j,6)=qqq     
        df=2.d0
        qqq=dchiin(qq,df)
        percas(j,7)=qqq
        percas(j,9)=qqq     
        percas(j,10)=qqq
        percas(j,12)=qqq
        df=4.d0
        qqq=dchiin(qq,df)
        percas(j,8)=qqq
        percas(j,11)=qqq 
c For Kuiper test VK asymptotic percentiles percas(j,13) are not computed:
c asymptotic test based on asymptotic p-value *pvalVKas* for the observed
c value of statistic VK (as for the output of DKSONE ==> see KS(6))
        percas(j,13)=0.d0
c The same as above for the Kolmogorov-Smirnov and Kuiper tests on the
c significand fractional part KSFR and VKFR (see asymptotic p-value *pvalVKFRas*)
        percas(j,14)=0.d0
        percas(j,15)=0.d0         
      end do
c A 64-bit Mersenne Twister generator is used (ioptMT=9)
      CALL RNOPT(ioptMT)
c Seed for random generation of Benford values is set automatically 
c (by the system clock: no CALL to RNSET)
c      open(unit=3,file='benf_iseed.txt')
c      read(3,*) iseed0
c      close(3)
c      CALL RNSET(iseed0)
      CALL RNGET(iseed0)
      write(*,*) iseed0
c      
C COMPUTES BENFORD PROBABILITIES AND RELATED QUANTITIES
c First digit and first-two digits
      do j=1,9
        dd=dfloat(j)
        benf1p(j)=dlog10(1.d0+1.d0/dd)
        benf1n(j)=benf1p(j)*fn
        do jj=1,10
          dd2=dfloat(jj)
          benf12p(j,jj)=dlog10(1.d0+1.d0/(10.d0*dd+dd2-1.d0))
          benf12n(j,jj)=benf12p(j,jj)*fn
        end do
      end do
c Moments for Hotelling test
      ee=DEXP(1.d0)    
      benfE=dlog10(ee)
      do j=1,9
        do jj=1,9
          benfVAR(j,jj)=-benfE**2
          if(j.eq.jj) benfVAR(j,j)=benfE*(dfloat(j)+0.5d0-benfE)
        end do
      end do
      nca=9
      nra=9
      lda=9
      ldginv=9
      CALL DLSGRR(nra,nca,benfVAR,lda,tol,irank,benfVARinv,ldginv)
C    
C IF *iread=\1*:
C * START OF SIMULATION LOOP OVER *nsim1* SIMULATIONS FOR COMPUTING TEST
C   STATISTICS UNDER THE NULL (BENFORD) HYPOTHESIS
C USE IMSL ROUTINES FOR RANDOM NUMBER GENERATION
C * SAVE SIMULATED TEST STATISTICS FOR SUBSEQUENT NULL DF-QUANTILE
C   ESTIMATION - Test statistics are saved in *single precision*
C OTHERWISE (iread=1): SKIP THIS STEP (READ *PERCENTILES* FROM FILE)
500   format(2f12.6)
      IF(iread.eq.1) THEN
        open(unit=51,file='benfX2_perc.txt')
        read(51,500) (percind(j,1),qprop(j),j=1,nqprop)
        close(51)
        open(unit=52,file='benfX2_2d_perc.txt')
        read(52,500) (percind(j,2),qprop(j),j=1,nqprop)
        close(52)
        open(unit=53,file='benfQ_perc.txt')
        read(53,500) (percind(j,3),qprop(j),j=1,nqprop)
        close(53)
        open(unit=54,file='benfKS_perc.txt')
        read(54,500) (percind(j,4),qprop(j),j=1,nqprop)
        close(54)
        open(unit=55,file='benfX2DIF_perc.txt')
        read(55,500) (percind(j,5),qprop(j),j=1,nqprop)
        close(55)
        open(unit=56,file='benfX2DIF2T_perc.txt')
        read(56,500) (percind(j,6),qprop(j),j=1,nqprop)
        close(56)
        open(unit=57,file='benfSCORE1_perc.txt')
        read(57,500) (percind(j,7),qprop(j),j=1,nqprop)
        close(57)
        open(unit=58,file='benfSCORE2_perc.txt')
        read(58,500) (percind(j,8),qprop(j),j=1,nqprop)
        close(58)
        open(unit=59,file='benfSCOREBIC_perc.txt')
        read(59,500) (percind(j,9),qprop(j),j=1,nqprop)
        close(59)
        open(unit=60,file='benfXLRT1_perc.txt')
        read(60,500) (percind(j,10),qprop(j),j=1,nqprop)
        close(60)
        open(unit=61,file='benfXLRT2_perc.txt')
        read(61,500) (percind(j,11),qprop(j),j=1,nqprop)
        close(61)
        open(unit=62,file='benfXLRTBIC_perc.txt')
        read(62,500) (percind(j,12),qprop(j),j=1,nqprop)
        close(62)
        open(unit=63,file='benfVK_perc.txt')
        read(63,500) (percind(j,13),qprop(j),j=1,nqprop)
        close(63)
        open(unit=64,file='benfKSFR_perc.txt')
        read(64,500) (percind(j,14),qprop(j),j=1,nqprop)
        close(64)
        open(unit=65,file='benfVKFR_perc.txt')
        read(65,500) (percind(j,15),qprop(j),j=1,nqprop)
        close(65)
c       
        open(unit=71,file='benfX2Q_perc.txt')
        read(71,500) (perccomb(j,1),testsize(j),j=1,nqprop)
        close(71)
        open(unit=72,file='benfQKS_perc.txt')
        read(72,500) (perccomb(j,2),testsize(j),j=1,nqprop)
        close(72)
        open(unit=73,file='benfX2QKS_perc.txt')
        read(73,500) (perccomb(j,3),testsize(j),j=1,nqprop)
        close(73)
        open(unit=74,file='benfX2DIFKS_perc.txt')
        read(74,500) (perccomb(j,4),testsize(j),j=1,nqprop)
        close(74)
        open(unit=75,file='benfX2DIF2TKS_perc.txt')
        read(75,500) (perccomb(j,5),testsize(j),j=1,nqprop)
        close(75)
        open(unit=76,file='benfX2DIFQSCORE_perc.txt')
        read(76,500) (perccomb(j,6),testsize(j),j=1,nqprop)
        close(76)
        open(unit=77,file='benfX2DIFQXLRT_perc.txt')
        read(77,500) (perccomb(j,7),testsize(j),j=1,nqprop)
        close(77)
        open(unit=78,file='benfX2DIFQVK_perc.txt')
        read(78,500) (perccomb(j,8),testsize(j),j=1,nqprop)
        close(78)
        open(unit=79,file='benfX2DIFKSFR_perc.txt')
        read(79,500) (perccomb(j,9),testsize(j),j=1,nqprop)
        close(79)
        open(unit=80,file='benfX2DIFVKFR_perc.txt')
        read(80,500) (perccomb(j,10),testsize(j),j=1,nqprop)
        close(80)
        GOTO 1000
      END IF
c Formats
20    format(i9)
30    Format(f15.7)
c increased number of digits for saving significands (isavesig=1)
300   format(f18.15)
301   format(f18.15,2x,f28.15)
302   format(f18.15,2x,f18.15,2x,f18.15,2x,f18.15)
630   format(2f20.9,8x,"Means")
631   format(2f20.9,8x,"Variances")
632   format(f20.9,8x,"Correlation")
633   format("Quadratic forms: X2 and Q",8x,"# sim =",i9,8x,
     + "n =",i6)
6322  format(8x,f12.9,8x,"Probability X2DIF < 0")
634   format("Quadratic forms: X2 and X2DIF = Q - X2",8x,"# sim =",i9,
     + 8x,"n =",i6)
63220 format(8x,f12.9,8x,"Probability X2DIF2T < 0")
6340  format("Quadratic forms: X2 and X2DIF2T = |Q - X2|",8x,"# sim =",
     + i9,8x,"n =",i6)
C      
      write(*,*) "NULL SIMULATION"
c Initialize some summary statistics
      checkX2DIF=0.d0
      checkX2DIF2T=0.d0
      avencompscoreb=0.d0
      varncompscoreb=0.d0
      maxncompscoreb=0.d0
      avencompxlrtb=0.d0
      varncompxlrtb=0.d0
      maxncompxlrtb=0.d0
      do jj=1,2
        X2Qmean(jj)=0.0
        X2DIFmean(jj)=0.0
        X2DIF2Tmean(jj)=0.0
        do j=1,2
          corX2Q(j,jj)=0.0
          corX2DIF(j,jj)=0.0
          corX2DIF2T(j,jj)=0.0
        end do
      end do
c Open files for saving test statistics
      OPEN(unit=31, file='benfX2.txt')
      OPEN(unit=32, file='benfX2_2d.txt')
      OPEN(unit=33, file='benfQ.txt')
      OPEN(unit=34, file='benfKS.txt')
      OPEN(unit=35, file='benfX2DIF.txt')
      OPEN(unit=36, file='benfX2DIF2T.txt')
      OPEN(unit=37, file='benfSCORE1.txt')
      OPEN(unit=38, file='benfSCORE2.txt')
      OPEN(unit=39, file='benfSCOREBIC.txt')
      OPEN(unit=40, file='benfXLRT1.txt')
      OPEN(unit=41, file='benfXLRT2.txt')
      OPEN(unit=42, file='benfXLRTBIC.txt')      
      OPEN(unit=43, file='benfVK.txt')      
      OPEN(unit=390, file='benfncompSCOREBIC.txt')
      OPEN(unit=420, file='benfncompXLRTBIC.txt')
      OPEN(unit=44, file='benfKSFR.txt')
      OPEN(unit=45, file='benfVKFR.txt')      
c    
c possibly save simulated data and significands - NOT CURRENT OPTION
c      open(unit=90,file='benf_simdata.txt')
c      open(unit=91,file='benf_simsig.txt')
c
C START OF SIMULATION LOOP OVER *nsim1* SIMULATIONS
      do jj=1,nsim1
        mm=mod(jj,loopn)
        if(mm.eq.0) write(*,*) "Benford sim # ",jj
c
C GENERATE A SAMPLE OF *n* OBSERVATIONS FROM A BENFORD RANDOM VARIABLE 
C AND COMPUTE SIGNIFICAND AND DIGITS FOR EACH OBSERVATION
        do i=1,n
          rn=DRNUNF()
          if(rn.le.0.d0.or.rn.ge.1.d0) then
            write(*,*) "error in rn =",rn," at unit #",i," in sim #",jj
            pause
          end if
          simdata(i)=rn
          s(i)=10.d0**simdata(i)
          d1(i)=floor(s(i))        
          if(d1(i).lt.1.or.d1(i).gt.9) then
             write(*,*) "error in d1"," at unit #",i," in sim #",jj
             pause
          end if
          ss=10.d0*s(i)
          d2(i)=floor(ss)-10.d0*floor(s(i))
          if(d2(i).lt.0.or.d2(i).gt.9) then
             write(*,*) "error in d2"," at unit #",i," in sim #",jj
             pause
           end if
          d12(i,1)=d1(i)
          d12(i,2)=d2(i)
c COMPUTE sFR(i) = <s(i)> = significand *fractional part*: to be used 
c                           in tests KSFR and VKFR
          sFR(i)=s(i)-floor(s(i))          
        end do
c possibly save data and significands
c        write(90,90) (simdata(i), i=1,n)
c        write(91,90) (s(i), i=1,n)
c90      format(f15.7) 
c       
C COMPUTE AND SAVE TEST STATISTICS ==> SAVE TO FILE IN ORDER TO REDUCE
C MEMORY REQUIREMENTS (only one dummy vector will be required for computing
C null percentiles) 
C EACH STATISTIC IS COMPUTED ON THE *SAME RANDOM SAMPLE* ==> THE SIMULATED 
C TEST VALUES CAN BE USED FOR COMBINATION OF TESTS
C Also compute summary statistics + correlation for X2, Q and X2DIF-X2DIF2T
c DOUBLE PRECISION IS USED FOR COMPUTING - SINGLE PRECISION FOR SAVING
        jsimul=jj
c Pearson X2 on first digit
        call X2ONE(jsimul,n,maxn,d1,benf1n,X2)
        write(31,30) X2
        X2Qmean(1)=X2Qmean(1)+X2/fn1
        corX2Q(1,1)=corX2Q(1,1)+X2*X2/fn1
        X2DIFmean(1)=X2DIFmean(1)+X2/fn1
        corX2DIF(1,1)=corX2DIF(1,1)+X2*X2/fn1
        X2DIF2Tmean(1)=X2DIFmean(1)
        corX2DIF2T(1,1)=corX2DIF(1,1)
c
c Pearson X2 on first-two digits
        call X2TWO(jsimul,n,maxn,d12,benf12n,X2_2d)
        write(32,30) X2_2d
c
c Hotelling Sum-invariance test Q
        call HOTSUM(jsimul,n,maxn,d1,s,benfE,benfVARinv,Q)
        write(33,30) Q 
        X2Qmean(2)=X2Qmean(2)+Q/fn1
        corX2Q(2,2)=corX2Q(2,2)+Q*Q/fn1
        corX2Q(1,2)=corX2Q(1,2)+X2*Q/fn1
c
c Kolmogorov-Smirnov test on the significand S(X)
c j=1,...,6 output values in vector KS(j) from IMSL routine
        call DKSONE(BENFCDF,n,s,KS,nmiss)
        write(34,34) (KS(j),j=1,6)
34      format(6f12.6)
c
c X2DIF = Q - X2
        X2DIF = Q - X2
        if(X2DIF.lt.0.d0) checkX2DIF = checkX2DIF + 1.d0/fn1
        write(35,30) X2DIF
        X2DIFmean(2)=X2DIFmean(2)+X2DIF/fn1
        corX2DIF(2,2)=corX2DIF(2,2)+X2DIF*X2DIF/fn1
        corX2DIF(1,2)=corX2DIF(1,2)+X2*X2DIF/fn1
c X2DIF2T = |Q - X2|
        X2DIF2T = dabs(Q - X2)
        if(X2DIF2T.lt.0.d0) checkX2DIF2T = checkX2DIF2T + 1.d0/fn1
        write(36,30) X2DIF2T
        X2DIF2Tmean(2)=X2DIF2Tmean(2)+X2DIF2T/fn1
        corX2DIF2T(2,2)=corX2DIF2T(2,2)+X2DIF2T*X2DIF2T/fn1
        corX2DIF2T(1,2)=corX2DIF2T(1,2)+X2*X2DIF2T/fn1
c        
c SCORE + LR TESTS ==> USE SUBROUTINE TEST() WITH 
c                      *ITERATIVE* BIC SELECTION ALGORITHM
        CALL TEST(jsimul,n,maxn,ncompmax,s,SCORE1,SCORE2,SCOREBIC,
     +    ncomps,XLRT1,XLRT2,XLRTBIC,ncompx)        
c        CALL TEST2(jsimul,n,maxn,ncompmax,s,SCORE1,SCORE2,SCOREBIC,
c     +    ncomps,XLRT1,XLRT2,XLRTBIC,ncompx)          
        avencompscoreb=avencompscoreb+dfloat(ncomps)/fn1
        varncompscoreb=varncompscoreb+dfloat(ncomps*ncomps)/fn1
        if(ncomps.eq.ncompmax) maxncompscoreb=maxncompscoreb+1.d0/fn1
        write(37,30) SCORE1
        write(38,30) SCORE2
        write(39,30) SCOREBIC
        write(390,20) ncomps
        avencompxlrtb=avencompxlrtb+dfloat(ncompx)/fn1
        varncompxlrtb=varncompxlrtb+dfloat(ncompx*ncompx)/fn1
        if(ncompx.eq.ncompmax) maxncompxlrtb=maxncompxlrtb+1.d0/fn1
        write(40,30) XLRT1
        write(41,30) XLRT2
        write(42,30) XLRTBIC
        write(420,20) ncompx
c        
c Kuiper test (based on output of DKSONE)
        VK=KS(2)+KS(3)
        write(43,30) VK
c
c Kolmogorov-Smirnov test on the significand *fractional part*
c j=1,...,6 output values in vector KS(j) from IMSL routine
        call DKSONE(BENFCDFFR,n,sFR,KSFR,nmiss)
        write(44,34) (KSFR(j),j=1,6)
c        
c Kuiper test on the significand *fractional part* (based on output of DKSONE)
        VKFR=KSFR(2)+KSFR(3)
        write(45,30) VKFR
c
C END OF SIMULATION LOOP OVER *nsim1* SIMULATIONS 
      end do 
c      close(90)
c      close(91)       
      CLOSE(31)
      CLOSE(32)
      CLOSE(33)
      CLOSE(34)
      CLOSE(35)  
      CLOSE(36)  
      CLOSE(37)
      CLOSE(38)
      CLOSE(39)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(390)
      CLOSE(420) 
      CLOSE(44)
      CLOSE(45)     
      write(*,*) "END OF NULL SIMULATION"
c
c Compute summary statistics + correlation beteween X2 and Q 
      do j=1,2
        corX2Q(j,j)=corX2Q(j,j)-(X2Qmean(j)*X2Qmean(j))
      end do
      corX2Q(1,2)=corX2Q(1,2)-(X2Qmean(1)*X2Qmean(2))
      corX2Q(1,2)=corX2Q(1,2)/((corX2Q(1,1)*corX2Q(2,2))**0.5)
      corX2Q(2,1)=corX2Q(1,2)
      open(unit=63,file='benfX2Q_summary.txt')
      write(63,630) (X2Qmean(j),j=1,2)
      write(63,*)
      write(63,631) (corX2Q(j,j),j=1,2)
      write(63,*)
      write(63,632) corX2Q(1,2)
      write(63,*)
      write(63,633) nsim1,n
      close(63)
c Compute summary statistics + correlation beteween X2 and X2DIF 
      do j=1,2
        corX2DIF(j,j)=corX2DIF(j,j)-(X2DIFmean(j)*X2DIFmean(j))
      end do
      corX2DIF(1,2)=corX2DIF(1,2)-(X2DIFmean(1)*X2DIFmean(2))
      corX2DIF(1,2)=corX2DIF(1,2)/((corX2DIF(1,1)*corX2DIF(2,2))**0.5)
      corX2DIF(2,1)=corX2DIF(1,2)
      open(unit=65,file='benfX2DIF_summary.txt')
      write(65,630) (X2DIFmean(j),j=1,2)
      write(65,*)
      write(65,631) (corX2DIF(j,j),j=1,2)
      write(65,*)
      write(65,632) corX2DIF(1,2)
      write(65,*)
      write(65,6322) checkX2DIF
      write(65,*)
      write(65,634) nsim1,n
      close(65)  
c Compute summary statistics + correlation beteween X2 and X2DIF2T 
      do j=1,2
        corX2DIF2T(j,j)=corX2DIF2T(j,j)-(X2DIF2Tmean(j)*X2DIF2Tmean(j))
      end do
      corX2DIF2T(1,2)=corX2DIF2T(1,2)-(X2DIF2Tmean(1)*X2DIF2Tmean(2))
      corX2DIF2T(1,2)=corX2DIF2T(1,2)/((corX2DIF2T(1,1)*
     + corX2DIF2T(2,2))**0.5)
      corX2DIF2T(2,1)=corX2DIF2T(1,2)
      open(unit=66,file='benfX2DIF2T_summary.txt')
      write(66,630) (X2DIF2Tmean(j),j=1,2)
      write(66,*)
      write(66,631) (corX2DIF2T(j,j),j=1,2)
      write(66,*)
      write(66,632) corX2DIF2T(1,2)
      write(66,*)
      write(66,6322) checkX2DIF2T
      write(66,*)
      write(66,634) nsim1,n
      close(66)  
c      
c
C
C READ SIMULATED TEST STATISTICS *FROM FILE* FOR NULL DF AND QUANTILE 
C ESTIMATION (if iread=\1)==> THIS STEP CAN BE AVOIDED OR IMPROVED IF 
C MORE MEMORY IS AVAILABLE AND FURTHER VECTORS CAN BE STORED
C
C *SORT* THE TEST STATISTICS + COMPUTE *ORDER STATISTICS* ==> COMPUTE AND SAVE 
C PERCENTILES FOR EACH TEST STATISTIC FOR SUBSEQUENT USE IN POWER COMPUTATION
C (see old versions for separate percentile computation using subroutine EQTIL)
C FOR COMPUTING THE EXACT P-VALUES:
C * USE ROUTINE *SVRGP* INSTEAD OF SVRGN
C * SAVE THE *SORTED* TEST STATISTICS IN FILES *sort.txt ==> SAVE THE
C   SORTED VALUES IN *DESCENDING* ORDER TO SAVE COMPUTING TIME 
C   (both here and in power computation)
C * SAVE THE *SORTING PERMUTATION* IN FILES *perm.txt ==> save in a separate 
C   file because the permutation is not (yet) required for power computation
C
C This read/save step can be avoided if further vectors can be stored in memory
C
C SINGLE PRECISION IS USED FOR PERCENTILE COMPUTATION IN ORDER TO ALLOW
C A LARGER VALUE OF *maxsim1*
C ONE LOOP FOR EACH TEST STATISTIC IN ORDER TO SAVE MEMORY SPACE: 
C use the dummy verctor *benfarray* for all test statistics 
C use the dummy verctor *ipermstat* for all test statistics 
      write(*,*) "START SORTING AND PERCENTILE COMPUTATION:"
c
c Start with Pearson X2 on first-two digits ==> do not save sorted values
c and do not compute the sorting permutation: not required (yet) for 
c combination of tests (use routine SVRGN instead of SVRGP)
      write(*,*) "FIRST-TWO DIGITS PEARSON TEST X2_2d"
      OPEN(unit=32, file='benfX2_2d.txt')
      do jj=1,nsim1
        read(32,30) X2_2ds
        benfarray(jj)=X2_2ds
      end do
      CLOSE(32)
      CALL SVRGN (nsim1, benfarray, benfarray)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)      
      do j=1,nqprop
        percind(j,2)=pp(j)
      end do
      open(unit=52,file='benfX2_2d_perc.txt')
      write(52,500) (percind(j,2),qprop(j),j=1,nqprop)
      write(52,*)
      write(52,502) nsim1,n
502   format("X2 First-two Digits",8x,"# sim =",i9,8x,"n =",i6)
      close(52)
c
c Pearson X2 on first digit
      write(*,*) "FIRST-DIGIT PEARSON TEST X2"
      OPEN(unit=31, file='benfX2.txt')
      do jj=1,nsim1
        read(31,30) X2s
        benfarray(jj)=X2s     
        ipermstat(jj)=jj        
      end do
      CLOSE(31)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=131,file='benfX2sort.txt')
      OPEN(unit=231,file='benfX2perm.txt')
      do jj=1,nsim1
        write(131,30) benfarray(nsim1-jj+1)
        write(231,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(131)
      CLOSE(231)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,1)=pp(j)
      end do
      open(unit=51,file='benfX2_perc.txt')
      write(51,500) (percind(j,1),qprop(j),j=1,nqprop)
      write(51,*)
      write(51,501) nsim1,n
501   format("X2 First Digit",8x,"# sim =",i9,8x,"n =",i6)
      close(51)
c
c Hotelling Sum-invariance test Q
      write(*,*) "HOTELLING SUM-INVARIANCE TEST Q"
      OPEN(unit=33, file='benfQ.txt')
      do jj=1,nsim1
        read(33,30) Qs
        benfarray(jj)=Qs
        ipermstat(jj)=jj
      end do
      CLOSE(33)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=133,file='benfQsort.txt')
      OPEN(unit=233,file='benfQperm.txt')      
      do jj=1,nsim1
        write(133,30) benfarray(nsim1-jj+1)
        write(233,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(133)
      CLOSE(233)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,3)=pp(j)
      end do
      open(unit=53,file='benfQ_perc.txt')
      write(53,500) (percind(j,3),qprop(j),j=1,nqprop)
      write(53,*)
      write(53,503) nsim1,n
503   format("Hotelling Sum-invariance test Q",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(53)
c
c Kolmogorov-Smirnov test on the significand S(X)
c j=1,...,6 output values in vector KS(j) from IMSL subroutine DKSONE
c ==> sorts the *first output* D_n = statistic for two-sided test
      write(*,*) "TWO-SIDED KOLMOGOROV-SMIRNOV TEST KS"
      OPEN(unit=34, file='benfKS.txt')
      do jj=1,nsim1
        read(34,34) (KSs(j),j=1,6)
        benfarray(jj)=KSs(1)
        ipermstat(jj)=jj        
      end do
      CLOSE(34)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=134,file='benfKSsort.txt')      
      OPEN(unit=234,file='benfKSperm.txt')
      do jj=1,nsim1
        write(134,30) benfarray(nsim1-jj+1)
        write(234,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(234)      
      CLOSE(134)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,4)=pp(j)
      end do
      open(unit=54,file='benfKS_perc.txt')
      write(54,500) (percind(j,4),qprop(j),j=1,nqprop)
      write(54,*)
      write(54,504) nsim1,n
504   format("Two-sided Kolmogorov-Smirnov test",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(54)
c
c X2DIF = Q - X2
      write(*,*) "X2DIF = Q - X2"
      OPEN(unit=35, file='benfX2DIF.txt')
      do jj=1,nsim1
        read(35,30) X2DIFs
        benfarray(jj)=X2DIFs
        ipermstat(jj)=jj        
      end do
      CLOSE(35)
      CALL SVRGP (nsim1, benfarray, benfarray,ipermstat)
      OPEN(unit=135,file='benfX2DIFsort.txt')
      OPEN(unit=235,file='benfX2DIFperm.txt')
      do jj=1,nsim1
        write(135,30) benfarray(nsim1-jj+1)
        write(235,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(135)
      CLOSE(235)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,5)=pp(j)
      end do
      open(unit=55,file='benfX2DIF_perc.txt')
      write(55,500) (percind(j,5),qprop(j),j=1,nqprop)
      write(55,*)
      write(55,505) nsim1,n
505   format("X2DIF = X2 - Q",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(55)
c
c X2DIF2T = |Q - X2|
c Use the previous vector *ipermX2DIF* to save memory space - NO 24/11/2021
      write(*,*) "X2DIF2T = |Q - X2|"
      OPEN(unit=36, file='benfX2DIF2T.txt')
      do jj=1,nsim1
        read(36,30) X2DIF2Ts
        benfarray(jj)=X2DIF2Ts
        ipermstat(jj)=jj        
      end do
      CLOSE(36)
      CALL SVRGP (nsim1, benfarray, benfarray,ipermstat)
      OPEN(unit=136,file='benfX2DIF2Tsort.txt')
      OPEN(unit=236,file='benfX2DIF2Tperm.txt')
      do jj=1,nsim1
        write(136,30) benfarray(nsim1-jj+1)
        write(236,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(136)
      CLOSE(236)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,6)=pp(j)
      end do
      open(unit=56,file='benfX2DIF2T_perc.txt')
      write(56,500) (percind(j,6),qprop(j),j=1,nqprop)
      write(56,*)
      write(56,506) nsim1,n
506   format("X2DIF = |Q - X2|",8x,"# sim =",i9,8x,"n =",i6)
      close(56) 
c 
c SCORE AND LIKELIHOOD RATIO (INDIVIDUAL) TESTS 
c NEW UNIT NUMBERS FOR PERCENTILE FILES
      write(*,*) "SCORE TEST WITH 1 COMPONENT"
      OPEN(unit=37, file='benfSCORE1.txt')
      do jj=1,nsim1
        read(37,30) SCORE1s
        benfarray(jj)=SCORE1s
        ipermstat(jj)=jj       
      end do
      CLOSE(37)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=137,file='benfSCORE1sort.txt')
      OPEN(unit=237,file='benfSCORE1perm.txt')      
      do jj=1,nsim1
        write(137,30) benfarray(nsim1-jj+1)
        write(237,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(137)
      CLOSE(237)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,7)=pp(j)
      end do
      open(unit=57,file='benfSCORE1_perc.txt')
      write(57,500) (percind(j,7),qprop(j),j=1,nqprop)
      write(57,*)
      write(57,507) nsim1,n
507   format("Score test with 1 component SCORE1",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(57)
c
      write(*,*) "SCORE TEST WITH 2 COMPONENTS"
      OPEN(unit=38, file='benfSCORE2.txt')
      do jj=1,nsim1
        read(38,30) SCORE2s
        benfarray(jj)=SCORE2s
        ipermstat(jj)=jj       
      end do
      CLOSE(38)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=138,file='benfSCORE2sort.txt')
      OPEN(unit=238,file='benfSCORE2perm.txt')      
      do jj=1,nsim1
        write(138,30) benfarray(nsim1-jj+1)
        write(238,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(138)
      CLOSE(238)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,8)=pp(j)
      end do
      open(unit=58,file='benfSCORE2_perc.txt')
      write(58,500) (percind(j,8),qprop(j),j=1,nqprop)
      write(58,*)
      write(58,508) nsim1,n
508   format("Score test with 2 components SCORE2",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(58)
c
      write(*,*) "SCORE TEST WITH BIC SELECTION OF COMPONENTS"
      OPEN(unit=39, file='benfSCOREBIC.txt')
      do jj=1,nsim1
        read(39,30) SCOREBICs
        benfarray(jj)=SCOREBICs
        ipermstat(jj)=jj       
      end do
      CLOSE(39)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=139,file='benfSCOREBICsort.txt')
      OPEN(unit=239,file='benfSCOREBICperm.txt')      
      do jj=1,nsim1
        write(139,30) benfarray(nsim1-jj+1)
        write(239,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(139)
      CLOSE(239)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,9)=pp(j)
      end do
      open(unit=59,file='benfSCOREBIC_perc.txt')
      write(59,500) (percind(j,9),qprop(j),j=1,nqprop)
      write(59,*)
      write(59,509) nsim1,n
509   format("Score test with BIC selection SCOREBIC",8x,"# sim =",i9,
     + 8x,"n =",i6)
      close(59)
c
      write(*,*) "LRT WITH 1 COMPONENT"
      OPEN(unit=40, file='benfXLRT1.txt')
      do jj=1,nsim1
        read(40,30) XLRT1s
        benfarray(jj)=XLRT1s
        ipermstat(jj)=jj       
      end do
      CLOSE(40)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=140,file='benfXLRT1sort.txt')
      OPEN(unit=240,file='benfXLRT1perm.txt')      
      do jj=1,nsim1
        write(140,30) benfarray(nsim1-jj+1)
        write(240,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(140)
      CLOSE(240)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,10)=pp(j)
      end do
      open(unit=60,file='benfXLRT1_perc.txt')
      write(60,500) (percind(j,10),qprop(j),j=1,nqprop)
      write(60,*)
      write(60,510) nsim1,n
510   format("LRT with 1 component XLRT1",8x,"# sim =",i9,8x,"n =",i6)
      close(60)
c
      write(*,*) "LRT WITH 2 COMPONENTS"
      OPEN(unit=41, file='benfXLRT2.txt')
      do jj=1,nsim1
        read(41,30) XLRT2s
        benfarray(jj)=XLRT2s
        ipermstat(jj)=jj       
      end do
      CLOSE(41)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=141,file='benfXLRT2sort.txt')
      OPEN(unit=241,file='benfXLRT2perm.txt')      
      do jj=1,nsim1
        write(141,30) benfarray(nsim1-jj+1)
        write(241,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(141)
      CLOSE(241)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,11)=pp(j)
      end do
      open(unit=61,file='benfXLRT2_perc.txt')
      write(61,500) (percind(j,11),qprop(j),j=1,nqprop)
      write(61,*)
      write(61,511) nsim1,n
511   format("LRT with 2 components XLRT2",8x,"# sim =",i9,8x,"n =",i6)
      close(61)
c
      write(*,*) "LRT WITH BIC SELECTION OF COMPONENTS"
      OPEN(unit=42, file='benfXLRTBIC.txt')
      do jj=1,nsim1
        read(42,30) XLRTBICs
        benfarray(jj)=XLRTBICs
        ipermstat(jj)=jj       
      end do
      CLOSE(42)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=142,file='benfXLRTBICsort.txt')
      OPEN(unit=242,file='benfXLRTBICperm.txt')      
      do jj=1,nsim1
        write(142,30) benfarray(nsim1-jj+1)
        write(242,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(142)
      CLOSE(242)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,12)=pp(j)
      end do
      open(unit=62,file='benfXLRTBIC_perc.txt')
      write(62,500) (percind(j,12),qprop(j),j=1,nqprop)
      write(62,*)
      write(62,512) nsim1,n
512   format("LRT with BIC selction XLRTBIC",8x,"# sim =",i9,8x,"n =",
     + i6)
      close(62)      
c
c Kuiper test VK
      write(*,*) "KUIPER TEST VK"
      OPEN(unit=43, file='benfVK.txt')
      do jj=1,nsim1
        read(43,30) VKs
        benfarray(jj)=VKs
        ipermstat(jj)=jj
      end do
      CLOSE(43)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=143,file='benfVKsort.txt')
      OPEN(unit=243,file='benfVKperm.txt')      
      do jj=1,nsim1
        write(143,30) benfarray(nsim1-jj+1)
        write(243,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(143)
      CLOSE(243)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,13)=pp(j)
      end do
      open(unit=63,file='benfVK_perc.txt')
      write(63,500) (percind(j,13),qprop(j),j=1,nqprop)
      write(63,*)
      write(63,603) nsim1,n
603   format("Kuiper test VK",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(63) 
c
c Kolmogorov-Smirnov test on the significand *fractional part*
c j=1,...,6 output values in vector KSFR(j) from IMSL subroutine DKSONE
c ==> sorts the *first output* D_n = statistic for two-sided test
      write(*,*) "TWO-SIDED KOLMOGOROV-SMIRNOV TEST KSFR ON THE ",
     + "SIGNIFICAND FRACTIONAL PART"
      OPEN(unit=44, file='benfKSFR.txt')
      do jj=1,nsim1
        read(44,34) (KSFRs(j),j=1,6)
        benfarray(jj)=KSFRs(1)
        ipermstat(jj)=jj        
      end do
      CLOSE(44)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=144,file='benfKSFRsort.txt')      
      OPEN(unit=244,file='benfKSFRperm.txt')
      do jj=1,nsim1
        write(144,30) benfarray(nsim1-jj+1)
        write(244,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(244)      
      CLOSE(144)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,14)=pp(j)
      end do
      open(unit=64,file='benfKSFR_perc.txt')
      write(64,500) (percind(j,14),qprop(j),j=1,nqprop)
      write(64,*)
      write(64,604) nsim1,n
604   format("Two-sided Kolmogorov-Smirnov test on the significand ",
     + "*fractional part*",8x,"# sim =",i9,8x,"n =",i6)
      close(64)
c
c Kuiper test VKFR on the significand *fractional part*
      write(*,*) "KUIPER TEST VKFR ON THE SIGNIFICAND FRACTIONAL PART"
      OPEN(unit=45, file='benfVKFR.txt')
      do jj=1,nsim1
        read(45,30) VKFRs
        benfarray(jj)=VKFRs
        ipermstat(jj)=jj
      end do
      CLOSE(45)
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=145,file='benfVKFRsort.txt')
      OPEN(unit=245,file='benfVKFRperm.txt')      
      do jj=1,nsim1
        write(145,30) benfarray(nsim1-jj+1)
        write(245,20) ipermstat(nsim1-jj+1)
      end do
      CLOSE(145)
      CLOSE(245)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc,pp, nmiss)
      do j=1,nqprop
        percind(j,15)=pp(j)
      end do
      open(unit=65,file='benfVKFR_perc.txt')
      write(65,500) (percind(j,15),qprop(j),j=1,nqprop)
      write(65,*)
      write(65,605) nsim1,n
605   format("Kuiper test VK on the significand *fractional part*",
     + 8x,"# sim =",i9,8x,"n =",i6)
      close(65) 
c 
c
C
C READ SIMULATED + SORTED TEST STATISTICS *FROM FILE* AGAIN FOR EACH TEST
C TO BE *COMBINED* ==> THIS STEP CAN BE AVOIDED OR IMPROVED IF MORE MEMORY 
C IS AVAILABLE AND FURTHER VECTORS CAN BE STORED 
C *NOT EFFICIENT*: COMPUTE THE QUANTILES OF EACH COMBINED TEST *SEPARATELY*
C IN ORDER TO AVOID STORAGE OF ADDITIONAL DUMMY VECTORS 
C COMPUTE QUANTILES OF P-VALUES ==> SORT AND SAVE THE P-VALUES (AS FOR THE
C OTHER TEST STATISTICS) BUT NOW IN *ASCENDING ORDER* 
C Use *ipermstat* as dummy vector for the sorting permutation
C Sorted test statistics are read but not (yet) required ==> use *iperm* only 
C N.B.: Sorting permutation of p-values not (yet) required for power computation
C USE *iperc2* INSTEAD OF *iperc* TO DEFINE ORDER STATISTICS OF P-VALUES
C (see old versions for separate percentile computation using subroutine EQTIL)
c Check unit numbers of files _perc.txt: TO DO
c
c Combined Test X2 + Q
      icomb(1)=0
      write(*,*) "COMBINED TEST X2 + Q"
      OPEN(unit=131, file='benfX2sort.txt')
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=231, file='benfX2perm.txt')
      OPEN(unit=233, file='benfQperm.txt')
      do jj=1,nsim1
        read(131,30) X2ss
        read(231,20) iperm1
        read(133,30) Qss
        read(233,20) iperm2
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(131)
      CLOSE(133)
      CLOSE(231)
      CLOSE(233)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=138,file='benfX2Qsort.txt')
      OPEN(unit=238,file='benfX2Qperm.txt')
      do jj=1,nsim1
        write(138,30) benfarray(jj)
        write(238,20) ipermstat(jj)        
      end do
      CLOSE(138)
      CLOSE(238)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,1)=pp(j)
      end do
      open(unit=71,file='benfX2Q_perc.txt')
      write(71,500) (perccomb(j,1),testsize(j),j=1,nqprop)
      write(71,*)
      write(71,580) nsim1,n
580   format("Percentiles of combined test X2 + Q",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(71)
c        
c Combined Test Q + KS
      write(*,*) "COMBINED TEST Q + KS"
      icomb(2)=0
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=134, file='benfKSsort.txt')
      OPEN(unit=233, file='benfQperm.txt')
      OPEN(unit=234, file='benfKSperm.txt')
      do jj=1,nsim1        
        read(133,30) Qss
        read(233,20) iperm1
        read(134,30) KSss
        read(234,20) iperm2
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)   
      end do
      CLOSE(233)
      CLOSE(234)
      CLOSE(133)
      CLOSE(134)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj        
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=139,file='benfQKSsort.txt')
      OPEN(unit=239,file='benfQKSperm.txt')
      do jj=1,nsim1
        write(139,30) benfarray(jj)
        write(239,20) ipermstat(jj)
      end do
      CLOSE(139)
      CLOSE(239)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,2)=pp(j)
      end do
      open(unit=72,file='benfQKS_perc.txt')
      write(72,500) (perccomb(j,2),testsize(j),j=1,nqprop)
      write(72,*)
      write(72,590) nsim1,n
590   format("Percentiles of combined test Q + KS",8x,"# sim =",i9,8x,
     + "n =",i6)
      close(72)
c
c Combined Test X2 + Q + KS
      icomb(3)=0
      write(*,*) "COMBINED TEST X2 + Q + KS"
      OPEN(unit=131, file='benfX2sort.txt')
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=134, file='benfKSsort.txt')
      OPEN(unit=231, file='benfX2perm.txt')
      OPEN(unit=233, file='benfQperm.txt')
      OPEN(unit=234, file='benfKSperm.txt')
      do jj=1,nsim1
        read(131,30) X2ss
        read(231,20) iperm1
        read(133,30) Qss
        read(233,20) iperm2
        read(134,30) KSss
        read(234,20) iperm3
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
        pval3(iperm3)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(131)
      CLOSE(133)
      CLOSE(134)
      CLOSE(231)
      CLOSE(233)
      CLOSE(234)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        if(pval3(jj).le.pvalmin) pvalmin=pval3(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=140,file='benfX2QKSsort.txt')
      OPEN(unit=240,file='benfX2QKSperm.txt')
      do jj=1,nsim1
        write(140,30) benfarray(jj)
        write(240,20) ipermstat(jj)
      end do
      CLOSE(140)
      CLOSE(240)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,3)=pp(j)
      end do
      open(unit=73,file='benfX2QKS_perc.txt')
      write(73,500) (perccomb(j,3),testsize(j),j=1,nqprop)
      write(73,*)
      write(73,610) nsim1,n
610   format("Percentiles of combined test X2 + Q + KS",8x,"# sim =",i9,
     + 8x,"n =",i6)
      close(73)
c
c Combined Test X2DIF + Q + KS 
      icomb(4)=0
      write(*,*) "COMBINED TEST X2DIF + Q + KS"
      OPEN(unit=135, file='benfX2DIFsort.txt')
      OPEN(unit=134, file='benfKSsort.txt')
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=235, file='benfX2DIFperm.txt')
      OPEN(unit=234, file='benfKSperm.txt')
      OPEN(unit=233, file='benfQperm.txt')
      do jj=1,nsim1
        read(135,30) X2DIFss
        read(235,20) iperm1
        read(134,30) KSss
        read(234,20) iperm2
        read(133,30) Qss
        read(233,20) iperm3
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
        pval3(iperm3)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(135)
      CLOSE(134)
      CLOSE(133)
      CLOSE(235)
      CLOSE(234)
      CLOSE(233)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        if(pval3(jj).le.pvalmin) pvalmin=pval3(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=141,file='benfX2DIFKSsort.txt')
      OPEN(unit=241,file='benfX2DIFKSperm.txt')
      do jj=1,nsim1
        write(141,30) benfarray(jj)
        write(241,20) ipermstat(jj)
      end do
      CLOSE(141)
      CLOSE(241)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,4)=pp(j)
      end do
      open(unit=74,file='benfX2DIFKS_perc.txt')
      write(74,500) (perccomb(j,4),testsize(j),j=1,nqprop)
      write(74,*)
      write(74,611) nsim1,n
611   format("Percentiles of combined test X2DIF + Q + KS",8x,
     + "# sim =",i9,8x,"n =",i6)
      close(74)  
c
c Combined Test X2DIF2T + Q + KS 
      icomb(5)=0      
      write(*,*) "COMBINED TEST X2DIF2T + Q + KS"
      OPEN(unit=136, file='benfX2DIF2Tsort.txt')
      OPEN(unit=134, file='benfKSsort.txt')
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=236, file='benfX2DIF2Tperm.txt')
      OPEN(unit=234, file='benfKSperm.txt')
      OPEN(unit=233, file='benfQperm.txt')
      do jj=1,nsim1
        read(136,30) X2DIF2Tss
        read(236,20) iperm1
        read(134,30) KSss
        read(234,20) iperm2
        read(133,30) Qss
        read(233,20) iperm3
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
        pval3(iperm3)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(136)
      CLOSE(134)
      CLOSE(133)
      CLOSE(236)
      CLOSE(234)
      CLOSE(233)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        if(pval3(jj).le.pvalmin) pvalmin=pval3(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj        
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=142,file='benfX2DIF2TKSsort.txt')
      OPEN(unit=242,file='benfX2DIF2TKSperm.txt')
      do jj=1,nsim1
        write(142,30) benfarray(jj)
        write(242,20) ipermstat(jj)
      end do
      CLOSE(142)
      CLOSE(242)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,5)=pp(j)
      end do
      open(unit=75,file='benfX2DIF2TKS_perc.txt')
      write(75,500) (perccomb(j,5),testsize(j),j=1,nqprop)
      write(75,*)
      write(75,612) nsim1,n
612   format("Percentiles of combined test X2DIF2T + Q + KS",8x,
     + "# sim =",i9,8x,"n =",i6)
      close(75)  
c
c Combined Test X2DIF + Q + SCOREBIC 
      icomb(6)=0
      write(*,*) "COMBINED TEST X2DIF + Q + SCOREBIC"
      OPEN(unit=136, file='benfX2DIFsort.txt')
      OPEN(unit=134, file='benfSCOREBICsort.txt')
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=236, file='benfX2DIFperm.txt')
      OPEN(unit=234, file='benfSCOREBICperm.txt')
      OPEN(unit=233, file='benfQperm.txt')
      do jj=1,nsim1
        read(136,30) X2DIFss
        read(236,20) iperm1
        read(134,30) SCOREBICss
        read(234,20) iperm2
        read(133,30) Qss
        read(233,20) iperm3
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
        pval3(iperm3)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(136)
      CLOSE(134)
      CLOSE(133)
      CLOSE(236)
      CLOSE(234)
      CLOSE(233)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        if(pval3(jj).le.pvalmin) pvalmin=pval3(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj        
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=143,file='benfX2DIFQSCOREsort.txt')
      OPEN(unit=243,file='benfX2DIFQSCOREperm.txt')
      do jj=1,nsim1
        write(143,30) benfarray(jj)
        write(243,20) ipermstat(jj)
      end do
      CLOSE(143)
      CLOSE(243)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,6)=pp(j)
      end do
      open(unit=76,file='benfX2DIFQSCORE_perc.txt')
      write(76,500) (perccomb(j,6),testsize(j),j=1,nqprop)
      write(76,*)
      write(76,613) nsim1,n
613   format("Percentiles of combined test X2DIF + Q + SCOREBIC",8x,
     + "# sim =",i9,8x,"n =",i6)
      close(76)  
c
c Combined Test X2DIF + Q + XLRTBIC 
      icomb(7)=0
      write(*,*) "COMBINED TEST X2DIF + Q + XLRTBIC"
      OPEN(unit=136, file='benfX2DIFsort.txt')
      OPEN(unit=134, file='benfXLRTBICsort.txt')
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=236, file='benfX2DIFperm.txt')
      OPEN(unit=234, file='benfXLRTBICperm.txt')
      OPEN(unit=233, file='benfQperm.txt')
      do jj=1,nsim1
        read(136,30) X2DIFss
        read(236,20) iperm1
        read(134,30) XLRTBICss
        read(234,20) iperm2
        read(133,30) Qss
        read(233,20) iperm3
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
        pval3(iperm3)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(136)
      CLOSE(134)
      CLOSE(133)
      CLOSE(236)
      CLOSE(234)
      CLOSE(233)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        if(pval3(jj).le.pvalmin) pvalmin=pval3(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj        
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=144,file='benfX2DIFQXLRTsort.txt')
      OPEN(unit=244,file='benfX2DIFQXLRTperm.txt')
      do jj=1,nsim1
        write(144,30) benfarray(jj)
        write(244,20) ipermstat(jj)
      end do
      CLOSE(144)
      CLOSE(244)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,7)=pp(j)
      end do
      open(unit=77,file='benfX2DIFQXLRT_perc.txt')
      write(77,500) (perccomb(j,7),testsize(j),j=1,nqprop)
      write(77,*)
      write(77,614) nsim1,n
614   format("Percentiles of combined test X2DIF + Q + XLRTBIC",8x,
     + "# sim =",i9,8x,"n =",i6)
      close(77)  
c
c Combined Test X2DIF + Q + KUIPER 
      icomb(8)=0      
      write(*,*) "COMBINED TEST X2DIF + Q + KUIPER"
      OPEN(unit=136, file='benfX2DIFsort.txt')
      OPEN(unit=134, file='benfVKsort.txt')
      OPEN(unit=133, file='benfQsort.txt')
      OPEN(unit=236, file='benfX2DIFperm.txt')
      OPEN(unit=234, file='benfVKperm.txt')
      OPEN(unit=233, file='benfQperm.txt')
      do jj=1,nsim1
        read(136,30) X2DIFss
        read(236,20) iperm1
        read(134,30) VKss
        read(234,20) iperm2
        read(133,30) Qss
        read(233,20) iperm3
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
        pval3(iperm3)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(136)
      CLOSE(134)
      CLOSE(133)
      CLOSE(236)
      CLOSE(234)
      CLOSE(233)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        if(pval3(jj).le.pvalmin) pvalmin=pval3(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj        
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=145,file='benfX2DIFQVKsort.txt')
      OPEN(unit=245,file='benfX2DIFQVKperm.txt')
      do jj=1,nsim1
        write(145,30) benfarray(jj)
        write(245,20) ipermstat(jj)
      end do
      CLOSE(145)
      CLOSE(245)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,8)=pp(j)
      end do
      open(unit=78,file='benfX2DIFQVK_perc.txt')
      write(78,500) (perccomb(j,8),testsize(j),j=1,nqprop)
      write(78,*)
      write(78,615) nsim1,n
615   format("Percentiles of combined test X2DIF + Q + KUIPER",8x,
     + "# sim =",i9,8x,"n =",i6)
      close(78) 
c
c Combined Test X2DIF + KSFR
      icomb(9)=0
      write(*,*) "COMBINED TEST X2DIF + KSFR"
      OPEN(unit=135, file='benfX2DIFsort.txt')
      OPEN(unit=134, file='benfKSFRsort.txt')
      OPEN(unit=235, file='benfX2DIFperm.txt')
      OPEN(unit=234, file='benfKSFRperm.txt')
      do jj=1,nsim1
        read(135,30) X2DIFss
        read(235,20) iperm1
        read(134,30) KSFRss
        read(234,20) iperm2
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(135)
      CLOSE(134)
      CLOSE(235)
      CLOSE(234)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=141,file='benfX2DIFKSFRsort.txt')
      OPEN(unit=241,file='benfX2DIFKSFRperm.txt')
      do jj=1,nsim1
        write(141,30) benfarray(jj)
        write(241,20) ipermstat(jj)
      end do
      CLOSE(141)
      CLOSE(241)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,9)=pp(j)
      end do
      open(unit=79,file='benfX2DIFKSFR_perc.txt')
      write(79,500) (perccomb(j,9),testsize(j),j=1,nqprop)
      write(79,*)
      write(79,619) nsim1,n
619   format("Percentiles of combined test X2DIF + KSFR",8x,
     + "# sim =",i9,8x,"n =",i6)
      close(79) 
c
c Combined Test X2DIF + VKFR
      icomb(10)=0
      write(*,*) "COMBINED TEST X2DIF + VKFR"
      OPEN(unit=135, file='benfX2DIFsort.txt')
      OPEN(unit=134, file='benfVKFRsort.txt')
      OPEN(unit=235, file='benfX2DIFperm.txt')
      OPEN(unit=234, file='benfVKFRperm.txt')
      do jj=1,nsim1
        read(135,30) X2DIFss
        read(235,20) iperm1
        read(134,30) VKFRss
        read(234,20) iperm2
        pval1(iperm1)=float(jj-1)/float(nsim1)
        pval2(iperm2)=float(jj-1)/float(nsim1)      
      end do
      CLOSE(135)
      CLOSE(134)
      CLOSE(235)
      CLOSE(234)
      do jj=1,nsim1
        pvalmin=pval1(jj)
        if(pval2(jj).le.pvalmin) pvalmin=pval2(jj)
        benfarray(jj)=pvalmin
        ipermstat(jj)=jj
      end do
      CALL SVRGP (nsim1, benfarray, benfarray, ipermstat)
      OPEN(unit=141,file='benfX2DIFVKFRsort.txt')
      OPEN(unit=241,file='benfX2DIFVKFRperm.txt')
      do jj=1,nsim1
        write(141,30) benfarray(jj)
        write(241,20) ipermstat(jj)
      end do
      CLOSE(141)
      CLOSE(241)
      CALL ORDST (nsim1,benfarray,nqprop,0,iperc2,pp,nmiss)
      do j=1,nqprop
        perccomb(j,10)=pp(j)
      end do
      open(unit=80,file='benfX2DIFVKFR_perc.txt')
      write(80,500) (perccomb(j,10),testsize(j),j=1,nqprop)
      write(80,*)
      write(80,620) nsim1,n
620   format("Percentiles of combined test X2DIF + VKFR",8x,
     + "# sim =",i9,8x,"n =",i6)
      close(80)      
c
      write(*,*) "END OF PERCENTILE COMPUTATION"
1000  CONTINUE
C
C
C START POWER COMPUTATION
      if(perccomb(1,1).le.testsize(1)) icomb(1)=1
      if(perccomb(1,2).le.testsize(1)) icomb(2)=1
      if(perccomb(1,3).le.testsize(1)) icomb(3)=1      
      if(perccomb(1,4).le.testsize(1)) icomb(4)=1
      if(perccomb(1,5).le.testsize(1)) icomb(5)=1
      if(perccomb(1,6).le.testsize(1)) icomb(6)=1
      if(perccomb(1,7).le.testsize(1)) icomb(7)=1
      if(perccomb(1,8).le.testsize(1)) icomb(8)=1
      ncont=floor(rcont*fn)
      ngood=n-ncont
      if(idistr.ge.0) then
        write(*,*) "POWER SIMULATION"
        write(*,*) "Contamination rate",rcont
        write(*,*) "# of Benford values in each sample",ngood,
     +   " out of",n
        write(*,*) "# of contaminated values in each sample",
     +   ncont," out of",n
      else
        write(*,*) "EMPIRICAL ANALYISIS - READ DATA FROM FILE"
      end if
c Initialize power and other arrays
      do j=1,npar
        do jj=1,nqprop
          do jjj=1,ntestind
            powind(j,jj,jjj)=0.d0
            powas(j,jj,jjj)=0.d0
            if(j.eq.1) pind(jj,jjj)=dfloat(percind(jj,jjj))      
          end do        
          do jjj=1,ntestcomb
            powcomb(j,jj,jjj)=0.d0
          end do        
        end do      
        avencompscores(j)=0.d0
        varncompscores(j)=0.d0
        maxncompscores(j)=0.d0
        avencompxlrts(j)=0.d0
        varncompxlrts(j)=0.d0
        maxncompxlrts(j)=0.d0
      end do
c
c Set parameter values: npar = max # of parameters to be read from file
c                      par() = vector of parameter values (read from file)
c                   npartrue = actual # of parameters (distribution-specific)
c                    partrue = specific parameter value used in simulation 
c Compute distribution-specific nominal parameters to be checked in 
c simulations ==> nominal parameters assume that *rcont=1*
c ==> files checkmeant(j) and checkvart(j) for mean and sigma of real
c     value distributions (except Stable)
c ==> files checksigt(9,j) and checkdigt(9,j) for significand and digits
C     of Benford-type distributions
c     j=1,...,npartrue
c Open distribution-specific files powX2* for non null-correlations 
      do j=1,npar
        par(j)=0.d0
      end do
      xloc=0.d0
      scal=0.d0
c idistr < 0 ==> Read data from file + jittering in the 
c                10th decimal digit to avoid ties in DKSONE
c                (write warning if negative data)
c                jittering in the 4th decimal digit if the observed
c                value is large
      if(idistr.lt.0) then
        npartrue=1
        iapprox=0
        ineg=0
c * if idistr < -1 keep *nsim2* as in file .in1 ==> do NOT 
c   set nsim2=1 for analysis of multiple (simulated) data sets
c   (less efficient but compatible with previous code)
c    * read data only if *idistr=-1* ==> less efficient but 
c    compatible with previous code  (if idistr<-1 read
c    data within the loop over *nsim2*  simulations)
         if(idistr.eq.(-1)) then
           nsim2=1
           OPEN(unit=11,file='benf_data.txt')
           OPEN(unit=12,file='benf_jitter.txt')        
           do i=1,n
             read(11,*) simdata(i)
             if(simdata(i).lt.0.d0) then
               ineg=1
               write(*,*) "Warning: observation i is negative;  i = ",i,
     +          simdata(i)
             end if
             if(simdata(i).eq.0.d0) then
               write(*,*)
             write(*,*) "Warning: observation i is 0 ==> NOT ALLOWED; ",
     +         " i = ",i,simdata(i)
               write(*,*)
               go to 9999
             end if
c Perform jittering in the data only if *ijit.gt.0*
             if(ijit.gt.0) then
               rng=DRNUNF()                           
c More jittering if the observed value is large to avoid ties in DKSONE
c more jittering for GM + COVID data set     
               if(simdata(i).lt.10.d0**3) then
!               if(simdata(i).lt.10.d0**1) then
!                 simdata(i)=simdata(i)+rng*(10.d0**(-9))   
!                 write(12,11) i,simdata(i),rng,rng*(10.d0**(-9))  
                 simdata(i)=simdata(i)+rng*(10.d0**(-6))   
                 write(12,11) i,simdata(i),rng,rng*(10.d0**(-6))
!                 simdata(i)=simdata(i)+rng*(10.d0**(-5))   
!                 write(12,11) i,simdata(i),rng,rng*(10.d0**(-5))                 
               else if(simdata(i).lt.10.d0**4) then
!               else if(simdata(i).lt.10.d0**2) then
                 simdata(i)=simdata(i)+rng*(10.d0**(-4))   
                 write(12,11) i,simdata(i),rng,rng*(10.d0**(-4)) 
c More jittering to reduce the probability of ties if 
c the observed value is very large
               else if(simdata(i).lt.10.d0**5) then
!               else if(simdata(i).lt.10.d0**3) then
!                 simdata(i)=simdata(i)+rng*(10.d0**(-2))   
!                 write(12,11) i,simdata(i),rng,rng*(10.d0**(-2))
                 simdata(i)=simdata(i)+rng*(10.d0**(-3))   
                 write(12,11) i,simdata(i),rng,rng*(10.d0**(-3))
               else
                 simdata(i)=simdata(i)+rng*(10.d0**(0))   
                 write(12,11) i,simdata(i),rng,rng*(10.d0**(0))                
               end if
11             format(i9,f36.16,10x,f18.16,10x,f18.16)
             end if    
           end do
c OPEN FILES FOR SAVING SIGNIFICANDS
c ALSO SAVE DATA IN THE SAME FILE
           open(unit=101, file='signif_DATA.txt') 
c       
           CLOSE(11)
           CLOSE(12)
         end if
      end if      
c idistr=0: Benford distribution
      if(idistr.eq.0) then
        npartrue=1
        do j=1,npartrue
          do jj=1,9
            checksigt(jj,j)=benfE
            checkdigt(jj,j)=benf1p(jj)
          end do
        end do     
c if *idig1=1* 1st digit is Benford ==> different file names
        if(idig1.eq.1) then
          open(unit=73,file='powX2Q_summary_BENFBENF.txt')
          open(unit=75,file='powX2DIF_summary_BENFBENF.txt')
          open(unit=77,file='powX2DIF2T_summary_BENFBENF.txt')
          open(unit=101, file='signif_BENFBENF.txt')   
        else      
          open(unit=73,file='powX2Q_summary_BENF.txt')
          open(unit=75,file='powX2DIF_summary_BENF.txt')
          open(unit=77,file='powX2DIF2T_summary_BENF.txt')
          open(unit=101, file='signif_BENF.txt')   
        end if
      end if
c idistr=1: Normal distribution (3 scales * 2 locations)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.1) THEN
        OPEN(unit=11,file='par_norm.txt')
        npartrue=npar
        npar2=npartrue/2
        if(floor(float(npar2)).ne.npar2) then
          write(*,*) "npar not even: the selected distribution is ",
     +     "invalid"
          go to 9999
        end if
        read(11,*) (par(j),j=1,npar2)
        read(11,*) par(npar2+1),par(npar2+2)
        CLOSE(11)
        do j=1,npar2
          xloc=par(npar2+1)
          partrue=par(j)
          checkmeant(j)=xloc
          checkvart(j)=partrue
        end do
        do j=npar2+1,npartrue
          xloc=par(npar2+2)
          partrue=par(j-npar2)
          checkmeant(j)=xloc
          checkvart(j)=partrue
        end do
c if *idig1=1* 1st digit is Benford ==> different file names
        if(idig1.eq.1) then
          open(unit=73,file='powX2Q_summary_BENFNORM.txt')
          open(unit=75,file='powX2DIF_summary_BENFNORM.txt')
          open(unit=77,file='powX2DIF2T_summary_BENFNORM.txt')
          open(unit=101, file='signif_BENFNORM1.txt') 
          open(unit=102, file='signif_BENFNORM2.txt') 
          open(unit=103, file='signif_BENFNORM3.txt') 
          open(unit=104, file='signif_BENFNORM4.txt') 
          open(unit=105, file='signif_BENFNORM5.txt') 
          open(unit=106, file='signif_BENFNORM6.txt') 
        else
          open(unit=73,file='powX2Q_summary_NORM.txt')
          open(unit=75,file='powX2DIF_summary_NORM.txt')
          open(unit=77,file='powX2DIF2T_summary_NORM.txt')
          open(unit=101, file='signif_NORM1.txt') 
          open(unit=102, file='signif_NORM2.txt') 
          open(unit=103, file='signif_NORM3.txt') 
          open(unit=104, file='signif_NORM4.txt') 
          open(unit=105, file='signif_NORM5.txt') 
          open(unit=106, file='signif_NORM6.txt')         
        end if
      END IF         
c idistr=2: Lognormal distribution (fixed scale read from file)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.2) THEN
        npartrue=npar-1
        OPEN(unit=11,file='par_lognorm.txt')
        read(11,*) (par(j),j=1,npartrue)
        read(11,*) par(npartrue+1)
        CLOSE(11)
        scal=par(npartrue+1)
        do j=1,npartrue
          partrue=par(j)
          checkmeant(j)=scal*dexp(0.5d0*partrue*partrue)
          checkvart(j)=scal*scal*dexp(partrue*partrue)*
     +     (dexp(partrue*partrue)-1.d0) 
          checkvart(j)=dsqrt(checkvart(j))  
        end do
c if *idig1=1* 1st digit is Benford ==> different file names
        if(idig1.eq.1) then
          open(unit=73,file='powX2Q_summary_BENFLOGNORM.txt')
          open(unit=75,file='powX2DIF_summary_BENFLOGNORM.txt')        
          open(unit=77,file='powX2DIF2T_summary_BENFLOGNORM.txt')  
          open(unit=101, file='signif_BENFLOGNORM1.txt') 
          open(unit=102, file='signif_BENFLOGNORM2.txt') 
          open(unit=103, file='signif_BENFLOGNORM3.txt') 
          open(unit=104, file='signif_BENFLOGNORM4.txt') 
          open(unit=105, file='signif_BENFLOGNORM5.txt')      
        else
          open(unit=73,file='powX2Q_summary_LOGNORM.txt')
          open(unit=75,file='powX2DIF_summary_LOGNORM.txt')        
          open(unit=77,file='powX2DIF2T_summary_LOGNORM.txt')  
          open(unit=101, file='signif_LOGNORM1.txt') 
          open(unit=102, file='signif_LOGNORM2.txt') 
          open(unit=103, file='signif_LOGNORM3.txt') 
          open(unit=104, file='signif_LOGNORM4.txt') 
          open(unit=105, file='signif_LOGNORM5.txt')      
        end if        
      END IF         
c idistr=3: Weibull distribution (fixed scale parameter set = 1 by default)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.3) THEN
        npartrue=npar
        OPEN(unit=11,file='par_weib.txt')
        read(11,*) (par(j),j=1,npartrue)
        CLOSE(11)
        scal=1.d0
        do j=1,npartrue
          partrue=par(j)
          checkmeant(j)=scal*dgamma((partrue+1.d0)/partrue)
          gg1=dgamma((partrue+2.d0)/partrue)
          gg2=dgamma((partrue+1.d0)/partrue)
          checkvart(j)=scal*scal*(gg1-gg2*gg2)
          checkvart(j)=dsqrt(checkvart(j))
        end do
c if *idig1=1* 1st digit is Benford ==> different file names
        if(idig1.eq.1) then     
          open(unit=73,file='powX2Q_summary_BENFWEIB.txt')
          open(unit=75,file='powX2DIF_summary_BENFWEIB.txt')
          open(unit=77,file='powX2DIF2T_summary_BENFWEIB.txt')
          open(unit=101, file='signif_BENFWEIB1.txt') 
          open(unit=102, file='signif_BENFWEIB2.txt') 
          open(unit=103, file='signif_BENFWEIB3.txt') 
          open(unit=104, file='signif_BENFWEIB4.txt') 
          open(unit=105, file='signif_BENFWEIB5.txt') 
          open(unit=106, file='signif_BENFWEIB6.txt')         
        else
          open(unit=73,file='powX2Q_summary_WEIB.txt')
          open(unit=75,file='powX2DIF_summary_WEIB.txt')
          open(unit=77,file='powX2DIF2T_summary_WEIB.txt')
          open(unit=101, file='signif_WEIB1.txt') 
          open(unit=102, file='signif_WEIB2.txt') 
          open(unit=103, file='signif_WEIB3.txt') 
          open(unit=104, file='signif_WEIB4.txt') 
          open(unit=105, file='signif_WEIB5.txt') 
          open(unit=106, file='signif_WEIB6.txt')                 
        end if
      END IF         
c idistr=4: Gamma distribution (fixed scale parameter set = 1 by default)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.4) THEN
        npartrue=npar
        OPEN(unit=11,file='par_gam.txt')
        read(11,*) (par(j),j=1,npartrue)
        CLOSE(11)
        scal=1.d0
        do j=1,npartrue
          partrue=par(j)
          checkmeant(j)=scal*partrue
          checkvart(j)=scal*scal*partrue
          checkvart(j)=dsqrt(checkvart(j))        
        end do
c if *idig1=1* 1st digit is Benford ==> different file names  
        if(idig1.eq.1) then      
          open(unit=73,file='powX2Q_summary_BENFGAM.txt')
          open(unit=75,file='powX2DIF_summary_BENFGAM.txt')
          open(unit=77,file='powX2DIF2T_summary_BENFGAM.txt')
          open(unit=101, file='signif_BENFGAM1.txt') 
          open(unit=102, file='signif_BENFGAM2.txt') 
          open(unit=103, file='signif_BENFGAM3.txt') 
          open(unit=104, file='signif_BENFGAM4.txt') 
          open(unit=105, file='signif_BENFGAM5.txt') 
          open(unit=106, file='signif_BENFGAM6.txt')          
        else
          open(unit=73,file='powX2Q_summary_GAM.txt')
          open(unit=75,file='powX2DIF_summary_GAM.txt')
          open(unit=77,file='powX2DIF2T_summary_GAM.txt')
          open(unit=101, file='signif_GAM1.txt') 
          open(unit=102, file='signif_GAM2.txt') 
          open(unit=103, file='signif_GAM3.txt') 
          open(unit=104, file='signif_GAM4.txt') 
          open(unit=105, file='signif_GAM5.txt') 
          open(unit=106, file='signif_GAM6.txt')                  
        end if
      END IF         
c idistr=5: Beta distribution (3 values of alpha * 2 values of beta)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.5) THEN
        OPEN(unit=11,file='par_beta.txt')
        npartrue=npar
        npar2=npartrue/2
        if(floor(float(npar2)).ne.npar2) then
          write(*,*) "npar not even: the selected distribution is ",
     +     "invalid"
          go to 9999
        end if
        read(11,*) (par(j),j=1,npar2)
        read(11,*) par(npar2+1),par(npar2+2)
        CLOSE(11)        
        do j=1,npar2
          xloc=par(npar2+1)
          partrue=par(j)
          checkmeant(j)=partrue/(partrue+xloc)
          checkvart(j)=partrue*xloc/((partrue+xloc)*(partrue+xloc)*
     +     (partrue+xloc+1.d0))
          checkvart(j)=dsqrt(checkvart(j))     
        end do
        do j=npar2+1,npartrue
          xloc=par(npar2+2)
          partrue=par(j-npar2)
          checkmeant(j)=partrue/(partrue+xloc)
          checkvart(j)=partrue*xloc/((partrue+xloc)*(partrue+xloc)*
     +     (partrue+xloc+1.d0))
          checkvart(j)=dsqrt(checkvart(j))
        end do
c if *idig1=1* 1st digit is Benford ==> different file names
        if(idig1.eq.1) then        
          open(unit=73,file='powX2Q_summary_BENFBETA.txt')
          open(unit=75,file='powX2DIF_summary_BENFBETA.txt')
          open(unit=77,file='powX2DIF2T_summary_BENFBETA.txt')
          open(unit=101, file='signif_BENFBETA1.txt') 
          open(unit=102, file='signif_BENFBETA2.txt') 
          open(unit=103, file='signif_BENFBETA3.txt') 
          open(unit=104, file='signif_BENFBETA4.txt') 
          open(unit=105, file='signif_BENFBETA5.txt') 
          open(unit=106, file='signif_BENFBETA6.txt')                  
        else
          open(unit=73,file='powX2Q_summary_BETA.txt')
          open(unit=75,file='powX2DIF_summary_BETA.txt')
          open(unit=77,file='powX2DIF2T_summary_BETA.txt')
          open(unit=101, file='signif_BETA1.txt') 
          open(unit=102, file='signif_BETA2.txt') 
          open(unit=103, file='signif_BETA3.txt') 
          open(unit=104, file='signif_BETA4.txt') 
          open(unit=105, file='signif_BETA5.txt') 
          open(unit=106, file='signif_BETA6.txt')                          
        end if
      END IF    
c idistr=6: Uniform distribution on (0,a)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.6) THEN
        npartrue=npar
        OPEN(unit=11,file='par_unif.txt')
        read(11,*) (par(j),j=1,npartrue)
        CLOSE(11)
        do j=1,npartrue
          partrue=par(j)
          checkmeant(j)=partrue/2.d0
          checkvart(j)=partrue*partrue/12.d0
        end do
c if *idig1=1* 1st digit is Benford ==> different file names
        if(idig1.eq.1) then        
          open(unit=73,file='powX2Q_summary_BENFUNIF.txt')
          open(unit=75,file='powX2DIF_summary_BENFUNIF.txt')
          open(unit=77,file='powX2DIF2T_summary_BENFUNIF.txt')
          open(unit=101, file='signif_BENFUNIF1.txt') 
          open(unit=102, file='signif_BENFUNIF2.txt') 
          open(unit=103, file='signif_BENFUNIF3.txt') 
          open(unit=104, file='signif_BENFUNIF4.txt') 
          open(unit=105, file='signif_BENFUNIF5.txt') 
          open(unit=106, file='signif_BENFUNIF6.txt')          
        else
          open(unit=73,file='powX2Q_summary_UNIF.txt')
          open(unit=75,file='powX2DIF_summary_UNIF.txt')
          open(unit=77,file='powX2DIF2T_summary_UNIF.txt')
          open(unit=101, file='signif_UNIF1.txt') 
          open(unit=102, file='signif_UNIF2.txt') 
          open(unit=103, file='signif_UNIF3.txt') 
          open(unit=104, file='signif_UNIF4.txt') 
          open(unit=105, file='signif_UNIF5.txt') 
          open(unit=106, file='signif_UNIF6.txt')                  
        end if
      END IF
c idistr=7: Stable distribution - 5 shape parameter values + 1 (fixed) 
c                                 skewness parameter value
c                                 (Assume location = 0 and scale = 1)
c Instead of nominal paramemeter values mean and sigma, compute the  
c nominal tail probability for tail value tolstab(j), j=1,...,npartrue
c ==> checkstabt(j)
c TO DO: compute nominal significand + digits 
c Choice of the simulation algorithm for Stable distribution:
c read istab ==> istab=1: IMSL routine DRNSTA from Chambers et al (JASA, 1976); 
c                         uses both alpha and beta
c                istab=2: Kanter algorithm for positive stable distrib. adapted by XXXX; 
c                         uses alpha only
c                istab=3: Kanter algorithm for positive stable distrib. adapted by Chambers et al.; 
c                         uses alpha only 
      IF(idistr.eq.7) THEN
        npartrue=npar-1
        OPEN(unit=11,file='par_stab.txt')
        read(11,*) (par(j),j=1,npartrue)
        read(11,*) par(npartrue+1)
        read(11,*) istab
        CLOSE(11)
        do j=1,npartrue
          if(par(j).le.0.d0.or.par(j).gt.2.d0) then
            write(*,*) "Error: invalid alpha of Stable distribution"
            go to 9999
          end if
        end do  
        if(par(npartrue+1).lt.(-1.d0).or.par(npartrue+1).gt.1.d0) then
          write(*,*) "Error: invalid beta of Stable distribution"
          go to 9999
        end if
        do j=1,npartrue
          partrue=par(j)
          xloc=par(npartrue+1)
          scal=1.d0
c default value for tolstab: 10**12 when alpha=0.1 ==> the values of tolstab(j)
c are chosen in order to keep approximately similar tail probabilities 
c as the default value (approx 0.05, but note the the proportionality constant 
c changes with alpha ==> the nominal tail probability ranges from 0.04 to 0.06)
          tolstab(j)=(1.d12)**(1.d-1/partrue)
          tprob=tailprob(partrue,xloc,scal,tolstab(j))
          checkstabt(j)=tprob
        end do
c if *idig1=1* 1st digit is Benford ==> different file names        
        if(idig1.eq.1) then
          open(unit=73,file='powX2Q_summary_BENFSTAB.txt')
          open(unit=75,file='powX2DIF_summary_BENFSTAB.txt')        
          open(unit=77,file='powX2DIF2T_summary_BENFSTAB.txt')        
          open(unit=101, file='signif_BENFSTAB1.txt') 
          open(unit=102, file='signif_BENFSTAB2.txt') 
          open(unit=103, file='signif_BENFSTAB3.txt') 
          open(unit=104, file='signif_BENFSTAB4.txt') 
          open(unit=105, file='signif_BENFSTAB5.txt') 
        else
          open(unit=73,file='powX2Q_summary_STAB.txt')
          open(unit=75,file='powX2DIF_summary_STAB.txt')        
          open(unit=77,file='powX2DIF2T_summary_STAB.txt')        
          open(unit=101, file='signif_STAB1.txt') 
          open(unit=102, file='signif_STAB2.txt') 
          open(unit=103, file='signif_STAB3.txt') 
          open(unit=104, file='signif_STAB4.txt') 
          open(unit=105, file='signif_STAB5.txt')         
        end if
      END IF    
c idistr=8: Generalized Benford distribution - 6 parameter values =\ 0
      IF(idistr.eq.8) THEN
        npartrue=npar
        OPEN(unit=11,file='par_gbl.txt')
        read(11,*) (par(j),j=1,npartrue)
        CLOSE(11)
        do j=1,npartrue
          if(par(j).eq.0.d0) then
            write(*,*) "Error: null parameter of GBL distribution"
            go to 9999
          end if
        end do
        do j=1,npartrue
          partrue=par(j)
          do jj=1,9
            dd2=dfloat(jj)
            gg1=partrue/(partrue+1.d0)
            gg2=(dd2+1.d0)**(partrue+1.d0)-dd2**(partrue+1.d0)
            checksigt(jj,j)=gg1*gg2/(10.d0**partrue-1.d0)
            gg1=(dd2+1.d0)**partrue-dd2**partrue
            gg2=10.d0**partrue-1.d0
            checkdigt(jj,j)=gg1/gg2
          end do
        end do
c if *idig1=1* 1st digit is Benford ==> *NO* for this distr.
        open(unit=73,file='powX2Q_summary_GBL.txt')
        open(unit=75,file='powX2DIF_summary_GBL.txt')
        open(unit=77,file='powX2DIF2T_summary_GBL.txt')
        open(unit=101, file='signif_GBL1.txt') 
        open(unit=102, file='signif_GBL2.txt') 
        open(unit=103, file='signif_GBL3.txt') 
        open(unit=104, file='signif_GBL4.txt') 
        open(unit=105, file='signif_GBL5.txt') 
        open(unit=106, file='signif_GBL6.txt')                  
      END IF
c idistr = 9: BenfGBL = First-Digit Benford + Generalized Benford - as above
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.9) THEN
        npartrue=npar
        OPEN(unit=11,file='par_benfgbl.txt')
        read(11,*) (par(j),j=1,npartrue)
        CLOSE(11)
        do j=1,npartrue
          if(par(j).eq.0.d0) then
            write(*,*) "Error: null parameter of GBL distribution"
            go to 9999
          end if
        end do
c if *idig1=1* 1st digit is Benford ==> *NO* for this distr.        
        open(unit=73,file='powX2Q_summary_BENFGBL.txt')
        open(unit=75,file='powX2DIF_summary_BENFGBL.txt')
        open(unit=77,file='powX2DIF2T_summary_BENFGBL.txt')
        open(unit=101, file='signif_BENFGBL1.txt') 
        open(unit=102, file='signif_BENFGBL2.txt') 
        open(unit=103, file='signif_BENFGBL3.txt') 
        open(unit=104, file='signif_BENFGBL4.txt') 
        open(unit=105, file='signif_BENFGBL5.txt') 
        open(unit=106, file='signif_BENFGBL6.txt')                  
      END IF
c idistr=10: 2-component Generalized Benford Mixture (3 parameter values with 
c different sign in the two components * 2 mixing proportions)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.10) THEN
        OPEN(unit=11,file='par_mixgbl.txt')
        npartrue=npar
        npar2=npartrue/2
        if(floor(float(npar2)).ne.npar2) then
          write(*,*) "npar not even: the selected distribution is ",
     +     "invalid"
          go to 9999
        end if
        read(11,*) (par(j),j=1,npar2)
        read(11,*) par(npar2+1),par(npar2+2)
        CLOSE(11)
c if *idig1=1* 1st digit is Benford ==> *NO* for this distr.        
        open(unit=73,file='powX2Q_summary_MIXGBL.txt')
        open(unit=75,file='powX2DIF_summary_MIXGBL.txt')
        open(unit=77,file='powX2DIF2T_summary_MIXGBL.txt')
        open(unit=101, file='signif_MIXGBL1.txt') 
        open(unit=102, file='signif_MIXGBL2.txt') 
        open(unit=103, file='signif_MIXGBL3.txt') 
        open(unit=104, file='signif_MIXGBL4.txt') 
        open(unit=105, file='signif_MIXGBL5.txt') 
        open(unit=106, file='signif_MIXGBL6.txt')                          
      END IF  
c idistr = 11: Dirac = First-Digit Dirac distribution (read npar values 
c                      of the fixed replacing digit ==> read as real*8)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.11) THEN
        npartrue=npar
        OPEN(unit=11,file='par_dirac.txt')
        read(11,*) (par(j),j=1,npartrue)
        CLOSE(11)
        do j=1,npartrue
          if(par(j).lt.1.d0.or.par(j).gt.9.d0) then
            write(*,*) "Error: invalid digit in Dirac distribution"
            go to 9999
          end if
        end do
c if *idig1=1* 1st digit is Benford ==> *NO* for this distr.        
        open(unit=73,file='powX2Q_summary_DIRAC.txt')
        open(unit=75,file='powX2DIF_summary_DIRAC.txt')
        open(unit=77,file='powX2DIF2T_summary_DIRAC.txt')
        open(unit=101, file='signif_DIRAC1.txt') 
        open(unit=102, file='signif_DIRAC2.txt') 
        open(unit=103, file='signif_DIRAC3.txt') 
        open(unit=104, file='signif_DIRAC4.txt') 
        open(unit=105, file='signif_DIRAC5.txt') 
        open(unit=106, file='signif_DIRAC6.txt')                  
      END IF        
c idistr = 12: Truncated Benford = Last-Digit(s) Dirac distribution 
c              (read 1 value of the fixed digit that replaces the last
c               j=1,...,npar digit values ==> read as real*8)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.12) THEN
        npartrue=npar
        OPEN(unit=11,file='par_trunc.txt')
        read(11,*) par(1)
        CLOSE(11)
        if(par(1).lt.0.d0.or.par(1).gt.9.d0) then
          write(*,*) "Error: invalid digit in Truncated distribution"
          go to 9999
        end if
c if *idig1=1* 1st digit is Benford ==> *NO* for this distr.        
        open(unit=73,file='powX2Q_summary_TRUNC.txt')
        open(unit=75,file='powX2DIF_summary_TRUNC.txt')
        open(unit=77,file='powX2DIF2T_summary_TRUNC.txt')
        open(unit=101, file='signif_TRUNC1.txt') 
        open(unit=102, file='signif_TRUNC2.txt') 
        open(unit=103, file='signif_TRUNC3.txt') 
        open(unit=104, file='signif_TRUNC4.txt') 
        open(unit=105, file='signif_TRUNC5.txt') 
        open(unit=106, file='signif_TRUNC6.txt')                          
      END IF    
c idistr =120: NEW Truncated Benford = Last-Digit(s) Dirac distribution 
c              with fixed digit 0 in a given proportion of cases
c              (read the given proportions of cases in which the last
c               npar=6 digit values are replaced by 0)
c par(1): proportion of cases where the fixed digit value 0 replaces the last 6 digits: 
c *only one significant digit remains in S(X)*
c par(2): proportion of cases where the fixed digit value 0 replaces the last 5 digits: 
c *only two significant digits remain in S(X)* 
c etc.
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.120) THEN
        npartrue=1
        OPEN(unit=11,file='par_truncNEW.txt')
        read(11,*) (par(j), j=1,npar)
        CLOSE(11)
        do j=1,npar
          ntrunc(j)=dnint(par(j)*fn)
          if(j.eq.1) then
            cumtrunc(j)=ntrunc(j)
          else
            cumtrunc(j)=cumtrunc(j-1)+ntrunc(j)
          end if
        end do
        if(cumtrunc(npar).gt.n) then
          write(*,*) "Error: invalid proportions in NEW Truncated ",
     +     "distribution"
          go to 9999
        end if
c if *idig1=1* 1st digit is Benford ==> *NO* for this distr.
        open(unit=73,file='powX2Q_summary_TRUNCNEW.txt')
        open(unit=75,file='powX2DIF_summary_TRUNCNEW.txt')
        open(unit=77,file='powX2DIF2T_summary_TRUNCNEW.txt')
        open(unit=101, file='signif_TRUNCNEW.txt') 
      END IF
c idistr=13: Cardiod distribution (3 values of the concentration parameter a
c            * 2 values of the location parameter b)
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.13) THEN
        OPEN(unit=11,file='par_card.txt')
        npartrue=npar
        npar2=npartrue/2
        if(floor(float(npar2)).ne.npar2) then
          write(*,*) "npar not even: the selected distribution is ",
     +     "invalid"
          go to 9999
        end if        
        read(11,*) (par(j),j=1,npar2)
        read(11,*) par(npar2+1),par(npar2+2)
        CLOSE(11)
        do j=1,npar2
          if(par(j).lt.(0.d0).or.par(j).gt.1.d0) then
            write(*,*) "Error: invalid concentration parameter a in ",
     +       "Cardioid distribution"
            go to 9999
          end if
        end do
        do j=(npar2+1),(npar2+2)
          if(par(j).lt.0.d0.or.par(j).gt.1.d0) then
            write(*,*) "Error: invalid location parameter b in ",
     +       "Cardioid distribution"
            go to 9999
          end if
          par(j)=par(j)*2.d0*pi
        end do
        do j=1,npar2
          xloc=par(npar2+1)
          partrue=par(j)
          checkmeant(j)=(partrue*dsin(xloc)+pi)/(2.d0*pi)
          checkvart(j)=(-3.d0*(partrue**2)*dsin(xloc)*dsin(xloc)+
     +     6.d0*partrue*dcos(xloc)+pi**2)/(12.d0*pi**2)
          checkvart(j)=dsqrt(checkvart(j))     
        end do
        do j=npar2+1,npartrue
          xloc=par(npar2+2)
          partrue=par(j-npar2)
          checkmeant(j)=(partrue*dsin(xloc)+pi)/(2.d0*pi)
          checkvart(j)=(-3.d0*(partrue**2)*dsin(xloc)*dsin(xloc)+
     +     6.d0*partrue*dcos(xloc)+pi**2)/(12.d0*pi**2)
          checkvart(j)=dsqrt(checkvart(j))     
        end do       
c if *idig1=1* 1st digit is Benford ==> different file names
c TO DO      
        open(unit=73,file='powX2Q_summary_CARD.txt')
        open(unit=75,file='powX2DIF_summary_CARD.txt')
        open(unit=77,file='powX2DIF2T_summary_CARD.txt')
        open(unit=101, file='signif_CARD1.txt') 
        open(unit=102, file='signif_CARD2.txt') 
        open(unit=103, file='signif_CARD3.txt') 
        open(unit=104, file='signif_CARD4.txt') 
        open(unit=105, file='signif_CARD5.txt') 
        open(unit=106, file='signif_CARD6.txt')                  
      END IF
C idistr=14: NNTS (Non Negative Trigonometric Sums: Fernandez-Duran, Biometrics, 2004)
C            npartrue=1 ==> read: *numNNTS <= npar - 1 = 5* 
C                                             (npar = max numNNTS + intercept) 
C                           read: *parNNTS = complex parameter vector* 
C                                  (intercept always real)  
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.14) THEN
        OPEN(unit=11,file='par_nnts.txt')
        npartrue=1
        read(11,*) numNNTS
        write(*,*) "numNNTS = ",numNNTS
        sumabsparNNTS=0.d0
        do j=1,numNNTS+1
          read(11,111) parNNTS(j)
111       format(f16.14,4x,f16.14) 
           absparNNTS(j)=cdabs(parNNTS(j))
           sumabsparNNTS=sumabsparNNTS+absparNNTS(j)*absparNNTS(j)
           write(*,*) "j = ",j,parNNTS(j),absparNNTS(j),sumabsparNNTS
        end do
        CLOSE(11)
        if(dimag(parNNTS(1)).ne.0.d0) then
          write(*,*) "Error: invalid intercept in NNTS distribution"
          go to 9999
        end if
c        write(*,*)
c        write(*,*) (parNNTS(j),j=1,numNNTS+1)
c        write(*,*) sumabsparNNTS,sngl(sumabsparNNTS)
c        pause 
!c Check the parameter constraint to single precision to avoid numerical approx
!        if(sngl(sumabsparNNTS).ne.1.0) then
!          write(*,*) "Error: invalid parameter constraint in NNTS ",
!     +     "distribution"
!          go to 9999
!        else
!          write(*,*) "Valid parameter constraint in NNTS distribution",
!     +     sngl(sumabsparNNTS)
!        end if       
c October 2022: Check the parameter constraint with tolerance to avoid 
c               numerical approx when parameter values are in single precision
        if(dabs(1.d0-sumabsparNNTS).gt.10.d-6) then
          write(*,*) "Error: invalid parameter constraint in NNTS ",
     +     "distribution"
          write(*,*) "sumabsparNNTS = ",sumabsparNNTS
          go to 9999
        else
          write(*,*) "Valid parameter constraint in NNTS distribution",
     +     sngl(sumabsparNNTS)
        end if         
c TO DO: compute nominal mean and variance for general numNNTS and parNNTS(j)
c          checkmeant(1)=
c          checkvart(1)=
c if *idig1=1* 1st digit is Benford ==> different file names
c TO DO
        open(unit=73,file='powX2Q_summary_NNTS.txt')
        open(unit=75,file='powX2DIF_summary_NNTS.txt')
        open(unit=77,file='powX2DIF2T_summary_NNTS.txt')
        open(unit=101, file='signif_NNTS.txt') 
      END IF
c NEW VERSION - OCTOBER 2022
c idistr=15: heteroschedastic 2-component Normal Mixture 
c (2 location parameters + 2 scale parameters) * 2 mixing proportions 
c TO DO: compute nominal significand + digits 
      IF(idistr.eq.15) THEN
        OPEN(unit=11,file='par_mixnorm.txt')
        npartrue=npar/3
        npar2=npar/3
        if(floor(float(npar2)).ne.npar2) then
          write(*,*) "npar not even: the selected distribution is ",
     +     "invalid"
          go to 9999
        end if
        read(11,*) (par(j),j=1,npar2)
        read(11,*) (par(j),j=npar2+1,npar2+2)
        read(11,*) par(npar2+3),par(npar2+4)
        CLOSE(11)
c if *idig1=1* 1st digit is Benford ==> different file names 
        if(idig1.eq.1) then       
          open(unit=73,file='powX2Q_summary_BENFMIXNORM.txt')
          open(unit=75,file='powX2DIF_summary_BENFMIXNORM.txt')
          open(unit=77,file='powX2DIF2T_summary_BENFMIXNORM.txt')
          open(unit=101, file='signif_BENFMIXNORM1.txt') 
          open(unit=102, file='signif_BENFMIXNORM2.txt') 
          open(unit=103, file='signif_BENFMIXNORM3.txt') 
          open(unit=104, file='signif_BENFMIXNORM4.txt') 
        else
          open(unit=73,file='powX2Q_summary_MIXNORM.txt')
          open(unit=75,file='powX2DIF_summary_MIXNORM.txt')
          open(unit=77,file='powX2DIF2T_summary_MIXNORM.txt')
          open(unit=101, file='signif_MIXNORM1.txt') 
          open(unit=102, file='signif_MIXNORM2.txt') 
          open(unit=103, file='signif_MIXNORM3.txt') 
          open(unit=104, file='signif_MIXNORM4.txt') 
        end if
      END IF 
c idistr =16: 3-component Normal Mixture: 3 fixed location parameters: par(1)-par(3)
c                                         3 fixed scale parameters: par(4)-par(6)
c                                         3 fixed mixing proportions in simulation: 1/3 
c TO DO: compute nominal significand + digits
      IF(idistr.eq.16) THEN
        OPEN(unit=11,file='par_3mixnorm.txt')
        npartrue=1
        npar2=npar/2
        if(floor(float(npar2)).ne.npar2) then
          write(*,*) "npar not even: the selected distribution is ",
     +     "invalid"
          go to 9999
        end if
        read(11,*) (par(j),j=1,npar2)
        read(11,*) (par(j),j=npar2+1,npar)
        CLOSE(11)
c if *idig1=1* 1st digit is Benford ==> different file names 
        if(idig1.eq.1) then       
          open(unit=73,file='powX2Q_summary_BENF3MIXNORM.txt')
          open(unit=75,file='powX2DIF_summary_BENF3MIXNORM.txt')
          open(unit=77,file='powX2DIF2T_summary_BENF3MIXNORM.txt')
          open(unit=101, file='signif_BENF3MIXNORM1.txt') 
        else
          open(unit=73,file='powX2Q_summary_3MIXNORM.txt')
          open(unit=75,file='powX2DIF_summary_3MIXNORM.txt')
          open(unit=77,file='powX2DIF2T_summary_3MIXNORM.txt')
          open(unit=101, file='signif_3MIXNORM1.txt') 
        end if
      END IF  
c      
c
c NEW OPTION: Open files for possibly saving test statistics with 
c             simulated data ==> *idistr >= 0*
c     save test statistics when isavetest=1 (old option: isavetest=isavesig)
c                                           new option: isavetest read from file
c     see files comp*.txt if *idistr < 0* (idistr.lt.(-1) and idistr.eq.(-1)) 
      if((idistr.ge.0).and.(isavetest.eq.1)) then
        OPEN(unit=410, file='simX2.txt')
        OPEN(unit=420, file='simX2_2d.txt')
        OPEN(unit=430, file='simQ.txt')
        OPEN(unit=440, file='simKS.txt')
        OPEN(unit=450, file='simX2DIF.txt')
        OPEN(unit=460, file='simX2DIF2T.txt')
        OPEN(unit=470, file='simSCORE1.txt')
        OPEN(unit=480, file='simSCORE2.txt')
        OPEN(unit=490, file='simSCOREBIC.txt')
        OPEN(unit=500, file='simXLRT1.txt')
        OPEN(unit=510, file='simXLRT2.txt')
        OPEN(unit=520, file='simXLRTBIC.txt')
        OPEN(unit=530, file='simVK.txt')
        OPEN(unit=540, file='simKSFR.txt')
        OPEN(unit=550, file='simVKFR.txt')
      end if
c    
cc possibly save simulated data and significands - NOT THE CURRENT OPTION
cc      open(unit=80,file='sim_simdata.txt')
cc      open(unit=81,file='sim_simsig.txt')
c
c      open(unit=73,file='powX2Q_summary.txt')
c      open(unit=75,file='powX2DIF_summary.txt')
c      open(unit=77,file='powX2DIF2T_summary.txt')
c
c Open files for saving test statistics and exact p-values
c                 if *idistr.lt.(-1)* (emprical analysis of multiple data sets)
c NEW OPTION: always save test statistics when *idistr < 0*
c             *idistr.lt.(-1)* (emprical analysis of multiple data sets)
c             *idistr.eq.(-1)* (emprical analysis of single data set)
c             isavetest not required
c For X2_2d + SCORE1 + SCORE2 + XLRT1 + XLRT2 the exact p-value is not (yet) computed
      if(idistr.lt.0) then
        OPEN(unit=51, file='compX2.txt')
        OPEN(unit=52, file='compX2_2d.txt')
        OPEN(unit=53, file='compQ.txt')
        OPEN(unit=54, file='compKS.txt')
        OPEN(unit=55, file='compX2DIF.txt')
        OPEN(unit=56, file='compX2DIF2T.txt')
        OPEN(unit=57, file='compSCORE1.txt')
        OPEN(unit=58, file='compSCORE2.txt')
        OPEN(unit=59, file='compSCOREBIC.txt')
        OPEN(unit=60, file='compXLRT1.txt')
        OPEN(unit=61, file='compXLRT2.txt')
        OPEN(unit=62, file='compXLRTBIC.txt')        
        OPEN(unit=63, file='compVK.txt')        
        OPEN(unit=64, file='compKSFR.txt')
        OPEN(unit=65, file='compVKFR.txt')        
      end if
c
C START LOOP OVER *j=1,npartrue* parameters 
C N.B: The loops over *npartrue* and *nsim2* could be exchanged ==> more 
C efficient code. However they are kept in this order to provide  
C easier control of the random simulation flow for each parameter value
      do j=1,npartrue
        checkX2DIF=0.d0
        checkX2DIF2T=0.d0
        checkstab(j)=0.d0
        checkmean(j)=0.d0
        checkvar(j)=0.d0
        do jj=1,9
          checksig(jj,j)=0.d0
          checkdig(jj,j)=0.d0
        end do  
        do jj=1,2
          X2Qmean(jj)=0.0
          X2DIFmean(jj)=0.0
          X2DIF2Tmean(jj)=0.0
          do jjj=1,2
            corX2Q(jjj,jj)=0.0
            corX2DIF(jjj,jj)=0.0
            corX2DIF2T(jjj,jj)=0.0
          end do
        end do
C  nsim2 > 1  if  *idistr.lt.(-1)*  
C  nsim2 = 1  if  *idistr.eq.(-1)*  
C  ==> save statistics and p-values in empirical analysis
        if(idistr.lt.0) then
          OPEN(unit=212,file='benf_datacheck.txt')
          OPEN(unit=11,file='benf_data.txt')
          OPEN(unit=12,file='benf_jitter.txt')        
        end if           
C START OF SIMULATION LOOP OVER *jj=1,nsim2* SIMULATIONS
        do jj=1,nsim2
          mm=mod(jj,loopp)
          if(mm.eq.0) write(*,*) "power sim # ",jj,"    parameter #",j
c
C GENERATE A SAMPLE OF *n* OBSERVATIONS FROM THE REQUIRED RANDOM VARIABLE 
C ACCORDING TO THE VALUE OF PARAMETER *idistr* (read from file)
C AND COMPUTE SIGNIFICAND AND DIGITS FOR EACH OBSERVATION
C USE IMSL ROUTINES FOR RANDOM NUMBER GENERATION
C N.B.: the first *n-ncont* observations in each sample are simulated from the
C       Benford random variable ==> contamination model with contamination
C       rate *rcont*
C       FULL POWER ANALYSIS IF *rcont=1* ==> *ncont=n*, *ngood=0*: *all*
C       the sample data are generated from the distribution defined 
C       by *idistr*
C       OUTLIER ANALYSIS FROM BENFORD DISTRIBUTION IF *0<rcont<1* and
C       *0<ncont<n* ==> *ngood>0* 
C       (rcont=0 is equivalent to rcont=1 + idistr=0)
c
c Start loop over sample elements
          do i=1,n
C BENFORD VALUE IN THE CONTAMINATION MODEL ==> if ncont<n and ngood>0
            if(i.le.ngood) then 
              rn=DRNUNF()
              if(rn.le.0.d0.or.rn.ge.1.d0) then
                write(*,*) "error in rn =",rn," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              simdata(i)=rn
              s(i)=10.d0**simdata(i)
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))
c Skip check on *idistr* and computation of s(i) ==> go to digit computation 
              goto 1001
            end if
c
c  if idistr < -1 read data within the loop over *nsim2* simulations
c  (if idistr = -1 data are read on input)
            if(idistr.lt.(-1)) then
              read(11,*) simdata(i)
              if(simdata(i).eq.0.d0) then
                write(*,*)
                write(*,*) "Warning: in simulation ",jj,
     +           "observation ",i,"   is 0 ==> NOT ALLOWED; ",simdata(i)
                write(*,*)
                go to 9999
              end if
c  Perform jittering in the data only if *ijit.gt.0*
c  Use the same jittering rule as in the empirical analysis
c                 of data (idistr.eq.(-1)) to avoid ties in DKSONE
              if(ijit.gt.0) then
                rng=DRNUNF()                           
                if(simdata(i).lt.10.d0**3) then
                  simdata(i)=simdata(i)+rng*(10.d0**(-6))   
                  write(12,11) i,simdata(i),rng,rng*(10.d0**(-6))
                else
                  simdata(i)=simdata(i)+rng*(10.d0**(-3))   
                  write(12,11) i,simdata(i),rng,rng*(10.d0**(-3))           
                end if
              end if    
            end if
C CONTAMINATED VALUE IN THE CONTAMINATION MODEL ACCORDING TO *idistr*
c idistr=0: Benford - No parameters
            if(idistr.eq.0) then 
              rn=DRNUNF()
              if(rn.le.0.d0.or.rn.ge.1.d0) then
                write(*,*) "error in rn =",rn," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              simdata(i)=rn
              s(i)=10.d0**simdata(i)
c if *idig1=1* 1st digit is Benford ==> compute the modified
c variable simdata(i) and the modified significand s(i) 
              if(idig1.eq.1) then 
                iss=floor(s(i))
                ss=s(i)-dfloat(iss) 
                rnb=DRNUNF()
                ssb=10.d0**rnb
                d1b=floor(ssb) 
                simdata(i)=ss+dfloat(d1b)                             
                ss=dlog10(dabs(simdata(i)))
                iss=floor(ss)
                ss=ss-dfloat(iss)
                s(i)=10.d0**ss              
c Now s(i) = simdata(i) ==> the previous 4 lines could be deleted: they  
c                           are kept for homogeneity with subsequent code
c                           for other non-Benford distributions
c          = fractional part of the original Benford number s(i) in (1,10)
c            + 1st digit d1b of the new Benford number ssb in (1,10)
              end if
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))              
c Skip the computation of s(i) for the other distributions
              goto 1001
            end if
c idistr=1: Normal distribution ==> npartrue = npar: 
c npar/2 scale parameter values (standard deviation) * 2 location parameter 
c values (npar must be an integer)
            if(idistr.eq.1) then
              if(j.le.npar2) then
                partrue=par(j)
                xloc=par(npar2+1)
              else
                partrue=par(j-npar2)
                xloc=par(npar2+2)
              end if
              rn=DRNNOF()
              simdata(i)=rn*partrue+xloc
            end if
c idistr=2: Lognormal distribution ==> npartrue = npar-1 shape parameter 
c values + (fixed) scale parameter 
c N.B.: in routine DRNLNL scal = Mean of the underlying normal distribution
c       ==> scal = log(par(npartrue+1)) if par(npartrue+1) = scale parameter
c                  of the Lognormal distribution
            if(idistr.eq.2) then
              partrue=par(j)
              scal=dlog(par(npartrue+1))              
              CALL DRNLNL(1,scal,partrue,rn)
              simdata(i)=rn
            end if
c idistr=3: Weibull distribution ==> npartrue = npar shape parameter values 
c the scale parameter is fixed and set = 1 (default of routine DRNWIB)
            if(idistr.eq.3) then
              partrue=par(j)
              CALL DRNWIB(1,partrue,rn)
              simdata(i)=rn
            end if
c idistr=4: Gamma distribution ==> npartrue = npar shape parameter values 
c the scale parameter is fixed and set = 1 (default of routine DRNGAM)
            if(idistr.eq.4) then
              partrue=par(j)
              CALL DRNGAM(1,partrue,rn)
              simdata(i)=rn
            end if
c idistr=5: Beta distribution ==> npartrue = npar: 
c npar/2 shape parameter values alpha * 2 shape parameter values beta
c (npar must be an integer)
            if(idistr.eq.5) then
              if(j.le.npar2) then
                partrue=par(j)
                xloc=par(npar2+1)
              else
                partrue=par(j-npar2)
                xloc=par(npar2+2)
              end if
              CALL DRNBET(1,partrue,xloc,rn)
              simdata(i)=rn
            end if
c idistr=6: Uniform distribution on (0,a) ==> npartrue = npar range parameter values 
            if(idistr.eq.6) then
              partrue=par(j)
              rn=DRNUNF()
              simdata(i)=rn*partrue
            end if     
c idistr=7: Stable distribution ==> npartrue = npar-1 shape parameter 
c values alpha + (fixed) skewness parameter beta
c N.B.: in routine DRNSTA scal = transformed parameter betaprime to be
c computed from beta=xloc and alpha=partrue ==> see Chambers et al. (JASA, 1976)
c If *beta=1* ==> *betaprime=1*   NEW 25/07/2021: always use *betaprime()*
c else        ==> use function *betaprime()* 
c Compute the empirical tail probability in vector checkstab(j) ==> n*nsim2
c simulations for each parameter value (j=1,...,npartrue)
            if(idistr.eq.7) then
              partrue=par(j)
              xloc=par(npartrue+1)
c always use *betaprime()*
              scal=betaprime(xloc,partrue)
c save betaprime
              betap(j)=scal 
c choice of the simulation algorithm for Stable distribution
c according to *istab*
              if(istab.eq.1) then
                CALL DRNSTA(1,partrue,scal,rn)
              end if
              if(istab.eq.2) then
                call positstab(partrue,rn)
              end if
              if(istab.eq.3) then
                call positstab2(partrue,rn)    
              end if 
              simdata(i)=rn
              if(dabs(simdata(i)).gt.tolstab(j)) 
     +         checkstab(j)=checkstab(j)+1.d0/(fn*fn2)
            end if
c idistr=8: Generalized Benford distribution ==> npartrue = npar
            if(idistr.eq.8) then 
              partrue=par(j)
              rn=DRNUNF()
              if(rn.le.0.d0.or.rn.ge.1.d0) then
                write(*,*) "error in rn =",rn," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              simdata(i)=rn*(10.d0**partrue-1.d0)
              s(i)=(1.d0+simdata(i))**(1/partrue)
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))
c Skip the computation of s(i) for the other distributions
              goto 1001
            end if
c idistr=9: BenfGBL distribution: First-Digit Benford + Generalized 
c           Benford ==> npartrue = npar
            if(idistr.eq.9) then 
              partrue=par(j)
              rnb=DRNUNF()
              if(rnb.le.0.d0.or.rnb.ge.1.d0) then
                write(*,*) "error in rnb =",rnb," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              ssb=10.d0**rnb
              d1b=floor(ssb)        
              if(d1b.lt.1.or.d1b.gt.9) then
                write(*,*) "error in d1b at unit #",i," in pow sim #",jj
                pause
              end if
              partrue=par(j)
              rng=DRNUNF()
              if(rng.le.0.d0.or.rng.ge.1.d0) then
                write(*,*) "error in rng =",rng," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              simdata(i)=rng*(10.d0**partrue-1.d0)
              ssg=(1.d0+simdata(i))**(1/partrue)
              d1g=floor(ssg)        
              if(d1g.lt.1.or.d1g.gt.9) then
                write(*,*) "error in d1g at unit #",i," in pow sim #",jj
                pause
              end if          
              d1dif=abs(d1g-d1b)
              if(ssg.le.ssb) then
                s(i)=ssg+dfloat(d1dif)
              else
                s(i)=ssg-dfloat(d1dif)
              end if
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))
c Skip the computation of s(i) for the other distributions
              goto 1001
            end if
c idistr=10: 2-component Generalized Benford mixture distribution 
c            ==> npartrue = npar: npar/2 values of the GBL parameter (with 
c            opposite sign in the two components) * 2 values of the mixing 
c            proportion (npar must be an integer)
            if(idistr.eq.10) then
              if(j.le.npar2) then
                partrue=par(j)
                xloc=par(npar2+1)
              else
                partrue=par(j-npar2)
                xloc=par(npar2+2)
              end if
              rn1=DRNUNF()
              if(rn1.le.0.d0.or.rn1.ge.1.d0) then
                write(*,*) "error in rn1 =",rn1," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              simdata(i)=rn1*(10.d0**partrue-1.d0)
              ss1=(1.d0+simdata(i))**(1/partrue)
              rn2=DRNUNF()
              if(rn2.le.0.d0.or.rn2.ge.1.d0) then
                write(*,*) "error in rn2 =",rn2," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              simdata(i)=rn2*(10.d0**(-1.d0*partrue)-1.d0)
              ss2=(1.d0+simdata(i))**(1/(-1.d0*partrue))
              rn=DRNUNF()
              if(rn.le.0.d0.or.rn.ge.1.d0) then
                write(*,*) "error in rn =",rn," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if
              if(rn.le.xloc) then
                s(i)=ss1
              else
                s(i)=ss2
              end if  
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))              
c Skip the computation of s(i) for the other distributions
              goto 1001
            end if
c idistr=11: Dirac distribution for the first digit
            if(idistr.eq.11) then 
              rnb=DRNUNF()
              if(rnb.le.0.d0.or.rnb.ge.1.d0) then
                write(*,*) "error in rnb =",rnb," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if            
              ssb=10.d0**rnb
              d1b=floor(ssb)        
              if(d1b.lt.1.or.d1b.gt.9) then
                write(*,*) "error in d1b at unit #",i," in pow sim #",jj
                pause
              end if
              partrue=par(j)
              ipartrue=floor(partrue)
              d1dif=abs(ipartrue-d1b)
              if(partrue.gt.ssb) then
                s(i)=ssb+dfloat(d1dif)
              else
                s(i)=ssb-dfloat(d1dif)
              end if
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))
c Skip the computation of s(i) for the other distributions
              goto 1001
            end if          
c idistr=12: Truncated Benford distribution = Dirac distribution for the 
c last j=1,...,npar digit(s) ==> the significand is truncated at the *7th* 
c decimal place (8 digits) and the *last* npar=6 digits are set to the
c fixed digit value par(1)   ==> if j=1: 6 decimal digits are retained, 
c                                the others are set to par(1);
c                                if j=6: 1 decimal digit is retained, 
c                                the others are set to par(1);
c                            ==> the significand is then *jittered* in the
c                                10th decimal place because of the continuity 
c                                requirement of DKSONE
c more jittering in large samples (n>200)
            if(idistr.eq.12) then 
              rnb=DRNUNF()
              if(rnb.le.0.d0.or.rnb.ge.1.d0) then
                write(*,*) "error in rnb =",rnb," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if  
              partrue=par(1)
              ipartrue=floor(partrue)          
              ssb=10.d0**rnb
              ss=ssb*(10.d0**(6-j+1))              
              iss8=floor(ss)
              ss=(dfloat(iss8))*(10.d0**(-6+j-1))                         
              do ij=1,j
                ss=ss+partrue*10.d0**(-7+ij-1)
              end do
              rng=DRNUNF()
c more jittering also when n=200
              if(n.lt.200) then
                s(i)=ss+rng*(10.d0**(-9))
              else
                s(i)=ss+rng*(10.d0**(-7))
              end if     
              if(s(i).ge.10.d0) then 
                write(*,*) "Error in jittering:",j,jj,i,s(i)
                go to 9999
              end if
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) then
c more jittering also when n=200
                if(n.lt.200) then                
                  write(jout,302) ss,rng*(10.d0**(-9)),
     +             s(i),dlog10(s(i))
                else
                  write(jout,302) ss,rng*(10.d0**(-7)),
     +             s(i),dlog10(s(i))
                end if                
              end if         
c Skip the computation of s(i) for the other distributions
              goto 1001
            end if  
c idistr =120: NEW Truncated Benford = Last-Digit(s) Dirac distribution 
c              with fixed digit 0 in a given proportion of cases
c              (read the given proportions of cases in which the last
c               npar=6 digit values are replaced by 0)
c * The Benford significand is first truncated at the *6th* decimal place 
c   (7 digits) ==> *j=1* in ss below (npartrue=1)
c * par(1): proportion of cases where the fixed digit value 0 replaces the last 6 digits: 
c   *only one significant digit remains in S(X)*
c * par(2): proportion of cases where the fixed digit value 0 replaces the last 5 digits: 
c   *only two significant digits remain in S(X)* 
c etc.
            if(idistr.eq.120) then 
              partrue=0.d0
              rnb=DRNUNF()
              if(rnb.le.0.d0.or.rnb.ge.1.d0) then
                write(*,*) "error in rnb =",rnb," at unit #",i,
     +           " in pow sim #",jj         
                pause
              end if  
              ssb=10.d0**rnb
              ss=ssb*(10.d0**(6-j+1))              
              iss8=floor(ss)
              ss=(dfloat(iss8))*(10.d0**(-6+j-1))
              if(i.le.cumtrunc(1)) then                                     
                ij=6
                iss8=floor(ss)
                ss=(dfloat(iss8))+partrue*10.d0**(-7+ij-1)
              else if(i.gt.cumtrunc(1).and.i.le.cumtrunc(2)) then
                ij=5                                     
                iss8=floor(ss*(10.d0**1))
                ss=(dfloat(iss8*(10.d0**(-1))))+
     +           partrue*10.d0**(-7+ij-1)
              else if(i.gt.cumtrunc(2).and.i.le.cumtrunc(3)) then
                ij=4                                     
                iss8=floor(ss*(10.d0**2))
                ss=(dfloat(iss8*(10.d0**(-2))))+
     +           partrue*10.d0**(-7+ij-1)
              else if(i.gt.cumtrunc(3).and.i.le.cumtrunc(4)) then
                ij=3                                     
                iss8=floor(ss*(10.d0**3))
                ss=(dfloat(iss8*(10.d0**(-3))))+
     +           partrue*10.d0**(-7+ij-1)
              else if(i.gt.cumtrunc(4).and.i.le.cumtrunc(5)) then
                ij=2                                     
                iss8=floor(ss*(10.d0**4))
                ss=(dfloat(iss8*(10.d0**(-4))))+
     +           partrue*10.d0**(-7+ij-1)
              else if(i.gt.cumtrunc(5).and.i.le.cumtrunc(6)) then
                ij=1                                     
                iss8=floor(ss*(10.d0**5))
                ss=(dfloat(iss8*(10.d0**(-5))))+
     +           partrue*10.d0**(-7+ij-1)
              end if             
              rng=DRNUNF()
c more jittering also when n=200
              if(n.lt.200) then              
                s(i)=ss+rng*(10.d0**(-9))
c even more jittering for big n and idistr.eq.120           
              else if(n.lt.500) then 
                s(i)=ss+rng*(10.d0**(-7))
              else                
!                s(i)=ss+rng*(10.d0**(-6))
                s(i)=ss+rng*(10.d0**(-5))
              end if     
              if(s(i).ge.10.d0) then 
                write(*,*) "Error in jittering:",j,jj,i,s(i)
                go to 9999
              end if
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) then
c more jittering also when n=200
                if(n.lt.200) then
                  write(jout,302) ss,rng*(10.d0**(-9)),
     +             s(i),dlog10(s(i))
c even more jittering for big n and idistr.eq.120
                else if(n.lt.500) then
                  write(jout,302) ss,rng*(10.d0**(-7)),
     +             s(i),dlog10(s(i))
                else
!                  write(jout,302) ss,rng*(10.d0**(-6)),
!     +             s(i),dlog10(s(i))                
                  write(jout,302) ss,rng*(10.d0**(-5)),
     +             s(i),dlog10(s(i))                
                end if                
              end if         
c Skip the computation of s(i) for the other distributions
              goto 1001
            end if  
c            
c idistr=13: Cardioid distribution ==> npartrue = npar 
c npar/2 concentration parameter values a * 2 location parameter values b
c (npar must be an integer)
            if(idistr.eq.13) then
              if(j.le.npar2) then
                partrue=par(j)
                xloc=par(npar2+1)
              else
                partrue=par(j-npar2)
                xloc=par(npar2+2)
              end if
1010          rn1=DRNUNF()
              rn2=DRNUNF()
              xx1=1.d0+partrue*dcos(2.d0*pi*rn1+xloc)
              xx2=(1.d0+partrue)*rn2
              if(xx1.lt.xx2) then
                goto 1010
              end if
              rn=rn1
              simdata(i)=rn
              s(i)=10.d0**simdata(i)
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))
c Skip the computation of s(i) for the other distributions
              goto 1001              
            end if            
c idistr=14: NNTS (Non Negative Trigonometric Sums: Fernandez-Duran, Biometrics, 2004)
c            ==> npartrue=1:  numNNTS trigonometric sums (complex parameters) 
c                             + real intercept  
            if(idistr.eq.14) then
              sumabsparNNTS=0.d0
              do jjj=1,numNNTS+1
                do jjjj=1,numNNTS+1
                  sumabsparNNTS=sumabsparNNTS+absparNNTS(jjj)*
     +             absparNNTS(jjjj)                
                end do
              end do
1011          rn1=DRNUNF()
              rn2=DRNUNF()
              sumNNTS=(0.d0,0.d0)
              do jjj=1,numNNTS+1
c exp(ix): computational alternative 1 (straightforward)
                eNNTS=2.d0*pi*dfloat(jjj-1)*rn1
                expNNTS=dcmplx(0.d0,eNNTS)
                sumNNTS=sumNNTS+parNNTS(jjj)*cdexp(expNNTS)          
c exp(ix): computational alternative 2 (using Euler's formula)
c                eNNTS=2.d0*pi*dfloat(jjj-1)*rn1
c                expNNTS=dcmplx(dcos(eNNTS),dsin(eNNTS))
cc if exp(-ix) is needed
cc                expNNTS=DCONJG(expNNTS)                              
c                sumNNTS=sumNNTS+parNNTS(jjj)*expNNTS                           
              end do
              abssumNNTS=(cdabs(sumNNTS))**2
              prodNNTS=sumabsparNNTS*rn2
              if(abssumNNTS.lt.prodNNTS) then
                goto 1011
              end if        
              rn=rn1
              simdata(i)=rn
              s(i)=10.d0**simdata(i)
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
              jout=100+j
              if(isavesig.eq.1) write(jout,300) dlog10(s(i))
c Skip the computation of s(i) for the other distributions
              goto 1001              
            end if          
c            
c NEW VERSION - OCTOBER 2022
c idistr=15: Heteroschedastic Normal Mixture distribution ==> npartrue = 2: 
c (2 location parameters + 2 scale parameters) * 2 mixing proportions 
c par(1): location parameter of the first mixture component
c par(2): location parameter of the second mixture component
c par(3): scale parameter of the first mixture component
c par(4): scale parameter of the second mixture component
c par(5): first value of the mixing proportion for *component 1*
c par(6): second value of the mixing proportion for *component 1*
c npar2 = npar/3 = 2
            if(idistr.eq.15) then
                partrue=par(npar2+2+j)
                xloc=par(1)
                shift=par(2)
                scal1=par(3)
                scal2=par(4)
              rn1=DRNNOF()
              ss1=rn1*scal1+xloc
              rn2=DRNNOF()
              ss2=rn2*scal2+shift
              rn=DRNUNF()
              if(rn.le.partrue) then
                simdata(i)=ss1
              else
                simdata(i)=ss2
              end if
            end if 
c            
c idistr=16: 3-component Normal Mixture distribution ==> npartrue = 1: 
c 3 fixed location parameters (xloc): par(1)-par(3)
c 3 fixed scale (standard deviation) parameters (scal): par(4)-par(6)
c 3 fixed mixing proportions (partrue): 1/3
            if(idistr.eq.16) then
              partrue=1.d0/3.d0
              rn=DRNUNF()
              rn1=DRNNOF()
              if(rn.le.partrue) then
                xloc=par(1)
                scal=par(4)                          
              else if(rn.le.(2.d0*partrue)) then
                xloc=par(2)
                scal=par(5)                          
              else 
                xloc=par(3)
                scal=par(6)                          
              end if            
              simdata(i)=rn1*scal+xloc            
            end if
c
c
c Compute significand s(i)
            ss=dlog10(dabs(simdata(i)))
            iss=floor(ss)
            ss=ss-dfloat(iss)
            s(i)=10.d0**ss            
c if *idig1=1* 1st digit is Benford ==> compute the modified
c variable simdata(i) and the modified significand s(i)             
            if(idig1.eq.1) then           
              iss=floor(s(i))
              ss=s(i)-dfloat(iss)              
              rnb=DRNUNF()
              ssb=10.d0**rnb
              d1b=floor(ssb)        
              simdata(i)=ss+dfloat(d1b)                      
              ss=dlog10(dabs(simdata(i)))
              iss=floor(ss)
              ss=ss-dfloat(iss)
              s(i)=10.d0**ss
            end if
c SAVE SIGNIFICAND s(i) in (0,1): fractional part of log10(|x|)
c (if dummy isavesig = 1)
c ALSO SAVE DATA IN THE SAME FILE
            jout=100+j
            if(isavesig.eq.1) write(jout,301) dlog10(s(i)),simdata(i)            
1001        continue
c COMPUTE sFR(i) = <s(i)> = significand *fractional part*: to be used 
c                           in tests KSFR and VKFR
            sFR(i)=s(i)-floor(s(i))
c
c Check empirical mean and standard deviation of simulated data 
c (exclude the Stable distribution)
            if(idistr.ne.7) then
              checkmean(j)=checkmean(j)+simdata(i)/(fn*fn2)
              checkvar(j)=checkvar(j)+simdata(i)*simdata(i)/(fn*fn2)
            end if
c Compute digits
            d1(i)=floor(s(i))        
            if(d1(i).lt.1.or.d1(i).gt.9) then
              write(*,*) "error in d1 at unit #",i," in pow sim #",jj
              pause
            end if
            ss=10.d0*s(i)
            d2(i)=floor(ss)-10.d0*floor(s(i))
            if(d2(i).lt.0.or.d2(i).gt.9) then
              write(*,*) "error in d2 at unit #",i," in pow sim #",jj
              pause
            end if
            d12(i,1)=d1(i)
            d12(i,2)=d2(i)
c Check mean significand and first digit probability
            checksig(d1(i),j)=checksig(d1(i),j)+s(i)/(fn*fn2)
            checkdig(d1(i),j)=checkdig(d1(i),j)+1.d0/(fn*fn2)
c save data, significands and digits for empirical data 
c (idistr < 0; npartrue=1) ==> check data
            if(idistr.lt.0) then
              write(212,12) i,simdata(i),s(i),d12(i,1),d12(i,2)
12            format(i9,f36.16,4x,f18.16,i4,i4)  
            end if
c
cc NO: *jittering for all the significands*
cc Sort the significands if idistr=12 to avoid ties in DKSONE ==> if the difference
cc between two significands is < tol, then jitter one of them at the 9th 
cc decimal place
c          if(idistr.eq.12) then
c            call DSVRGN(n,s,s)
c            do i=2,n
c              xdif=s(i)-s(i-1)
c              if(xdif.lt.tol) then
c                rng=DRNUNF()
c                s(i-1)=s(i-1)+rng*(10.d0**(-9))
c                if(s(i-1).ge.10.d0) then 
c                  write(*,*) "Error in jittering",s(i-1)
c                  go to 9999
c                end if
c              end if
c            end do
c          end if
c       
c End loop over sample elements
          end do
c possibly save data and significands
c          write(80,80) (simdata(i), i=1,n),j
c          write(81,80) (s(i), i=1,n),j
c80        format(f15.7,8x,"parameter",i5) 
c          if(idistr.lt.0) then
c            CLOSE(12)
c          end if
C
C COMPUTE TEST STATISTICS ==> NEW OPTION: possibly save to file if 
C                             isavetest=1 (old option isavetest=isavesig)
C                                         new option: isavetest read from file
!40        format(f12.6,8x,"parameter",i5)
!41        format(f12.6,4x,"ncomp",i4,8x,"parameter",i5)
!44        format(6f12.6,8x,"parameter",i5)
c reduce output to save space in files
40        format(f12.6,6x,i3)
41        format(f12.6,4x,"ncomp",i4,6x,i3)
44        format(6f12.6,6x,i3)
C COMPUTE (exact + asymptotic) POWER OF TESTS
          jsimul=jj
c
c Pearson X2 on first digit
          call X2ONE(jsimul,n,maxn,d1,benf1n,X2)
          X2Qmean(1)=X2Qmean(1)+X2/fn2
          corX2Q(1,1)=corX2Q(1,1)+X2*X2/fn2
          X2DIFmean(1)=X2DIFmean(1)+X2/fn2
          corX2DIF(1,1)=corX2DIF(1,1)+X2*X2/fn2
          X2DIF2Tmean(1)=X2DIFmean(1)
          corX2DIF2T(1,1)=corX2DIF(1,1)
c
c Pearson X2 on first-two digits
          call X2TWO(jsimul,n,maxn,d12,benf12n,X2_2d)
c
c Hotelling Sum-invariance test Q
          call HOTSUM(jsimul,n,maxn,d1,s,benfE,benfVARinv,Q)
          X2Qmean(2)=X2Qmean(2)+Q/fn2
          corX2Q(2,2)=corX2Q(2,2)+Q*Q/fn2
          corX2Q(1,2)=corX2Q(1,2)+X2*Q/fn2
c
c Kolmogorov-Smirnov test on the significand S(X)
c j=1,...,6 output values in vector KS(j) from IMSL subroutine
          call DKSONE(BENFCDF,n,s,KS,nmiss)
c
c X2DIF = Q - X2
          X2DIF = Q - X2
          if(X2DIF.lt.0.d0) checkX2DIF = checkX2DIF + 1.d0/fn2
          X2DIFmean(2)=X2DIFmean(2)+X2DIF/fn2
          corX2DIF(2,2)=corX2DIF(2,2)+X2DIF*X2DIF/fn2
          corX2DIF(1,2)=corX2DIF(1,2)+X2*X2DIF/fn2
c
c X2DIF2T = |Q - X2|
          X2DIF2T = dabs(Q - X2)
          if(X2DIF2T.lt.0.d0) checkX2DIF2T = checkX2DIF2T + 1.d0/fn2
          X2DIF2Tmean(2)=X2DIF2Tmean(2)+X2DIF2T/fn2
          corX2DIF2T(2,2)=corX2DIF2T(2,2)+X2DIF2T*X2DIF2T/fn2
          corX2DIF2T(1,2)=corX2DIF2T(1,2)+X2*X2DIF2T/fn2
c
c SCORE + LR TESTS ==> USE SUBROUTINE TEST() WITH 
c                      *ITERATIVE* BIC SELECTION ALGORITHM
        CALL TEST(jsimul,n,maxn,ncompmax,s,SCORE1,SCORE2,SCOREBIC,
     +    ncomps,XLRT1,XLRT2,XLRTBIC,ncompx)
c        CALL TEST2(jsimul,n,maxn,ncompmax,s,SCORE1,SCORE2,SCOREBIC,
c     +    ncomps,XLRT1,XLRT2,XLRTBIC,ncompx)
          avencompscores(j)=avencompscores(j)+dfloat(ncomps)/fn2
          varncompscores(j)=varncompscores(j)+dfloat(ncomps*ncomps)/fn2
          if(ncomps.eq.ncompmax) maxncompscores(j)=maxncompscores(j)+
     +     1.d0/fn2
          avencompxlrts(j)=avencompxlrts(j)+dfloat(ncompx)/fn2
          varncompxlrts(j)=varncompxlrts(j)+dfloat(ncompx*ncompx)/fn2
          if(ncompx.eq.ncompmax) maxncompxlrts(j)=maxncompxlrts(j)+
     +     1.d0/fn2
c
c Kuiper test VK (based on output of DKSONE) + asymptotic p-value
          VK=KS(2)+KS(3)
          CALL PVK(VK,pvalVKas,n)          
c
c Kolmogorov-Smirnov test on the significand *fractional part*
c j=1,...,6 output values in vector KSFR(j) from IMSL subroutine
          call DKSONE(BENFCDFFR,n,sFR,KSFR,nmiss)
c
c Kuiper test VKFR (based on output of DKSONE) + asymptotic p-value
          VKFR=KSFR(2)+KSFR(3)
          CALL PVK(VKFR,pvalVKFRas,n)        
c
c NEW OPTION with simulated data (*idistr > 0*): save test statistics if
c            isavetest=1 (isavetest read from file)
c see files comp*.txt if idistr.lt.0 (NEW OPTION for saving test statistics +
c            exact p-values in empirical analysis)
          if((idistr.ge.0).and.(isavetest.eq.1)) then
            write(410,40) X2,j
            write(420,40) X2_2d,j
            write(430,40) Q,j
            write(440,44) (KS(jjjj),jjjj=1,6),j
            write(450,40) X2DIF,j
            write(460,40) X2DIF2T,j          
            write(470,40) SCORE1,j
            write(480,40) SCORE2,j
            write(490,41) SCOREBIC,ncomps,j
            write(500,40) XLRT1,j
            write(510,40) XLRT2,j
            write(520,41) XLRTBIC,ncompx,j
            write(530,40) VK,j
            write(540,44) (KSFR(jjjj),jjjj=1,6),j
            write(550,40) VKFR,j
          end if
c
c
c skip the loop for computing exact p-values of combined 
c test statistics if dummy *itest=1* ==> perform 
c individual tests only (but for simplicity retain 
c subsequent computations)
          if(itest.eq.1) goto 11000
c
c Exact *combined* tests (min-pvalue principle): need to read *SORTED*
c null test statistics from file for exact p-value computation ==> statistics
c are sorted in *DESCENDING* order to reduce computing time under (distant)
c alternatives (this read step could be avoided if more memory is available 
c and further verctors can be stored)  
c An *approximate* but *faster* procedure under close alternatives can
c also be adopted (according to the value of *iapprox*):
c * compute the individual p-values if at least one of the statistic 
c   to be combined is > perc(1,*)
c * set all individual p-values = 1 otherwise: the individual p-values are larger
c   than testsize(1) ==> do not obtain the full distribution of p-values 
c                        under the alternative
c                    ==> power is available only for test sizes <= testsize(1)
c                    ==> in general the approximation is expected to be good if the 
c                        the alternative is far from Benford (but the reduction
c                        in computing time will be limited) 
c                    ==> the approximation is expected to be *exact* if the 
c                        exact percentile of p-values of the combined test 
c                        is < testsize(1): see the check value *icomb* for each 
c                        combined test
c
c * N.B.: the speed of the approximate procedure is greatly reduced when
c         combined tests based on X2DIF and the significand fractional part
c         (KSFR + VKFR) are computed (February 2022 + April 2022) ==> they
c         are expected to be negatively correlated with the other tests 
c 
c * More efficient algorithms could be adopted: efficient selection algorithm
c                                               or computational trick by
c                                               YYYY - ** TO DO **
c
c For a single test statistic a *faster* procedure could possibly be adopted:
c * simulate all the nsim2 random numbers under the alternative at once
c * compute the corresponding nsim2 test statistics
c * sort the nsim2 values of the test statistics
c * compute the exact p-value for each *sorted* test statistic under the alternative 
c   by comparing it with the *sorted* test statistics under the null ==> the
c   nsim1 loop for p-value computation will reduce to one or few steps 
c 
c Exact p-values of X2_2d + SCORE1 + SCORE2 + XLRT1 + XLRT2 are not (yet)
c computed ==> not (yet) used in combined tests - TO DO
c
          if(iapprox.eq.1) then
            ipval=0
            if(X2.gt.pind(1,1)) ipval=1
            if(Q.gt.pind(1,3)) ipval=1
            if(KS(1).gt.pind(1,4)) ipval=1
            if(X2DIF.gt.pind(1,5)) ipval=1
            if(X2DIF2T.gt.pind(1,6)) ipval=1
            if(SCOREBIC.gt.pind(1,9)) ipval=1
            if(XLRTBIC.gt.pind(1,12)) ipval=1
            if(VK.gt.pind(1,13)) ipval=1   
            if(KSFR(1).gt.pind(1,14)) ipval=1
            if(VKFR.gt.pind(1,15)) ipval=1
            if(ipval.eq.0) then
              pvalX2=1.0
              pvalQ=1.0
              pvalKS=1.0
              pvalX2DIF=1.0
              pvalX2DIF2T=1.0
              pvalSCOREBIC=1.0
              pvalXLRTBIC=1.0
              pvalVK=1.0
              pvalKSFR=1.0
              pvalVKFR=1.0
              go to 1002
            end if
          end if  
          OPEN(unit=131,file='benfX2sort.txt')
          OPEN(unit=133,file='benfQsort.txt')
          OPEN(unit=134,file='benfKSsort.txt')
          OPEN(unit=135,file='benfX2DIFsort.txt') 
          OPEN(unit=136,file='benfX2DIF2Tsort.txt') 
          OPEN(unit=139,file='benfSCOREBICsort.txt') 
          OPEN(unit=142,file='benfXLRTBICsort.txt') 
          OPEN(unit=143,file='benfVKsort.txt') 
          OPEN(unit=144,file='benfKSFRsort.txt')
          OPEN(unit=145,file='benfVKFRsort.txt')
          pvalX2=0.0
          pvalQ=0.0
          pvalKS=0.0
          pvalX2DIF=0.0
          pvalX2DIF2T=0.0
          pvalSCOREBIC=0.0
          pvalXLRTBIC=0.0
          pvalVK=0.0          
          pvalKSFR=0.0          
          pvalVKFR=0.0          
          do ij=1,nsim1
            istopX2=1
            istopQ=1
            istopKS=1
            istopX2DIF=1
            istopX2DIF2T=1
            istopSCOREBIC=1          
            istopXLRTBIC=1          
            istopKSFR=1          
            istopVKFR=1          
            read(131,30) X2ss
            if(dfloat(X2ss).gt.X2) then
              pvalX2=pvalX2+1.0/float(nsim1)
              istopX2=0
            end if
            read(133,30) Qss
            if(dfloat(Qss).gt.Q) then 
              pvalQ=pvalQ+1.0/float(nsim1)
              istopQ=0
            end if
            read(134,30) KSss
            if(dfloat(KSss).gt.KS(1)) then
              pvalKS=pvalKS+1.0/float(nsim1)
              istopKS=0
            end if
            read(135,30) X2DIFss
            if(dfloat(X2DIFss).gt.X2DIF) then
              pvalX2DIF=pvalX2DIF+1.0/float(nsim1)
              istopX2DIF=0
            end if
            read(136,30) X2DIF2Tss
            if(dfloat(X2DIF2Tss).gt.X2DIF2T) then
              pvalX2DIF2T=pvalX2DIF2T+1.0/float(nsim1)
              istopX2DIF2T=0
            end if            
            read(139,30) SCOREBICss
            if(dfloat(SCOREBICss).gt.SCOREBIC) then
              pvalSCOREBIC=pvalSCOREBIC+1.0/float(nsim1)
              istopSCOREBIC=0
            end if
            read(142,30) XLRTBICss
            if(dfloat(XLRTBICss).gt.XLRTBIC) then
              pvalXLRTBIC=pvalXLRTBIC+1.0/float(nsim1)
              istopXLRTBIC=0
            end if
            read(143,30) VKss
            if(dfloat(VKss).gt.VK) then
              pvalVK=pvalVK+1.0/float(nsim1)
              istopVK=0
            end if          
            read(144,30) KSFRss
            if(dfloat(KSFRss).gt.KSFR(1)) then
              pvalKSFR=pvalKSFR+1.0/float(nsim1)
              istopKSFR=0
            end if          
            read(145,30) VKFRss
            if(dfloat(VKFRss).gt.VKFR) then
              pvalVKFR=pvalVKFR+1.0/float(nsim1)
              istopVKFR=0
            end if          
            istop=istopX2*istopQ*istopKS*istopX2DIF*istopX2DIF2T
            istop=istop*istopSCOREBIC*istopXLRTBIC*istopVK
            istop=istop*istopKSFR*istopVKFR
            if(istop.eq.1) goto 1002
          end do
1002      continue
          CLOSE(131)
          CLOSE(133)
          CLOSE(134)
          CLOSE(135)
          CLOSE(136)
          CLOSE(139)
          CLOSE(142)
          CLOSE(143)
          CLOSE(144)
          CLOSE(145)
c Skip computation of all combined test statistics
c if dummy *itest=1* 
11000     continue
c
C COMPUTE POWER (exact + asymptotic) ==> *jjj=1,nqprop* significance levels
C for each test ==> See intial comments for the list of avaialable tests
c Exact *individual* tests ==> *nqprop* significance levels
c Asymptotic *individual* tests ==> *nqprop* significance levels
c Exact *combined* tests ==> based on percentiles of *p-values* instead of 
c test statistics (*single precision* is used)
c * NO asymptotic combined tests *
          do jjj=1,nqprop
c Individual tests          
            if(X2.gt.pind(jjj,1)) powind(j,jjj,1)=powind(j,jjj,1)+
     +       1.d0/fn2
            if(X2.gt.percas(jjj,1)) powas(j,jjj,1)=powas(j,jjj,1)+
     +       1.d0/fn2
            if(X2_2d.gt.pind(jjj,2)) powind(j,jjj,2)=powind(j,jjj,2)+
     +       1.d0/fn2
            if(X2_2d.gt.percas(jjj,2)) powas(j,jjj,2)=powas(j,jjj,2)+
     +       1.d0/fn2
            if(Q.gt.pind(jjj,3)) powind(j,jjj,3)=powind(j,jjj,3)+
     +       1.d0/fn2
            if(Q.gt.percas(jjj,3)) powas(j,jjj,3)=powas(j,jjj,3)+
     +       1.d0/fn2
            if(KS(1).gt.pind(jjj,4)) powind(j,jjj,4)=powind(j,jjj,4)+
     +       1.d0/fn2
            if(KS(6).lt.testsize2(jjj)) powas(j,jjj,4)=powas(j,jjj,4)+
     +       1.d0/fn2            
            if(X2DIF.gt.pind(jjj,5)) powind(j,jjj,5)=powind(j,jjj,5)+
     +       1.d0/fn2
            if(X2DIF.gt.percas(jjj,5)) powas(j,jjj,5)=powas(j,jjj,5)+
     +       1.d0/fn2
            if(X2DIF2T.gt.pind(jjj,6)) powind(j,jjj,6)=powind(j,jjj,6)+
     +       1.d0/fn2
            if(X2DIF2T.gt.percas(jjj,6)) powas(j,jjj,6)=powas(j,jjj,6)+
     +       1.d0/fn2
            if(SCORE1.gt.pind(jjj,7)) powind(j,jjj,7)=powind(j,jjj,7)+
     +       1.d0/fn2
            if(SCORE1.gt.percas(jjj,7)) powas(j,jjj,7)=powas(j,jjj,7)+
     +       1.d0/fn2
            if(SCORE2.gt.pind(jjj,8)) powind(j,jjj,8)=powind(j,jjj,8)+
     +       1.d0/fn2
            if(SCORE2.gt.percas(jjj,8)) powas(j,jjj,8)=powas(j,jjj,8)+
     +       1.d0/fn2
            if(SCOREBIC.gt.pind(jjj,9)) powind(j,jjj,9)=powind(j,jjj,9)+
     +       1.d0/fn2
            if(SCOREBIC.gt.percas(jjj,9)) powas(j,jjj,9)=powas(j,jjj,9)+
     +       1.d0/fn2
            if(XLRT1.gt.pind(jjj,10)) powind(j,jjj,10)=powind(j,jjj,10)+
     +       1.d0/fn2
            if(XLRT1.gt.percas(jjj,10)) powas(j,jjj,10)=powas(j,jjj,10)+
     +       1.d0/fn2
            if(XLRT2.gt.pind(jjj,11)) powind(j,jjj,11)=powind(j,jjj,11)+
     +       1.d0/fn2
            if(XLRT2.gt.percas(jjj,11)) powas(j,jjj,11)=powas(j,jjj,11)+
     +       1.d0/fn2
            if(XLRTBIC.gt.pind(jjj,12)) powind(j,jjj,12)=
     +       powind(j,jjj,12)+1.d0/fn2
            if(XLRTBIC.gt.percas(jjj,12)) powas(j,jjj,12)=
     +       powas(j,jjj,12)+1.d0/fn2
            if(VK.gt.pind(jjj,13)) powind(j,jjj,13)=powind(j,jjj,13)+
     +       1.d0/fn2
            if(pvalVKas.lt.testsize2(jjj)) powas(j,jjj,13)=
     +       powas(j,jjj,13)+1.d0/fn2
            if(KSFR(1).gt.pind(jjj,14)) powind(j,jjj,14)=
     +       powind(j,jjj,14)+1.d0/fn2
            if(KSFR(6).lt.testsize2(jjj)) powas(j,jjj,14)=
     +       powas(j,jjj,14)+1.d0/fn2
            if(VKFR.gt.pind(jjj,15)) powind(j,jjj,15)=powind(j,jjj,15)+
     +       1.d0/fn2
            if(pvalVKFRas.lt.testsize2(jjj)) powas(j,jjj,15)=
     +       powas(j,jjj,15)+1.d0/fn2
c Combined test X2-Q
            pvalmin=pvalX2
            if(pvalQ.le.pvalmin) pvalmin=pvalQ            
            if(pvalmin.lt.perccomb(jjj,1)) powcomb(j,jjj,1)=
     +       powcomb(j,jjj,1)+1.d0/fn2
            pvalX2Q=pvalmin        
c Combined test Q-KS
            pvalmin=pvalKS
            if(pvalQ.le.pvalmin) pvalmin=pvalQ                        
            if(pvalmin.lt.perccomb(jjj,2)) powcomb(j,jjj,2)=
     +       powcomb(j,jjj,2)+1.d0/fn2
            pvalQKS=pvalmin        
c Combined test X2-Q-KS
            pvalmin=pvalX2
            if(pvalQ.le.pvalmin) pvalmin=pvalQ            
            if(pvalKS.le.pvalmin) pvalmin=pvalKS            
            if(pvalmin.lt.perccomb(jjj,3)) powcomb(j,jjj,3)=
     +       powcomb(j,jjj,3)+1.d0/fn2
            pvalX2QKS=pvalmin        
c Combined test X2DIF-Q-KS
            pvalmin=pvalX2DIF
            if(pvalKS.le.pvalmin) pvalmin=pvalKS
            if(pvalQ.le.pvalmin) pvalmin=pvalQ           
            if(pvalmin.lt.perccomb(jjj,4)) powcomb(j,jjj,4)=
     +       powcomb(j,jjj,4)+1.d0/fn2
            pvalX2DIFKS=pvalmin        
c Combined test X2DIF2T-Q-KS   
            pvalmin=pvalX2DIF2T
            if(pvalKS.le.pvalmin) pvalmin=pvalKS
            if(pvalQ.le.pvalmin) pvalmin=pvalQ           
            if(pvalmin.lt.perccomb(jjj,5)) powcomb(j,jjj,5)=
     +       powcomb(j,jjj,5)+1.d0/fn2
            pvalX2DIF2TKS=pvalmin        
c Combined test X2DIF-Q-SCOREBIC
            pvalmin=pvalX2DIF
            if(pvalQ.le.pvalmin) pvalmin=pvalQ           
            if(pvalSCOREBIC.le.pvalmin) pvalmin=pvalSCOREBIC           
            if(pvalmin.lt.perccomb(jjj,6)) powcomb(j,jjj,6)=
     +       powcomb(j,jjj,6)+1.d0/fn2
            pvalX2DIFQSCORE=pvalmin        
c Combined test X2DIF-Q-XLRTBIC
            pvalmin=pvalX2DIF
            if(pvalQ.le.pvalmin) pvalmin=pvalQ           
            if(pvalXLRTBIC.le.pvalmin) pvalmin=pvalXLRTBIC           
            if(pvalmin.lt.perccomb(jjj,7)) powcomb(j,jjj,7)=
     +       powcomb(j,jjj,7)+1.d0/fn2
            pvalX2DIFQXLRT=pvalmin        
c Combined test X2DIF-Q-Kuiper
            pvalmin=pvalX2DIF
            if(pvalQ.le.pvalmin) pvalmin=pvalQ           
            if(pvalVK.le.pvalmin) pvalmin=pvalVK           
            if(pvalmin.lt.perccomb(jjj,8)) powcomb(j,jjj,8)=
     +       powcomb(j,jjj,8)+1.d0/fn2
            pvalX2DIFQVK=pvalmin 
c Combined test X2DIF-KSFR
            pvalmin=pvalX2DIF
            if(pvalKSFR.le.pvalmin) pvalmin=pvalKSFR
            if(pvalmin.lt.perccomb(jjj,9)) powcomb(j,jjj,9)=
     +       powcomb(j,jjj,9)+1.d0/fn2
            pvalX2DIFKSFR=pvalmin
c Combined test X2DIF-VKFR
            pvalmin=pvalX2DIF
            if(pvalVKFR.le.pvalmin) pvalmin=pvalVKFR
            if(pvalmin.lt.perccomb(jjj,10)) powcomb(j,jjj,10)=
     +       powcomb(j,jjj,10)+1.d0/fn2
            pvalX2DIFVKFR=pvalmin        
c
          end do
c NEW OPTION: always save test statistics + exact p-values if *idistr < 0*
c     if *idistr.lt.(-1)* (emprical analysis of multiple data sets)
c     if *idistr.eq.(-1)* (emprical analysis of single data set)
c     isavetest is not required
c
c Exact p-values for X2_2d + SCORE1 + SCORE2 + XLRT1 + XLRT2 are not (yet)
c computed in empirical analysis of multiple data sets ==> the p-values 
c are computed in the loop for combined tests - TO DO
c for the moment save test statistics only ==> compute exact p-values 
c for a *single data set* from output files of test statistics
c (X2_2d is not sorted: see use routine SVRGN instead of SVRGP)
c
          if(idistr.lt.0) then          
            write(51,50) X2,pvalX2
c            write(52,50) X2_2d,pvalX2_2d
            write(52,51) X2_2d
            write(53,50) Q,pvalQ
            write(54,54) (KS(jjjj),jjjj=1,6),pvalKS
            write(55,50) X2DIF,pvalX2DIF
            write(56,50) X2DIF2T,pvalX2DIF2T
c            write(57,50) SCORE1,pvalSCORE1
c            write(58,50) SCORE2,pvalSCORE2
            write(57,51) SCORE1
            write(58,51) SCORE2
            write(59,52) SCOREBIC,ncomps,pvalSCOREBIC
c            write(60,50) XLRT1,pvalXLRT1
c            write(61,50) XLRT2,pvalXLRT2
            write(60,51) XLRT1
            write(61,51) XLRT2
            write(62,52) XLRTBIC,ncompx,pvalXLRTBIC            
            write(63,55) VK,pvalVKas,pvalVK           
            write(64,54) (KSFR(jjjj),jjjj=1,6),pvalKSFR
            write(65,55) VKFR,pvalVKFRas,pvalVKFR  
          end if
50        format(f12.6,8x,f9.6)
55        format(f12.6,8x,f9.6,8x,f9.6)
54        format(6f12.6,8x,f9.6)            
51        format(f12.6)
52        format(f12.6,4x,i4,8x,f9.6)
c
C
C END OF SIMULATION LOOP OVER *nsim2* SIMULATIONS 
        end do 
cc      close(80)
cc      close(81)       
c Close files for saving test statistics and p-values
c if *idistr < -1* (emprical analysis of multiple data sets) 
      if(idistr.lt.(-1)) then
        CLOSE(51)
        CLOSE(52)
        CLOSE(53)
        CLOSE(54)
        CLOSE(55)  
        CLOSE(56)
        CLOSE(57)
        CLOSE(58)
        CLOSE(59)
        CLOSE(60)
        CLOSE(61)
        CLOSE(62)  
        CLOSE(63) 
        CLOSE(64)
        CLOSE(65)     
      end if
      CLOSE(11)
      CLOSE(12)
c
c Compute summary statistics + correlation beteween X2 and Q 
c skip this step when idistr < 0 (empirical analysis)
        if(idistr.ge.0) then
         do jjj=1,2
          corX2Q(jjj,jjj)=corX2Q(jjj,jjj)-(X2Qmean(jjj)*X2Qmean(jjj))
         end do
         corX2Q(1,2)=corX2Q(1,2)-(X2Qmean(1)*X2Qmean(2))
         corX2Q(1,2)=corX2Q(1,2)/((corX2Q(1,1)*corX2Q(2,2))**0.5)
         corX2Q(2,1)=corX2Q(1,2)
         write(73,630) (X2Qmean(jjj),jjj=1,2)
         write(73,*)
         write(73,631) (corX2Q(jjj,jjj),jjj=1,2)
         write(73,*)
         write(73,632) corX2Q(1,2)
         write(73,*)
         write(73,733) nsim2,n,j,idistr
733      format("Quadratic forms: X2 and Q",4x,"# power sim =",i7,4x,
     +   "n =",i6,4x,"parameter # =",i3,4x,"Distribution # =",i3)
         write(73,*)
         write(73,*)
c Compute summary statistics + correlation beteween X2 and X2DIF 
         do jjj=1,2
           corX2DIF(jjj,jjj)=corX2DIF(jjj,jjj)-
     +     (X2DIFmean(jjj)*X2DIFmean(jjj))
         end do
         corX2DIF(1,2)=corX2DIF(1,2)-(X2DIFmean(1)*X2DIFmean(2))
        corX2DIF(1,2)=corX2DIF(1,2)/((corX2DIF(1,1)*corX2DIF(2,2))**0.5)
         corX2DIF(2,1)=corX2DIF(1,2)
         write(75,630) (X2DIFmean(jjj),jjj=1,2)
         write(75,*)
         write(75,631) (corX2DIF(jjj,jjj),jjj=1,2)
         write(75,*)
         write(75,632) corX2DIF(1,2)
         write(75,*)
         write(75,6322) checkX2DIF
         write(75,*)
         write(75,734) nsim2,n,j,idistr
734      format("Quadratic forms: X2 and X2DIF = Q - X2",4x,"# power ",
     +   "sim =",i7,
     +   4x,"n =",i6,4x,"parameter # =",i3,4x,"Distribution # =",i3)
         write(75,*)
         write(75,*)
c Compute summary statistics + correlation beteween X2 and X2DIF2T 
         do jjj=1,2
           corX2DIF2T(jjj,jjj)=corX2DIF2T(jjj,jjj)-
     +     (X2DIF2Tmean(jjj)*X2DIF2Tmean(jjj))
         end do
         corX2DIF2T(1,2)=corX2DIF2T(1,2)-(X2DIF2Tmean(1)*X2DIF2Tmean(2))
         corX2DIF2T(1,2)=corX2DIF2T(1,2)/((corX2DIF2T(1,1)*
     +   corX2DIF2T(2,2))**0.5)
         corX2DIF2T(2,1)=corX2DIF2T(1,2)
         write(77,630) (X2DIF2Tmean(jjj),jjj=1,2)
         write(77,*)
         write(77,631) (corX2DIF2T(jjj,jjj),jjj=1,2)
         write(77,*)
         write(77,632) corX2DIF2T(1,2)
         write(77,*)
         write(77,6322) checkX2DIF2T
         write(77,*)
         write(77,737) nsim2,n,j,idistr
737      format("Quadratic forms: X2 and X2DIF2T = |Q - X2|",4x,
     +   "# power sim =",i7,
     +   4x,"n =",i6,4x,"parameter # =",i3,4x,"Distribution # =",i3)
         write(77,*)
         write(77,*)
        end if
c        
c Compute empirical standard deviation of values
c (exclude the Stable distribution)
        if(idistr.ne.7) then
          checkvar(j)=checkvar(j)-checkmean(j)*checkmean(j)
          checkvar(j)=dsqrt(checkvar(j))        
        end if
C      
C END OF LOOP OVER *npartrue* PARAMETER VALUES
      end do
      close(73)
      close(75)  
      close(77)
c
      CLOSE(410)
      CLOSE(420)
      CLOSE(430)
      CLOSE(440)
      CLOSE(450)  
      CLOSE(460)
      CLOSE(470)
      CLOSE(480)
      CLOSE(490)
      CLOSE(500)
      CLOSE(510)
      CLOSE(520)
      CLOSE(530)  
      CLOSE(540)
      CLOSE(550)
c CLOSE FILES FOR SAVING SIGNIFICANDS      
      close(101)
      if(idistr.eq.14) goto 1020
      if(idistr.eq.16) goto 1020
      if(idistr.gt.0) then
        close(102)
        close(103)
        close(104)
        if(idistr.eq.15) goto 1020      
        close(105)
      end if
      if(idistr.eq.2) goto 1020
      if(idistr.eq.7) goto 1020
      close(106)
1020  continue      
c      
      if(idistr.ge.0) then
        write(*,*) "END OF POWER SIMULATION"
      else 
        write(*,*) "END OF EMPIRICAL ANALYSIS"
      end if     
C
C
C WRITE POWER RESULTS IN CORRESPONDING FILE
C Empirical analysis of data read from file (idistr < 0)
C Write output file of tests statistics + p-values
C Include information about idistr = -1 or idistr < -1
      if(idistr.lt.0) then
        OPEN(unit=9,file='outDATA.txt')
        write(9,991) n,nsim1,ineg,ijit
991     format(1x,"** ANALYSIS OF EMPIRICAL DATA: see input file **",//
     +   8x,"n =",i9//
     +   1x,"# sim for null quantile estimation",i12,//
     +  1x,"Check on negative values in data (if ineg=1): ineg =",i3,//,
     +   1x,"Jittering in the data (if ijit=1): ijit =",i3,//
     +   1x,"(Not valid if nsim2>1)")
        if(idistr.eq.(-1)) then
          write(9,*)
          write(9,9911) nsim2
9911    format(1x,"SINGLE DATA SET: nsim2 = ",i9)
        else
          write(9,*)
          write(9,9912) nsim2
9912    format(1x,"MULTIPLE DATA SETS: nsim2 = ",i9,3x,"(see output ",
     +   "files for individual statistics and p-values)")
        end if
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300)
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
        OPEN(unit=99,file='outSTAT.txt')
        write(99,991) n,nsim1,ineg,ijit
        write(99,*)
        if(idistr.eq.(-1)) then
          write(99,*)
          write(99,9911) nsim2
        else
          write(99,*)
          write(99,9912) nsim2
        end if
        write(99,*)
        write(99,*) "** TEST STATISTICS AND P-VALUES **"
        write(99,*) "ONLY * LAST SAMPLE *  IF * nsim2 > 1 *"
        write(99,*) "SEE  * OUTPUT FILES *  FOR ALL nsim2 VALUES"
        write(99,*) "NO  * COMBINED TESTS *  IF *itest=1*"
        write(99,*) "NO  * EXACT P-VALUES *  IF *itest=1*"
        write(99,*) "itest = ",itest
        write(99,*)
        df=8.d0
        qqq=1.d0-dchidf(X2,df)
        write(99,9981)X2,qqq,pvalX2
9981    format(1x,"X2 =",f12.6,9x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)
        df=89.d0
        qqq=1.d0-dchidf(X2_2d,df)
!        write(99,9982)X2_2d,qqq,pvalX2_2d
!9982    format(1x,"X2_2d =",f12.6,6x,"Asymptotic p-value =",f9.6,6x,
!     +   "Exact p-value =",f9.6) 
       write(99,9982)X2_2d,qqq
9982    format(1x,"X2_2d =",f12.6,6x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value = TO DO")
        df=9.d0
        qqq=1.d0-dchidf(Q,df)
        write(99,9983)Q,qqq,pvalQ
9983    format(1x,"Q =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)
        write(99,9984)KS(1),KS(6),pvalKS
9984    format(1x,"KS =",f12.6,9x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)
        df=1.d0
        qqq=1.d0-dchidf(X2DIF,df)
        write(99,9985)X2DIF,qqq,pvalX2DIF
9985    format(1x,"X2DIF =",f12.6,6x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)
        df=1.d0
        qqq=1.d0-dchidf(X2DIF2T,df)
        write(99,9986)X2DIF2T,qqq,pvalX2DIF2T
9986    format(1x,"X2DIF2T =",f12.6,4x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)
        df=2.d0
        qqq=1.d0-dchidf(SCORE1,df)
!        write(99,9987)SCORE1,qqq,pvalSCORE1
!9987    format(1x,"SCORE1 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
!     +   "Exact p-value =",f9.6)
        write(99,9987)SCORE1,qqq
9987    format(1x,"SCORE1 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value = TO DO")
        df=4.d0
        qqq=1.d0-dchidf(SCORE2,df)
!        write(99,9988)SCORE2,qqq,pvalSCORE2
!9988    format(1x,"SCORE2 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
!     +   "Exact p-value =",f9.6)
        write(99,9988)SCORE2,qqq
9988    format(1x,"SCORE2 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value = TO DO")
        df=2.d0
        qqq=1.d0-dchidf(SCOREBIC,df)
        write(99,9979)SCOREBIC,qqq,pvalSCOREBIC
9979    format(1x,"SCOREBIC =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)
        df=2.d0
        qqq=1.d0-dchidf(XLRT1,df)
!        write(99,9980)XLRT1,qqq,pvalXLRT1
!9980    format(1x,"XLRT1 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
!     +   "Exact p-value =",f9.6)
        write(99,9980)XLRT1,qqq
9980    format(1x,"XLRT1 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value = TO DO")
        df=4.d0
        qqq=1.d0-dchidf(XLRT2,df)
!        write(99,9918)XLRT2,qqq,pvalXLRT2
!9918    format(1x,"XLRT2 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
!     +   "Exact p-value =",f9.6)    
        write(99,9918)XLRT2,qqq
9918    format(1x,"XLRT2 =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value = TO DO")      
        df=2.d0
        qqq=1.d0-dchidf(XLRTBIC,df)
        write(99,9928)XLRTBIC,qqq,pvalXLRTBIC
9928    format(1x,"XLRTBIC =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)
!        qqq=0.d0
        write(99,9938)VK,pvalVKas,pvalVK
9938    format(1x,"VK =",f12.6,10x,"Asymptotic p-value =",f9.6,6x,
     +   "Exact p-value =",f9.6)   
        write(99,99841)KSFR(1),KSFR(6),pvalKSFR
99841   format(1x,"KS on the significand *fractional part* =",f12.6,9x,
     +   "Asymptotic p-value =",f9.6,6x,"Exact p-value =",f9.6)   
        write(99,99381)VKFR,pvalVKFRas,pvalVKFR
99381   format(1x,"VK on the significand *fractional part* =",f12.6,10x,
     +   "Asymptotic p-value =",f9.6,6x,"Exact p-value =",f9.6)     
        write(99,9948)
9948   format(1x,"NB: exact p-values for SCORE1-SCORE2-XLRT1-XLRT2 not", 
     +  " yet computed in empirical analysis of multiple data sets",/,
     +  " (they are computed in the loop for combined tests)",/,
     +  "compute exact p-values from *output files* of test statistics")  
        write(99,*)     
        write(99,*)     
        write(99,9991)pvalX2Q
9991    format(1x,"Combined test X2-Q:          Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,1),jjj=1,nqprop)
9990    format(1x,"Quantiles of Min p-value of combined test:",6f9.6)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=138,file='benfX2Qsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(138,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2Q) pvalcom=pvalcom+1.0/float(nsim1)
        end do
        CLOSE(138)
        write(99,9989) pvalcom
9989    format(1x,"Exact p-value of Min p-value test statistics:",f9.6)      
        write(99,*)     
        write(99,9992)pvalQKS
9992    format(1x,"Combined test Q-KS:          Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,2),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=139,file='benfQKSsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(139,30) benfarray(jj)
          if(benfarray(jj).lt.pvalQKS) pvalcom=pvalcom+1.0/float(nsim1)
        end do
        CLOSE(139)
        write(99,9989) pvalcom
        write(99,*)     
        write(99,9993)pvalX2QKS
9993    format(1x,"Combined test X2-Q-KS:       Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,3),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=140,file='benfX2QKSsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(140,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2QKS) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(140)
        write(99,9989) pvalcom        
        write(99,*)     
        write(99,9994)pvalX2DIFKS
9994    format(1x,"Combined test X2DIF-Q-KS:    Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,4),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=141,file='benfX2DIFKSsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(141,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2DIFKS) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(141)
        write(99,9989) pvalcom
        write(99,*)     
        write(99,9995)pvalX2DIF2TKS
9995    format(1x,"Combined test X2DIF2T-Q-KS:  Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,5),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=142,file='benfX2DIF2TKSsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(142,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2DIF2TKS) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(142)
        write(99,9989) pvalcom
        write(99,*)    
        write(99,9996)pvalX2DIFQSCORE
9996    format(1x,"Combined test X2DIF-Q-SCOREBIC:  Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,6),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=143,file='benfX2DIFQSCOREsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(143,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2DIFQSCORE) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(143)
        write(99,9989) pvalcom
        write(99,*)    
        write(99,9997)pvalX2DIFQXLRT
9997    format(1x,"Combined test X2DIF-Q-XLRTBIC:  Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,7),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=144,file='benfX2DIFQXLRTsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(144,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2DIFQXLRT) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(144)
        write(99,9989) pvalcom
        write(99,*)    
        write(99,9998)pvalX2DIFQVK
9998    format(1x,"Combined test X2DIF-Q-KUIPER:  Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,8),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=145,file='benfX2DIFQVKsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(145,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2DIFQVK) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(145)
        write(99,9989) pvalcom
        write(99,*)     
        write(99,99941)pvalX2DIFKSFR
99941    format(1x,"Combined test X2DIF-KSFR:    Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,9),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=141,file='benfX2DIFKSFRsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(141,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2DIFKSFR) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(141)
        write(99,9989) pvalcom     
        write(99,*)     
        write(99,99942)pvalX2DIFVKFR
99942   format(1x,"Combined test X2DIF-VKFR:    Min p-value =",f9.6)
        write(99,9990) (perccomb(jjj,10),jjj=1,nqprop)
c Compute exact p-value of Min p-value statistic
        OPEN(unit=141,file='benfX2DIFVKFRsort.txt')
        pvalcom=0.0
        do jj=1,nsim1
          read(141,30) benfarray(jj)
          if(benfarray(jj).lt.pvalX2DIFVKFR) pvalcom=pvalcom+
     +     1.0/float(nsim1)
        end do
        CLOSE(141)
        write(99,9989) pvalcom
        write(99,*)    
        write(99,*)
        write(99,899) (testsize(j),j=1,nqprop)
899     format(30x,"Test sizes: ",6f10.6)
        CLOSE(99)
      end if
c      
      if(idistr.eq.0) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFBENF.txt')      
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outBENF.txt')
        end if
        write(9,900) nsim2,n,npartrue,idistr,nsim1
900     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* BENFORD RANDOM VARIABLE: check size *",//,
     +   1x,"# sim for null quantile estimation",i12)      
        write(9,*)
        write(9,9000) 
9000    format(1x,"Check of significands")
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
          write(9,9002) (checksigt(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9100) 
9100    format(1x,"Check of first-digit probabilities")
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
          write(9,9002) (checkdigt(jj,j),jj=1,9)
        end do
9010    format(1x,"Parameter value j = ",i3)
9001    format(1x,"Empirical mean:               ",9f12.6)
9002    format(1x,"Nominal (assuming rcont=1.0): ",9f12.6)
      end if
c
      if(idistr.eq.1) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFNORM.txt')      
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outNORM.txt')              
        end if
961     format(1x,"** WARNING **  **MANIPULATED BENFORD DISTRIBUTION **    
     +   ==> 1ST DIGIT IS BENFORD")
        write(9,901) nsim2,n,npartrue,idistr,(par(j),j=1,npar2),
     +   par(npar/2+1),par(npar/2+2),nsim1
901     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* NORMAL RANDOM VARIABLE *",//,
     +   1x,"Values of the scale parameter (standard deviation):",
     +   3f10.3,//,
     +   1x,"Values of the location parameter:",2f10.3,//
     +   1x,"Parameter setting 1: scale1*location1",/,
     +   1x,"Parameter setting 2: scale2*location1",/,
     +   1x,"Parameter setting 3: scale3*location1",/,
     +   1x,"Parameter setting 4: scale1*location2",/,
     +   1x,"Parameter setting 5: scale2*location2",/,
     +   1x,"Parameter setting 6: scale3*location2",//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
8000    format(1x,"Check of mean value ")
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
8100    format(1x,"Check of standard deviation ")
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,8012) (checkvart(j), j=1,npartrue)
8001    format(1x,"Empirical means:                ",6f15.6)
8002    format(1x,"Nominal means (assume rcont=1): ",6f15.6)
8011    format(1x,"Empirical std dev:                ",6f15.6)
8012    format(1x,"Nominal std dev (assume rcont=1): ",6f15.6)
        write(9,*)
        write(9,9200) 
9200    format(1x,"Check of significands (empirical only)")        
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300)
9300    format(1x,"Check of first-digit probabilities",
     +   " (empirical only)") 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.2) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFLOGNORM.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outLOGNORM.txt')        
        end if
        write(9,902) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   par(npartrue+1),nsim1
902     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* LOGNORMAL RANDOM VARIABLE *",//,
     +   1x,"Values of the shape parameter:",
     +   5f10.3,//,
     +   1x,"Fixed value of the scale parameter:",f10.3,//
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do     
      end if
c      
      if(idistr.eq.3) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFWEIB.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outWEIB.txt')        
        end if        
        write(9,903) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   nsim1
903     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* (Standard) WEIBULL RANDOM VARIABLE *",//,
     +   1x,"Values of the shape parameter:",
     +   6f10.3,//,
     +   1x,"The scale parameter is fixed and set = 1",//
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.4) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFGAM.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outGAM.txt')        
        end if        
        write(9,904) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   nsim1
904     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* (Standard) GAMMA RANDOM VARIABLE *",//,
     +   1x,"Values of the shape parameter:",
     +   6f10.3,//,
     +   1x,"The scale parameter is fixed and set = 1",//
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.5) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFBETA.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outBETA.txt')          
        end if        
        write(9,905) nsim2,n,npartrue,idistr,(par(j),j=1,npar2),
     +   par(npar/2+1),par(npar/2+2),nsim1
905     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* BETA RANDOM VARIABLE *",//,
     +   1x,"Values of the shape parameter alpha:",
     +   3f10.3,//,
     +   1x,"Values of the shape parameter beta: ",2f10.3,//
     +   1x,"Parameter setting 1: alpha1*beta1",/,
     +   1x,"Parameter setting 2: alpha2*beta1",/,
     +   1x,"Parameter setting 3: alpha3*beta1",/,
     +   1x,"Parameter setting 4: alpha1*beta2",/,
     +   1x,"Parameter setting 5: alpha2*beta2",/,
     +   1x,"Parameter setting 6: alpha3*beta2",//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.6) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFUNIF.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outUNIF.txt')        
        end if      
        write(9,906) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   nsim1
906     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* UNIFORM RANDOM VARIABLE on (0,a) *",//,
     +   1x,"Values of the range parameter:",
     +   6f10.3,//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.7) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFSTAB.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outSTAB.txt')        
        end if        
        write(9,907) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   par(npartrue+1),(betap(j),j=1,npartrue),nsim1
907     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* STABLE RANDOM VARIABLE *",//,
     +   1x,"Values of the shape parameter:",
     +   5f10.3,//,
     +   1x,"Fixed value of the skewness parameter ",
     +    " (used when istab=1):",f10.3,//
     +   1x,"Betaprime = skewness after reparametrization:",5f10.3,//
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,970) (tolstab(j),j=1,npartrue)
970     format(1x,"Check of tail probabilities for parameter-dependent",
     +   " tail values: (suitable when istab=1)",/
     +   22x,5f16.2)
        write(9,971) (checkstab(j),j=1,npartrue)
971     format(1x,"Empirical probabilities:  ",5f15.12)
        write(9,972) (checkstabt(j),j=1,npartrue)
972     format(1x,"Nominal probabilities:    ",5f15.12)
        write(9,*)
        if(istab.eq.1) write(9,973) istab
973     format("istab =",i4,4x,"routine DRNSTA")
        if(istab.eq.2) write(9,974) istab
974     format("istab =",i4,4x,"Kanter algorithm for positive Stable")
        if(istab.eq.3) write(9,975) istab
975     format("istab =",i4,4x,"Kanter algorithm for positive Stable,"
     +   " as presented by Chambers et al. (JASA, 1976, Eq. (2.2))")
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.8) then
        OPEN(unit=9,file='outGBL.txt')
        write(9,908) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   nsim1
908     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* GENERALIZED BENFORD RANDOM VARIABLE *",//,
     +   1x,"Values of the parameter:",
     +   6f10.3,//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,9000) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
          write(9,9002) (checksigt(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9100) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
          write(9,9002) (checkdigt(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.9) then
        OPEN(unit=9,file='outBENFGBL.txt')
        write(9,909) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   nsim1
909     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* BENFGBL RANDOM VARIABLE: First-Digit Benford + ",
     +   "Generalized Benford *",//,
     +   1x,"Values of the parameter:",
     +   6f10.3,//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.10) then
        OPEN(unit=9,file='outMIXGBL.txt')
        write(9,910) nsim2,n,npartrue,idistr,(par(j),j=1,npar2),
     +   par(npar/2+1),par(npar/2+2),nsim1
910     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* 2-COMPONENT GENERALIZED BENFORD MIXTURE DISTRIBUTION *",
     +   //,
     +   1x,"Values of the GBL parameter",/
     +   1x,"> 0 in the 1st component and < 0 in the second component:",
     +   3f10.3,//,
     +   1x,"Values of the mixing proportion for *component 1*: ",
     +   2f10.3,//
     +   1x,"Parameter setting 1: par1*mix1",/,
     +   1x,"Parameter setting 2: par2*mix1",/,
     +   1x,"Parameter setting 3: par3*mix1",/,
     +   1x,"Parameter setting 4: par1*mix2",/,
     +   1x,"Parameter setting 5: par2*mix2",/,
     +   1x,"Parameter setting 6: par3*mix2",//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.11) then
        OPEN(unit=9,file='outDIRAC.txt')
        write(9,911) nsim2,n,npartrue,idistr,(par(j),j=1,npartrue),
     +   nsim1
911     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* FIRST-DIGIT DIRAC DISTRIBUTION: *",//,
     +   1x,"Fixed values of the replacing first digit:",
     +   6f10.3,//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.12) then
        OPEN(unit=9,file='outTRUNC.txt')
        write(9,912) nsim2,n,npartrue,idistr,par(1),nsim1
912     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     + 25x,"* TRUNCATED BENFORD ==> LAST-DIGIT DIRAC DISTRIBUTION: *",/,
     + 25x,"  + JITTERING IN THE 10th DECIMAL PLACE TO AVOID ",
     +  "TIES IN KS TEST (more jittering when n>=200)",//
     +   1x,"Fixed value of the replacing last digit(s):",
     +   f10.3,//, 
     +   1x,"Parameter setting 1: 6 decimal digits are retained, the",
     +   " others are set to the fixed value",/, 
     +   1x,"Parameter setting 2: 5 decimal digits are retained, the",
     +   " others are set to the fixed value",/, 
     +   1x,"Parameter setting 3: 4 decimal digits are retained, the",
     +   " others are set to the fixed value",/, 
     +   1x,"Parameter setting 4: 3 decimal digits are retained, the",
     +   " others are set to the fixed value",/, 
     +   1x,"Parameter setting 5: 2 decimal digits are retained, the",
     +   " others are set to the fixed value",/, 
     +   1x,"Parameter setting 6: 1 decimal digit is retained, the",
     +   " others are set to the fixed value",//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c
      if(idistr.eq.120) then
        OPEN(unit=9,file='outTRUNCNEW.txt')
        write(9,9120) nsim2,n,npartrue,idistr,(j,j=1,npar),
     +   (par(j),j=1,npar),(ntrunc(j),j=1,npar),(cumtrunc(j),j=1,npar),
     +   nsim1
9120    format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +25x,"*NEW TRUNCATED BENFORD* ==> DISTRIBUTION WITH LAST DIGIT(S)",
     + " 0 IN FIXED PROPORTIONS OF CASES (TO BE READ AS INPUT)",/,
     + 25x,"  + JITTERING IN THE 10th DECIMAL PLACE TO AVOID ",
     +  "TIES IN KS TEST (more jittering when n>=200)",//,
     +   1x,"# of significant digits j = ",9x,6i18,//,
     +   1x,"Proportion with j significant digits",3x,6f18.3,//, 
     +   1x,"# of cases with j significant digits",1x,6i18,//, 
     +   1x,"Cumulative # of cases",16x,6i18,//, 
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.13) then
        OPEN(unit=9,file='outCARD.txt')
        xx1=(par(npar/2+1))/(2.d0*pi)
        xx2=(par(npar/2+2))/(2.d0*pi)
        write(9,913) nsim2,n,npartrue,idistr,(par(j),j=1,npar2),
     +   xx1,xx2,par(npar2+1),par(npar2+2),nsim1
913     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* CARDIOID RANDOM VARIABLE *",//,
     +   1x,"Values of the concentration parameter 0<a<1:",
     +   3f10.3,//,
     +   1x,"Values of the location parameter 0<b<1",2f10.3,/
     +   1x,"Values of the location parameter b (after multiplication ",
     +    "by 2*pi)",2f10.3,//
     +   1x,"Parameter setting 1: concentration1*location1",/,
     +   1x,"Parameter setting 2: concentration2*location1",/,
     +   1x,"Parameter setting 3: concentration3*location1",/,
     +   1x,"Parameter setting 4: concentration1*location2",/,
     +   1x,"Parameter setting 5: concentration2*location2",/,
     +   1x,"Parameter setting 6: concentration3*location2",//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300)
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if     
c      
      if(idistr.eq.14) then
        OPEN(unit=9,file='outNNTS.txt')
        write(9,914) nsim2,n,npartrue,idistr,numNNTS,nsim1
914     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* NON NEGATIVE TRIGONOMETRIC SUMS (NNTS) VARIABLE *",//,
     +   1x,"Number of trigonometric sum (excluding intercept): ",
     +   "numNNTS =",i3,//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,9141) 
9141    format(10x,"Intercept of NNTS")
        write(9,*) parNNTS(1)
        write(9,9142) 
9142    format(10x,"Complex parameters of NNTS")
        do j=2,numNNTS+1
          write(9,*) parNNTS(j)
        end do
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
c        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
c        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,*)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300)
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if      
c
c NEW VERSION - OCTOBER 2022      
      if(idistr.eq.15) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENFMIXNORM.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='outMIXNORM.txt')        
        end if        
        write(9,915) nsim2,n,npartrue,idistr,(par(j),j=1,npar2),
     +   (par(j),j=npar2+1,npar2+2),par(npar2+3),par(npar2+4),nsim1
915     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"*2-COMPONENT NORMAL MIXTURE DISTRIBUTION (New version)*",
     +   //,
     +   1x,"Location component 1:",f10.3,6x,"Location component 2:",
     +   f10.3,/
     +   1x,"Scale component 1:",f10.3,6x,"Scale component 2:",
     +   f10.3,/
     +   1x,"Values of the mixing proportion for *component 1*: ",
     +   2f10.3,//
     +   1x,"Parameter setting 1: heteroschedastic mixture * mix1",/,
     +   1x,"Parameter setting 2: heteroschedastic mixture * mix2",/,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
c        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
c        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c      
      if(idistr.eq.16) then
c consider the case idig1=1 ==> 1st digit Benford
        if(idig1.eq.1) then
          OPEN(unit=9,file='outBENF3MIXNORM.txt')
          write(9,961)
          write(9,*)
          write(9,*)
        else
          OPEN(unit=9,file='out3MIXNORM.txt')        
        end if        
        write(9,916) nsim2,n,npartrue,idistr,(par(j),j=1,npar2),
     +   (par(j),j=npar2+1,npar),partrue,nsim1
916     format(1x,"POWER OF TESTS - EXACT + ASYMPTOTIC",//,1x,"# sim =",
     +   i9,
     +   8x,"n =",i6,8x,"# param =",i3,8x,"Distribution # =",i3,//,
     +   25x,"* 3-COMPONENT NORMAL MIXTURE DISTRIBUTION WITH FIXED ",
     +    "WEIGHTS *",
     +   //,
     +   1x,"Mixture location parameters       ",3f16.6/,
     +   1x,"Mixture scale (std dev) parameters",3f16.6/,
     +   1x,"Fixed mixing proportion",f12.6,//,
     +   1x,"# sim for null quantile estimation",i12)
        write(9,*)
        write(9,8000) 
        write(9,8001) (checkmean(j), j=1,npartrue)
c        write(9,8002) (checkmeant(j), j=1,npartrue)
        write(9,*)
        write(9,8100) 
        write(9,8011) (checkvar(j), j=1,npartrue)
c        write(9,8012) (checkvart(j), j=1,npartrue)
        write(9,9200) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checksig(jj,j),jj=1,9)
        end do
        write(9,*)
        write(9,9300) 
        do j=1,npartrue
          write(9,9010) j
          write(9,9001) (checkdig(jj,j),jj=1,9)
        end do
      end if
c    
c
      write(9,*)
      write(9,80) rcont,ncont,n
80    format(1x,"Contamination rate =",f10.6,8x,/
     + 1x,"# of contaminated observations in each sample",i8,
     + "   out of",i8)
      write(9,*)
      if(rcont.eq.1.d0) then
        write(9,81) 
      else
        write(9,82)
      end if
81    format(25x,"** FULL POWER ANALYIS **")
82    format(25x,"** CONTAMINATION MODEL FOR OUTLIER ANALYIS FROM ",
     + "BENFORD DISTRIBUTION **")
      write(9,*)
      write(9,*)
      if(iapprox.eq.1) then      
        write(9,83) (icomb(j),j=1,ntestcomb)
      else
        write(9,84)
      end if
83    format(1x,"Use faster but possibly *approximated* algorithm for",
     + " computing p-values of *combined* tests",/,
     + 1x,"Check approximations:",/
     + 1x,"icombX2Q     =",i4,/
     + 1x,"icombQKS     =",i4,/ 
     + 1x,"icombX2QKS   =",i4,/ 
     + 1x,"icombX2DIFKS =",i4,/
     + 1x,"icombX2DIF2TKS =",i4,/
     + 1x,"icombX2DIFQSCORE =",i4,/
     + 1x,"icombX2DIFQXLRT =",i4,/
     + 1x,"icombX2DIFQVK =",i4,/
     + 1x,"icombX2DIFKSFR =",i4,/
     + 1x,"icombX2DIFVKFR =",i4)
84    format(10x,"Use *full* algorithm for computing p-values of",
     + " *combined* tests")
      write(9,*)
      write(9,*)
      write(9,*)
      write(9,*) "                        Test sizes:"
      write(9,90) (testsize(j),j=1,nqprop)
90    format(23x,6f10.6)
      write(9,*)   
      write(9,*) 
      if(itest.eq.1) then
        write(9,88)itest
88      format(10x,"EXACT TESTS  -  * Combined tests are NOT performed *
     + :  itest =",i3,4x,"(see power of individual tests only)")
      else
        write(9,888)itest
888      format(20x,"EXACT TESTS  -  * Combined tests are performed *
     + :  itest =",i3)
      end if
      write(9,*)
      if(iread.eq.1) then 
        write(9,881) iread 
      else
        write(9,882) iread 
      end if
881   format(1x,"iread =",i4,6x,"Exact percentiles are read from file")
882   format(1x,"iread =",i4,6x,"Exact percentiles are simulated")      
      do jj=1,ntestind
        write(9,890) jj,(percind(j,jj),j=1,nqprop)
      end do
      do jj=1,ntestcomb
        write(9,891) jj,(perccomb(j,jj),j=1,nqprop)
      end do
890   format(1x,"Percentiles of individual test",i4,6x,6f12.6)    
891   format(1x,"Percentiles of combined test  ",i4,6x,6f12.6)    
      write(9,*)
      if(iread.eq.1) then
        write(9,8877)
8877    format(1x,"See saved files for results on BIC selection of # of" 
     +   " components under Benford's law")
      else
        write(9,8887) ncompmax
8887    format(1x,"Max # of components for BIC selection =",i6)        
        write(9,8888) avencompscoreb,varncompscoreb-avencompscoreb**2,
     +   maxncompscoreb
8888    format(1x,"Benford' law",/
     +   1x,"BIC selection of # of components for *SCORE TEST*",/
     +   1x,"AVE(ncomp) =",f12.6,8x,"VAR(ncomp) =",f12.6,/
     +   1x,"Prob(ncomp=ncompmax) =",f10.6)
        write(9,8889) avencompxlrtb,varncompxlrtb-avencompxlrtb**2,
     +   maxncompxlrtb   
8889    format(1x,"BIC selection of # of components for *LR TEST*",/
     +   1x,"AVE(ncomp) =",f12.6,8x,"VAR(ncomp) =",f12.6,/
     +   1x,"Prob(ncomp=ncompmax) =",f10.6)
      end if      
      write(9,*)      
      do j=1,npartrue
        write(9,91) j
91      format(1x,"Parameter setting",i3)
        do jj=1,ntestind
          write(9,92) jj,(powind(j,jjj,jj),jjj=1,nqprop)
92        format(1x,"Power for  *individual*  test ",i3,4x,6f10.6)
        end do
        do jj=1,ntestcomb
          write(9,922) jj,(powcomb(j,jjj,jj),jjj=1,nqprop)
922       format(1x,"Power for  *combined*    test ",i3,4x,6f10.6)
        end do
        write(9,*)
      end do
      write(9,*)   
      write(9,*) 
      write(9,89)
89    format(25x,"ASYMPTOTIC TESTS (individual tests only)",/
     + "(For KS the asymptotic p-value from DKSONE is used instead of ",
     +  "asymptotic percentiles)",/,
     + "(For Kuiper test the asymptotic p-value is computed instead of",
     +  "asymptotic percentiles)")
      write(9,*)
c Write percentiles
c For the Kolmogorov-Smirnov test do not compute percentiles but use 
c the asymptotic p-value from DKSONE
      do jj=1,ntestind
        write(9,890) jj,(percas(j,jj),j=1,nqprop)
      end do
      write(9,*)  
      do j=1,npartrue
        write(9,91) j
        do jj=1,ntestind
          write(9,92) jj,(powas(j,jjj,jj),jjj=1,nqprop)
        end do
        write(9,*)
      end do
      write(9,*)
      write(9,*)
      write(9,8887) ncompmax
      write(9,8688) 
      write(9,8788) (avencompscores(j),j=1,npar)
      write(9,8988) (varncompscores(j)-(avencompscores(j))**2,j=1,npar)
      write(9,8798) (maxncompscores(j),j=1,npar)
8688  format(1x,"Selected distribution",/
     + 1x,"BIC selection of # of components for *SCORE TEST*")
8788  format(1x,"AVE(ncomp) =",6f12.6)
8988  format(1x,"VAR(ncomp) =",6f12.6)
8798  format(1x,"Prob(ncomp=ncompmax) =",6f10.6)
      write(9,8789) 
8789  format(1x,"BIC selection of # of components for *LR TEST*")
      write(9,8788) (avencompxlrts(j),j=1,npar)
      write(9,8988) (varncompxlrts(j)-(avencompxlrts(j))**2,j=1,npar)
      write(9,8798) (maxncompxlrts(j),j=1,npar)
      write(9,*)
      write(9,*)
      write(9,*)
      write(9,*)
      write(9,*)" List of available *individual* tests:"
      write(9,*)" 1:   X2      = 1st digit Chi-square test"
      write(9,*)" 2:   X2_2d   = 2-digit Chi-square test"
      write(9,*)" 3:   Q       = Hotelling test based on sum invariance"
      write(9,*)" 4:   KS      = Kolmogorov-Smirnov two-sided test"
      write(9,*)" 5:   X2DIF   = Q - X2 test (one tail)"
      write(9,*)" 6:   X2DIF2T = |Q - X2| test (two tails: 09/06/21)"
      write(9,*)" 7:   SCORE1  = Score test with 1 component"
      write(9,*)" 8:   SCORE2  = Score test with 2 components"
      write(9,*)" 9:   SCOREBIC= Score test with BIC selection of ",
     +           "components"
      write(9,*)" 10:  XLRT1  = Likelihood ratio test with 1 component"
      write(9,*)" 11:  XLRT2  = Likelihood ratio test with 2 components"
      write(9,*)" 12:  XLRTBIC= Likelihood ratio test with BIC ",
     +           "selection of components"
      write(9,*)" 13:  VK     = Kuiper test"   
      write(9,*)" 14:  ZZZZZ - TO DO"
      write(9,*)" 15:  ZZZZZ - TO DO"      
      write(9,*)
      write(9,*)" List of available *combined* tests:"
      write(9,*)" 1:   X2Q     = Combined X2-Q test"
      write(9,*)" 2:   QKS     = Combined Q-KS test"
      write(9,*)" 3:   X2QKS   = Combined X2-Q-KS test"
      write(9,*)" 4:   X2DIFKS = Combined X2DIF-Q-KS test (new version)"
      write(9,*)" 5:   X2DIF2TKS= Combined X2DIF2T-Q-KS test"
      write(9,*)" 6:   X2DIFQSCORE= Combined X2DIF-Q-SCOREBIC test"
      write(9,*)" 7:   X2DIFQXLRT = Combined X2DIF-Q-XLRTBIC test"
      write(9,*)" 8:   X2DIFQVK   = Combined X2DIF-Q-KUIPER test"
      write(9,*)" 9:   X2DIFKSFR  = ZZZZZ - TO DO"
      write(9,*)" 10:  X2DIFVKFR  = ZZZZZ - TO DO"
      write(9,*)
      write(9,*)" Asymptotic approximations (for individual tests):"
      write(9,*)" 1:  X2      = Chi-square with df = 8"
      write(9,*)" 2:  X2_2d   = Chi-square with df = 89"
      write(9,*)" 3:  Q       = Chi-square with df = 9"
      write(9,*)" 4:  KS      = p-value from routine DKSONE"
      write(9,*)" 5:  X2DIF   = Chi-square with df = 1"
      write(9,*)" 6:  X2DIF2T = Chi-square with df = 1"
      write(9,*)" 7:  SCORE1  = Score test with 1 component  df = 2"
      write(9,*)" 8:  SCORE2  = Score test with 2 components df = 4"
      write(9,*)" 9:  SCOREBIC= Score test with BIC selection of ",
     +           "components  df = 2"
      write(9,*)" 10: XLRT1  = Likelihood ratio test with 1 component",
     +           " df = 2"
      write(9,*)" 11: XLRT2  = Likelihood ratio test with 2 components",
     *           " df = 4" 
      write(9,*)" 12: XLRTBIC= Likelihood ratio test with BIC ",
     +           "selection of components df = 2"
      write(9,*)" 13: VK       = Kuiper test (see references)"
      write(9,*)" 14: ZZZZZ - TO DO" 
      write(9,*)" 15: ZZZZZ - TO DO"
      write(9,*)
      write(9,*)
      write(9,*) "iseed0 =",iseed0
C
C
C COPY AND PASTE OUTPUT ==> (parameter setting) * (test)
C table of power values for each significance level
C CHECK THE VALUES OF *ntestind* and *ntestcomb* in format 98-980 
C                                               and format 99-990
C
97    format("Test size",f10.6)
980   format(1x,15i10)
98    format(i1,15f10.6)
990   format(1x,10i10)
99    format(i1,10f10.6)      
      write(9,*)
      write(9,*)
      write(9,*)
      write(9,*)
      write(9,*) "COPY AND PASTE OUTPUT ==> (parameter setting) * (test) 
     + table of power values for each significance level"
      write(9,*)
      write(9,*)
      do jjj=1,nqprop
        write(9,97) testsize(jjj)
        write(9,*)
        write(9,*) "Individual *exact* tests"
        write(9,980) (jj,jj=1,ntestind)
        do j=1,npartrue
          write(9,98) j,(powind(j,jjj,jj),jj=1,ntestind)
        end do
        write(9,*)
        write(9,*) "Individual *asymptotic* tests"
        write(9,980) (jj,jj=1,ntestind)
        do j=1,npartrue
          write(9,98) j,(powas(j,jjj,jj),jj=1,ntestind)
        end do
        write(9,*)
        write(9,*) "Combined *exact* tests"
        write(9,990) (jj,jj=1,ntestcomb)
        do j=1,npartrue
          write(9,99) j,(powcomb(j,jjj,jj),jj=1,ntestcomb)
        end do
        write(9,*)
        write(9,*)               
      end do      
C
C      
      CLOSE(9)
C
C
C
C
      call RNGET(iseed)
      write(*,*) iseed
c Saves input seed and current seed for random generation 
c of Benford values 
      open(unit=2,file='benf_iseed0.txt')
      write(2,*) iseed0
      close(2)
      open(unit=2,file='benf_iseed.txt')
      write(2,*) iseed
      close(2)
C
C
C
C
9999  CONTINUE
      write(*,*) "END OF PROGRAM"
      END
C
C 
C *************************************************************
C
      SUBROUTINE X2ONE(jsimul,n,maxn,d1,benf1n,X2)
C COMPUTES PEARSON X2 ON FIRST DIGIT
      IMPLICIT NONE
      integer*4 i,j,jj,jj1,jsimul,n,maxn,d1(maxn),d1freq(9)
      real*8 benf1n(9),X2,xx
c computes first-digit frequencies   
      do j=1,9
        d1freq(j)=0
      end do
      X2=0.d0   
      do i=1,n
        jj1=d1(i)
        if(jj1.lt.1.or.jj1.gt.9) then
          write(*,*) "X2ONE: error in 1st digit of sample unit",i,
     +    " d1 =",jj1
          pause
        end if
        d1freq(jj1)=d1freq(jj1)+1
      end do    
      do j=1,9
        xx = (dfloat(d1freq(j))-benf1n(j))**2
        X2 = X2 + xx/benf1n(j)
      end do
      RETURN
      END
C
C *************************************************************
C
      SUBROUTINE X2TWO(jsimul,n,maxn,d12,benf12n,X2_2d)
C COMPUTES PEARSON X2 ON FIRST-TWO DIGITS
      IMPLICIT NONE
      integer*4 i,j,jj,jj1,jj2,jsimul,n,maxn,d12(maxn,2),d12freq(9,10)
      real*8 benf12n(9,10),X2_2d,xx
c computes first-two digits frequencies   
      do j=1,9
        do jj=1,10
          d12freq(j,jj)=0
        end do
      end do
      X2_2d=0.d0   
      do i=1,n
        jj1=d12(i,1)
        jj2=d12(i,2)
        if(jj1.lt.1.or.jj1.gt.9) then
          write(*,*) "X2TWO: error in 1st digit of sample unit",i,
     +    " d1 =",jj1
          pause
        end if
        if(jj2.lt.0.or.jj1.gt.9) then
          write(*,*) "X2TWO: error in 2nd digit of sample unit",i,
     +    " d2 =",jj2
          pause
        end if
        d12freq(jj1,jj2+1)=d12freq(jj1,jj2+1)+1
      end do    
      do j=1,9
        do jj=1,10
          xx = (dfloat(d12freq(j,jj))-benf12n(j,jj))**2
          X2_2d = X2_2d + xx/benf12n(j,jj)
        end do
      end do
      RETURN
      END
C
C ******************************************************************
C
      SUBROUTINE HOTSUM(jsimul,n,maxn,d1,s,benfE,benfVARinv,Q)
C COMPUTES HOTELLING SUM-INVARIANCE TEST Q
      IMPLICIT NONE
c      integer*8 i,j,jj,jj1,jsimul,n,maxn,d1(maxn)
      integer*4 i,j,jj,jj1,jsimul,n,maxn,d1(maxn)
      real*8 s(maxn),benfE,benfVARinv(9,9),ave(9),zave(9),Q,a,b,v,fn
      fn=dfloat(n)
      do j=1,9
        zave(j)=0.d0
        ave(j)=benfE
      end do
c Compute average significand for sum-invariance test
      do i=1,n
        jj1=d1(i)
        if(jj1.lt.1.or.jj1.gt.9) then
          write(*,*) "HOTSUM: error in 1st digit of sample unit",i,
     +    " d1 =",jj1
          pause
        end if
        zave(jj1)=zave(jj1)+s(i)/fn
      end do
c Compute Hotelling statistic
      Q=0.d0
      do j=1,9
        do jj=1,9
          a=zave(j)-ave(j)
          b=zave(jj)-ave(jj)
          v=benfVARinv(j,jj)
          Q=Q+a*b*v
        end do
      end do
      Q=Q*fn
      RETURN
      END
C
C *******************************************************************
C
      SUBROUTINE TEST(jsimul,n,maxn,ncompmax,s,SCORE1,SCORE2,SCOREBIC,
     +    ncomps,XLRT1,XLRT2,XLRTBIC,ncompx)
C COMPUTES SCORE + LR TESTS WITH 1 COMPONENT, 2 COMPONENTS 
C AND BIC SELECTION OF COMPONENTS
C *ITERATIVE* BIC SELECTION ALGORITHM: N = SMALLEST # OF COMPONENTS SUCH 
C THAT THE OBJECTIVE FUNCTION INCREASES ==> see Bogdan et al. (AISM, 2002)
C                                           Eq. (2.6): TO BE PREFERRED 
C (the value of ncompmax is not essential == see alternative version *TEST2* 
C  for  * ncompmax = fixed # of iterations *)
C 
      IMPLICIT NONE
      integer*4 i,j,jj,jjj,ij,jsimul,n,maxn,ncompmax,ncomps,ncompx
      integer iscore,ixlrt,istop,niter
      real*8 s(maxn),fn,SCORE1,SCORE2,SCOREBIC,XLRT1,XLRT2,XLRTBIC
      real*8 funs,funsmax,funx,funxmax,score,scorej,xlrt(ncompmax)
      real*8 pi,arg,sum1,sum2,ss,d,xj,dj,dij  
      real*8 eNNTS,r,c1,rr,check,tol,gss
      complex*16 expNNTS,sumNNTS 
      complex*16 evec(maxn,ncompmax+1),cold(ncompmax+1),
     + coldc(ncompmax+1),cnew(ncompmax+1),cml(ncompmax+1),
     + cmlc(ncompmax+1)
      complex*16 frac,den
c      
      pi=3.14159265358979d0
      d=2.d0
      tol=(10.d0)**(-9)
      fn=dfloat(n)     
      score=0.d0
      funsmax=-1.d0/tol
      funxmax=-1.d0/tol
      iscore=0
      ixlrt=0
      SCORE1=0.d0
      SCORE2=0.d0
      SCOREBIC=0.d0
      XLRT1=0.d0
      XLRT2=0.d0
      XLRTBIC=0.d0
      ncomps=0
      ncompx=0
c
      cold(1)=dcmplx(1.d0,0.d0)
      cml(1)=cold(1)
      do j=1,ncompmax
        cold(j+1)=dcmplx(0.d0,0.d0)
        cml(j+1)=cold(j+1)
      end do
      r=(cdabs(cold(1)))**2
c
      do j=1,ncompmax
c   
c SCORE computations for # of components = j 
        sum1=0.d0
        sum2=0.d0  
c Penalty as in Bogdan et al: # of parameters = xj = 2*j = d*dj
c computation of d*dj is also performed in arg=d*pi*dj*ss for clarity 
        dj=dfloat(j)
        xj=d*dj
        do i=1,n
          ss=dlog10(s(i))
          arg=d*pi*dj*ss
          sum1=sum1+dcos(arg)
          sum2=sum2+dsin(arg)
c preliminary computations for ML estimation and LR Test
c exp(-ix): computational alternative 1 (straightforward)
          evec(i,1)=dcmplx(1.d0,0.d0)
          eNNTS=-2.d0*pi*dj*ss
          expNNTS=dcmplx(0.d0,eNNTS)
          sumNNTS=cdexp(expNNTS)
          evec(i,j+1)=sumNNTS          
!c exp(ix): computational alternative 2 (using Euler's formula)
!          evec(i,1)=dcmplx(1.d0,0.d0)
!          eNNTS=2.d0*pi*dj*ss
!          expNNTS=dcmplx(dcos(eNNTS),dsin(eNNTS))
!c if exp(-ix) is needed
!          sumNNTS=DCONJG(expNNTS)                              
!          evec(i,j+1)=sumNNTS       
          cold(j+1)=cold(j+1)+evec(i,j+1)/fn
        end do
        sum1=sum1/fn
        sum2=sum2/fn
        scorej=d*fn*(sum1*sum1+sum2*sum2)
        score=score+scorej        
        funs=score-xj*dlog(fn)
        if(j.eq.1) SCORE1=score
        if(j.eq.2) SCORE2=score   
c skip if iscore=1 
c (the statement is placed here because LRT computation must be performed)
        if(iscore.eq.1) goto 9991         
        if(funs.lt.funsmax) then
          ncomps=j-1
          SCOREBIC=score-scorej
          iscore=1
        else 
          ncomps=j
          funsmax=funs
          SCOREBIC=score
        end if                  
c      
9991    continue
c LRT computations for # of components = j ==> skip if ixlrt=1
        if(ixlrt.eq.1) goto 9992
c Maximum Likelihood estimation        
        r=r+(cdabs(cold(j+1)))**2
        do jj=1,j+1
          cold(jj)=cold(jj)/dsqrt(r)
          coldc(jj)=DCONJG(cold(jj))
          cnew(jj)=dcmplx(0.d0,0.d0)       
        end do  
        niter=1     
7001    continue
        do i=1,n
          den=dcmplx(0.d0,0.d0)
          do jj=1,j+1
            den=den+coldc(jj)*evec(i,jj)            
          end do
          do jj=1,j+1
            frac=evec(i,jj)/den
            cnew(jj)=cnew(jj)+frac/fn 
          end do
        end do
        c1=dble(cnew(1))
        cnew(1)=dcmplx(c1,0.d0)     
        rr=0.d0
        do jj=1,j+1
          rr=rr+(cdabs(cnew(jj)))**2
        end do
        check=0.d0
        do jj=1,j+1
          cnew(jj)=cnew(jj)/dsqrt(rr)
          check=check+cdabs(cold(jj)-cnew(jj))
        end do 
        if(check.lt.tol) then
          do jj=1,j+1
            cml(jj)=cnew(jj)
          end do
        else
          niter=niter+1
          do jj=1,j+1
            cold(jj)=cnew(jj)
            coldc(jj)=DCONJG(cold(jj))
          end do
          goto 7001
        end if
c Likelihood computations ==> to improve the algorithm some calculations
c                             may be anticipated in previous loops
c calculations in terms of log-likelihood        
!        xlrt(j)=1.d0
        xlrt(j)=0.d0
        do i=1,n
          ss=dlog10(s(i))
          gss=0.d0
          do jj=1,j+1
            do jjj=1,j+1
              ij=(jj-1)-(jjj-1)
c             ij=jj-jjj
              dij=dfloat(ij)
              eNNTS=2.d0*pi*dij*ss
              expNNTS=dcmplx(0.d0,eNNTS)
              sumNNTS=cdexp(expNNTS)
              cmlc(jjj)=DCONJG(cml(jjj))
              gss=gss+cml(jj)*cmlc(jjj)*sumNNTS
            end do
          end do       
!          xlrt(j)=xlrt(j)/gss
          xlrt(j)=xlrt(j)+dlog(gss)
        end do     
!        xlrt(j)=-2.d0*dlog(xlrt(j))
        xlrt(j)=2.d0*xlrt(j)
        funx=xlrt(j)-xj*dlog(fn)
        if(j.eq.1) XLRT1=xlrt(j)
        if(j.eq.2) XLRT2=xlrt(j) 
        if(funx.lt.funxmax) then
          ncompx=j-1
          XLRTBIC=xlrt(j-1)
          ixlrt=1
        else 
          ncompx=j
          funxmax=funx
          XLRTBIC=xlrt(j)
        end if            
c
9992    continue       
c Skip out of the loop over the # of components if:
c * both SCOREBIC and XLRT are selected (iscore=1 and ixlrt=1) 
c * both SCORE2 and XLRT2 are computed (j>2)
        istop=iscore*ixlrt*j
        if(istop.gt.2) goto 9999
c         
      end do 
c           
9999  continue
      RETURN
      END          
C
C *******************************************************************
C
      SUBROUTINE TEST2(jsimul,n,maxn,ncompmax,s,SCORE1,SCORE2,SCOREBIC,
     +    ncomps,XLRT1,XLRT2,XLRTBIC,ncompx)
C NEW 24/11/2011: COMPUTES SCORE + LR TESTS WITH 1 COMPONENT, 2 COMPONENTS 
C AND BIC SELECTION OF COMPONENTS
C *NON ITERATIVE* BIC SELECTION ALGORITHM: N = SMALLEST # OF COMPONENTS SUCH 
C THAT THE OBJECTIVE FUNCTION IS MAXIMIZED ==> see Bogdan et al. (AISM, 2002)
C                                              §5 where argmax is used
C (ncompmax = fixed # of iterations == see alternative version *TEST* 
C  for an iterative algorithm in which ncompmax is not essential)
C 
      IMPLICIT NONE
      integer*4 i,j,jsimul,n,maxn,ncompmax,ncomps,ncompx
      real*8 s(maxn),fn,SCORE1,SCORE2,SCOREBIC,XLRT1,XLRT2,XLRTBIC
      real*8 funs,funsmax,funx,funxmax,score,scorej
      real*8 pi,arg,sum1,sum2,ss,d,xj,dj     
      pi=3.14159265358979d0
      d=2.d0
      fn=dfloat(n)     
      score=0.d0
      funsmax=-1.d0*dfloat(maxn)
      SCORE1=0.d0
      SCORE2=0.d0
      SCOREBIC=0.d0
      XLRT1=0.d0
      XLRT2=0.d0
      XLRTBIC=0.d0
      ncomps=0
      ncompx=0
c
      do j=1,ncompmax
c   
c SCORE computations 
        sum1=0.d0
        sum2=0.d0  
c Penalty as in Bogdan et al: # of parameters = xj = 2*j = d*dj
c computation of d*dj is also performed in arg=d*pi*dj*ss for clarity 
        dj=dfloat(j)
        xj=d*dj
        xj=d*dfloat(j)
        do i=1,n
          ss=dlog10(s(i))
          arg=d*pi*dj*ss
          sum1=sum1+dcos(arg)
          sum2=sum2+dsin(arg)
        end do
        sum1=sum1/fn
        sum2=sum2/fn
        scorej=d*fn*(sum1*sum1+sum2*sum2)
        score=score+scorej        
        funs=score-xj*dlog(fn)
        if(j.eq.1) SCORE1=score
        if(j.eq.2) SCORE2=score   
        if(funs.gt.funsmax) then
          ncomps=j
          SCOREBIC=score
          funsmax=funs
        end if                
c
c LRT computations for # of components = j ==> skip if ixlrt=1
        XLRT1=1.d0
        XLRT2=1.D0
        XLRTBIC=1.d0
        ncompx=1 
c         
      end do 
c           
      RETURN
      END
C
C *******************************************************************
C
      DOUBLE PRECISION FUNCTION BENFCDF(X)
C RETURNS THE FUNCTION LOG10(X) FOR X in [1,10)
      IMPLICIT NONE
      real*8 X,a,b
      a=1.d0
      b=10.d0
      if(X.lt.a) then
        BENFCDF=0.d0
        write(*,*) "Warning: input significand LT 1"
      else if(X.ge.b) then
        BENFCDF=1.d0
        write(*,*) "Warning: input significand GE 10"
      else
        BENFCDF=DLOG10(X)
      end if
      RETURN
      END
C
C *******************************************************************
C
      DOUBLE PRECISION FUNCTION BENFCDFFR(U)
C ZZZZZ - TO DO
C                                          
      IMPLICIT NONE
      real*8 U,a,b,x,xj,cum
      integer j
      a=0.d0
      b=1.d0
      if(U.lt.a) then
        BENFCDFFR=0.d0
        write(*,*) "Warning: input *fractional part* LT 0"
      else if(U.ge.b) then
        BENFCDFFR=1.d0
        write(*,*) "Warning: input *fractional part* GE 1"
      else
        BENFCDFFR=0.d0
      end if
      RETURN
      END
C
C *******************************************************************
C
      DOUBLE PRECISION FUNCTION BETAPRIME(beta,alpha)
C TRANSFORMS THE SKEWNESS PARAMETER BETA OF A STABLE DISTRIBUTION (WITH 
C SHAPE PARAMETER ALPHA) INTO BETAPRIME AS REQUIRED BY ROUTINE DRNSTA 
C (see Chambers et al., JASA 1976)
      IMPLICIT NONE
      real*8 alpha,beta,a,b,c1,c2,pi
      pi=3.14159265358979d0
      if(alpha.eq.1.d0) then
        betaprime=beta
        RETURN
      end if
      c1=pi*(1.d0-alpha)/2.d0
      a=-dtan(c1)
      c2=-1.d0*pi*beta*(1.d0-dabs(1.d0-alpha))/2.d0
      b=dtan(c2)
      betaprime=a*b
      RETURN
      END
C
C *******************************************************************
C
      DOUBLE PRECISION FUNCTION TAILPROB(alpha,beta,sig,xtail)
C COMPUTES THE NOMINAL TAIL PROBABILITY FOR TAIL VALUE *XTAIL* IN A 
C STABLE DISTRIBUTION WITH SHAPE PARAMETER ALPHA, SKEWNESS PARAMETER BETA 
C AND SCALE PARAMETER SIG (see Rachev et al., p. 61)
      IMPLICIT NONE
      real*8 alpha,beta,sig,xtail,a1,a2,b,c1,c2,c,gg,p1,p2,prob,pi
      real*8 dgamma
      pi=3.14159265358979d0
      if(alpha.eq.1.d0) then
        c=2.d0/pi
      else
        c1=1.d0-alpha
        c2=dcos(pi*alpha/2.d0)
        gg=dgamma(2.d0-alpha)
        c=c1/(gg*c2)
      end if
      a1=(1.d0+beta)/2.0
      a2=(1.d0-beta)/2.0
      b=sig**alpha
      p1=c*a1*b
      p2=c*a2*b
      prob=(p1+p2)/(xtail**alpha)
      tailprob=prob
      RETURN
      END 
C
C *******************************************************************
C
      subroutine positstab(alpha,rn)
C SPECIFIC SIMULATION ALGORITHM FOR A *POSITIVE STABLE* 
C DISTRIBUTION - KANTER ALGORITHM AS ADAPTED BY XXXX
C (ignores betaprime = beta = xscal)
      IMPLICIT NONE
      real*8 alpha,rn,pi,rnu,rne,theta,rne2,drnunf,c1,c2,c3,c4
      pi=3.14159265358979d0
      rnu=DRNUNF()
      theta=rnu*pi
      call drnexp(1,rne)
      rne2=rne**(1.d0-alpha)
      c1=dsin(alpha*theta)
      c1=c1**alpha
      c2=dsin((1.d0-alpha)*theta)
      c2=c2**(1.d0-alpha)
      c3=dsin(theta)
      c4=c1*c2*(1.d0/c3)*(1.d0/rne2)
      rn=c4**(1.d0/alpha)
      RETURN
      END
C
C *******************************************************************
C
      subroutine positstab2(alpha,rn)
C SPECIFIC SIMULATION ALGORITHM FOR A *POSITIVE STABLE* 
C DISTRIBUTION - KANTER ALGORITHM AS ADAPTED BY Chambers et al. (1976)
C (ignores betaprime = beta = xscal)
      IMPLICIT NONE
      real*8 alpha,rn,pi,rnu,rne,theta,drnunf,a1,a2,a3,a,rat
      pi=3.14159265358979d0
      rnu=DRNUNF()
      theta=rnu*pi
      call drnexp(1,rne)
      a1=dsin((1.d0-alpha)*theta)
      a2=(dsin(alpha*theta))**(alpha/(1.d0-alpha))
      a3=(dsin(theta))**(1.d0/(1.d0-alpha))
      a=a1*a2/a3
      rat=a/rne
      rn=rat**((1.d0-alpha)/alpha)
      RETURN
      END     
C
C *******************************************************************
C
      subroutine PVK(VK0,pvalVKas,n)
C COMPUTES THE ASYMPTOTIC P-VALUE OF KUIPER TEST: SEE
C * MARDIA (1972), p. 178
C * JAMMALAMADAKA & SENGUPTA (2001), p. 155
C The series is approximated with nterm=100 terms
C VK0 = input value of the Kuiper statistic
C VK  = VK0*dsqrt(n) for asymptotic approximation of the distribution function
      IMPLICIT NONE
      real*8 VK0,VK,pvalVKas
      real*8 fn,fi,a,b,sum1,sum2,t,f
      integer n,nterm,i      
c      
      nterm=100
      fn=dfloat(n)
      VK=VK0*dsqrt(fn)
      t=2.d0
      f=4.d0
      sum1=0.d0
      sum2=0.d0
      do i=1,nterm
        fi=dfloat(i)
        a=f*fi*fi*VK*VK-1.d0
        b=dexp(-t*fi*fi*VK*VK)
        sum1=sum1+t*a*b
        a=a-t
        sum2=sum2+fi*fi*a*b
      end do
      a=t*f*VK
      b=3.d0*dsqrt(fn)
      pvalVKas=sum1-sum2*(a/b)
      if(pvalVKas.lt.0d0) then
        pvalVKas=0.d0
      else if(pvalVKas.gt.1.d0) then
        pvalVKas=1.d0
      end if      
      RETURN
      END
