Module RPtools

Use BranchingRatios
Use Control
Use EplusEminusProduction
Use Experiment
Use LoopCouplings
Use LoopMasses
Use MSSM_data
Use RP_data
Use StandardModel
Use SugraRuns

 Interface Fit_Neutrino_Data
  Module Procedure Fit_Neutrino_Data_sp
 End Interface

 Interface ParMin
  Module Procedure ParMin_sp
 End Interface

Contains


 Subroutine Calculate_Bi(mu, eps, vevL, vevSM, gp, g, M2L, B)
 !----------------------------------------------------------------------
 ! calculates the R-parity violating B-parameters [GeV^2]
 ! written by Werner Porod, 18.08.05
 !----------------------------------------------------------------------
  Implicit None
   Complex(dp), Intent(in) :: mu, eps(3) ! bilinear parameters of Superpotential
   Complex(dp), Intent(in) :: M2L(3,3)   ! left slepton masses
   Real(dp), Intent(in) :: gp, g         ! U(1) & SU(2) gauge couplings
   Real(dp), Intent(in) :: vevSM(2)      ! (v_d, v_u)
   Real(dp), Intent(in) :: vevL(3)       ! (v_1, v_2, v_3)
   Complex(dp), Intent(out) :: B(3)      ! B-parameters [GeV^2]

   Integer :: i1, i2
   Real(dp) :: DTerm
   Complex(dp) :: SumC

   DTerm = 0.125_dp * (g**2+gp**2)  &
       & * ( vevSM(1)**2 - vevSM(2)**2 + Dot_product(vevL,vevL))

   Do i1=1,3
    sumC = - (M2L(i1,i1)+Dterm) * vevL(i1) + vevSM(1) * Conjg(mu) * eps(i1)
    Do i2=1,3
     sumC = sumC - Conjg( eps(i1) ) * eps(i2) * vevL(i2)
     If (i1.Ne.i2) sumC = sumC - vevL(i2) * M2_L(i2,i1)
    End Do
    B(i1) = Real( sumC,dp ) / vevSM(2) 
  End Do
   
 End Subroutine Calculate_Bi


 Subroutine Calculate_RP_Parameters(delta, eps, vevL, B_i, Lambda, RSpm, RS0 &
                                 & , RP0, U, V, mN7, N7, kont)
 Implicit None
  Real(dp), Intent(in) :: delta     ! numerical precision for RGE running
  Integer , Intent (inout) :: kont  ! checks if everyhint is o.k. = 0

  Real(dp), Intent(out) :: vevL(3)      ! sneutrino vevs
  Complex(dp), Intent(out) :: eps(3)    ! superpotential parameters epsilon_i
  Complex(dp), Intent(out) :: B_i(3)    ! soft SUSY parameters B_i [GeV^2]
  Real(dp), Intent(out) :: Lambda(3)    ! RP lambda vector = mu v_i + eps_i v_d
  Complex(dp), Intent(out) :: RSpm(8,8) ! charged scalar mixing matrix
  Real(dp), Intent(out) :: RS0(5,5)     ! neutral scalar mixing matrix
  Real(dp), Intent(out) :: RP0(5,5)     ! neutral pseudoscalar mixing matrix
  Complex(dp), Intent(out) :: U(5,5), V(5,5) ! chargino mixing matrices
  Real(dp), Intent(out) :: mN7(7)    ! loop corrected neutrino/neutralino masses
  Complex(dp), Intent(out) :: N7(7,7)   ! neutralino/neutrino mixing matrix

  Integer :: count, isol
  Real(dp) :: g0(213), mudim, tz, dt, gauge_mZ(3), M2H_mZ(3)
  Complex(dp), Dimension(3,3) :: y_l_mZ,  y_d_mZ, y_u_mZ , Al_mZ, Ad_mZ &
           &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ
  Complex(dp) :: Mi_mZ(3), B_mZ, B_4(4), bi(4)
 ! neutrino constraints
  Real(dp) :: Lam_Sq, m_sq, m2_atm_rp, m2_sol_rp, tan2_sol, tan2_atm &
   & , Ue32, tan2_sol_opt, tan2_atm_opt, epsT12, epsT22
  Logical :: check
 ! RP-masses + mixing
  Real(dp) :: mGlu, mC(5), mC2(5), mSdown(6), mSdown2(6), mSup(6), mSup2(6)   &
    & , mP0(5), mP02(5), mS0(5), mS02(5), mSpm(8), mSpm2(8) &
    & , mN(7), mN1L(7), mN2(7), mN1L2(7)
  Complex(dp) :: PhaseGlu, RSdown(6,6), RSup(6,6), N(7,7), N1L(7,7), lam(3), epst(3)

 !-------------------------------------------------
 ! evolve the parameters down to m_Z if necessary 
 !-------------------------------------------------
  mudim = GetRenormalizationScale()
  Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
                    &,M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g0)

  If (mudim.Ne.mZ2) Then
 
   tz = 0.5_dp * Log(mZ2/mudim)
   dt = tz / 100._dp
   g0(1) = Sqrt(5._dp / 3._dp ) * g0(1)
 
   Call odeint(g0, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)
   g0(1) = Sqrt(3._dp / 5._dp ) * g0(1)
 
  End If
  
  Call GToParameters(g0,gauge_mZ, y_l_mZ,  y_d_mZ, y_u_mZ, Mi_mZ, Al_mZ, Ad_mZ &
          &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ, M2H_mZ, mu_mZ, B_mZ)

  kont = 0

  !--------------------------------------------
  ! try to find a consist set of RP-parameters
  !--------------------------------------------
  count = 0
  !-------------------------------------------------
  ! first set of masses to start with
  ! needed to get a first estimate of RP parameters
  !-------------------------------------------------
  tan2_atm_opt = 0.5_dp * (tan2_atm_min + tan2_atm_max)
  tan2_sol_opt = 0.5_dp * (tan2_sol_min + tan2_sol_max)

  eps = 1.e-4_dp * abs(mu_mZ)
  eps(2) = -eps(2)

  bi(1) = mu_mZ
  bi(2:4) = eps
  b_4 = B_mZ

  Lam_sq = 4._dp * Sqrt(m2_atm)                                               &
         & * Real(Mi_mZ(1)*Mi_mZ(2)*mu_mZ**2                                  &
         &    -0.5*vevSM(1)*vevSM(2)*Mi_mZ(2)*mu_mZ*(g0(1)**2+g0(2)**2),dp)   &
         &       / (g0(1)**2*Real(Mi_mZ(2),dp) + g0(2)**2*Real(Mi_mZ(1),dp))
  Lam_Sq = Abs(Lam_Sq)
  Lambda(2:3) = Sqrt(0.4995*Lam_sq)
  Lambda(1) = Sqrt(Lam_sq-2._dp*Lambda(2)**2)

  vevL = (Lambda - eps*vevSM(1)) / mu_mZ
  !--------------------------------------
  ! tree level masses
  !--------------------------------------
  Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4) )

  Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                 &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ   &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N              &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup             &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm      &
                    &, GenerationMixing, kont)
  !----------------------------
  ! approx. squared solar mass
  !-----------------------------
   epsT12 = ( eps(1) * (lambda(2)**2 + lambda(3)**2)                 &
          & - lambda(1) * (lambda(2)*eps(2)+lambda(3)*eps(3)) )**2   &
          & / ( (lambda(2)**2 + lambda(3)**2)                        &
          &   *(lambda(1)**2+lambda(2)**2 + lambda(3)**2))
   epsT22 = (lambda(3)*eps(2)-eps(3)*lambda(2))**2  &
          & /  (lambda(2)**2 + lambda(3)**2)

   m_Sq = 2._dp * oo16pi2 * (epsT12+epsT22) / Abs(mu_mZ)**2                   &
        &       * ( 3._dp*mf_d_mZ(3)*Rsdown(5,5)*Rsdown(5,6)*Y_d(3,3)**2      &
        &                *Log(mSdown2(6)/mSdown2(5))                          &
        &         + mf_l_mZ(3)* Rslepton(5,5)*Rslepton(5,6)*Y_l(3,3)**2       &
        &                     * Log(Slepton(6)%m2/Slepton(5)%m2))

   m_Sq = m_Sq**2

   If (m_Sq.Lt.m2_sol_min) eps = (m2_sol/m_sq)**0.25_dp * eps
   If (m_Sq.Gt.m2_sol_max) eps = (m2_sol/m_sq)**0.25_dp * eps

   bi(2:4) = eps
   vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
   !-------------------------------------------------------------
   ! recalculating B_i to fulfill tadpoles
   !-------------------------------------------------------------
   Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4) )

  !--------------------------------------
  !recalculation of tree level masses
  !--------------------------------------
   Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                    &, GenerationMixing, kont)
  !-----------------------------------------------------
  ! the iteration to get a consistent set of parameters
  !-----------------------------------------------------
  isol=100
  Do count = 0,isol

   Call NeutralinoMass_Loop_RP(g0(1), g0(2), Y_d_mZ, Y_l_mZ, Y_u_mZ, vevSM     &
          & , vevL, Mi_mZ(1), Mi_mZ(2), mu_mZ, eps, mC, mC2, U, V, mSup2, RSup &
          & , mSdown2, RSdown, mS02, RS0, mP02, RP0, mSpm2, RSpm, uD_L, uD_R   &
          & , uU_L, uU_R, mN, mN2, N, mN1L, mN1L2, N1L, kont, .False.)

   m2_sol_rp = mN1L2(2)-mN1L2(1)
   m2_atm_rp = mN1L2(3)-mN1L2(1)

   !------------------------------------------------
   ! checking experimental data, first the masses
   !------------------------------------------------
   check = .True. 
   If ((m2_atm_rp.Lt.m2_atm_min).or.(m2_atm_rp.Gt.m2_atm_max)) Then
    check = .False.
    Lambda = (m2_atm/m2_atm_rp)**0.25_dp * Lambda
   End If

   If ((m2_sol_rp.Lt.m2_sol_min).or.(m2_sol_rp.Gt.m2_sol_max)) Then
    check = .False.
    eps = (m2_sol/m2_sol_rp)**0.25_dp * eps
    bi(2:4) = eps
   End If

   If (.Not.check) Then
    vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                     &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                     &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                     &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                     &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                     &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                     &, GenerationMixing, kont)
    Cycle  !stop here and try the next iteration
   End If
   !-------------------------------------
   ! and now atmospheric and solar angle
   !-------------------------------------
   tan2_atm = Abs(N1L(3,6)/N1L(3,7))**2
   tan2_sol = Abs(N1L(2,5)/N1L(1,5))**2

   If ((tan2_atm.Lt.tan2_atm_min).or.(tan2_atm.Gt.tan2_atm_max)) Then
    check = .False.
    Lambda(2) = Sqrt(tan2_atm/(tan2_atm+tan2_atm_opt)) * Lambda(2)
    Lambda(3) = Sqrt(tan2_atm_opt/(tan2_atm+tan2_atm_opt)) * Lambda(3)
   End If

   If ((tan2_sol.Lt.tan2_sol_min).or.(tan2_sol.Gt.tan2_sol_max)) Then
    check = .False.

    lam = Cmplx(lambda,0._dp,dp)
    Call CalculateEpsTilde_from_Eps(epsT,eps,lam)

    epsT(1) = Sqrt(tan2_sol_opt/(tan2_sol+tan2_sol_opt)) * epsT(1)
    epsT(2) = Sqrt(tan2_sol/(tan2_sol+tan2_sol_opt)) * epsT(2)
    Call CalculateEps_from_EpsTilde(eps,epsT,lam)

    bi(2:4) = eps
   End If

   If (.Not.check) Then
    vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                     &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                     &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                     &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                     &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                     &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                     &, GenerationMixing, kont)
    Cycle
   End If
   !-----------------------------------
   ! and now the reactor constraint
   !-----------------------------------
   Ue32 = Abs(N1l(3,5))**2

   If (Ue32.Gt.Ue32_max) Then
    check = .False.
    Lambda(1) = Lambda(1) / 2._dp
    Lam_Sq = Dot_product(Lambda, lambda)

    vevL = (Lambda - eps*vevSM(1)) / mu_mZ 
    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, B_4(2:4))
   End If
   !-----------------------------------
   ! leave loop if everything is fine
   !-----------------------------------
   If (check) Exit
   !---------------------------------------------------
   ! else recalculate tree level masses
   !---------------------------------------------------
   Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3) &
                    &, bi, B_4, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ                  &
                    &, M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ    &
                    &, mGlu, PhaseGlu, mC, mC2, U, V, mN, mN2, N               &
                    &, mSdown, mSdown2, RSdown, mSup, mSup2, RSup              &
                    &, mP0, mP02, RP0, mS0, mS02, RS0, mSpm, mSpm2, RSpm       &
                    &, GenerationMixing, kont)

  End Do

  If(.Not.check)Then 
    Write(ErrCan,*)'Error, no solution found',isol
!    Write(*,*)'Error, no solution found',isol
!  Else
!    Write(*,*)'solution found after:',count,Sqrt(Sum(Abs(eps)**2))/Real(mu_mZ,dp) 
  End If

  mN7 = mN1L
  N7 = N1L
  B_i = B_4(2:4)

 End Subroutine Calculate_RP_Parameters

 Subroutine CalculateEpsTilde_from_Eps(epsT, eps, Lambda)
 Implicit None
  Complex(dp), Intent(in), Dimension(3) :: eps, Lambda 
  Complex(dp), Intent(out), Dimension(3) :: epsT

  Real(dp) :: AbsLam2, Abs23, AbsLam22, Abs232


  Abs232 = Abs(Lambda(2))**2 + Abs(Lambda(3))**2
  AbsLam22 = Abs232 + Abs(Lambda(1))**2

  Abs23 = Sqrt(Abs232)
  AbsLam2 = Sqrt(AbsLam22)

  epsT(1) = ( eps(1) * Abs232                                   &
          & - Lambda(1) * (Lambda(2)*eps(2)+Lambda(3)*eps(3))   &
          & ) / (Abs23*AbsLam2)
  epsT(2) = (Lambda(3)*eps(2)-Lambda(2)*eps(3)) / Abs23
  epsT(3) = (Lambda(1)*eps(1)+Lambda(2)*eps(2)+Lambda(3)*eps(3)) / AbsLam2

 End Subroutine CalculateEpsTilde_from_Eps

 Subroutine CalculateEps_from_EpsTilde(eps, epsT, Lambda)
 Implicit None
  Complex(dp), Intent(in), Dimension(3) :: epsT, Lambda 
  Complex(dp), Intent(out), Dimension(3) :: eps

  Real(dp) :: AbsLam2, Abs23, AbsLam22, Abs232

  Abs232 = Abs(Lambda(2))**2 + Abs(Lambda(3))**2
  AbsLam22 = Abs232 + Abs(Lambda(1))**2

  Abs23 = Sqrt(Abs232)
  AbsLam2 = Sqrt(AbsLam22)

  eps(1) = ( epsT(3) * Lambda(1) + epsT(1) * Abs23  )  / AbsLam2

  eps(2) = (-epsT(1) * Lambda(1) * Lambda(2) + epsT(3) * Lambda(2) * Abs23  &
         & + epsT(2) * Lambda(3) * AbsLam2 ) / (Abs23 * AbsLam2)

  eps(3) = (-epsT(1) * Lambda(1) * Lambda(3) + epsT(3) * Lambda(3) * Abs23  &
         & - epsT(2) * Lambda(2) * AbsLam2 ) / (Abs23 * AbsLam2)

 End Subroutine CalculateEps_from_EpsTilde

  
 Subroutine Fit_Neutrino_Data_sp(Nfit, m2_atm_min, m2_atm_max, tan2_atm_min   &
       & , tan2_atm_max, m2_sol_min, m2_sol_max, tan2_sol_min, tan2_sol_max   &
       & , Ue32_max, g, gp, rh0, hpncs, lmbd, vevSM, vR, vS, vP, mPhi, MR, mu &
       & , M1, M2, hnu, vL)
  !--------------------------------------------------------------------------
  ! adjustes h_nu and v_L such that neutrino data are most likely fullfilled
  ! based on the routine FitNeu by Martin Hirsch
  !     options:
  !     Nfit:   Comment:
  !       <0    returns without doing anything, added by Werner
  !       0     No fit at all, Eps and Lam Random numbers, see below
  !       1     (Lam=Atm,Eps=Sol,sign condition = yes),
  !              28% success rate with current NuData
  !
  !       2     (Lam=Atm,Eps=Sol,sign condition = no),
  !              12% success rate with current NuDat
  !
  !       3     (Lam=Sol,Eps=Atm,sign condition = yes),
  !              28% success rate with current NuData
  !
  !       4     (Lam=Sol,Eps=Atm,sign condition = no),
  !               9% success rate with current NuData
  !
  !     Note: Success rate can be enhanced by narrowing the random regions,
  !           but at the expense of not covering uniformly the (currently)
  !           allowed neutrino parameter space
  !--------------------------------------------------------------------------
  Implicit none
   Integer, intent(in) :: Nfit
   Real(dp), intent(in) :: m2_atm_min, m2_atm_max, tan2_atm_min, tan2_atm_max &
       & , m2_sol_min, m2_sol_max, tan2_sol_min, tan2_sol_max, Ue32_max       &
       & , g, gp, rh0, hpncs, lmbd, vevSM(2), vR, vS, vP, mPhi, MR
   Complex(dp), Intent(in) :: M1, M2, mu
   Real(dp), Dimension(3), Intent(out) :: hnu, vL

   Real(dp) :: epssq, lambda, eps(3), lam(3), mptot, mRtot, mphot, x1, x2, x3 &
       & , det7, c1, c2, c3, vu, vd, sgnc
   Real(dp) :: vec8(8), vec13(13) ! vectors of random numbers
   integer :: i1
   !----------------------------------------------
   ! no adjustment if Nfit < 0
   !----------------------------------------------
   if (Nfit.lt.0) return 

   Iname = Iname + 1
   NameOfUnit(Iname) = "Fit_Neutrino_Data_sp"

   !--------------------------------------------
   ! initialization
   !--------------------------------------------
   epssq = 0._dp
   lambda = 0._dp
   eps = 0._dp
   lam = 0._dp

   If (Nfit.gt.0) then
    vd = vevSM(1)
    vu = vevSM(2)
    mptot=mPhi + lmbd * vP * oosqrt2
    mRtot=MR+hpncs*vP * oosqrt2
    mphot=g**2*m1+gp**2*m2
    x1=8.0_dp*m1*m2*mu
    x1=x1*(mptot*mRtot*mu-hpncs**2*mu*vR*vS+rh0**2*mRtot*vd*vu)
    x2=4.0_dp*mu*vd*(mptot*mRtot-hpncs**2*vR*vS)*vu
    x3=rh0**2*mRtot*(vd**2+vu**2)**2
    det7=mRtot*(x1-mphot*(x2+x3))/8.0_dp
  !      det4=mphot/(4.0_dp*m1*m2*mu**2-2.0_dp*mphot*mu*vd*vu)
  !
  !1.2  === c1,c2 and c3 ===
  !
    c1=-hpncs**2*vR*vS*mu+mptot*mRtot*mu
    c1=(c1+rh0**2*mRtot*vd*vu)*mphot*mRtot/4.0_dp/mu
    c2=rh0*mphot*mRtot*(rh0*mRtot+hpncs*mu)*vu*(vu**2-vd**2)/8.0_dp/mu
    c3=(rh0*mRtot+hpncs*mu)**2*vu**2*(2.0_dp*m1*m2*mu-mphot*vd*vu)
    c3=c3/4.0_dp/mu
    c1=-c1/det7
    c2=-c2/det7
    c3=-c3/det7
   end if

   Select Case(Nfit)
   Case(0)
    Call Random_Number(vec8)
    lambda=10.0_dp**(1.0_dp - 4.0_dp*vec8(1))
    epssq =10.0_dp**(0.0_dp - 6.0_dp*vec8(2))
    lam = vec8(3:5)
    lam = lam * lambda / Sqrt(Dot_product(lam,lam)) 
    eps = vec8(6:8)
    eps = eps * Sqrt(epssq/Dot_product(eps,eps)) 

   Case(1)
    Call Random_Number(vec13)
    lambda = Sqrt(sqrt(m2_atm_min+(m2_atm_max-m2_atm_min)*vec13(1))/abs(c1))
    epssq =  abs(sqrt(m2_sol_min+(m2_sol_max-m2_sol_min)*vec13(2))/abs(c3))

    lam(1) = sqrt(Ue32_max * vec13(3))
    lam(2) = Sqrt(tan2_atm_min + (tan2_atm_max - tan2_atm_min) * vec13(4))
    lam(3) = 1._dp
    lam = lam * lambda / Sqrt(Dot_product(lam,lam))

    eps(1) = sqrt(tan2_sol_min + (tan2_sol_max - tan2_sol_min) * vec13(5))
    eps(2) = 1._dp
    eps(3) = sqrt(tan2_sol_min + 0.8_dp * (tan2_sol_max - tan2_sol_min) * vec13(6))
    eps = eps * Sqrt(epssq/Dot_product(eps,eps))

    Do i1=1,3
     If (vec13(6+i1).Lt.0.5_dp) lam(i1) = - lam(i1)
     If (vec13(9+i1).Lt.0.5_dp) eps(i1) = - eps(i1)
    end do

    sgnc=(eps(2)/eps(3))*(lam(2)/lam(3))
    If (sgnc.Gt.0._dp) lam(2) = - lam(2)

   Case(3)
    Call Random_Number(vec13)
    lambda = Sqrt(sqrt(m2_sol_min+(m2_sol_max-m2_sol_min)*vec13(2))/Abs(c1))
    epssq = Abs(sqrt(m2_atm_min+(m2_atm_max-m2_atm_min)*vec13(1))/c3)

    lam(1) = sqrt(tan2_sol_min + (tan2_sol_max - tan2_sol_min) * vec13(5))
    lam(2) = 1._dp
    lam(3) = sqrt(tan2_sol_min + 0.8_dp * (tan2_sol_max - tan2_sol_min) * vec13(6))
    lam = lam * lambda / Sqrt(Dot_product(lam,lam))

    eps(1) = sqrt(Ue32_max * vec13(3))
    eps(2) = sqrt(tan2_atm_min + (tan2_atm_max - tan2_atm_min) * vec13(4))
    eps(3) = 1._dp
    eps = eps * Sqrt(epssq/Dot_product(eps,eps))

    Do i1=1,3
     If (vec13(6+i1).Lt.0.5_dp) lam(i1) = - lam(i1)
     If (vec13(9+i1).Lt.0.5_dp) eps(i1) = - eps(i1)
    end do

    sgnc=(eps(2)/eps(3))*(lam(2)/lam(3))
    If (sgnc.Gt.0._dp) lam(2) = - lam(2)

   Case default
    Write(ErrCan,*) "Problem in routine "//NameOfUnit(Iname)
    Write(ErrCan,*) "Option Nfit=",Nfit," does not exist."
    If (ErrorLevel.Gt.-2) call TerminateProgram
   end select
    vL = (lam - eps*vd) / mu
    hnu = sqrt2 * eps / vR

    Iname = Iname - 1
 
 End Subroutine Fit_Neutrino_Data_sp


 Subroutine Model_bilinear_Rparity(add_Rparity, HighScaleModel, delta, epsI   &
       & , deltaM, ratioWoM, m32, grav_fac, CalcTBD, Ecms, Pm, Pp, ISR, Beam  &
       & , SigSup , SigSdown, SigC, SigChi0, SigS0, SigSP, SigHp, M_GUT, kont)
  Implicit None

  !-----------------------------------------------------------
  ! input / ouput
  !-----------------------------------------------------------
  Logical, Intent(in) :: add_Rparity
  Character(len=15), Intent(inout) :: HighScaleModel
  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: delta   ! required precision for spectrum calculation
  Real(dp), Intent(out) :: M_GUT  ! scale of SUSY boundary conditions,
                                  ! usually the GUT scale
 !-------------------------------
 ! widths and branching ratios
 !-------------------------------
  Real(dp), Intent(in) ::  m32, grav_fac, epsI, deltaM, ratioWoM
  Logical, Intent(in) :: CalcTBD
 !--------------------------------
 ! cross section calculation
 !--------------------------------
  Real(dp), Intent(in) :: Ecms(:), Pm(:), Pp(:)
  Logical, Intent(in) :: ISR(:), Beam(:)
  Real(dp), Intent(out) :: SigSup(:,:,:) , SigSdown(:,:,:), SigC(:,:,:)   &
     & , SigChi0(:,:,:), SigS0(:,:), SigSP(:,:,:), SigHp(:,:,:)

  !-----------------------------------------------------------
  ! local variables
  !-----------------------------------------------------------
  Integer :: i1, i_min(3)
  Real(dp) :: sinW2, gp, g,vev
  Real(dp) :: g0(213), mudim, tz, dt, gauge_mZ(3), M2H_mZ(3), mN1L(7), mN1L2(7)
  Real(dp) :: mGlu_T, mSdown_T(6), mSdown2_T(6), mSup_T(6), mSup2_T(6)
  Complex(dp), Dimension(3,3) :: y_l_mZ,  y_d_mZ, y_u_mZ , Al_mZ, Ad_mZ &
           &, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ
  Complex(dp) :: Mi_mZ(3), B_mZ, bi(4), BiEpsi(4), N1L(7,7), PhaseGlu_T &
           & , RSdown_T(6,6), RSup_T(6,6)
  !--------------------------------
  ! cross section calculation
  !--------------------------------
  Integer :: p_max

 !------------------------------------------------------
 ! internal information on particle identities
 !------------------------------------------------------
  Integer :: id_gl, id_ph, id_grav
  Integer, Dimension(1) :: id_Z, id_W
  Integer, Dimension(3) :: id_nu, id_l, id_d, id_u

  Iname = Iname + 1
  NameOfUnit(Iname) = "Model_bilinear_Rparity"

  kont = 0

  Call Initialize_RPexplicit(GenerationMixing, id_gl, id_ph, id_Z, id_W &
               & , id_nu, id_l, id_d, id_u, id_grav)

  If (add_Rparity) Then ! calculate parameters from a high scale model
                       ! assuming conserved R-parity

   Call FirstGuess(phase_mu, tanb, Mi, M2_E, M2_L, A_l, M2_D, M2_Q, M2_U, A_d &
           & , A_u, mu, B, M2_H, gp, g, Y_l, Y_d, Y_u, vevSM, mP02, mP0, kont)
   
   If (kont.Ne.0) Return

   Call TreeMasses(gp, g, vevSM, Mi(1), Mi(2), Mi(3), mu, B                &
        &, tanb, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d               &
        &, A_u, Y_d, Y_u, Glu%m, PhaseGlu, ChiPm%m, ChiPm%m2, U, V, Chi0%m &
        &, Chi0%m2, N, Sneut%m, Sneut%m2, Rsneut, Slepton%m, Slepton%m2    &
        &, RSlepton, Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup        &
        &, P0%m, P0%m2, RP0, S0%m, S0%m2, RS0, Spm%m, Spm%m2, RSpm         &
        &, GenerationMixing, kont, .False.) ! tree-level Higgs mass

   If (kont.Ne.0) Return

  !-------------------------------------------------------------------------
  ! Iterative procedure to get the 1-loop masses and mixing matrices;
  ! delta is the relative precision required for the masses after two runs
  ! The procedure is stopped after the i-th iteration 
  ! if delta > |m(i) - m(i-1)|/m(i) for all SUSY. In the case that more than
  ! 20 iterations are necessary the iteration loop is left and a warning is 
  ! given.
  !--------------------------------------------------------------------------
   ! In the SPA convention the the renormalization scale is fixed with 1 TeV
   If (SPA_Convention) Call SetRGEScale(1.e3_dp**2) ! 

   Call Sugra(delta, vevSM, ChiPm%m, U, V, Chi0%m, N, S0%m, S0%m2, RS0    &
     & , P0%m, P0%m2, RP0, Spm%m, Spm%m2, RSpm, Sdown%m, Sdown%m2, RSdown &
     & , Sup%m, Sup%m2, RSup, Slepton%m, Slepton%m2, RSlepton, Sneut%m    &
     & , Sneut%m2, RSneut, Glu%m, PhaseGlu, gauge, uL_L, uL_R, uD_L, uD_R &
     & , uU_L, uU_R, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D   &
     & , M2_Q, M2_U, M2_H, mu, B, m_GUT, kont, WriteOut, n_run)

   If (l_fit_RP_parameters) Then

    Call Calculate_RP_Parameters(delta_mass, eps, vevL, Beps, Lam_ex, RSpm8 &
                     & , RS05, RP05, U5, V5, Chi07%m, N7, kont)

   Else ! .not.l_fit_RP_parameters
   !-------------------------------------------------
   ! evolve the parameters down to m_Z if necessary 
   !-------------------------------------------------
    mudim = GetRenormalizationScale()
    Call ParametersToG(gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
                    &,M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g0)

    If (mudim.Ne.mZ2) Then
 
     tz = 0.5_dp * Log(mZ2/mudim)
     dt = tz / 100._dp
     g0(1) = Sqrt(5._dp / 3._dp ) * g0(1)
 
     Call odeint(g0, 213, 0._dp, tz, delta, dt, 0._dp, rge213, kont)
     g0(1) = Sqrt(3._dp / 5._dp ) * g0(1)
 
    End If
  
    Call GToParameters(g0, gauge_mZ, y_l_mZ,  y_d_mZ, y_u_mZ, Mi_mZ, Al_mZ  &
          & , Ad_mZ, Au_mZ,M2E_mZ, M2L_mZ, M2D_mZ, M2Q_mZ, M2U_mZ, M2H_mZ   &
          & , mu_mZ, B_mZ)

    Call Calculate_Bi(mu_mZ, eps, vevL, vevSM, g0(1), g0(2), M2L_mZ, Beps)
    BiEpsi(1) = B_mZ
    BiEpsi(2:4) = Beps
    bi(1) = mu_mZ
    bi(2:4) = eps
 
    Call TreeMassesEps3(g0(1), g0(2), vevSM, vevL, Mi_mZ(1), Mi_mZ(2)        &
           & , Mi_mZ(3), bi, BiEpsi, M2E_mZ, M2L_mZ, Al_mZ, Y_l_mZ           &
           & , M2D_mZ, M2U_mZ, M2Q_mZ, Ad_mZ, Au_mZ, Y_d_mZ, Y_u_mZ          &
           & , mGlu_T, PhaseGlu_T, ChiPm5%m, ChiPm5%m2, U5, V5, Chi07%m      &
           & , Chi07%m2, N7, mSdown_T, mSdown2_T, RSdown_T, mSup_T, mSup2_T  &
           & , RSup_T, P05%m, P05%m2, RP05, S05%m, S05%m2, RS05              &
           & , Spm8%m, Spm8%m2, RSpm8, GenerationMixing, kont)

    If (kont.Ne.0) Then
     Iname = Iname - 1
     Return
    End If
  
    Call NeutralinoMass_Loop_RP(g0(1), g0(2), Y_d_mZ, Y_l_mZ, Y_u_mZ, vevSM    &
         & , vevL, Mi_mZ(1), Mi_mZ(2), mu_mZ, eps, ChiPm5%m, ChiPm5%m2, U5, V5 &
         & , mSup2_T, RSup_T, mSdown2_T, RSdown_T, S05%m2, RS05, P05%m2, RP05  &
         & , Spm8%m2, RSpm8, uD_L, uD_R, uU_L, uU_R, Chi07%m, Chi07%m2, N7     &
         & , mN1L, mN1L2, N1L, kont, .False.)

    Chi07%m = mN1L
    Chi07%m2 = mN1L2
    N7 = N1L

   End If  ! l_fit_RP_parameters

   !-------------------------------------------------------------------
   ! replacing the tree level masses by the 1-loop masses of the MSSM
   ! assuming that the effect of the RP parameter is tiny
   ! this is justified in the region of parameter space where neutrino
   ! data are correctely explained
   !-------------------------------------------------------------------
    ChiPm5(1:3)%m = mf_l
    ChiPm5(4:5)%m = ChiPm%m
    Chi07(4:7)%m = Chi0%m

    P05(1)%m = mZ
    Do i1=2,5
     If ( (RP05(i1,1)**2 +RP05(i1,2)**2).Gt.0.5_dp) P05(i1)%m = P0(2)%m
     If ( RP05(i1,3)**2.Gt.0.5_dp) P05(i1)%m = Sneut(1)%m
     If ( RP05(i1,4)**2.Gt.0.5_dp) P05(i1)%m = Sneut(2)%m
     If ( RP05(i1,5)**2.Gt.0.5_dp) P05(i1)%m = Sneut(3)%m
    End Do
    P05%m2 = P05%m**2

    i_min = 0
    Do i1=1,5
     If ( (RS05(i1,1)**2 +RS05(i1,2)**2).Gt.0.5_dp) Then
      S05(i1)%m = S0(1 + i_min(1) )%m
      i_min(1) = i_min(1) + 1
     End If
     If ( RS05(i1,3)**2.Gt.0.5_dp) S05(i1)%m = Sneut(1)%m
     If ( RS05(i1,4)**2.Gt.0.5_dp) S05(i1)%m = Sneut(2)%m
     If ( RS05(i1,5)**2.Gt.0.5_dp) S05(i1)%m = Sneut(3)%m
    End Do
    S05%m2 = S05%m**2

    Spm8(1)%m = mW
    i_min = 0
    Do i1=2,8
     If ( (Abs(RSpm8(i1,1))**2 +Abs(RSpm8(i1,2))**2).Gt.0.5_dp) &
        &  Spm8(i1)%m = Spm(2)%m
     If ( (Abs(RSpm8(i1,3))**2 +Abs(RSpm8(i1,6))**2).Gt.0.5_dp) Then
      Spm8(i1)%m = Slepton(1 + i_min(1) )%m
      i_min(1) = i_min(1) + 1
     End If
     If ( (Abs(RSpm8(i1,4))**2 +Abs(RSpm8(i1,7))**2).Gt.0.5_dp) Then
      Spm8(i1)%m = Slepton(3 + i_min(2) )%m
      i_min(2) = i_min(2) + 1
     End If
     If ( (Abs(RSpm8(i1,5))**2 +Abs(RSpm8(i1,8))**2).Gt.0.5_dp) Then
      Spm8(i1)%m = Slepton(5 + i_min(3) )%m
      i_min(3) = i_min(3) + 1
     End If
    End Do
    Spm8%m2 = Spm8%m**2

   HighScaleModel = "RPexplicit"

  Else ! .not.add_Rparity
   sinW2 = 1._dp - mW2/mZ2

   alpha_mZ = Alpha_MSbar(mZ, mW)
   gauge(1) = Sqrt( 4._dp*pi*alpha_mZ/(1._dp-sinW2) )
   gauge(2) = Sqrt( 4._dp*pi*alpha_mZ/sinW2)
   gauge(3) = Sqrt( 4._dp*pi*alphas_mZ)
   gp = gauge(1)
   g = gauge(2)

   vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
   vevSM(1) = vev / Sqrt(1._dp + tanb**2)
   vevSM(2) = tanb * vevSM(1)

   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp

   Do i1=1,3
    Y_l(i1,i1) = sqrt2 * mf_l_mZ(i1) / vevSM(1)
    Y_d(i1,i1) = sqrt2 * mf_d_mZ(i1) / vevSM(1)
    Y_u(i1,i1) = sqrt2 * mf_u_mZ(i1) / vevSM(2)
   End Do
   If (GenerationMixing) Then
    i1 = GetYukawaScheme()
    If (i1.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM),Y_d) 
    End If
   End If

   bi(1) = mu
   bi(2:4) = eps
   BiEpsi(1) = B * mu
   BiEpsi(2:4) = Beps * eps 
   Call TreeMassesEps3(gp, g, vevSM, vevL, Mi(1), Mi(2), Mi(3), bi, BiEpsi    &
      & , M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d, A_u, Y_d, Y_u          &
      & , Glu%m, PhaseGlu, ChiPm5%m, ChiPm5%m2, U5, V5, Chi07%m, Chi07%m2, N7 &
      & , Sdown%m, Sdown%m2, RSdown, Sup%m, Sup%m2, RSup, P05%m, P05%m2, RP05 &
      & , S05%m, S05%m2, RS05, Spm8%m, Spm8%m2, RSpm8, GenerationMixing, kont)

  End If ! add_Rparity

   !----------------------------------------------
   ! reorder state identification if necessary
   !----------------------------------------------
   If (.Not.GenerationMixing) Then
    Call Swap_Order_Sf(RSdown(1,1), Sdown(1)%id, Sdown(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSdown(3,3), Sdown(3)%id, Sdown(4)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(1,1), Sup(1)%id, Sup(2)%id, id_p, c_name)    
    Call Swap_Order_Sf(RSup(3,3), Sup(3)%id, Sup(4)%id, id_p, c_name)    
   End If
  !-----------------------------------------------------------------------
  ! as there is a large hierachy between m_nu and m_chi, there might have
  ! been numerical problems. However, in all known cases this has been
  ! the mass of the lightest neutrino, which is of no importance for the
  ! following. Therefore, the error flag is reset to zero. 
  ! This is of course dangerous and has to be improved.
  !-----------------------------------------------------------------------
  If (kont.Eq.-14) kont = 0

  
  If ((L_BR).And.(kont.Eq.0)) Then

   Call CalculateBR(0, id_nu, 0, id_l, 3, id_d, 3, id_u, 1, id_Z, 1, id_W, 0   &
    & , 0, 6, 6, 7, 5, 1, 5, 5, 8, id_grav, id_gl, id_ph, gauge, Glu, PhaseGlu &
    & , ChiPm5, U5, V5, Chi07, N7, Sup, RSup, Sdown, RSdown, uD_L, uD_R, uU_L  &
    & , uU_R, S05, RS05, P05, RP05, Spm8, RSpm8, epsI, deltaM, CalcTBD         &
    & , ratioWoM, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, eps, vevSM, vevL, F_Gmsb   &
    & , m32, grav_fac)

  End If


 !---------------------------------------------------------------------------
 ! Calculation of the cross sections in e+ e- annihilation provided L_Cs is
 ! set .TRUE. (default) and that the routine Sugra has finished
 ! correctly (kont.eq.0) 
 ! The following input quantities can be specified in the file 
 ! CrossSections.in: Ecms .... c.m.s enerergy in GeV
 !                   Pm ...... degree of longitudinal polarization of incoming
 !                             electron
 !                   Pp ...... degree of longitudinal polarization of incoming
 !                             positron
 !                   ISR ..... if .TRUE. then the effect of initial state
 !                             radiation will be included
 ! In the case that the file CrossSections.in does not exist, the following
 ! default values are used: Ecms = 500 GeV, Pm = Pp = 0, ISR = .TRUE.
 !----------------------------------------------------------------------------
  If ((L_CS).And.(kont.Eq.0)) Then
   p_max = Size(Pm)

   Do i1=1,p_max
    If (Ecms(i1).Eq.0._dp) Exit
    Call CalculateCrossSections(Ecms(i1), Pm(i1), Pp(i1), ISR(i1), Beam(i1)    &
           & , "Tesla800", Sup%m, RSup, mf_u, Sdown%m, RSdown, mf_d, Glu%m     &
           & , SigSup(i1,:,:), SigSdown(i1,:,:), ChiPm5%m, U5, V5, Chi07%m, N7 &
           & , SigC(i1,:,:), SigChi0(i1,:,:), S05%m, RS05, vevSM, vevL         &
           & , P05%m, RP05, Spm8%m, RSpm8, SigS0(i1,:), SigSP(i1,:,:)          &
           & , SigHp(i1,:,:) )
   End Do

  End If

  lam_ex = mu_mZ * vevL + vevSM(1) * eps
  mf_nu = Chi07(1:3)%m

  Iname = Iname - 1

 Contains

  Subroutine Swap_Order_Sf(Rij, id1, id2, id_p, names)
  Implicit None
   Complex(dp), Intent(in) :: Rij
   Integer, Intent(in) :: id1, id2
   Integer, Intent(inout) :: id_p(:)
   Character(len=12), Intent(inout) :: names(:)

   Integer :: k(2)
   Character(len=12) :: nam(2)
   !--------------------------------------------------------------------
   ! changes name and particle id, in case that the lighter of two
   ! sfermions is a right-sfermion
   !--------------------------------------------------------------------
   If (Abs(Rij).Lt.0.5_dp) Then
    k = id_p(id1:id1+1)
    id_p(id1:id1+1) = id_p(id2:id2+1)
    id_p(id2:id2+1) = k

    nam = names(id1:id1+1)
    names(id1:id1+1) = names(id2:id2+1)
    names(id2:id2+1) = nam
   End If

  End Subroutine Swap_Order_Sf


 End Subroutine Model_bilinear_Rparity

 Subroutine ParMin_sp(g, gp, hn, rh0, h, lmbd, Ahn, Ah, RAh0, Almbd, vevSM, vL &
      & , vR, vS, vP, mu, MPhi, MR, B, bmr, BMphi, mH2, mL2, mR2, mS2, mP2)
 !------------------------------------------------------------------------
 ! this routine solves the tad-pol equations in terms of the soft susy masses
 ! squared in the spontaneous model.
 ! taken from J.Romao and adjusted for SPheno
 !------------------------------------------------------------------------
 Implicit None

   Real(dp), Intent(in) :: g, gp, rh0, h, lmbd, Ah, RAh0, Almbd, vevSM(2) &
      & , vL(3) , vR, vS, vP, MPhi, MR, bmr, BMphi, hn(3), Ahn(3)
   Complex(dp), Intent(in) :: B, mu
   Real(dp), Intent(out) :: mH2(2), mL2(3), mR2, mS2, mP2

  Real(dp) :: g2, gp2, vd, vu, t1, t2, t3

  Iname = Iname + 1
  NameOfUnit(Iname) = "ParMin_sp"

  g2=g*g
  gp2=gp*gp
  vd = vevSM(1)
  vu = vevSM(2)
  
  mH2(1)=(B*vu+rh0*vu*( (lmbd*vp**2)/4._dp+(h*vR*vS)/2._dp        &
    & -(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)+(RAh0*vp*vu)/sqrt2           &
    & +(gp2*vd*(-vd**2/2._dp+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp  &
    & -vL(3)**2/2._dp))/4._dp-(- mu -(rh0*vp)/sqrt2)*sqrt2*(-(rh0*vd   &
    & *vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1)*vL(1))/sqrt2         &
    & +(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3))/sqrt2))/sqrt2)-(g2     &
    & *(2*vd*sqrt2*(vd**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2 &
    & /2._dp)-vd*sqrt2*(vd**2/2._dp+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2 &
    & /2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vd

  t1=1/vu
  t2=B*vd-(vR**2*vu*hn(1)**2)/2._dp-(vR**2*vu*hn(2)**2)/2._dp&
    &-(vR**2*vu*hn(3)**2)/2._dp+rh0*vd*( (lmbd*vp**2)/4._dp+(h&
    &*vR*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)
  t3=(RAh0*vd*vp)/sqrt2-(- mu -(rh0*vp)/sqrt2)*(-(rh0*vp*vu) &
    /2._dp-( mu *vu)/sqrt2)*sqrt2-vR*((Ahn(1)*vL(1))/sqrt2&
    &+(Ahn(2)*vL(2))/sqrt2+(Ahn(3)*vL(3))/sqrt2)-sqrt2*((h&
    &*vp*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu&
    &*hn(2)*vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)*((hn(1)*vL(1)) &
    /sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3))/sqrt2)-(gp2&
    &*vu*(-vd**2/2._dp+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp&
    &-vL(3)**2/2._dp))/4._dp-(g2*(vu**3*sqrt2-vu*sqrt2*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2)
  mH2(2)=t1*(t2+t3)
  mL2(1)=(-((vR*vu*Ahn(1))/sqrt2)-vu*hn(1)*((h*vp&
    &*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2) &
    *vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+(gp2*vL(1)*(-vd**2/2._dp&
    &+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp-vL(3)**2/2._dp))/4._dp-vR&
    &*hn(1)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2))/sqrt2)-(g2*(2*sqrt2*vL(1)*(vd**2/2._dp+vL(1)**2 &
    /2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)-sqrt2*vL(1)*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vL(1)
  mL2(2)=(-((vR*vu*Ahn(2))/sqrt2)-vu*hn(2)*((h*vp&
    &*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2) &
    *vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+(gp2*vL(2)*(-vd**2/2._dp&
    &+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp-vL(3)**2/2._dp))/4._dp-vR&
    &*hn(2)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2))/sqrt2)-(g2*(2*sqrt2*vL(2)*(vd**2/2._dp+vL(1)**2 &
    /2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)-sqrt2*vL(2)*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vL(2)
  mL2(3)=(-((vR*vu*Ahn(3))/sqrt2)-vu*hn(3)*((h*vp&
    &*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2) &
    *vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+(gp2*vL(3)*(-vd**2/2._dp&
    &+vu**2/2._dp-vL(1)**2/2._dp-vL(2)**2/2._dp-vL(3)**2/2._dp))/4._dp-vR&
    &*hn(3)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2))/sqrt2)-(g2*(2*sqrt2*vL(3)*(vd**2/2._dp+vL(1)**2 &
    /2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)-sqrt2*vL(3)*(vd**2/2._dp&
    &+vu**2/2._dp+vL(1)**2/2._dp+vL(2)**2/2._dp+vL(3)**2/2._dp)))/(4._dp*sqrt2))/vL(3)

  mP2=(-(BMphi*vp)-h*vR*((h*vp*vR)/2._dp+(MR*vR) &
    /sqrt2)+rh0*vu*(-(rh0*vp*vu)/2._dp-( mu *vu)/sqrt2)-(Almbd&
    &*vp**2)/(2._dp*sqrt2)-(Ah*vR*vS)/sqrt2+(RAh0*vd*vu)/sqrt2&
    &-(Mphi+(lmbd*vp)/sqrt2)*( (lmbd*vp**2) &
    /4._dp+(h*vR*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)*sqrt2&
    &-h*vS*((h*vp*vS)/2._dp+(MR*vS)/sqrt2+(vu*hn(1)*vL(1))/2._dp&
    &+(vu*hn(2)*vL(2))/2._dp+(vu*hn(3)*vL(3))/2._dp)+rh0*vd*(-(rh0&
    &*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1)*vL(1))/sqrt2&
    &+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3))/sqrt2))/sqrt2))/vp
  mR2=(-(BMR*vR)-h*vR*( (lmbd*vp**2)/4._dp+(h*vR&
    &*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2)-(Ah*vp*vR)/sqrt2&
    &-(MR+(h*vp)/sqrt2)*sqrt2*((h*vp*vS)/2._dp+(MR*vS)/sqrt2&
    &+(vu*hn(1)*vL(1))/2._dp+(vu*hn(2)*vL(2))/2._dp+(vu*hn(3)&
    & *vL(3))/2._dp))/vS
  mS2=(-(BMR*vS)-(vR*vu**2*hn(1)**2)/2._dp-(vR*vu**2&
    &*hn(2)**2)/2._dp-(vR*vu**2*hn(3)**2)/2._dp-h*vS*( (lmbd&
    &*vp**2)/4._dp+(h*vR*vS)/2._dp-(rh0*vd*vu)/2._dp+(Mphi*vp)/sqrt2) &
    -(Ah*vp*vS)/sqrt2-(MR+(h*vp)/sqrt2)*((h*vp*vR)/2._dp+(MR*vR) &
    /sqrt2)*sqrt2-(vu*Ahn(1)*vL(1))/sqrt2-(vu*Ahn(2) &
    *vL(2))/sqrt2-(vu*Ahn(3)*vL(3))/sqrt2-sqrt2*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) &
    /sqrt2)*(-(rh0*vd*vp)/2._dp-( mu *vd)/sqrt2+(vR*((hn(1) &
    *vL(1))/sqrt2+(hn(2)*vL(2))/sqrt2+(hn(3)*vL(3)) /sqrt2))/sqrt2))/vR

  Iname = Iname - 1

 End Subroutine ParMin_sp


End Module RPtools
