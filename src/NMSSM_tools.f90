Module NMSSM_tools

!---------------------------------
! loading necessary modules
!---------------------------------
Use BranchingRatios
Use Control
Use EplusEminusProduction
Use LoopCouplings
Use LoopFunctions 
Use LoopMasses
Use Model_Data
Use RGEs
Use StandardModel
Use SugraRuns

Contains


 Subroutine Model_NMSSM(m32, Grav_fac, F_GMSB, Ecms, Pm, Pp, ISR, Beam      &
   & , SigSup , SigSdown, SigSle, SigSn, SigC, SigChi0, SigS0, SigSP, SigHp &
   & , kont)
 Implicit None

  !---------------------------------------
  ! input / output
  !---------------------------------------
  Integer, Intent (inout) :: kont
  !-------------------------------------------------------------------
  ! are needed in case of light gravitinos, e.g. in GMSB like models
  !-------------------------------------------------------------------
  Real(dp), Intent(in) :: m32, Grav_fac, F_GMSB
  !--------------------------------
  ! cross section calculation
  !--------------------------------
  Real(dp), Intent(in) :: Ecms(:), Pm(:), Pp(:)
  Logical, Intent(in) :: ISR(:), Beam(:)
  Real(dp), Intent(out) :: SigSup(:,:,:) , SigSdown(:,:,:), SigSle(:,:,:) &
         & , SigSn(:,:,:), SigC(:,:,:), SigChi0(:,:,:), SigS0(:,:)        &
         & , SigSP(:,:,:), SigHp(:,:,:)

  !---------------------------------------
  ! local variables
  !---------------------------------------
    Real(dp) :: gp, g, vev, vev2, vevs_DR(2), sinW2, mP0_T(3), mP02_T(3) &
           &, RP0_T(3,3), mS0_T(3), mS02_T(3), RS0_T(3,3), mSpm_T(2)   &
           &, mSpm2_T(2), tad(3), mZ2_run, mW2_run, US(3,3), UP(3,3)   &
           &, cosb2, sinb2, matMSSM(2,2)
  Real(dp) :: gauge_mZ(3), ht, hb, htau, lam_mZ, kappa_mZ, Scale_Q, tz, dt &
           &, g8(8), g9(9), g6(6), mf_l_save(3), mf_d_save(3), mf_u_save(3)
  Real(dp) :: delta = 1.e-5_dp, logQ, alphas_DR, p2_in(3), pp2_in(3)
  Complex(dp) :: mu_in, dmZ2, RSpm_T(2,2)
  Integer :: i1 
  Real(dp), Parameter :: &
    & as2loop = 1._dp / 24._dp + 2011._dp * oo32Pi2 / 12._dp          &
    &         + Log2 / 12._dp - oo8Pi2 * Zeta3                        &
    & , log2loop_a = 123._dp * oo32Pi2, log2loop_b = 33._dp * oo32Pi2
 !-------------------------------
 ! widths and branching ratios
 !-------------------------------
  Real(dp) ::  epsI, deltaM, ratioWoM
  Logical :: CalcTBD
  !--------------------------------
  ! cross section calculation
  !--------------------------------
  Integer :: p_max
  Complex(dp) :: Ylp(3,3)

  Iname = Iname + 1 
  NameOfUnit(Iname) = "Model_NMSSM"

  kont = 0

  sinW2 = 1._dp - mW2/mZ2 
  
  alpha_mZ = Alpha_MSbar(mZ, mW)
  gauge_mZ(1) = Sqrt(4._dp * pi * alpha_mZ/(1._dp - sinW2))
  gauge_mZ(2) = Sqrt(4._dp * pi * alpha_mZ/sinW2) 
  alphas_DR =  AlphaS_mZ / (1._dp - oo4pi *AlphaS_mZ )
  gauge_mZ(3) = Sqrt(4._dp * pi * alphas_DR)

  tanb_mZ = tanb

  vev = Sqrt(mZ2 * (1._dp - sinW2) *sinW2 / (pi * alpha_mZ)) 
  vevSM(1) = vev / Sqrt(1._dp + tanb**2)
  vevSM(2) = tanb * vevSM(1)
! wird in InputOutput definiert zu vP = sqrt2 * mu_eff/lambda

  Y_l = 0._dp
  Y_d = 0._dp
  Y_u = 0._dp
  A_l = 0._dp
  A_d = 0._dp 
  A_u = 0._dp 

  Scale_Q = Sqrt(GetRenormalizationScale())

  htau = Sqrt2 * mf_l_mZ(3) / vevSM(1)
  hb = Sqrt2 * mf_d_mZ(3) / vevSM(1)
  ht = Sqrt2 * mf_u_mZ(3) / vevSM(2)
  logQ = 2._dp * Log(mZ/mf_u(3))
  ht = ht * (1._dp - alphas_DR / (3._dp*pi) * (5._dp +3._dp * LogQ  &
       &   + (as2loop + log2loop_a * logQ                      &
       &                         + log2loop_b * logQ**2) *4*pi*alphas_DR))

  !-------------------------------------------------------------
  ! running up within MSSM to get first approximation at Q_EWSB
  !-------------------------------------------------------------
   g6(1:3) = gauge_mZ
   g6(1) = Sqrt(5._dp/3._dp) * g6(1)
   g6(4) = htau
   g6(5) = hb
   g6(6) = ht

   tz = Log(Scale_Q/mZ)
   dt = tz / 50._dp
   Call odeint(g6, 6, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge6, kont)
   !-----------------------------------------------------------------------
   ! adding NMSSM couplings, run down to m_Z to get the NMSSM couplings at 
   ! m_Z, repeating this a second time including tan(beta)
   !-----------------------------------------------------------------------
   g8(1:6) = g6
   g8(7) = h0
   g8(8) = 0.5_dp * lam ! using convention of hep-ph/9505326

   Call odeint(g8, 8, tz, 0._dp, 0.1_dp*delta, - dt, 0._dp, rge8_NMSSM, kont)

   g8(1:3) = gauge_mZ
   g8(1) = Sqrt(5._dp/3._dp) * g8(1) 
   g8(4) = htau
   g8(5) = hb
   g8(6) = ht
   Call odeint(g8, 8, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge8_NMSSM, kont)

   g8(7) = h0
   g8(8) = 0.5_dp * lam ! using convention of hep-ph/9505326

   Call odeint(g8, 8, tz, 0._dp, 0.1_dp*delta, - dt, 0._dp, rge8_NMSSM, kont)

   g9(1:3) = gauge_mZ
   g9(1) = Sqrt(5._dp/3._dp) * g9(1)
   g9(4) = htau
   g9(5) = hb
   g9(6) = ht
   g9(7:8) = g8(7:8) 
   g9(9) = Log(tanb)
   Call odeint(g9, 9, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge9_NMSSM, kont)

   gp = Sqrt(3._dp/5._dp) * g9(1)
   g = g9(2)
   gauge(2:3) = g9(2:3)
   gauge(1) = gp
   Y_l(3,3) = g9(4)
   Y_d(3,3) = g9(5)
   Y_u(3,3) = g9(6)
 
   A_u(3,3) = Y_u(3,3) * AoY_u(3,3)
   A_d(3,3) = Y_d(3,3) * AoY_d(3,3)
   A_l(3,3) = Y_l(3,3) * AoY_l(3,3)

   tanb_Q = Exp(g9(9))

   vevSM(1) = vev / Sqrt(1._dp + tanb_Q**2)
   vevSM(2) = tanb_Q * vevSM(1)
 
   mu_in = 0._dp

   If (GenerationMixing) Then
    i1 = GetYukawaScheme()
    If (i1.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM),Y_d) 
    End If
   End If

 ! calculates all SusyMasses in the NMSSM at tree level except for
 ! neutral Higgs bosons where the 1-loop effective potential is used
 ! false = Higgs bosons treelevel 
 ! true = Higgs bosons 1loop effective potential

 If (.Not.External_Spectrum) &
  & Call TreeMassesNMSSM(gp, g, vevSM, vP, Mi(1), Mi(2), Mi(3), mu_in, B, h0  &
            &, lam, A_h0, A_lam, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q    &
            &, A_d, A_u, Y_d, Y_u, mGlu, PhaseGlu, mC, mC2, U, V, mN5      &
            &, mN52, N5, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2      &
            &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2, RSup        &
            &, mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T, mSpm_T, mSpm2_T &
            &, RSpm_T, GenerationMixing, kont, .True.)

  If (kont.Ne.0) Then
   Iname = Iname - 1
   Return 
  End If

    If (GenerationMixing) Then
     Call FermionMass(Y_l,vevSM(1),mf_l,uL_L,uL_R,kont)
     Call FermionMass(Y_d,vevSM(1),mf_d,uD_L,uD_R,kont)
     Call FermionMass(Y_u,vevSM(2),mf_u,uU_L,uU_R,kont)
    Else
     uL_L = id3C
     uL_R = id3C
     uD_L = id3C
     uD_R = id3C
     uU_L = id3C
     uU_R = id3C
    End If

  !----------------------------------------
  ! decays of SUSY and Higgs particles
  !----------------------------------------
  If ((L_BR).And.(kont.Eq.0)) Then
   ! relative precision for the calculation of phase space integrals
   epsI = 1.e-5_dp
   deltaM = 1.e-3_dp ! if mass/phase space is smaller, than mass is treated as 0
   CalcTBD = .False. ! if .True. than calculation of 3-body decays is enforced
   ratioWoM = 0._dp ! 1.e-4_dp

   Call CalculateBR(gauge, mGlu, PhaseGlu, mC, U, V, mN5, N5, mSneut, RSneut &
     & , mSlepton, RSlepton, mSup, RSup, mSdown, RSdown, uL_L, uL_R          &
     & , uD_L, uD_R, uU_L, uU_R, mS03, RS03, mP03, RP03, mSpm, RSpm, epsI    &
     & , deltaM, CalcTBD, kont, ratioWoM, Y_d, A_d, Y_l, A_l, Y_u, A_u       &
     & , mu, h0, A_h0, lam, A_lam, vevSM, vP, F_Gmsb, m32, grav_fac          &
     & , gP_Sl, gT_Sl, BR_Sl, gP_Sn, gT_Sn, BR_Sn, gP_Sd, gT_Sd, BR_Sd       &
     & , gP_Su, gT_Su, BR_Su, gP_C, gT_C, BR_C, gP_N5, gT_N5, BR_N5          &
     & , gP_Glu, gT_Glu, BR_Glu, gP_P03, gT_P03, BR_P03                      &
     & , gP_S03, gT_S03, BR_S03, gP_Spm, gT_Spm, BR_Spm)

  End If

  !-------------------------------------------------------------
  ! cross section for e+ e- annihilation
  ! The following input quantities can be specified either in the file 
  ! CrossSections.in or in the file LesHouches.in:
  !                    Ecms .... c.m.s enerergy in GeV
  !                   Pm ...... degree of longitudinal polarization of incoming
  !                             electron
  !                   Pp ...... degree of longitudinal polarization of incoming
  !                             positron
  !                   ISR ..... if .TRUE. then the effect of initial state
  !                             radiation will be included
  ! In the case that none of these files exist, the following
  ! default values are used: Ecms = 500 GeV, Pm = Pp = 0, ISR = .TRUE.
  !-------------------------------------------------------------
  If ((L_CS).And.(kont.Eq.0)) Then
    Ylp = Y_l / gauge(2)
    p_max = Size(Pm)

    Do i1=1,p_max
     If (Ecms(i1).Eq.0._dp) Exit
     Call CalculateCrossSections(Ecms(i1), Pm(i1), Pp(i1), ISR(i1), Beam(i1)   &
            & , "Tesla800", mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu       &
            & , SigSup(i1,:,:), SigSdown(i1,:,:), mSlepton, RSlepton, Ylp      &
            & , mSneut, RSneut, SigSle(i1,:,:), SigSn(i1,:,:), mC, U, V        &
            & , mN5, N5, SigC(i1,:,:), SigChi0(i1,1:5,1:5), mS03, RS03, vevSM  &
            & , mP03, RP03, mSpm, RSpm, SigS0(i1,1:3), SigSP(i1,1:3,1:2)       &
            & , SigHp(i1,1,1) )
    End Do
   End If

  Iname = Iname - 1
 
 End Subroutine Model_NMSSM

End Module NMSSM_tools
