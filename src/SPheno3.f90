!------------------------------------------------------------------------
! This is the main program for the package SPheno, version 3.0beta
! SPheno stands for S(upersymmetric) Pheno(menology) 
! A documentation is provided in the file Spheno.ps
! Please send any comments and/or bugs to porod@physik.uni-wuerzburg.de
! In case of a bug:
!   - please check that the integer in the first line of Control.in
!     is 0 or larger. If it is negative, please set it to 0. This
!     enables the writing of several warnings. Afterwards rerun the 
!     program.
!   - please send the file Messages.out  and all input files (*.in) to
!     the email above, if possible with a copy of the screen output 
! written by Werner Porod
!-----------------------------------------------------------------
Program SPheno

!---------------------------------
! loading necessary modules
!---------------------------------
Use BranchingRatios
Use Chargino3Decays
Use Control
Use EplusEminusProduction
Use InputOutput
Use LoopFunctions
Use LoopMasses
Use LowEnergy
Use Mathematics
Use Model_Data
Use NMSSM_tools 
Use RPtools
Use RGEs
Use StandardModel
Use SugraRuns

! Use f90_unix_proc
 Implicit None
 !-------------------------------
 ! widths and branching ratios
 !-------------------------------
  Real(dp) ::  epsI, deltaM
 !--------------------------------
 ! cross section calculation
 !--------------------------------
 Integer, Parameter :: p_max=100
 Complex(dp) :: Ylp(3,3)
 Real(dp) :: Ecms(p_max), Pm(p_max), Pp(p_max), SigSup(p_max,6,6)        &
         & , SigSdown(p_max,6,6), SigSle(p_max,6,6), SigSn(p_max,3,3)    &
         & , SigC(p_max,5,5), SigChi0(p_max,7,7), SigS0(p_max,5)         &
         & , SigSP(p_max,5,4), SigHp(p_max,7,7)
 Logical :: ISR(p_max)=.False., Beam(p_max)=.False.
 !----------------------------------
 ! low energy constraints
 !----------------------------------
 Real(dp) :: BRbtosgamma, Bs_mumu, BrBToSLL, BR_Bu_TauNu, a_e, a_mu, a_tau &
   & , d_e, d_mu, d_tau, BrMutoEGamma, BrTautoEGamma, BrTautoMuGamma       &
   & , BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau, BtoSNuNu
 Complex(dp) :: DeltaMBd, DeltaMBs
 !------------------------------------
 ! auxiliary parameters
 !------------------------------------
 Logical :: CalcTBD
 Integer :: kont, i1, i_min(3)
 Real(dp) :: m_Gut, ratioWoM

 !--------------------------------------------------------------------------
 ! set all parameters to zero
 !--------------------------------------------------------------------------
 Call Set_All_Parameters_0()

 !--------------------------------------------------------------------------
 ! This routine call routines to
 !   - initializ the error system
 !   - to calculate the constants which are required for the 
 !     calculation of the loop functions
 !   - to get the standard model parameters
 !   - to get SUSY parameters
 ! The following steps are performed to get the parameters and flags
 !   (i) The file LesHouches.in exists containing all necessary information.
 !       In this case the remaining steps are skipped
 !   (ii) Check if Control.in exists to set the error level and
 !        to check if widths and cross sections shall be calculated
 !   (iii) Check if StandardModel.in exists to change the default SM values
 !   (iv) Read the information concerning the SUSY model from the file
 !        HighScale.in
 ! Note that either the file LesHouches.in or the file HighScale.in
 ! must exist.
 !--------------------------------------------------------------------------
 kont = 0
 Call ReadingData(kont)

 !---------------------------------------------
 ! parameters for branching ratio calculations
 !---------------------------------------------
  ! relative precision for the calculation of phase space integrals
  epsI = 1.e-5_dp
  deltaM = 1.e-3_dp ! if mass/phase space is smaller, than mass is treated as 0
  CalcTBD = .False. ! if .True. than calculation of 3-body decays is enforced
  ratioWoM = 0._dp ! 1.e-4_dp
 

 If (HighScaleModel.Eq."NMSSM") Then ! NMSSM model

  Call Model_NMSSM(m32, Grav_fac, F_GMSB, Ecms, Pm, Pp, ISR, Beam           &
   & , SigSup , SigSdown, SigSle, SigSn, SigC, SigChi0, SigS0, SigSP, SigHp &
   & , kont)

 Elseif ((HighScaleModel.Eq."RPexplicit").Or.(Add_Rparity)) Then ! bilinear RP

  Call Model_bilinear_Rparity(add_Rparity, HighScaleModel, delta_mass, epsI     &
       & , deltaM, ratioWoM, m32, grav_fac, CalcTBD, Ecms, Pm, Pp, ISR, Beam    &
       & , SigSup , SigSdown, SigC, SigChi0, SigS0, SigSP, SigHp, M_GUT, kont)
  
 Else If (kont.Eq.0) Then  ! models with MSSM particle content
                           ! at the electroweak scale 
  
  !---------------------------------------------------------------------------
  ! calculation of the spectrum, the following parameters can be changed
  ! with the help of the SLHA input file LesHouches.in or the file Control.in
  ! 
  !---------------------------------------------------------------------------
   Call CalculateSpectrum(n_run, delta_mass, WriteOut, kont, tanb &
    & , vevSM, mC, U, V, mN, N, mS0, mS02, RS0, mP0, mP02, RP0, mSpm, mSpm2    &
    & , RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup, mSlepton, mSlepton2  &
    & , RSlepton, mSneut, mSneut2, RSneut, mGlu, PhaseGlu, gauge               &
    & , uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, Y_l, Y_d, Y_u                      &
    & , Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B           &
    & , m_GUT)
 !-------------------------------------------------------------------
 ! Calculation of the branching ratios and widths provided L_BR is
 ! set .TRUE. (default) and that the routine Sugra has finished
 ! correctly (kont.eq.0) 
 !-------------------------------------------------------------------

  If ((L_BR).And.(kont.Eq.0)) Then

   Call CalculateBR(gauge, mGlu, PhaseGlu, mC, U, V, mN, N, mSneut, RSneut  &
     & , mSlepton, RSlepton, mSup, RSup, mSdown, RSdown, uL_L, uL_R         &
     & , uD_L, uD_R, uU_L, uU_R, mS0, RS0, mP0, RP0, mSpm, RSpm, epsI       &
     & , deltaM, CalcTBD, kont, ratioWoM, Y_d, A_d, Y_l, A_l, Y_u, A_u, mu  &
     & , vevSM, F_Gmsb, m32, grav_fac, gP_Sl, gT_Sl, BR_Sl                  &
     & , gP_Sn, gT_Sn, BR_Sn , gP_Sd, gT_Sd, BR_Sd, gP_Su, gT_Su, BR_Su     &
     & , gP_C, gT_C, BR_C, gT_N, gP_N4_2, BR_N4_2, gP_N4_3, BR_N4_3         &
     & , gP_Glu, gT_Glu, BR_Glu, gP_P0, gT_P0, BR_P0      &
     & , gP_S0, gT_S0, BR_S0, gP_Spm, gT_Spm, BR_Spm)

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
   Ylp = Y_l / gauge(2)
   Do i1=1,p_max
    If (Ecms(i1).Eq.0._dp) Exit
    Call CalculateCrossSections(Ecms(i1), Pm(i1), Pp(i1), ISR(i1), Beam(i1)  &
           & , "Tesla800", mSup, RSup, mf_u, mSdown, RSdown, mf_d, mglu      &
           & , SigSup(i1,:,:), SigSdown(i1,:,:), mSlepton, RSlepton, Ylp     &
           & , mSneut, RSneut, SigSle(i1,:,:), SigSn(i1,:,:), mC, U, V       &
           & , mN, N, SigC(i1,1:2,1:2), SigChi0(i1,1:4,1:4), mS0, RS0, vevSM &
           & , mP0, RP0, mSpm, RSpm, SigS0(i1,1:2), SigSP(i1,1:2,1)          &
           & , SigHp(i1,1,1) )
   End Do

  End If

 End If

 If ((kont.Eq.0).And.(HighScaleModel.Ne."NMSSM")) Then

  Call CalculateLowEnergyConstraints(gauge, Y_l, Y_d, Y_u      &
    & , mSpm2, RSpm, mC, U, V, mN, N , mSup2, RSup, mSdown2    &
    & , RSdown, mSlepton2, RSlepton, mSneut2, RSneut           &
    & , BRbtosgamma, DeltaMBd, DeltaMBs, BrBToSLL, BtoSNuNu, BR_Bu_TauNu &
    & , a_e, a_mu, a_tau, d_e, d_mu, d_tau                     &
    & , BrMutoEGamma, BrTautoEGamma, BrTautoMuGamma, BrMu3e, BrTau3e, BrTau3Mu &
    & , BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau)

 Else
  BRbtosgamma = 0._dp
  BToSNuNu = 0._dp
  BrBToSLL = 0._dp
  DeltaMBd = 0._dp
  DeltaMBs = 0._dp
  bs_mumu = 0._dp
  BR_Bu_TauNu = 0._dp
  a_e = 0._dp
  a_mu = 0._dp
  a_tau = 0._dp
  d_e = 0._dp
  d_mu = 0._dp
  d_tau = 0._dp
  BrMutoEGamma = 0._dp
  BrTautoEGamma = 0._dp
  BrTautoMuGamma = 0._dp
  BrMu3e = 0._dp
  BrTau3e = 0._dp
  BrTau3Mu = 0._dp
  BR_Z_e_mu = 0._dp
  BR_Z_e_tau = 0._dp
  BR_Z_mu_tau = 0._dp
  rho_parameter = 0._dp
 End If


! If (kont.Ne.0) Then
!  Write(*,*) "there has been a problem, kont = ",kont
! Else

  Call LesHouches_Out(67, 11, kont, HighScaleModel, M_GUT                   &
      & , BRbtosgamma, Bs_MuMu, DeltaMBd, DeltaMBs, BrBToSLL, BtoSNuNu      &
      & , BR_Bu_TauNu  &
      & , a_e, a_mu, a_tau, d_e, d_mu, d_tau, BrMuToEGamma, BrTauToEGamma  &
      & , BrTauToMuGamma, BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau &
      & , BR_Z_mu_tau                       &
      & , Rho_parameter, Ecms, Pm, Pp, ISR, SigSup, SigSdown, SigSle      &
      & , SigSn, SigChi0, SigC, SigS0, SigSP, SigHp)
  !------------------------------------------------------------
  ! programs like micrOmegas do not yet use flavour mixing, in
  ! this case a modified SLHA file is needed
  !------------------------------------------------------------
  If (Write_SLHA1) Call WriteSPhenoOutputLHA1(35, M_GUT)

! End If

 Call closing() ! closes the files
 If ((kont.Ne.0).And.Non_Zero_Exit) Stop 99

Contains


 Subroutine CalculateLowEnergyConstraints(gauge, Y_l, Y_d, Y_u, mSpm2, RSpm   &
   & , mC, U, V, mN, N , mSup2, RSup, mSdown2, RSdown, mSlepton2              &
   & , RSlepton, mSneutrino2, RSneutrino, BRbtosgamma, DeltaMBd, DeltaMBs     &
   & , BrBToSLL, BtoSNuNu, BR_Bu_TauNu, a_e, a_mu, a_tau, d_e, d_mu, d_tau    &
   & , BrMutoEGamma, BrTautoEGamma, BrTautoMuGamma, BrMu3e, BrTau3e, BrTau3Mu &
   & , BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau)

 Implicit None
  Real(dp), Intent(in) :: gauge(3), mSpm2(2), mC(2), mN(4), mSneutrino2(3) &
     & , mSup2(6), mSdown2(6), mSlepton2(6)
  Complex(dp), Intent(in), Dimension(3,3) :: Y_d, Y_u, Y_l, RSneutrino 
  Complex(dp), Intent(in) ::  RSpm(2,2), U(2,2), V(2,2), RSup(6,6)   &
     & , RSdown(6,6), RSlepton(6,6), N(4,4)
  Real(dp), Intent(out) :: BRbtosgamma, BrBToSLL, BR_Bu_TauNu, a_mu, a_e      &
   & , a_tau, d_e, d_mu, d_tau, BrMutoEGamma, BrTautoEGamma, BrTautoMuGamma   &
   & , BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau, BtoSNuNu
  Complex(dp), Intent(out) :: DeltaMBd, DeltaMBs
  Real(dp) :: MuE_conv_Ti

  Real(dp) :: mGlu_T, mC_T(2), mC2_T(2), mN_T(4), mN2_T(4), mSneutrino_T(3)   &
     & , mSneutrino2_T(3), mSlepton_T(6), mSlepton2_T(6), mSdown_T(6)         &
     & , mSdown2_T(6), mSup_T(6), mSup2_T(6), mP0_T(2), mP02_T(2), RP0_T(2,2) &
     & , mS0_T(2), mS02_T(2), RS0_T(2,2), mSpm_T(2), mSpm2_T(2),mZ2_run, mW2_run
  Complex(dp) :: Phase_Glu_T, U_T(2,2), V_T(2,2), N_T(4,4), Rsneut_T(3,3)  &
     & , RSlepton_T(6,6), RSdown_T(6,6), RSup_T(6,6), RSpm_T(2,2)
  Real(dp) :: mudim, dt, tz, g2(213), vev2, vevs_DR(2), sinW2_DR
  Complex(dp) :: CKMad(3,3), mu_mZ, B_mZ, Mi_mZ(3), dmZ2
  Complex(dp) :: cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), cpl_CLSn_L(2,3,3) &
    & , cpl_CLSn_R(2,3,3)
  Integer :: i1,i2,i3, scheme
 Real(dp) :: GMutoEGamma, GTautoEGamma, GTautoMuGamma
 Complex(dp) :: cpl_LNSl_L(3,4,6), cpl_LNSl_R(3,4,6), cpl_NuNSn_R(3,4,3)       &
    & , cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6), cpl_NNZ_L(4,4), cpl_NNZ_R(4,4) &
    & , cpl_DGSd_L(3,6), cpl_DGSd_R(3,6), cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6) &
    & , cpl_CCZ_L(2,2), cpl_CCZ_R(2,2), cpl_SnSnZ(3,3), cpl_SlSlZ(6,6)
  Real(dp) :: mf_t(3), BtoSEE, EDM_e(3), EDM_mu(3), EDM_tau(3), gU1, gSU2  &
     & , cosW, sinW2, KtoPiNuNu
  Complex(dp) :: c7(7), c7p(6), c8(7), c8p(6)

  Y_l_mZ = Transpose(Y_l) 
  Y_d_mZ = Transpose(Y_d) 
  Y_u_mZ = Transpose(Y_u) 
  A_l_mZ = Transpose(A_l) 
  A_d_mZ = Transpose(A_d) 
  A_u_mZ = Transpose(A_u) 
  M2_E_mZ = M2_E
  M2_L_mZ = M2_L
  M2_D_mZ = M2_D
  M2_U_mZ = M2_U
  M2_Q_mZ = M2_Q
  Mi_mZ = Mi
  mu_mZ = mu
  B_mZ = B
  M2_H_mZ = M2_H

  Call ParametersToG(gauge, y_l_mZ, y_d_mZ, y_u_mZ, Mi_mZ, A_l_mZ, A_d_mZ   &
          & , A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ, M2_H_mZ &
          & , mu_mZ, B_mZ, g2)

  mudim = GetRenormalizationScale()
  tz = 0.5_dp * Log(mZ**2/mudim)

  If (tz.Ne.0._dp) Then
   dt = tz / 100._dp
   g2(1) = Sqrt(5._dp / 3._dp ) * g2(1)

   Call odeint(g2, 213, 0._dp, tz, delta_mass, dt, 0._dp, rge213, kont)
   g2(1) = Sqrt(3._dp / 5._dp ) * g2(1)

   Call GToParameters(g2, gauge_mZ, Y_l_mZ, Y_d_mZ, Y_u_mZ, Mi_mZ, A_l_mZ &
                  & , A_d_mZ, A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ  &
                  & , M2_U_mZ, M2_H_mZ, mu_mZ, B_mZ)
   Else
    gauge_mZ = gauge
   End If

   Y_l_mZ = Transpose(Y_l_mZ) 
   Y_d_mZ = Transpose(Y_d_mZ) 
   Y_u_mZ = Transpose(Y_u_mZ) 
   A_l_mZ = Transpose(A_l_mZ) 
   A_d_mZ = Transpose(A_d_mZ) 
   A_u_mZ = Transpose(A_u_mZ) 

  !-------------------------------------
  ! calculate running masses at m_Z
  !-------------------------------------
  sinW2 = 1._dp - mW2 / mZ2
  vev2 =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
  vevSM(1) = vev2 / Sqrt(1._dp + tanb_mZ**2)
  vevSM(2) = tanb_mZ * vevSM(1)

  mZ2_run = (gauge_mZ(1)**2+gauge_mZ(2))**2*0.25*(vevSM(1)**2+vevSM(2)**2)
  mW2_run = gauge_mZ(2)**2*0.25*(vevSM(1)**2+vevSM(2)**2)

  Call TreeMassesMSSM2(gauge_mZ(1), gauge_mZ(2), vevSM, Mi_mZ(1), Mi_mZ(2)   &
     & , Mi_mZ(3), mu_mZ, B_mZ, tanb_mZ, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ       &
     & , M2_D_mZ, M2_U_mZ, M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ           &
     & , uU_L, uU_R ,uD_L, uD_R, uL_L, uL_R                                  &
     & , mGlu_T, Phase_Glu_T, mC_T, mC2_T, U_T, V_T, mN_T, mN2_T, N_T        &
     & , mSneutrino_T, mSneutrino2_T, Rsneut_T, mSlepton_T, mSlepton2_T      &
     & , RSlepton_T, mSdown_T, mSdown2_T, RSdown_T, mSup_T, mSup2_T, RSup_T  &
     & , mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T, mSpm_T, mSpm2_T, RSpm_T &
     & , mZ2_run, mW2_run, GenerationMixing, kont, .False., .False.)

  if (kont.ne.0) then ! there is a problem with the running masses, use pole masses instead
   mGlu_T = mglu
   Phase_Glu_T = PhaseGlu
   mC_T = mC
   mC2_T = mC2
   U_T = U
   V_T = V
   mN_T = mN
   mN2_T = mN2
   N_T = N
   mSneutrino2_T = mSneutrino2
   Rsneut_T = Rsneut
   mSlepton2_T = mSlepton2
   RSlepton_T = RSlepton
   mSdown2_T = mSdown2
   RSdown_T = Rsdown
   mSup2_T = mSup2 
   RSup_T = RSup
   mP0_T = mP0
   mP02_T = mP02 
   RP0_T = RP0
   mS0_T = mS0 
   mS02_T = mS02 
   RS0_T = RS0 
   mSpm_T = mSpm
   mSpm2_T = mSpm2 
   RSpm_T = RSpm
  End If 

  If (.Not.GenerationMixing) Then ! need to add quark mixing for the following
   If (scheme.Eq.1) Then
    uU_L = CKM
    Y_u_mZ = Matmul(Transpose(CKM),Y_u_mZ)
   Else
    Call Adjungate(CKM, CKMad)
    uD_L = CKMad
    Y_d_mZ = Matmul(CKM,Y_d_mZ)
   End If
  End If

  !---------------------------------------
  ! BR(b-> s gamma)
  !---------------------------------------
!  Call BtoQGamma(2, mf_d_mZ, gauge_mZ, mf_u, mW, Y_d_mZ, uD_L, uD_R, Y_u_mZ  &
!   & , uU_L, uU_R, mSpm2, RSpm, mC, U, V, mSup2, RSup, A_u_mZ  &
!   & , mSdown2, RSdown, A_d_mZ, mglu, phaseGlu, mN, N, mu_mZ    &
!   & , mS02, RS0, mP02, RP0, vevSm, BRbtosgamma, c7, c7p, c8, c8p)
  Call BtoQGamma(2, mf_d_mZ, gauge_mZ, mf_u, mW, Y_d_mZ, uD_L, uD_R, Y_u_mZ  &
   & , uU_L, uU_R, mSpm2_T, RSpm_T, mC_T, U_T, V_T, mSup2_T, RSup_T, A_u_mZ  &
   & , mSdown2_T, RSdown_T, A_d_mZ, mglu_T, phase_Glu_T, mN_T, N_T, mu_mZ    &
   & , mS02_T, RS0_T, mP02_T, RP0_T, vevSm, BRbtosgamma, c7, c7p, c8, c8p)

!  Call Delta_MB(1, mf_u, gauge_mZ, Y_u_mZ, uU_L, uU_R, Y_d_mZ, uD_L, uD_R     &
!      & , mC, U, V, mN, N, mGlu, phaseGlu, mS02, RS0, mP02, RP0, mSpm2, RSpm  &
!      & , mSup2, RSup, A_u_mZ, mu_mZ, mSdown2, RSdown, A_d_mZ, vevSM, DeltaMBd)
!  Call Delta_MB(2, mf_u, gauge_mZ, Y_u_mZ, uU_L, uU_R, Y_d_mZ, uD_L, uD_R     &
!      & , mC, U, V, mN, N, mGlu, phaseGlu, mS02, RS0, mP02, RP0, mSpm2, RSpm  &
!      & , mSup2, RSup, A_u_mZ, mu_mZ, mSdown2, RSdown, A_d_mZ, vevSM, DeltaMBs)
  Call Delta_MB(1, mf_u, gauge_mZ, Y_u_mZ, uU_L, uU_R, Y_d_mZ, uD_L, uD_R  &
      & , mC_T, U_T, V_T, mN_T, N_T, mGlu_T, phase_Glu_T, mS02_T, RS0_T    &
      & , mP02_T, RP0_T, mSpm2_T, RSpm_T, mSup2_T, RSup_T, A_u_mZ, mu_mZ   &
      & , mSdown2_T, RSdown_T, A_d_mZ, vevSM, DeltaMBd)
  Call Delta_MB(2, mf_u, gauge_mZ, Y_u_mZ, uU_L, uU_R, Y_d_mZ, uD_L, uD_R  &
      & , mC_T, U_T, V_T, mN_T, N_T, mGlu_T, phase_Glu_T, mS02_T, RS0_T    &
      & , mP02_T, RP0_T, mSpm2_T, RSpm_T, mSup2_T, RSup_T, A_u_mZ, mu_mZ   &
      & , mSdown2_T, RSdown_T, A_d_mZ, vevSM, DeltaMBs)
  ! conversion to pico-seconds
  DeltaMBd = 1.e-12_dp*DeltaMBd/hbar
  DeltaMBs = 1.e-12_dp*DeltaMBs/hbar

  Do i1=1,2
   Do i2=1,3
    Do i3=1,3     
     Call CoupCharginoSfermion(i1, i2, i3, gauge_mZ(2), -0.5_dp, RSneut_T  &
             & , Y_l_mZ, Zero33C, uL_L, uL_R, U_T, V_T, cpl_CLSn_L(i1,i2,i3) &
             & , cpl_CLSn_R(i1,i2,i3) )
    End Do
    Do i3=1,6
     Call CoupCharginoSfermion(i1, i2, i3, gauge_mZ(2), 0.5_dp, RSlepton_T &
      & , Y_l_mZ, Zero33C, uL_L, uL_R, U_T, V_T, cpl_CNuSl_L(i1,i2,i3)         &
      & , cpl_CNuSl_R(i1,i2,i3) )
     Call CoupCharginoSfermion(i1, i2, i3, gauge_mZ(2), -0.5_dp, RSup_T   &
           & , Y_D_mZ, Y_U_mZ, id3C, id3C, U_T, V_T, cpl_CDSu_L(i1,i2,i3)      &
           & , cpl_CDSu_R(i1,i2,i3))
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i2=1,4
    Do i3=1,6  
     Call CoupNeutralinoSlepton(i1, i2, i3, gauge_mZ(1), gauge_mZ(2), RSlepton_T &
         & , uL_L, uL_R, Y_l_mZ, N_T, cpl_LNSl_L(i1,i2,i3), cpl_LNSl_R(i1,i2,i3) )
     Call CoupNeutralinoSdown(i1, i2, i3, gauge_mZ(1), gauge_mZ(2), RSdown_T &
         & , uD_L, uD_R, Y_d_mZ, N_T, cpl_DNSd_L(i1,i2,i3), cpl_DNSd_R(i1,i2,i3) )
    End Do
    Do i3=1,3  
     Call CoupNeutralinoSneutrino(i1, i2, i3, gauge_mZ(1), gauge_mZ(2), N &
           & , RSneut_T, uL_R, cpl_nuNSn_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  Do i1=1,3
    Do i3=1,6  
     Call CoupGluinoSquark3(gauge_mZ(3), phase_Glu_T, i1, i3, RSdown_T, uD_L, uD_R &
           & , cpl_DGSd_L(i1,i3), cpl_DGSd_R(i1,i3) )
    End Do
  End Do

  !-------------------------
  ! B_s -> mu+ mu-
  !-------------------------
!  Call Bs_to_MuMu(2, mf_u, gauge_mZ, Y_u_mZ, uU_L, uU_R, Y_d_mZ, uD_L, uD_R    &
!      & , mC, U, V, mN, N, mGlu, phaseGlu, mS02, RS0, mP02, RP0, mSpm2, RSpm   &
!      & , mSup2, RSup, A_u, mu, mSdown2, RSdown, A_d, vevSM, mSneutrino2       &
!      & , mSlepton2, cpl_CDSu_L, cpl_CDSu_R, cpl_CLSn_L, cpl_CLSn_R, cpl_LNSl_L&
!      & , cpl_LNSl_R, cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R, Bs_mumu )
  Call Bs_to_MuMu(2, mf_u, gauge_mZ, Y_u_mZ, uU_L, uU_R, Y_d_mZ, uD_L, uD_R   &
       & , mC_T, U_T, V_T, mN_T, N_T, mGlu_T, phase_Glu_T, mS02_T, RS0_T      &
       & , mP02_T, RP0_T, mSpm2_T, RSpm_T, mSup2_T, RSup_T, A_u_mZ, mu_mZ     &
       & , mSdown2_T, RSdown_T, A_d, vevSM, mSneutrino2_T, mSlepton2_T        &
       & , cpl_CDSu_L, cpl_CDSu_R, cpl_CLSn_L, cpl_CLSn_R, cpl_LNSl_L         &
       & , cpl_LNSl_R, cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R, Bs_mumu )

  !-------------------------
  ! b -> s l+ l-
  !-------------------------
!  Call BToSLL(gauge_mZ, mf_d_mZ, mf_u_mZ, mW, Y_d_mZ, uD_L, uD_R, Y_u_mZ       &
!       & , uU_L, uU_R, Y_l_mZ, uL_L, UL_R, mSneutrino2, Rsneut, mSlepton2, Rslepton&
!       & , mSpm2, RSpm, mC, U, V, mSup2, RSup, A_u_mZ, mSdown2, RSdown, A_d_mZ &
!       & , mglu, phaseglu, mN, N, mu_mZ, mS02, RS0, mP02, RP0, vevSM           &
!       & , BtoSEE, BrBtoSLL)
  Call BToSLL(gauge_mZ, mf_d_mZ, mf_u_mZ, mW, Y_d_mZ, uD_L, uD_R, Y_u_mZ      &
        & , uU_L, uU_R, Y_l_mZ, uL_L, UL_R, mSneutrino2_T, Rsneut_T           &
        & , mSlepton2_T, Rslepton_T, mSpm2_T, RSpm_T, mC_T, U_T, V_T, mSup2_T &
        & , RSup_T, A_u_mZ, mSdown2_T, RSdown_T, A_d_mZ, mglu_T, phase_glu_T  &
        & , mN_T, N_T, mu_mZ, mS02_T, RS0_T, mP02_T, RP0_T, vevSM             &
        & , BtoSEE, BrBtoSLL)
  !---------------------------------
  ! b -> s nu nu, no QCD corrections
  !---------------------------------
!   Call B_To_SNuNu(gauge_mZ, mf_d_mZ, mf_u_mZ, mW, Y_d_mZ, uD_L, uD_R, Y_u_mZ &
!        & , uU_L, uU_R, Y_l_mZ, mSneut2, Rsneut, mSlepton2, Rslepton, mSpm2   &
!        & , RSpm, mC, U, V, mSup2, RSup, mSdown2, RSdown, mglu, phaseGlu      &
!        & , mN, N, vevSM, .False., BtoSNuNu)
   Call B_To_SNuNu(gauge_mZ, mf_d_mZ, mf_u_mZ, mW, Y_d_mZ, uD_L, uD_R, Y_u_mZ  &
     & , uU_L, uU_R, Y_l_mZ, mSneutrino2_T, Rsneut_T, mSlepton2_T, Rslepton_T  &
     & , mSpm2_T, RSpm_T, mC_T, U_T, V_T, mSup2_T, RSup_T, mSdown2_T, RSdown_T &
     & , mglu_T, phase_Glu_T , mN_T, N_T, vevSM, .False., BtoSNuNu)
  !-------------------
  ! B^-_u -> tau nu
  !-------------------
!  BR_Bu_TauNu = Bm_to_l_nu(3,1, mSpm2(2), tanb, RSpm, Y_d_mZ, uU_L &
!              &           , uD_R , Y_l_mZ, vevSM)
  BR_Bu_TauNu = Bm_to_l_nu(3,1, mSpm2_T(2), tanb_mZ, RSpm_T, Y_d_mZ, uU_L &
              &           , uD_R , Y_l_mZ, vevSM)

  !------------------------
  ! K -> pi nu nu
  !------------------------
   Call K_To_PiNuNu(gauge_mZ, mf_d_mZ, mf_u_mZ, mW, mZ, Y_d_mZ, uD_L, uD_R, Y_u_mZ &
   & , uU_L, uU_R, Y_l_mZ, mSneutrino2_T, Rsneut_T, mSlepton2_T, Rslepton_T        &
   & , mSpm2_T, RSpm_T, mC_T, U_T, V_T, mSup2_T, RSup_T, mSdown2_T, RSdown_T       &
   & , mglu_T, phase_Glu_T, mN_T, N_T, vevSM, .False., KtoPiNuNu)

  !------------------------------------------------------------------
  ! leptonic electric dipole moments
  !------------------------------------------------------------------
!  Call Lepton_EDM3(1, mN, mSlepton2, cpl_LNSl_L, cpl_LNSl_R      &
!                & , mC, mSneut2, cpl_CLSn_L, cpl_CLSn_R, EDM_e   )
!  Call Lepton_EDM3(2, mN, mSlepton2, cpl_LNSl_L, cpl_LNSl_R      &
!                & , mC, mSneut2, cpl_CLSn_L, cpl_CLSn_R, EDM_mu  )
!  Call Lepton_EDM3(3, mN, mSlepton2, cpl_LNSl_L, cpl_LNSl_R      &
!                 & , mC, mSneut2, cpl_CLSn_L, cpl_CLSn_R, EDM_tau )
  Call Lepton_EDM3(1, mN_T, mSlepton2_T, cpl_LNSl_L, cpl_LNSl_R        &
                & , mC_T, mSneutrino2_T, cpl_CLSn_L, cpl_CLSn_R, EDM_e   )
  Call Lepton_EDM3(2, mN_T, mSlepton2_T, cpl_LNSl_L, cpl_LNSl_R        &
                & , mC_T, mSneutrino2_T, cpl_CLSn_L, cpl_CLSn_R, EDM_mu  )
  Call Lepton_EDM3(3, mN_T, mSlepton2_T, cpl_LNSl_L, cpl_LNSl_R        &
                & , mC_T, mSneutrino2_T, cpl_CLSn_L, cpl_CLSn_R, EDM_tau )
  d_e = EDM_e(3)
  d_mu = EDM_mu(3)
  d_tau = EDM_tau(3)
  !------------------------------------------------------------------
  ! leptonic anomalous magnetic moments
  !------------------------------------------------------------------
!  Call Gminus2(1, mSneut2, mC, cpl_CLSn_L, cpl_CLSn_R, mSlepton2, mN &
!             &, cpl_LNSl_L, cpl_LNSl_R, a_e, GenerationMixing)
!  Call Gminus2(2, mSneut2, mC, cpl_CLSn_L, cpl_CLSn_R, mSlepton2, mN &
!             &, cpl_LNSl_L, cpl_LNSl_R, a_mu, GenerationMixing)
!  Call Gminus2(3, mSneut2, mC, cpl_CLSn_L, cpl_CLSn_R, mSlepton2, mN &
!             &, cpl_LNSl_L, cpl_LNSl_R, a_tau, GenerationMixing)
  Call Gminus2(1, mSneutrino2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlepton2_T &
             & , mN_T, cpl_LNSl_L, cpl_LNSl_R, a_e, GenerationMixing)
  Call Gminus2(2, mSneutrino2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlepton2_T &
             & , mN_T, cpl_LNSl_L, cpl_LNSl_R, a_mu, GenerationMixing)
  Call Gminus2(3, mSneutrino2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlepton2_T &
             & , mN_T, cpl_LNSl_L, cpl_LNSl_R, a_tau, GenerationMixing)
  !------------------------------------------------------------------
  ! rare decays of leptons: l -> l' gamma
  !------------------------------------------------------------------
  BrMutoEGamma = 0._dp
  BrTautoEGamma = 0._dp
  BrTautoMuGamma = 0._dp
  If (GenerationMixing) Then
!!$  Do i1=1,2
!!$   Do i2=1,3
!!$    Do i3=1,3     
!!$     Call CoupCharginoSfermion(i1, i2, i3, gauge_mZ(2), -0.5_dp, RSneut  &
!!$             & , Y_l_mZ, Zero33C, id3c, id3c, U, V, cpl_CLSn_L(i1,i2,i3) &
!!$             & , cpl_CLSn_R(i1,i2,i3) )
!!$    End Do
!!$   End Do
!!$  End Do
!!$  Do i1=1,3
!!$   Do i2=1,4
!!$    Do i3=1,6  
!!$     Call CoupNeutralinoSlepton(i1, i2, i3, gauge_mZ(1), gauge_mZ(2), RSlepton &
!!$         & , id3C, id3C, Y_l_mZ, N, cpl_LNSl_L(i1,i2,i3), cpl_LNSl_R(i1,i2,i3) )
!!$    End Do
!!$   End Do
!!$  End Do
!!$   Call LtoLpGamma(2, 1, mSneutrino2, mC, cpl_CLSn_L, cpl_CLSn_R, mSlepton2, mN &
!!$                  &, cpl_LNSl_L, cpl_LNSl_R, GMutoEGamma, BrMutoEGamma)
!!$   Call LtoLpGamma(3, 1, mSneutrino2, mC, cpl_CLSn_L, cpl_CLSn_R, mSlepton2, mN &
!!$                 &, cpl_LNSl_L, cpl_LNSl_R, GTautoEGamma, BrTautoEGamma)
!!$   Call LtoLpGamma(3, 2, mSneutrino2, mC, cpl_CLSn_L, cpl_CLSn_R, mSlepton2, mN &
!!$                 &, cpl_LNSl_L, cpl_LNSl_R, GTautoMuGamma, BrTautoMuGamma)
!write(*,*) "a",BrMutoEGamma,BrTautoEGamma, BrTautoMuGamma

  Do i1=1,2
   Do i2=1,3
    Do i3=1,3     
     Call CoupCharginoSfermion(i1, i2, i3, gauge_mZ(2), -0.5_dp, RSneut_T  &
             & , Y_l_mZ, Zero33C, uL_L, uL_R, U_T, V_T, cpl_CLSn_L(i1,i2,i3) &
             & , cpl_CLSn_R(i1,i2,i3) )
    End Do
   End Do
  End Do
  Do i1=1,3
   Do i2=1,4
    Do i3=1,6  
     Call CoupNeutralinoSlepton(i1, i2, i3, gauge_mZ(1), gauge_mZ(2), RSlepton_T &
         & , uL_L, uL_R, Y_l_mZ, N_T, cpl_LNSl_L(i1,i2,i3), cpl_LNSl_R(i1,i2,i3) )
    End Do
   End Do
  End Do
   Call LtoLpGamma(2, 1, mSneutrino2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlepton2_T, mN_T &
                  &, cpl_LNSl_L, cpl_LNSl_R, GMutoEGamma, BrMutoEGamma)
   Call LtoLpGamma(3, 1, mSneutrino2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlepton2_T, mN_T &
                 &, cpl_LNSl_L, cpl_LNSl_R, GTautoEGamma, BrTautoEGamma)
   Call LtoLpGamma(3, 2, mSneutrino2_T, mC_T, cpl_CLSn_L, cpl_CLSn_R, mSlepton2_T, mN_T &
                 &, cpl_LNSl_L, cpl_LNSl_R, GTautoMuGamma, BrTautoMuGamma)
!write(*,*) "b",BrMutoEGamma,BrTautoEGamma, BrTautoMuGamma
  End If
  !------------------------------------------------------------------
  ! rare decays of leptons: l -> 3 l' 
  !------------------------------------------------------------------
  BrMu3E = 0._dp
  BrTau3E = 0._dp
  BrTau3Mu = 0._dp
  gU1 = gauge_mZ(1)
  gSU2 = gauge_mZ(2)
  If (GenerationMixing) Then
!   Call BR_lj_to_3li(2, 1, gU1, gSU2, Y_l, mSlepton2, RSlepton, mN, N &
!          & , mSneut2, RSneut, mC, U, V, mS02, RS0, mP02, RP0, A_l, mu, vevSM  &
!          & , BrMu3e)
!   Call BR_lj_to_3li(3, 1, gU1, gSU2, Y_l, mSlepton2, RSlepton, mN, N &
!          & , mSneut2, RSneut, mC, U, V, mS02, RS0, mP02, RP0, A_l, mu, vevSM  &
!          & , BrTau3e)
!   Call BR_lj_to_3li(3, 2, gU1, gSU2, Y_l, mSlepton2, RSlepton, mN, N &
!          & , mSneut2, RSneut, mC, U, V, mS02, RS0, mP02, RP0, A_l, mu, vevSM  &
!          & , BrTau3Mu)
   Call BR_lj_to_3li(2, 1, gU1, gSU2, Y_l_mZ, uL_L, uL_R, mSlepton2_T, RSlepton_T   &
          & , mN_T, N_T, mSneutrino2_T, RSneut_T, mC_T, U_T, V_T, mS02_T, RS0_T     &
          & , mP02_T, RP0_T, A_l_mZ, mu_mZ, vevSM, BrMu3e)
   Call BR_lj_to_3li(3, 1, gU1, gSU2, Y_l_mZ, uL_L, uL_R, mSlepton2_T, RSlepton_T   &
          & , mN_T, N_T, mSneutrino2_T, RSneut_T, mC_T, U_T, V_T, mS02_T, RS0_T     &
          & , mP02_T, RP0_T, A_l_mZ, mu_mZ, vevSM, BrTau3e)
   Call BR_lj_to_3li(3, 2, gU1, gSU2, Y_l_mZ, uL_L, uL_R, mSlepton2_T, RSlepton_T   &
          & , mN_T, N_T, mSneutrino2_T, RSneut_T, mC_T, U_T, V_T, mS02_T, RS0_T     &
          & , mP02_T, RP0_T, A_l_mZ, mu_mZ, vevSM, BrTau3Mu)
  End If
  !-------------
  ! delta(rho)
  !-------------
  rho_parameter = DeltaRho(mZ2, mW2, mP02, RP0, mSneutrino2, RSneutrino  &
                &         , mSlepton2, RSlepton, mSup2, RSup, mSdown2    &
                &         , RSdown, mC, U, V, mN, N)
  !----------------------------------
  ! rare Z-boson decays into leptons
  !----------------------------------
  BR_Z_e_mu = 0._dp
  BR_Z_e_tau = 0._dp
  BR_Z_mu_tau = 0._dp

  If (GenerationMixing) Then
   cosW = gSU2 / Sqrt(gU1**2 + gSU2**2)
   sinW2 = gU1**2 / (gU1**2 + gSU2**2)

   cpl_SnSnZ = 0._dp
   Call CoupSneutrinoZ(gSU2, sinW2, cpl_SnSnZ(1,1))
   cpl_SnSnZ(2,2) = cpl_SnSnZ(1,1)
   cpl_SnSnZ(3,3) = cpl_SnSnZ(1,1)

   Do i1=1,6
    Do i2=1,6
     Call CoupSleptonZ(i1, i2, gSU2, sinW2, RSlepton_T, cpl_SlSlZ(i1,i2))
    End Do
   End Do

   Do i1=1,4
    Do i2=1,4
     Call CoupNeutralinoZ(i1, i2, N_T, gSU2, cosW &
                       & , cpl_NNZ_L(i1,i2), cpl_NNZ_R(i1,i2))
    End Do
   End Do

   Do i1=1,2
    Do i2=1,2
     Call CoupCharginoZ(i1, i2, U_T, V_T, gSU2, cosW &
                     & , cpl_CCZ_L(i1,i2), cpl_CCZ_R(i1,i2))
    End Do
   End Do

   Call ZtoLiLj(1, 2, .False.  &
       & , mC_T, mC2_T, cpl_CLSn_L, cpl_CLSn_R, mSneutrino2_T, cpl_SnSnZ &
       & , mN_T, mN2_T, cpl_LNSl_L, cpl_LNSl_R, mSlepton2_T, cpl_SlSlZ   &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R        &
       & , BR_Z_e_mu)

   Call ZtoLiLj(1, 3, .True.  &
       & , mC_T, mC2_T, cpl_CLSn_L, cpl_CLSn_R, mSneutrino2_T, cpl_SnSnZ &
       & , mN_T, mN2_T, cpl_LNSl_L, cpl_LNSl_R, mSlepton2_T, cpl_SlSlZ   &
 !      & , mC, mC2, cpl_CLSn_L, cpl_CLSn_R, mSneut2, cpl_SnSnZ &
 !      & , mN, mN2, cpl_LNSl_L, cpl_LNSl_R, mSlepton2, cpl_SlSlZ &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R        &
       & , BR_Z_e_tau)

   Call ZtoLiLj(2, 3, .True.  &
       & , mC_T, mC2_T, cpl_CLSn_L, cpl_CLSn_R, mSneutrino2_T, cpl_SnSnZ &
       & , mN_T, mN2_T, cpl_LNSl_L, cpl_LNSl_R, mSlepton2_T, cpl_SlSlZ   &
!       & , mC, mC2, cpl_CLSn_L, cpl_CLSn_R, mSneut2, cpl_SnSnZ &
!       & , mN, mN2, cpl_LNSl_L, cpl_LNSl_R, mSlepton2, cpl_SlSlZ &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R        &
       & , BR_Z_mu_tau)

  End If ! GenerationMixing

  !-----------------------------
  ! this still needs to be done
  !-----------------------------
  MuE_conv_ti = 0._dp

  !-------------------------------------------------------------------------
  ! neutrino masses and mixings, if the dim-5 operator has non-zero entries
  !-------------------------------------------------------------------------
  If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
   MnuL5 = MnuL5 * vevSM(2)**2
   Call NeutrinoMass_1L(MnuL5, gauge(1), gauge(2), Y_l_mZ, mC2_T, U_T, V_T, mN2_T, N_T &
           & , mSlepton2_T, Rslepton_T, mSneutrino2_T, Rsneut_T, mf_nu, Unu, kont)
  End If

 End Subroutine CalculateLowEnergyConstraints


 Subroutine CalculateSpectrum(n_run, delta, WriteOut, kont, tanb, vevSM     &
     & , mC, U, V, mN, N, mS0, mS02, RS0, mP0, mP02, RP0, mSpm, mSpm2, RSpm &
     & , mSdown, mSdown2, RSdown, mSup, mSup2, RSup, mSlepton, mSlepton2    &
     & , RSlepton, mSneut, mSneut2, RSneut, mGlu, PhaseGlu, gauge           &
     & , uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, Y_l, Y_d, Y_u                  &
     & , Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B       &
     & , m_GUT)
 Implicit None
  Integer, Intent(in) :: n_run
  Logical, Intent(in) :: WriteOut
  Real(dp), Intent(in) :: delta
  Integer, Intent(inout) :: kont
  Real(dp), Intent(inout) :: gauge(3), M2_H(2), tanb, vevSM(2), mP0(2) &
     & , mP02(2), RP0(2,2), mGlu, mSpm(2), mSpm2(2)
  Complex(dp), Intent(inout) :: Y_l(3,3), Y_d(3,3), Y_u(3,3), mu, B, Mi(3)  &
     & , A_l(3,3), A_d(3,3), A_u(3,3), M2_E(3,3), M2_L(3,3), M2_D(3,3)      &
     & , M2_U(3,3), M2_Q(3,3), PhaseGlu, RSpm(2,2)
  Real(dp), Intent(inout) :: mS0(2), mS02(2), RS0(2,2) &
     & , mC(2), mN(4), mSneut(3), mSneut2(3), mSlepton(6), mSlepton2(6) &
     & , mSdown(6), mSdown2(6), mSup(6), mSup2(6), m_GUT
  Complex(dp), Intent(inout) :: U(2,2), V(2,2), N(4,4), Rsneut(3,3)  &
     & , RSlepton(6,6), RSdown(6,6), RSup(6,6)
  Complex(dp), Intent(inout), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L &
     & , uU_R

  Real(dp) :: mZ2_t, mW2_t, Scale_Q, g, gp, dt, tz, g1(57), g7(7), vev      &
       & , mN2(4), mC2(2), g2(213), mglu_T, mP0_T(2), mP02_T(2), RP0_T(2,2)
  Real(dp) :: mC_T(2), mN_T(4), mS0_T(2), mSpm_T(2), mSup_T(6), mSdown_T(6)  &
    & , mSlepton_T(6), mSneut_T(3), mSup2_T(6), mSdown2_T(6), mSlepton2_T(6) &
    & , mSneut2_T(3), mS02_T(2), mSpm2_T(2), RS0_T(2,2), mC2_T(2), mN2_T(4)  &
    & , mass_new(32), mass_old(32), diff_m(32), tanb_in, sinW2, mf3(3)
  Real(dp) :: g6(6), g8(8), g9(9)
  Complex(dp) :: U_T(2,2), V_T(2,2), N_T(4,4), RSpm_T(2,2), RSdown_T(6,6)  &
    & , RSup_T(6,6), RSlepton_T(6,6), RSneut_T(3,3), CKM_Q(3,3)
  Integer :: i1, i2, ierr
  Logical :: Converge

  !--------------
  Real(dp) :: amz, as_5, mf_nuT(3),mf_ut(3),mf_lt(3),mf_dt(3)
  complex(dp) :: ckmt(3,3), pmnst(3,3)
  !------------------------------------------------------------------
  ! Performing a first, very rough calculation of the parameters
  ! using 1-loop RGEs and tree-level boundary conditions for gauge and
  ! Yukawa couplings at m_Z for the parameters.
  ! These parameters are used to get a first set of tree-level masses
  ! which are needed for the first of computation of the SUSY thresholds
  ! to gauge and Yukawa couplings. In case of the general MSSM the
  ! parameters are already given.
  !------------------------------------------------------------------
  kont = 0
  If ( (HighScaleModel(1:4).Eq."MSSM").Or.(HighScaleModel.Eq."pMSSM") )  Then
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

  Else ! high scale model

   If (len(trim(Old_data)).ne.0) then
    Call FirstGuess(phase_mu, tanb, Mi, M2_E, M2_L, A_l, M2_D   &
           & , M2_Q, M2_U, A_d, A_u, mu, B, M2_H, gp, g, Y_l  &
           & , Y_d, Y_u, vevSM, mP02, mP0, kont, .True., delta, Trim(Old_data) )
   Else 
    Call FirstGuess(phase_mu, tanb, Mi, M2_E, M2_L, A_l, M2_D   &
           & , M2_Q, M2_U, A_d, A_u, mu, B, M2_H, gp, g, Y_l  &
           & , Y_d, Y_u, vevSM, mP02, mP0, kont)
   End If
  End If

  If (External_Spectrum) Return ! using the externaly given spectrum

  If (HighScaleModel.Eq."MSSMtree") Then 
   If (GenerationMixing) Then
    Ad_sckm = A_d
    Au_sckm = A_u
    M2D_sckm = M2_D
    M2Q_sckm = M2_Q
    M2U_sckm = M2_U
    Call Switch_from_superCKM(Y_d, Y_u, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm &
              & , M2U_sckm, A_d, A_u, M2_D, M2_Q, M2_U, .False. )
   End If
   !------------------------------------------------------------------
   ! last flag implies that Higgs masses at 1-loop effective potential
   !------------------------------------------------------------------
   Call TreeMasses(gp, g, vevSM, Mi(1), Mi(2), Mi(3), mu, B        &
        &, tanb, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d       &
        &, A_u, Y_d, Y_u, mGlu_T, PhaseGlu, mC, mC2, U, V, mN      &
        &, mN2, N, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2    &
        &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2, RSup    &
        &, mP0_T, mP02_T, RP0_T, mS0, mS02, RS0, mSpm, mSpm2, RSpm &
        &, GenerationMixing, kont, .True.) 

   tanb_Q = tanb
   mA2_Q = mP02(2)
   vev_Q = Sqrt(vevSM(1)**2 + vevSM(2)**2)
   rho_parameter = 0

  Else 

   Call TreeMasses(gp, g, vevSM, Mi(1), Mi(2), Mi(3), mu, B        &
        &, tanb, M2_E, M2_L, A_l, Y_l, M2_D, M2_U, M2_Q, A_d       &
        &, A_u, Y_d, Y_u, mGlu_T, PhaseGlu, mC, mC2, U, V, mN      &
        &, mN2, N, mSneut, mSneut2, Rsneut, mSlepton, mSlepton2    &
        &, RSlepton, mSdown, mSdown2, RSdown, mSup, mSup2, RSup    &
        &, mP0_T, mP02_T, RP0_T, mS0, mS02, RS0, mSpm, mSpm2, RSpm &
        &, GenerationMixing, kont, .False.) ! tree-level Higgs mass

   mass_old(1) = Abs(mglu_T)
   mass_old(2:3) = Abs(mC)
   mass_old(4:7) = Abs(mN)
   mass_old(8:9) = mS0
   mass_old(10) = mP0(2)  
   mass_old(11) = mSpm(2)  
   mass_old(12:17) = mSup   
   mass_old(18:23) = mSdown   
   mass_old(24:29) = mSlepton 
   mass_old(30:32) = mSneut
  End If

  If (kont.Ne.0) Return

   ! In the SPA convention the the renormalization scale is fixed with 1 TeV
  If (SPA_Convention) Call SetRGEScale(1.e3_dp**2)
  Scale_Q = Sqrt(GetRenormalizationScale())

  If (HighScaleModel.Eq."MSSM") Then
  ! MSSM parameters and masses at loop level, all parameters are given at
  ! scale Q, mu and m_A(tree) serve as input in the Higgs sector
   tanb_Q = tanb

   ! calculate first gauge and Yukawa in DR-scheme at m_Z
   kont = 0
   Call BoundaryEW(1, vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T, mSpm  &
    & , mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup            &
    & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut                     &
    & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R, mGlu_T, PhaseGlu, mZ2_t, mW2_t &
    & , delta, g1, kont)

   converge = .False.
   Do i1=1,n_run
    tz = Log(Scale_Q/mZ)
    dt = tz / 50._dp
    Call odeint(g1, 57, 0._dp, tz, 0.1_dp*delta, dt, 0._dp, rge57, kont)
    !----------------------------------------------
    ! evolve first parameters up to Q
    !----------------------------------------------
    Call GToCouplings(g1, gauge, Y_l, Y_d, Y_u)
    gauge(1) = Sqrt(3._dp/5._dp) * gauge(1)

    A_l = ZeroC
    A_u = ZeroC
    A_d = ZeroC

    If (GenerationMixing) Then
     Y_l = Transpose(Y_l) 
     Y_d = Transpose(Y_d) 
     Y_u = Transpose(Y_u) 

     A_l = AoY_l * Y_l
     A_d = AoY_d * Y_d
     A_u = AoY_u * Y_u

    
     Call FermionMass(Y_u, sqrt2, mf3, uU_L, uU_R, ierr)
     Call FermionMass(Y_d, sqrt2, mf3, uD_L, uD_R, ierr)
     CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

     M2_D = Matmul( Matmul( Transpose(Conjg(uD_R)), M2D_sckm), uD_R)
     M2_U = Matmul( Matmul( Transpose(Conjg(uU_R)), M2U_sckm), uU_R)
     M2_Q = Matmul( Matmul( Transpose(Conjg(uU_L)), M2Q_sckm), uU_L)

    Else

     Do i2=1,3
      A_l(i2,i2) = AoY_l(i2,i2) * Y_l(i2,i2)
      A_d(i2,i2) = AoY_d(i2,i2) * Y_d(i2,i2)
      A_u(i2,i2) = AoY_u(i2,i2) * Y_u(i2,i2)
     End Do
    End If

    Call LoopMassesMSSM_3(tanb_Q, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d, A_u &
       & , M2_E, M2_L, M2_D, M2_Q, M2_U, mu, B, 0.1_dp*delta              &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0      &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2        &
       & , RSup, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2           &
       & , RSneut, mGlu, PhaseGlu, M2_H, kont)

    mass_new(1) = Abs(mglu)
    mass_new(2:3) = Abs(mC)
    mass_new(4:7) = Abs(mN)
    mass_new(8:9) = mS0
    mass_new(10) = mP0(2)  
    mass_new(11) = mSpm(2)  
    mass_new(12:17) = mSup   
    mass_new(18:23) = mSdown   
    mass_new(24:29) = mSlepton 
    mass_new(30:32) = mSneut

    diff_m = Abs(mass_new - mass_old)
    Where (mass_old.Ne.0._dp) diff_m = diff_m / mass_old

    If (Maxval(diff_m).Lt.0.1_dp*delta) Then
     converge = .True.
     Exit
    Else
     mass_old = mass_new
    End If

    gauge(1) = Sqrt(5._dp/3._dp) * gauge(1)
    Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u       &
                  & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2)

    g7(1:3) = g2(1:3) ! gauge couplings
    g7(4) = g2(20) ! tau Yukawa
    g7(5) = g2(38) ! b Yukawa
    g7(6) = g2(56) ! t Yukawa
    g7(7) = Log( tanb_Q )

    !----------------------------------------------
    ! evolve next tan(beta) down to m_Z
    !----------------------------------------------
    Call odeint(g2, 213, tz, 0._dp, 0.1_dp*delta, dt, 0._dp, rge213, kont)

    Call GToParameters(g2, gauge_mZ, Y_l_mZ, Y_d_mZ, Y_u_mZ, Mi_mZ, A_l_mZ &
       & , A_d_mZ, A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ     &
       & , M2_H_mZ, mu_mZ, B_mZ)

    Call odeint(g7, 7, tz, 0._dp, 0.1_dp*delta, dt, 0._dp, rge7, kont)

    tanb_mZ = Exp( g7(7) )

    sinW2 = 1._dp - mW2 / mZ2
    vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
    vevSM(1) = vev / Sqrt(1._dp + tanb_mZ**2)
    vevSM(2) = tanb_mZ * vevSM(1)
    g2(1) = Sqrt(3._dp/5._dp) * g2(1)

    Call TreeMasses(g2(1), g2(2), vevSM, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3), mu_mZ &
      & , B_mZ, tanb_mZ, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ, M2_D_mZ, M2_U_mZ  &
      & , M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ, mGlu_T, PhaseGlu, mC_T    &
      & , mC2_T, U_T, V_T, mN_T, mN2_T, N_T, mSneut_T, mSneut2_T, Rsneut_T   &
      & , mSlepton_T, mSlepton2_T, RSlepton_T, mSdown_T, mSdown2_T, RSdown_T &
      & , mSup_T, mSup2_T, RSup_T, mP0_T, mP02_T, RP0_T, mS0_T, mS02_T       &
      & , RS0_T, mSpm_T, mSpm2_T, RSpm_T, GenerationMixing, kont, .False.)

    If (Min(Minval(mSup2_T), Minval(mSdown2_T), Minval(mSlepton2_T)   &
       &   , Minval(mSneut2_T), Minval(mS02_T), Minval(mP02_T)       &
       &   , Minval(mSpm2_T)).Gt. 0._dp ) Then
     mC = mC_T
     mC2 = mC2_T
     U = U_T
     V = V_T
     mN = mN_T
     mN2 = mN2_T
     N = N_T
     mSneut = mSneut_T
     mSneut2 = mSneut2_T
     Rsneut = Rsneut_T
     mSlepton = mSlepton_T
     mSlepton2 = mSlepton2_T
     RSlepton = RSlepton_T
     mSDown = mSDown_T
     mSDown2 = mSDown2_T
     RSDown = RSDown_T
     mSup = mSup_T
     mSup2 = mSup2_T
     RSup = RSup_T
     mS0 = mS0_T
     mS02 = mS02_T
     RS0 = RS0_T
     mSpm = mSpm_T
     mSpm2 = mSpm2_T
     RSpm = RSpm_T
     YukScen = 1 ! using running masses for boundary conditions at mZ
    Else
     mGlu_T = mGlu
     mP0_T = mP0
     mP02_T = mP02
     RP0_T = RP0
     YukScen = 2 ! using pole masses for boundary conditions at mZ
    End If

    kont = 0
    Call BoundaryEW(i1+1,vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
      & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup   &
      & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R      &
      & , uD_L, uD_R , uL_L, uL_R, mGlu_T, PhaseGlu, mZ2_t, mW2_t, delta  &
      & , g1, kont)

   End Do

   If ((kont.eq.0).and.(.not.converge)) then
    Write (ErrCan,*) 'Problem in subroutine CalculateSpectrum!!'
    Write (ErrCan,*) "After",n_run,"iterations no convergence found"
    kont = -1200
   End If

  Else If (HighScaleModel.Eq."MSSM1") Then
  ! MSSM parameters and masses at loop level, all parameters are given at
  ! scale Q, m^2_(H_i) serve as input in the Higgs sector
   tanb_Q = tanb

   ! calculate first gauge and Yukawa in DR-scheme at m_Z
   kont = 0
   Call BoundaryEW(1, vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T, mSpm  &
    & , mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup            &
    & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut                     &
    & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R, mGlu_T, PhaseGlu, mZ2_t, mW2_t &
    & , delta, g1, kont)

   converge = .False.
   tanb_mZ = tanb ! first gues, justified because tanb runs weakly
   Do i1=1,n_run
    !----------------------------------------------
    ! evolve first parameters up to Q
    !----------------------------------------------
    tz = Log(Scale_Q/mZ)
    dt = tz / 50._dp
    Call odeint(g1, 57, -tz, 0._dp, 0.1_dp*delta, -dt, 0._dp, rge57, kont)
    Call GToCouplings(g1, gauge, Y_l, Y_d, Y_u)
    gauge(1) = Sqrt(3._dp/5._dp) * gauge(1)

    If (SPA_Convention) Then
     tanb_in = tanb_Q
    Else
     tanb_in = tanb_mZ
    End If

    If (GenerationMixing) Then
     Y_l = Transpose(Y_l) 
     Y_d = Transpose(Y_d) 
     Y_u = Transpose(Y_u) 
     
     Call Switch_from_superCKM(Y_d, Y_u, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm &
               & , M2U_sckm, A_d, A_u, M2_D, M2_Q, M2_U, .False. )
    Else ! .not. GenerationMixing
     A_l = ZeroC
     A_u = ZeroC
     A_d = ZeroC
     Do i2=1,3
      A_l(i2,i2) = AoY_l(i2,i2) * Y_l(i2,i2)
      A_d(i2,i2) = AoY_d(i2,i2) * Y_d(i2,i2)
      A_u(i2,i2) = AoY_u(i2,i2) * Y_u(i2,i2)
     End Do
    End If

    Call LoopMassesMSSM(delta, tanb_in, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d  &
       & , A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, phase_mu, mu, B, i1   &
       & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R                             &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP0, mP02, RP0      &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2        &
       & , RSup, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2           &
       & , RSneut, mGlu, PhaseGlu, kont)

    mass_new(1) = Abs(mglu)
    mass_new(2:3) = Abs(mC)
    mass_new(4:7) = Abs(mN)
    mass_new(8:9) = mS0
    mass_new(10) = mP0(2)  
    mass_new(11) = mSpm(2)  
    mass_new(12:17) = mSup   
    mass_new(18:23) = mSdown   
    mass_new(24:29) = mSlepton 
    mass_new(30:32) = mSneut

    diff_m = Abs(mass_new - mass_old)
    Where (mass_old.Ne.0._dp) diff_m = diff_m / mass_old

    If (Maxval(diff_m).Lt.0.1_dp*delta) Then
     converge = .True.
     Exit
    Else
     mass_old = mass_new
    End If

    gauge(1) = Sqrt(5._dp/3._dp) * gauge(1)
    Y_l = Transpose(Y_l)
    Y_d = Transpose(Y_d)
    Y_u = Transpose(Y_u)
    A_l = Transpose(A_l)
    A_d = Transpose(A_d)
    A_u = Transpose(A_u)
    Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u       &
                  & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2)
    !------------------------------
    ! now the running
    !------------------------------
    Call odeint(g2, 213, tz, 0._dp,  0.1_dp*delta, dt, 0._dp, rge213, kont)
    Call GToParameters(g2, gauge_mZ, Y_l_mZ, Y_d_mZ, Y_u_mZ, Mi_mZ, A_l_mZ &
       & , A_d_mZ, A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ     &
       & , M2_H_mZ, mu_mZ, B_mZ)
    Y_l = Transpose(Y_l)
    Y_d = Transpose(Y_d)
    Y_u = Transpose(Y_u)
    A_l = Transpose(A_l)
    A_d = Transpose(A_d)
    A_u = Transpose(A_u)
    Y_l_mZ = Transpose(Y_l_mZ)
    Y_d_mZ = Transpose(Y_d_mZ)
    Y_u_mZ = Transpose(Y_u_mZ)
    A_l_mZ = Transpose(A_l_mZ)
    A_d_mZ = Transpose(A_d_mZ)
    A_u_mZ = Transpose(A_u_mZ)
    !-----------------------------------------------------
    ! tanb_mZ is has been calculated in routine LoopMasses
    !-----------------------------------------------------
    sinW2 = 1._dp - mW2 / mZ2
    vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
    vevSM(1) = vev / Sqrt(1._dp + tanb_mZ**2)
    vevSM(2) = tanb_mZ * vevSM(1)
    g2(1) = Sqrt(3._dp/5._dp) * g2(1)
    
   Call TreeMassesMSSM2(g2(1), g2(2), vevSM, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3)   &
      & , mu_mZ, B_mZ , tanb_mZ, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ, M2_D_mZ   &
      & , M2_U_mZ, M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ, uU_L, uU_R       &
      & , uD_L, uD_R, uL_L, uL_R, mGlu_T, PhaseGlu, mC_T, mC2_T, U_T, V_T    &
      & , mN_T, mN2_T, N_T, mSneut_T, mSneut2_T, Rsneut_T, mSlepton_T        &
      & , mSlepton2_T, RSlepton_T, mSdown_T, mSdown2_T, RSdown_T, mSup_T     &
      & , mSup2_T, RSup_T, mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T        &
      & , mSpm_T, mSpm2_T, RSpm_T, mZ2_t, mW2_t, GenerationMixing, kont      &
      & , .False., .False.)

    If  (kont.Ne.0)  Return

    If (Min(Minval(mSup2_T), Minval(mSdown2_T), Minval(mSlepton2_T)  &
       &   , Minval(mSneut2_T), Minval(mS02_T), Minval(mP02_T)       &
       &   , Minval(mSpm2_T)).Gt. 0._dp ) Then
     mC = mC_T
     mC2 = mC2_T
     U = U_T
     V = V_T
     mN = mN_T
     mN2 = mN2_T
     N = N_T
     mSneut = mSneut_T
     mSneut2 = mSneut2_T
     Rsneut = Rsneut_T
     mSlepton = mSlepton_T
     mSlepton2 = mSlepton2_T
     RSlepton = RSlepton_T
     mSDown = mSDown_T
     mSDown2 = mSDown2_T
     RSDown = RSDown_T
     mSup = mSup_T
     mSup2 = mSup2_T
     RSup = RSup_T
     mS0 = mS0_T
     mS02 = mS02_T
     RS0 = RS0_T
     mSpm = mSpm_T
     mSpm2 = mSpm2_T
     RSpm = RSpm_T
     YukScen = 1 ! using running masses for boundary conditions at mZ
    Else
     YukScen = 2 ! using pole masses for boundary conditions at mZ
     mglu_T = mglu
     mP0_T = mP0
     mP02_T = mP02
     RP0_T = RP0     
    End If

    kont = 0
    Call BoundaryEW(i1+1,vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
      & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup   &
      & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R      &
      & , uD_L, uD_R , uL_L, uL_R, mGlu_T, PhaseGlu, mZ2_t, mW2_t, delta  &
      & , g1, kont)

   End Do

   If ((kont.eq.0).and.(.not.converge)) then
    Write (ErrCan,*) 'Problem in subroutine CalculateSpectrum, model ' &
                              //trim(HighScaleModel)//'!!'
    Write (ErrCan,*) "After",n_run,"iterations no convergence found"
    kont = -1200
   End If

  Else If (HighScaleModel.Eq."pMSSM") Then

   ! calculate first gauge and Yukawa in DR-scheme at m_Z
   Call BoundaryEW(1, vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T, mSpm  &
    & , mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup            &
    & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut                     &
    & , uU_L, uU_R, uD_L, uD_R, uL_L, uL_R, mGlu_T, PhaseGlu, mZ2_t, mW2_t &
    & , delta, g1, kont)

   converge = .False.
   mP0(1) = mZ
   mP02(1) = mZ2
   Do i1=1,n_run
    tz = Log(Scale_Q/mZ)
    dt = tz / 50._dp
    kont = 0
    Call odeint(g1, 57, 0._dp, tz,  0.1_dp*delta, dt, 0._dp, rge57, kont)

    !----------------------------------------------
    ! evolve first parameters up to Q
    !----------------------------------------------
    Call GToCouplings(g1, gauge, Y_l, Y_d, Y_u)
    gauge(1) = Sqrt(3._dp/5._dp) * gauge(1)

    If (GenerationMixing) Then
     Y_l = Transpose(Y_l) 
     Y_d = Transpose(Y_d) 
     Y_u = Transpose(Y_u) 
     
     Call Switch_from_superCKM(Y_d, Y_u, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm &
               & , M2U_sckm, A_d, A_u, M2_D, M2_Q, M2_U, .False. )
    Else ! .not. GenerationMixing
     A_l = ZeroC
     A_u = ZeroC
     A_d = ZeroC
     Do i2=1,3
      A_l(i2,i2) = AoY_l(i2,i2) * Y_l(i2,i2)
      A_d(i2,i2) = AoY_d(i2,i2) * Y_d(i2,i2)
      A_u(i2,i2) = AoY_u(i2,i2) * Y_u(i2,i2)
     End Do
    End If

    kont = 0

    Call LoopMassesMSSM_2(delta, tanb, gauge, Y_l, Y_d, Y_u, Mi, A_l, A_d &
       & , A_u, M2_E, M2_L, M2_D, M2_Q, M2_U, mu                          &
       & , mC, mC2, U, V, mN, mN2, N, mS0, mS02, RS0, mP02, RP0           &
       & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2        &
       & , RSup, mSlepton, mSlepton2, RSlepton, mSneut, mSneut2           &
       & , RSneut, mGlu, PhaseGlu, M2_H, B, kont)

    mass_new(1) = Abs(mglu)
    mass_new(2:3) = Abs(mC)
    mass_new(4:7) = Abs(mN)
    mass_new(8:9) = mS0
    mass_new(10) = mP0(2)  
    mass_new(11) = mSpm(2)  
    mass_new(12:17) = mSup   
    mass_new(18:23) = mSdown   
    mass_new(24:29) = mSlepton 
    mass_new(30:32) = mSneut
    diff_m = Abs(mass_new - mass_old)
    Where (mass_old.Ne.0._dp) diff_m = diff_m / mass_old

    If (Maxval(diff_m).Lt.delta) Then
     converge = .True.
     Exit
    Else
     mass_old = mass_new
    End If

    If (i1.Eq.50) Then
     Write(*,*) "Problem with accuracy (pMSSM)",diff_m,delta,i1
     Exit
    End If

    gauge(1) = Sqrt(5._dp/3._dp) * gauge(1)

    Y_l = Transpose(Y_l) 
    Y_d = Transpose(Y_d) 
    Y_u = Transpose(Y_u) 
    A_l = Transpose(A_l) 
    A_d = Transpose(A_d) 
    a_u = Transpose(A_u) 
  
   Call ParametersToG(gauge, y_l, y_d, y_u, Mi, A_l, A_d, A_u       &
                  & , M2_E, M2_L, M2_D, M2_Q, M2_U, M2_H, mu, B, g2)
    Call odeint(g2, 213, tz, 0._dp, 0.1_dp*delta, dt, 0._dp, rge213, kont)

    Call GToParameters(g2, gauge_mZ, Y_l_mZ, Y_d_mZ, Y_u_mZ, Mi_mZ, A_l_mZ &
       & , A_d_mZ, A_u_mZ, M2_E_mZ, M2_L_mZ, M2_D_mZ, M2_Q_mZ, M2_U_mZ     &
       & , M2_H_mZ, mu_mZ, B_mZ)

    Y_l_mZ = Transpose(Y_l_mZ) 
    Y_d_mZ = Transpose(Y_d_mZ) 
    Y_u_mZ = Transpose(Y_u_mZ) 
    A_l_mZ = Transpose(A_l_mZ) 
    A_d_mZ = Transpose(A_d_mZ) 
    a_u_mZ = Transpose(A_u_mZ)
 
    sinW2 = 1._dp - mW2 / mZ2
    vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
    vevSM(1) = vev / Sqrt(1._dp + tanb**2)
    vevSM(2) = tanb * vevSM(1)
    g2(1) = Sqrt(3._dp/5._dp) * g2(1)
    
    Call TreeMasses(g2(1), g2(2), vevSM, Mi_mZ(1), Mi_mZ(2), Mi_mZ(3), mu_mZ  &
      & , B_mZ , tanb, M2_E_mZ, M2_L_mZ, A_l_mZ, Y_l_mZ, M2_D_mZ, M2_U_mZ     &
      & , M2_Q_mZ, A_d_mZ, A_u_mZ, Y_d_mZ, Y_u_mZ, mGlu_T, PhaseGlu, mC_T     &
      & , mC2_T, U_T, V_T, mN_T, mN2_T, N_T, mSneut_T, mSneut2_T, Rsneut_T    &
      & , mSlepton_T, mSlepton2_T, RSlepton_T, mSdown_T, mSdown2_T, RSdown_T  &
      & , mSup_T, mSup2_T, RSup_T, mP0_T, mP02_T, RP0_T, mS0_T, mS02_T, RS0_T &
      & , mSpm_T, mSpm2_T, RSpm_T, GenerationMixing, kont, .False.)

    If (Min(Minval(mSup2_T), Minval(mSdown2_T), Minval(mSlepton2_T)   &
       &   , Minval(mSneut2_T), Minval(mS02_T), Minval(mP02_T)       &
       &   , Minval(mSpm2_T)).Gt. 0._dp ) Then
     mC = mC_T
     mC2 = mC2_T
     U = U_T
     V = V_T
     mN = mN_T
     mN2 = mN2_T
     N = N_T
     mSneut = mSneut_T
     mSneut2 = mSneut2_T
     Rsneut = Rsneut_T
     mSlepton = mSlepton_T
     mSlepton2 = mSlepton2_T
     RSlepton = RSlepton_T
     mSDown = mSDown_T
     mSDown2 = mSDown2_T
     RSDown = RSDown_T
     mSup = mSup_T
     mSup2 = mSup2_T
     RSup = RSup_T
     mS0 = mS0_T
     mS02 = mS02_T
     RS0 = RS0_T
     mSpm = mSpm_T
     mSpm2 = mSpm2_T
     RSpm = RSpm_T
     YukScen = 1 ! using running masses for boundary conditions at mZ
    Else
     YukScen = 2 ! using pole masses for boundary conditions at mZ
     mglu_T = mglu
     mP0_T = mP0
     mP02_T = mP02
     RP0_T = RP0     
    End If

    kont = 0
    Call BoundaryEW(i1+1,vevSM, mC, U, V, mN, N, mS02, RS0, mP02_T, RP0_T &
      & , mSpm, mSpm2, RSpm, mSdown, mSdown2, RSdown, mSup, mSup2, RSup   &
      & , mSlepton, mSlepton2, RSlepton, mSneut2, RSneut, uU_L, uU_R      &
      & , uD_L, uD_R , uL_L, uL_R, mGlu_T, PhaseGlu, mZ2_t, mW2_t, delta  &
      & , g1, kont)

   End Do

   If ((kont.eq.0).and.(.not.converge)) then
    Write (ErrCan,*) 'Problem in subroutine CalculateSpectrum, model ' &
                              //trim(HighScaleModel)//'!!'
    Write (ErrCan,*) "After",n_run,"iterations no convergence found"
    kont = -1200
   End If

  Else If (HighScaleModel.Eq."MSSMtree") Then
   mP0 = mP0_T
   mP02 = mP02_T
   RP0 = RP0_T
   mglu = mGlu_T
   If (GenerationMixing) Then

    Call FermionMass(Y_l,vevSM(1),mf_l_mZ,uL_L,uL_R,kont)
    Call FermionMass(Y_d,vevSM(1),mf_d_mZ,uD_L,uD_R,kont)
    Call FermionMass(Y_u,vevSM(2),mf_u_mZ,uU_L,uU_R,kont)

   Else
    uL_L = id3C
    uL_R = id3C
    uD_L = id3C
    uD_R = id3C
    uU_L = id3C
    uU_R = id3C
   End If
   
  Else ! high scale models
   mP0 = mP0_T
   mP02 = mP02_T
   RP0 = RP0_T
   mglu = mGlu_T
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

   Call Sugra(delta, vevSM, mC, U, V, mN, N, mS0, mS02, RS0 &
     & , mP0, mP02, RP0, mSpm, mSpm2, RSpm, mSdown, mSdown2 &
     & , RSdown, mSup, mSup2, RSup, mSlepton, mSlepton2     &
     & , RSlepton, mSneut, mSneut2, RSneut, mGlu, PhaseGlu  &
     & , gauge, uL_L, uL_R, uD_L, uD_R, uU_L, uU_R, Y_l     &
     & , Y_d, Y_u, Mi, A_l, A_d, A_u, M2_E, M2_L, M2_D      &
     & , M2_Q, M2_U, M2_H, mu, B, m_GUT, kont, WriteOut, n_run)
  End If

 End Subroutine CalculateSpectrum


 Subroutine ReadingData(kont)
 !--------------------------------------------------------
 ! reading the input, all routines used can be found in
 ! InputOutput.f90 
 !--------------------------------------------------------
 Implicit None
  Integer, Intent(out) :: kont

  Logical :: file_exists

  kont = -123456

 !------------------------------------------------------------------
 ! Checked, if the file LesHouches.in exists
 !------------------------------------------------------------------
  Inquire(file="LesHouches.in",exist=file_exists)
 !---------------------------------------
 !   if yes, use the data from this file
 !---------------------------------------
  If (file_exists) Then
   kont = 1

   Call LesHouches_Input(kont, HighScaleModel, Ecms, Pm, Pp, ISR, F_GMSB)

   LesHouches_Format = .True.
  Else
   LesHouches_Format = .False.
   !-----------------------------------
   ! perform the original squence
   !-----------------------------------
   !--------------------------------------------------------------------------
   ! Intialisation routines
   ! InitialzeControl: initialze the control system and opens the output file
   !           can be influenced using the file Control.in 
   ! InitializeLoopFunctions: calculates the constants needed for the 
   !            calculation of the loop functions
   ! InitializeStandardModel: sets the Standard Model parameters; in the case
   !            that the file StandardModel.in is present the data will be read
   !            in, otherwise default values will be used as described in the
   !            manual
   !--------------------------------------------------------------------------
   Call InitializeControl(11, "SPheno3.out", "SPheno3")
   Call InitializeLoopFunctions
   Call InitializeStandardModel
   If (l_CS) Call InitializeCrossSections(Ecms, Pm, Pp, ISR)
   !------------------------------------------------------------------
   ! reading the data for the specification of the high scale model,
   ! source code of the subroutine is in the module InputOutput
   !------------------------------------------------------------------
   Inquire(file="HighScale.in",exist=file_exists)

   If (file_exists) Then
    kont = 0
    Call HighScaleInput()
   Else
    Write(*,*) &
     & "Neither the file 'HighScale.in' nor the file 'LesHouches.in'"
    Write(*,*) "has been found. Please provide an input file."
    Call TerminateProgram 
   End If
  End If

 End Subroutine ReadingData

End Program SPheno

