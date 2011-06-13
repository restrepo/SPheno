Module InputOutput
!--------------------------------------------------------------
! this module contains routines for input/output
!--------------------------------------------------------------

Use Control
Use Experiment
Use LHC_observables
Use Model_Data
Use SugraRuns

Interface WriteMatrixBlock
 Module Procedure WriteMatrixBlockC, WriteMatrixBlockR
End Interface

Character(len=15) :: HighScaleModel
! input/output according to SUSY Les Houches accord
Logical, Save :: LesHouches_Format
! if R-parity is added at low energies
Logical, Save ::  Add_Rparity = .False.  
! transfer of GMSB info
Real(dp), Save :: grav_fac = 1._dp
! used in combination with Fittino
character(len=80), save :: Old_Data=""
! using 1st SLHA2 output with flavour ordered states
Logical, Save, Private :: Use_Flavour_States = .False.
! branching ratios larger than BrMin are written out
Real(dp), Save, Private :: BrMin=1.e-4_dp
! cross sections larger than SigMin [in fb] are written out
Real(dp), Save, Private :: SigMin=1.e-3_dp
! contains information on possible inconsitencies in the input
Integer, Save, Private :: in_kont(2)
! version number
Character(len=8), Save, Private :: version="v3beta48"
! name of 'input-program'
Character(len=40), Private :: sp_info 
! tempory variables for Higgs mixing in case of NMSSM
Real(dp), Private, Dimension(3,3) :: RS03_save, RP03_save
! for Les Houches input/output
Logical :: Write_SLHA1 = .False. ! write a second SLHA output file
                                 ! using SLHA1 standard only called SPheno_1.spc
Integer, Private :: i_cpv=0
Logical, Private :: l_RP_Pythia = .False. ! Pythia only takes 4x4 matrix 
                                         ! for neutralinos and 2x2 for charginos
Logical, Private :: LWrite_LHC_Observables = .False. ! give LHC observables in the output

Contains



 Subroutine HighScaleInput() 
 !--------------------------------------------------------------
 ! specifies high scale model, reads the data from HighScale.in
 ! implemented models: mSUGRA, 
 !--------------------------------------------------------------
 Implicit None

  Integer :: i1, YukawaScheme
  Real(dp) :: M0_in, M0_ii_in(3), t_in, M_hlf_in, A0_in, sign_mu, vev, test &
    &  , Scale_Q, sinW2, sinb2, cosb2, cos2b
  Logical :: check

  Iname = Iname + 1
  NameOfUnit(Iname) = "HighScaleInput"

  Open(93,file="HighScale.in",status="old")

  Read(93,*) HighScaleModel
  HighScaleModel = Adjustl(HighScaleModel)

  AoY_d_0 = 0._dp
  AoY_l_0 = 0._dp
  AoY_nu_0 = 0._dp
  AoY_u_0 = 0._dp
  Mi_0 = 0._dp
  M2_D_0 = 0._dp
  M2_E_0 = 0._dp
  M2_L_0 = 0._dp
  M2_R_0 = 0._dp
  M2_Q_0 = 0._dp
  M2_U_0 = 0._dp
  M2_H_0 = 0._dp

  !-------------------------------------------------------
  ! necessary to exclude right handed neutrinos from RGEs
  ! is set positive in the corresponding model
  !-------------------------------------------------------
  MNuR = - 1.e-9_dp

  !-----------------------------------------------------------------------
  ! these variables are only used in GMSB and will be set correctly below
  !-----------------------------------------------------------------------
  F_GMSB = 1.e12_dp
  m32 = 1.e20_dp ! set an abitrary large gravitino mass

  If (HighScaleModel.Eq."mSugra") Then
   !-----------------------------------------------------------
   ! needed for the setting of the correct boundary conditions
   !-----------------------------------------------------------
   check = SetHighScaleModel("SUGRA")

   If (.Not.check) Then
    Write(ErrCan,*) "Strange error in subroutine HighScaleInput."
    Write(ErrCan,*) "Unable to put name of high scale model, aborting."
    Call TerminateProgram
   End If

   Read(93,*) M_hlf_in
   Read(93,*) M0_in
   Read(93,*) A0_in
   Read(93,*) tanb
   Read(93,*) sign_mu
   Read(93,*) TwoLoopRGE

   Mi_0 = M_hlf_in
   Do i1=1,3
    AoY_d_0(i1,i1) = A0_in
    M2_D_0(i1,i1) = M0_in**2
   End Do
   AoY_l_0 = AoY_d_0
   AoY_u_0 = AoY_d_0
   M2_E_0 = M2_D_0
   M2_L_0 = M2_D_0
   M2_Q_0 = M2_D_0
   M2_U_0 = M2_D_0

   M2_H_0 = M0_in**2

  Else If (HighScaleModel.Eq."mSugraNuR") Then
   !-----------------------------------------------------------
   ! needed for the setting of the correct boundary conditions
   !-----------------------------------------------------------
   check = SetHighScaleModel("SUGRA_NuR")

   If (.Not.check) Then
    Write(ErrCan,*) "Strange error in subroutine HighScaleInput."
    Write(ErrCan,*) "Unable to put name of high scale model, aborting."
    Call TerminateProgram
   End If

   Read(93,*) M_hlf_in
   Read(93,*) M0_in
   Read(93,*) A0_in
   Read(93,*) tanb
   Read(93,*) sign_mu
   Read(93,*) MNuR
   Read(93,*) mf_nu
   Read(93,*) TwoLoopRGE

   Mi_0 = M_hlf_in
   Do i1=1,3
    AoY_d_0(i1,i1) = A0_in
    M2_D_0(i1,i1) = M0_in**2
   End Do
   AoY_l_0 = AoY_d_0
   AoY_nu_0 = AoY_d_0
   AoY_u_0 = AoY_d_0
   M2_E_0 = M2_D_0
   M2_R_0 = M2_D_0
   M2_L_0 = M2_D_0
   M2_Q_0 = M2_D_0
   M2_U_0 = M2_D_0

   M2_H_0 = M0_in**2

  Else If (HighScaleModel.Eq."SU5") Then
   HighScaleModel = "SUGRA_SU5"
   !-----------------------------------------------------------
   ! needed for the setting of the correct boundary conditions
   !-----------------------------------------------------------
   check = SetHighScaleModel("SUGRA_SU5")

   If (.Not.check) Then
    Write(ErrCan,*) "Strange error in subroutine HighScaleInput."
    Write(ErrCan,*) "Unable to put name of high scale model, aborting."
    Call TerminateProgram
   End If

   Read(93,*) M_hlf_in
   Read(93,*) M0_in
   Read(93,*) A0_in
   Read(93,*) tanb
   Read(93,*) sign_mu
   Read(93,*) MNuR
   Read(93,*) mf_nu
   Read(93,*) lam_0
   Read(93,*) lamp_0
   Read(93,*) M_SO_10
   Read(93,*) D_SO_10
   Read(93,*) TwoLoopRGE

   Mi_0 = M_hlf_in

   Do i1=1,3
    AoY_d_0(i1,i1) = A0_in
    M2_D_0(i1,i1) = M0_in**2
   End Do
   AoY_l_0 = AoY_d_0
   AoY_nu_0 = AoY_d_0
   AoY_u_0 = AoY_d_0
   M2_E_0 = M2_D_0
   M2_R_0 = M2_D_0
   M2_L_0 = M2_D_0
   M2_Q_0 = M2_D_0
   M2_U_0 = M2_D_0

   M2_H_0 = M0_in**2

  Elseif (HighScaleModel.Eq."AMSB") Then
   !-----------------------------------------------------------
   ! needed for the setting of the correct boundary conditions
   !-----------------------------------------------------------
   check = SetHighScaleModel("AMSB")

   If (.Not.check) Then
    Write(ErrCan,*) "Strange error in subroutine HighScaleInput."
    Write(ErrCan,*) "Unable to put name of high scale model, aborting."
    Call TerminateProgram
   End If

   Read(93,*) M_32
   Read(93,*) M0_amsb
   Read(93,*) tanb
   Read(93,*) sign_mu
   Read(93,*) TwoLoopRGE

  Else If (HighScaleModel.Eq."GMSB") Then
   !-----------------------------------------------------------
   ! needed for the setting of the correct boundary conditions
   !-----------------------------------------------------------
   check = SetHighScaleModel("GMSB")
   Read(93,*) Lambda
   Read(93,*) MlambdaS
   Read(93,*) n5plets
   Read(93,*) n10plets
   Read(93,*) A0_in
   Read(93,*) tanb
   Read(93,*) sign_mu
   Read(93,*) TwoLoopRGE
   Do i1=1,3
    AoY_d_0(i1,i1) = A0_in
   End Do
   AoY_l_0 = AoY_d_0
   AoY_u_0 = AoY_d_0
 
   F_GMSB = Lambda * MlambdaS ! needed for the calculation of NLSP
   m32 = 2.4e-10_dp * F_GMSB  ! gravitino mass in eV  

  Else If (     (HighScaleModel.Eq."String_OIIa")  &
          & .Or.(HighScaleModel.Eq."String_OIIb")  &
          & .Or.(HighScaleModel.Eq."String_OI")    ) Then

   check = .False.
   If (HighScaleModel.Eq."String_OIIa") check = SetHighScaleModel("Str_A")
   If (HighScaleModel.Eq."String_OIIb") check = SetHighScaleModel("Str_B")
   If (HighScaleModel.Eq."String_OI") check = SetHighScaleModel("Str_C")

   If (.Not.check) Then
    Write(ErrCan,*) "Strange error in subroutine HighScaleInput."
    Write(ErrCan,*) "Unable to put  name of high scale model, aborting."
    Call TerminateProgram
   End If

   Read(93,*) m32      !   gravitino mass
   Read(93,*) t_in     ! vacuum value of the moduli
   Read(93,*) g_s2     ! string coupling squared
   Read(93,*) sinT2    ! sin(theta) squared
   Read(93,*) delta_GS ! delta_GS

   If (HighScaleModel.Eq."String_OI") Then
    ! first index if for t_i, second for generation
    Read(93,*) nE_ai(1,1), nL_ai(1,1)                
    Read(93,*) nD_ai(1,1), nU_ai(1,1), nQ_ai(1,1)
    Read(93,*) nH1_ai(1), nH2_ai(1)
    nE_ai(1,2:3) = nE_ai(1,1)
    nL_ai(1,2:3) = nL_ai(1,1)
    nD_ai(1,2:3) = nD_ai(1,1)
    nU_ai(1,2:3) = nU_ai(1,1)
    nQ_ai(1,2:3) = nQ_ai(1,1)
   Else
    nE_ai(1,1:3) = -1
    nL_ai(1,1:3) = -1
    nD_ai(1,1:3) = -1
    nU_ai(1,1:3) = -1
    nQ_ai(1,1:3) = -1
    nH1_ai(1) = -1
    nH2_ai(1) = -1
   End If

   Read(93,*) tanb
   Read(93,*) sign_mu
   Read(93,*) TwoLoopRGE

   ! with some elder compilers this gives numerically more stable results
   If (SinT2.Eq.1._dp) Then 
    cosT2 = 0._dp 
    cosT = 0._dp
    SinT = 1._dp
   Else If (SinT2.Eq.0._dp) Then
    cosT2 = 1._dp 
    cosT = 1._dp
    SinT = 0._dp
   Else
    cosT2 = 1._dp - sinT2
    cosT = Sqrt(cosT2)
    SinT = Sqrt(SinT2)
   End If

   ! dilaton field, assuming k=-ln(s+\bar{s}) and 2 g_s^2 = < s+\bar{s} >
   ! for the moment everything is reel
   g_s = Sqrt(g_s2)
   k_s = - 0.5_dp * g_s2
   k_sb = k_s
   k_ss = k_s**2
   oosqrt_k_ss = 1._dp / Abs(k_s)
   phase_s = 1._dp

   ! moduli fields, assuming for the moment that all moduli fields couple
   ! equally, thus I can express everything with the help of one field
   ! for the moment everything is reel
   num_t = 1
   Do i1=1,num_t ! for a latter extension to multiple moduli fields
    t(i1) = t_in
    ThetaA(i1) = 1._dp
    ReT(i1) = Real(t(i1), dp)
    Phase_T(i1) = t(i1) / Abs( t(i1) )
    G2t(i1) = Eisenstein( t(i1) )
    ReG2ThetaT(i1) = 2._dp * ReT(i1) * G2t(i1) * ThetaA(i1)
    LnReDedekind(i1) = Log( ReT(i1) * Abs( Dedekind( t(i1) ) )**4 )
   End Do

   ! sums needed in the calculation of the gaugino masses
   ! U(1)
   SumC_O1(1) = -0.6_dp * (3 + Sum( nE_ai(1,:)) )                            &
              & -0.3_dp * (3 + Sum( nL_ai(1,:)) + nH1_ai(1) + nH2_ai(1) )    &
              & -0.2_dp * (3 + Sum( nD_ai(1,:)) )                            &
              & -0.8_dp * (3 + Sum( nU_ai(1,:)) )                            &
              & -0.1_dp * (3 + Sum( nQ_ai(1,:)) )
   SumC1(1) = -0.6_dp * (9 + Sum( nE_ai(1,:)) )                            &
            & -0.3_dp * (9 + Sum( nL_ai(1,:)) + nH1_ai(1) + nH2_ai(1) )    &
            & -0.2_dp * (9 + Sum( nD_ai(1,:)) )                            &
            & -0.8_dp * (9 + Sum( nU_ai(1,:)) )                            &
            & -0.1_dp * (9 + Sum( nQ_ai(1,:)) )
   SumC2(1) = 6.6_dp
   ! SU(2)
   SumC_O1(2) = -0.5_dp * (3 + Sum( nL_ai(1,:)) + nH1_ai(1) + nH2_ai(1) )    &
              & -1.5_dp * (3 + Sum( nQ_ai(1,:)) )
   SumC1(2) = -0.5_dp * (9 + Sum( nL_ai(1,:)) + nH1_ai(1) + nH2_ai(1) )    &
            & -1.5_dp * (9 + Sum( nQ_ai(1,:)) )
   SumC2(2) = 5._dp
   ! SU(3)
   SumC_O1(3) = -0.5_dp * (3 + Sum( nU_ai(1,:)) + Sum( nD_ai(1,:)) )    &
              & -1._dp * (3 + Sum( nQ_ai(1,:)) )
   SumC1(3) = -0.5_dp * (9 + Sum( nU_ai(1,:)) + Sum( nD_ai(1,:)) )    &
            & -1._dp * (9 + Sum( nQ_ai(1,:)) )
   SumC2(3) = 3._dp

  Else If (HighScaleModel.Eq."Oscar") Then
   !-----------------------------------------------------------
   ! needed for the setting of the correct boundary conditions
   !-----------------------------------------------------------
   check = SetHighScaleModel("Oscar")

   If (.Not.check) Then
    Write(ErrCan,*) "Strange error in subroutine HighScaleInput."
    Write(ErrCan,*) "Unable to put name of high scale model, aborting."
    Call TerminateProgram
   End If

   Read(93,*) Mi_0(1)
   Read(93,*) M2_E_0(1,1),M2_E_0(2,2),M2_E_0(3,3)
   Read(93,*) M2_L_0(1,1),M2_L_0(2,2),M2_L_0(3,3)
   Read(93,*) M2_D_0(1,1),M2_D_0(2,2),M2_D_0(3,3)
   Read(93,*) M2_Q_0(1,1),M2_Q_0(2,2),M2_Q_0(3,3)
   Read(93,*) M2_U_0(1,1),M2_U_0(2,2),M2_U_0(3,3)
   Read(93,*) M2_H_0(1),M2_H_0(2)
   Read(93,*) AoY_q_0(1,1), AoY_q_0(2,2), AoY_q_0(3,3)
   Read(93,*) AoY_u_0(1,1), AoY_u_0(2,2), AoY_u_0(3,3)
   Read(93,*) AoY_d_0(1,1), AoY_d_0(2,2), AoY_d_0(3,3)
   Read(93,*) AoY_l_0(1,1), AoY_l_0(2,2), AoY_l_0(3,3)
   Read(93,*) tanb
   Read(93,*) phase_mu
   Read(93,*) TwoLoopRGE

   If (Abs(phase_mu).Ne.1._dp) phase_mu = phase_mu/abs(phase_mu)

   Mi_0 = Mi_0(1)
  Else If (HighScaleModel.Eq."Sugra") Then
   !-----------------------------------------------------------
   ! needed for the setting of the correct boundary conditions
   !-----------------------------------------------------------
   check = SetHighScaleModel("SUGRA")

   If (.Not.check) Then
    Write(ErrCan,*) "Strange error in subroutine HighScaleInput."
    Write(ErrCan,*) "Unable to put name of high scale model, aborting."
    Call TerminateProgram
   End If

   Read(93,*) M0_ii_in
   Mi_0 = M0_ii_in
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_E_0(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_L_0(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_D_0(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_Q_0(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_U_0(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in(1:2)
   M2_H_0 = M0_ii_in(1:2)**2
   Read(93,*) M0_ii_in
   Do i1=1,3
    AoY_u_0(i1,i1) = M0_ii_in(i1)
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    AoY_d_0(i1,i1) = M0_ii_in(i1)
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    AoY_l_0(i1,i1) = M0_ii_in(i1)
   End Do
   Read(93,*) tanb
   Read(93,*) sign_mu
   Read(93,*) TwoLoopRGE

 
  Else If (HighScaleModel.Eq."MSSMtree") Then

   Mi = 0._dp
   M2_E = 0._dp
   M2_L = 0._dp
   M2_D = 0._dp
   M2_U = 0._dp
   M2_Q = 0._dp
   AoY_u = 0._dp
   AoY_d = 0._dp
   AoY_l = 0._dp

   Read(93,*) M0_ii_in
   Mi = M0_ii_in
   If (GenerationMixing) then
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_E(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_L(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_D(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_Q(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_U(i1,:) = M0_ii_in
    End Do

    Do i1=1,3
     Read(93,*) M0_ii_in
     A_u(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     A_d(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     A_l(i1,:) = M0_ii_in
    End Do

   Else 
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_E(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_L(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_D(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_Q(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_U(i1,i1) = M0_ii_in(i1)
    End Do

    Read(93,*) M0_ii_in
    Do i1=1,3
     A_u(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     A_d(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     A_l(i1,i1) = M0_ii_in(i1)
    End Do
   end if ! Generation Mixing
   Read(93,*) tanb
   Read(93,*) M0_ii_in(1:2)
   mu = M0_ii_in(1)
   sign_mu = Real(mu,dp)/Abs(mu)
   mP0(2) = M0_ii_in(2)
   mP02(2) = mP0(2)**2

   test = SetRenormalizationScale(mZ2)
   ! this parameter is needed for consistency
   B = mP02(2) * tanb / (1._dp + tanb**2)
   ! calculating Yukawa couplings + rescaling A-Parameters
   sinW2 = 1._dp - mW2/mZ2
   vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
   vevSM(1) = vev / Sqrt(1._dp + tanb**2)
   vevSM(2) = tanb * vevSM(1)
   cosb2 = 1._dp / (1._dp + tanb**2)
   sinb2 = tanb**2 * cosb2
   cos2b = cosb2 - sinb2
   M2_H(1) = ( cos2b*sinb2*( mP02(2) + mZ2)          &
              & - cos2b*(Abs(mu)**2 + 0.5_dp * mZ2) ) / (cosb2 - sinb2)
   M2_H(2) = ( -cos2b*cosb2*(mP02(2) + mZ2)           &
              & + cos2b*(Abs(mu)**2 + 0.5_dp * mZ2) ) / (-cosb2 + sinb2)

   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp
   Do i1=1,3
    y_l(i1,i1) = sqrt2 * mf_L(i1) / vevSM(1)
    y_d(i1,i1) = sqrt2 * mf_D(i1) / vevSM(1)
    y_u(i1,i1) = sqrt2 * mf_U(i1) / vevSM(2)
   End Do
   If (GenerationMixing) then
    YukawaScheme = GetYukawaScheme()
    If (YukawaScheme.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM),Y_d) 
    End If
   end if

  Else If (HighScaleModel.Eq."MSSM") Then

   Mi = 0._dp
   M2_E = 0._dp
   M2_L = 0._dp
   M2_D = 0._dp
   M2_U = 0._dp
   M2_Q = 0._dp
   AoY_u = 0._dp
   AoY_d = 0._dp
   AoY_l = 0._dp
   A_u = 0._dp
   A_d = 0._dp
   A_l = 0._dp

   Read(93,*) M0_ii_in
   Mi = M0_ii_in

   if (GenerationMixing) then
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_E(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_L(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_D(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_Q(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     M2_U(i1,:) = M0_ii_in
    End Do

    Do i1=1,3
     Read(93,*) M0_ii_in
     A_u(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     A_d(i1,:) = M0_ii_in
    End Do
    Do i1=1,3
     Read(93,*) M0_ii_in
     A_l(i1,:) = M0_ii_in
    End Do

   else  ! .not.GenerationMixing
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_E(i1,i1) = M0_ii_in(i1)**2
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_L(i1,i1) = M0_ii_in(i1)**2
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_D(i1,i1) = M0_ii_in(i1)**2
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_Q(i1,i1) = M0_ii_in(i1)**2
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     M2_U(i1,i1) = M0_ii_in(i1)**2
    End Do

    Read(93,*) M0_ii_in
    Do i1=1,3
     AoY_u(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     AoY_d(i1,i1) = M0_ii_in(i1)
    End Do
    Read(93,*) M0_ii_in
    Do i1=1,3
     AoY_l(i1,i1) = M0_ii_in(i1)
    End Do
   end if ! GenerationMixing
   Read(93,*) tanb,Scale_Q

   Read(93,*) M0_ii_in(1:2)

   mu = M0_ii_in(1)
   sign_mu = Real(mu,dp)/Abs(mu)
   mP0(2) = M0_ii_in(2)   ! tree level pseudoscalar mass
   mP02(2) = mP0(2)**2
   B = mP02(2) * tanb / (1._dp + tanb**2)

   test = SetRenormalizationScale(Scale_Q**2)
   ! calculating Yukawa couplings + rescaling A-Parameters
   sinW2 = 1._dp - mW2/mZ2
   vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
   vevSM(1) = vev / Sqrt(1._dp + tanb**2)
   vevSM(2) = tanb * vevSM(1)

   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp
   If (GenerationMixing) then
    YukawaScheme = GetYukawaScheme()
    If (YukawaScheme.Eq.1) Then
     Y_u = Matmul(Transpose(CKM),Y_u) 
    Else
     Y_d = Matmul(Conjg(CKM),Y_d) 
    End If

   else

    Do i1=1,3
     y_l(i1,i1) = sqrt2 * mf_L(i1) / vevSM(1)
     A_l(i1,i1) = AoY_l(i1,i1) * y_l(i1,i1)
     y_d(i1,i1) = sqrt2 * mf_D(i1) / vevSM(1)
     A_d(i1,i1) = AoY_d(i1,i1) * y_d(i1,i1)
     y_u(i1,i1) = sqrt2 * mf_U(i1) / vevSM(2)
     A_u(i1,i1) = AoY_u(i1,i1) * y_u(i1,i1)
    End Do
   end if

  Else If (HighScaleModel.Eq."pMSSM") Then

   Mi = 0._dp
   M2_E = 0._dp
   M2_L = 0._dp
   M2_D = 0._dp
   M2_U = 0._dp
   M2_Q = 0._dp
   AoY_u = 0._dp
   AoY_d = 0._dp
   AoY_l = 0._dp

   Read(93,*) M0_ii_in
   Mi = M0_ii_in

   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_E(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_L(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_D(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_Q(i1,i1) = M0_ii_in(i1)**2
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    M2_U(i1,i1) = M0_ii_in(i1)**2
   End Do

   Read(93,*) M0_ii_in
   Do i1=1,3
    AoY_u(i1,i1) = M0_ii_in(i1)
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    AoY_d(i1,i1) = M0_ii_in(i1)
   End Do
   Read(93,*) M0_ii_in
   Do i1=1,3
    AoY_l(i1,i1) = M0_ii_in(i1)
   End Do
   Read(93,*) tanb,Scale_Q
   Read(93,*) M0_ii_in(1:2)
   mu = M0_ii_in(1)
   sign_mu = Real(mu,dp)/Abs(mu)
   mP0(2) = M0_ii_in(2)
   mP02(2) = mP0(2)**2

   test = SetRenormalizationScale(Scale_Q**2)
    ! this parameter is needed for consistency
    B = mP02(2) * tanb / (1._dp + tanb**2)
    ! calculating Yukawa couplings + rescaling A-Parameters
   vev =  Sqrt( mZ2 * (1._dp - sinW2) * SinW2 / (pi * alpha_mZ) )
   vevSM(1) = vev / Sqrt(1._dp + tanb**2)
   vevSM(2) = tanb * vevSM(1)

   Y_l = 0._dp
   Y_d = 0._dp
   Y_u = 0._dp
   Do i1=1,3
    y_l(i1,i1) = sqrt2 * mf_L(i1) / vevSM(1)
    A_l(i1,i1) = AoY_l(i1,i1) * y_l(i1,i1)
    y_d(i1,i1) = sqrt2 * mf_D(i1) / vevSM(1)
    A_d(i1,i1) = AoY_d(i1,i1) * y_d(i1,i1)
    y_u(i1,i1) = sqrt2 * mf_U(i1) / vevSM(2)
    A_u(i1,i1) = AoY_u(i1,i1) * y_u(i1,i1)
   End Do


  Else
   Write(ErrCan,*) "Error in subroutine HighScaleInput. Model "//HighScaleModel
   Write(ErrCan,*) "is not defined, aborting."
   Write(*,*) "Error in subroutine HighScaleInput. Model "//HighScaleModel
   Write(*,*) "is not defined, aborting."
   Call TerminateProgram
  End If

  if (HighScaleModel.Ne."Oscar") phase_mu = sign_mu

  Close(93)

  Iname = Iname - 1

 End Subroutine HighScaleInput 


 Subroutine LesHouches_Input(kont, HighScaleModel, Ecms, Pm, Pp, l_ISR, Fgmsb)
 !--------------------------------------------------------------------
 ! reads in data using the Les Houches standard as defined in 
 ! hep-ph/03111231
 ! output:
 !   - kont .................. is 0 if the input is consistent, non-zero if
 !                             there is a problem
 !   - HighScaleModel ........ string specifiying the model
 !   - Ecms .................. center of mass energy in GeV
 !   - Pm .................... degree of polarisation of incoming electrons
 !   - Pp .................... degree of polarisation of incoming positrons
 !   - l_ISR ................. if .true. then calculate initial state rediation
 !                             stemming from the incoming elctron/positron beams
 !   - Fgmsb ................. the vev of the F-component in the GMSB model
 !--------------------------------------------------------------------
 Implicit None
  Integer, Intent(out) :: kont
  Real(dp), Intent(out) :: Fgmsb, Ecms(:), Pm(:), Pp(:)
  Character(len=15), Intent(out) :: HighScaleModel
  Logical, Intent(out) :: l_ISR(:)
  
  Character(len=80) :: read_line
  Integer :: i_mod=-1, i_sm=-1, i_par=-1, set_mod_par(27)=-1 &
    & , i1, p_max, p_act, i_sp, i_model=-1, i_particles=-1, i_rp=0
  Real(dp) :: wert, Abs_Mu2, cosb2, cos2b, sinb2, RG0(3,3) &
    & , mat_D(3,3), R2(2,2), s12, s13, s23, c12, c13, c23
  Logical :: check, calc_ferm, check_alpha(2), test_l
  Complex(dp) :: lam_vS
  Logical, Save :: l_open = .True.

  Iname = Iname + 1
  NameOfUnit(Iname) = "LesHouches_Input"

  check_alpha = .False. ! used to check consistency of alpha(mZ) calculation
  in_kont = 0

  Call InitializeStandardModel
  Call InitializeLoopFunctions
  !-------------------------------------------
  ! this has been shifted in case of a loop
  !-------------------------------------------
  i_mod = -1
  i_sm = -1
  i_par = -1
  set_mod_par = -1

  ErrorLevel = -1
  GenerationMixing=.False.
  L_BR=.True.
  L_CS=.False.
  L_ISR = .False. 
  If (l_open) Then
   Open(ErrCan,file="Messages.out",status="unknown")
   Open(11,file="SPheno.out",status="unknown")
   l_open = .False.
  End If

  !-------------------------------------------------------
  ! set all model parameters to zero
  !-------------------------------------------------------
  Call Set_All_Parameters_0()
  lam_vs = 0._dp
  sp_info = " "
  !-------------------------------------------------------
  ! necessary to exclude right handed neutrinos from RGEs
  ! is set positive in the corresponding model
  !-------------------------------------------------------
  MNuR = - 1.e-9_dp
  !-------------------------------------------------------
  !take highest precision, will be changed at a later stage
  !-------------------------------------------------------
  TwoLoopRGE = .True.
  !-----------------------------------------------------------------------
  ! these variables are only used in GMSB and will be set correctly below
  !-----------------------------------------------------------------------
  Fgmsb = 1.e12_dp
  m32 = 1.e20_dp ! set an abitrary large gravitino mass in eV

  kont = 0

  Open(99,file="LesHouches.in",status="old",err=200)

  Do ! reading file
   Read(99,"(a80)",End=200,err=200) read_line
! Write(*,*) trim(read_line)
   If (read_line(1:1).Eq."#") Cycle ! ignore comments for the moment
   If (read_line.Eq." ") Cycle ! ignore empty lines for the moment

   Call PutUpperCase(read_line)
   If (read_line(1:5).Eq."BLOCK") Then ! assigning values for the select case
    If (read_line(7:12).Eq."MODSEL") Then
     Call Read_MODSEL(99, i_particles, i_model, i_cpv, i_rp, kont)
     If (i_cpv.Eq.0) Then ! one has to recalculated the CKM to the real case
                                     ! because InitializeStandardModel assumes a non-zero phase
     s12 = lam_wolf
     s23 = s12**2 * A_wolf
     s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2) 
     c12 = Sqrt(1._dp-s12*s12)
     c23 = Sqrt(1._dp-s23*s23)
     c13 = Sqrt(1._dp-s13*s13)

     CKM(1,1) = c12 * c13
     CKM(1,2) = s12 * c13
     CKM(1,3) = s13      
     CKM(2,1) = -s12*c23 -c12*s23*s13 
     CKM(2,2) = c12*c23 -s12*s23*s13 
     CKM(2,3) = s23 * c13
     CKM(3,1) = s12*s23 -c12*c23*s13 
     CKM(3,2) = -c12*s23 - s12*c23*s13 
     CKM(3,3) = c23 * c13
    End If

    Else If (read_line(7:14).Eq."SMINPUTS") Then
     Call Read_SMinput(99)

    Else If (read_line(7:12).Eq."MINPAR") Then
     Call Read_MINPAR(99, 0, i_model, set_mod_par, kont)

    Else If (read_line(7:14).Eq."IMMINPAR") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMINPAR") 
      Cycle
     End If
     Call Read_MINPAR(99, 1, i_model, set_mod_par, kont)

    Else If (read_line(7:12).Eq."EXTPAR") Then
     Call Read_EXTPAR(99, 0, i_model, set_mod_par, kont)

    Else If (read_line(7:14).Eq."IMEXTPAR") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMEXTPAR") 
      Cycle
     End If
     Call Read_EXTPAR(99, 1, i_model, set_mod_par, kont)

    Else If (read_line(7:17).Eq."SPHENOINPUT") Then
     Call Read_SPhenoInput(99)

    Else If (read_line(7:12).Eq."SPINFO") Then
     Call  Read_SPINFO(99, kont)

    Elseif ((read_line(7:12).Eq."HIGMIX").Or.(read_line(7:12).Eq."NMHMIX")) Then
     Call ReadMatrixR(99, 3, RS03_save, "RS03_save", kont)

    Else If ((read_line(7:10).Eq."AMIX").Or.(read_line(7:12).Eq."NMAMIX")) Then
     Call ReadMatrixR(99, 3, RP03_save, "RP03_save", kont)

    Else If (read_line(7:15).Eq."RVKAPPAIN") Then
     Call ReadVectorC(99, 3, eps, 0, "Re(epsilon)", kont)

    Else If (read_line(7:17).Eq."IMRVKAPPAIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMRVKAPPAIN") 
      Cycle
     End If
     Call ReadVectorC(99, 3, eps, 1, "Im(epsilon)", kont)
     set_mod_par(26) = 1

    Else If (read_line(7:15).Eq."RVSNVEVIN") Then
     Call ReadVectorR(99, 3, vevL, "v_L", kont)  
     set_mod_par(27) = 1

!     Else If (read_line(7:17).Eq."RVLAMPBDAIN") Then
!      Call ReadTensorC(99, 3, Rp_lam, 0, "Re(lambda_ijk)", kont)
! 
!     Else If (read_line(7:19).Eq."IMRVLAMPBDAIN") Then
!      if (i_cpv.lt.2) then
!       call Warn_CPV(i_cpv, "IMRVLAMPBDAIN") 
!       cycle
!      end if
!      Call ReadTensorC(99, 3, Rp_lam, 1, "Im(lambda_ijk)", kont)
! 
!     Else If (read_line(7:18).Eq."RVLAMPBDAPIN") Then
!      Call ReadTensorC(99, 3, Rp_lamp, 0, "Re(lambda'_ijk)", kont)
! 
!     Else If (read_line(7:20).Eq."IMRVLAMPBDAPIN") Then
!      if (i_cpv.lt.2) then
!       call Warn_CPV(i_cpv, "IMRVLAMPBDAPIN") 
!       cycle
!      end if
!      Call ReadTensorC(99, 3, Rp_lamp, 1, "Im(lambda'_ijk)", kont)

    Else If (read_line(7:12).Eq."VCKMIN") Then
     Call Read_CKM(99,i_cpv)

    Else If (read_line(7:10).Eq."MSL2") Then
     Call ReadMatrixC(99, 3, M2_L, 0, "Re(M2_L)", kont, 1)
     M2_L_0 = M2_L
     l_ML = .True.
     set_mod_par(11:13) = 1

    Else If (read_line(7:12).Eq."IMMSL2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSL2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2_L, 1, "Im(M2_L)", kont, 1)
     M2_L_0 = M2_L
     l_ML = .True.
     set_mod_par(11:13) = 1

    Else If (read_line(7:10).Eq."MSE2") Then
     Call ReadMatrixC(99, 3, M2_E, 0, "Re(M2_E)", kont, 1)
     M2_E_0 = M2_E
     l_ME = .True.
     set_mod_par(14:16) = 1

    Else If (read_line(7:12).Eq."IMMSE2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSE2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2_E, 1, "Im(M2_E)", kont, 1)
     M2_E_0 = M2_E
     l_ME = .True.
     set_mod_par(14:16) = 1

    Else If (read_line(7:10).Eq."MSQ2") Then
     Call ReadMatrixC(99, 3, M2_Q, 0, "Re(M2_Q)", kont, 1)
     M2_Q_0 = M2_Q
     l_MQ = .True.
     set_mod_par(17:19) = 1

    Else If (read_line(7:12).Eq."IMMSQ2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSQ2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2_Q, 1, "Im(M2_Q)", kont, 1)
     M2_Q_0 = M2_Q
     l_MQ = .True.
     set_mod_par(17:19) = 1

    Else If (read_line(7:10).Eq."MSU2") Then
     Call ReadMatrixC(99, 3, M2_U, 0, "Re(M2_U)", kont, 1)
     M2_U_0 = M2_U
     l_MU = .True.
     set_mod_par(20:22) = 1

    Else If (read_line(7:12).Eq."IMMSU2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSU2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2_U, 1, "Im(M2_U)", kont, 1)
     M2_U_0 = M2_U
     l_MU = .True.
     set_mod_par(20:22) = 1

    Else If (read_line(7:10).Eq."MSD2") Then
     Call ReadMatrixC(99, 3, M2_D, 0, "Re(M2_D)", kont, 1)
     M2_D_0 = M2_D
     l_MD = .True.
     set_mod_par(23:25) = 1

    Else If (read_line(7:12).Eq."IMMSD2") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMMSD2") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, M2_D, 1, "Im(M2_D)", kont, 1)
     M2_D_0 = M2_D
     l_MD = .True.
     set_mod_par(23:25) = 1

    Else If (read_line(7:10).Eq."TUIN") Then
     Call ReadMatrixC(99, 3, A_U_0, 0, "Re(T_U)", kont)
     A_U = A_U_0
     l_Au = .True.
    Else If (read_line(7:12).Eq."IMTUIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMTU") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, A_U_0, 1, "Im(T_U)", kont)
     A_U = A_U_0
     l_Au = .True.

    Else If (read_line(7:10).Eq."TDIN") Then
     Call ReadMatrixC(99, 3, A_D_0, 0, "Re(T_D)", kont)
     A_D = A_D_0
     l_Ad = .True.
    Else If (read_line(7:12).Eq."IMTDIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMTD") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, A_D_0, 1, "Im(T_D)", kont)
     A_D = A_D_0
     l_Ad = .True.

    Else If (read_line(7:10).Eq."TEIN") Then
     Call ReadMatrixC(99, 3, A_l_0, 0, "Re(T_E)", kont)
     A_l = A_l_0
     l_Al = .True.
    Else If (read_line(7:12).Eq."IMTEIN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMTE") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, A_L_0, 1, "Im(T_E)", kont)
     A_l = A_l_0
     l_Al = .True.

    Else If (read_line(7:10).Eq."YNU0")Then
     Fixed_Nu_Yukawas = .True.
     Call ReadMatrixC(99, 3, Y_nu_0, 0, "Re(Y_nu_0)", kont)

    Else If (read_line(7:12).Eq."IMYNU0") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMYNU0") 
      Cycle
     End If
     Fixed_Nu_Yukawas = .True.
     Call ReadMatrixC(99, 3, Y_nu_0, 1, "Im(Y_nu_0)", kont)

    Else If (read_line(7:10).Eq."MNUR") Then
     Call ReadVectorR(99, 3, MnuR, "M_nu_R", kont)

    Else If (read_line(7:13).Eq."VPMNSIN") Then
     Call Read_PMNS(99)

    Else If (read_line(7:9).Eq."YT0") Then
     Call ReadMatrixC(99, 3, Y_T_0, 0, "Re(Y_T_0)", kont, 1)

    Else If (read_line(7:11).Eq."IMYT0") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMYT0") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, Y_T_0, 1, "Im(Y_T_0)", kont, 2)

! Florian Staub
    Else If (read_line(7:11).Eq."YB3IN") Then

     Call ReadMatrixC(99, 3, Yb3_H24_gut,0, "Yb3", kont)
     Yb30_H24(3,:,:) = Yb3_H24_gut
     Yb30_H24(2,:,:) = Yb30_H24(3,:,:)
     Yb30_H24(2,:,3) = 0._dp
     Yb30_H24(1,:,:) = Yb30_H24(2,:,:)
     Yb30_H24(1,:,2) = 0._dp
     Yw30_H24 = Yb30_H24
     Yx30_H24 = Yb30_H24


    Else If (read_line(7:10).Eq."YTIN") Then

     Call ReadMatrixC(99, 3, YT_H15_gut,0, "Yt", kont)
     YT_H15_GUT(2,1) = YT_H15_GUT(1,2)
     YT_H15_GUT(3,1) = YT_H15_GUT(1,3)
     YT_H15_GUT(3,2) = YT_H15_GUT(2,3)
     YT0_H15 = YT_H15_GUT
     YZ0_H15 = YT_H15_GUT
     YS0_H15 = YT_H15_GUT
     Y_T_0 = YT_H15_gut

    Else If (read_line(7:11).Eq."MWMIN") Then
     Call ReadMatrixC(99, 3, MWM30,0, "MWM3", kont)
     MGM3 = MWM30
     MWM3 = MWM30
     MBM3 = MWM30
     MXM3 = MWM30

    Else If (read_line(7:13).Eq."IMYB3IN") Then
     If (i_cpv.Lt.2) Then
      Call Warn_CPV(i_cpv, "IMYB3IN") 
      Cycle
     End If
     Call ReadMatrixC(99, 3, Yb30_H24, 1, "Im(Yb3)", kont)
! Florian Staub

    Else If (read_line(7:12).Eq."HIGGS3") Then
     Call Read_Higgs3(99)

    Else If (read_line(7:10).Eq."MASS") Then
     Call Read_MASS(99)

    Else If (read_line(7:10).Eq."UMIX") Then
     Call ReadMatrixC(99, 2, U, 0, "Re(U)", kont)

    Else If (read_line(7:12).Eq."IMUMIX") Then
     Call ReadMatrixC(99, 2, U, 1, "Im(U)", kont)

    Else If (read_line(7:10).Eq."VMIX") Then
     Call ReadMatrixC(99, 2, V, 0, "Re(V)", kont)

    Else If (read_line(7:12).Eq."IMVMIX") Then
     Call ReadMatrixC(99, 2, V, 1, "Im(V)", kont)

    Else If ((read_line(7:10).Eq."NMIX").Or. (read_line(7:12).Eq."NMNMIX")) Then
     If (HighScaleModel.Eq."NMSSM") Then
      Call ReadMatrixC(99, 5, N5, 0, "Re(N)", kont)
     Else
      Call ReadMatrixC(99, 4, N, 0, "Re(N)", kont)
     End If

    Else If (read_line(7:12).Eq."IMNMIX") Then
     If (HighScaleModel.Eq."NMSSM") Then
      Call ReadMatrixC(99, 5, N5, 1, "Im(N)", kont)
     Else
      Call ReadMatrixC(99, 4, N, 1, "Im(N)", kont)
     End If

    Else If (read_line(7:13).Eq."STOPMIX") Then
     Call ReadMatrixR(99, 2, R2, "Re(R_~t)", kont)
     If (Rsup(1,1).eq.0._dp) Rsup = Id6C             ! first initialization
     RSup(5:6,5:6) = Cmplx(R2, Aimag(RSup(5:6,5:6)), dp)

    Else If (read_line(7:15).Eq."IMSTOPMIX") Then
     Call ReadMatrixR(99, 2, R2, "Im(R_~t)", kont)
     If (Rsup(1,1).eq.0._dp) Rsup = Id6C             ! first initialization
     RSup(5:6,5:6) = Cmplx(Real(RSup(5:6,5:6),dp), R2, dp)

    Else If (read_line(7:13).Eq."SBOTMIX") Then
     Call ReadMatrixR(99, 2, R2, "Re(R_~b)", kont)
     If (Rsdown(1,1).eq.0._dp) Rsdown = Id6C
     RSdown(5:6,5:6) = Cmplx(R2, Aimag(RSdown(5:6,5:6)), dp)

    Else If (read_line(7:15).Eq."IMSBOTMIX") Then
     Call ReadMatrixR(99, 2, R2, "Im(R_~b)", kont)
     If (Rsdown(1,1).eq.0._dp) Rsdown = Id6C
     RSdown(5:6,5:6) = Cmplx(Real(RSdown(5:6,5:6),dp), R2, dp)

    Else If (read_line(7:13).Eq."STAUMIX") Then
     Call ReadMatrixR(99, 2, R2, "Re(R_~tau)", kont)
     If (Rslepton(1,1).eq.0._dp) Rslepton = Id6C
     RSlepton(5:6,5:6) = Cmplx(R2, Aimag(RSlepton(5:6,5:6)), dp)

    Else If (read_line(7:15).Eq."IMSTAUMIX") Then
     Call ReadMatrixR(99, 2, R2, "Im(R_~tau)", kont)
     If (Rslepton(1,1).eq.0._dp) Rslepton = Id6C
     RSlepton(5:6,5:6) = Cmplx(Real(RSlepton(5:6,5:6),dp), R2, dp)

    Else If (read_line(7:24).Eq."NEUTRINOBOUNDSIN") Then
     Call Read_Neutrino_Bounds(99)

    Else If (read_line(7:19).Eq."STARTDATAFILE") Then
     Read(99,*) Old_Data
     Old_data= Trim(Old_data) ! to avoid trailing blanks

    Else
     If (output_screen) Write(*,*) "Warning, the following block is ignored"
     If (output_screen) Write(*,*) Trim(read_line)
     Write(ErrCan,*) "Warning, the following block is ignored"
     Write(ErrCan,*) Trim(read_line)

    End If
   End If

   If (kont.Ne.0) Return

  End Do 
200 Close(99)

  !---------------------------------------------
  ! saving parameters in the super CKm basis
  !---------------------------------------------
  If (l_Ad) Then
   A_D_0 = Transpose(A_D_0)   ! SLHA convention contain a transpose
   A_D = A_D_0
   Ad_sckm = A_d_0
  End If
  If (l_Au) Then
   A_U_0 = Transpose(A_U_0)   ! SLHA convention contain a transpose
   A_U = A_U_0
   Au_sckm = A_u_0
  End If
  If (l_MD) M2D_sckm = M2_D_0
  If (l_MQ) M2D_sckm = M2_Q_0
  If (l_MU) M2D_sckm = M2_U_0
!-----------------------------------------------
! now some checks and additional settings
!-----------------------------------------------
  If (i_particles.Eq.1) Then  ! MSSM particle content
   If (i_model.Eq.0) Then 
    If ((set_mod_par(7).Eq.1).And.(set_mod_par(8).Eq.1)) Then
     If (i_rp.Eq.0) HighScaleModel = "MSSM1"
     If (set_mod_par(9).Eq.1) Then
      Write(ErrCan,*)  "m^2_H1 and m^2_H2 have been specified together with mu"
      Write(ErrCan,*)  "mu will be ignored"
      set_mod_par(9) = 0
      in_kont(2) = 1
     End If
     If (set_mod_par(9).Eq.-1) set_mod_par(9) = 0 ! to avoid problems with the sum
     If (set_mod_par(10).Eq.1) Then
      Write(ErrCan,*) &
           &  "m^2_H1 and m^2_H2 have been specified together with mA^2"
      Write(ErrCan,*)  "mA^2 will be ignored"
      in_kont(2) = 1
      set_mod_par(10) = 0
     End If
     If (set_mod_par(10).Eq.-1) set_mod_par(10) = 0 ! to avoid problems with the sum
     !-------------------------------
     ! first guess of mu and B, mA
     !-------------------------------
     If (i_rp.Eq.0) Then
      cosb2 = 1._dp / (1._dp + tanb**2)
      sinb2 = tanb**2 * cosb2
      cos2b = cosb2 - sinb2
      Abs_Mu2 = (M2_H(2) * sinb2 - M2_H(1) * cosb2 )/ cos2b - 0.5_dp * mZ2
      If (Abs_mu2.Lt.0._dp) Abs_mu2 = 1.e4_dp
      mu = Sqrt(abs_mu2) * phase_mu
      B = (M2_H(1) + M2_H(2) + 2._dp *  Abs_Mu2) * tanb / (1+tanb**2)
      mP02(2) = Abs(B) * (1._dp/tanb + tanb)
      mP0(2) = Sqrt(mp02(2))
     Else
     End If
    Else If ((set_mod_par(9).Eq.1).And.(set_mod_par(10).Eq.1)) Then
!     HighScaleModel = "MSSM"
     If (set_mod_par(7).Eq.1) Then
      Write(ErrCan,*)  "mu and m_A0 have been specified together with M^2_Hd"
      Write(ErrCan,*)  "M^2_Hd will be ignored"
      set_mod_par(7) = 0
      in_kont(2) = 1
     End If
     If (set_mod_par(7).Eq.-1) set_mod_par(7) = 0 ! to avoid problems with the sum
     If (set_mod_par(8).Eq.1) Then
      Write(ErrCan,*)  "mu and m_A0 have been specified together with M^2_Hu"
      Write(ErrCan,*)  "M^2_Hu will be ignored"
      set_mod_par(8) = 0
      in_kont(2) = 1
     End If
     If (set_mod_par(8).Eq.-1) set_mod_par(8) = 0 ! to avoid problems with the sum
     !-------------------------------
     ! first guess of B, M^2_H_i
     !-------------------------------
     mA2_Q =  mP02(2)
     B = mP02(2) * tanb / (1._dp + tanb**2)

    Else 
     Write(ErrCan,*) "Higgs sector has not been specified consistently. Aborting"
     kont = -307
     Call AddError(307)
     Return
    End If

   End If ! i_mod.eq.0

     kont = -305 ! model has not specified completly
     If ((i_model.Eq.0).And.(Sum(set_mod_par(1:25)).Eq.23)) kont = 0 ! MSSM
     If ((i_model.Eq.1).And.(Sum(set_mod_par(1:5)).Eq.5)) kont = 0 ! mSugra 
     If ((i_model.Eq.2).And.(Sum(set_mod_par(1:5)).Eq.5)) Then ! GMSB 
      kont = 0
      AoY_d_0 = 0._dp
      AoY_l_0 = 0._dp
      AoY_nu_0 = 0._dp
      AoY_u_0 = 0._dp
      Fgmsb = Lambda * MlambdaS ! needed for the calculation of NLSP
      If (grav_fac.Ge.0) Then
       m32 = grav_fac * 2.4e-10_dp * Fgmsb  ! gravitino mass in eV  
      Else
       m32 = 2.4e-10_dp * Fgmsb  ! gravitino mass in eV  
      End If

     End If
     If ((i_model.Eq.3).And.(Sum(set_mod_par(1:4)).Eq.4)) kont = 0 ! AMSB
     If (kont.Ne.0) Call AddError(Abs(kont))

  Else If (i_particles.Eq.2) Then  ! MSSM + nu_R particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam

   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2_L(1,1) = M2_L(3,3)
   If (set_mod_par(12).Ne.1) M2_L(2,2) = M2_L(3,3)
   If (set_mod_par(14).Ne.1) M2_E(1,1) = M2_E(3,3)
   If (set_mod_par(15).Ne.1) M2_E(2,2) = M2_E(3,3)
   If (set_mod_par(17).Ne.1) M2_Q(1,1) = M2_Q(3,3)
   If (set_mod_par(18).Ne.1) M2_Q(2,2) = M2_Q(3,3)
   If (set_mod_par(20).Ne.1) M2_U(1,1) = M2_U(3,3)
   If (set_mod_par(21).Ne.1) M2_U(2,2) = M2_U(3,3)
   If (set_mod_par(23).Ne.1) M2_D(1,1) = M2_D(3,3)
   If (set_mod_par(24).Ne.1) M2_D(2,2) = M2_D(3,3)
   !----------------------------------------------------
   ! check if Higgs masses and mixings are given
   !----------------------------------------------------
   A_nu(1,:) = AoY_nu(1,:) * Y_nu(1,:) 

  Else If (i_particles.Eq.5) Then  ! MSSM + 2 nu_R particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam

   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2_L(1,1) = M2_L(3,3)
   If (set_mod_par(12).Ne.1) M2_L(2,2) = M2_L(3,3)
   If (set_mod_par(14).Ne.1) M2_E(1,1) = M2_E(3,3)
   If (set_mod_par(15).Ne.1) M2_E(2,2) = M2_E(3,3)
   If (set_mod_par(17).Ne.1) M2_Q(1,1) = M2_Q(3,3)
   If (set_mod_par(18).Ne.1) M2_Q(2,2) = M2_Q(3,3)
   If (set_mod_par(20).Ne.1) M2_U(1,1) = M2_U(3,3)
   If (set_mod_par(21).Ne.1) M2_U(2,2) = M2_U(3,3)
   If (set_mod_par(23).Ne.1) M2_D(1,1) = M2_D(3,3)
   If (set_mod_par(24).Ne.1) M2_D(2,2) = M2_D(3,3)
   !----------------------------------------------------
   ! check if Higgs masses and mixings are given
   !----------------------------------------------------
   A_nu(1:2,:) = AoY_nu(1:2,:) * Y_nu(1:2,:) 
   A_h02 = Ao_h02 * h02
   A_lam222 = Ao_lam222 * lam2
   A_lam112 = Ao_lam112 * lam112
   A_lam122 = Ao_lam122 * lam122

  Else If (i_particles.Eq.3) Then  ! NMSSM particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam

   !--------------------------------------
   ! the input is term of an effective mu
   !--------------------------------------
   vP = sqrt2 * lam_vS / h0
   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2_L(1,1) = M2_L(3,3)
   If (set_mod_par(12).Ne.1) M2_L(2,2) = M2_L(3,3)
   If (set_mod_par(14).Ne.1) M2_E(1,1) = M2_E(3,3)
   If (set_mod_par(15).Ne.1) M2_E(2,2) = M2_E(3,3)
   If (set_mod_par(17).Ne.1) M2_Q(1,1) = M2_Q(3,3)
   If (set_mod_par(18).Ne.1) M2_Q(2,2) = M2_Q(3,3)
   If (set_mod_par(20).Ne.1) M2_U(1,1) = M2_U(3,3)
   If (set_mod_par(21).Ne.1) M2_U(2,2) = M2_U(3,3)
   If (set_mod_par(23).Ne.1) M2_D(1,1) = M2_D(3,3)
   If (set_mod_par(24).Ne.1) M2_D(2,2) = M2_D(3,3)
   !----------------------------------------------------
   ! check if Higgs masses and mixings are given
   !----------------------------------------------------
   If (      (Maxval(mSpm).Gt.0._dp).And.(Maxval(mP03).Gt.0._dp)            &
      & .And.(Maxval(mS03).Gt.0._dp).And.(Maxval(Abs(RP03_save)).Gt.0._dp) &
      & .And.(Maxval(Abs(RS03_save)).Gt.0._dp) ) Then
!    External_Higgs = .True.
    mP03(1) = mZ
    mP032 = mP03**2
    mS032 = mS03**2
    mSpm(1) = mW
    mSpm2 = mSpm**2

    RS03 = RS03_save
    RG0 = 0._dp
    RG0(1,1) = - 1._dp / Sqrt(1._dp + tanb**2)
    RG0(2,2) = RG0(1,1)
    RG0(3,3) = 1._dp
    RG0(1,2) = tanb * RG0(1,1)
    RG0(2,1) = RG0(1,2)
!    RP03 = Matmul(RP03_save, RG0)
    RP03 = RP03_save
    RSpm = RG0(1:2,1:2)
   End If

  Else If (i_particles.Eq.4) Then  ! NMSSM + nu_R + S  particle content
   Call SetRGEscale(mZ2)
   tanb_mZ = tanb
   A_h0 = Ao_h0 * h0
   A_lam = Ao_lam * lam
   A_nu(1,:) = AoY_nu(1,:) * Y_nu(1,:) 
   A_pns = Ao_hpns * h_pns
   If (lam_vs.Ne.0._dp) Then ! calulate either h0 or v_P
    If ((h0.Ne.0._dp).And.(vP.Eq.0._dp)) vp = sqrt2 * lam_vs / h0
    If ((h0.Eq.0._dp).And.(vP.Ne.0._dp)) h0 = sqrt2 * lam_vs / vp
   End If 
   !----------------------------------------------------
   ! check if 1st and 2nd generation parameters are set
   !----------------------------------------------------
   If (set_mod_par(11).Ne.1) M2_L(1,1) = M2_L(3,3)
   If (set_mod_par(12).Ne.1) M2_L(2,2) = M2_L(3,3)
   If (set_mod_par(14).Ne.1) M2_E(1,1) = M2_E(3,3)
   If (set_mod_par(15).Ne.1) M2_E(2,2) = M2_E(3,3)
   If (set_mod_par(17).Ne.1) M2_Q(1,1) = M2_Q(3,3)
   If (set_mod_par(18).Ne.1) M2_Q(2,2) = M2_Q(3,3)
   If (set_mod_par(20).Ne.1) M2_U(1,1) = M2_U(3,3)
   If (set_mod_par(21).Ne.1) M2_U(2,2) = M2_U(3,3)
   If (set_mod_par(23).Ne.1) M2_D(1,1) = M2_D(3,3)
   If (set_mod_par(24).Ne.1) M2_D(2,2) = M2_D(3,3)

  End If
  !---------------------------------------------
  ! warning if alpha(mZ) and alpha(0) are given
  !---------------------------------------------
  If (check_alpha(1).And.check_alpha(2)) Then
   Write(ErrCan,*) "Warning: alpha(0) and alpha(mZ) have been specified,"
   Write(ErrCan,*) "the consisteny has not been checked!"
   in_kont(1) = 1 
  End If
  !------------------------------------------
  ! recalculate quantities to be sure
  !------------------------------------------
  gmZ = gamZ * mZ
  gmZ2 = gmZ**2

 !---------
 ! W-boson, first rough estimate
 !---------
  mW2 = mZ2 * (0.5_dp + Sqrt(0.25_dp-Alpha_Mz*pi / (sqrt2*G_F*mZ2))) / 0.985_dp

  mW = Sqrt(mW2)      ! mass
  gamW = 2.06_dp     ! width
  gamW2 = gamW**2
  gmW = gamW * mW
  gmW2 = gmW**2
  !------------------------------------------------------------------
  ! the running fermion masses at m_Z need to be recalculated
  !------------------------------------------------------------------
  Alpha_mZ = Alpha_MSbar(mZ, mW)
  If (calc_ferm) Call CalculateRunningMasses(mf_l, mf_d, mf_u            &
                           &  , Q_light_quarks, alpha_mZ, alphas_mZ, mZ  &
                           &  , mf_l_mZ, mf_d_mZ, mf_u_mZ, kont)

  If (GenerationMixing) Then
   M2Q_sckm = M2_Q
   M2D_sckm = M2_D
   M2U_sckm = M2_U
  End If

  If (HighScaleModel.Eq."SUGRA_NuR") Then
   Y_nu_0 = Transpose(Y_nu_0) ! in the RGEs the transposed Yukawas are used
   Y_nu_mR(1,:,1) = Y_nu_0(:,1)
   Y_nu_mR(1,:,1:2) = Y_nu_0(:,1:2)
   Y_nu_mR(3,:,:) = Y_nu_0
  End If

  If ((HighScaleModel(1:9).Eq."SUGRA_NuR").And.(D_SO_10.Ne.0._dp)) Then
    mat_D = 0._dp
    mat_D(1,1) = D_SO_10
    mat_D(2,2) = mat_D(1,1)
    mat_D(3,3) = mat_D(1,1)
    M2_E_0 = M2_E_0 + mat_D
    M2_L_0 = M2_L_0 - 3._dp * mat_D
    M2_R_0 = M2_R_0 + 5._dp * mat_D
    M2_D_0 = M2_D_0 - 3._dp * mat_D
    M2_Q_0 = M2_Q_0 + mat_D
    M2_U_0 = M2_U_0 + mat_D
    M2_H_0(1) =  M2_H_0(1) + 2._dp * D_SO_10
    M2_H_0(2) =  M2_H_0(2) - 2._dp * D_SO_10
  End If

  If (External_spectrum) Then
   If (GenerationMixing) Then
    Write(*,*) "Sorry, but the use of an external spectrum with generationmixing"
    Write(*,*) "is not implemented"
    Stop 99
   Else
    RSneut = id3C
    phaseGlu = 1._dp
    !----------------------------------------------------
    ! check for the order of sfermions
    !----------------------------------------------------
    If (mSup(1).Gt.mSup(2)) Then
     wert = mSup(2)
     mSup(2) = mSup(1)
     mSup(1) = wert
     RSup(1,1) = 0._dp
     RSup(2,2) = 0._dp
     Rsup(1,2) = 1
     Rsup(2,1) = -1._dp
    End If
    If (mSup(3).Gt.mSup(4)) Then
     wert = mSup(4)
     mSup(4) = mSup(3)
     mSup(3) = wert
     RSup(3,3) = 0._dp
     RSup(4,4) = 0._dp
     Rsup(3,4) = 1
     Rsup(4,3) = -1._dp
    End If
    If (mSdown(1).Gt.mSdown(2)) Then
     wert = mSdown(2)
     mSdown(2) = mSdown(1)
     mSdown(1) = wert
     RSdown(1,1) = 0._dp
     RSdown(2,2) = 0._dp
     Rsdown(1,2) = 1
     Rsdown(2,1) = -1._dp
    End If
    If (mSdown(3).Gt.mSdown(4)) Then
     wert = mSdown(4)
     mSdown(4) = mSdown(3)
     mSdown(3) = wert
     RSdown(3,3) = 0._dp
     RSdown(4,4) = 0._dp
     Rsdown(3,4) = 1
     Rsdown(4,3) = -1._dp
    End If
    If (mSlepton(1).Gt.mSlepton(2)) Then
     wert = mSlepton(2)
     mSlepton(2) = mSlepton(1)
     mSlepton(1) = wert
     RSlepton(1,1) = 0._dp
     RSlepton(2,2) = 0._dp
     Rslepton(1,2) = 1
     Rslepton(2,1) = -1._dp
    End If
    If (mSlepton(3).Gt.mSlepton(4)) Then
     wert = mSlepton(4)
     mSlepton(4) = mSlepton(3)
     mSlepton(3) = wert
     RSlepton(3,3) = 0._dp
     RSlepton(4,4) = 0._dp
     Rslepton(3,4) = 1
     Rslepton(4,3) = -1._dp
    End If
   ! changing neutralino masses to positive values
    If (HighScaleModel.Eq."NMSSM") Then
     Do i1=1,5
      If (mN5(i1).Lt.0._dp) Then
       mN5(i1) = - mN5(i1)
       N5(i1,:) = (0._dp,1._dp) * N5(i1,:)
      End If
     End Do
    Else
     Do i1=1,4
      If (mN(i1).Lt.0._dp) Then
       mN(i1) = - mN(i1)
       N(i1,:) = (0._dp,1._dp) * N(i1,:)
      End If
     End Do
    End If
   End If
  End If

  Iname = Iname - 1

 Contains

  Subroutine Read_Neutrino_Bounds(io)
  Implicit None
   Integer, Intent(in) :: io

    Do 
     Read(io,*,End=200) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)
     Case(1)
      m2_atm_min = wert
     Case(2)
      m2_atm_max = wert
     Case(3)
      tan2_atm_min = wert
     Case(4)
      tan2_atm_max = wert
     Case(5)
      m2_sol_min = wert
     Case(6)
      m2_sol_max = wert
     Case(7)
      tan2_sol_min = wert
     Case(8)
      tan2_sol_max = wert
     Case(9)
      Ue32_min = wert
     Case(10)
      Ue32_max = wert
     Case default
      Write(ErrCan,*) "Reading block NeutrinoBoundsIn"
      Write(ErrCan,*) "Particle with id=",i1," is unknown"
      Write(ErrCan,*) "The assigned value is",wert
     End Select

    End Do

  200 Return

  End Subroutine Read_Neutrino_Bounds


  Subroutine Warn_CPV(i_cpv, name)
  Implicit None 
   Integer, Intent(in) :: i_cpv
   Character(len=*), Intent(in) :: name
   If (i_cpv.Eq.0) Write(ErrCan,*) "CP violation is switched off"
   If (i_cpv.Eq.1) Write(ErrCan,*) "CP violation beyond CKM is switched off"
   Write(ErrCan,*) "Ignoring block "//Trim(name)
   If (ErrorLevel.Eq.2) Call TerminateProgram
  End Subroutine Warn_CPV


  Subroutine Read_MASS(io)
  Implicit None
   Integer, Intent(in) :: io

    Do 
     Read(io,*,End=200) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)
     Case(25)
      mS0(1) = wert
      mS02(1) = mS0(1)**2
      mS03(1) = wert
      mS032(1) = mS02(1)
      mS05(1) = wert
      mS052(1) = mS02(1)
     Case(35)
      mS0(2) = wert
      mS02(2) = mS0(2)**2
      mS03(2) = wert
      mS032(2) = mS02(2)
      mS05(2) = wert
      mS052(2) = mS02(2)
     Case(36)
      mP0(2) = wert
      mP02(2) = mP0(2)**2
      mP03(2) = wert
      mP032(2) = mP02(2)
      mP05(2) = wert
      mP052(2) = mP02(2)
     Case(37)
      mSpm(2) = wert
      mSpm2(2) = mSpm(2)**2
     Case(45)
      mS03(3) = wert
      mS032(3) = wert**2
     Case(46)
      mP03(3) = wert
      mP032(3) = wert**2
     Case(1000001)
      mSdown(1) = wert
      mSdown2(1) = wert**2
     Case(1000002)
      mSup(1) = wert
      mSup2(1) = wert**2
     Case(1000003)
      mSdown(3) = wert
      mSdown2(3) = wert**2
     Case(1000004)
      mSup(3) = wert
      mSup2(3) = wert**2
     Case(1000005)
      mSdown(5) = wert
      mSdown2(5) = wert**2
     Case(1000006)
      mSup(5) = wert
      mSup2(5) = wert**2
     Case(1000011)
      mSlepton(1) = wert
      mSlepton2(1) = wert**2
     Case(1000012)
      mSneut(1) = wert
      mSneut2(1) = wert**2
     Case(1000013)
      mSlepton(3) = wert
      mSlepton2(3) = wert**2
     Case(1000014)
      mSneut(2) = wert
      mSneut2(2) = wert**2
     Case(1000015)
      mSlepton(5) = wert
      mSlepton2(5) = wert**2
     Case(1000016)
      mSneut(3) = wert
      mSneut2(3) = wert**2
     Case(1000022)
      mN(1) = wert
      mN2(1) = wert**2
      mN5(1) = wert
      mN52(1) = mN2(1)
     Case(1000023)
      mN(2) = wert
      mN2(2) = wert**2
      mN5(2) = wert
      mN52(2) = mN2(2)
     Case(1000024)
      mC(1) = wert
      mC2(1) = wert**2
     Case(1000025)
      mN(3) = wert
      mN2(3) = wert**2
      mN5(3) = wert
      mN52(3) = mN2(3)
     Case(1000035)
      mN(4) = wert
      mN2(4) = wert**2
      mN5(4) = wert
      mN52(4) = mN2(4)
     Case(1000037)
      mC(2) = wert
      mC2(2) = wert**2
     Case(1000045)
      mN5(5) = wert
      mN52(5) = wert**2
     Case(1000021)
      mGlu = wert
     Case(2000001)
      mSdown(2) = wert
      mSdown2(2) = wert**2
     Case(2000002)
      mSup(2) = wert
      mSup2(2) = wert**2
     Case(2000003)
      mSdown(4) = wert
      mSdown2(4) = wert**2
     Case(2000004)
      mSup(4) = wert
      mSup2(4) = wert**2
     Case(2000005)
      mSdown(6) = wert
      mSdown2(6) = wert**2
     Case(2000006)
      mSup(6) = wert
      mSup2(6) = wert**2
     Case(2000011)
      mSlepton(2) = wert
      mSlepton2(2) = wert**2
     Case(2000013)
      mSlepton(4) = wert
      mSlepton2(4) = wert**2
     Case(2000015)
      mSlepton(6) = wert
      mSlepton2(6) = wert**2
     Case default
      Write(ErrCan,*) "Particle with id=",i1," is unknown"
      Write(ErrCan,*) "The assigned mass is",wert," GeV"
     End Select

    End Do

  200 Return

  End Subroutine Read_MASS


  Subroutine Read_Higgs3(io)
  Implicit None
   Integer, Intent(in) :: io

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."#").Or.(read_line(1:1).Eq."B")  &
                                .Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)
     Case(1)
      M_H3 = wert
      ! as initalization, will be computed more precisely later
      MS15_mH3 = wert
      MT15_mH3 = wert
      MZ15_mH3 = wert
     Case(2)
      lam12_0(1) = Cmplx(wert, Aimag(lam12_0(1) ), dp)
     Case(3)
      lam12_0(1) = Cmplx(Real(lam12_0(1), dp ), wert, dp) 
     Case(4)
      lam12_0(2) = Cmplx(wert, Aimag(lam12_0(2) ), dp)
     Case(5)
      lam12_0(2) = Cmplx(Real(lam12_0(2), dp ), wert, dp) 
     Case(6)
      If (wert.Eq.1) Then
       Fifteen_plet = .True.
      Else
       Fifteen_plet = .False.
      End If
     End Select

    End Do

   200 Return

  End Subroutine Read_Higgs3

  Subroutine Read_PMNS(io)
  Implicit None
   Integer, Intent(in) :: io

   Real(dp) :: s12, s13, s23, c12, c13, c23, phase, alpha, beta

    s12 = 0._dp
    s13 = 0._dp
    s23 = 0._dp
    phase = 0._dp
    alpha = 0._dp
    beta = 0._dp

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."#").Or.(read_line(1:1).Eq."B")  &
                                .Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)     
     Case(1)
      theta_12 = wert
      s12 = Sin(wert)
     Case(2)
      theta_23 = wert
      s23 = Sin(wert)
     Case(3)
      theta_13 = wert
      s13 = Sin(wert)
     Case(4)
      delta_nu = wert
      phase = wert
     Case(5)
      alpha_nu1 = wert
      alpha = wert
     Case(6)
      alpha_nu2 = wert
      beta = wert
     Case default
     End Select

    End Do

 200 c12 = Sqrt(1._dp-s12*s12)
    c23 = Sqrt(1._dp-s23*s23)
    c13 = Sqrt(1._dp-s13*s13)

    Unu(1,1) = c12 * c13
    Unu(1,2) = s12 * c13
    Unu(2,3) = s23 * c13
    Unu(3,3) = c23 * c13
    If (phase.Ne.0._dp) Then
     Unu(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
     Unu(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     Unu(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     Unu(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
     Unu(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    Else
     Unu(1,3) = s13
     Unu(2,1) = -s12*c23 -c12*s23*s13
     Unu(2,2) = c12*c23 -s12*s23*s13
     Unu(3,1) = s12*s23 -c12*c23*s13
     Unu(3,2) = -c12*s23 - s12*c23*s13
    End If
    ! Majorana phases
    If (alpha.Ne.0._dp) Unu(1,:) =  Unu(1,:) * Exp( (0._dp,-0.5_dp) * alpha )
    If (beta.Ne.0._dp) Unu(2,:) =  Unu(2,:) * Exp( (0._dp,-0.5_dp) * beta )

  End Subroutine Read_PMNS


  Subroutine Read_CKM(io, i_cpv)
  Implicit None
   Integer, Intent(in) :: io, i_cpv

   Real(dp) :: s12, s13, s23, c12, c13, c23, phase
    
    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."#").Or.(read_line(1:1).Eq."B")  &
                                .Or.(read_line(1:1).Eq."b") ) Exit ! this loop
     Read(io,*) i1, wert, read_line
     Select Case(i1)     
     Case(1)
      lam_wolf = wert
     Case(2)
      A_wolf = wert
     Case(3)
      rho_wolf = wert
     Case(4)
      eta_wolf = wert
     Case default
     End Select

    End Do

 200   s12 = lam_wolf
    s23 = s12**2 * A_wolf
    s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2)
    If (i_cpv.Eq.0) Then
     Write(ErrCan,*) "Warning: CP violation is switched of, ignoring CKM phase."
     phase = 0._dp
    Else
     phase = Atan(eta_wolf/rho_wolf)
    End If


    c12 = Sqrt(1._dp-s12*s12)
    c23 = Sqrt(1._dp-s23*s23)
    c13 = Sqrt(1._dp-s13*s13)

    CKM(1,1) = c12 * c13
    CKM(1,2) = s12 * c13
    CKM(2,3) = s23 * c13
    CKM(3,3) = c23 * c13
    If (phase.Ne.0._dp) Then
     CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
     CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
     CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
    Else
     CKM(1,3) = s13
     CKM(2,1) = -s12*c23 -c12*s23*s13
     CKM(2,2) = c12*c23 -s12*s23*s13
     CKM(3,1) = s12*s23 -c12*c23*s13
     CKM(3,2) = -c12*s23 - s12*c23*s13
    End If

  End Subroutine Read_CKM

  Subroutine Read_SPINFO(io, kont)
  Implicit None
   Integer, Intent(in) :: io
   Integer, Intent(inout) :: kont

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line

     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_sp, read_line

     If (i_sp.Eq.1) Then
      sp_info = Trim(read_line)//" "//Trim(sp_info)
     Else If (i_sp.Eq.2) Then
      sp_info = Trim(sp_info)//" "//Trim(read_line)
     Else If (i_sp.Eq.4) Then ! there is some inconsistency, exit
      kont = -306
      Call AddError(306)
      Iname = Iname - 1
      Return
     Else
      Write(ErrCan,*) "Unknown entry in BLOCK SPINFO, is ignored:"
      Write(ErrCan,*) i_sp, read_line
     End If
    End Do

   200 Return

  End Subroutine Read_SPINFO

  Subroutine Read_MODSEL(io, i_particles, i_model, i_cpv, i_rp, kont)
  Implicit None
   Integer, Intent(in) :: io
   Integer, Intent(out) :: i_particles, i_model, i_cpv, i_rp
   Integer, Intent(inout) :: kont

   Integer :: i_mod, i_test
   Real(dp) :: r_mod
   Character(len=80) :: read_line

   i_cpv = 0
   i_rp = 0

    Do 
     Read(io,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_test,r_mod ! ,read_line
     If (i_test.Ne.12) Then
      Backspace(io)
      Read(io,*) i_test,i_mod ! ,read_line
     End If

     If (i_test.Eq.1) Then
      i_particles = i_test
      i_model = i_mod
      If ((i_mod.Lt.0).Or.(i_mod.Gt.3)) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "MSSM, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      Else If (i_mod.Eq.0) Then
       HighScaleModel = "MSSM"
      Else If (i_mod.Eq.1) Then
       HighScaleModel = "mSugra"
       check = SetHighScaleModel("SUGRA")
      Else If (i_mod.Eq.2) Then
       HighScaleModel = "GMSB"
       check = SetHighScaleModel("GMSB")
      Else If (i_mod.Eq.3) Then
       HighScaleModel = "AMSB"
       check = SetHighScaleModel("AMSB")
      End If

     Else If (i_test.Eq.3) Then
      External_Higgs = .False.
      mP03 = 0._dp
      mS03 = 0._dp
      mSpm = 0._dp
      RP03_save = 0._dp
      RS03_save = 0._dp
      If (i_mod.Eq.1) Then
       i_particles = i_test
       i_model = i_mod
       HighScaleModel = "NMSSM"
      Else If (i_mod.Eq.2) Then  ! adding SU(5)
       HighScaleModel = "SUGRA_SU5"
       check = SetHighScaleModel("SUGRA_SU5")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.3) Then  ! adding nu_R, one scale only
       HighScaleModel = "SUGRA_NuR1"
       check = SetHighScaleModel("SUGRA_NuR1")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.4) Then  ! adding nu_R, three scales
       HighScaleModel = "SUGRA_NuR"
       check = SetHighScaleModel("SUGRA_NuR")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.5) Then  ! adding Higgs triplet 
       HighScaleModel = "SEESAW_II"
       check = SetHighScaleModel("SEESAW_II")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.6) Then  ! adding one nu_R
       HighScaleModel = "NURRP1"
       check = SetHighScaleModel("NURRP1")
       i_particles = 2
       i_model = 4
      Else If (i_mod.Eq.7) Then  ! adding one nu_R, S, Phi
       HighScaleModel = "RPspon"
       check = SetHighScaleModel("RPspon")
       i_particles = 4
       i_model = 6
      Else If (i_mod.Eq.8) Then  ! adding two nu_R
       HighScaleModel = "NURRP2"
       check = SetHighScaleModel("NURRP2")
       i_particles = 5
       i_model = 5
! Florian Staub
      Else If (i_mod.Eq.10) Then  ! adding three 24-plet 
       HighScaleModel = "SEESAW_III_3G"
       check = SetHighScaleModel("SEESAW_III_3G")
       i_particles = 1
       i_model = 1
      Else If (i_mod.Eq.11) Then  ! adding one 15-plet
       HighScaleModel = "SEESAW_II_SARAH"
       check = SetHighScaleModel("SEESAW_II_SARAH")
       i_particles = 1
       i_model = 1
! Florian Staub
      Else If (i_mod.Ne.0) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "NMSSM, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.4) Then
      If (i_mod.Eq.1) Then
       i_rp = 1

      Else If (i_mod.Ne.0) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "Unknown entry for Block MODSEL ",i_test,i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.5) Then
      i_cpv = i_mod
      If ((i_mod.Lt.0).Or.(i_mod.Gt.2)) Then
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "Unknown entry for Block MODSEL ",i_test,i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.6) Then
      If (i_mod.Eq.0) Then
       GenerationMixing = .False.
      Else If ((i_mod.Ge.1).And.(i_mod.Le.3)) Then
       GenerationMixing = .True.
      Else
       Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
       Write(ErrCan,*) "GenerationMixing, Unknown entry for Block MODSEL ",i_mod
       kont = -302
       Call AddError(-kont)
       Return
      End If

     Else If (i_test.Eq.12) Then
      Call SetRGEScale(r_mod**2)  ! set Q_EWSB

     End If
    End Do ! i_mod

    If ((i_rp.Eq.1).And.(i_model.Eq.0)) Then
     HighScaleModel = "RPexplicit"
     check = SetHighScaleModel("RPexplicit")
    Else If ((i_rp.Eq.1).And.(i_model.Ge.1).And.(i_model.Le.3)) Then
     Add_Rparity = .True.
    End If

  End Subroutine Read_MODSEL

  Subroutine Read_SPhenoInput(io)
  Implicit None
   Integer, Intent(in) :: io

   Integer :: i_par
   Real(dp) :: wert
   Character(len=80) :: read_line

    ! This initialization is necessary for the arrar of production infos
    p_max = Size(Ecms)
    p_act = 0
    Ecms = 0._dp
    Pm = 0._dp
    Pp = 0._dp
    l_ISR = .False.
    Do 
     Read(io,*,End=200,err=200) read_line
!     Write(*,*) trim(read_line)
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*,End=200) i_par,wert ! ,read_line
!     write(*,*) i_par,wert,trim(read_line)
     Select Case(i_par)
     Case(1)
      ErrorLevel = Int(wert)

     Case(2)
      If (Int(wert).Ne.0) Then
       SPA_convention = .True.
       Call SetRGEScale(1.e3_dp**2)
      End If

     Case(3)
      If (Int(wert).Ne.0) Then 
       External_Spectrum = .True.
       External_Higgs = .True.
      End If

     Case(4)
      If (Int(wert).Ne.0) Use_Flavour_States = .True.

     Case(5)
      If (Int(wert).Ne.0) FermionMassResummation = .False.

     Case(11)  ! whether to calculate  branching ratios or not
      If (Int(wert).Eq.1) L_BR = .True.
      If (Int(wert).Eq.0) L_BR = .False.

     Case(12) ! minimal value such that a branching ratio is written out
      Call SetWriteMinBR(wert)

     Case(21)  ! whether to calculate cross sections or not
      If (Int(wert).Eq.1) L_CS = .True.
      If (Int(wert).Eq.0) L_CS = .False.

     Case(22) ! cms energy
      p_act = p_act + 1
      ! this test is necessary to avoid a memory violation
      If (p_act.Le.p_max) Then
       Ecms(p_act) = wert
      Else
       If (output_screen) &
           & Write(*,*) "The number of required points for the calculation"// &
           &  " of cross sections exceeds",p_max
       If (output_screen) &
           & Write(*,*) "Ignoring this information"
       If (output_screen) &
     &  Write(*,*) "Please enlarge the corresponding arrays in the main program."
       Write(ErrCan,*) "The number of required points for the calculation"// &
               &   " of cross sections exceeds",p_max
       Write(ErrCan,*) "Ignoring this information"
       Write(ErrCan,*) &
         &"Please enlarge the corresponding arrays in the main program."
      End If

     Case (23) ! polarisation of incoming e- beam
      If (Abs(wert).Gt.1._dp) Then
       If (output_screen) Write(*,*) &
           & "e- beam polarisation has to between -1 and 1 and not",wert
       If (output_screen) &
           & Write(*,*) "using now unpolarised e- beam"
       Write(ErrCan,*) &
          & "e- beam polarisation has to between -1 and 1 and not",wert
       Write(ErrCan,*) "using now unpolarised e- beam"
       If (p_act.Le.p_max) Pm(p_act) = 0
      Else
       If (p_act.Le.p_max) Pm(p_act) = wert
      End If

     Case (24) ! polarisation of incoming e+ beam
      If (Abs(wert).Gt.1._dp) Then
       If (output_screen) Write(*,*) &
           & "e+ beam polarisation has to between -1 and 1 and not",wert
       If (output_screen) &
           & Write(*,*) "using now unpolarised e+ beam"
       Write(ErrCan,*) &
          & "e+ beam polarisation has to between -1 and 1 and not",wert
       Write(ErrCan,*) "using now unpolarised e+ beam"
       If (p_act.Le.p_max) Pp(p_act) = 0
      Else
       If (p_act.Le.p_max) Pp(p_act) = wert
      End If

     Case(25)
      If ((wert.Eq.1._dp).And.(p_act.Le.p_max)) L_ISR(p_act) = .True.

     Case(26) ! minimal value such that a cross section is written out
      Call SetWriteMinSig(wert)

     Case(31) ! setting a fixed GUT scale if wert > 0
      If (wert.Gt.0._dp) Call SetGUTScale(wert)

     Case(32) ! requires strict unification
      If (Int(wert).Ne.0) check = SetStrictUnification(.True.)

     Case(34) ! precision of mass calculation
      delta_mass = wert

     Case(35) ! maximal number of iterations
      n_run = Int(wert)

     Case(36) ! write out debug information
      If (wert.Eq.0) Then
       WriteOut = .False.
      Else
       WriteOut = .True.
      End If

     Case(37) ! if =1 -> CKM thourgh V_u, if =2 CKM through V_d 
      If ((wert.Eq.1._dp).Or.(wert.Eq.2._dp)) i1 =  SetYukawaScheme(Int(wert))

     Case(38) ! set looplevel of RGEs
      If (wert.Ne.2._dp) Then
       TwoLoopRGE=.False.
      Else
       TwoLoopRGE=.True.
      End If

     Case(39) ! write additional SLHA1 file
      If (wert.Eq.1._dp) Write_SLHA1 = .True.

     Case(40) ! alpha(0)
      check_alpha(2) = .True.
      Alpha = 1._dp / wert

     Case(41) ! Z-boson width
      gamZ = wert

     Case(42) ! W-boson width
      gamW = wert


     Case(80) ! exit for sure with non-zero value if a problem occurs
      If (wert.Eq.1) Non_Zero_Exit = .True.      

     Case(90) ! add R-parity at low energies
      If (wert.Eq.1) Add_Rparity = .True.      

     Case(91) ! fit RP parameters such, that neutrino data are o.k.
      If (wert.Eq.1) l_fit_RP_parameters = .True.      

     Case(92) ! for Pythia input
      If (wert.Eq.1) l_RP_Pythia = .True.      

     Case(100) ! use bsstep instead of 
      If (wert.Eq.1) test_l = Set_Use_bsstep_instead_of_rkqs(.True.)

     Case(101) ! use bsstep instead of 
      If (wert.Eq.1) test_l = Set_Use_rzextr_instead_of_pzextr(.True.)

     Case(110) ! write output for LHC observables
      If (wert.Eq.1) Then
       LWrite_LHC_Observables = .True.
      Else
       LWrite_LHC_Observables = .False.
      End If

     Case Default
      If (output_screen) Write(*,*) &
           & "Problem while reading SPhenoInput, ignoring unknown entry" &
           & ,i_par,wert
      Write(Errcan,*) &
          & "Problem while reading  SPhenoInput, ignoring unknown entry" &
               & ,i_par,wert
     End Select ! i_par

    End Do  ! i_par 

   200 Return

  End Subroutine Read_SPhenoInput

  Subroutine Read_EXTPAR(io, i_c, i_model, set_mod_par, kont)
  Implicit None
   Integer, Intent(in) :: io, i_c, i_model
   Integer, Intent(inout) :: kont, set_mod_par(:)

   Integer :: i_par
   Real(dp) :: wert
   Character(len=80) :: read_line

    If (i_model.Lt.0) Then ! check if model is already defined
     Write(ErrCan,*) &
     "You must first specify the model before the model parameters can be set."
     kont = -303
     Call AddError(-kont)
     Return
    End If
!    Write(ErrCan,*) "Reading EXTPAR"

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) trim(read_line)
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_par,wert ! ,read_line

     If ((i_par.Eq.0).And.(i_c.Eq.0)) Then
      If (i_model.Eq.0) Call SetRGEScale(wert**2)  ! in case of MSSM
      If (i_model.Eq.1) Call SetGUTScale(wert)     ! Sugra
      If (i_model.Eq.3) Call SetGUTScale(wert)     ! AMSB
     Else If (i_par.Eq.1) Then 
      If (i_c.Eq.0) Mi(1) = Cmplx(wert, Aimag(Mi(1)), dp) 
      If (i_c.Eq.0) Mi_0(1) = Cmplx(wert, Aimag(Mi_0(1)), dp) 
      If (i_c.Eq.1) Mi(1) = Cmplx(Real(Mi(1),dp), wert, dp)
      If (i_c.Eq.1) Mi_0(1) = Cmplx(Real(Mi_0(1),dp), wert, dp)  
      set_mod_par(1) = 1
     Else If (i_par.Eq.2) Then 
      If (i_c.Eq.0) Mi(2) = Cmplx(wert, Aimag(Mi(2)), dp)
      If (i_c.Eq.0) Mi_0(2) = Cmplx(wert, Aimag(Mi_0(2)), dp)
      If (i_c.Eq.1) Mi(2) = Cmplx(Real(Mi(2),dp), wert, dp)
      If (i_c.Eq.1) Mi_0(2) = Cmplx(Real(Mi_0(2),dp), wert, dp)  
      set_mod_par(2) = 1
     Else If (i_par.Eq.3) Then 
      If (i_c.Eq.0) Mi(3) = Cmplx(wert, Aimag(Mi(3)), dp) 
      If (i_c.Eq.0) Mi_0(3) = Cmplx(wert, Aimag(Mi_0(3)), dp) 
      If (i_c.Eq.1) Mi(3) = Cmplx(Real(Mi(3),dp), wert, dp)
      If (i_c.Eq.1) Mi_0(3) = Cmplx(Real(Mi_0(3),dp), wert, dp)  
      set_mod_par(3) = 1
     Else If (i_par.Eq.11) Then 
      If (i_c.Eq.0) AoY_u(3,3) = Cmplx(wert, Aimag(AoY_u(3,3)), dp) 
      If (i_c.Eq.1) AoY_u(3,3) = Cmplx(Real(AoY_u(3,3),dp), wert, dp) 
      At_save = AoY_u(3,3)
      AoY_u_0 = AoY_u 
      set_mod_par(4) = 1
     Else If (i_par.Eq.12) Then 
      If (i_c.Eq.0) AoY_d(3,3) = Cmplx(wert, Aimag(AoY_d(3,3)), dp)
      If (i_c.Eq.1) AoY_d(3,3) = Cmplx(Real(AoY_d(3,3),dp), wert, dp) 
      Ab_save = AoY_d(3,3)
      AoY_d_0 = AoY_d 
      set_mod_par(5) = 1
     Else If (i_par.Eq.13) Then 
      If (i_c.Eq.0) AoY_l = Cmplx(wert, Aimag(AoY_l(3,3)), dp) 
      If (i_c.Eq.1) AoY_l = Cmplx(Real(AoY_l(3,3),dp), wert, dp) 
      Atau_save = AoY_l(3,3)
      AoY_l_0 = AoY_l 
      set_mod_par(6) = 1
     Else If ((i_par.Eq.21).And.(i_c.Eq.0)) Then 
      M2_H(1) = wert
      M2_H_0(1) = wert
      set_mod_par(7) = 1

     Else If ((i_par.Eq.22).And.(i_c.Eq.0)) Then 
      M2_H(2) = wert
      M2_H_0(2) = wert
      set_mod_par(8) = 1

     Else If (i_par.Eq.23) Then
      If ((i_model.Eq.0).Or.(HighScaleModel.Eq."NMSSM")) Then 
       If (i_c.Eq.0) mu = Cmplx(wert, Aimag(mu), dp)
       If (i_c.Eq.1) mu = Cmplx(Real(mu,dp), wert, dp)
       set_mod_par(9) = 1
      Else
       Write(ErrCan,*) "mu can only be specified in the general MSSM and is"
       Write(ErrCan,*) "ignored for the high scale model SUGRA, GMSB and AMSB"
      End If

     Else If ((i_par.Eq.24).And.(i_c.Eq.0)) Then 
      If (i_model.Eq.0) Then 
       mP02(2) = wert
       If (mP02(2).Ge.0._dp) mP0(2) = Sqrt(mP02(2))
       set_mod_par(10) = 1
      Else
       Write(ErrCan,*) "m_A0(Q) can only be specified in the general MSSM and is"
       Write(ErrCan,*) "ignored for high scale models like SUGRA, GMSB or AMSB"
      End If
 
     Else If ((i_par.Eq.25).And.(i_c.Eq.0)) Then 
      tanb = wert

     Else If ((i_par.Eq.26).And.(i_c.Eq.0)) Then 
      If (i_model.Eq.0) Then 
       HighScaleModel = "pMSSM"
       check = SetHighScaleModel("pMSSM")
       mP0(2) = wert
       mP02(2) = wert**2
       set_mod_par(10) = 1
      Else If (HighScaleModel.Eq."RPspon") Then 
       mP0(2) = wert
       mP02(2) = wert**2
      Else
       Write(ErrCan,*) "m_A0 can only be specified in the general MSSM and is"
       Write(ErrCan,*) "ignored for high scale models like SUGRA, GMSB or AMSB"
      End If

     Else If ((i_par.Eq.31).And.(i_c.Eq.0)) Then 
      M2_L(1,1) = wert**2
      M2_L_0(1,1) = M2_L(1,1)
      set_mod_par(11) = 1
     Else If ((i_par.Eq.32).And.(i_c.Eq.0)) Then 
      M2_L(2,2) = wert**2
      M2_L_0(2,2) = M2_L(2,2)
      set_mod_par(12) = 1
     Else If ((i_par.Eq.33).And.(i_c.Eq.0)) Then 
      M2_L(3,3) = wert**2
      M2_L_0(3,3) = M2_L(3,3)
      set_mod_par(13) = 1
     Else If ((i_par.Eq.34).And.(i_c.Eq.0)) Then 
      M2_E(1,1) = wert**2
      M2_E_0(1,1) = M2_E(1,1)
      set_mod_par(14) = 1
     Else If ((i_par.Eq.35).And.(i_c.Eq.0)) Then 
      M2_E(2,2) = wert**2
      M2_E_0(2,2) = M2_E(2,2)
      set_mod_par(15) = 1
     Else If ((i_par.Eq.36).And.(i_c.Eq.0)) Then 
      M2_E(3,3) = wert**2
      M2_E_0(3,3) = M2_E(3,3)
      set_mod_par(16) = 1
     Else If ((i_par.Eq.41).And.(i_c.Eq.0)) Then 
      M2_Q(1,1) = wert**2
      M2_Q_0(1,1) = M2_Q(1,1)
      set_mod_par(17) = 1
     Else If ((i_par.Eq.42).And.(i_c.Eq.0)) Then 
      M2_Q(2,2) = wert**2
      M2_Q_0(2,2) = M2_Q(2,2)
      set_mod_par(18) = 1
     Else If ((i_par.Eq.43).And.(i_c.Eq.0)) Then 
      M2_Q(3,3) = wert**2
      M2_Q_0(3,3) = M2_Q(3,3)
      set_mod_par(19) = 1
     Else If ((i_par.Eq.44).And.(i_c.Eq.0)) Then 
      M2_U(1,1) = wert**2
      M2_U_0(1,1) = M2_U(1,1)
      set_mod_par(20) = 1
     Else If ((i_par.Eq.45).And.(i_c.Eq.0)) Then 
      M2_U(2,2) = wert**2
      M2_U_0(2,2) = M2_U(2,2)
      set_mod_par(21) = 1
     Else If ((i_par.Eq.46).And.(i_c.Eq.0)) Then 
      M2_U(3,3) = wert**2
      M2_U_0(3,3) = M2_U(3,3)
      set_mod_par(22) = 1
     Else If ((i_par.Eq.47).And.(i_c.Eq.0)) Then 
      M2_D(1,1) = wert**2
      M2_D_0(1,1) = M2_D(1,1)
      set_mod_par(23) = 1
     Else If ((i_par.Eq.48).And.(i_c.Eq.0)) Then 
      M2_D(2,2) = wert**2
      M2_D_0(2,2) = M2_D(2,2)
      set_mod_par(24) = 1
     Else If ((i_par.Eq.49).And.(i_c.Eq.0)) Then 
      M2_D(3,3) = wert**2
      M2_D_0(3,3) = M2_D(3,3)
      set_mod_par(25) = 1

     !------------------------------
     ! NMSSM
     !------------------------------
     Else If (i_par.Eq.61) Then 
      If (i_c.Eq.0) h0 = Cmplx(wert, Aimag(h0),dp) 
      If (i_c.Eq.1) h0 = Cmplx(Real(h0,dp), wert, dp)
     Else If (i_par.Eq.62) Then 
      ! different convention with respect to Cyril
      If (i_c.Eq.0) lam = Cmplx(2._dp * wert, Aimag(lam),dp)
      If (i_c.Eq.1) lam = Cmplx(Real(lam,dp), 2._dp * wert, dp)
     Else If (i_par.Eq.63) Then 
      If (i_c.Eq.0) Ao_h0 = Cmplx(wert, Aimag(Ao_h0),dp) 
      If (i_c.Eq.1) Ao_lam = Cmplx(Real(Ao_lam,dp), wert, dp)
     Else If (i_par.Eq.64) Then 
      If (i_c.Eq.0) Ao_lam = Cmplx(wert, Aimag(Ao_lam),dp) 
      If (i_c.Eq.1) Ao_lam = Cmplx(Real(Ao_lam,dp), wert, dp)
     Else If (i_par.Eq.65) Then
      If ((HighScaleModel.Eq."NMSSM").Or.(HighScaleModel.Eq."RPspon")) Then 
       If (i_c.Eq.0) lam_vS = Cmplx(wert, Aimag(lam_vS),dp)
       If (i_c.Eq.1) lam_vS = Cmplx(Real(lam_vS,dp), wert, dp)
       set_mod_par(9) = 1
      Else
       Write(ErrCan,*) "Attempt to use i_par == 65"
       Write(ErrCan,*) "mu can only be specified in this way in the NMSSM and"
       Write(ErrCan,*) "is ignored in model "//Trim(HighScaleModel)
      End If
     
     !---------------------------------------
     ! explicit R-parity breaking a la Munoz
     ! and spontaneous R-parity breaking
     !---------------------------------------
      Else If (i_par.Eq.71) Then
        M2_P = wert
      Else If (i_par.Eq.72) Then 
       vP = wert ! has to be replaced by v_R later in case of Munoz
      Else If (i_par.Eq.73) Then 
       Y_nu(1,1) = wert ! has to be replaced by hnu(1) later
      Else If (i_par.Eq.74) Then 
       Y_nu(1,2) = wert ! has to be replaced by hnu(1) later
      Else If (i_par.Eq.75) Then 
       Y_nu(1,3) = wert ! has to be replaced by hnu(1) later
      Else If (i_par.Eq.76) Then 
       AoY_nu(1,1) = wert ! has to be replaced by AoY_nu(1) later
      Else If (i_par.Eq.77) Then 
       AoY_nu(1,2) = wert ! has to be replaced by AoY_nu(1) later
      Else If (i_par.Eq.78) Then 
       AoY_nu(1,3) = wert ! has to be replaced by AoY_nu(1) later
     !---------------------------------------
     ! spontaneous R-parity breaking
     !---------------------------------------
      Else If (i_par.Eq.81) Then
        M2_S = wert
      Else If (i_par.Eq.82) Then
        M2_R(1,1) = wert
      Else If (i_par.Eq.83) Then 
       vS = wert 
      Else If (i_par.Eq.84) Then 
       vR = wert  
      Else If (i_par.Eq.85) Then 
       h_pns = wert  
      Else If (i_par.Eq.86) Then 
       Ao_hpns = wert  
      Else If (i_par.Eq.87) Then 
       M_rs = wert  
      Else If (i_par.Eq.88) Then 
       M_phi = wert  
      Else If (i_par.Eq.89) Then 
       BM_rs = wert  
      Else If (i_par.Eq.90) Then 
       BM_phi = wert  
     !---------------------------------------
     ! explicit R-parity breaking a la Munoz
     ! and spontaneous R-parity breaking
     !---------------------------------------
      Else If (i_par.Eq.91) Then
       If (i_c.Eq.0) h02 = Cmplx(wert, Aimag(h02),dp) 
       If (i_c.Eq.1) h02 = Cmplx(Real(h02,dp), wert, dp)
      Else If (i_par.Eq.92) Then 
        ! different convention  with respect to Cyril
       If (i_c.Eq.0) lam2 = Cmplx(2._dp * wert, Aimag(lam2),dp) 
       If (i_c.Eq.1) lam2 = Cmplx(Real(lam2,dp), 2._dp * wert , dp)
      Else If (i_par.Eq.93) Then 
       If (i_c.Eq.0) lam112 = Cmplx(wert, Aimag(lam112),dp) 
       If (i_c.Eq.1) lam112 = Cmplx(Real(lam112,dp), wert, dp)
      Else If (i_par.Eq.94) Then 
       If (i_c.Eq.0) lam122 = Cmplx(wert, Aimag(lam122),dp) 
       If (i_c.Eq.1) lam122 = Cmplx(Real(lam122,dp), wert, dp)
      Else If (i_par.Eq.95) Then 
       If (i_c.Eq.0) Ao_h02 = Cmplx(wert, Aimag(Ao_h02),dp)
       If (i_c.Eq.1) Ao_h02 = Cmplx(Real(Ao_h02,dp), wert, dp)
      Else If (i_par.Eq.96) Then 
       If (i_c.Eq.0) Ao_lam222 = Cmplx(wert, Aimag(Ao_lam222),dp) 
       If (i_c.Eq.1) Ao_lam222 = Cmplx(Real(Ao_lam222,dp), wert, dp)
      Else If (i_par.Eq.97) Then 
       If (i_c.Eq.0) Ao_lam112 = Cmplx(wert, Aimag(Ao_lam112),dp)
       If (i_c.Eq.1) Ao_lam112 = Cmplx(Real(Ao_lam112,dp), wert, dp) 
      Else If (i_par.Eq.98) Then 
       If (i_c.Eq.0) Ao_lam122 = Cmplx(wert, Aimag(Ao_lam122),dp)
       If (i_c.Eq.1) Ao_lam122 = Cmplx(Real(Ao_lam122,dp), wert, dp)
      Else If (i_par.Eq.99) Then 
       M2_R(2,2) = wert**2
      Else If (i_par.Eq.100) Then 
       vP2 = wert ! has to be replaced by v_R2 later
      Else If (i_par.Eq.101) Then 
       Y_nu(2,1) = wert ! has to be replaced by hnu2(1) later
      Else If (i_par.Eq.102) Then 
       Y_nu(2,2) = wert ! has to be replaced by hnu2(2) later
      Else If (i_par.Eq.103) Then 
       Y_nu(2,3) = wert ! has to be replaced by hnu2(3) later
      Else If (i_par.Eq.104) Then 
       AoY_nu(2,1) = wert ! has to be replaced by AoY_nu2(1) later
      Else If (i_par.Eq.105) Then 
       AoY_nu(2,2) = wert ! has to be replaced by AoY_nu2(2) later
      Else If (i_par.Eq.106) Then 
       AoY_nu(2,3) = wert ! has to be replaced by AoY_nu2(3) later

! Florian Staub, Seesaw II+III
      Else If (i_par.Eq.200) Then 
       If (i_c.Eq.0) MTM0 = Cmplx(wert, Aimag(MTM0),dp)
       If (i_c.Eq.1) MTM0 = Cmplx(Real(MTM0,dp), wert,dp)
       If (i_c.Eq.0) Then
        m_H3 = wert
       ! as initalization, will be computed more precisely later
        MS15_mH3 = wert
        MT15_mH3 = wert
        MZ15_mH3 = wert
       End If
      Else If (i_par.Eq.202) Then 
       If (i_c.Eq.0) Lambda1_gut = Cmplx(wert, Aimag(Lambda1_gut),dp)
       If (i_c.Eq.1) Lambda1_gut = Cmplx(Real(Lambda1_gut,dp), wert,dp)
       lam12_0(1) = Lambda1_gut
      Else If (i_par.Eq.203) Then 
       If (i_c.Eq.0) Lambda2_gut = Cmplx(wert, Aimag(Lambda2_gut),dp)
       If (i_c.Eq.1) Lambda2_gut = Cmplx(Real(Lambda2_gut,dp), wert,dp)
       lam12_0(2) = Lambda2_gut

! Florian Staub, Seesaw II+III

     Else
      If (i_c.Eq.0) Then
       If (output_screen)  Write(*,*) &
        & "Problem while reading EXTPAR, ignoring unknown (unsupported) entry"&
        & ,i_par,wert
       Write(Errcan,*) &
       & "Problem while reading EXTPAR, ignoring unknown  (unsupported) entry"&
       & ,i_par,wert
      Else If (i_c.Eq.1) Then
       If (output_screen)  Write(*,*) &
      & "Problem while reading IMEXTPAR, ignoring unknown (unsupported) entry"&
        & ,i_par,wert
       Write(Errcan,*) &
     & "Problem while reading IMEXTPAR, ignoring unknown  (unsupported) entry"&
       & ,i_par,wert
      End If
     End If
!     Write(errcan,*) i_par,wert
    End Do  ! i_par 

    !----------------------------------------------------------------
    ! check if T_f and A_f given, if yes, then A_f gets overwritten
    !----------------------------------------------------------------
    If (A_u(3,3).ne.ZeroC) At_save = ZeroC
    If (A_d(3,3).ne.ZeroC) Ab_save = ZeroC
    If (A_l(3,3).ne.ZeroC) Atau_save = ZeroC
    200 Return

  End Subroutine Read_EXTPAR 


  Subroutine Read_MINPAR(io, i_c, i_model, set_mod_par, kont)
  Implicit None
   Integer, Intent(in) :: io, i_c, i_model
   Integer, Intent(inout) :: kont, set_mod_par(:)

   Integer :: i_par
   Real(dp) :: wert
   Character(len=80) :: read_line

    If (i_model.Lt.0) Then ! check if model is already defined
     Write(ErrCan,*) &
     "You must first specify the model before the model parameters can be set."
     kont = -303
     Call AddError(-kont)
     Return
    End If

    Do 
     Read(io,*,End=200) read_line
!     Write(*,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_par,wert ! ,read_line
! Write(*,*) i_c,i_par,wert ! ,read_line
! Write(*,*) i_model
     If ((i_par.Eq.1).And.(i_c.Eq.0)) Then 
      If (i_model.Eq.1) Then ! mSugra, M_0
       set_mod_par(1) = 1
       M2_D_0 = 0._dp
       Do i1=1,3
        M2_D_0(i1,i1) = wert**2
       End Do
       M2_E_0 = M2_D_0
       M2_L_0 = M2_D_0
       M2_R_0 = M2_D_0
       M2_Q_0 = M2_D_0
       M2_U_0 = M2_D_0
       M2_H_0 = wert**2
       M2_T_0 = M2_H_0
      Else If (i_model.Eq.2) Then ! GMSB, Lambda
       set_mod_par(1) = 1
       Lambda = wert
      Else If (i_model.Eq.3) Then ! AMSB, common scalar mass
       set_mod_par(1) = 1
       M0_amsb = wert
      End If

     Else If (i_par.Eq.2) Then 
      If (i_model.Eq.1) Then ! mSugra, M_1/2
       set_mod_par(2) = 1
       If (i_c.Eq.0) Mi_0 =  Cmplx(wert, Aimag(Mi_0),dp) 
       If (i_c.Eq.1) Mi_0 =  Cmplx(Real(Mi_0,dp),  wert, dp)
      Else If ((i_model.Eq.2).And.(i_c.Eq.0)) Then ! GMSB, M_M
       set_mod_par(2) = 1
       MlambdaS = wert
      Else If ((i_model.Eq.3).And.(i_c.Eq.0)) Then ! AMSB, Gravitino mass
       set_mod_par(2) = 1
       M_32 = wert
      End If

     Else If ((i_par.Eq.3).And.(i_c.Eq.0)) Then 
      If ((i_model.Ge.0).And.(i_model.Le.3)) Then ! MSSM, mSugra, GMSB, AMSB
       set_mod_par(3) = 1
       tanb = wert
       tanb_mZ = wert
      End If

     Else If (i_par.Eq.4) Then 
      If ((i_model.Ge.0).And.(i_model.Le.3)) Then ! MSSM, mSugra, GMSB, AMSB, sign_mu
       set_mod_par(4) = 1
       If (i_c.Eq.0) phase_mu = Cmplx(wert, Aimag(phase_mu),dp)
       If (i_c.Eq.1) phase_mu = Cmplx(Real(phase_mu, dp),  wert, dp)
      End If

     Else If (i_par.Eq.5) Then 
      If (i_model.Eq.1) Then ! mSugra, A_0
       set_mod_par(5) = 1
       If (i_c.Eq.0) AoY_d_0 = Cmplx(wert, Aimag(AoY_d_0),dp) 
       If (i_c.Eq.1) AoY_d_0 = Cmplx( Real(AoY_d_0,dp) , wert, dp)
       AoY_l_0 = AoY_d_0
       AoY_u_0 = AoY_d_0
       AoY_nu_0 = AoY_d_0
       AoT_0 = AoY_d_0
       Aolam12_0 = AoY_d_0(1,1)
!       If (i_c.Eq.0) Alam12_0 =  Cmplx(wert, Aimag(Alam12_0),dp)
!       If (i_c.Eq.1) Alam12_0 =  Cmplx(Real(Alam12_0,dp) , wert, dp)
      Else If ((i_model.Eq.2).And.(i_c.Eq.0)) Then ! GMSB, n_5
       set_mod_par(5) = 1
       n5plets = wert
       n10plets = 0
      End If

     Else If ((i_par.Eq.6).And.(i_c.Eq.0)) Then 
      If (i_model.Eq.2) Then ! GMSB, Gravitino mass factor
       set_mod_par(6) = 1
       grav_fac = wert
      End If

     Else If ((i_par.Eq.7).And.(i_c.Eq.0)) Then ! SUGRA_SU5, SO(10) scale
      M_SO_10 = wert

     Else If ((i_par.Eq.8).And.(i_c.Eq.0)) Then ! SUGRA_SU5 or SUGRA_NuR, D-term
      D_SO_10 = wert

     Else If (i_par.Eq.9) Then ! SUGRA_SU5, real(lambda(m_GUT))
      If (i_c.Eq.0) lam_0 = Cmplx(wert, Aimag(lam_0),dp) 
      If (i_c.Eq.1) lam_0 = Cmplx( Real(lam_0,dp), wert, dp)
 
     Else If (i_par.Eq.10) Then ! SUGRA_SU5, real(lambda'(m_GUT))
      If (i_c.Eq.0) lamp_0 = Cmplx(wert, Aimag(lamp_0), dp) 
      If (i_c.Eq.1) lamp_0 = Cmplx( Real(lamp_0,dp), wert, dp) 
 
     Else 
      Write(ErrCan,*) "Error in routine "//NameOfUnit(Iname)
      If (i_c.Eq.0) Write(ErrCan,*) "Unknown entry for Block MINPAR ",i_par 
      If (i_c.Eq.1) Write(ErrCan,*) "Unknown entry for Block IMMINPAR ",i_par 
      Call AddError(304)
      If (ErrorLevel.Eq.2) Call TerminateProgram
     End If 

    End Do ! i_par

  200  Return

  End Subroutine Read_MINPAR

  Subroutine Read_SMinput(io)
  Implicit None
   Integer, Intent(in) :: io
   
   Integer :: i_sm
   Real(dp) :: wert
   Character(len=80) :: read_line

    Do 
     Read(io,*) read_line
     If (read_line(1:1).Eq."#") Cycle ! this loop
     Backspace(io) ! resetting to the beginning of the line
     If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Exit ! this loop

     Read(io,*) i_sm,wert ! ,read_line

     Select Case(i_sm)
     Case(1)
      check_alpha(1) = .True.
      MZ_input = .True.
      Alpha_MZ_MS = 1._dp / wert

     Case(2) ! G_F
      G_F = wert

     Case(3) ! alpha_s(m_Z)
      alphaS_mZ = wert

     Case(4) ! m_Z
      mZ = wert
      mZ2 = mZ**2
      calc_ferm = .True.

     Case(5) ! m_b(m_b)^MSbar
      mf_d(3) = wert
      mf_d2(3) = mf_d(3)**2
      calc_ferm = .True.

     Case(6) ! m_t^pole
      mf_u(3) = wert
      mf_u2(3) = mf_u(3)**2

     Case(7) ! m_tau^pole
      mf_l(3) = wert
      mf_l2(3) = mf_l(3)**2
      calc_ferm = .True.

     Case(8) ! m_nu_3, input is in GeV
      Mf_nu(3) = wert

     Case(11) ! electron mass
      mf_l(1) = wert
      mf_l2(1) = wert**2
      calc_ferm = .True.

     Case(12) ! m_nu_1, input is in GeV
      Mf_nu(1) = wert 

     Case(13) ! muon mass
      mf_l(2) = wert
      mf_l2(2) = wert**2
      calc_ferm = .True.

     Case(14) ! m_nu_2, input is in eV, transform to GeV
      Mf_nu(2) = wert 

     Case(21) ! d-quark mass at 2 GeV
      mf_d(1) = wert
      mf_d2(1) = wert**2
      calc_ferm = .True.

     Case(22) ! u-quark mass at 2 GeV
      mf_u(1) = wert
      mf_u2(1) = wert**2
      calc_ferm = .True.

     Case(23) ! s-quark mass at 2 GeV
      mf_d(2) = wert
      mf_d2(2) = wert**2
      calc_ferm = .True.

     Case(24) ! c-quark mass at Q=m_c
      mf_u(2) = wert
      mf_u2(2) = wert**2
      calc_ferm = .True.

     Case Default
      If (output_screen) &
           & Write(*,*) "Ignoring unknown entry for Block SMINPUTS ",i_sm 
      Write(ErrCan,*) "Ignoring unknown entry for Block SMINPUTS ",i_sm 
     End Select

    End Do ! i_sm

  End Subroutine Read_SMinput
 
 End  Subroutine LesHouches_Input


 Subroutine LesHouches_Out(io_L, io, kont, HighScaleModel, M_GUT           &
      & , BRbtosgamma, Bs_MuMu, DeltaMBd, DeltaMBs, BrBToSLL, BtoSNuNu     &
      & , BR_Bu_TauNu  &
      & , a_e, a_mu, a_tau, d_e, d_mu, d_tau, BrMuToEGamma, BrTauToEGamma  &
      & , BrTauToMuGamma, BrMu3e, BrTau3e, BrTau3Mu, BR_Z_e_mu, BR_Z_e_tau &
      & , BR_Z_mu_tau  &
      & , Rho_parameter, Ecms, Pm, Pp, ISR, SigSup, SigSdown, SigSle &
      & , SigSn, SigChi0, SigC, SigS0, SigSP, SigHp, omega, f_name)
 !--------------------------------------------------------------------
 ! writes out data using the Les Houches standard as defined in 
 ! hep-ph/0311123
 ! input:
 !   - HighScaleModel ........ string specifiying the model 
 !   - M_GUT ................. scale of high scale theory, in general a GUT theory
 !   - BRbtosgamma ........... 10^4 BR(b -> s gamma)
 !   - a_mu .................. SUSY contribution to (g-2)_muon
 !   - Rho_parameter ......... rho parameter
 !   - Ecms .................. center of mass energy in GeV
 !   - Pm .................... degree of polarisation of incoming electrons
 !   - Pp .................... degree of polarisation of incoming positrons
 !   - ISR ................... if .true. then calculate initial state rediation
 !                             stemming from the incoming elctron/positron beams
 !   - Sigsup ................ cross sections of u-type squarks
 !   - SigSdown .............. cross sections of d-type squarks
 !   - SigSle ................ cross sections of sleptons
 !   - SigSn ................. cross sections of sneutrinos
 !   - SigChi0 ............... cross sections of neutralinos
 !   - SigC .................. cross sections of charginos
 !   - SigSP ................. cross sections of neutral Higgs bosons
 !   - SigHp ................. cross section of charged Higgs boson
 !   - omega ................. relic density omega h^2, optional
 !   - f_name ................ alternative name of output file, optional
 !--------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: io_L, io, kont
  Real(dp), Intent(in) :: M_GUT, BRbtosgamma, Bs_MuMu, BrBToSLL, BR_Bu_TauNu &
      & , BtoSNuNu, a_e, a_mu, a_tau, d_e, d_mu, d_tau, BrMuToEGamma         &
      & , BrTauToEGamma, BrTauToMuGamma, BrMu3e, BrTau3e, BrTau3Mu        &
      & , BR_Z_e_mu, BR_Z_e_tau, BR_Z_mu_tau, Rho_parameter   &
      & , Ecms(:), Pm(:) &
      & , Pp(:), SigSup(:,:,:), SigSdown(:,:,:), SigSle(:,:,:), SigSn(:,:,:) &
      & , SigChi0(:,:,:), SigC(:,:,:), SigS0(:,:), SigSP(:,:,:), SigHp(:,:,:)
  Complex(dp), Intent(in) :: DeltaMBd, DeltaMBs
  Character(len=15), Intent(in) :: HighScaleModel
  Logical, Intent(in) :: ISR(:)
  Real(dp), Intent(in), Optional :: Omega
  Character(len=*), Intent(in), Optional :: f_name

  Integer :: i1, i2, i3, i4, id_sle(6), id_f, id_fp, id_su(6) &
      & , id_sd(6), i_zaehl, n_min, c_min, ii, jj
  Integer :: id_d_2(200,2), id_d_3(400,3), i_c2
  Complex(dp) :: nr(10,10)
  Real(dp) :: mnr(10), Q, BRtot, RG0(3,3), Brmin100, mat6R(6,6), mat3R(3,3) &
      & , mat5R(5,5), mat8R(8,8)
  Integer, Parameter :: n_max=500
  Real(dp), Dimension(n_max) :: Br, gP
  Real(dp) :: gT, MaxCont
  Character(len=10) :: name, c_m, c_snu(3), c_sle(6)     &
      & , c_su(6), c_sd(6), c_slep(6), c_grav, zeit
  Character(len=6) :: c_lm(3), c_lp(3), c_nu(3), c_gl, c_cp(5) &
     & , c_cm(5), c_c0(8), c_phot, c_sp(7), c_sm(7)
  Character(len=1) :: c_u(3), c_d(3)
  Character(len=4) :: c_S0(6), c_P0(5)
  Character(len=8) :: Datum
  Character(len=13) :: c_sfermion 
  Character(len=30), Dimension(n_max) :: Fnames, Lnames
  Integer :: n_n, n_c, n_s0, n_p0, n_spm, n_sl, n_sn, n_sd, n_su 
  Integer :: id_n(8), id_sp(7), id_sm(7), id_S0(6), id_P0(5), id_c(5) &
      & , id_snu(3)
  Integer, Parameter ::id_A0 = 36, id_Hp = 37                            &
      & , id_W = 24, id_Z = 23, id_ph = 22, id_gl = 21                   &
      & , id_l(3) = (/ 11, 13, 15 /), id_nu(3) = (/ 12, 14, 16 /)        &
      & , id_u(3) = (/ 2, 4, 6 /), id_d(3) = (/ 1, 3, 5 /)               &
      & , id_grav = 1000039, id_glu = 1000021
  Logical :: non_minimal ! checks if there is deviations from the minimal models
  Logical :: l_3body     ! checks if 3-body modes have been calculated
  Logical, Save :: l_open = .True. ! in case of a loop I want to open the
                                   ! output only once 
  !-------------------------------------------------------------------
  ! these variables are useful for the case of R-parity violation
  !------------------------------------------------------------------
  Integer :: delta_n_rp, delta_c_rp, i_eff
  Logical :: file_exists,  use_new_RP
  !--------------------------------------------------------------------- 
  ! mixing matrices for shifts to super-CKM and super-PMNS basis 
  !--------------------------------------------------------------------- 
  Integer :: ierr, i_errors(1100), io_L1, id_check(2)
  Real(dp) :: Yu(3), Yd(3), Yl(3), m_save
  Complex(dp) :: Rsave(6)
  Complex(dp), Dimension(3,3) :: CKM_Q, PMNS_Q, RSn_pmns
  Complex(dp), Dimension(6,6) :: RUsq_ckm, RDsq_ckm, RSl_pmns
  !--------------------------------------------------------------------- 
  ! LHC edge variables
  !--------------------------------------------------------------------- 
  Real(dp) :: LHC_observ(50), mSle(6), mSu(6), mSd(6)

  c_phot = "photon"
  c_grav = "~Gravitino"
  c_lm(1) = "e-"
  c_lp(1) = "e+"
  c_lm(2) = "mu-"
  c_lp(2) = "mu+"
  c_lm(3) = "tau-"
  c_lp(3) = "tau+"
  If (GenerationMixing) Then
   c_nu(1) = "nu_1"
   c_nu(2) = "nu_2"
   c_nu(3) = "nu_3"
  Else
   c_nu(1) = "nu_e"
   c_nu(2) = "nu_mu"
   c_nu(3) = "nu_tau"
  End If
  c_d(1) = "d"
  c_d(2) = "s"
  c_d(3) = "b"
  c_u(1) = "u"
  c_u(2) = "c"
  c_u(3) = "t"
  c_gl = "~g"
  c_cp(1) = "chi_1+"
  c_cp(2) = "chi_2+"
  c_cm(1) = "chi_1-"
  c_cm(2) = "chi_2-"
  c_c0(1) = "chi_10"
  c_c0(2) = "chi_20"
  c_c0(3) = "chi_30"
  c_c0(4) = "chi_40"
  id_n(1) = 1000022
  id_n(2) = 1000023
  id_n(3) = 1000025
  id_n(4) = 1000035

  c_sp(1) = "H+"
  c_sm(1) = "H-"
  id_sp = 0
  id_sp(1) = 37
  id_sm = - id_sp  

  n_sd = 6
  n_su = 6
 
  delta_n_rp = 0
  delta_c_rp = 0

  id_c(1) = 1000024
  id_c(2) = 1000037

  If (HighScaleModel.Eq."NMSSM") Then
   c_c0(5) = "chi_50"
   id_n(5) = 1000045
   c_s0(1) = "S0_1"
   c_s0(2) = "S0_2"
   c_s0(3) = "S0_3"
   id_S0(1) = 25
   id_S0(2) = 35
   id_S0(3) = 45
   c_P0(1) = "P0_1"
   c_P0(2) = "P0_2"
   id_P0(1) = 36
   id_P0(2) = 46
   n_n = 5
   n_c = 2
   n_s0 = 3
   n_p0 = 2 
   n_spm = 1 
   n_sl = 6
   n_sn = 3
   n_min = 1   ! first relevant neutralino index
   c_min = 1   ! first relevant chargino index

  Else If (HighScaleModel.Eq."RPexplicit") Then
   id_c(4:5) = id_c(1:2)
   id_c(1:3) = - id_l
   c_cp(4:5) = c_cp(1:2)
   c_cp(1:3) = c_lp
   c_cm(4:5) = c_cm(1:2)
   c_cm(1:3) = c_lm

   id_n(4:7) = id_n(1:4)
   id_n(1:3) = id_nu
   c_c0(4:7) = c_c0(1:4)
   c_c0(1) = "nu_1"
   c_c0(2) = "nu_2"
   c_c0(3) = "nu_3"
!   l_CS = .False.
!   L_BR = .False.
   n_n = 4
   n_c = 2
   n_s0 = 5
   n_p0 = 4 
   n_spm = 7 
   n_sl = 0
   n_sn = 0 
   delta_n_rp = 3
   delta_c_rp = 3
   Inquire(file="oldRP.dat",exist=file_exists)
   If (file_exists) Then
    Open(65,file="oldRP.dat")
    Read(65,*) use_new_RP
    Close(65)
   Else
    use_new_RP = .True.
   End If

   If (use_new_RP) Then
    mat5R = Abs(RP05)
    mat5R(1,:) = 0._dp ! this is the Goldstone boson
    Do i1=1,4
     c_P0(i1) = "P0_"//Bu(i1)
     MaxCont = Maxval(mat5R)
     Call FindPosition(5, mat5R, MaxCont, ii, jj)
     Select Case(jj)
     Case(3)
      id_p0(ii-1) = 2000012   ! Im(snu_e)
     Case(4)
      id_p0(ii-1) = 2000014   ! Im(snu_mu)
     Case(5)
      id_p0(ii-1) = 2000016   ! Im(snu_tau)
     Case default
      id_p0(ii-1) = 36        ! A_0
     End Select
     mat5R(ii,:) = 0._dp
     mat5R(:,jj) = 0._dp
    End Do

    i_zaehl = 1
    mat5R = Abs(RS05)
    Do i1=1,5
     c_s0(i1) = "S0_"//Bu(i1)
     MaxCont = Maxval(mat5R)
     Call FindPosition(5, mat5R, MaxCont, ii, jj)
     Select Case(jj)
     Case(3)
      id_s0(ii) = 1000012   ! Re(snu_e)
     Case(4)
      id_s0(ii) = 1000014   ! Re(snu_mu)
     Case(5)
      id_s0(ii) = 1000016   ! Re(snu_tau)
     Case default
      id_s0(ii) = 15 + i_zaehl*10
      i_zaehl = i_zaehl + 1
     End Select
     mat5R(ii,:) = 0._dp
     mat5R(:,jj) = 0._dp
    End Do

    i_zaehl = 1 
    mat8R = Abs(RSpm8)
    mat8R(1,:) = 0._dp ! this is the Goldstone boson
    Do i1=1,7
     c_sp(i1) = "S^+_"//Bu(i1)
     c_sm(i1) = "S^-_"//Bu(i1)
     MaxCont = Maxval(mat8R)
     Call FindPosition(8, mat8R, MaxCont, ii, jj)
     Select Case(jj)
     Case(1,2)
      id_Sp(ii-1) = 37           ! H^+
     Case(3)
      id_sp(ii-1) = -1000011     ! ~e_L
     Case(4)
      id_sp(ii-1) = -1000013     ! ~mu_L
     Case(6)
      id_sp(ii-1) = -2000011     ! ~e_R
     Case(7)
      id_sp(ii-1) = -2000013     ! ~mu_R
     Case default              ! stau_(i_zaehl)
      id_sp(ii-1) = -(1000000 * i_zaehl + 15)
      i_zaehl = i_zaehl + 1
     End Select
     mat8R(ii,:) = 0._dp
     mat8R(:,jj) = 0._dp
    End Do
    
   Else

    Do i1=1,4
     c_s0(i1) = "S0_"//Bu(i1)
     c_P0(i1) = "P0_"//Bu(i1)
     id_s0(i1) = 15 + 10 * i1
     id_p0(i1) = 26 + 10 * i1
    End Do
    c_s0(5) = "S0_5"
    id_s0(5) = 65
    Do i1=1,7
     c_sp(i1) = "S^+_"//Bu(i1)
     c_sm(i1) = "S^-_"//Bu(i1)
     id_sp(i1) = 27 + 10 * i1
    End Do
   End If

   id_sm = - id_sp
   n_min = 4   ! first relevant neutralino index
   c_min = 4   ! first relevant chargino index

  Else If (HighScaleModel.Eq."NURRP1") Then
   id_c(4:5) = id_c(1:2)
   id_c(1:3) = - id_l
   c_cp(4:5) = c_cp(1:2)
   c_cp(1:3) = c_lp
   c_cm(4:5) = c_cm(1:2)
   c_cm(1:3) = c_lm

   id_n(4:7) = id_n(1:4)
   id_n(1:3) = id_nu
   c_c0(4:7) = c_c0(1:4)
   c_c0(1) = "nu_1"
   c_c0(2) = "nu_2"
   c_c0(3) = "nu_3"

   n_n = 5
   n_c = 2
   n_s0 = 6
   n_p0 = 5 
   n_spm = 7 
   n_sl = 0
   n_sn = 0 
   delta_n_rp = 3
   delta_c_rp = 3

   i_zaehl = 1 
   mat8R = Abs(RSpm8)
   mat8R(1,:) = 0._dp ! this is the Goldstone boson
   Do i1=1,7
    c_sp(i1) = "S^+_"//Bu(i1)
    c_sm(i1) = "S^-_"//Bu(i1)
    MaxCont = Maxval(mat8R)
    Call FindPosition(8, mat8R, MaxCont, ii, jj)
    Select Case(jj)
    Case(1,2)
     id_Sp(ii) = 37           ! H^+
    Case(3)
     id_sp(i1) = -1000011     ! ~e_L
    Case(4)
     id_sp(i1) = -1000013     ! ~mu_L
    Case(6)
     id_sp(i1) = -2000011     ! ~e_R
    Case(7)
     id_sp(i1) = -2000013     ! ~mu_R
    Case default              ! stau_(i_zaehl)
     id_sp(i1) = -(1000000 * i_zaehl + 15)
     i_zaehl = i_zaehl + 1
    End Select
    mat8R(ii,:) = 0._dp
    mat8R(:,jj) = 0._dp
   End Do

   id_sm = - id_sp
   n_min = 4   ! first relevant neutralino index
   c_min = 4   ! first relevant chargino index

  Else
   c_s0(1) = "h0"
   c_s0(2) = "H0"
   id_S0(1) = 25
   id_S0(2) = 35
   c_P0(1) = "A0"
   id_P0(1) = 36
   n_n = 4
   n_c = 2
   n_s0 = 2
   n_p0 = 1 
   n_spm = 1 
   n_sl = 6
   n_sn = 3
   n_min = 1   ! first relevant neutralino index
   c_min = 1   ! first relevant chargino index
  End If

  Q = Sqrt( GetRenormalizationScale() )

  Call Date_and_time(datum,zeit)
  If (l_open) Then
   If (Present(f_name)) Then
    Open(io_L,file=Trim(f_name),status="unknown")
   Else
    Open(io_L,file="SPheno.spc",status="unknown")
   End If
   l_open = .False.
  End If
  !--------------------------------------------------------
  ! General information
  !--------------------------------------------------------
  ! Les Houches standard
  !-----------------------
   Write(io_L,100) "# SUSY Les Houches Accord 2 - MSSM spectrum + Decays"
   Write(io_L,100) "# SPheno "//version
   Write(io_L,100) &
     & "# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101"
   Write(io_L,100) "# in case of problems send email to porod@physik.uni-wuerzburg.de"
   Write(io_L,100) "# Created: "//Datum(7:8)//"."//Datum(5:6)//"."//Datum(1:4) &
     & //",  "//Zeit(1:2)//":"//Zeit(3:4)
   Write(io_L,100) "Block SPINFO         # Program information"
   Write(io_L,100) "     1   SPheno      # spectrum calculator"
   Write(io_L,100) "     2   "//version//"    # version number"
   !-----------------------------------------------
   ! check if somewhere a problem has had happened
   !-----------------------------------------------
   Call GetError(i_errors)
   !--------------------------------------
   ! a numerical problem might have happen
   !--------------------------------------
   If ((i_errors(1)+i_errors(3)+i_errors(5)+i_errors(7)+i_errors(8) &
     & + i_errors(10) + i_errors(12)+ Sum(i_errors(14:19))).Gt.0)   &
     & Write(io_L,100) &
 & "     3               # potential numerical problem, check file Messages.out"
   If (in_kont(1).Eq.1) Write(io_L,99) 3, &
    & "alpha(0) and alpha(mZ) have both been specified without check for"// &
    &  " consistency"
   If (in_kont(2).Eq.1) Write(io_L,99) 3, &
    & "redundant specification in Higgs sector"
   If (kont.Ne.0)   Write(io_L,100)  &
     "     4               # internal problem, see Messages.out for infos"
   Write(io_L,100) "#"
   Write(io_L,100) "Block SPhenoINFO     # SPheno specific information"
   If (TwoLoopRGE) Then
    Write(io_L,100) "    1      2         # using 2-loop RGEs"
   Else 
    Write(io_L,100) "    1      1         # using 1-loop RGEs"
   End If
   If (YukScen.Eq.1) Then
    Write(io_L,100) &
     &"    2      1         # using running masses for boundary conditions at mZ"
   Else
    Write(io_L,100) &
      &"    2      2         # using pole masses for boundary conditions at mZ"
   End If
   !------------------------------
   ! SPheno.out
   !------------------------------
  Write(io,123) "SPheno output file"
  Write(io,123) "Version "//version//" ,  "//"created: "//Datum(7:8)//"."// &
     & Datum(5:6)//"."//Datum(1:4)//",  "//Zeit(1:2)//":"//Zeit(3:4)
  Write(io,*) " "
123 Format(a)

  If (kont.Ne.0) Then
   Write(io,*) "There has been a problem during the run."
   Write(io,*) "Please check the file Messages.out for further information."
   Write(*,*) "There has been a problem during the run."
   Write(*,*) "Please check the file Messages.out for further information."
   Write(io,*) " "
   Return
  End If

  Write(io,*) " "
  If (YukScen.Eq.1) Then
   Write(io,*) &
     & "Running masses have been used for the boundary conditions at mZ" 
  Else If (YukScen.Eq.2) Then
   Write(io,*) "Pole masses have been used for the boundary conditions at mZ" 
  End If
  Write(io,*) "  Using Yukawa scheme :",GetYukawaScheme()
  Write(io,*) " "


   !--------------------------------------
   ! model information
   !--------------------------------------
   If (HighScaleModel.Eq."mSugra") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (i_cpv.Gt.0) Write(io_L,110) 5,i_cpv,"switching on CP violation"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# cos(phase_mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    Write(io_L,100) "#"
    If ((Aimag(mu).Ne.0._dp).Or.(Aimag(AoY_l_0(1,1)).Ne.0._dp)) Then
     Write(io_L,100) "Block IMMINPAR  # Input parameters, imaginary part"
     Write(io_L,101) 4, Aimag(phase_mu),"# sin(phase_mu)"
     Write(io_L,101) 5, Aimag(AoY_l_0(1,1)),"# Im(A0)"
    End If

    Write(io,*) "mSugra input at the GUT scale", m_GUT
    If (Aimag(Mi_0(1)).Eq.0._dp) Then
     Write(io,*) "M_1/2     : ",Real(Mi_0(1))
    Else
     Write(io,*) "M_1/2     : ",Mi_0(1)
    End If
    Write(io,*) "M_0       : ",Real(Sqrt(M2_D_0(1,1)))
    If (Aimag(AoY_d_0(1,1)).Eq.0._dp) Then
     Write(io,*) "A_0       : ",Real(AoY_d_0(1,1))
    Else
     Write(io,*) "A_0       : ",AoY_d_0(1,1)
    End If
    Write(io,*) "tan(beta) at m_Z : ", tanb_mZ
    If (Aimag(phase_mu).Eq.0._dp) Then
     Write(io,*) "phase(mu) : ", Real(phase_mu)
    Else
     Write(io,*) "phase(mu) : ", phase_mu
    End If

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
   Else If ((HighScaleModel(1:5).Eq."SUGRA").Or.   &
          & (HighScaleModel(1:5).Eq."SEESA")) Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,100) "    3    2    # mSUGRA model + SU(5)"
    Else If (HighScaleModel.Eq."SUGRA_NuR1") Then
     Write(io_L,100) "    3    3    # mSUGRA model + nu_R at a common scale"
     Write(io_L,100) "Block MnuR  # mass scale of the right handed neutrinos"
     Write(io_L,101) 1,MnuR(1),"# m_nu_R      "     
    Else If (HighScaleModel.Eq."SUGRA_NuR") Then
     Y_nu_0 = Transpose(Y_nu_0) ! in the RGEs the transposed Yukawas are used
     Write(io_L,100) "    3    4    # mSUGRA model + three  nu_R, ~nu_R"

     Call WriteMatrixBlockC(io_L,3,Y_nu_0,m_GUT &
                         & ,"Ynu0","GUT scale","Y_(nu,")

    Else If (HighScaleModel.Eq."SEESAW_II") Then
     Write(io_L,100) "    3    5    # mSUGRA model + Higgs triplett"

     Call WriteMatrixBlockC(io_L,3,Y_T_0,m_GUT &                ! symmetric
                         & ,"YT0","GUT scale","Y_(T,", .True.)  ! matrix

     Write(io_L,106) "Block Higgs3 Q=",m_GUT,"# (GUT scale)"
     Write(io_L,101) 1,M_H3(1),"# m_H3      "
     Write(io_L,101) 2,Real(lam12_0(1),dp),"# Re(lambda_1)"
     If (Aimag(lam12_0(1)).Ne.0._dp) &
         Write(io_L,101) 3,Aimag(lam12_0(1)),"# Im(lambda_1)"
     Write(io_L,101) 4,Real(lam12_0(2),dp),"# Re(lambda_2)"
     If (Aimag(lam12_0(2)).Ne.0._dp) &
         Write(io_L,101) 5,Aimag(lam12_0(2)),"# Im(lambda_2)"
     If (Fifteen_plet) Then
      Write(io_L,100) "    6    1               # using RGEs for 15-plet"
     Else
      Write(io_L,100) "    6    0               # using RGEs for 3-plet"
     End If
! Florian SARAH, Seesaw II+III
    Else If (HighScaleModel.Eq."SEESAW_III_3G") Then
     Write(io_L,100) "    3    10    # mSUGRA model + 3 24-plets"

     Call WriteMatrixBlockC(io_L,3,Yb3_h24_gut,m_GUT &
                         & ,"YB","GUT scale","Y_(b,")

     Write(io_L,106) "Block Higgs3 Q=",m_GUT,"# (GUT scale)"
     Write(io_L,101) 1,Real(MWM3_gut(1,1),dp),"# MWM      "
     Write(io_L,106) "Block Higgs3 Q=",Abs(MWM30(1,1)),"# (Triplet Scale 1 )"
     Write(io_L,101) 1,Abs(MWM30(1,1)),"# MWM   "
     Write(io_L,101) 2,Abs(MGM3(1,1)),"# MGM   "
     Write(io_L,101) 3,Abs(MBM3(1,1)),"# MBM   "                    
     Write(io_L,101) 4,Abs(MXM3(1,1)),"# MXM   "          
     Write(io_L,106) "Block Higgs3 Q=",Abs(MWM30(2,2)),"# (Triplet Scale 2 )"
     Write(io_L,101) 1,Abs(MWM30(2,2)),"# MWM   "
     Write(io_L,101) 2,Abs(MGM3(2,2)),"# MGM   "
     Write(io_L,101) 3,Abs(MBM3(2,2)),"# MBM   "                    
     Write(io_L,101) 4,Abs(MXM3(2,2)),"# MXM   "          
     Write(io_L,106) "Block Higgs3 Q=",Abs(MWM30(3,3)),"# (Triplet Scale 3 )"
     Write(io_L,101) 1,Abs(MWM30(3,3)),"# MWM   "
     Write(io_L,101) 2,Abs(MGM3(3,3)),"# MGM   "
     Write(io_L,101) 3,Abs(MBM3(3,3)),"# MBM   "                    
     Write(io_L,101) 4,Abs(MXM3(3,3)),"# MXM   "          

     Call WriteMatrixBlockC(io_L,3,Yb30_h24(1,:,:),Abs(MWM30(1,1)) &
                         & ,"YB","Triplet scale 1","Y_(b,")
     Call WriteMatrixBlockC(io_L,3,Yw30_h24(1,:,:),Abs(MWM30(1,1)) &
                         & ,"YW","Triplet scale 1","Y_(w,")
     Call WriteMatrixBlockC(io_L,3,Yx30_h24(1,:,:),Abs(MWM30(1,1)) &
                         & ,"YX","Triplet scale 1","Y_(x,")

     Call WriteMatrixBlockC(io_L,3,Yb30_h24(2,:,:),Abs(MWM30(2,2)) &
                         & ,"YB","Triplet scale 2","Y_(b,")
     Call WriteMatrixBlockC(io_L,3,Yw30_h24(2,:,:),Abs(MWM30(2,2)) &
                         & ,"YW","Triplet scale 2","Y_(w,")
     Call WriteMatrixBlockC(io_L,3,Yx30_h24(2,:,:),Abs(MWM30(2,2)) &
                         & ,"YX","Triplet scale 2","Y_(x,")

     Call WriteMatrixBlockC(io_L,3,Yb30_h24(3,:,:),Abs(MWM30(3,3)) &
                         & ,"YB","Triplet scale 3","Y_(b,")
     Call WriteMatrixBlockC(io_L,3,Yw30_h24(3,:,:),Abs(MWM30(3,3)) &
                         & ,"YW","Triplet scale 3","Y_(w,")
     Call WriteMatrixBlockC(io_L,3,Yx30_h24(3,:,:),Abs(MWM30(3,3)) &
                         & ,"YX","Triplet scale 3","Y_(x,")


    Else If (HighScaleModel.Eq."SEESAW_II_SARAH") Then
     Write(io_L,100) "    3    11    # mSUGRA model + 1 15-plet"

     Call WriteMatrixBlockC(io_L,3,YT_h15_gut,m_GUT &          ! symmetric
                         & ,"YT","GUT scale","Y_(T,", .True.)  ! matrix

     Write(io_L,106) "Block Higgs3 Q=",m_GUT,"# (GUT scale)"
     Write(io_L,101) 1,Abs(MTM_gut),"# MWM      "
     Write(io_L,106) "Block Higgs3 Q=",Abs(MTM0),"# (Triplet Scale)"
     Write(io_L,101) 1,Abs(MTM),"# MTM   "
     Write(io_L,101) 2,Abs(MZM),"# MZM   "
     Write(io_L,101) 3,Abs(MSM),"# MSM   "                    
     
     Call WriteMatrixBlockC(io_L,3,YT_h15,Abs(MTM0) &            ! symmetric
                         & ,"YT","triplet scale","Y_(T,", .True.)! matrix
     Call WriteMatrixBlockC(io_L,3,YS_h15,Abs(MTM0) &            ! symmetric
                         & ,"YS","triplet scale","Y_(S,", .True.)! matrix
     Call WriteMatrixBlockC(io_L,3,YZ_h15,Abs(MTM0) &
                         & ,"YZ","triplet scale","Y_(Z,")

     Write(io_L,106) "Block Lambda Q=",Abs(MTM0),"# (Triplet Scale)"
     Write(io_L,101) 1,Real(Lambda10,dp),"# Lambda1"
     Write(io_L,101) 2,Real(Lambda20,dp),"# Lambda2"
     If ((Aimag(Lambda10).Ne.0._dp).Or.(Aimag(Lambda20).Ne.0._dp)) Then
      Write(io_L,106) "Block IMLambda Q=",Abs(MTM0),"# (Triplet Scale)"
      Write(io_L,101) 1,Aimag(Lambda10),"# Lambda1"
      Write(io_L,101) 2,Aimag(Lambda20),"# Lambda2"
     End If
! Florian SARAH, Seesaw II+III      
    End If
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,101) 7,Max(M_SO_10,m_GUT),"# SO(10) scale"
     Write(io_L,101) 8,D_SO_10,"# D-terms at SO(10) scale"
    End If
    Write(io_L,100) "#"

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"

    Write(io,*) "mSugra input at the GUT scale", m_GUT
    If (Aimag(Mi_0(1)).Eq.0._dp) Then
     Write(io,*) "M_1/2     : ",Real(Mi_0(1))
    Else
     Write(io,*) "M_1/2     : ",Mi_0(1)
    End If
    Write(io,*) "M_0       : ",Real(Sqrt(M2_D_0(1,1)))
    If (Aimag(AoY_d_0(1,1)).Eq.0._dp) Then
     Write(io,*) "A_0       : ",Real(AoY_d_0(1,1))
    Else
     Write(io,*) "A_0       : ",AoY_d_0(1,1)
    End If
    Write(io,*) "tan(beta) at m_Z : ", tanb_mZ
    If (Aimag(phase_mu).Eq.0._dp) Then
     Write(io,*) "phase(mu) : ", Real(phase_mu)
    Else
     Write(io,*) "phase(mu) : ", phase_mu
    End If
    
   Else If (HighScaleModel.Eq."GMSB") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    2    # mGMSB model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Lambda,"# Lambda, scale of SUSY breaking"
    Write(io_L,101) 2,MlambdaS,"# M_M, messenger mass scale"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(n5plets+3*n10plets,dp) , &
                      & "# N_5 number of effective five-plets"
    Write(io_L,101) 6,grav_fac,"# c_grav, Gravitino mass factor"

    Write(io,*) "Lambda    : ",Lambda
    Write(io,*) "M_M       : ",MlambdaS
    Write(io,*) "n_5, n_10 : ",n5plets,n10plets
    If (Aimag(AoY_l_0(1,1)).Eq.0._dp) Then
     Write(io,*) "A_0       : ",Real(AoY_l_0(1,1))
    Else
     Write(io,*) "A_0       : ",AoY_l_0(1,1)
    End If
    Write(io,*) "tan(beta) : ", tanb
    If (Aimag(phase_mu).Eq.0._dp) Then
     Write(io,*) "phase(mu) : ", Real(phase_mu)
    Else
     Write(io,*) "phase(mu) : ", phase_mu
    End If
    Write(io,*) " "
    Write(io,*) "gravitino mass",m32," eV"

    Write(io_L,106) "Block gauge Q=",MlambdaS,"# (GMSB scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
   Else If (HighScaleModel.Eq."AMSB") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    3    # mAMSB model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,m0_amsb,"# M_0" 
    Write(io_L,101) 2,m_32,"# m_3/2, gravitino mass"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"

    Write(io,*) "AMSB input at the GUT scale", m_GUT
    Write(io,*) "M_3/2     : ",m_32
    Write(io,*) "M_0       : ",m0_amsb
    Write(io,*) "tan(beta) : ", tanb
    If (Aimag(phase_mu).Eq.0._dp) Then
     Write(io,*) "phase(mu) : ", Real(phase_mu)
    Else
     Write(io,*) "phase(mu) : ", phase_mu
    End If

   Else If (HighScaleModel.Eq."NMSSM") Then
    non_minimal = .True.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    3    1    # NMSSM model"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
    If (phase_mu.Ne.0._dp) Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"

    Write(io,*) "General NMSSM, tree level masses"
    Write(io,*) "tan(beta) at m_Z : ", tanb_mZ
    Write(io,*) "lambda           : ",Real(h0,dp)
    Write(io,*) "kappa            : ",Real(2._dp*lam,dp)
    Write(io,*) "A_lambda         : ",Real(Ao_h0,dp)
    Write(io,*) "A_kappa          : ",Real(Ao_lam,dp)
    
   Else If (HighScaleModel.Eq."RPexplicit") Then
    non_minimal = .True.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    4    1    # MSSM with explicit R-parity violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
    If (phase_mu.Ne.0._dp) Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"

    Write(io,*) "MSSM, explicit R-parity violation, tree level masses"
    Write(io,*) "tan(beta) at m_Z : ", tanb_mZ
    
   Else
    non_minimal = .True.
    Write(io_L,100) "# Either the general MSSM or a model has been used"
    Write(io_L,100) &
      & "# which has not yet been implemented in the LesHouches standard"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
   End If 

   !---------------------------------------------------
   ! parameters + masses for SPheno.spc
   !---------------------------------------------------
   Write(io_L,100) "Block SMINPUTS  # SM parameters"
   Write(io_L,102) 1, 1._dp / alpha_MSbar(mZ, mW),"# alpha_em^-1(MZ)^MSbar"
   Write(io_L,102) 2,G_F,"# G_mu [GeV^-2]"
   Write(io_L,102) 3,alphaS_MZ,"# alpha_s(MZ)^MSbar"
   Write(io_L,102) 4,mZ,"# m_Z(pole)"
   Write(io_L,102) 5,mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) 6,mf_u(3),"# m_t(pole)"
   Write(io_L,102) 7,mf_l(3),"# m_tau(pole)"

   !----------------------------------------------------------------
   ! in the case of GenerationMixing all parameters and mixings
   ! are given in the SuperCKM basis for squarks
   !----------------------------------------------------------------
   If (GenerationMixing) Then
    Write(io_L,100) "Block VCKMIN  # CKM matrix, Wolfenstein parameterization"
    Write(io_L,102) 1, lam_wolf,"# lambda"
    Write(io_L,102) 2,A_wolf,"# A"
    Write(io_L,102) 3,rho_wolf,"# rho bar"
    Write(io_L,102) 4,eta_wolf,"# eta bar"
    Call Switch_to_superCKM(Y_d, Y_u, A_d, A_u, M2_D, M2_Q, M2_U         &
              &, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm, M2U_sckm, .False. &
              &, RSdown, RSup, Rdsq_ckm, RUsq_ckm, CKM_Q, Yd, Yu )

    If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
     Call Switch_to_superPMNS(Y_l, MnuL5, A_l, M2_E, M2_L, Al_pmns, M2E_pmns  &
        &, M2L_pmns, .False., RSlepton, RSneut, Rsl_pmns, RSn_pmns, PMNS_Q, Yl)

    Else
     Call Switch_to_superPMNS(Y_l, id3C, A_l, M2_E, M2_L, Al_pmns, M2E_pmns  &
        &, M2L_pmns, .False., RSlepton, RSneut, Rsl_pmns, RSn_pmns, PMNS_Q, Yl)

    End If

   Else ! .non.GenerationMixing

    Do i1=1,3
     Yu(i1) = Real(Y_u(i1,i1),dp)
     Yd(i1) = Real(Y_d(i1,i1),dp)
     Yl(i1) = Real(Y_l(i1,i1),dp)
    End Do
    Al_pmns = A_l
    Ad_sckm = A_d
    Au_sckm = A_u

    M2D_SCKM = M2_D
    M2U_SCKM = M2_U
    M2Q_SCKM = M2_Q
    M2E_pmns = M2_E
    M2L_pmns = M2_L

    RUsq_ckm = RSup
    RDsq_ckm = RSdown

    RSn_pmns = RSneut
    RSl_pmns = RSlepton

   End If

   If (non_minimal) Then
    Write(io_L,100) "Block EXTPAR  # "
    Write(io_L,104) 0,Q , "# scale Q where the parameters below are defined"
    Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
    Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
    Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
    If (At_save.ne.0._dp) then
     Write(io_L,104) 11,Real(At_save,dp), "# A_t"
    Else
     Write(io_L,104) 11,Real(Au_sckm(3,3)/yu(3),dp), "# A_t"
    End If
    If (Ab_save.ne.0._dp) then
     Write(io_L,104) 12,Real(Ab_save,dp), "# A_b"
    Else
     Write(io_L,104) 12,Real(Ad_sckm(3,3)/yd(3),dp), "# A_b"
    End If
    If (Atau_save.ne.0._dp) then
     Write(io_L,104) 13,Real(Atau_save,dp), "# A_l"
    Else
     Write(io_L,104) 13,Real(Al_pmns(3,3)/yl(3),dp), "# A_l"
    End If
    Write(io_L,104) 23,Real(mu ,dp), "# mu "
    Write(io_L,104) 31,Sqrt(Real(M2L_pmns(1,1),dp)),"# M_(L,11)"
    Write(io_L,104) 32,Sqrt(Real(M2L_pmns(2,2),dp)),"# M_(L,22)"
    Write(io_L,104) 33,Sqrt(Real(M2L_pmns(3,3),dp)),"# M_(L,33)"
    Write(io_L,104) 34,Sqrt(Real(M2E_pmns(1,1),dp)),"# M_(E,11)"
    Write(io_L,104) 35,Sqrt(Real(M2E_pmns(2,2),dp)),"# M_(E,22)"
    Write(io_L,104) 36,Sqrt(Real(M2E_pmns(3,3),dp)),"# M_(E,33)"
    Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
    Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
    Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
    Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
    Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
    Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
    Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
    Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
    Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"
    If (HighScaleModel.Eq."NMSSM") Then
     Write(io_L,104) 61,Real(h0,dp),"# lambda"
     Write(io_L,104) 62,0.5_dp*Real(lam,dp),"# kappa"
     Write(io_L,104) 63,Real(ao_h0,dp),"# A_lambda"
     Write(io_L,104) 64,Real(Ao_lam,dp),"# A_kappa"
     Write(io_L,104) 65,oosqrt2*Real(vP*h0,dp),"# mu_eff = lambda * <S> "
    End If
   End If
      
! couplings
  Write(io_L,106) "Block gauge Q=",Q,"# (SUSY scale)"
  Write(io_L,104) 1,gauge(1),"# g'(Q)^DRbar"
  Write(io_L,104) 2,gauge(2),"# g(Q)^DRbar"
  Write(io_L,104) 3,gauge(3),"# g3(Q)^DRbar"

  Write(io_L,106) "Block Yu Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yu(1), "# Y_u(Q)^DRbar"
  Write(io_L,107) 2,2,Yu(2), "# Y_c(Q)^DRbar"
  Write(io_L,107) 3,3,Yu(3), "# Y_t(Q)^DRbar"

  Write(io_L,106) "Block Yd Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yd(1), "# Y_d(Q)^DRbar"
  Write(io_L,107) 2,2,Yd(2), "# Y_s(Q)^DRbar"
  Write(io_L,107) 3,3,Yd(3), "# Y_b(Q)^DRbar"

  Write(io_L,106) "Block Yl Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yl(1), "# Y_e(Q)^DRbar"
  Write(io_L,107) 2,2,Yl(2), "# Y_mu(Q)^DRbar"
  Write(io_L,107) 3,3,Yl(3), "# Y_tau(Q)^DRbar"

  If (GenerationMixing) Then 
                             
   Call WriteMatrixBlockC(io_L,3,CKM_Q,Q &
                         & ,"VCKM","V_CKM at the SUSY scale","V_(")

   Call WriteMatrixBlockC(io_L,3,PMNS_Q,Q &
                         & ,"VPMNS","V_PMNS at the SUSY scale","V_(")

  End If ! generationmixing


  If (GenerationMixing) Then

   Call WriteMatrixBlockC(io_L,3,Au_sckm,Q,"Tu","SUSY scale","T_(u,",tr=.True.)

  Else 
   Write(io_L,106) "Block Au Q=",Q,"# (SUSY scale)"
   If (Abs(y_u(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Au_sckm(1,1)/y_u(1,1),dp), "# A_u(Q)^DRbar"
   If (Abs(y_u(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Au_sckm(2,2)/y_u(2,2),dp), "# A_c(Q)^DRbar"
   If (Abs(y_u(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t(Q)^DRbar"
   If (Maxval(Abs(Aimag(Au_sckm))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAu Q=",Q,"# (SUSY scale)"
    If (Abs(y_u(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Au_sckm(1,1)/y_u(1,1)), "# Im(A_u)(Q)^DRbar"
    If (Abs(y_u(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Au_sckm(2,2)/y_u(2,2)), "# Im(A_c)(Q)^DRbar"
    If (Abs(y_u(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Au_sckm(3,3)/y_u(3,3)), "# Im(A_t)(Q)^DRbar"
   End If
  End If

  If (GenerationMixing) Then

   Call WriteMatrixBlockC(io_L,3,Ad_sckm,Q,"Td","SUSY scale","T_(d,",tr=.True.)

  Else 
   Write(io_L,106) "Block Ad Q=",Q,"# (SUSY scale)"
   If (Abs(y_d(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Ad_sckm(1,1)/y_d(1,1),dp), "# A_d(Q)^DRbar"
   If (Abs(y_d(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Ad_sckm(2,2)/y_d(2,2),dp), "# A_s(Q)^DRbar"
   If (Abs(y_d(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b(Q)^DRbar"
   If (Maxval(Abs(Aimag(Ad_sckm))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAd Q=",Q,"# (SUSY scale)"
    If (Abs(y_d(1,1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Ad_sckm(1,1)/y_d(1,1)), "# Im(A_d)(Q)^DRbar"
    If (Abs(y_d(2,2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Ad_sckm(2,2)/y_d(2,2)), "# Im(A_s)(Q)^DRbar"
    If (Abs(y_d(3,3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Ad_sckm(3,3)/y_d(3,3)), "# Im(A_b)(Q)^DRbar"
   End If
  End If

  ierr = 0
  If (GenerationMixing) Then
   !-------------------------------------------
   ! check if any off-diagonal term is non-zero
   !-------------------------------------------
   Do i1=1,3
    Do i2=1,3
     If ((i1.Ne.i2).And.(Abs(Al_pmns(i2,i1)).Ne.0._dp)) ierr = ierr + 1
    End Do
   End Do
  End If

  If (ierr.Ne.0) Then

   Call WriteMatrixBlockC(io_L,3,Al_pmns,Q,"Te","SUSY scale","T_(l,") ! ,tr=.True.) checken

  Else 
   Write(io_L,106) "Block Ae Q=",Q,"# (SUSY scale)"
   If (Abs(Yl(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Real(Al_pmns(1,1)/Yl(1),dp), "# A_e(Q)^DRbar"
   If (Abs(Yl(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Real(Al_pmns(2,2)/Yl(2),dp), "# A_mu(Q)^DRbar"
   If (Abs(Yl(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Real(Al_pmns(3,3)/Yl(3),dp), "# A_tau(Q)^DRbar"
   If (Maxval(Abs(Aimag(Al_pmns))).Gt.0._dp) Then
    Write(io_L,106) "Block IMAe Q=",Q,"# (SUSY scale)"
    If (Abs(Yl(1)).Gt.0._dp) &
        & Write(io_L,107) 1,1,Aimag(Al_pmns(1,1)/Yl(1)), "# Im(A_e)(Q)^DRbar"
    If (Abs(Yl(2)).Gt.0._dp) &
        & Write(io_L,107) 2,2,Aimag(Al_pmns(2,2)/Yl(2)), "# Im(A_mu)(Q)^DRbar"
    If (Abs(Yl(3)).Gt.0._dp) &
        & Write(io_L,107) 3,3,Aimag(Al_pmns(3,3)/Yl(3)), "# Im(A_tau)(Q)^DRbar"
   End If
  End If


  Write(io_L,106) "Block MSOFT Q=",Q,"# soft SUSY breaking masses at Q"
  Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
  Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
  Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
  Write(io_L,104) 21,M2_H(1),"# M^2_(H,d)"
  Write(io_L,104) 22,M2_H(2),"# M^2_(H,u)"

  Write(io_L,104) 31,Sqrt(Real(M2L_pmns(1,1),dp)),"# M_(L,11)"
  Write(io_L,104) 32,Sqrt(Real(M2L_pmns(2,2),dp)),"# M_(L,22)"
  Write(io_L,104) 33,Sqrt(Real(M2L_pmns(3,3),dp)),"# M_(L,33)"
  Write(io_L,104) 34,Sqrt(Real(M2E_pmns(1,1),dp)),"# M_(E,11)"
  Write(io_L,104) 35,Sqrt(Real(M2E_pmns(2,2),dp)),"# M_(E,22)"
  Write(io_L,104) 36,Sqrt(Real(M2E_pmns(3,3),dp)),"# M_(E,33)"
  Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
  Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
  Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
  Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
  Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
  Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
  Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
  Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
  Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"

  If (GenerationMixing) Then

   Call WriteMatrixBlockC(io_L,3,M2L_pmns,Q,"MSL2" &
           & ,"M^2_L soft SUSY breaking masses","M^2_(L,")

   Call WriteMatrixBlockC(io_L,3,M2E_pmns,Q,"MSE2" &
           & ,"M^2_E soft SUSY breaking masses","M^2_(E,")

   Call WriteMatrixBlockC(io_L,3,M2Q_SCKM,Q,"MSQ2" &
           & ,"M^2_Q soft SUSY breaking masses","M^2_(Q,")

   Call WriteMatrixBlockC(io_L,3,M2U_SCKM,Q,"MSU2" &
           & ,"M^2_U soft SUSY breaking masses","M^2_(U,")

   Call WriteMatrixBlockC(io_L,3,M2D_SCKM,Q,"MSD2" &
           & ,"M^2_D soft SUSY breaking masses","M^2_(D,")

  End If

  If (HighScaleModel.Eq."NMSSM") Then
   Write(io_L,104) 61,Real(h0,dp),"# lambda"
   Write(io_L,104) 62,0.5_dp*Real(lam,dp),"# kappa"
   Write(io_L,104) 63,Real(ao_h0,dp),"# A_lambda"
   Write(io_L,104) 64,Real(Ao_lam,dp),"# A_kappa"
   Write(io_L,104) 65,Real(oosqrt2*vP*h0,dp),"# mu_eff"

  Else If (HighScaleModel.Eq."RPexplicit") Then
   Write(io_L,106) "Block RVKAPPA Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(eps(1),dp),"# epsilon_1"
   Write(io_L,102) 2,Real(eps(2),dp),"# epsilon_2"
   Write(io_L,102) 3,Real(eps(3),dp),"# epsilon_3"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVKAPPA Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(eps(1)),"# Im(epsilon_1)"
    Write(io_L,102) 2,Aimag(eps(2)),"# Im(epsilon_2)"
    Write(io_L,102) 3,Aimag(eps(3)),"# Im(epsilon_3)"
   End If

   Write(io_L,106) "Block RVD Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(Beps(1),dp),"# Re( B_1 epsilon_1)"
   Write(io_L,102) 2,Real(Beps(2),dp),"# Re( B_2 epsilon_2)"
   Write(io_L,102) 3,Real(Beps(3),dp),"# Re( B_3 epsilon_3)"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVD Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(Beps(1)),"# Im( B_1 epsilon_1)"
    Write(io_L,102) 2,Aimag(Beps(2)),"# Im( B_2 epsilon_2)"
    Write(io_L,102) 3,Aimag(Beps(3)),"# Im( B_3 epsilon_3)"
   End If

   Write(io_L,106) "Block RVSNVEV Q=",Q,"# sneutrino vevs at Q"
   Write(io_L,102) 1,vevL(1),"# v_L_1"
   Write(io_L,102) 2,vevL(2),"# v_L_2"
   Write(io_L,102) 3,vevL(3),"# v_L_3"

   If (Maxval(Abs(Rp_lam)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDA Q=",Q,"# lambda_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lam(i1,i2,i3),dp) &
                     & ,"# lambda_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDA Q=",Q,"# Im(lambda_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lam(i1,i2,i3)) &
            & ,"# Im(lambda_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   If (Maxval(Abs(Rp_lamp)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDAP Q=",Q,"# lambda'_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lamp(i1,i2,i3),dp) &
                     & ,"# lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDAP Q=",Q,"# Im(lambda'_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lamp(i1,i2,i3)) &
            & ,"# Im(lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   Write(io_L,100) "Block SPhenoRP  # additional RP parameters"
   Write(io_L,102) 4,Lam_Ex(1),"# Lambda_1 = v_d epsilon_1 + mu v_L1"
   Write(io_L,102) 5,Lam_Ex(2),"# Lambda_2 = v_d epsilon_2 + mu v_L2"
   Write(io_L,102) 6,Lam_Ex(3),"# Lambda_3 = v_d epsilon_3 + mu v_L3"
   Write(io_L,102) 7,(mN7(3)**2-mN7(2)**2)*1.e18_dp,"# m^2_atm [eV^2]"
   Write(io_L,102) 8,(mN7(2)**2-mN7(1)**2)*1.e18_dp,"# m^2_sol [eV^2]"
   Write(io_L,102) 9,Abs(N7(3,6)/N7(3,7))**2,"# tan^2 theta_atm"
   Write(io_L,102) 10,Abs(N7(2,5)/N7(1,5))**2,"# tan^2 theta_sol"
   Write(io_L,102) 11,Abs(N7(3,5))**2,"# U_e3^2"
   Write(io_L,102) 15,vevSM(1),"# v_d"
   Write(io_L,102) 16,vevSM(2),"# v_u"
  End If


   Write(io_L,100) "Block MASS  # Mass spectrum"
   Write(io_L,100) "#   PDG code      mass          particle"
   Write(io_L,102) id_u(2),mf_u(2),"# m_c(m_c), MSbar"
   Write(io_L,102) id_d(3),mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) id_u(3),mf_u(3),"# m_t(pole)"
   Write(io_L,102) id_Z,mZ,"# m_Z(pole)"
   Write(io_L,102) id_W,mW,"# W+"
   If (HighScaleModel.Ne."RPexplicit") &
        & Write(io_L,102) id_l(3),mf_l(3),"# m_tau(pole)"
   If (HighScaleModel.Eq."NMSSM") Then
    Write(io_L,102) 25,mS03(1),"# leightest neutral scalar" 
    Write(io_L,102) 35,mS03(2),"# second neutral scalar" 
    Write(io_L,102) 45,mS03(3),"# third neutral scalar" 
    Write(io_L,102) 36,mP03(2),"# leighter pseudoscalar" 
    Write(io_L,102) 46,mP03(3),"# heavier pseudoscalar" 
    Write(io_L,102) 37,mSpm(2),"# H+"
   Else If (HighScaleModel.Eq."RPexplicit") Then
!    Write(io_L,102) 12,mN7(1),"# lightest neutrino"
!    Write(io_L,102) 14,mN7(2),"# second lightest neutrino"
!    Write(io_L,102) 16,mN7(3),"# heaviest neutrino"
    Write(io_L,102) id_s0(1),mS05(1),"# leightest neutral scalar" 
    Write(io_L,102) id_s0(2),mS05(2),"# 2nd neutral scalar" 
    Write(io_L,102) id_s0(3),mS05(3),"# 3rd neutral scalar" 
    Write(io_L,102) id_s0(4),mS05(4),"# 4th neutral scalar" 
    Write(io_L,102) id_s0(5),mS05(5),"# 5th neutral scalar" 
    Write(io_L,102) id_P0(1),mP05(2),"# leightest pseudoscalar" 
    Write(io_L,102) id_P0(2),mP05(3),"# 2nd pseudoscalar" 
    Write(io_L,102) id_P0(3),mP05(4),"# 3rd pseudoscalar" 
    Write(io_L,102) id_p0(4),mP05(5),"# 4th pseudoscalar"
    Write(io_L,102) Abs(id_sp(1)),mSpm8(2),"# leightest charged scalar" 
    Write(io_L,102) Abs(id_sp(2)),mSpm8(3),"# 2nd charged scalar" 
    Write(io_L,102) Abs(id_sp(3)),mSpm8(4),"# 3rd charged scalar" 
    Write(io_L,102) Abs(id_sp(4)),mSpm8(5),"# 4th charged scalar" 
    Write(io_L,102) Abs(id_sp(5)),mSpm8(6),"# 5th charged scalar" 
    Write(io_L,102) Abs(id_sp(6)),mSpm8(7),"# 6th charged scalar" 
    Write(io_L,102) Abs(id_sp(7)),mSpm8(8),"# 7th charged scalar" 
    
   Else
    Write(io_L,102) 25,mS0(1),"# h0" 
    Write(io_L,102) 35,mS0(2),"# H0" 
    Write(io_L,102) 36,mP0(2),"# A0" 
    Write(io_L,102) 37,mSpm(2),"# H+"
   End If
! squarks

  If (GenerationMixing) Then
   If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
     i_zaehl = 1
     mat6R = Abs(RDsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_sd(ii) = 1000001
       c_sd(ii) = "~d_L"
      Case(2)
       id_sd(ii) = 1000003
       c_sd(ii) = "~s_L-"
      Case(4)
       id_sd(ii) = 2000001
       c_sd(ii) = "~d_R"
      Case(5)
       id_sd(ii) = 2000003
       c_sd(ii) = "~s_R"
      Case default
       id_sd(ii) = 1000000 * i_zaehl + 5
       c_sd(ii) = "~b_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6 ! check ordering of sbottoms
      If (id_sd(ii).Eq.1000005)  id_check(1) = ii
      If (id_sd(ii).Eq.2000005)  id_check(2) = ii
     End Do
     If (id_check(1).Gt.id_check(2)) Then ! switch ordering
      ii = id_check(2)  ! the lighter one
      id_sd(ii) = 1000005
      c_sd(ii) = "~b_1"
      ii = id_check(1)  ! the heavier one
      id_sd(ii) = 2000005
      c_sd(ii) = "~b_2"
     End If
     Do ii=1,6
      Write(io_L,102) id_sd(ii),msdown(ii),"# "//Trim(c_sd(ii))
     End Do

     i_zaehl = 1
     mat6R = Abs(RUsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_su(ii) = 1000002
       c_su(ii) = "~u_L"
      Case(2)
       id_su(ii) = 1000004
       c_su(ii) = "~c_L-"
      Case(4)
       id_su(ii) = 2000002
       c_su(ii) = "~u_R"
      Case(5)
       id_su(ii) = 2000004
       c_su(ii) = "~c_R"
      Case default
       id_su(ii) = 1000000 * i_zaehl + 6
       c_su(ii) = "~t_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6 ! check ordering of stops
      If (id_su(ii).Eq.1000006)  id_check(1) = ii
      If (id_su(ii).Eq.2000006)  id_check(2) = ii
     End Do
     If (id_check(1).Gt.id_check(2)) Then ! switch ordering
      ii = id_check(2)  ! the lighter one
      id_su(ii) = 1000006
      c_su(ii) = "~t_1"
      ii = id_check(1)  ! the heavier one
      id_su(ii) = 2000006
      c_su(ii) = "~t_2"
     End If

     Do ii=1,6
      Write(io_L,102) id_su(ii),msup(ii),"# "//Trim(c_su(ii))
     End Do
      
    Else ! use mass ordering

     id_sd(1) = 1000001
     id_sd(2) = 1000003
     id_sd(3) = 1000005
     id_sd(4) = 2000001
     id_sd(5) = 2000003
     id_sd(6) = 2000005
     Do i1=1,6
      c_sd(i1) = "~d_"//Bu(i1)
      Write(io_L,102) id_sd(i1),msdown(i1),"# "//Trim(c_sd(i1))
     End Do

     id_su(1) = 1000002
     id_su(2) = 1000004
     id_su(3) = 1000006
     id_su(4) = 2000002
     id_su(5) = 2000004
     id_su(6) = 2000006
     Do i1=1,6
      c_su(i1) = "~u_"//Bu(i1)
      Write(io_L,102) id_su(i1),msup(i1),"# "//Trim(c_su(i1))
     End Do
    End If
  Else ! .not.GenerationMixing

   If (Abs(rsdown(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000001,msdown(1),"# ~d_L"
    Write(io_L,102) 2000001,msdown(2),"# ~d_R"
    id_sd(1) = 1000001
    id_sd(2) = 2000001
    c_sd(1) = "~d_L"
    c_sd(2) = "~d_R"
   Else
    Write(io_L,102) 1000001,msdown(2),"# ~d_L"
    Write(io_L,102) 2000001,msdown(1),"# ~d_R"
    id_sd(2) = 1000001
    id_sd(1) = 2000001
    c_sd(2) = "~d_L"
    c_sd(1) = "~d_R"
   End If
   If (Abs(rsup(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000002,msup(1),"# ~u_L"
    Write(io_L,102) 2000002,msup(2),"# ~u_R"
    id_su(1) = 1000002
    id_su(2) = 2000002
    c_su(1) = "~u_L"
    c_su(2) = "~u_R"
   Else
    Write(io_L,102) 1000002,msup(2),"# ~u_L"
    Write(io_L,102) 2000002,msup(1),"# ~u_R"
    id_su(2) = 1000002
    id_su(1) = 2000002
    c_su(2) = "~u_L"
    c_su(1) = "~u_R"
   End If
   If (Abs(rsdown(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000003,msdown(3),"# ~s_L"
    Write(io_L,102) 2000003,msdown(4),"# ~s_R"
    id_sd(3) = 1000003
    id_sd(4) = 2000003
    c_sd(3) = "~s_L"
    c_sd(4) = "~s_R"
   Else
    Write(io_L,102) 1000003,msdown(4),"# ~s_L"
    Write(io_L,102) 2000003,msdown(3),"# ~s_R"
    id_sd(4) = 1000003
    id_sd(3) = 2000003
    c_sd(4) = "~s_L"
    c_sd(3) = "~s_R"
   End If
   If (Abs(rsup(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000004,msup(3),"# ~c_L"
    Write(io_L,102) 2000004,msup(4),"# ~c_R"
    id_su(3) = 1000004
    id_su(4) = 2000004
    c_su(3) = "~c_L"
    c_su(4) = "~c_R"
   Else
    Write(io_L,102) 1000004,msup(4),"# ~c_L"
    Write(io_L,102) 2000004,msup(3),"# ~c_R"
    id_su(4) = 1000004
    id_su(3) = 2000004
    c_su(4) = "~c_L"
    c_su(3) = "~c_R"
   End If
   Write(io_L,102) 1000005,msdown(5),"# ~b_1"
   Write(io_L,102) 2000005,msdown(6),"# ~b_2"
   id_sd(5) = 1000005
   id_sd(6) = 2000005
   c_sd(5) = "~b_1"
   c_sd(6) = "~b_2"
   Write(io_L,102) 1000006,msup(5),"# ~t_1"
   Write(io_L,102) 2000006,msup(6),"# ~t_2"
   id_su(5) = 1000006
   id_su(6) = 2000006
   c_su(5) = "~t_1"
   c_su(6) = "~t_2"
  End If ! GenerationMixing
   
   If (HighScaleModel.Ne."RPexplicit") Then

! sleptons
    If (GenerationMixing) Then
     If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
 
      mat3R = Abs(RSn_pmns)
      Do i1=1,3
       MaxCont = Maxval(mat3R)
       Call FindPosition(3, mat3R, MaxCont, ii, jj)
       Select Case (jj)
       Case(1)
        id_snu(ii) = 1000012
        c_snu(ii) = "~nu_eL"
       Case(2)
        id_snu(ii) = 1000014
        c_snu(ii) = "~nu_muL"
       Case(3)
        id_snu(ii) = 1000016
        c_snu(ii) = "~nu_tauL"
       End Select
       mat3R(ii,:) = 0._dp
       mat3R(:,jj) = 0._dp
      End Do
      Do ii=1,3
       Write(io_L,102) id_snu(ii),msneut(ii),"# "//Trim(c_snu(ii))
      End Do

      i_zaehl = 1
      mat6R = Abs(RSl_pmns)
      Do i1=1,6
       MaxCont = Maxval(mat6R)
       Call FindPosition(6, mat6R, MaxCont, ii, jj)
       Select Case (jj)
       Case(1)
        id_sle(ii) = 1000011
        c_sle(ii) = "~e_L-"
        c_slep(ii) = "~e_L+"
       Case(2)
        id_sle(ii) = 1000013
        c_sle(ii) = "~mu_L-"
        c_slep(ii) = "~mu_L+"
       Case(4)
        id_sle(ii) = 2000011
        c_sle(ii) = "~e_R-"
        c_slep(ii) = "~e_R+"
       Case(5)
        id_sle(ii) = 2000013
        c_sle(ii) = "~mu_R-"
        c_slep(ii) = "~mu_R+"
       Case default
        id_sle(ii) = 1000000 * i_zaehl + 15
        c_sle(ii) = "~tau_"//bu(i_zaehl)//"-"
        c_slep(ii) = "~tau_"//bu(i_zaehl)//"+"
        i_zaehl = I_zaehl + 1
       End Select
       mat6R(ii,:) = 0._dp
       mat6R(:,jj) = 0._dp
      End Do
      Do ii=1,6 ! check ordering of staus
       If (id_sle(ii).Eq.1000015)  id_check(1) = ii
       If (id_sle(ii).Eq.2000015)  id_check(2) = ii
      End Do
      If (id_check(1).Gt.id_check(2)) Then ! switch ordering
       ii = id_check(2)  ! the lighter one
       id_sle(ii) = 1000015
       c_sle(ii) = "~tau_1-"
       c_slep(ii) = "~tau_1+"
       ii = id_check(1)  ! the heavier one
       id_sle(ii) = 2000015
       c_sle(ii) = "~tau_2-"
       c_slep(ii) = "~tau_2+"
      End If

      Do ii=1,6
       Write(io_L,102) id_sle(ii),mslepton(ii),"# "//Trim(c_sle(ii))
      End Do

     Else ! using mass ordering

      id_snu(1) = 1000012
      id_snu(2) = 1000014
      id_snu(3) = 1000016
      Do i1=1,3
       c_snu(i1) = "~nu_"//Bu(i1)
       Write(io_L,102) id_snu(i1),mSneut(i1),"# "//Trim(c_snu(i1))
      End Do

      id_sle(1) = 1000011
      id_sle(2) = 1000013
      id_sle(3) = 1000015
      id_sle(4) = 2000011
      id_sle(5) = 2000013
      id_sle(6) = 2000015
      Do i1=1,6
       c_sle(i1) = "~l_"//Bu(i1)
       c_slep(i1) = "~l+_"//Bu(i1)
       Write(io_L,102) id_sle(i1),mSlepton(i1),"# "//Trim(c_sle(i1))
      End Do
     End If

    Else ! .not.GenerationMixing

     id_snu = (/ 1000012, 1000014, 1000016 /)
     If (Abs(RSl_pmns(1,1)).Gt.0.5_dp) Then
      Write(io_L,102) 1000011,mslepton(1),"# ~e_L-"
      Write(io_L,102) 2000011,mslepton(2),"# ~e_R-"
      id_sle(1) = 1000011
      id_sle(2) = 2000011
      c_sle(1) = "~e_L-"
      c_sle(2) = "~e_R-"
      c_slep(1) = "~e_L+"
      c_slep(2) = "~e_R+"
     Else
      Write(io_L,102) 1000011,mslepton(2),"# ~e_L-"
      Write(io_L,102) 2000011,mslepton(1),"# ~e_R-"
      id_sle(2) = 1000011
      id_sle(1) = 2000011
      c_sle(2) = "~e_L-"
      c_sle(1) = "~e_R-"
      c_slep(2) = "~e_L+"
      c_slep(1) = "~e_R+"
     End If
     Write(io_L,102) 1000012,msneut(1),"# ~nu_eL"
     c_snu(1) = "~nu_eL"
     If (Abs(RSl_pmns(3,3)).Gt.0.5_dp) Then
      Write(io_L,102) 1000013,mslepton(3),"# ~mu_L-"
      Write(io_L,102) 2000013,mslepton(4),"# ~mu_R-"
      id_sle(3) = 1000013
      id_sle(4) = 2000013
      c_sle(3) = "~mu_L-"
      c_sle(4) = "~mu_R-"
      c_slep(3) = "~mu_L+"
      c_slep(4) = "~mu_R+"
     Else
      Write(io_L,102) 1000013,mslepton(4),"# ~mu_L-"
      Write(io_L,102) 2000013,mslepton(3),"# ~mu_R-"
      id_sle(4) = 1000013
      id_sle(3) = 2000013
      c_sle(4) = "~mu_L-"
      c_sle(3) = "~mu_R-"
      c_slep(4) = "~mu_L+"
      c_slep(3) = "~mu_R+"
     End If
     Write(io_L,102) 1000014,msneut(2),"# ~nu_muL"
     c_snu(2) = "~nu_muL"
     Write(io_L,102) 1000015,mslepton(5),"# ~tau_1-"
     Write(io_L,102) 2000015,mslepton(6),"# ~tau_2-"
     id_sle(5) = 1000015
     id_sle(6) = 2000015
     c_sle(5) = "~tau_1-"
     c_sle(6) = "~tau_2-"
     c_slep(5) = "~tau_1+"
     c_slep(6) = "~tau_2+"
     Write(io_L,102) 1000016,msneut(3),"# ~nu_tauL"
     c_snu(3) = "~nu_tauL"
    End If
   End If ! GenerationMixing

 ! gauginos/higgsinos
   Write(io_L,102) 1000021,mglu,"# ~g"
   ! checking for negative sign
   nr = ZeroC
   If (HighScaleModel.Eq."NMSSM") Then
    Do i1=1,5
     If (Sum(Abs(Real(N5(i1,:)))).Eq.0._dp) Then
      mNr(i1) = - mN5(i1)
      nr(i1,1:5) = (0._dp,-1._dp) * n5(i1,:)
     Else   
      mNr(i1) =  mN5(i1)
      nr(i1,1:5) = n5(i1,:)
     End If
    End Do

   Else If (HighScaleModel.Eq."RPexplicit") Then

    If (l_RP_Pythia) Then  ! Pythia only takes 4x4 matrix for neutralinos
                           ! and 2x2 for charginos
     mNr(1:7) = Abs(mN7(1:7))
     Do i1=1,4
      If (Sum(Abs(Real(N(i1,:)))).Lt.0.1_dp) Then
       nr(i1,1:4) = Aimag(n(i1,:))
      Else   
       nr(i1,1:4) = n(i1,:)
      End If
     End Do

    Else
     Do i1=1,7
      If (Sum(Abs(Real(N7(i1,:)))).Eq.0._dp) Then
       mNr(i1) = - mN7(i1)
       nr(i1,1:7) = Aimag(n7(i1,:))
      Else   
       mNr(i1) =  mN7(i1)
       nr(i1,1:7) = n7(i1,:)
      End If
     End Do
    End If

   Else
    Do i1=1,4
     If (Sum(Abs(Real(N(i1,:)))).Eq.0._dp) Then
      mNr(i1) = - mN(i1)
      nr(i1,1:4) = Aimag(n(i1,:))
     Else   
      mNr(i1) =  mN(i1)
      nr(i1,1:4) = n(i1,:)
     End If
    End Do
   End If

   If (HighScaleModel.Eq."RPexplicit") Then
    Do i1=1,7
     Write(io_L,102) id_n(i1),mnr(i1),"# "//Trim(c_c0(i1))
    End Do
    Do i1=1,5
     Write(io_L,102) id_c(i1),mc5(i1),"# "//Trim(c_cp(i1))
    End Do

   Else
    Write(io_L,102) 1000022,mnr(1),"# ~chi_10" 
    Write(io_L,102) 1000023,mnr(2),"# ~chi_20" 
    Write(io_L,102) 1000025,mnr(3),"# ~chi_30" 
    Write(io_L,102) 1000035,mnr(4),"# ~chi_40" 
    If (HighScaleModel.Eq."NMSSM") Write(io_L,102) 1000045,mnr(5),"# ~chi_50"
    Write(io_L,102) 1000024,mc(1),"# ~chi_1+" 
    Write(io_L,102) 1000037,mc(2),"# ~chi_2+"
   End If

   If (Maxval(MnuR).Gt.0._dp) Then
    Write(io_L,100) "# masses of right handed neutrinos"
    Write(io_L,100) "Block MnuR"
    Write(io_L,102) 1,mnur(1),"# m_nu_R_1"
    Write(io_L,102) 2,mnur(2),"# m_nu_R_2"
    Write(io_L,102) 3,mnur(3),"# m_nu_R_3"
   End If

! Mixing matrices 
  Write(io_L,100) "# Higgs mixing"
  If (HighScaleModel.Eq."NMSSM") Then
   Write(io_L,103) "Block HIGMIX Q=",Q, "# neutral scalar Higgs mixing"
   If (.Not.External_Higgs) Then
    RS03_save = RS03
    RG0 = 0._dp
    RG0(1,1) = - 1._dp / Sqrt(1._dp + tanb**2)
    RG0(2,2) = - RG0(1,1)
    RG0(3,3) = 1._dp
    RG0(1,2) = tanb * RG0(1,1)
    RG0(2,1) = RG0(1,2)
    RSpm = RG0(1:2,1:2)
!    RP03_save = Matmul(Transpose(RG0),RP03)
   End If

   Do i1=1,n_s0
    Do i2=1,n_s0
     Write(io_L,105) i1,i2,RS03_save(i1,i2),"# R_S0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
   Write(io_L,103) "Block AMIX Q=",Q, "# pseudoscalar Higgs mixing"
   Do i1=1,n_p0+1
    Do i2=1,n_p0+1
     Write(io_L,105) i1,i2,RP03_save(i1,i2) &
             & ,"# R_P0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
  Else If (HighScaleModel.Eq."RPexplicit") Then
   If (l_RP_Pythia) Then
    Write(io_L,100) "Block alpha # Effective Higgs mixing angle"
    Write(io_L,108) -Asin(RS0(1,1)),"# alpha"
    Write(io_L,103) "Block Hmix Q=",Q, "# Higgs mixing parameters"
    Write(io_L,104) 1,Real(mu,dp),"# mu"
    Write(io_L,104) 2,tanb_Q,"# tan[beta](Q)"
    Write(io_L,104) 3,vev_Q,"# v(Q)"
    Write(io_L,104) 4,mA2_Q,"# m^2_A(Q)"
    Write(io_L,100) "Block staumix  # stau mixing matrix"
    Write(io_L,105) 1,1,Real(RSl_pmns(5,5),dp),"# R_sta(1,1)"
    Write(io_L,105) 1,2,Real(RSl_pmns(5,6),dp),"# R_sta(1,2)"
    Write(io_L,105) 2,1,Real(RSl_pmns(6,5),dp),"# R_sta(2,1)"
    Write(io_L,105) 2,2,Real(RSl_pmns(6,6),dp),"# R_sta(2,2)"

   Else

   Write(io_L,103) "Block RVHMIX  Q=",Q, "# neutral scalar Higgs mixing"

   Do i1=1,n_s0
    Do i2=1,n_s0
     Write(io_L,105) i1,i2,RS05(i1,i2),"# R_S0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
   Write(io_L,103) "Block RVAMIX  Q=",Q, "# pseudoscalar Higgs mixing"
   Do i1=1,n_p0+1
    Do i2=1,n_p0+1
     Write(io_L,105) i1,i2,RP05(i1,i2) &
             & ,"# R_P0("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do 
   Write(io_L,103) "Block RVLMIX Q=",Q, "# charged Higgs mixing"
   Do i1=1,8
    Do i2=1,8
     Write(io_L,105) i1,i2,Real(Rspm8(i1,i2),dp) &
             & ,"# R_Spm("//Bu(i1)//","//Bu(i2)//")"
    End Do
   End Do
   End If

  Else
   Write(io_L,100) "Block alpha # Effective Higgs mixing angle"
   Write(io_L,108) -Asin(RS0(1,1)),"# alpha"
   Write(io_L,103) "Block Hmix Q=",Q, "# Higgs mixing parameters"
   Write(io_L,104) 1,Real(mu,dp),"# mu"
   Write(io_L,104) 2,tanb_Q,"# tan[beta](Q)"
   Write(io_L,104) 3,vev_Q,"# v(Q)"
   Write(io_L,104) 4,mA2_Q,"# m^2_A(Q)"
  End If

  If (generationmixing) Then
   Call WriteMatrixBlockC2(io_L, 6, RUsq_ckm, "USQmix" &
                         &, "u-sqark mixing matrix", "R_Su(")

   Call WriteMatrixBlockC2(io_L, 6, RDsq_ckm, "DSQmix" &
                         &, "d-sqark mixing matrix", "R_Sd(")

   If (HighScaleModel.Ne."RPexplicit") Then

    Call WriteMatrixBlockC2(io_L, 6, RSl_pmns, "SELmix" &
                         &, "slepton mixing matrix", "R_Sl(")

    Call WriteMatrixBlockC2(io_L, 3, RSn_pmns, "SNUmix" &
                         &, "sneutrino mixing matrix", "R_Sn(")

   End If

  Else ! .not.GenerationMixing

   Call WriteMatrixBlockC2(io_L, 2, RUsq_ckm(5:6,5:6), "stopmix" &
                         &, "stop mixing matrix", "R_st(")

   Call WriteMatrixBlockC2(io_L, 2, RDsq_ckm(5:6,5:6), "sbotmix" &
                         &, "sbottom mixing matrix", "R_sb(")

   If (HighScaleModel.Ne."RPexplicit") & 
    & Call WriteMatrixBlockC2(io_L, 2, RSl_pmns(5:6,5:6), "staumix" &
                            &, "stau mixing matrix", "R_sta(")
  End If

  If ((HighScaleModel.Eq."RPexplicit").And.(.Not.l_RP_Pythia)) Then
   Call WriteMatrixBlockC2(io_L, n_n+delta_n_rp, nr(1:n_n+delta_n_rp,1:n_n+delta_n_rp) &
                         &, "RVNmix", "/neutrino/neutralino mixing matrix", "N(")
   Call WriteMatrixBlockC2(io_L, 5, U5, "RVUmix" &
                         &, "lepton/chargino mixing matrix", "U(")
   Call WriteMatrixBlockC2(io_L, 5, V5, "RVVmix" &
                         &, "lepton/chargino mixing matrix", "V(")
  Else   
   Call WriteMatrixBlockC2(io_L, n_n, nr(1:n_n,1:n_n), "Nmix" &
                         &, "neutralino mixing matrix", "N(")
   Call WriteMatrixBlockC2(io_L, n_c, U, "Umix" &
                         &, "chargino mixing matrix", "U(")
   Call WriteMatrixBlockC2(io_L, n_c, V, "Vmix" &
                         &, "chargino mixing matrix", "V(")
  End If

  !---------------------------------------------------
  ! parameters + masses for SPheno.out
  !---------------------------------------------------
  If ( (HighScaleModel(2:5).Ne."MSSM").And.(HighScaleModel(1:4).Ne."MSSM")  &
     & .And. (HighScaleModel(1:2).Ne."RP") ) Then
   Write(io,*) "       g'             g             g_3"
   Write(io,5103) gauge_0
   Write(io,*) " "

   Write(io,*) "      Y_e            Y_mu          Y_tau"
   Call WriteComplexMatrix(io, Transpose(Y_l_0), .Not.GenerationMixing)
   Write(io,*) " "

   Write(io,*) "      Y_u            Y_c            Y_t"
   Call WriteComplexMatrix(io, Transpose(Y_u_0), .Not.GenerationMixing)
   Write(io,*) " "

   Write(io,*) "      Y_d            Y_s            Y_b"
   Call WriteComplexMatrix(io, Transpose(Y_d_0), .Not.GenerationMixing)
   Write(io,*) " "
  End If

  Write(io,*) " "
  Write(io,5101) "Parameters at the scale ", Q
  Write(io,*) " "

  If (HighScaleModel.Ne."RPexplicit") Then
   Call WriteMSSMParameters(io, CKM_Q, Yd, Yu, .True.)

   If (HighScaleModel.Eq."NMSSM") Then
    Call WriteMassesNMSSM(io, .True., mGlu, PhaseGlu, mC, U, V, mN5, N5  &
           &, mSneut, RSn_pmns, mSlepton, RSl_pmns, mSdown, RDsq_ckm, mSup &
           & , RUsq_ckm, mP03, RP03, mS03, RS03, mSpm, GenerationMixing)
   Else
    Call WriteMassesMSSM(io, .True., mGlu, PhaseGlu, mC, U, V, mN, N     &
           &, mSneut, RSn_pmns, mSlepton, RSl_pmns, mSdown, RDsq_ckm, mSup &
           & , RUsq_ckm, mP0, mS0, RS0, mSpm, GenerationMixing)
   End If
  End If
 !------------------
 ! branching ratios
 !------------------
 If (L_BR) Then
  BRmin100 = 100._dp * BRmin

  Write(io,*) " "
  Write(io,*) " Anti particles are marked with a * in case of"
  Write(io,*) " (s)neutrinos and (s)quarks in the decay section."
  Write(io,*) "                    Decay widths (GeV) and branching ratios"
  Write(io,*) " "

  !---------------------------------
  ! sleptons
  !---------------------------------
  Do i1=1,n_sl
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n + delta_n_rp
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    id_fp = (i1+1)/2
    Do i2=1,n_c + delta_c_rp
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_nu(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_nu(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_nu(id_fp)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" W-"
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_spm
     Do i2=1,n_sn
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_P0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,n_S0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Do i3 = 1,3
     Fnames(i_zaehl) = "gravitino "//Trim(c_lm(i3))
     Lnames(i_zaehl) = Trim(c_grav)//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_grav
     id_d_2(i_zaehl,2) = id_l(i3)
     i_zaehl = i_zaehl+1
    End Do

    Call WriteDecays2(io, " slepton_"//Bu(i1) , Fnames &
                     &, gP_Sl(i1,:), 100*BR_Sl(i1,:), gT_Sl(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n + delta_n_rp
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_lm(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_l(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = (i1+1)/2
    Do i2=1,n_c + delta_c_rp
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_nu(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_nu(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_nu(id_fp)
     i_zaehl = i_zaehl+1
    End Do
    i2=id_fp
    If (i2.Eq.1) Then
     c_sfermion = "e-sneutrino"
    Else If (i2.Eq.2) Then
     c_sfermion = "mu-sneutrino"
    Else If (i2.Eq.3) Then
     c_sfermion = "tau-sneutrino"
    End If
    Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W-"
    Lnames(i_zaehl) =  Trim(c_snu(i2))//" W-"
    id_d_2(i_zaehl,1) = id_snu(i2)
    id_d_2(i_zaehl,2) = -id_W
    i_zaehl = i_zaehl+1
    Do i3=1,n_spm
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sm(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_sm(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     i_zaehl = i_zaehl+1
    End Do
    If (i1.Le.2) Then
     c_sfermion = "selectron"
    Else If (i1.Le.4) Then
     c_sfermion = "smuon"
    Else 
     c_sfermion = "stau"
    End If

    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    i3 = (i1+1)/2
    Fnames(i_zaehl) = "gravitino "//Trim(c_lm(i3))
    Lnames(i_zaehl) = Trim(c_grav)//" "//Trim(c_lm(i3))
    id_d_2(i_zaehl,1) = id_grav
    id_d_2(i_zaehl,2) = id_l(i3)
    i_zaehl = i_zaehl+1

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2
    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Sl(i1,:), 100*BR_Sl(i1,:), gT_Sl(i1), BrMin100)
   End If ! GenerationMixing

   c_m = c_sle(i1)
   Write(io_L,200) id_sle(i1),gT_sl(i1),Trim(c_m)
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sl(i1,i2).Gt.BrMin) Write(io_L,201) BR_Sl(i1,i2),2,id_d_2(i2,:), &
            &                  Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! sneutrinos
  !---------------------------------
  Do i1=1,n_sn
   If (GenerationMixing) Then
    i_zaehl = 1
    i3 = i1
    Do i2=1,n_n + delta_n_rp
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_nu(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_nu(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_nu(i3)
      i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_c + delta_c_rp
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_c(i2)
      id_d_2(i_zaehl,2) = id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_sl
     Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" W+"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_spm
     Do i2=1,n_sl
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" Z"
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_P0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,n_S0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " Sneutrino_"//Bu(i1) , Fnames &
                     &, gP_Sn(i1,:), 100*BR_Sn(i1,:), gT_Sn(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = i1
    Do i2=1,n_n + delta_n_rp
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i1
    Do i2=1,n_c + delta_c_rp
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_lm(id_fp))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_lm(id_fp))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_l(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    If (i1.Eq.1) Then
     c_sfermion = "selectron_"
    Else If (i1.Eq.2) Then
     c_sfermion = "smuon_"
    Else If (i1.Eq.3) Then
     c_sfermion = "stau_"
    End If
    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sle(i2+(i1-1)*2))//" W+"
     id_d_2(i_zaehl,1) = id_sle(i2+(i1-1)*2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2+(i1-1)*2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sle(i2+(i1-1)*2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Eq.1) c_sfermion="e-sneutrino"
    If (i1.Eq.2) c_sfermion="mu-sneutrino"
    If (i1.Eq.3) c_sfermion="tau-sneutrino"

    Call WriteDecays2(io, Trim(c_sfermion) , Fnames &
                     &, gP_Sn(i1,:), 100*BR_Sn(i1,:), gT_Sn(i1), BrMin100)
   End If ! GenerationMixing

   c_m = c_snu(i1)
   Write(io_L,200) id_snu(i1),gT_sn(i1),Trim(c_m)
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sn(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Sn(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! d-squarks
  !---------------------------------
  Do i1=1,n_sd
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n + delta_n_rp
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_c + delta_c_rp
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = -id_c(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,3
     Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_glu
     id_d_2(i_zaehl,2) = id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" W-"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_spm
     Do i2=1,n_su
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_P0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,n_S0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " s-down_"//Bu(i1) , Fnames &
                     &, gP_Sd(i1,:), 100*BR_Sd(i1,:), gT_Sd(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n + delta_n_rp
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i3
    Do i2=1,n_c + delta_c_rp
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_u(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_u(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_u(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))
    Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_d(i3))
    id_d_2(i_zaehl,1) = id_glu
    id_d_2(i_zaehl,2) = id_d(i3)
    i_zaehl = i_zaehl+1

    If (i1.Le.2) Then
     c_sfermion = "s-up_"
     id_fp = 0
    Else If (i1.Le.4) Then
     c_sfermion = "s-charm_"
     id_fp = 2
    Else 
     c_sfermion = "stop_"
     id_fp = 4
    End If

    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_su(i2 + id_fp))//" W-"
     id_d_2(i_zaehl,1) = id_su(i2 + id_fp)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2+id_fp))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_su(i2+id_fp)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Le.2) Then
     c_sfermion = "s-down_"
    Else If (i1.Le.4) Then
     c_sfermion = "s-strange_"
    Else 
     c_sfermion = "sbottom_"
    End If


    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2
    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Sd(i1,:), 100*BR_Sd(i1,:), gT_Sd(i1), BrMin100)
   End If ! GenerationMixing

   c_m = c_sd(i1)
   Write(io_L,200) id_sd(i1),gT_sd(i1),Trim(c_m)
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sd(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Sd(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! u-squarks
  !---------------------------------
  Do i1=1,n_su
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n + delta_n_rp
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_c + delta_c_rp
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_c(i2)
      id_d_2(i_zaehl,2) = id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,3
     Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_glu
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" W+"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_spm
     Do i2=1,n_sd
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" Z"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i3=1,n_P0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,n_S0
     Do i2=1,i1-1
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " s-up_"//Bu(i1) , Fnames &
                     &, gP_Su(i1,:), 100*BR_Su(i1,:), gT_Su(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n + delta_n_rp
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i3
    Do i2=1,n_c + delta_c_rp
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_d(id_fp))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_d(id_fp))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_d(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    Fnames(i_zaehl) =  "gluino "//Trim(c_u(id_fp))
    Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_u(id_fp))
    id_d_2(i_zaehl,1) = id_glu
    id_d_2(i_zaehl,2) = id_u(id_fp)
    i_zaehl = i_zaehl+1

    If (i1.Le.2) Then
     c_sfermion = "s-down_"
     id_fp = 0
    Else If (i1.Le.4) Then
     c_sfermion = "s-strange_"
     id_fp = 2
    Else 
     c_sfermion = "sbottom_"
     id_fp = 4
    End If

    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sd(i2 + id_fp))//" W+"
     id_d_2(i_zaehl,1) = id_sd(i2 + id_fp)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2+id_fp))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_su(i2+id_fp)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Le.2) Then
     c_sfermion = "s-up_"
    Else If (i1.Le.4) Then
     c_sfermion = "s-charm_"
    Else 
     c_sfermion = "stop_"
    End If


    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" Z"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2

    If (i1.Eq.5) Then ! 3-body decays of lighter stop
     id_d_3 = 0
     Fnames(55) = "neutralino_1 c"
     Lnames(55) = Trim(c_c0(1))//" c"
     id_d_2(55,1) = id_n(1)
     id_d_2(55,2) = id_u(2)
     Fnames(56) = "neutralino_2 c"
     Lnames(56) = Trim(c_c0(2))//" c"
     id_d_2(56,1) = id_n(2)
     id_d_2(56,2) = id_u(2)
     i_zaehl = 57
     Fnames(57) = "neutralino_1 W b"
     Lnames(57) = Trim(c_c0(1))//" W+ b"
     id_d_3(1,1) = id_n(1)
     id_d_3(1,2) = id_W
     id_d_3(1,3) = id_d(3)
     Fnames(58) = "e-sneutrino e+ b"
     Lnames(58) = Trim(c_snu(1))//" e+ b"
     id_d_3(2,1) = id_snu(1)
     id_d_3(2,2) = -id_l(1)
     id_d_3(2,3) = id_d(3)
     Fnames(59) = "mu-sneutrino mu+ b"
     Lnames(59) = Trim(c_snu(2))//" mu+ b"
     id_d_3(3,1) = id_snu(2)
     id_d_3(3,2) = -id_l(2)
     id_d_3(3,3) = id_d(3)
     Fnames(60) = "tau-sneutrino tau+ b"
     Lnames(60) = Trim(c_snu(3))//" tau+ b"
     id_d_3(4,1) = id_snu(3)
     id_d_3(4,2) = -id_l(3)
     id_d_3(4,3) = id_d(3)
     Fnames(61) = "selectron_1 nu_e b"
     Lnames(61) = Trim(c_slep(1))//" nu_e b"
     id_d_3(5,1) = -id_sle(1)
     id_d_3(5,2) = id_nu(1)
     id_d_3(5,3) = id_d(3)
     Fnames(62) = "selectron_2 nu_e b"
     Lnames(62) = Trim(c_slep(2))//" nu_e b"
     id_d_3(6,1) = -id_sle(2)
     id_d_3(6,2) = id_nu(1)
     id_d_3(6,3) = id_d(3)
     Fnames(63) = "smuon_1 nu_mu b"
     Lnames(63) = Trim(c_slep(3))//" nu_mu b"
     id_d_3(7,1) = -id_sle(3)
     id_d_3(7,2) = id_nu(2)
     id_d_3(7,3) = id_d(3)
     Fnames(64) = "smuon_2 nu_mu b"
     Lnames(64) = Trim(c_slep(4))//" nu_mu b"
     id_d_3(8,1) = -id_sle(4)
     id_d_3(8,2) = id_nu(2)
     id_d_3(8,3) = id_d(3)
     Fnames(65) = "stau_1 nu_tau b"
     Lnames(65) = Trim(c_slep(5))//" nu_tau b"
     id_d_3(9,1) = -id_sle(5)
     id_d_3(9,2) = id_nu(3)
     id_d_3(9,3) = id_d(3)
     Fnames(66) = "stau_2 nu_tau b"
     Lnames(66) = Trim(c_slep(6))//" nu_tau b"
     id_d_3(10,1) = -id_sle(6)
     id_d_3(10,2) = id_nu(3)
     id_d_3(10,3) = id_d(3)
    End If

    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Su(i1,:), 100*BR_Su(i1,:), gT_Su(i1), BrMin100)
   End If ! GenerationMixing

   c_m = c_su(i1)
   Write(io_L,200) id_su(i1),gT_su(i1),Trim(c_m)
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Su(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Su(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do

   If ((i1.Eq.5).And.(Maxval(BR_Su(i1,57:66)).Gt.BRmin)) Then
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,10
      If  (BR_Su(i1,56+i2).Gt.BrMin) Write(io_L,202) BR_Su(i1,56+i2),3 &
             & ,id_d_3(i2,:), Trim(c_m)//" -> "//Trim(Lnames(56+i2))//")"
    End Do
   End If

  End Do   

  !--------------
  ! charginos
  !--------------

  Do i1=c_min,n_c+delta_c_rp
   i_zaehl = 1
   If (GenerationMixing) Then
    Do i2=1,n_sl
     i3 = (i2+1)/2
     Fnames(i_zaehl) = "slepton_"//Bu(i2)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_slep(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = -id_sle(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Do i3 = 1,3
      Fnames(i_zaehl) = "sneutrino_"//Bu(i2)//" "//Trim(c_lp(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_lp(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = -id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_su
     Do i3 = 1,3
      Fnames(i_zaehl) = "s-up_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = -id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_sd
     Do i3 = 1,3
      Fnames(i_zaehl) = "s-down_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = -id_sd(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

   Else ! GenerationMixing

    Do i2=1,n_sl
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="selectron_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="smuon_"
      id_f = i2 - 2
     Else 
      c_sfermion="stau_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_slep(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = -id_sle(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     i3 = i2
     If (i2.Eq.1) c_sfermion = "e-sneutrino"
     If (i2.Eq.2) c_sfermion = "mu-sneutrino"
     If (i2.Eq.3) c_sfermion = "tau-sneutrino"
     Fnames(i_zaehl) = Trim(c_sfermion)//" "//Trim(c_lp(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_lp(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_l(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-up_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-charm_"
      id_f = i2 - 2
     Else 
      c_sfermion="stop_"
      id_f = i2 - 2
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sd
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-down_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-strange_"
      id_f = i2 - 2
     Else 
      c_sfermion="sbottom_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = -id_sd(i2)
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do

   End If ! GenerationMixing

   Do i2=1,n_n+delta_n_rp
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" W+"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" W+"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_W
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_spm
    Do i2=1,n_n+delta_n_rp
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_sp(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i2=1,i1-1
    Fnames(i_zaehl) = "chargino_"//Bu(i2)//" Z"
    Lnames(i_zaehl) =  Trim(c_cp(i2))//" Z"
    id_d_2(i_zaehl,1) = id_c(i2)
    id_d_2(i_zaehl,2) = id_Z
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_p0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_p0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i3=1,n_s0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   i_c2 = i_zaehl-1
   gP = 0._dp
   BR = 0._dp

   i_eff = i1 - delta_c_rp

   gP(1:i_c2) = gP_C2(i_eff,1:i_c2)
   BR(1:i_c2) = BR_C2(i_eff,1:i_c2)

   c_m = c_cp(i1)
   Write(io_L,200) id_c(i1),gT_C(i_eff),Trim(c_cp(i1))
   If (Sum(gp_C2(i_eff,:)).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,i_zaehl - 1
     If (BR_C2(i_eff,i2).Gt.BrMin) Write(io_L,201) BR_C2(i_eff,i2) &
            & ,2,id_d_2(i2,:), Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
    End Do
   End If

   !-------------------------------------------------
   ! 3-body decays
   !-------------------------------------------------
   If (Maxval(BR_C3(i_eff,:)).Gt.BRmin) Then
    If (GenerationMixing) Then
     Do i2=1,n_n+delta_n_rp
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl-i_c2,1) = id_n(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
        id_d_3(i_zaehl-i_c2,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3-delta_c_rp
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_nu(i3)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_glu
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
        id_d_3(i_zaehl-i_c2,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
        id_d_3(i_zaehl-i_c2,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3-delta_c_rp
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
        id_d_3(i_zaehl-i_c2,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      If (delta_c_rp.Lt.3) Then
       Fnames(i_zaehl) = "chargino_"//Bu(i2)//"nu_i nu_j"
       Lnames(i_zaehl) = Trim(c_cp(i2))//"nu_i nu_j"
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_nu(1)
       id_d_3(i_zaehl-i_c2,3) = id_nu(3)
       i_zaehl = i_zaehl+1    
      End If
     End Do
  
    Else ! GenerationMixing

     Do i2=1,n_n+delta_n_rp
      Do i3=1,3
       Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,3-delta_c_rp
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_nu(i3)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_u(i3))
      Lnames(i_zaehl) = Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
      id_d_3(i_zaehl-i_c2,1) = id_glu
      id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = id_u(i3)
      i_zaehl = i_zaehl+1    
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = Trim(c_cp(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i3))
       Lnames(i_zaehl) = Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,3-delta_c_rp
       Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i3)
       i_zaehl = i_zaehl+1    
      End Do
      If (delta_c_rp.Lt.3) Then
       Fnames(i_zaehl) = "chargino_"//Bu(i2)//"nu_i nu_i"
       Lnames(i_zaehl) = Trim(c_cp(i2))//"nu_i nu_i"
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_nu(1)
       id_d_3(i_zaehl-i_c2,3) = id_nu(1)
       i_zaehl = i_zaehl+1    
      End If
     End Do

    End If ! GenerationMixing

    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,i_zaehl - i_c2 - 1
     If (BR_C3(i_eff,i2).Gt.BrMin) Write(io_L,202) BR_C3(i_eff,i2),3 &
             & ,id_d_3(i2,:) , Trim(c_m)//" -> "//Trim(Lnames(i2+i_c2))//")"
    End Do

    BR(i_c2+1:i_zaehl-1) = BR_C3(i_eff,1:i_zaehl-i_c2-1)
    gP(i_c2+1:i_zaehl-1) = gP_C3(i_eff,1:i_zaehl-i_c2-1)
   End If ! check for maxval(BR)

   Call WriteDecays2(io, " chargino_"//Bu(i1) , Fnames &
                       &, gP, 100*BR, gT_C(i_eff), BrMin100)

  End Do ! i1


  !--------------
  ! neutralinos
  !--------------
  Do i1=n_min,n_n+delta_n_rp
   i_zaehl = 1
   If (GenerationMixing) Then
    Do i2=1,n_sl
     Do i3 = 1,3
      Fnames(i_zaehl) = "slepton_"//Bu(i2)//" "//Trim(c_lp(i3))
      Fnames(i_zaehl+1) = "slepton_"//Bu(i2)//"^+ "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_lp(i3))
      Lnames(i_zaehl+1) =  Trim(c_slep(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = -id_l(i3)
      id_d_2(i_zaehl+1,1) = -id_sle(i2)
      id_d_2(i_zaehl+1,2) = id_l(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
    Do i2=1,n_sn
     i3 = i2
     Fnames(i_zaehl) = "sneutrino_"//Bu(i2)//" "//Trim(c_snu(i3))//"^*"
     Fnames(i_zaehl+1) = "sneutrino_"//Bu(i2)//"^* "//Trim(c_snu(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_nu(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_snu(i2))//"^* "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_nu(i3)
     id_d_2(i_zaehl+1,1) = -id_snu(i2)
     id_d_2(i_zaehl+1,2) = id_nu(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_su
     Do i3 = 1,3
      Fnames(i_zaehl) = "sup_"//Bu(i2)//" "//Trim(c_u(i3))//"^*"
      Fnames(i_zaehl+1) = "sup_"//Bu(i2)//"^* "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
      Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = -id_u(i3)
      id_d_2(i_zaehl+1,1) = -id_su(i2)
      id_d_2(i_zaehl+1,2) = id_u(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
    Do i2=1,n_sd
     Do i3 = 1,3
      Fnames(i_zaehl) = "sdown_"//Bu(i2)//" "//Trim(c_d(i3))//"^*"
      Fnames(i_zaehl+1) = "sdown_"//Bu(i2)//"^* "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
      Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_d(i3)
      id_d_2(i_zaehl+1,1) = -id_sd(i2)
      id_d_2(i_zaehl+1,2) = id_d(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
   Else ! GenerationMixing

    Do i2=1,n_sl
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="selectron_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="smuon_"
      id_f = i2 - 2
     Else 
      c_sfermion="stau_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_lp(i3))
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^+ "//Trim(c_lm(i3))
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_lp(i3))
     Lnames(i_zaehl+1) =  Trim(c_slep(i2))//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = -id_l(i3)
     id_d_2(i_zaehl+1,1) = -id_sle(i2)
     id_d_2(i_zaehl+1,2) = id_l(i3)
     i_zaehl = i_zaehl+2
    End Do

    Do i2=1,n_sn
     i3 = i2
     If (i2.Eq.1) c_sfermion = "e-sneutrino"
     If (i2.Eq.2) c_sfermion = "mu-sneutrino"
     If (i2.Eq.3) c_sfermion = "tau-sneutrino"
     Fnames(i_zaehl) = Trim(c_sfermion)//" "//Trim(c_snu(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//"^* "//Trim(c_snu(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_nu(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_snu(i2))//"^* "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_nu(i3)
     id_d_2(i_zaehl+1,1) = -id_snu(i2)
     id_d_2(i_zaehl+1,2) = id_nu(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_su
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-up_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-charm_"
      id_f = i2 - 2
     Else 
      c_sfermion="stop_"
      id_f = i2 - 2
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_u(i3)
     id_d_2(i_zaehl+1,1) = -id_su(i2)
     id_d_2(i_zaehl+1,2) = id_u(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_sd
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-down_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-strange_"
      id_f = i2 - 2
     Else 
      c_sfermion="sbottom_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     id_d_2(i_zaehl+1,1) = -id_sd(i2)
     id_d_2(i_zaehl+1,2) = id_d(i3)
     i_zaehl = i_zaehl+2
    End Do
    
   End If ! GenerationMixing

   Do i2=1,n_c+delta_c_rp
    Fnames(i_zaehl) = "chargino_"//Bu(i2)//"^+ W-"
    Fnames(i_zaehl+1) = "chargino_"//Bu(i2)//"^- W+"
    Lnames(i_zaehl) =  Trim(c_cp(i2))//" W-"
    Lnames(i_zaehl+1) =  Trim(c_cm(i2))//" W+"
    id_d_2(i_zaehl,1) = id_c(i2)
    id_d_2(i_zaehl,2) = -id_W
    id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl+2 
   End Do

   Do i3=1,n_spm
    Do i2=1,n_c+delta_c_rp
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//"^+ "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = "chargino_"//Bu(i2)//"^- "//Trim(c_sm(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) =  Trim(c_cm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = -id_sp(i3)
     id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl+2   
    End Do
   End Do


   Do i2=1,i1-1
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" Z"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" Z"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_Z
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_p0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_p0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i3=1,n_s0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Fnames(i_zaehl) = "Gravitino photon"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_ph
   i_zaehl = i_zaehl + 1
   Fnames(i_zaehl) = "Gravitino Z"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_Z
   i_zaehl = i_zaehl + 1
   Fnames(i_zaehl) = "Gravitino h0"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_s0(1)
   i_zaehl = i_zaehl + 1

   If (HighScaleModel.Eq."NMSSM") Then
    i_zaehl = 150 
   Else If (HighScaleModel.Eq."RPexplicit") Then
    i_zaehl = 186
   Else
    i_zaehl = 150 
   End If

   Do i2=1,i1-1
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" photon"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" photon"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_ph
    i_zaehl = i_zaehl+1    
   End Do

   i_c2 = i_zaehl-1
   gP = 0._dp
   BR = 0._dp
   If (HighScaleModel.Eq."NMSSM") Then
    gT = gT_N5(i1)
    gP(1:i_c2) = gP_N2(i1,1:i_c2)
    BR(1:i_c2) = BR_N2(i1,1:i_c2)

   Else If (HighScaleModel.Eq."RPexplicit") Then
    i_eff = i1 - 3
    gT = gT_N(i_eff)
    gP(1:i_c2) = gP_N4_2(i_eff,1:i_c2)
    BR(1:i_c2) = BR_N4_2(i_eff,1:i_c2)

   Else
    gT = gT_N(i1)
    gP(1:i_c2) = gP_N4_2(i1,1:i_c2)
    BR(1:i_c2) = BR_N4_2(i1,1:i_c2)
   End If
   c_m = c_c0(i1)

    Write(io_L,200) id_n(i1),gT,Trim(c_c0(i1))
   If (Sum(gP).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,i_zaehl - 1
     If (BR(i2).Gt.BrMin) Write(io_L,201) BR(i2),2,id_d_2(i2,:), &
            &                   Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
    End Do
   End If

   !-------------------------------------------------
   ! 3-body decays
   !-------------------------------------------------
   l_3body = .False.
   If (HighScaleModel.Eq."NMSSM") Then
    If (Maxval(BR_N3(i1,:)).Gt.BRmin) l_3body = .True.
   Else If (HighScaleModel.Eq."RPexplicit") Then
    If (Maxval(BR_N4_3(i1-3,:)).Gt.BRmin) l_3body = .True.
   Else
    If (Maxval(BR_N4_3(i1,:)).Gt.BRmin) l_3body = .True.
   End If

   If (l_3body) Then
    If (GenerationMixing) Then
     Do i2=1,n_c+delta_c_rp
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_c(i2)
        id_d_3(i_zaehl,2) = id_d(i3)
        id_d_3(i_zaehl,3) = -id_u(i4)
        id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
        i_zaehl = i_zaehl+2
       End Do
      End Do
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_l(i3)
       id_d_3(i_zaehl,3) = -id_nu(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do   

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_u(i3)
       id_d_3(i_zaehl,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_d(i3)
       id_d_3(i_zaehl,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i2=1,i1-1+delta_n_rp
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_u(i3)
        id_d_3(i_zaehl,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_d(i3)
        id_d_3(i_zaehl,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" nu_i nu_j"
      Lnames(i_zaehl) = Trim(c_c0(i2))//" nu_i nu_j"
      id_d_3(i_zaehl,1) = id_n(i2)
      id_d_3(i_zaehl,2) = - id_nu(1)
      id_d_3(i_zaehl,3) = id_nu(3)
      i_zaehl = i_zaehl+1    
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_l(i3)
        id_d_3(i_zaehl,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
     
     End Do
  
    Else ! GenerationMixing

     Do i2=1,n_c+delta_c_rp
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_d(i3)
       id_d_3(i_zaehl,3) = -id_u(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
      Do i3=1,3-delta_c_rp
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_l(i3)
       id_d_3(i_zaehl,3) = -id_nu(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do   

      Do i3=1,3
       i4 = i3
       Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_u(i3)
       id_d_3(i_zaehl,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do

      Do i3=1,3
       i4=i3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_d(i3)
       id_d_3(i_zaehl,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do

      Do i2=1,i1-1
       Do i3=1,3
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_u(i3)
        id_d_3(i_zaehl,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
       Do i3=1,3
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_d(i3)
        id_d_3(i_zaehl,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
       If (delta_n_rp.Lt.3) Then
        Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" nu_i nu_i"
        Lnames(i_zaehl) = Trim(c_c0(i2))//" nu_i nu_i"
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_nu(1)
        id_d_3(i_zaehl,3) = id_nu(3)
        i_zaehl = i_zaehl+1    
       End If
       Do i3=1,3 - delta_n_rp 
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_l(i3)
        id_d_3(i_zaehl,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do    
      End Do
 
     If (HighScaleModel.Eq."RPexplicit") Then
      Do i2=1,n_c+delta_c_rp
       Do i3=1,Min(i2,3)
        If (i2.Eq.i3) Then
         Fnames(i_zaehl) = "neutrino_e"//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Lnames(i_zaehl) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         id_d_3(i_zaehl,1) = id_n(1)
         id_d_3(i_zaehl,2) = - id_l(i2)
         id_d_3(i_zaehl,3) = id_l(i3)
         i_zaehl = i_zaehl+1
        Else
         Fnames(i_zaehl) = "neutrino_e"//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Fnames(i_zaehl+1) = &
            & "neutrino_e"//" "//Trim(c_cp(i3))//" "//Trim(c_cm(i2))
         Lnames(i_zaehl) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Lnames(i_zaehl+1) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i3))//" "//Trim(c_cm(i2))
         id_d_3(i_zaehl,1) = id_n(1)
         id_d_3(i_zaehl,2) = - id_c(i2)
         id_d_3(i_zaehl,3) = id_c(i3)
         id_d_3(i_zaehl+1,1) = id_n(1)
         id_d_3(i_zaehl+1,2:3) = -id_d_3(i_zaehl,2:3) 
         i_zaehl = i_zaehl+2
        End If
       End Do
      End Do
      Do i2=4,i1-1
       Do i3=1,3
        Do i4=1,i3
         If (i3.Eq.i4) Then
          Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          id_d_3(i_zaehl,1) = id_n(i2)
          id_d_3(i_zaehl,2) = - id_l(i3)
          id_d_3(i_zaehl,3) = id_l(i4)
          i_zaehl = i_zaehl+1
         Else
          Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Fnames(i_zaehl+1) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lm(i3))//" "//Trim(c_lp(i4))
          Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Lnames(i_zaehl+1) = &
            & Trim(c_c0(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_lp(i4))
          id_d_3(i_zaehl,1) = id_n(i2)
          id_d_3(i_zaehl,2) = - id_l(i3)
          id_d_3(i_zaehl,3) = id_l(i4)
          id_d_3(i_zaehl+1,1) = id_n(i2)
          id_d_3(i_zaehl+1,2:3) = -id_d_3(i_zaehl,2:3) 
          i_zaehl = i_zaehl+2
         End If
        End Do
       End Do
      End Do
      Do i2=4,i1-1
       Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_nu(1))//" "//Trim(c_nu(1))
       Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_c0(1))//" "//Trim(c_c0(1))
       id_d_3(i_zaehl,1) = id_n(i2)
       id_d_3(i_zaehl,2) = id_n(1)
       id_d_3(i_zaehl,3) = id_n(1)
       i_zaehl = i_zaehl+1
      End Do
      Fnames(i_zaehl) = Trim(c_nu(1))//" "//Trim(c_nu(1))//" "//Trim(c_nu(1))
      Lnames(i_zaehl) =  Trim(c_c0(1))//" "//Trim(c_c0(1))//" "//Trim(c_c0(1))
      id_d_3(i_zaehl,1) = id_n(1)
      id_d_3(i_zaehl,2) = id_n(1)
      id_d_3(i_zaehl,3) = id_n(1)
      i_zaehl = i_zaehl+1
     End If
    End If ! GenerationMixing
   
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"

    If (HighScaleModel.Eq."NMSSM") Then
     BR(i_c2+1:i_zaehl-1) = BR_N3(i1,1:i_zaehl-i_c2-1)
     gP(i_c2+1:i_zaehl-1) = gP_N3(i1,1:i_zaehl-i_c2-1)
    Else If (HighScaleModel.Eq."RPexplicit") Then
     BR(i_c2+1:i_zaehl-1) = BR_N4_3(i1-3,1:i_zaehl-i_c2-1)
     gP(i_c2+1:i_zaehl-1) = gP_N4_3(i1-3,1:i_zaehl-i_c2-1)

    Else
     BR(i_c2+1:i_zaehl-1) = BR_N4_3(i1,1:i_zaehl-i_c2-1)
     gP(i_c2+1:i_zaehl-1) = gP_N4_3(i1,1:i_zaehl-i_c2-1)
    End If

    Do i2=i_c2 + 1,i_zaehl - 1
     If (BR(i2).Gt.BrMin) Write(io_L,202) BR(i2),3,id_d_3(i2,:) &
             &                , Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
    End Do

   End If

  Call WriteDecays2(io, " neutralino_"//Bu(i1), Fnames, gP, 100*BR, gT, BrMin100)
  End Do

  !--------------
  ! gluino
  !--------------
  i_zaehl = 1
  If (GenerationMixing) Then
   Do i2=1,n_su
    Do i3 = 1,3
     Fnames(i_zaehl) = "sup_"//Bu(i2)//" "//Trim(c_u(i3))//"^*"
     Fnames(i_zaehl+1) = "sup_"//Bu(i2)//"^* "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_u(i3)
     id_d_2(i_zaehl+1,1) = -id_su(i2)
     id_d_2(i_zaehl+1,2) = id_u(i3)
     i_zaehl = i_zaehl+2
    End Do
   End Do
    
   Do i2=1,n_sd
    Do i3 = 1,3
     Fnames(i_zaehl) = "sdown_"//Bu(i2)//" "//Trim(c_d(i3))//"^*"
     Fnames(i_zaehl+1) = "sdown_"//Bu(i2)//"^* "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     id_d_2(i_zaehl+1,1) = -id_sd(i2)
     id_d_2(i_zaehl+1,2) = id_d(i3)
     i_zaehl = i_zaehl+2
    End Do
   End Do
 
    
  Else ! GenerationMixing

   Do i2=1,n_su
    i3 = (i2+1)/2
    If (i2.Le.2) Then
     c_sfermion="s-up_"
     id_f = i2
    Else If (i2.Le.4) Then
     c_sfermion="s-charm_"
     id_f = i2 - 2
    Else 
     c_sfermion="stop_"
     id_f = i2 - 2
    End If
    Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))//"^*"
    Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(i3))
    Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
    Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
    id_d_2(i_zaehl,1) = id_su(i2)
    id_d_2(i_zaehl,2) = -id_u(i3)
    id_d_2(i_zaehl+1,1) = -id_su(i2)
    id_d_2(i_zaehl+1,2) = id_u(i3)
    i_zaehl = i_zaehl+2
   End Do
    
   Do i2=1,n_sd
    i3 = (i2+1)/2
    If (i2.Le.2) Then
     c_sfermion="s-down_"
     id_f = i2
    Else If (i2.Le.4) Then
     c_sfermion="s-strange_"
     id_f = i2 - 2
    Else 
     c_sfermion="sbottom_"
     id_f = i2 - 4
    End If
    Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))//"^*"
    Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_d(i3))
    Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
    Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
    id_d_2(i_zaehl,1) = id_sd(i2)
    id_d_2(i_zaehl,2) = -id_d(i3)
    id_d_2(i_zaehl+1,1) = -id_sd(i2)
    id_d_2(i_zaehl+1,2) = id_d(i3)
    i_zaehl = i_zaehl+2
   End Do
    
   id_f = 1
   c_sfermion="stop_"
   Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(2))//"^*"
   Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(2))
   Lnames(i_zaehl) =  Trim(c_su(5))//" "//Trim(c_u(2))//"^*"
   Lnames(i_zaehl+1) =  Trim(c_su(5))//" "//Trim(c_u(2))
   id_d_2(i_zaehl,1) = id_su(5)
   id_d_2(i_zaehl,2) = -id_u(2)
   id_d_2(i_zaehl+1,1) = -id_su(5)
   id_d_2(i_zaehl+1,2) = id_u(2)
   i_zaehl = i_zaehl+2

  End If ! GenerationMixing

  Do i2=1,n_n+delta_n_rp
   Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" gluon"
   Lnames(i_zaehl) =  Trim(c_c0(i2))//" gluon"
   id_d_2(i_zaehl,1) = id_n(i2)
   id_d_2(i_zaehl,2) = id_gl
   i_zaehl = i_zaehl+1    
  End Do

  i_c2 = i_zaehl-1
  gP = 0._dp
  BR = 0._dp
  gP(1:i_c2) = gP_Glu2(1:i_c2)
  BR(1:i_c2) = BR_Glu2(1:i_c2)

  c_m = c_gl

  Write(io_L,200) id_glu,gT_Glu,Trim(c_gl)
  If (Sum(gp_Glu2).Gt.0._dp) Then ! 2-body decays
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1 
    If (BR_Glu2(i2).Gt.BrMin) Write(io_L,201) BR_Glu2(i2),2,id_d_2(i2,:), &
            &                   Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do
  End If

  !-------------------------------------------------
  ! 3-body decays
  !-------------------------------------------------
  If (Maxval(BR_Glu3).Gt.BRmin) Then
   If (GenerationMixing) Then
    Do i2=1,n_n+delta_n_rp
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
       Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do
    End Do

    Do i2=1,n_c+delta_c_rp
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = -id_u(i4)
       id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do
    End Do 

    Do i2=1,6
     Do i3=1,3
      Fnames(i_zaehl) = "sup_"//Bu(i2)//"^* W+ "//Trim(c_d(i3))
      Fnames(i_zaehl+1) = "sup_"//Bu(i2)//" W- "//Trim(c_d(i3))//"^*"
      Lnames(i_zaehl) = Trim(c_su(i2))//"^* W+ "//Trim(c_d(i3))
      Lnames(i_zaehl+1) = Trim(c_su(i2))//" W- "//Trim(c_d(i3))//"^*"
      id_d_3(i_zaehl-i_c2,1) = -id_su(i2)
      id_d_3(i_zaehl-i_c2,2) = id_W
      id_d_3(i_zaehl-i_c2,3) = id_d(i3)
      id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
      i_zaehl = i_zaehl+2
     End Do
    End Do  

   Else ! GenerationMixing
    Do i2=1,n_n+delta_n_rp
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
      id_d_3(i_zaehl-i_c2,1) = id_n(i2)
      id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
      id_d_3(i_zaehl-i_c2,3) = id_u(i4)
      i_zaehl = i_zaehl+1    
     End Do
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
      Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
      id_d_3(i_zaehl-i_c2,1) = id_n(i2)
      id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = id_d(i4)
      i_zaehl = i_zaehl+1    
     End Do
    End Do

    Do i2=1,n_c+delta_c_rp
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      id_d_3(i_zaehl-i_c2,1) = id_c(i2)
      id_d_3(i_zaehl-i_c2,2) = id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = -id_u(i4)
      id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
      i_zaehl = i_zaehl+2
     End Do
    End Do 

    Do i2=5,6
     i3=3
     Fnames(i_zaehl) = "stop_"//Bu(i2-4)//"^* W+ "//Trim(c_d(i3))
     Fnames(i_zaehl+1) = "stop_"//Bu(i2-4)//" W- "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl) = Trim(c_su(i2))//"^* W+ "//Trim(c_d(i3))
     Lnames(i_zaehl+1) = Trim(c_su(i2))//" W- "//Trim(c_d(i3))//"^*"
     id_d_3(i_zaehl-i_c2,1) = -id_su(i2)
     id_d_3(i_zaehl-i_c2,2) = id_W
     id_d_3(i_zaehl-i_c2,3) = id_d(i3)
     id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
     i_zaehl = i_zaehl+2
    End Do  
   End If ! GenerationMixing

   Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
   Do i2=1,i_zaehl - i_c2 - 2
     If (BR_Glu3(i2).Gt.BrMin) Write(io_L,202) BR_glu3(i2),3,id_d_3(i2,:) &
             &                , Trim(c_m)//" -> "//Trim(Lnames(i2+i_c2))//")"
   End Do

   BR(i_c2+1:i_zaehl-2) = BR_Glu3(1:i_zaehl-i_c2-2)
   gP(i_c2+1:i_zaehl-2) = gP_glu3(1:i_zaehl-i_c2-2)

  End If

  Call WriteDecays2(io, " gluino_" , Fnames, gP, 100*BR, gT_Glu, BrMin100)


  !-----------------
  ! neutral scalars
  !-----------------
  Do i1=1,n_s0
   c_m = c_s0(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    Do i2 = 1, 3 - delta_c_rp
     Do i3 = i2, 3 - delta_c_rp
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Fnames(i_zaehl+1) = c_lp(i2)//" "//c_lm(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl+1) = Trim(c_lm(i2))//" "//Trim(c_lP(i3))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,3
     Do i3 = i2,3
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i2)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+1) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Fnames(i_zaehl+10) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+1) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+10) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i3)
       id_d_2(i_zaehl+10,:) = - id_d_2(i_zaehl+9,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 9
    Do i2 = 1, n_sl
     Do i3 = i2, n_sl
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'slepton^+_'//Bu(i2)//' slepton^-_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       Lnames(i_zaehl+1) = Trim(c_sle(i2))//" "//Trim(c_slep(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1, n_sn
     Do i3 = i2, n_sn
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'sneutrino_'//Bu(i2)//' sneutrino^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       id_d_2(i_zaehl,1) = -id_snu(i2)
       id_d_2(i_zaehl,2) = id_snu(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'sneutrino_'//Bu(i2)//' sneutrino^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'sneutrino^*_'//Bu(i2)//' sneutrino_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       id_d_2(i_zaehl,1) = -id_snu(i2)
       id_d_2(i_zaehl,2) = id_snu(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,6
     Do i3 = i2,6
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'd-squark^*_'//Bu(i2)//' d-squark_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       Lnames(i_zaehl+1) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Fnames(i_zaehl+37) = &
                    & 'u-squark^*_'//Bu(i2)//' u-squark_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       Lnames(i_zaehl+37) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       id_d_2(i_zaehl+37,:) = - id_d_2(i_zaehl+36,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 36

   Else ! .not.GenerationMixing
    ! leptons
    Do i2=1,3 - delta_c_rp
     Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = -id_l(i2)
     id_d_2(i_zaehl,2) = id_l(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = -id_d(i2)
     id_d_2(i_zaehl,2) = id_d(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! u-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_u(i2))
     id_d_2(i_zaehl,1) = -id_u(i2)
     id_d_2(i_zaehl,2) = id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! sleptons
    Do i2=1,3 - delta_c_rp
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        Lnames(i_zaehl+1) = Trim(c_sle(i3))//" "//Trim(c_slep(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! sneutrinos 
    Do i2=1,n_sn
     Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i2))
     id_d_2(i_zaehl,1) = -id_snu(i2)
     id_d_2(i_zaehl,2) = id_snu(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        Lnames(i_zaehl+1) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! u-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        Lnames(i_zaehl+1) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2

    If (delta_c_rp.Eq.1) Then
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     i2 = 2
    Else If (delta_c_rp.Eq.3) Then
     i2 = 0
    Else
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     Fnames(3) = 'taus '
     i2 = 3
    End If
    Fnames(i2+1) = 'd-quark '
    Fnames(i2+2) = 's-quark '
    Fnames(i2+3) = 'b-quark '
    Fnames(i2+4) = 'u-quark '
    Fnames(i2+5) = 'c-quark '
    Fnames(i2+6) = 't-quark '

    If (delta_c_rp.Eq.1) Then
     Fnames(9) = 'Selectron 1 1 '
     Fnames(10) = 'Selectron 1 2 '
     Fnames(11) = 'Selectron 2 1 '
     Fnames(12) = 'Selectron 2 2 '
     Fnames(13) = 'Smuon 1 1 '
     Fnames(14) = 'Smuon 1 2 '
     Fnames(15) = 'Smuon 2 1 '
     Fnames(16) = 'Smuon 2 2 '
     Fnames(17) = 'e-Sneutrino'
     Fnames(18) = 'mu-Sneutrino'
     i2 = 18
    Else If (delta_c_rp.Eq.3) Then
     i2 = 6
    Else
     Fnames(10) = 'Selectron 1 1 '
     Fnames(11) = 'Selectron 1 2 '
     Fnames(12) = 'Selectron 2 1 '
     Fnames(13) = 'Selectron 2 2 '
     Fnames(14) = 'Smuon 1 1 '
     Fnames(15) = 'Smuon 1 2 '
     Fnames(16) = 'Smuon 2 1 '
     Fnames(17) = 'Smuon 2 2 '
     Fnames(18) = 'Stau 1 1 '
     Fnames(19) = 'Stau 1 2 '
     Fnames(20) = 'Stau 2 1 '
     Fnames(21) = 'Stau 2 2 '
     Fnames(22) = 'e-Sneutrino'
     Fnames(23) = 'mu-Sneutrino'
     Fnames(24) = 'tau-Sneutrino'
     i2 = 24
    End If
    Fnames(i2+1) = 'Sdown 1 1 '
    Fnames(i2+2) = 'Sdown 1 2 '
    Fnames(i2+3) = 'Sdown 2 1 '
    Fnames(i2+4) = 'Sdown 2 2 '
    Fnames(i2+5) = 'S-strange 1 1 '
    Fnames(i2+6) = 'S-strange 1 2 '
    Fnames(i2+7) = 'S-strange 2 1 '
    Fnames(i2+8) = 'S-strange 2 2 '
    Fnames(i2+9) = 'Sbottom 1 1 '
    Fnames(i2+10) = 'Sbottom 1 2 '
    Fnames(i1+11) = 'Sbottom 2 1 '
    Fnames(i1+12) = 'Sbottom 2 2 '
    Fnames(i1+13) = 'Sup 1 1 '
    Fnames(i1+14) = 'Sup 1 2 '
    Fnames(i1+15) = 'Sup 2 1 '
    Fnames(i1+16) = 'Sup 2 2 '
    Fnames(i1+17) = 'S-charm 1 1 '
    Fnames(i1+18) = 'S-charm 1 2 '
    Fnames(i1+19) = 'S-charm 2 1 '
    Fnames(i1+20) = 'S-charm 2 2 '
    Fnames(i1+21) = 'Stop 1 1 '
    Fnames(i1+22) = 'Stop 1 2 '
    Fnames(i1+23) = 'Stop 2 1 '
    Fnames(i1+24) = 'Stop 2 2 '

   End If

   ! neutralinos 
   Do i2=1,n_n + delta_n_rp
    Do i3=i2,n_n + delta_n_rp
     Fnames(i_zaehl) = 'neutralino_'//Bu(i2)//' neutralino_'//Bu(i3)
     Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   ! charginos
   Do i2=1, n_c + delta_c_rp
    Do i3=i2, n_c + delta_c_rp
     If (i2.Eq.i3) Then
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i2))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      i_zaehl = i_zaehl + 1
     Else
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Fnames(i_zaehl+1) = 'chargino^-_'//Bu(i2)//' chargino^+_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i3))
      Lnames(i_zaehl+1) = Trim(c_cp(i2))//" "//Trim(c_cm(i3))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
      i_zaehl = i_zaehl + 2
     End If
    End Do
   End Do

   Fnames(i_zaehl) = 'Z0 Z0'
   Lnames(i_zaehl) = 'Z0 Z0'
   id_d_2(i_zaehl,1) = id_Z
   id_d_2(i_zaehl,2) = id_Z
   i_zaehl = i_zaehl + 1

   Fnames(i_zaehl) = 'W+ W-'
   Lnames(i_zaehl) = 'W- W+'
   id_d_2(i_zaehl,1) = id_W
   id_d_2(i_zaehl,2) = - id_W
   i_zaehl = i_zaehl + 1

   Do i2=1,n_p0
    Fnames(i_zaehl) = 'Z0 '//Trim(c_p0(i2))
    Lnames(i_zaehl) = 'Z0 '//Trim(c_p0(i2))
    id_d_2(i_zaehl,1) = id_Z
    id_d_2(i_zaehl,2) = id_P0(i2)
    i_zaehl = i_zaehl + 1
   End Do
   Do i2=1,n_p0
    Do i3=i2,n_p0
     Fnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_P0(i2)
     id_d_2(i_zaehl,2) = id_P0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i2=1,i1-1
    Do i3=i2,i1-1
     Fnames(i_zaehl) = Trim(c_s0(i2))//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) = Trim(c_s0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_S0(i2)
     id_d_2(i_zaehl,2) = id_S0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Fnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    Lnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Lnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_W
    id_d_2(i_zaehl,2) = id_sm(i2)
    id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl + 2
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    Lnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_sp(i2)
    id_d_2(i_zaehl,2) = id_sm(i2)
    i_zaehl = i_zaehl + 1
    Do i3=i2+1,n_spm
     Fnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_sp(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Fnames(i_zaehl) = "g g"
   Fnames(i_zaehl + 1) = Trim(c_phot)//" "//Trim(c_phot)
   Lnames(i_zaehl) = Fnames(i_zaehl)
   Lnames(i_zaehl+1) = Fnames(i_zaehl+1)
   id_d_2(i_zaehl,:) = id_gl
   id_d_2(i_zaehl+1,:) = id_ph
   i_zaehl = i_zaehl + 2
   
   If ( HighScaleModel.Eq."NMSSM") Then
    gT = gT_S03(i1)
    BR(1:200) = BR_S03(i1,:)
    gP(1:200) = gP_S03(i1,:)
   Else If ( HighScaleModel.Eq."RPexplicit") Then
    gT = gT_S05(i1)
    BR(1:200) = BR_S05(i1,:)
    gP(1:200) = gP_S05(i1,:)
   Else
    gT = gT_S0(i1)
    BR(1:200) = BR_S0(i1,:)
    gP(1:200) = gP_S0(i1,:)
   End If    

   Write(io_L,200) id_S0(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    
   
   Fnames(i_zaehl) = "W^+ W^-*" 
   Fnames(i_zaehl+1) = "W^- W^+*"
   Fnames(i_zaehl+2) = "Z Z^*" 

   Call WriteDecays2(io,c_m , Fnames, gP, 100*BR, gT, BrMin100)
   If ((BR(i_zaehl).Gt.BrMin).Or.(BR(i_zaehl+2).Gt.BrMin)) Then
    Write(io_L,100) "# writing decays into V V* as 3-body decays"
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    If (BR(i_zaehl).Gt.BrMin) Then
     Do i2=1,3
      Write(io_L,202) BrWln(i2) * BR(i_zaehl),3,id_w,id_l(i2),-id_nu(i2) &
      & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_lm(i2))//" "//Trim(c_nu(i2))//")"
      Write(io_L,202) BrWln(i2) * BR(i_zaehl),3,-id_w,-id_l(i2),id_nu(i2) &
      & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_lp(i2))//" "//Trim(c_nu(i2))//")"
     End Do
     BRtot = 0.5_dp * Sum(BrWqq)
     Do i2=1,3
      Write(io_L,202) Brtot * Abs(CKM(1,i2))**2 * BR(i_zaehl) &
        & , 3, id_w,-id_u(1),id_d(i2) &
        & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_u(1))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(1,i2))**2 * BR(i_zaehl) &
        & ,3,-id_w,id_u(1),-id_d(i2) &
        & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_u(1))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(2,i2))**2 * BR(i_zaehl) &
        & , 3, id_w,-id_u(2),id_d(i2) &
        & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_u(2))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(2,i2))**2 * BR(i_zaehl) &
       & ,3,-id_w,id_u(2),-id_d(i2) &
       & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_u(2))//" "//Trim(c_d(i2))//")"
     End Do
    End If
    i_zaehl = i_zaehl + 2
    If (BR(i_zaehl).Gt.BrMin) Then 
     Write(io_L,202) BrZinv*BR(i_zaehl),3,id_Z,id_nu(1),-id_nu(1), &
        & Trim(c_m)//" -> Z0 nu_i nu_i)"
     Do i2=1,3  
      Write(io_L,202) BrZll(i2)*BR(i_zaehl),3,id_Z,id_l(i2),-id_l(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_lm(i2))//" "//Trim(c_lp(i2))//")"
     End Do
     Do i2=1,3  
      Write(io_L,202) BrZqq(i2)*BR(i_zaehl),3,id_Z,id_d(i2),-id_d(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_d(i2))//" "//Trim(c_d(i2))//")"
     End Do
     Do i2=1,2 
      Write(io_L,202) BrZqq(i2)*BR(i_zaehl),3,id_Z,id_u(i2),-id_u(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_u(i2))//" "//Trim(c_u(i2))//")"
     End Do
    End If
   End If

  End Do


  !----------------------
  ! neutral pseudoscalar
  !----------------------
  Do i1=1,n_p0
   c_m = c_p0(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    Do i2 = 1, 3 - delta_c_rp
     Do i3 = i2, 3 - delta_c_rp
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Fnames(i_zaehl+1) = c_lp(i2)//" "//c_lm(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl+1) = Trim(c_lm(i2))//" "//Trim(c_lP(i3))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,3
     Do i3 = i2,3
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i2)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+1) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Fnames(i_zaehl+10) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+1) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+10) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i3)
       id_d_2(i_zaehl+10,:) = - id_d_2(i_zaehl+9,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 9
    Do i2 = 1, n_sl
     Do i3 = i2, n_sl
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'slepton^+_'//Bu(i2)//' slepton^-_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       Lnames(i_zaehl+1) = Trim(c_sle(i2))//" "//Trim(c_slep(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,6
     Do i3 = i2,6
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'd-squark^*_'//Bu(i2)//' d-squark_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       Lnames(i_zaehl+1) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Fnames(i_zaehl+37) = &
                    & 'u-squark^*_'//Bu(i2)//' u-squark_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       Lnames(i_zaehl+37) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       id_d_2(i_zaehl+37,:) = - id_d_2(i_zaehl+36,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 36

   Else ! .not.GenerationMixing
    ! leptons
    Do i2=1,3 - delta_c_rp
     Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = -id_l(i2)
     id_d_2(i_zaehl,2) = id_l(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = -id_d(i2)
     id_d_2(i_zaehl,2) = id_d(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! u-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_u(i2))
     id_d_2(i_zaehl,1) = -id_u(i2)
     id_d_2(i_zaehl,2) = id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! sleptons
    Do i2=1,3 - delta_c_rp
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        Lnames(i_zaehl+1) = Trim(c_sle(i3))//" "//Trim(c_slep(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! d-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        Lnames(i_zaehl+1) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! u-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        Lnames(i_zaehl+1) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2

    If (delta_c_rp.Eq.1) Then
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     i2 = 2
    Else If (delta_c_rp.Eq.3) Then
     i2 = 0
    Else
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     Fnames(3) = 'taus '
     i2 = 3
    End If
    Fnames(i2+1) = 'd-quark '
    Fnames(i2+2) = 's-quark '
    Fnames(i2+3) = 'b-quark '
    Fnames(i2+4) = 'u-quark '
    Fnames(i2+5) = 'c-quark '
    Fnames(i2+6) = 't-quark '

    If (delta_c_rp.Eq.1) Then
     Fnames(9) = 'Selectron 1 1 '
     Fnames(10) = 'Selectron 1 2 '
     Fnames(11) = 'Selectron 2 1 '
     Fnames(12) = 'Selectron 2 2 '
     Fnames(13) = 'Smuon 1 1 '
     Fnames(14) = 'Smuon 1 2 '
     Fnames(15) = 'Smuon 2 1 '
     Fnames(16) = 'Smuon 2 2 '
     i2 = 16
    Else If (delta_c_rp.Eq.3) Then
     i2 = 6
    Else
     Fnames(10) = 'Selectron 1 1 '
     Fnames(11) = 'Selectron 1 2 '
     Fnames(12) = 'Selectron 2 1 '
     Fnames(13) = 'Selectron 2 2 '
     Fnames(14) = 'Smuon 1 1 '
     Fnames(15) = 'Smuon 1 2 '
     Fnames(16) = 'Smuon 2 1 '
     Fnames(17) = 'Smuon 2 2 '
     Fnames(18) = 'Stau 1 1 '
     Fnames(19) = 'Stau 1 2 '
     Fnames(20) = 'Stau 2 1 '
     Fnames(21) = 'Stau 2 2 '
     i2 = 21
    End If
    Fnames(i2+1) = 'Sdown 1 1 '
    Fnames(i2+2) = 'Sdown 1 2 '
    Fnames(i2+3) = 'Sdown 2 1 '
    Fnames(i2+4) = 'Sdown 2 2 '
    Fnames(i2+5) = 'S-strange 1 1 '
    Fnames(i2+6) = 'S-strange 1 2 '
    Fnames(i2+7) = 'S-strange 2 1 '
    Fnames(i2+8) = 'S-strange 2 2 '
    Fnames(i2+9) = 'Sbottom 1 1 '
    Fnames(i2+10) = 'Sbottom 1 2 '
    Fnames(i1+11) = 'Sbottom 2 1 '
    Fnames(i1+12) = 'Sbottom 2 2 '
    Fnames(i1+13) = 'Sup 1 1 '
    Fnames(i1+14) = 'Sup 1 2 '
    Fnames(i1+15) = 'Sup 2 1 '
    Fnames(i1+16) = 'Sup 2 2 '
    Fnames(i1+17) = 'S-charm 1 1 '
    Fnames(i1+18) = 'S-charm 1 2 '
    Fnames(i1+19) = 'S-charm 2 1 '
    Fnames(i1+20) = 'S-charm 2 2 '
    Fnames(i1+21) = 'Stop 1 1 '
    Fnames(i1+22) = 'Stop 1 2 '
    Fnames(i1+23) = 'Stop 2 1 '
    Fnames(i1+24) = 'Stop 2 2 '

   End If

   ! neutralinos 
   Do i2=1,n_n + delta_n_rp
    Do i3=i2,n_n + delta_n_rp
     Fnames(i_zaehl) = 'neutralino_'//Bu(i2)//' neutralino_'//Bu(i3)
     Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   ! charginos
   Do i2=1, n_c + delta_c_rp
    Do i3=i2, n_c + delta_c_rp
     If (i2.Eq.i3) Then
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i2))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      i_zaehl = i_zaehl + 1
     Else
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Fnames(i_zaehl+1) = 'chargino^-_'//Bu(i2)//' chargino^+_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i3))
      Lnames(i_zaehl+1) = Trim(c_cp(i2))//" "//Trim(c_cm(i3))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
      i_zaehl = i_zaehl + 2
     End If
    End Do
   End Do

   Do i2=1,n_spm
    Fnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Fnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    Lnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Lnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_W
    id_d_2(i_zaehl,2) = id_sm(i2)
    id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl + 2
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    Lnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_sp(i2)
    id_d_2(i_zaehl,2) = id_sm(i2)
    i_zaehl = i_zaehl + 1
    Do i3=i2+1,n_spm
     Fnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_sp(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl + 2
    End Do
   End Do

   Do i2=1,n_s0
    Fnames(i_zaehl) = 'Z0 '//Trim(c_s0(i2))
    Lnames(i_zaehl) = 'Z0 '//Trim(c_s0(i2))
    id_d_2(i_zaehl,1) = id_Z
    id_d_2(i_zaehl,2) = id_s0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   Do i2=1,i1-1
    Do i3=1,n_s0
     Fnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_P0(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do

   If ( HighScaleModel.Eq."NMSSM") Then
    gT = gT_P03(i1+1)
    BR(1:200) = BR_P03(i1+1,:)
    gP(1:200) = gP_P03(i1+1,:)
   Else If ( HighScaleModel.Eq."RPexplicit") Then
    gT = gT_P05(i1+1)
    BR(1:200) = BR_P05(i1+1,:)
    gP(1:200) = gP_P05(i1+1,:)
   Else
    gT = gT_P0(i1+1)
    BR(1:200) = BR_P0(i1+1,:)
    gP(1:200) = gP_P0(i1+1,:)
   End If    

   Write(io_L,200) id_p0(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    
   If (BR(i_zaehl+3).Gt.BrMin)  Write(io_L,201) BR(i_zaehl+3),2,id_gl  &
           &  ,id_gl,Trim(c_m)//" -> g g)"
 
   Call WriteDecays2(io,c_m , Fnames, gP, 100*BR, gT, BrMin100)

  End Do


  !-----------------
  ! charged scalars
  !-----------------
  Do i1=1,n_spm
   c_m = c_sm(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    !-----------------------------------------------------------------
    ! leptons, if RP is violated, these BRs are included in the
    ! chargino/neutralino final  states
    !-----------------------------------------------------------------
    Do i2=1,3 - delta_c_rp
     Lnames(i_zaehl) = Trim(c_nu(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = id_l(i2)
     id_d_2(i_zaehl,2) = -id_nu(i2)
     i_zaehl = i_zaehl + 1
    End Do

    !--------
    ! quarks
    !--------
    Do i2 = 1, 3
     Do i3 = 1,3
      Lnames(i_zaehl) = Trim(c_u(i3))//" "//Trim(c_d(i2))
      id_d_2(i_zaehl,1) = id_d(i2)
      id_d_2(i_zaehl,2) = -id_u(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do

    !-----------------------------------------------------------------
    ! sleptons, in case of RP violation -> S0 S-, P0 S- final states
    !-----------------------------------------------------------------
    Do i2 = 1,2*(3-delta_c_rp)
     Do i3 = 1,3-delta_c_rp
      Lnames(i_zaehl) = Trim(c_sle(i2))//" "//Trim(c_snu(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = -id_snu(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do
    !------------------
    ! into squarks
    !------------------
    Do i2 = 1,6
     Do i3 = 1,6
      Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_su(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_su(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do

   Else ! no Generation Mixing
    Do i2=1,3 - delta_c_rp
     Lnames(i_zaehl) = Trim(c_nu(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = id_l(i2)
     id_d_2(i_zaehl,2) = -id_nu(i2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = id_d(i2)
     id_d_2(i_zaehl,2) = -id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1, 2*(3 - delta_c_rp)
     Lnames(i_zaehl) = Trim(c_sle(i2))//" "//Trim(c_snu((i2+1)/2))
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = -id_snu((i2+1)/2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1,6
     If (i2.Le.2) i3=1
     If (i2.Le.4) i3=3
     If (i2.Le.6) i3=5
     Do i4=i3,i3+1
      Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_su(i4))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_su(i4)
      i_zaehl = i_zaehl + 1
     End Do
    End Do
   End If

   Do i2=1,n_c+delta_c_rp
    Do i3=1,n_n+delta_n_rp
     Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = - id_c(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do

   Do i2=1,n_p0
    Lnames(i_zaehl) = " W- "//Trim(c_p0(i2))
    id_d_2(i_zaehl,1) = -id_W
    id_d_2(i_zaehl,2) = id_p0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   Do i2=1,n_s0
    Lnames(i_zaehl) = " W- "//Trim(c_s0(i2))
    id_d_2(i_zaehl,1) = -id_W
    id_d_2(i_zaehl,2) = id_s0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   Do i3=1,3
    Lnames(177+i3) = Trim(c_grav)//" "//Trim(c_lm(i3))
    id_d_2(177+i3,1) = id_grav
    id_d_2(177+i3,2) = id_l(i3)
   End Do

   If ( HighScaleModel.Eq."NMSSM") Then
    gT = gT_Spm(i1+1)
    BR(1:200) = BR_Spm(i1+1,:)
    gP(1:200) = gP_Spm(i1+1,:)
   Else If ( HighScaleModel.Eq."RPexplicit") Then
    gT = gT_Spm8(i1+1)
    BR(1:200) = BR_Spm8(i1+1,:)
    gP(1:200) = gP_Spm8(i1+1,:)
   Else
    gT = gT_Spm(i1+1)
    BR(1:200) = BR_Spm(i1+1,:)
    gP(1:200) = gP_Spm(i1+1,:)
   End If    

   Write(io_L,200) id_sm(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    
   Do i2=178,180
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    

   If (Maxval(BR(181:200)).Gt.BrMin) Then
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    i_zaehl= 181
    Do i2=1,3
     Do i3=i2,3
      If (BR(i_zaehl).Gt.BrMin)  Write(io_L,202) BR(i_zaehl),3,id_l(i2),-id_l(i3),id_sm(1) &
      & ,Trim(c_m)//" ->  "//Trim(c_lm(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_sm(1))//")"
       If ((i2.Ne.i3).And.(BR(i_zaehl).Gt.BrMin) ) &
      & Write(io_L,202) BR(i_zaehl+1),3,id_l(i3),-id_l(i2),id_sm(1) &
      & ,Trim(c_m)//" ->  "//Trim(c_lm(i3))//" "//Trim(c_lp(i2))//" "//Trim(c_sm(1))//")"
       i_zaehl = i_zaehl + 1
       If (i2.Ne.i3) i_zaehl = i_zaehl + 1
      End Do
     End Do

    i_zaehl= 191
    Do i2=1,3
     Do i3=i2,3
      If (BR(i_zaehl).Gt.BrMin)  Write(io_L,202) BR(i_zaehl),3,id_l(i2),id_l(i3),id_sp(1) &
      & ,Trim(c_m)//" ->  "//Trim(c_lm(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_sp(1))//")"
       i_zaehl = i_zaehl + 1
      End Do
     End Do
   End If

  End Do

 End If ! L_BR  

  If (L_CS) Then
   Write(io_L,100) "Block SPhenoCrossSections  # cross sections"
   Do i3=1,Size(Ecms)
    If (Ecms(i3).Eq.0._dp) Exit
   If (ISR(i3)) Then
    Write(io_L,4712) Ecms(i3),Pm(i3),Pp(i3)," 1  # e+ e- XS, Pe-, Pe+,  including ISR"
   Else
    Write(io_L,4712) Ecms(i3),Pm(i3),Pp(i3)," 0  # e+ e- XS"
   End If  
   Write(io_L,100) "#      Sigma [fb]    NDA        ID1     ID2"
   ! u-squarks
   Do i1=1,6
    Do i2=1,6
     If (SigSup(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigSup(i3,i1,i2), 2 &
        &      , id_su(i1), -id_su(i2), c_su(i1)//" "//Trim(c_su(i2))//"*"
    End Do
   End Do
   ! d-squarks
   Do i1=1,6
    Do i2=1,6
     If (SigSdown(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigSdown(i3,i1,i2), 2 &
        &      , id_sd(i1), -id_sd(i2), c_sd(i1)//" "//Trim(c_sd(i2))//"*"
    End Do
   End Do
   If (HighScaleModel.Eq."RPexplicit") Then
   Else
    ! sleptons
    Do i1=1,6
     Do i2=1,6
      If (SigSle(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigSle(i3,i1,i2), 2 &
        &      , id_sle(i1), -id_sle(i2), c_sle(i1)//" "//c_slep(i2)
     End Do
    End Do
    ! sneutrinos
    Do i1=1,3
     If (SigSn(i3,i1,i1).Gt.SigMin) Write(io_L,4711) SigSn(i3,i1,i1), 2 &
        &      , id_snu(i1), -id_snu(i1), c_snu(i1)//" "//Trim(c_snu(i1))//"*"
    End Do
   End If
   ! neutralinos
   Do i1=1,n_n + delta_n_rp
    Do i2=i1,n_n + delta_n_rp
     If (SigChi0(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigChi0(i3,i1,i2), 2 &
        &      , id_n(i1), id_n(i2), c_c0(i1)//" "//c_c0(i2)
    End Do
   End Do
   ! charginos
   Do i1=1,n_c + delta_c_rp
    Do i2=1,n_c + delta_c_rp
     If (SigC(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigC(i3,i1,i2), 2 &
        &      , id_c(i1), -id_c(i2), c_cm(i1)//" "//c_cp(i2)
    End Do
   End Do
   
   If (HighScaleModel.Eq."RPexplicit") Then
    Do i1=1,n_s0
     If (SigS0(i3,i1).Gt.SigMin) Write(io_L,4711) SigS0(i3,i1),2, id_S0(i1) &
                                                & , id_Z , Trim(c_s0(i1))//" Z"
    End Do
    Do i1=1,n_s0
     Do i2=1,n_p0
      If (SigSP(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigSP(i3,i1,i2),2 &
                & , id_S0(i1), id_p0(i2) , Trim(c_s0(i1))//" "//Trim(c_p0(i2))
     End Do
    End Do
    Do i1=1,n_spm
     Do i2=1,n_spm
      If (SigHP(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigHP(i3,i1,i2),2 &
                & , id_Sp(i1), id_sm(i2) , Trim(c_sp(i1))//" "//Trim(c_sm(i2))
     End Do
    End Do

   Else If (HighScaleModel.Eq."NMSSM") Then
    Do i1=1,n_s0
     If (SigS0(i3,i1).Gt.SigMin) Write(io_L,4711) SigS0(i3,i1),2, id_S0(i1) &
                                                & , id_Z , Trim(c_s0(i1))//" Z"
    End Do
    Do i1=1,n_s0
     Do i2=1,n_p0
      If (SigSP(i3,i1,i2).Gt.SigMin) Write(io_L,4711) SigSP(i3,i1,i2),2 &
                & , id_S0(i1), id_p0(i2) , Trim(c_s0(i1))//" "//Trim(c_p0(i2))
     End Do
    End Do
    If (SigHP(i3,1,1).Gt.SigMin) &
       & Write(io_L,4711) SigHP(i3,1,1), 2, id_Hp, -id_Hp , "H+ H-"

   Else
    If (SigS0(i3,1).Gt.SigMin) Write(io_L,4711) SigS0(i3,1),2, id_S0(1), id_Z &
                                                & , "h0 Z"
    If (SigS0(i3,2).Gt.SigMin) &
      & Write(io_L,4711) SigS0(i3,2), 2, id_S0(2), id_Z, "H0 Z"
    If (SigSP(i3,1,1).Gt.SigMin) &
      & Write(io_L,4711) SigSP(i3,1,1),2, id_S0(1), id_A0, "h0 A0"
    If (SigSP(i3,2,1).Gt.SigMin) &
       &  Write(io_L,4711) SigSP(i3,2,1),2, id_S0(2), id_A0 , "H0 A0"
    If (SigHP(i3,1,1).Gt.SigMin) &
       & Write(io_L,4711) SigHP(i3,1,1), 2, id_Hp, -id_Hp , "H+ H-"
    End If
   End Do
  End If ! L_CS

  Write(io,*) " "
  Write(io,*) "Low energy constraints "
  Write(io,5410) " 10^4 BR(b -> s gamma) :",BRbtosgamma
  Write(io,5410) " BR(b -> s mu+ mu-)    :", BrBToSLL
  Write(io,5410) " BR(b -> s nu nu)      :", BToSNuNu
  Write(io,5410) " BR(B_s -> mu+ mu)     :",bs_mumu
  Write(io,5410) " BR(B_u -> tau nu)     :", BR_Bu_TauNu
  Write(io,5411) " Delta(M_B_d)          :",deltambd
  Write(io,5411) " Delta(M_B_s)          :",deltambs
  Write(io,5410) " Delta(a_e)            :",a_e
  Write(io,5410) " Delta(a_mu)           :",a_mu
  Write(io,5410) " Delta(a_tau)          :",a_tau
  Write(io,5410) " D_e                   :",d_e
  Write(io,5410) " D_mu                  :",d_mu
  Write(io,5410) " D_tau                 :",d_tau
  Write(io,5410) " Br(mu -> e gamma)     :",BrMutoEGamma
  Write(io,5410) " Br(mu -> 3 e)         :",BrMu3e
  Write(io,5410) " Br(tau -> e gamma)    :",BrTautoEGamma
  Write(io,5410) " Br(tau -> mu gamma)   :",BrTautoMuGamma
  Write(io,5410) " Br(tau -> 3 e)        :",BrTau3e
  Write(io,5410) " Br(tau -> 3 mu)       :",BrTau3mu
  Write(io,5410) " rho parameter         :",rho_parameter
  If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
   Write(io,5410) " m_nu_1 (in eV)       :",mf_nu(1)*1.e9_dp
   Write(io,5410) " m_nu_2 (in eV)       :",mf_nu(2)*1.e9_dp
   Write(io,5410) " m_nu_3 (in eV)       :",mf_nu(3)*1.e9_dp
  End If
  Write(io,*) " "

  Write(io_L,100) "Block SPhenoLowEnergy  # low energy observables"
  Write(io_L,101) 1,1.e-4*BrBToSGamma," # BR(b -> s gamma)"
  Write(io_L,101) 2,BrBToSLL," # BR(b -> s mu+ mu-)"
  Write(io_L,101) 3,BToSNuNu," # BR(b -> s nu nu)"
  Write(io_L,101) 4,Bs_mumu," # BR(Bs -> mu+ mu-)"
  Write(io_L,101) 5,BR_Bu_TauNu," # BR(B_u -> tau nu)"
  Write(io_L,101) 6,Abs(DeltaMbd)," # |Delta(M_Bd)| [ps^-1] "
  Write(io_L,101) 7,Abs(DeltaMbs)," # |Delta(M_Bs)| [ps^-1] "

  Write(io_L,101) 10,A_e," # Delta(g-2)_electron"
  Write(io_L,101) 11,A_mu," # Delta(g-2)_muon"
  Write(io_L,101) 12,A_tau," # Delta(g-2)_tau"
  Write(io_L,101) 13,d_e," # electric dipole moment of the electron"
  Write(io_L,101) 14,d_mu," # electric dipole moment of the muon"
  Write(io_L,101) 15,d_tau," # electric dipole moment of the tau"
  Write(io_L,101) 16,BrMutoEGamma," # Br(mu -> e gamma)"
  Write(io_L,101) 17,BrTautoEGamma," # Br(tau -> e gamma)"
  Write(io_L,101) 18,BrTautoMuGamma," # Br(tau -> mu gamma)"
  Write(io_L,101) 19,BrMu3E," # Br(mu -> 3 e)"
  Write(io_L,101) 20,BrTau3E," # Br(tau -> 3 e)"
  Write(io_L,101) 21,BrTau3Mu," # Br(tau -> 3 mu)"

  Write(io_L,101) 30,rho_parameter," # Delta(rho_parameter)"
  Write(io_L,101) 40,BR_Z_e_mu," # BR(Z -> e mu)"
  Write(io_L,101) 41,BR_Z_e_tau," # BR(Z -> e tau)"
  Write(io_L,101) 42,BR_Z_mu_tau," # BR(Z -> mu tau)"
  If (Maxval(Abs(MnuL5)).Gt.0._dp) Then
   Write(io_L,100) "Block NuMass # neutrino masses"
   Write(io_L,101) 1,mf_nu(1)*1.e9_dp," # m_nu_1 in eV"
   Write(io_L,101) 2,mf_nu(2)*1.e9_dp," # m_nu_2 in eV"
   Write(io_L,101) 3,mf_nu(3)*1.e9_dp," # m_nu_3 in eV"
   Write(io_L,101) 4,(mf_nu(3)**2-mf_nu(2)**2)*1.e18_dp &
                & ," # Delta(m^2_atm) in eV^2"
   Write(io_L,101) 5,(mf_nu(2)**2-mf_nu(1)**2)*1.e18_dp &
                & ," # Delta(m^2_sol) in eV^2"
!   Call WriteMatrixBlockC2(io_L, 3, Unu, "Unu" &
!                         &, "neutrino mixing matrix", "U_(nu,")
  End If
  
  If (Present(omega)) Then
   Write(io_L,100) "Block Omega # omega h^2"
   Write(io_L,101) 1,omega," # omega h^2"   
  End If

  If (LWrite_LHC_Observables) Then
   Do i1=1,6
    if (id_sle(i1).eq.1000011) then
     mSle(2) = mSlepton(i1)
    else if (id_sle(i1).eq.2000011) then
     mSle(1) = mSlepton(i1)
    else if (id_sle(i1).eq.1000013) then
     mSle(4) = mSlepton(i1)
    else if (id_sle(i1).eq.2000013) then
     mSle(3) = mSlepton(i1)
    else if (id_sle(i1).eq.1000015) then
     mSle(5) = mSlepton(i1)
    else if (id_sle(i1).eq.2000015) then
     mSle(6) = mSlepton(i1)
    end if

    if (id_su(i1).eq.1000002) then
     mSu(2) = mSup(i1)
    else if (id_su(i1).eq.2000002) then
     mSu(1) = mSup(i1)
    else if (id_su(i1).eq.1000004) then
     mSu(4) = mSup(i1)
    else if (id_su(i1).eq.2000004) then
     mSu(3) = mSup(i1)
    else if (id_su(i1).eq.1000006) then
     mSu(5) = mSup(i1)
    else if (id_su(i1).eq.2000006) then
     mSu(6) = mSup(i1)
    end if

    if (id_sd(i1).eq.1000001) then
     mSd(2) = mSdown(i1)
    else if (id_sd(i1).eq.2000001) then
     mSd(1) = mSdown(i1)
    else if (id_sd(i1).eq.1000003) then
     mSd(4) = mSdown(i1)
    else if (id_sd(i1).eq.2000003) then
     mSd(3) = mSdown(i1)
    else if (id_sd(i1).eq.1000005) then
     mSd(5) = mSdown(i1)
    else if (id_sd(i1).eq.2000005) then
     mSd(6) = mSdown(i1)
    end if
   end do

   call Calc_LHC_observables(mN, N, mC, U, V, mSle, Rslepton        &
      & , mSd, RSdown, mSu, RSup, mGlu, PhaseGlu, mS0, RS0, mP0, RP0     &
      & , mSpm, RSpm, gauge, Y_u, Y_d, A_u, A_d, mu, vevSM, .False. & ! GenerationMixing &
      & , LHC_observ)

   Write(io_L,100) "Block LHCobservables # edge observables for LHC"
   i_zaehl = 1
   Do i1=2,n_n
    Do i2=1,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # e+e- edge with right selectron, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=2,n_n
    Do i2=1,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # e+e- edge with left selectron, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=2,n_n
    Do i2=1,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # mu+mu- edge with right smuon, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=2,n_n
    Do i2=1,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # mu+mu- edge with left smuon, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=2,n_n
    Do i2=1,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # tau+tau- edge with lighter stau, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i1=2,n_n
    Do i2=1,i1-1
     Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # tau+tau- edge with heavier stau, chi^0_"//Bu(i1)//", chi^0_"//Bu(i2)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+ l- q edge, averaging over d_L, s_L, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+ l- q threshold, averaging over d_L, s_l, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+-_near q edge, averaging over d_L, s_l, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+-_far q edge, averaging over d_L, s_l, u_L, c_L"
   i_zaehl = i_zaehl + 1
   Write(io_L,101) i_zaehl,LHC_observ(i_zaehl) &
       & ," # l+ l- b threshold"
   i_zaehl = i_zaehl + 1


  End If
!  Close(io_L)

 99 Format(1x,i5,3x,a)
100 Format(a)
101 Format(2x,i3,2x,1P,e16.8,2x,a) 
102 Format(1x,i9,3x,1P,e16.8,2x,a)
103 Format(a13,1P,e16.8,2x,a)
104 Format(i4,2x,1P,e16.8,2x,a)
105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)
107 Format(2i3,3x,1P,e16.8,3x,a)
108 Format(9x,1P,E16.8,0P,3x,a)
109 Format(1x,3i3,3x,1P,e16.8,3x,a)
110 Format(3x,2i3,3x,"# ",a)
200 Format("DECAY",1x,I9,3x,1P,E16.8,0P,3x,"# ",a)
201 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i9,1x),2x,"# BR(",a)
202 Format(3x,1P,e16.8,0p,3x,I2,3x,3(i9,1x),2x,"# BR(",a)

4711 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i9,1x),2x," # ",A)
4712 Format("XS 11 -11 ",F7.1," ",F5.2," ",F5.2," ",A)

!--------------------------------
! from WriteOutput
!--------------------------------
5101 Format(a,f15.7)
5103 Format(1P,3e16.8)
5410 Format(a25,1p,e16.7)
5411 Format(a25,1p,"(",e16.7,",",e16.7,")")
5807 Format(a9,f16.7,a5)
5808 Format("  Degree of polarization: P_e- =",f9.6," P_e+ =",f9.6)

 End Subroutine LesHouches_Out


 Subroutine LesHouches_out_start(io_L, io, kont, use_rge, Q)
 Implicit None
  Integer, Intent(in) :: io_L, io, kont, use_rge
  Real(dp), Intent(out) :: Q

  Character(len=8) :: Datum
  Character(len=10) :: Zeit
  Logical, Save :: l_open = .True. ! in case of a loop I want to open the
                                   ! output only once 
  Integer :: i_errors(1100)

  Q = Sqrt( GetRenormalizationScale() )

  Call Date_and_time(datum,zeit)
  If (l_open) Then
   Open(io_L,file="SPheno.spc",status="unknown")
   l_open = .False.
  End If
  !--------------------------------------------------------
  ! General information
  !--------------------------------------------------------
  ! Les Houches standard
  !-----------------------
   Write(io_L,100) "# SUSY Les Houches Accord 2.beta - MSSM spectrum + Decays"
   Write(io_L,100) "# SPheno "//version
   Write(io_L,100) &
     & "# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101"
   Write(io_L,100) &
     & "# in case of problems send email to porod@physik.uni-wuerzburg.de"
   Write(io_L,100) "# Created: "//Datum(7:8)//"."//Datum(5:6)//"."//Datum(1:4) &
     & //",  "//Zeit(1:2)//":"//Zeit(3:4)
   Write(io_L,100) "Block SPINFO         # Program information"
   Write(io_L,100) "     1   SPheno      # spectrum calculator"
   Write(io_L,100) "     2   "//version//"    # version number"
   !-----------------------------------------------
   ! check if somewhere a problem has had happened
   !-----------------------------------------------
   Call GetError(i_errors)
   !--------------------------------------
   ! a numerical problem might have happen
   !--------------------------------------
   If ((i_errors(1)+i_errors(3)+i_errors(5)+i_errors(7)+i_errors(8) &
     & + i_errors(10) + i_errors(12)+ Sum(i_errors(14:19))).Gt.0)   &
     & Write(io_L,100) &
 & "     3               # potential numerical problem, check file Messages.out"
   If (in_kont(1).Eq.1) Write(io_L,99) 3, &
    & "alpha(0) and alpha(mZ) have both been specified without check for"// &
    &  " consistency"
   If (in_kont(2).Eq.1) Write(io_L,99) 3, &
    & "redundant specification in Higgs sector"
   If (use_rge.Eq.0) Then
    Write(io_L,100) "#"
    Write(io_L,100) "Block SPhenoINFO     # SPheno specific information"
    If (TwoLoopRGE) Then
     Write(io_L,100) "    1      2         # using 2-loop RGEs"
    Else 
     Write(io_L,100) "    1      1         # using 1-loop RGEs"
    End If
    If (YukScen.Eq.1) Then
     Write(io_L,100) &
    &"    2      1         # using running masses for boundary conditions at mZ"
    Else
     Write(io_L,100) &
       &"    2      2         # using pole masses for boundary conditions at mZ"
    End If
   End If
   !------------------------------
   ! SPheno.out
   !------------------------------
  Write(io,123) "SPheno output file"
  Write(io,123) "Version "//version//" ,  "//"created: "//Datum(7:8)//"."// &
     & Datum(5:6)//"."//Datum(1:4)//",  "//Zeit(1:2)//":"//Zeit(3:4)
  Write(io,*) " "
123 Format(a)

  If (kont.Ne.0) Then
   Write(io,*) "There has been a problem during the run."
   Write(io,*) "Please check the file Messages.out for further information."
   Write(*,*) "There has been a problem during the run."
   Write(*,*) "Please check the file Messages.out for further information."
   Write(io,*) " "
  End If

  Write(io,*) " "
  If (YukScen.Eq.1) Then
   Write(io,*) &
     & "Running masses have been used for the boundary conditions at mZ" 
  Else If (YukScen.Eq.2) Then
   Write(io,*) "Pole masses have been used for the boundary conditions at mZ" 
  End If
  Write(io,*) "  Using Yukawa scheme :",GetYukawaScheme()
  Write(io,*) " "

   !---------------------------------------------------
   ! parameters + masses for SPheno.spc
   !---------------------------------------------------
   Write(io_L,100) "Block SMINPUTS  # SM parameters"
   Write(io_L,102) 1, 1._dp / alpha_MSbar(mZ, mW),"# alpha_em^-1(MZ)^MSbar"
   Write(io_L,102) 2,G_F,"# G_mu [GeV^-2]"
   Write(io_L,102) 3,alphaS_MZ,"# alpha_s(MZ)^MSbar"
   Write(io_L,102) 4,mZ,"# m_Z(pole)"
   Write(io_L,102) 5,mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) 6,mf_u(3),"# m_t(pole)"
   Write(io_L,102) 7,mf_l(3),"# m_tau(pole)"
   !------------------------
   ! input CKM
   !------------------------
   If (GenerationMixing) Then
    Write(io_L,100) "Block VCKMIN  # CKM matrix, Wolfenstein parameterization"
    Write(io_L,102) 1, lam_wolf,"# lambda"
    Write(io_L,102) 2,A_wolf,"# A"
    Write(io_L,102) 3,rho_wolf,"# rho bar"
    Write(io_L,102) 4,eta_wolf,"# eta bar"
   End If

  99 Format(1x,i5,3x,a)
  100 Format(a)
  102 Format(1x,i5,3x,1P,e16.8,2x,a)

 End Subroutine LesHouches_out_start

 Subroutine LH_Write_GY(io_L, Q, g, Yd, Yu, Yl, Ynu, YT)
 implicit none
  Integer, intent(in) :: io_L
  Real(dp), Intent(in) :: Q, g(3), Yd(3), Yu(3)
  Complex(dp), Optional, Intent(in) :: Yl(3,3), Ynu(3,3), YT(3,3)

  Integer :: i1, i2, ierr

  Write(io_L,106) "Block gauge Q=",Q,"# (SUSY scale)"
  Write(io_L,104) 1,g(1),"# g'(Q)^DRbar"
  Write(io_L,104) 2,g(2),"# g(Q)^DRbar"
  Write(io_L,104) 3,g(3),"# g3(Q)^DRbar"

  Write(io_L,106) "Block Yu Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yu(1), "# Y_u(Q)^DRbar"
  Write(io_L,107) 2,2,Yu(2), "# Y_c(Q)^DRbar"
  Write(io_L,107) 3,3,Yu(3), "# Y_t(Q)^DRbar"

  Write(io_L,106) "Block Yd Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yd(1), "# Y_d(Q)^DRbar"
  Write(io_L,107) 2,2,Yd(2), "# Y_s(Q)^DRbar"
  Write(io_L,107) 3,3,Yd(3), "# Y_b(Q)^DRbar"

  Write(io_L,106) "Block Ye Q=",Q,"# (SUSY scale)"

  ierr = 0
  If (GenerationMixing) Then
   !-------------------------------------------
   ! check if any off-diagonal term is non-zero
   !-------------------------------------------
   Do i1=1,3
    Do i2=1,3
     If ((i1.Ne.i2).And.(Abs(yl(i2,i1)).Ne.0._dp)) ierr = ierr + 1
    End Do
   End Do
  End If

  If (ierr.Ne.0) Then
   Do i1=1,3
    Do i2=1,3
     Write(io_L,105) i2,i1,Real(yl(i2,i1),dp),"# Y_(l,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
  Else 
   Write(io_L,107) 1,1,Real(Yl(1,1),dp), "# Y_e(Q)^DRbar"
   Write(io_L,107) 2,2,Real(Yl(2,2),dp), "# Y_mu(Q)^DRbar"
   Write(io_L,107) 3,3,Real(Yl(3,3),dp), "# Y_tau(Q)^DRbar"
  End If
  If (Maxval(Abs(Aimag(yl))).Gt.0._dp) Then
   Write(io_L,106) "Block IMYe Q=",Q,"# (SUSY scale)"
   If (GenerationMixing) Then
    Do i1=1,3
      Do i2=1,3
      Write(io_L,105) i2,i1,Aimag(yl(i2,i1)),"# Im(Y_(l,"//bu(i2)//bu(i1)//"))"
     End Do
    End Do
   Else 
    Write(io_L,107) 1,1,Aimag(Yl(1,1)), "# Im(Y_e)(Q)^DRbar"
    Write(io_L,107) 2,2,Aimag(Yl(2,2)), "# Im(Y_mu)(Q)^DRbar"
    Write(io_L,107) 3,3,Aimag(Yl(3,3)), "# Im(Y_tau)(Q)^DRbar"
   End If
  End If

104 Format(i4,2x,1P,e16.8,2x,a)
105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)
107 Format(2i3,3x,1P,e16.8,4x,a)

 end Subroutine LH_Write_GY

 Subroutine LH_write_decays(io_L, io, n_c, n_n, n_sl, n_sn, n_sd, n_su, n_S0 &
          & , n_P0, n_Spm, id_l, id_nu, id_d, id_u, id_W, id_Z, id_ph, id_gl &
          & , id_glu, id_grav, id_c, id_n, id_sle, id_snu, c_min, n_min      &
          & , id_sd, id_su, id_S0, id_P0, id_Sm, id_Sp, c_gl, c_lm           &
          & , c_lp, c_nu, c_d, c_u, c_Cm, c_Cp, c_c0, c_sle, c_slep, c_Snu   &
          & , c_Sd, c_Su, c_S0, c_P0, c_Sm, c_Sp, c_grav, c_phot             &
          & , gT_C, gP_C2, gP_C3, BR_C2, BR_C3                               &
          & , gT_N, gP_N2, gP_N3, BR_N2, BR_N3, gT_S0, gP_S0, BR_S0, gT_P0   &
          & , gP_P0, BR_P0, gT_Spm, gP_Spm, BR_Spm )

 implicit none
  Integer, Intent(in) :: io_L, io, n_c, n_n, n_sl, n_sn, n_sd, n_su, n_S0     &
   & , n_P0, n_Spm, id_l(3), id_nu(3), id_d(3), id_u(3), id_c(:), id_n(:)     &
   & , id_sle(:), id_snu(:), id_sd(:), id_su(:), id_S0(:), id_P0(:), id_Sm(:) &
   & , id_Sp(:), id_W, id_Z, id_glu, id_grav, id_gl, id_ph, c_min, n_min 
  Real(dp), Intent(in) :: gT_N(:), gP_N2(:,:), gP_N3(:,:), BR_N2(:,:)       &
   & , BR_N3(:,:), gT_S0(:), gP_S0(:,:), BR_S0(:,:), gT_P0(:), gP_P0(:,:)   &
   & , BR_P0(:,:), gT_Spm(:), gP_Spm(:,:), BR_Spm(:,:), gT_C(:), gP_C2(:,:) &
   & , gP_C3(:,:), BR_C2(:,:), BR_C3(:,:)
  Character(len=10), intent(in) :: c_snu(:), c_sle(:), c_slep(:)
  Character(len=6), intent(in) :: c_lm(3), c_lp(3), c_nu(3), c_sp(:), c_sm(:) &
   & , c_phot
  Character(len=1), intent(in) :: c_u(3), c_d(3)
  Character(len=2), intent(in) :: c_gl, c_grav
  Character(len=4), intent(in) :: c_S0(:), c_P0(:), c_su(6), c_sd(6)
  Character(len=5), intent(in) :: c_cp(:), c_cm(:), c_c0(:)

  Character(len=13) :: c_sfermion 

  Integer :: i_zaehl, id_d_2(200,2), id_d_3(400,3), i1, i2, i3, i4, id_f &
     & , id_fp, i_c2, i_eff
  Real(dp) :: BRmin100, BRtot, gT
  Logical :: l_3body
  Integer, Parameter :: n_max=500
  Real(dp), Dimension(n_max) :: Br, gP
  Character(len=30), Dimension(n_max) :: Fnames, Lnames
  Character(len=10) :: c_m

  Write(io,*) " "
  Write(io,*) " Anti particles are marked with a * in case of"
  Write(io,*) " (s)neutrinos and (s)quarks in the decay section."
  Write(io,*) "                    Decay widths (GeV) and branching ratios"
  Write(io,*) " "

  BRmin100 = 100._dp * BRmin
  !---------------------------------
  ! sleptons
  !---------------------------------
  Do i1=1,n_sl
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    id_fp = (i1+1)/2
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_nu(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_nu(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_nu(id_fp)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" W-"
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Do i3 = 1,3
     Fnames(i_zaehl) = "gravitino "//Trim(c_lm(i3))
     Lnames(i_zaehl) = Trim(c_grav)//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_grav
     id_d_2(i_zaehl,2) = id_l(i3)
     i_zaehl = i_zaehl+1
    End Do

    Call WriteDecays2(io, " slepton_"//Bu(i1) , Fnames &
                     &, gP_Sl(i1,:), 100*BR_Sl(i1,:), gT_Sl(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_lm(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_l(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = (i1+1)/2
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_nu(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_nu(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_nu(id_fp)
     i_zaehl = i_zaehl+1
    End Do
    i2=id_fp
    If (i2.Eq.1) Then
     c_sfermion = "e-sneutrino"
    Else If (i2.Eq.2) Then
     c_sfermion = "mu-sneutrino"
    Else If (i2.Eq.3) Then
     c_sfermion = "tau-sneutrino"
    End If
    Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W-"
    Lnames(i_zaehl) =  Trim(c_snu(i2))//" W-"
    id_d_2(i_zaehl,1) = id_snu(i2)
    id_d_2(i_zaehl,2) = -id_W
    i_zaehl = i_zaehl+1
    Do i3=1,n_spm
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sm(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_sm(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     i_zaehl = i_zaehl+1
    End Do
    If (i1.Le.2) Then
     c_sfermion = "selectron"
    Else If (i1.Le.4) Then
     c_sfermion = "smuon"
    Else 
     c_sfermion = "stau"
    End If

    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    i3 = (i1+1)/2
    Fnames(i_zaehl) = "gravitino "//Trim(c_lm(i3))
    Lnames(i_zaehl) = Trim(c_grav)//" "//Trim(c_lm(i3))
    id_d_2(i_zaehl,1) = id_grav
    id_d_2(i_zaehl,2) = id_l(i3)
    i_zaehl = i_zaehl+1

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2
    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Sl(i1,:), 100*BR_Sl(i1,:), gT_Sl(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_sle(i1),gT_sl(i1),Trim( c_sle(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sl(i1,i2).Gt.BrMin) Write(io_L,201) BR_Sl(i1,i2),2,id_d_2(i2,:), &
            &                  Trim( c_sle(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! sneutrinos
  !---------------------------------
  Do i1=1,n_sn
   If (GenerationMixing) Then
    i_zaehl = 1
    i3 = i1
    Do i2=1,n_n
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_nu(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_nu(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_nu(i3)
      i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_c
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_c(i2)
      id_d_2(i_zaehl,2) = id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_sl
     Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" W+"
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sl
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "slepton_"//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" Z"
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "sneutrino_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " Sneutrino_"//Bu(i1) , Fnames &
                     &, gP_Sn(i1,:), 100*BR_Sn(i1,:), gT_Sn(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = i1
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i1
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_lm(id_fp))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_lm(id_fp))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_l(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    If (i1.Eq.1) Then
     c_sfermion = "selectron_"
    Else If (i1.Eq.2) Then
     c_sfermion = "smuon_"
    Else If (i1.Eq.3) Then
     c_sfermion = "stau_"
    End If
    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sle(i2+(i1-1)*2))//" W+"
     id_d_2(i_zaehl,1) = id_sle(i2+(i1-1)*2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2+(i1-1)*2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sle(i2+(i1-1)*2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Eq.1) c_sfermion="e-sneutrino"
    If (i1.Eq.2) c_sfermion="mu-sneutrino"
    If (i1.Eq.3) c_sfermion="tau-sneutrino"

    Call WriteDecays2(io, Trim(c_sfermion) , Fnames &
                     &, gP_Sn(i1,:), 100*BR_Sn(i1,:), gT_Sn(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_snu(i1),gT_sn(i1),Trim(c_snu(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sn(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Sn(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_snu(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! d-squarks
  !---------------------------------
  Do i1=1,n_sd
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_c
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = -id_c(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,3
     Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_glu
     id_d_2(i_zaehl,2) = id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" W-"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " s-down_"//Bu(i1) , Fnames &
                     &, gP_Sd(i1,:), 100*BR_Sd(i1,:), gT_Sd(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i3
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_u(id_fp))
     Lnames(i_zaehl) =  Trim(c_cm(i2))//" "//Trim(c_u(id_fp))
     id_d_2(i_zaehl,1) = -id_c(i2)
     id_d_2(i_zaehl,2) = id_u(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))
    Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_d(i3))
    id_d_2(i_zaehl,1) = id_glu
    id_d_2(i_zaehl,2) = id_d(i3)
    i_zaehl = i_zaehl+1

    If (i1.Le.2) Then
     c_sfermion = "s-up_"
     id_fp = 0
    Else If (i1.Le.4) Then
     c_sfermion = "s-charm_"
     id_fp = 2
    Else 
     c_sfermion = "stop_"
     id_fp = 4
    End If

    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W-"
     Lnames(i_zaehl) =  Trim(c_su(i2 + id_fp))//" W-"
     id_d_2(i_zaehl,1) = id_su(i2 + id_fp)
     id_d_2(i_zaehl,2) = -id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sm(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2+id_fp))//" "//Trim(c_sm(i3))
      id_d_2(i_zaehl,1) = id_su(i2+id_fp)
      id_d_2(i_zaehl,2) = id_sm(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Le.2) Then
     c_sfermion = "s-down_"
    Else If (i1.Le.4) Then
     c_sfermion = "s-strange_"
    Else 
     c_sfermion = "sbottom_"
    End If


    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" Z"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2
    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Sd(i1,:), 100*BR_Sd(i1,:), gT_Sd(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_sd(i1),gT_sd(i1),Trim(c_sd(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Sd(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Sd(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_sd(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do


  End Do   

  !---------------------------------
  ! u-squarks
  !---------------------------------
  Do i1=1,n_su
   If (GenerationMixing) Then
    i_zaehl = 1
    Do i2=1,n_n
     Do i3=1,3
      Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = id_n(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_c
     Do i3=1,3
      Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_c(i2)
      id_d_2(i_zaehl,2) = id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i3=1,3
     Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_glu
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" W+"
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     Do i3=1,n_spm
      Fnames(i_zaehl) =  "s-down_"//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" Z"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" Z"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,i1-1
     Do i3=1,n_P0
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,i1-1
     Do i3=1,n_S0
      Fnames(i_zaehl) =  "s-up_"//Bu(i2)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    Call WriteDecays2(io, " s-up_"//Bu(i1) , Fnames &
                     &, gP_Su(i1,:), 100*BR_Su(i1,:), gT_Su(i1), BrMin100)

   Else ! GenerationMixing

    i_zaehl = 1
    i3 = (i1+1)/2
    Do i2=1,n_n
     Fnames(i_zaehl) =  "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do
    id_fp = i3
    Do i2=1,n_c
     Fnames(i_zaehl) =  "chargino_"//Bu(i2)//" "//Trim(c_d(id_fp))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_d(id_fp))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_d(id_fp)
     i_zaehl = i_zaehl+1
    End Do

    Fnames(i_zaehl) =  "gluino "//Trim(c_u(id_fp))
    Lnames(i_zaehl) =  Trim(c_gl)//" "//Trim(c_u(id_fp))
    id_d_2(i_zaehl,1) = id_glu
    id_d_2(i_zaehl,2) = id_u(id_fp)
    i_zaehl = i_zaehl+1

    If (i1.Le.2) Then
     c_sfermion = "s-down_"
     id_fp = 0
    Else If (i1.Le.4) Then
     c_sfermion = "s-strange_"
     id_fp = 2
    Else 
     c_sfermion = "sbottom_"
     id_fp = 4
    End If

    Do i2=1,2
     Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" W+"
     Lnames(i_zaehl) =  Trim(c_sd(i2 + id_fp))//" W+"
     id_d_2(i_zaehl,1) = id_sd(i2 + id_fp)
     id_d_2(i_zaehl,2) = id_W
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,2
     Do i3=1,n_spm
      Fnames(i_zaehl) =  Trim(c_sfermion)//Bu(i2)//" "//Trim(c_sp(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2+id_fp))//" "//Trim(c_sp(i3))
      id_d_2(i_zaehl,1) = id_su(i2+id_fp)
      id_d_2(i_zaehl,2) = id_sp(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

    If (i1.Le.2) Then
     c_sfermion = "s-up_"
    Else If (i1.Le.4) Then
     c_sfermion = "s-charm_"
    Else 
     c_sfermion = "stop_"
    End If


    If ((i1.Eq.2).Or.(i1.Eq.4).Or.(i1.Eq.6)) Then
     i2 = i1-1
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" Z"
     Lnames(i_zaehl) =  Trim(c_su(i2))//" Z"
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = id_Z
     i_zaehl = i_zaehl+1

     Do i3=1,n_P0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_p0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_p0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_p0(i3)
      i_zaehl = i_zaehl+1
     End Do

     Do i3=1,n_S0
      Fnames(i_zaehl) = Trim(c_sfermion)//Bu(1)//" "//Trim(c_S0(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_S0(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = id_S0(i3)
      i_zaehl = i_zaehl+1
     End Do
    End If ! i1 is even

    id_f = i1
    If (i1.Gt.2) id_f = id_f - 2
    If (i1.Gt.4) id_f = id_f - 2

    If (i1.Eq.5) Then ! 3-body decays of lighter stop
     id_d_3 = 0
     Fnames(55) = "neutralino_1 c"
     Lnames(55) = Trim(c_c0(1))//" c"
     id_d_2(55,1) = id_n(1)
     id_d_2(55,2) = id_u(2)
     Fnames(56) = "neutralino_2 c"
     Lnames(56) = Trim(c_c0(2))//" c"
     id_d_2(56,1) = id_n(2)
     id_d_2(56,2) = id_u(2)
     i_zaehl = 57
     Fnames(57) = "neutralino_1 W b"
     Lnames(57) = Trim(c_c0(1))//" W+ b"
     id_d_3(1,1) = id_n(1)
     id_d_3(1,2) = id_W
     id_d_3(1,3) = id_d(3)
     Fnames(58) = "e-sneutrino e+ b"
     Lnames(58) = Trim(c_snu(1))//" e+ b"
     id_d_3(2,1) = id_snu(1)
     id_d_3(2,2) = -id_l(1)
     id_d_3(2,3) = id_d(3)
     Fnames(59) = "mu-sneutrino mu+ b"
     Lnames(59) = Trim(c_snu(2))//" mu+ b"
     id_d_3(3,1) = id_snu(2)
     id_d_3(3,2) = -id_l(2)
     id_d_3(3,3) = id_d(3)
     Fnames(60) = "tau-sneutrino tau+ b"
     Lnames(60) = Trim(c_snu(3))//" tau+ b"
     id_d_3(4,1) = id_snu(3)
     id_d_3(4,2) = -id_l(3)
     id_d_3(4,3) = id_d(3)
     Fnames(61) = "selectron_1 nu_e b"
     Lnames(61) = Trim(c_slep(1))//" nu_e b"
     id_d_3(5,1) = -id_sle(1)
     id_d_3(5,2) = id_nu(1)
     id_d_3(5,3) = id_d(3)
     Fnames(62) = "selectron_2 nu_e b"
     Lnames(62) = Trim(c_slep(2))//" nu_e b"
     id_d_3(6,1) = -id_sle(2)
     id_d_3(6,2) = id_nu(1)
     id_d_3(6,3) = id_d(3)
     Fnames(63) = "smuon_1 nu_mu b"
     Lnames(63) = Trim(c_slep(3))//" nu_mu b"
     id_d_3(7,1) = -id_sle(3)
     id_d_3(7,2) = id_nu(2)
     id_d_3(7,3) = id_d(3)
     Fnames(64) = "smuon_2 nu_mu b"
     Lnames(64) = Trim(c_slep(4))//" nu_mu b"
     id_d_3(8,1) = -id_sle(4)
     id_d_3(8,2) = id_nu(2)
     id_d_3(8,3) = id_d(3)
     Fnames(65) = "stau_1 nu_tau b"
     Lnames(65) = Trim(c_slep(5))//" nu_tau b"
     id_d_3(9,1) = -id_sle(5)
     id_d_3(9,2) = id_nu(3)
     id_d_3(9,3) = id_d(3)
     Fnames(66) = "stau_2 nu_tau b"
     Lnames(66) = Trim(c_slep(6))//" nu_tau b"
     id_d_3(10,1) = -id_sle(6)
     id_d_3(10,2) = id_nu(3)
     id_d_3(10,3) = id_d(3)
    End If

    Call WriteDecays2(io, Trim(c_sfermion)//Bu(id_f) , Fnames &
                     &, gP_Su(i1,:), 100*BR_Su(i1,:), gT_Su(i1), BrMin100)
   End If ! GenerationMixing

   Write(io_L,200) id_su(i1),gT_su(i1),Trim(c_su(i1))
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl - 1
    If (BR_Su(i1,i2).Gt.BrMin) &
      & Write(io_L,201) BR_Su(i1,i2),2,id_d_2(i2,:), &
            &          Trim(c_su(i1))//" -> "//Trim(Lnames(i2))//")"
   End Do

   If ((i1.Eq.5).And.(Maxval(BR_Su(i1,57:66)).Gt.BRmin)) Then
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,10
      If  (BR_Su(i1,56+i2).Gt.BrMin) Write(io_L,202) BR_Su(i1,56+i2),3 &
             & ,id_d_3(i2,:), Trim(c_su(i1))//" -> "//Trim(Lnames(56+i2))//")"
    End Do
   End If

  End Do   

  !--------------
  ! charginos
  !--------------

  Do i1=c_min,n_c
   i_zaehl = 1
   If (GenerationMixing) Then
    Do i2=1,n_sl
     i3 = (i2+1)/2
     Fnames(i_zaehl) = "slepton_"//Bu(i2)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_slep(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = -id_sle(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     Do i3 = 1,3
      Fnames(i_zaehl) = "sneutrino_"//Bu(i2)//" "//Trim(c_lp(i3))
      Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_lp(i3))
      id_d_2(i_zaehl,1) = id_snu(i2)
      id_d_2(i_zaehl,2) = -id_l(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_su
     Do i3 = 1,3
      Fnames(i_zaehl) = "s-up_"//Bu(i2)//" "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = -id_d(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do
    Do i2=1,n_sd
     Do i3 = 1,3
      Fnames(i_zaehl) = "s-down_"//Bu(i2)//" "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = -id_sd(i2)
      id_d_2(i_zaehl,2) = id_u(i3)
      i_zaehl = i_zaehl+1
     End Do
    End Do

   Else ! GenerationMixing

    Do i2=1,n_sl
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="selectron_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="smuon_"
      id_f = i2 - 2
     Else 
      c_sfermion="stau_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_nu(i3))
     Lnames(i_zaehl) =  Trim(c_slep(i2))//" "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = -id_sle(i2)
     id_d_2(i_zaehl,2) = id_nu(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sn
     i3 = i2
     If (i2.Eq.1) c_sfermion = "e-sneutrino"
     If (i2.Eq.2) c_sfermion = "mu-sneutrino"
     If (i2.Eq.3) c_sfermion = "tau-sneutrino"
     Fnames(i_zaehl) = Trim(c_sfermion)//" "//Trim(c_lp(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_lp(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_l(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_su
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-up_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-charm_"
      id_f = i2 - 2
     Else 
      c_sfermion="stop_"
      id_f = i2 - 2
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     i_zaehl = i_zaehl+1
    End Do
    Do i2=1,n_sd
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-down_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-strange_"
      id_f = i2 - 2
     Else 
      c_sfermion="sbottom_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = -id_sd(i2)
     id_d_2(i_zaehl,2) = id_u(i3)
     i_zaehl = i_zaehl+1
    End Do

   End If ! GenerationMixing

   Do i2=1,n_n
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" W+"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" W+"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_W
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_spm
    Do i2=1,n_n
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_sp(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i2=1,i1-1
    Fnames(i_zaehl) = "chargino_"//Bu(i2)//" Z"
    Lnames(i_zaehl) =  Trim(c_cp(i2))//" Z"
    id_d_2(i_zaehl,1) = id_c(i2)
    id_d_2(i_zaehl,2) = id_Z
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_p0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_p0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i3=1,n_s0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   i_c2 = i_zaehl-1
   gP = 0._dp
   BR = 0._dp

   i_eff = i1 - (n_c-2)

   gP(1:i_c2) = gP_C2(i_eff,1:i_c2)
   BR(1:i_c2) = BR_C2(i_eff,1:i_c2)

   Write(io_L,200) id_c(i1),gT_C(i_eff),Trim(c_cp(i1))
   If (Sum(gp_C2(i_eff,:)).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,i_zaehl - 1
     If (BR_C2(i_eff,i2).Gt.BrMin) Write(io_L,201) BR_C2(i_eff,i2) &
            & ,2,id_d_2(i2,:), Trim(c_cp(i1))//" -> "//Trim(Lnames(i2))//")"
    End Do
   End If

   !-------------------------------------------------
   ! 3-body decays
   !-------------------------------------------------
   If (Maxval(BR_C3(i_eff,:)).Gt.BRmin) Then
    If (GenerationMixing) Then
     Do i2=1,n_n
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl-i_c2,1) = id_n(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
        id_d_3(i_zaehl-i_c2,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3-(n_c-2)
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_nu(i3)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_glu
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
        id_d_3(i_zaehl-i_c2,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
        id_d_3(i_zaehl-i_c2,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,5-n_c
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl-i_c2,1) = id_c(i2)
        id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
        id_d_3(i_zaehl-i_c2,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      If ((n_c-2).Lt.3) Then
       Fnames(i_zaehl) = "chargino_"//Bu(i2)//"nu_i nu_j"
       Lnames(i_zaehl) = Trim(c_cp(i2))//"nu_i nu_j"
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_nu(1)
       id_d_3(i_zaehl-i_c2,3) = id_nu(3)
       i_zaehl = i_zaehl+1    
      End If
     End Do
  
    Else ! GenerationMixing

     Do i2=1,n_n
      Do i3=1,3
       Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,5-n_c
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_nu(i3)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_u(i3))
      Lnames(i_zaehl) = Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
      id_d_3(i_zaehl-i_c2,1) = id_glu
      id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = id_u(i3)
      i_zaehl = i_zaehl+1    
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = Trim(c_cp(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i3))
       Lnames(i_zaehl) = Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i3)
       i_zaehl = i_zaehl+1    
      End Do
      Do i3=1,5-n_c
       Fnames(i_zaehl) = &
         & "chargino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl) = &
         & Trim(c_cp(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i3))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_l(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i3)
       i_zaehl = i_zaehl+1    
      End Do
      If ((n_c-2).Lt.3) Then
       Fnames(i_zaehl) = "chargino_"//Bu(i2)//"nu_i nu_i"
       Lnames(i_zaehl) = Trim(c_cp(i2))//"nu_i nu_i"
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_nu(1)
       id_d_3(i_zaehl-i_c2,3) = id_nu(1)
       i_zaehl = i_zaehl+1    
      End If
     End Do

    End If ! GenerationMixing

    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    Do i2=1,i_zaehl - i_c2 - 1
     If (BR_C3(i_eff,i2).Gt.BrMin) Write(io_L,202) BR_C3(i_eff,i2),3 &
           & ,id_d_3(i2,:) , Trim(c_cp(i1))//" -> "//Trim(Lnames(i2+i_c2))//")"
    End Do

    BR(i_c2+1:i_zaehl-1) = BR_C3(i_eff,1:i_zaehl-i_c2-1)
    gP(i_c2+1:i_zaehl-1) = gP_C3(i_eff,1:i_zaehl-i_c2-1)
   End If ! check for maxval(BR)

   Call WriteDecays2(io, " chargino_"//Bu(i1) , Fnames &
                       &, gP, 100*BR, gT_C(i_eff), BrMin100)

  End Do ! i1


  !--------------
  ! neutralinos
  !--------------
  Do i1=n_min,n_n
   i_zaehl = 1
   If (GenerationMixing) Then
    Do i2=1,n_sl
     Do i3 = 1,3
      Fnames(i_zaehl) = "slepton_"//Bu(i2)//" "//Trim(c_lp(i3))
      Fnames(i_zaehl+1) = "slepton_"//Bu(i2)//"^+ "//Trim(c_lm(i3))
      Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_lp(i3))
      Lnames(i_zaehl+1) =  Trim(c_slep(i2))//" "//Trim(c_lm(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = -id_l(i3)
      id_d_2(i_zaehl+1,1) = -id_sle(i2)
      id_d_2(i_zaehl+1,2) = id_l(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
    Do i2=1,n_sn
     i3 = i2
     Fnames(i_zaehl) = "sneutrino_"//Bu(i2)//" "//Trim(c_snu(i3))//"^*"
     Fnames(i_zaehl+1) = "sneutrino_"//Bu(i2)//"^* "//Trim(c_snu(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_nu(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_snu(i2))//"^* "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_nu(i3)
     id_d_2(i_zaehl+1,1) = -id_snu(i2)
     id_d_2(i_zaehl+1,2) = id_nu(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_su
     Do i3 = 1,3
      Fnames(i_zaehl) = "sup_"//Bu(i2)//" "//Trim(c_u(i3))//"^*"
      Fnames(i_zaehl+1) = "sup_"//Bu(i2)//"^* "//Trim(c_u(i3))
      Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
      Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
      id_d_2(i_zaehl,1) = id_su(i2)
      id_d_2(i_zaehl,2) = -id_u(i3)
      id_d_2(i_zaehl+1,1) = -id_su(i2)
      id_d_2(i_zaehl+1,2) = id_u(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
    Do i2=1,n_sd
     Do i3 = 1,3
      Fnames(i_zaehl) = "sdown_"//Bu(i2)//" "//Trim(c_d(i3))//"^*"
      Fnames(i_zaehl+1) = "sdown_"//Bu(i2)//"^* "//Trim(c_d(i3))
      Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
      Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_d(i3)
      id_d_2(i_zaehl+1,1) = -id_sd(i2)
      id_d_2(i_zaehl+1,2) = id_d(i3)
      i_zaehl = i_zaehl+2
     End Do
    End Do
    
   Else ! GenerationMixing

    Do i2=1,n_sl
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="selectron_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="smuon_"
      id_f = i2 - 2
     Else 
      c_sfermion="stau_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_lp(i3))
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^+ "//Trim(c_lm(i3))
     Lnames(i_zaehl) =  Trim(c_sle(i2))//" "//Trim(c_lp(i3))
     Lnames(i_zaehl+1) =  Trim(c_slep(i2))//" "//Trim(c_lm(i3))
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = -id_l(i3)
     id_d_2(i_zaehl+1,1) = -id_sle(i2)
     id_d_2(i_zaehl+1,2) = id_l(i3)
     i_zaehl = i_zaehl+2
    End Do

    Do i2=1,n_sn
     i3 = i2
     If (i2.Eq.1) c_sfermion = "e-sneutrino"
     If (i2.Eq.2) c_sfermion = "mu-sneutrino"
     If (i2.Eq.3) c_sfermion = "tau-sneutrino"
     Fnames(i_zaehl) = Trim(c_sfermion)//" "//Trim(c_snu(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//"^* "//Trim(c_snu(i3))
     Lnames(i_zaehl) =  Trim(c_snu(i2))//" "//Trim(c_nu(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_snu(i2))//"^* "//Trim(c_nu(i3))
     id_d_2(i_zaehl,1) = id_snu(i2)
     id_d_2(i_zaehl,2) = -id_nu(i3)
     id_d_2(i_zaehl+1,1) = -id_snu(i2)
     id_d_2(i_zaehl+1,2) = id_nu(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_su
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-up_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-charm_"
      id_f = i2 - 2
     Else 
      c_sfermion="stop_"
      id_f = i2 - 2
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_u(i3)
     id_d_2(i_zaehl+1,1) = -id_su(i2)
     id_d_2(i_zaehl+1,2) = id_u(i3)
     i_zaehl = i_zaehl+2
    End Do
    
    Do i2=1,n_sd
     i3 = (i2+1)/2
     If (i2.Le.2) Then
      c_sfermion="s-down_"
      id_f = i2
     Else If (i2.Le.4) Then
      c_sfermion="s-strange_"
      id_f = i2 - 2
     Else 
      c_sfermion="sbottom_"
      id_f = i2 - 4
     End If
     Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))//"^*"
     Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     id_d_2(i_zaehl+1,1) = -id_sd(i2)
     id_d_2(i_zaehl+1,2) = id_d(i3)
     i_zaehl = i_zaehl+2
    End Do
    
   End If ! GenerationMixing

   Do i2=1,n_c
    Fnames(i_zaehl) = "chargino_"//Bu(i2)//"^+ W-"
    Fnames(i_zaehl+1) = "chargino_"//Bu(i2)//"^- W+"
    Lnames(i_zaehl) =  Trim(c_cp(i2))//" W-"
    Lnames(i_zaehl+1) =  Trim(c_cm(i2))//" W+"
    id_d_2(i_zaehl,1) = id_c(i2)
    id_d_2(i_zaehl,2) = -id_W
    id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl+2 
   End Do

   Do i3=1,n_spm
    Do i2=1,n_c
     Fnames(i_zaehl) = "chargino_"//Bu(i2)//"^+ "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = "chargino_"//Bu(i2)//"^- "//Trim(c_sm(i3))
     Lnames(i_zaehl) =  Trim(c_cp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) =  Trim(c_cm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_c(i2)
     id_d_2(i_zaehl,2) = -id_sp(i3)
     id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl+2   
    End Do
   End Do


   Do i2=1,i1-1
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" Z"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" Z"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_Z
    i_zaehl = i_zaehl+1    
   End Do

   Do i3=1,n_p0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_p0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Do i3=1,n_s0
    Do i2=1,i1-1
     Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) =  Trim(c_c0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl+1    
    End Do
   End Do

   Fnames(i_zaehl) = "Gravitino photon"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_ph
   i_zaehl = i_zaehl + 1
   Fnames(i_zaehl) = "Gravitino Z"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_Z
   i_zaehl = i_zaehl + 1
   Fnames(i_zaehl) = "Gravitino h0"
   Lnames(i_zaehl) = Fnames(i_zaehl)
   id_d_2(i_zaehl,1) = id_grav
   id_d_2(i_zaehl,2) = id_s0(1)
   i_zaehl = i_zaehl + 1

   If (HighScaleModel.Eq."NMSSM") Then
    i_zaehl = 150 
   Else If (HighScaleModel.Eq."RPexplicit") Then
    i_zaehl = 186
   Else If (HighScaleModel.Eq."NURRP1") Then
    i_zaehl = 190
   Else
    i_zaehl = 150 
   End If

   Do i2=1,i1-1
    Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" photon"
    Lnames(i_zaehl) =  Trim(c_c0(i2))//" photon"
    id_d_2(i_zaehl,1) = id_n(i2)
    id_d_2(i_zaehl,2) = id_ph
    i_zaehl = i_zaehl+1    
   End Do

   i_c2 = i_zaehl-1
   gP = 0._dp
   BR = 0._dp
   
   i_eff = i1 + 1 - n_min

   gT = gT_N(i_eff)
   gP(1:i_c2) = gP_N2(i_eff,1:i_c2)
   BR(1:i_c2) = BR_N2(i_eff,1:i_c2)

   Write(io_L,200) id_n(i1),gT,Trim(c_c0(i1))
   If (Sum(gP).Gt.0._dp) Then ! 2-body decays
    Write(io_L,100) "#    BR                NDA      ID1      ID2"
    Do i2=1,i_zaehl - 1
     If (BR(i2).Gt.BrMin) Write(io_L,201) BR(i2),2,id_d_2(i2,:), &
            &                   Trim(c_c0(i1))//" -> "//Trim(Lnames(i2))//")"
    End Do
   End If

   !-------------------------------------------------
   ! 3-body decays
   !-------------------------------------------------
   l_3body = (Maxval(BR_N3(i_eff,:)).Gt.BRmin)
   i_zaehl = i_zaehl + (n_n-i1+1)
   i_c2 = i_zaehl - 1

   If (l_3body) Then
    If (GenerationMixing) Then
     Do i2=1,n_c
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_c(i2)
        id_d_3(i_zaehl,2) = id_d(i3)
        id_d_3(i_zaehl,3) = -id_u(i4)
        id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
        i_zaehl = i_zaehl+2
       End Do
      End Do
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_l(i3)
       id_d_3(i_zaehl,3) = -id_nu(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do   

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_u(i3)
       id_d_3(i_zaehl,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_d(i3)
       id_d_3(i_zaehl,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do

     Do i2=1,i1-1+(n_n-4)
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_u(i3)
        id_d_3(i_zaehl,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_d(i3)
        id_d_3(i_zaehl,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
      Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" nu_i nu_j"
      Lnames(i_zaehl) = Trim(c_c0(i2))//" nu_i nu_j"
      id_d_3(i_zaehl,1) = id_n(i2)
      id_d_3(i_zaehl,2) = - id_nu(1)
      id_d_3(i_zaehl,3) = id_nu(3)
      i_zaehl = i_zaehl+1    
      Do i3=1,3
       Do i4=1,3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_l(i3)
        id_d_3(i_zaehl,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do
      End Do
     
     End Do
  
    Else ! GenerationMixing

     Do i2=1,n_c
      Do i3=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_d(i3)
       id_d_3(i_zaehl,3) = -id_u(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
      Do i3=1,5-n_c
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_nu(i3))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_nu(i3))
       id_d_3(i_zaehl,1) = id_c(i2)
       id_d_3(i_zaehl,2) = id_l(i3)
       id_d_3(i_zaehl,3) = -id_nu(i3)
       id_d_3(i_zaehl+1,:) = -id_d_3(i_zaehl,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do   

      Do i3=1,3
       i4 = i3
       Fnames(i_zaehl) =  "gluino "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_u(i3)
       id_d_3(i_zaehl,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do

      Do i3=1,3
       i4=i3
       Fnames(i_zaehl) =  "gluino "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       Lnames(i_zaehl) = &
         & Trim(c_gl)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))//"^*"
       id_d_3(i_zaehl,1) = id_glu
       id_d_3(i_zaehl,2) = - id_d(i3)
       id_d_3(i_zaehl,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do

      Do i2=1,i1-1
       Do i3=1,3
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_u(i3)
        id_d_3(i_zaehl,3) = id_u(i4)
        i_zaehl = i_zaehl+1    
       End Do
       Do i3=1,3
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_d(i3)
        id_d_3(i_zaehl,3) = id_d(i4)
        i_zaehl = i_zaehl+1    
       End Do
       If ((n_n-4).Lt.3) Then
        Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" nu_i nu_i"
        Lnames(i_zaehl) = Trim(c_c0(i2))//" nu_i nu_i"
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_nu(1)
        id_d_3(i_zaehl,3) = id_nu(3)
        i_zaehl = i_zaehl+1    
       End If
       Do i3=1,7-n_n
        i4 = i3
        Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
        id_d_3(i_zaehl,1) = id_n(i2)
        id_d_3(i_zaehl,2) = - id_l(i3)
        id_d_3(i_zaehl,3) = id_l(i4)
        i_zaehl = i_zaehl+1    
       End Do    
      End Do
 
     If (     (HighScaleModel.Eq."RPexplicit") &
        & .or.(HighScaleModel(1:5).Eq."NURRP") ) Then
      Do i2=1,n_c
       Do i3=1,Min(i2,3)
        If (i2.Eq.i3) Then
         Fnames(i_zaehl) = &
            & "neutrino_e"//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Lnames(i_zaehl) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         id_d_3(i_zaehl,1) = id_n(1)
         id_d_3(i_zaehl,2) = - id_l(i2)
         id_d_3(i_zaehl,3) = id_l(i3)
         i_zaehl = i_zaehl+1
        Else
         Fnames(i_zaehl) = "neutrino_e"//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Fnames(i_zaehl+1) = &
            & "neutrino_e"//" "//Trim(c_cp(i3))//" "//Trim(c_cm(i2))
         Lnames(i_zaehl) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i2))//" "//Trim(c_cm(i3))
         Lnames(i_zaehl+1) = &
            & Trim(c_c0(1))//" "//Trim(c_cp(i3))//" "//Trim(c_cm(i2))
         id_d_3(i_zaehl,1) = id_n(1)
         id_d_3(i_zaehl,2) = - id_c(i2)
         id_d_3(i_zaehl,3) = id_c(i3)
         id_d_3(i_zaehl+1,1) = id_n(1)
         id_d_3(i_zaehl+1,2:3) = -id_d_3(i_zaehl,2:3) 
         i_zaehl = i_zaehl+2
        End If
       End Do
      End Do
      Do i2=4,i1-1
       Do i3=1,3
        Do i4=1,i3
         If (i3.Eq.i4) Then
          Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          id_d_3(i_zaehl,1) = id_n(i2)
          id_d_3(i_zaehl,2) = - id_l(i3)
          id_d_3(i_zaehl,3) = id_l(i4)
          i_zaehl = i_zaehl+1
         Else
          Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Fnames(i_zaehl+1) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_lm(i3))//" "//Trim(c_lp(i4))
          Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_lp(i3))//" "//Trim(c_lm(i4))
          Lnames(i_zaehl+1) = &
            & Trim(c_c0(i2))//" "//Trim(c_lm(i3))//" "//Trim(c_lp(i4))
          id_d_3(i_zaehl,1) = id_n(i2)
          id_d_3(i_zaehl,2) = - id_l(i3)
          id_d_3(i_zaehl,3) = id_l(i4)
          id_d_3(i_zaehl+1,1) = id_n(i2)
          id_d_3(i_zaehl+1,2:3) = -id_d_3(i_zaehl,2:3) 
          i_zaehl = i_zaehl+2
         End If
        End Do
       End Do
      End Do
      Do i2=4,i1-1
       Fnames(i_zaehl) = &
           & "neutralino_"//Bu(i2-3)//" "//Trim(c_nu(1))//" "//Trim(c_nu(1))
       Lnames(i_zaehl) = &
            & Trim(c_c0(i2))//" "//Trim(c_c0(1))//" "//Trim(c_c0(1))
       id_d_3(i_zaehl,1) = id_n(i2)
       id_d_3(i_zaehl,2) = id_n(1)
       id_d_3(i_zaehl,3) = id_n(1)
       i_zaehl = i_zaehl+1
      End Do
      Fnames(i_zaehl) = Trim(c_nu(1))//" "//Trim(c_nu(1))//" "//Trim(c_nu(1))
      Lnames(i_zaehl) =  Trim(c_c0(1))//" "//Trim(c_c0(1))//" "//Trim(c_c0(1))
      id_d_3(i_zaehl,1) = id_n(1)
      id_d_3(i_zaehl,2) = id_n(1)
      id_d_3(i_zaehl,3) = id_n(1)
      i_zaehl = i_zaehl+1
     End If
    End If ! GenerationMixing
   
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"

    BR(i_c2+1:i_zaehl-1) = BR_N3(i_eff,1:i_zaehl-i_c2-1)
    gP(i_c2+1:i_zaehl-1) = gP_N3(i_eff,1:i_zaehl-i_c2-1)

    Do i2=i_c2 + 1,i_zaehl - 1
     If (BR(i2).Gt.BrMin) Write(io_L,202) BR(i2),3,id_d_3(i2,:) &
             &                , Trim(c_c0(i1))//" -> "//Trim(Lnames(i2))//")"
    End Do

   End If

  Call WriteDecays2(io, " neutralino_"//Bu(i1), Fnames, gP, 100*BR, gT, BrMin100)
  End Do

  !--------------
  ! gluino
  !--------------
  i_zaehl = 1
  If (GenerationMixing) Then
   Do i2=1,n_su
    Do i3 = 1,3
     Fnames(i_zaehl) = "sup_"//Bu(i2)//" "//Trim(c_u(i3))//"^*"
     Fnames(i_zaehl+1) = "sup_"//Bu(i2)//"^* "//Trim(c_u(i3))
     Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
     id_d_2(i_zaehl,1) = id_su(i2)
     id_d_2(i_zaehl,2) = -id_u(i3)
     id_d_2(i_zaehl+1,1) = -id_su(i2)
     id_d_2(i_zaehl+1,2) = id_u(i3)
     i_zaehl = i_zaehl+2
    End Do
   End Do
    
   Do i2=1,n_sd
    Do i3 = 1,3
     Fnames(i_zaehl) = "sdown_"//Bu(i2)//" "//Trim(c_d(i3))//"^*"
     Fnames(i_zaehl+1) = "sdown_"//Bu(i2)//"^* "//Trim(c_d(i3))
     Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
     id_d_2(i_zaehl,1) = id_sd(i2)
     id_d_2(i_zaehl,2) = -id_d(i3)
     id_d_2(i_zaehl+1,1) = -id_sd(i2)
     id_d_2(i_zaehl+1,2) = id_d(i3)
     i_zaehl = i_zaehl+2
    End Do
   End Do
 
    
  Else ! GenerationMixing

   Do i2=1,n_su
    i3 = (i2+1)/2
    If (i2.Le.2) Then
     c_sfermion="s-up_"
     id_f = i2
    Else If (i2.Le.4) Then
     c_sfermion="s-charm_"
     id_f = i2 - 2
    Else 
     c_sfermion="stop_"
     id_f = i2 - 2
    End If
    Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(i3))//"^*"
    Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(i3))
    Lnames(i_zaehl) =  Trim(c_su(i2))//" "//Trim(c_u(i3))//"^*"
    Lnames(i_zaehl+1) =  Trim(c_su(i2))//" "//Trim(c_u(i3))
    id_d_2(i_zaehl,1) = id_su(i2)
    id_d_2(i_zaehl,2) = -id_u(i3)
    id_d_2(i_zaehl+1,1) = -id_su(i2)
    id_d_2(i_zaehl+1,2) = id_u(i3)
    i_zaehl = i_zaehl+2
   End Do
    
   Do i2=1,n_sd
    i3 = (i2+1)/2
    If (i2.Le.2) Then
     c_sfermion="s-down_"
     id_f = i2
    Else If (i2.Le.4) Then
     c_sfermion="s-strange_"
     id_f = i2 - 2
    Else 
     c_sfermion="sbottom_"
     id_f = i2 - 4
    End If
    Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_d(i3))//"^*"
    Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_d(i3))
    Lnames(i_zaehl) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))//"^*"
    Lnames(i_zaehl+1) =  Trim(c_sd(i2))//" "//Trim(c_d(i3))
    id_d_2(i_zaehl,1) = id_sd(i2)
    id_d_2(i_zaehl,2) = -id_d(i3)
    id_d_2(i_zaehl+1,1) = -id_sd(i2)
    id_d_2(i_zaehl+1,2) = id_d(i3)
    i_zaehl = i_zaehl+2
   End Do
    
   id_f = 1
   c_sfermion="stop_"
   Fnames(i_zaehl) = Trim(c_sfermion)//Bu(id_f)//" "//Trim(c_u(2))//"^*"
   Fnames(i_zaehl+1) = Trim(c_sfermion)//Bu(id_f)//"^* "//Trim(c_u(2))
   Lnames(i_zaehl) =  Trim(c_su(5))//" "//Trim(c_u(2))//"^*"
   Lnames(i_zaehl+1) =  Trim(c_su(5))//" "//Trim(c_u(2))
   id_d_2(i_zaehl,1) = id_su(5)
   id_d_2(i_zaehl,2) = -id_u(2)
   id_d_2(i_zaehl+1,1) = -id_su(5)
   id_d_2(i_zaehl+1,2) = id_u(2)
   i_zaehl = i_zaehl+2

  End If ! GenerationMixing

  Do i2=1,n_n
   Fnames(i_zaehl) = "neutralino_"//Bu(i2)//" gluon"
   Lnames(i_zaehl) =  Trim(c_c0(i2))//" gluon"
   id_d_2(i_zaehl,1) = id_n(i2)
   id_d_2(i_zaehl,2) = id_gl
   i_zaehl = i_zaehl+1    
  End Do

  i_c2 = i_zaehl-1
  gP = 0._dp
  BR = 0._dp
  gP(1:i_c2) = gP_Glu2(1:i_c2)
  BR(1:i_c2) = BR_Glu2(1:i_c2)

  Write(io_L,200) id_glu,gT_Glu,Trim(c_gl)
  If (Sum(gp_Glu2).Gt.0._dp) Then ! 2-body decays
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1 
    If (BR_Glu2(i2).Gt.BrMin) Write(io_L,201) BR_Glu2(i2),2,id_d_2(i2,:), &
            &                   Trim(c_gl)//" -> "//Trim(Lnames(i2))//")"
   End Do
  End If

  !-------------------------------------------------
  ! 3-body decays
  !-------------------------------------------------
  If (Maxval(BR_Glu3).Gt.BRmin) Then
   If (GenerationMixing) Then
    Do i2=1,n_n
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
       id_d_3(i_zaehl-i_c2,3) = id_u(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
       Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
       id_d_3(i_zaehl-i_c2,1) = id_n(i2)
       id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = id_d(i4)
       i_zaehl = i_zaehl+1    
      End Do
     End Do
    End Do

    Do i2=1,n_c
     Do i3=1,3
      Do i4=1,3
       Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
       id_d_3(i_zaehl-i_c2,1) = id_c(i2)
       id_d_3(i_zaehl-i_c2,2) = id_d(i3)
       id_d_3(i_zaehl-i_c2,3) = -id_u(i4)
       id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
       i_zaehl = i_zaehl+2
      End Do
     End Do
    End Do 

    Do i2=1,6
     Do i3=1,3
      Fnames(i_zaehl) = "sup_"//Bu(i2)//"^* W+ "//Trim(c_d(i3))
      Fnames(i_zaehl+1) = "sup_"//Bu(i2)//" W- "//Trim(c_d(i3))//"^*"
      Lnames(i_zaehl) = Trim(c_su(i2))//"^* W+ "//Trim(c_d(i3))
      Lnames(i_zaehl+1) = Trim(c_su(i2))//" W- "//Trim(c_d(i3))//"^*"
      id_d_3(i_zaehl-i_c2,1) = -id_su(i2)
      id_d_3(i_zaehl-i_c2,2) = id_W
      id_d_3(i_zaehl-i_c2,3) = id_d(i3)
      id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
      i_zaehl = i_zaehl+2
     End Do
    End Do  

   Else ! GenerationMixing
    Do i2=1,n_n
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
        & "neutralino_"//Bu(i2)//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_u(i3))//" "//Trim(c_u(i4))
      id_d_3(i_zaehl-i_c2,1) = id_n(i2)
      id_d_3(i_zaehl-i_c2,2) = - id_u(i3)
      id_d_3(i_zaehl-i_c2,3) = id_u(i4)
      i_zaehl = i_zaehl+1    
     End Do
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
         & "neutralino_"//Bu(i2)//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
      Lnames(i_zaehl) = &
         & Trim(c_c0(i2))//" "//Trim(c_d(i3))//" "//Trim(c_d(i4))
      id_d_3(i_zaehl-i_c2,1) = id_n(i2)
      id_d_3(i_zaehl-i_c2,2) = - id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = id_d(i4)
      i_zaehl = i_zaehl+1    
     End Do
    End Do

    Do i2=1,n_c
     Do i3=1,3
      i4=i3
      Fnames(i_zaehl) = &
        & "chargino_"//Bu(i2)//"^+ "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Fnames(i_zaehl+1) = &
        & "chargino_"//Bu(i2)//"^- "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl) = &
        & Trim(c_cp(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      Lnames(i_zaehl+1) = &
        & Trim(c_cm(i2))//" "//Trim(c_d(i3))//" "//Trim(c_u(i4))
      id_d_3(i_zaehl-i_c2,1) = id_c(i2)
      id_d_3(i_zaehl-i_c2,2) = id_d(i3)
      id_d_3(i_zaehl-i_c2,3) = -id_u(i4)
      id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
      i_zaehl = i_zaehl+2
     End Do
    End Do 

    Do i2=5,6
     i3=3
     Fnames(i_zaehl) = "stop_"//Bu(i2-4)//"^* W+ "//Trim(c_d(i3))
     Fnames(i_zaehl+1) = "stop_"//Bu(i2-4)//" W- "//Trim(c_d(i3))//"^*"
     Lnames(i_zaehl) = Trim(c_su(i2))//"^* W+ "//Trim(c_d(i3))
     Lnames(i_zaehl+1) = Trim(c_su(i2))//" W- "//Trim(c_d(i3))//"^*"
     id_d_3(i_zaehl-i_c2,1) = -id_su(i2)
     id_d_3(i_zaehl-i_c2,2) = id_W
     id_d_3(i_zaehl-i_c2,3) = id_d(i3)
     id_d_3(i_zaehl-i_c2+1,:) = -id_d_3(i_zaehl-i_c2,:)
     i_zaehl = i_zaehl+2
    End Do  
   End If ! GenerationMixing

   Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
   Do i2=1,i_zaehl - i_c2 - 2
     If (BR_Glu3(i2).Gt.BrMin) Write(io_L,202) BR_glu3(i2),3,id_d_3(i2,:) &
             &                , Trim(c_gl)//" -> "//Trim(Lnames(i2+i_c2))//")"
   End Do

   BR(i_c2+1:i_zaehl-2) = BR_Glu3(1:i_zaehl-i_c2-2)
   gP(i_c2+1:i_zaehl-2) = gP_glu3(1:i_zaehl-i_c2-2)

  End If

  Call WriteDecays2(io, " gluino_" , Fnames, gP, 100*BR, gT_Glu, BrMin100)


  !-----------------
  ! neutral scalars
  !-----------------
  Do i1=1,n_s0
   c_m = c_s0(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    Do i2 = 1, 5 - n_c
     Do i3 = i2, 5 - n_c
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Fnames(i_zaehl+1) = c_lp(i2)//" "//c_lm(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl+1) = Trim(c_lm(i2))//" "//Trim(c_lP(i3))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,3
     Do i3 = i2,3
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i2)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+1) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Fnames(i_zaehl+10) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+1) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+10) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i3)
       id_d_2(i_zaehl+10,:) = - id_d_2(i_zaehl+9,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 9
    Do i2 = 1, n_sl
     Do i3 = i2, n_sl
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'slepton^+_'//Bu(i2)//' slepton^-_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       Lnames(i_zaehl+1) = Trim(c_sle(i2))//" "//Trim(c_slep(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1, n_sn
     Do i3 = i2, n_sn
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'sneutrino_'//Bu(i2)//' sneutrino^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       id_d_2(i_zaehl,1) = -id_snu(i2)
       id_d_2(i_zaehl,2) = id_snu(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'sneutrino_'//Bu(i2)//' sneutrino^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'sneutrino^*_'//Bu(i2)//' sneutrino_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i3))
       id_d_2(i_zaehl,1) = -id_snu(i2)
       id_d_2(i_zaehl,2) = id_snu(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,6
     Do i3 = i2,6
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'd-squark^*_'//Bu(i2)//' d-squark_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       Lnames(i_zaehl+1) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Fnames(i_zaehl+37) = &
                    & 'u-squark^*_'//Bu(i2)//' u-squark_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       Lnames(i_zaehl+37) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       id_d_2(i_zaehl+37,:) = - id_d_2(i_zaehl+36,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 36

   Else ! .not.GenerationMixing
    ! leptons
    Do i2=1,5 - n_c
     Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = -id_l(i2)
     id_d_2(i_zaehl,2) = id_l(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = -id_d(i2)
     id_d_2(i_zaehl,2) = id_d(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! u-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_u(i2))
     id_d_2(i_zaehl,1) = -id_u(i2)
     id_d_2(i_zaehl,2) = id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! sleptons
    Do i2=1,5 - n_c
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        Lnames(i_zaehl+1) = Trim(c_sle(i3))//" "//Trim(c_slep(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! sneutrinos 
    Do i2=1,n_sn
     Lnames(i_zaehl) = Trim(c_snu(i2))//" "//Trim(c_snu(i2))
     id_d_2(i_zaehl,1) = -id_snu(i2)
     id_d_2(i_zaehl,2) = id_snu(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        Lnames(i_zaehl+1) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! u-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        Lnames(i_zaehl+1) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2

    If ((n_c-2).Eq.1) Then
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     i2 = 2
    Else If ((n_c-2).Eq.3) Then
     i2 = 0
    Else
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     Fnames(3) = 'taus '
     i2 = 3
    End If
    Fnames(i2+1) = 'd-quark '
    Fnames(i2+2) = 's-quark '
    Fnames(i2+3) = 'b-quark '
    Fnames(i2+4) = 'u-quark '
    Fnames(i2+5) = 'c-quark '
    Fnames(i2+6) = 't-quark '

    If ((n_c-2).Eq.1) Then
     Fnames(9) = 'Selectron 1 1 '
     Fnames(10) = 'Selectron 1 2 '
     Fnames(11) = 'Selectron 2 1 '
     Fnames(12) = 'Selectron 2 2 '
     Fnames(13) = 'Smuon 1 1 '
     Fnames(14) = 'Smuon 1 2 '
     Fnames(15) = 'Smuon 2 1 '
     Fnames(16) = 'Smuon 2 2 '
     Fnames(17) = 'e-Sneutrino'
     Fnames(18) = 'mu-Sneutrino'
     i2 = 18
    Else If ((n_c-2).Eq.3) Then
     i2 = 6
    Else
     Fnames(10) = 'Selectron 1 1 '
     Fnames(11) = 'Selectron 1 2 '
     Fnames(12) = 'Selectron 2 1 '
     Fnames(13) = 'Selectron 2 2 '
     Fnames(14) = 'Smuon 1 1 '
     Fnames(15) = 'Smuon 1 2 '
     Fnames(16) = 'Smuon 2 1 '
     Fnames(17) = 'Smuon 2 2 '
     Fnames(18) = 'Stau 1 1 '
     Fnames(19) = 'Stau 1 2 '
     Fnames(20) = 'Stau 2 1 '
     Fnames(21) = 'Stau 2 2 '
     Fnames(22) = 'e-Sneutrino'
     Fnames(23) = 'mu-Sneutrino'
     Fnames(24) = 'tau-Sneutrino'
     i2 = 24
    End If
    Fnames(i2+1) = 'Sdown 1 1 '
    Fnames(i2+2) = 'Sdown 1 2 '
    Fnames(i2+3) = 'Sdown 2 1 '
    Fnames(i2+4) = 'Sdown 2 2 '
    Fnames(i2+5) = 'S-strange 1 1 '
    Fnames(i2+6) = 'S-strange 1 2 '
    Fnames(i2+7) = 'S-strange 2 1 '
    Fnames(i2+8) = 'S-strange 2 2 '
    Fnames(i2+9) = 'Sbottom 1 1 '
    Fnames(i2+10) = 'Sbottom 1 2 '
    Fnames(i1+11) = 'Sbottom 2 1 '
    Fnames(i1+12) = 'Sbottom 2 2 '
    Fnames(i1+13) = 'Sup 1 1 '
    Fnames(i1+14) = 'Sup 1 2 '
    Fnames(i1+15) = 'Sup 2 1 '
    Fnames(i1+16) = 'Sup 2 2 '
    Fnames(i1+17) = 'S-charm 1 1 '
    Fnames(i1+18) = 'S-charm 1 2 '
    Fnames(i1+19) = 'S-charm 2 1 '
    Fnames(i1+20) = 'S-charm 2 2 '
    Fnames(i1+21) = 'Stop 1 1 '
    Fnames(i1+22) = 'Stop 1 2 '
    Fnames(i1+23) = 'Stop 2 1 '
    Fnames(i1+24) = 'Stop 2 2 '

   End If

   ! neutralinos 
   Do i2=1,n_n
    Do i3=i2,n_n
     Fnames(i_zaehl) = 'neutralino_'//Bu(i2)//' neutralino_'//Bu(i3)
     Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   ! charginos
   Do i2=1, n_c
    Do i3=i2, n_c
     If (i2.Eq.i3) Then
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i2))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      i_zaehl = i_zaehl + 1
     Else
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Fnames(i_zaehl+1) = 'chargino^-_'//Bu(i2)//' chargino^+_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i3))
      Lnames(i_zaehl+1) = Trim(c_cp(i2))//" "//Trim(c_cm(i3))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
      i_zaehl = i_zaehl + 2
     End If
    End Do
   End Do

   Fnames(i_zaehl) = 'Z0 Z0'
   Lnames(i_zaehl) = 'Z0 Z0'
   id_d_2(i_zaehl,1) = id_Z
   id_d_2(i_zaehl,2) = id_Z
   i_zaehl = i_zaehl + 1

   Fnames(i_zaehl) = 'W+ W-'
   Lnames(i_zaehl) = 'W- W+'
   id_d_2(i_zaehl,1) = id_W
   id_d_2(i_zaehl,2) = - id_W
   i_zaehl = i_zaehl + 1

   Do i2=1,n_p0
    Fnames(i_zaehl) = 'Z0 '//Trim(c_p0(i2))
    Lnames(i_zaehl) = 'Z0 '//Trim(c_p0(i2))
    id_d_2(i_zaehl,1) = id_Z
    id_d_2(i_zaehl,2) = id_P0(i2)
    i_zaehl = i_zaehl + 1
   End Do
   Do i2=1,n_p0
    Do i3=i2,n_p0
     Fnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_p0(i3))
     Lnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_p0(i3))
     id_d_2(i_zaehl,1) = id_P0(i2)
     id_d_2(i_zaehl,2) = id_P0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i2=1,i1-1
    Do i3=i2,i1-1
     Fnames(i_zaehl) = Trim(c_s0(i2))//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) = Trim(c_s0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_S0(i2)
     id_d_2(i_zaehl,2) = id_S0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Fnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    Lnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Lnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_W
    id_d_2(i_zaehl,2) = id_sm(i2)
    id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl + 2
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    Lnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_sp(i2)
    id_d_2(i_zaehl,2) = id_sm(i2)
    i_zaehl = i_zaehl + 1
    Do i3=i2+1,n_spm
     Fnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_sp(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Fnames(i_zaehl) = "g g"
   Fnames(i_zaehl + 1) = Trim(c_phot)//" "//Trim(c_phot)
   Lnames(i_zaehl) = Fnames(i_zaehl)
   Lnames(i_zaehl+1) = Fnames(i_zaehl+1)
   id_d_2(i_zaehl,:) = id_gl
   id_d_2(i_zaehl+1,:) = id_ph
   i_zaehl = i_zaehl + 2

   gT = gT_S0(i1)
   BR(1:250) = BR_S0(i1,:)
   gP(1:250) = gP_S0(i1,:)

   Write(io_L,200) id_S0(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    
   
   Fnames(i_zaehl) = "W^+ W^-*" 
   Fnames(i_zaehl+1) = "W^- W^+*"
   Fnames(i_zaehl+2) = "Z Z^*" 

   Call WriteDecays2(io,c_m , Fnames, gP, 100*BR, gT, BrMin100)
   If ((BR(i_zaehl).Gt.BrMin).Or.(BR(i_zaehl+2).Gt.BrMin)) Then
    Write(io_L,100) "# writing decays into V V* as 3-body decays"
    Write(io_L,100) "#    BR                NDA      ID1      ID2       ID3"
    If (BR(i_zaehl).Gt.BrMin) Then
     Do i2=1,3
      Write(io_L,202) BrWln(i2) * BR(i_zaehl),3,id_w,id_l(i2),-id_nu(i2) &
      & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_lm(i2))//" "//Trim(c_nu(i2))//")"
      Write(io_L,202) BrWln(i2) * BR(i_zaehl),3,-id_w,-id_l(i2),id_nu(i2) &
      & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_lp(i2))//" "//Trim(c_nu(i2))//")"
     End Do
     BRtot = 0.5_dp * Sum(BrWqq)
     Do i2=1,3
      Write(io_L,202) Brtot * Abs(CKM(1,i2))**2 * BR(i_zaehl) &
        & , 3, id_w,-id_u(1),id_d(i2) &
        & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_u(1))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(1,i2))**2 * BR(i_zaehl) &
        & ,3,-id_w,id_u(1),-id_d(i2) &
        & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_u(1))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(2,i2))**2 * BR(i_zaehl) &
        & , 3, id_w,-id_u(2),id_d(i2) &
        & ,Trim(c_m)//" -> W+ W-* -> W+ "//Trim(c_u(2))//" "//Trim(c_d(i2))//")"
      Write(io_L,202) Brtot * Abs(CKM(2,i2))**2 * BR(i_zaehl) &
       & ,3,-id_w,id_u(2),-id_d(i2) &
       & ,Trim(c_m)//" -> W- W+* -> W- "//Trim(c_u(2))//" "//Trim(c_d(i2))//")"
     End Do
    End If
    i_zaehl = i_zaehl + 2
    If (BR(i_zaehl).Gt.BrMin) Then 
     Write(io_L,202) BrZinv*BR(i_zaehl),3,id_Z,id_nu(1),-id_nu(1), &
        & Trim(c_m)//" -> Z0 nu_i nu_i)"
     Do i2=1,3  
      Write(io_L,202) BrZll(i2)*BR(i_zaehl),3,id_Z,id_l(i2),-id_l(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_lm(i2))//" "//Trim(c_lp(i2))//")"
     End Do
     Do i2=1,3  
      Write(io_L,202) BrZqq(i2)*BR(i_zaehl),3,id_Z,id_d(i2),-id_d(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_d(i2))//" "//Trim(c_d(i2))//")"
     End Do
     Do i2=1,2 
      Write(io_L,202) BrZqq(i2)*BR(i_zaehl),3,id_Z,id_u(i2),-id_u(i2), &
        & Trim(c_m)//" -> Z0 "//Trim(c_u(i2))//" "//Trim(c_u(i2))//")"
     End Do
    End If
   End If

  End Do


  !----------------------
  ! neutral pseudoscalar
  !----------------------
  Do i1=1,n_p0
   c_m = c_p0(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    Do i2 = 1, 5 - n_c
     Do i3 = i2, 5 - n_c
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_lm(i2)//" "//c_lp(i3)
       Fnames(i_zaehl+1) = c_lp(i2)//" "//c_lm(i3)
       Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i3))
       Lnames(i_zaehl+1) = Trim(c_lm(i2))//" "//Trim(c_lP(i3))
       id_d_2(i_zaehl,1) = -id_l(i2)
       id_d_2(i_zaehl,2) = id_l(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,3
     Do i3 = i2,3
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i2)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i2)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+1) = c_d(i2)//" "//c_d(i3)
       Fnames(i_zaehl+9) = c_u(i2)//" "//c_u(i3)
       Fnames(i_zaehl+10) = c_u(i2)//" "//c_u(i3)
       Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+1) = Trim(c_d(i2))//" "//Trim(c_d(i3))
       Lnames(i_zaehl+9) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       Lnames(i_zaehl+10) = Trim(c_u(i2))//" "//Trim(c_u(i3))
       id_d_2(i_zaehl,1) = -id_d(i2)
       id_d_2(i_zaehl,2) = id_d(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       id_d_2(i_zaehl+9,1) = -id_u(i2)
       id_d_2(i_zaehl+9,2) = id_u(i3)
       id_d_2(i_zaehl+10,:) = - id_d_2(i_zaehl+9,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 9
    Do i2 = 1, n_sl
     Do i3 = i2, n_sl
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'slepton^-_'//Bu(i2)//' slepton^+_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'slepton^+_'//Bu(i2)//' slepton^-_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_slep(i2))//" "//Trim(c_sle(i3))
       Lnames(i_zaehl+1) = Trim(c_sle(i2))//" "//Trim(c_slep(i3))
       id_d_2(i_zaehl,1) = -id_sle(i2)
       id_d_2(i_zaehl,2) = id_sle(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    Do i2 = 1,6
     Do i3 = i2,6
      If (i2.Eq.i3) Then
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       i_zaehl = i_zaehl + 1
      Else
       Fnames(i_zaehl) = &
                    & 'd-squark_'//Bu(i2)//' d-squark^*_'//Bu(i3)
       Fnames(i_zaehl+1) = &
                    & 'd-squark^*_'//Bu(i2)//' d-squark_'//Bu(i3)
       Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       Lnames(i_zaehl+1) = Trim(c_sd(i2))//" "//Trim(c_sd(i3))
       id_d_2(i_zaehl,1) = -id_sd(i2)
       id_d_2(i_zaehl,2) = id_sd(i3)
       id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
       Fnames(i_zaehl+36) = &
                    & 'u-squark_'//Bu(i2)//' u-squark^*_'//Bu(i3)
       Fnames(i_zaehl+37) = &
                    & 'u-squark^*_'//Bu(i2)//' u-squark_'//Bu(i3)
       Lnames(i_zaehl+36) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       Lnames(i_zaehl+37) = Trim(c_su(i2))//" "//Trim(c_su(i3))
       id_d_2(i_zaehl+36,1) = -id_su(i2)
       id_d_2(i_zaehl+36,2) = id_su(i3)
       id_d_2(i_zaehl+37,:) = - id_d_2(i_zaehl+36,:)
       i_zaehl = i_zaehl + 2
      End If
     End Do
    End Do
    i_zaehl = i_zaehl + 36

   Else ! .not.GenerationMixing
    ! leptons
    Do i2=1,5 - n_c
     Lnames(i_zaehl) = Trim(c_lp(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = -id_l(i2)
     id_d_2(i_zaehl,2) = id_l(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! d-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_d(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = -id_d(i2)
     id_d_2(i_zaehl,2) = id_d(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! u-quarks
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_u(i2))
     id_d_2(i_zaehl,1) = -id_u(i2)
     id_d_2(i_zaehl,2) = id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    ! sleptons
    Do i2=1,5 - n_c
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_slep(i3))//" "//Trim(c_sle(i4))
        Lnames(i_zaehl+1) = Trim(c_sle(i3))//" "//Trim(c_slep(i4))
        id_d_2(i_zaehl,1) = -id_sle(i3)
        id_d_2(i_zaehl,2) = id_sle(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! d-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        Lnames(i_zaehl+1) = Trim(c_sd(i3))//" "//Trim(c_sd(i4))
        id_d_2(i_zaehl,1) = -id_sd(i3)
        id_d_2(i_zaehl,2) = id_sd(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2
    ! u-squarks
    Do i2=1,3
     Do i3=2*i2-1,2*i2
      Do i4=i3,2*i2
       If (i3.Eq.i4) Then
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        i_zaehl = i_zaehl + 1
       Else
        Lnames(i_zaehl) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        Lnames(i_zaehl+1) = Trim(c_su(i3))//" "//Trim(c_su(i4))
        id_d_2(i_zaehl,1) = -id_su(i3)
        id_d_2(i_zaehl,2) = id_su(i4)
        id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
        i_zaehl = i_zaehl + 2
       End If
      End Do ! i4
     End Do ! i3
    End Do ! i2

    If ((n_c-2).Eq.1) Then
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     i2 = 2
    Else If ((n_c-2).Eq.3) Then
     i2 = 0
    Else
     Fnames(1) = 'electrons '
     Fnames(2) = 'muons '
     Fnames(3) = 'taus '
     i2 = 3
    End If
    Fnames(i2+1) = 'd-quark '
    Fnames(i2+2) = 's-quark '
    Fnames(i2+3) = 'b-quark '
    Fnames(i2+4) = 'u-quark '
    Fnames(i2+5) = 'c-quark '
    Fnames(i2+6) = 't-quark '

    If ((n_c-2).Eq.1) Then
     Fnames(9) = 'Selectron 1 1 '
     Fnames(10) = 'Selectron 1 2 '
     Fnames(11) = 'Selectron 2 1 '
     Fnames(12) = 'Selectron 2 2 '
     Fnames(13) = 'Smuon 1 1 '
     Fnames(14) = 'Smuon 1 2 '
     Fnames(15) = 'Smuon 2 1 '
     Fnames(16) = 'Smuon 2 2 '
     i2 = 16
    Else If ((n_c-2).Eq.3) Then
     i2 = 6
    Else
     Fnames(10) = 'Selectron 1 1 '
     Fnames(11) = 'Selectron 1 2 '
     Fnames(12) = 'Selectron 2 1 '
     Fnames(13) = 'Selectron 2 2 '
     Fnames(14) = 'Smuon 1 1 '
     Fnames(15) = 'Smuon 1 2 '
     Fnames(16) = 'Smuon 2 1 '
     Fnames(17) = 'Smuon 2 2 '
     Fnames(18) = 'Stau 1 1 '
     Fnames(19) = 'Stau 1 2 '
     Fnames(20) = 'Stau 2 1 '
     Fnames(21) = 'Stau 2 2 '
     i2 = 21
    End If
    Fnames(i2+1) = 'Sdown 1 1 '
    Fnames(i2+2) = 'Sdown 1 2 '
    Fnames(i2+3) = 'Sdown 2 1 '
    Fnames(i2+4) = 'Sdown 2 2 '
    Fnames(i2+5) = 'S-strange 1 1 '
    Fnames(i2+6) = 'S-strange 1 2 '
    Fnames(i2+7) = 'S-strange 2 1 '
    Fnames(i2+8) = 'S-strange 2 2 '
    Fnames(i2+9) = 'Sbottom 1 1 '
    Fnames(i2+10) = 'Sbottom 1 2 '
    Fnames(i1+11) = 'Sbottom 2 1 '
    Fnames(i1+12) = 'Sbottom 2 2 '
    Fnames(i1+13) = 'Sup 1 1 '
    Fnames(i1+14) = 'Sup 1 2 '
    Fnames(i1+15) = 'Sup 2 1 '
    Fnames(i1+16) = 'Sup 2 2 '
    Fnames(i1+17) = 'S-charm 1 1 '
    Fnames(i1+18) = 'S-charm 1 2 '
    Fnames(i1+19) = 'S-charm 2 1 '
    Fnames(i1+20) = 'S-charm 2 2 '
    Fnames(i1+21) = 'Stop 1 1 '
    Fnames(i1+22) = 'Stop 1 2 '
    Fnames(i1+23) = 'Stop 2 1 '
    Fnames(i1+24) = 'Stop 2 2 '

   End If

   ! neutralinos 
   Do i2=1,n_n
    Do i3=i2,n_n
     Fnames(i_zaehl) = 'neutralino_'//Bu(i2)//' neutralino_'//Bu(i3)
     Lnames(i_zaehl) = Trim(c_c0(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = id_n(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do
   ! charginos
   Do i2=1, n_c
    Do i3=i2, n_c
     If (i2.Eq.i3) Then
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i2))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      i_zaehl = i_zaehl + 1
     Else
      Fnames(i_zaehl) = 'chargino^+_'//Bu(i2)//' chargino^-_'//Bu(i3)
      Fnames(i_zaehl+1) = 'chargino^-_'//Bu(i2)//' chargino^+_'//Bu(i3)
      Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_cp(i3))
      Lnames(i_zaehl+1) = Trim(c_cp(i2))//" "//Trim(c_cm(i3))
      id_d_2(i_zaehl,1) = - id_c(i2)
      id_d_2(i_zaehl,2) = id_c(i3)
      id_d_2(i_zaehl+1,:) = - id_d_2(i_zaehl,:)
      i_zaehl = i_zaehl + 2
     End If
    End Do
   End Do

   Do i2=1,n_spm
    Fnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Fnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    Lnames(i_zaehl+1) = 'W- '//Trim(c_sp(i2))
    Lnames(i_zaehl) = 'W+ '//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_W
    id_d_2(i_zaehl,2) = id_sm(i2)
    id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
    i_zaehl = i_zaehl + 2
   End Do
   Do i2=1,n_spm
    Fnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    Lnames(i_zaehl) = Trim(c_sp(i2))//Trim(c_sm(i2))
    id_d_2(i_zaehl,1) = id_sp(i2)
    id_d_2(i_zaehl,2) = id_sm(i2)
    i_zaehl = i_zaehl + 1
    Do i3=i2+1,n_spm
     Fnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Fnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     Lnames(i_zaehl) = Trim(c_sp(i2))//" "//Trim(c_sm(i3))
     Lnames(i_zaehl+1) = Trim(c_sm(i2))//" "//Trim(c_sp(i3))
     id_d_2(i_zaehl,1) = id_sp(i2)
     id_d_2(i_zaehl,2) = id_sm(i3)
     id_d_2(i_zaehl+1,:) = -id_d_2(i_zaehl,:)
     i_zaehl = i_zaehl + 2
    End Do
   End Do

   Do i2=1,n_s0
    Fnames(i_zaehl) = 'Z0 '//Trim(c_s0(i2))
    Lnames(i_zaehl) = 'Z0 '//Trim(c_s0(i2))
    id_d_2(i_zaehl,1) = id_Z
    id_d_2(i_zaehl,2) = id_s0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   Do i2=1,i1-1
    Do i3=1,n_s0
     Fnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_s0(i3))
     Lnames(i_zaehl) = Trim(c_p0(i2))//" "//Trim(c_s0(i3))
     id_d_2(i_zaehl,1) = id_P0(i2)
     id_d_2(i_zaehl,2) = id_s0(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do

   gT = gT_P0(i1+1)
   BR(1:250) = BR_P0(i1+1,:)
   gP(1:250) = gP_P0(i1+1,:)

   Write(io_L,200) id_p0(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    
   If (BR(i_zaehl+3).Gt.BrMin)  Write(io_L,201) BR(i_zaehl+3),2,id_gl  &
           &  ,id_gl,Trim(c_m)//" -> g g)"
 
   Call WriteDecays2(io,c_m , Fnames, gP, 100*BR, gT, BrMin100)

  End Do


  !-----------------
  ! charged scalars
  !-----------------
  Do i1=1,n_spm
   c_m = c_sm(i1)
   i_zaehl = 1

   If (GenerationMixing) Then
    !-----------------------------------------------------------------
    ! leptons, if RP is violated, these BRs are included in the
    ! chargino/neutralino final  states
    !-----------------------------------------------------------------
    Do i2=1,5 - n_c
     Lnames(i_zaehl) = Trim(c_nu(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = id_l(i2)
     id_d_2(i_zaehl,2) = -id_nu(i2)
     i_zaehl = i_zaehl + 1
    End Do

    !--------
    ! quarks
    !--------
    Do i2 = 1, 3
     Do i3 = 1,3
      Lnames(i_zaehl) = Trim(c_u(i3))//" "//Trim(c_d(i2))
      id_d_2(i_zaehl,1) = id_d(i2)
      id_d_2(i_zaehl,2) = -id_u(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do

    !-----------------------------------------------------------------
    ! sleptons, in case of RP violation -> S0 S-, P0 S- final states
    !-----------------------------------------------------------------
    Do i2 = 1,2*(5-n_c)
     Do i3 = 1,5-n_c
      Lnames(i_zaehl) = Trim(c_sle(i2))//" "//Trim(c_snu(i3))
      id_d_2(i_zaehl,1) = id_sle(i2)
      id_d_2(i_zaehl,2) = -id_snu(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do
    !------------------
    ! into squarks
    !------------------
    Do i2 = 1,6
     Do i3 = 1,6
      Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_su(i3))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_su(i3)
      i_zaehl = i_zaehl + 1
     End Do
    End Do

   Else ! no Generation Mixing
    Do i2=1,5-n_c
     Lnames(i_zaehl) = Trim(c_nu(i2))//" "//Trim(c_lm(i2))
     id_d_2(i_zaehl,1) = id_l(i2)
     id_d_2(i_zaehl,2) = -id_nu(i2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1,3
     Lnames(i_zaehl) = Trim(c_u(i2))//" "//Trim(c_d(i2))
     id_d_2(i_zaehl,1) = id_d(i2)
     id_d_2(i_zaehl,2) = -id_u(i2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1, 2*(5 - n_c)
     Lnames(i_zaehl) = Trim(c_sle(i2))//" "//Trim(c_snu((i2+1)/2))
     id_d_2(i_zaehl,1) = id_sle(i2)
     id_d_2(i_zaehl,2) = -id_snu((i2+1)/2)
     i_zaehl = i_zaehl + 1
    End Do
    Do i2=1,6
     If (i2.Le.2) i3=1
     If (i2.Le.4) i3=3
     If (i2.Le.6) i3=5
     Do i4=i3,i3+1
      Lnames(i_zaehl) = Trim(c_sd(i2))//" "//Trim(c_su(i4))
      id_d_2(i_zaehl,1) = id_sd(i2)
      id_d_2(i_zaehl,2) = -id_su(i4)
      i_zaehl = i_zaehl + 1
     End Do
    End Do
   End If

   Do i2=1,n_c
    Do i3=1,n_n
     Lnames(i_zaehl) = Trim(c_cm(i2))//" "//Trim(c_c0(i3))
     id_d_2(i_zaehl,1) = - id_c(i2)
     id_d_2(i_zaehl,2) = id_n(i3)
     i_zaehl = i_zaehl + 1
    End Do
   End Do

   Do i2=1,n_p0
    Lnames(i_zaehl) = " W- "//Trim(c_p0(i2))
    id_d_2(i_zaehl,1) = -id_W
    id_d_2(i_zaehl,2) = id_p0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   Do i2=1,n_s0
    Lnames(i_zaehl) = " W- "//Trim(c_s0(i2))
    id_d_2(i_zaehl,1) = -id_W
    id_d_2(i_zaehl,2) = id_s0(i2)
    i_zaehl = i_zaehl + 1
   End Do

   gT = gT_Spm(i1+1)
   BR(1:200) = BR_Spm(i1+1,:)
   gP(1:200) = gP_Spm(i1+1,:)

   Write(io_L,200) id_sm(i1),gT,c_m
   Write(io_L,100) "#    BR                NDA      ID1      ID2"
   Do i2=1,i_zaehl-1
    If (BR(i2).Gt.BrMin)  Write(io_L,201) BR(i2),2,id_d_2(i2,:)  &
           &  ,Trim(c_m)//" -> "//Trim(Lnames(i2))//")"
   End Do    

  End Do

 100 Format(a)
 200 Format("DECAY",1x,I9,3x,1P,E16.8,0P,3x,"# ",a)
 201 Format(3x,1P,e16.8,0p,3x,I2,3x,2(i9,1x),2x,"# BR(",a)
 202 Format(3x,1P,e16.8,0p,3x,I2,3x,3(i9,1x),2x,"# BR(",a)
 
 end subroutine LH_write_decays

 Subroutine LH_write_ext(io_L, Q, bi)
 Implicit None
  Integer, Intent(in) :: io_L
  Logical, Intent(in) :: bi
  Real(dp), intent(in) :: Q

   Write(io_L,100) "Block EXTPAR  # "
   Write(io_L,104) 0,Q , "# scale Q where the parameters below are defined"
   Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
   Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
   Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
   Write(io_L,104) 11,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t"
   Write(io_L,104) 12,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b"
   Write(io_L,104) 13,Real(A_l(3,3)/y_l(3,3),dp), "# A_l"
   If (bi) Write(io_L,104) 23,Real(mu ,dp), "# mu "
   Write(io_L,104) 31,Sqrt(Real(M2_L(1,1),dp)),"# M_(L,11)"
   Write(io_L,104) 32,Sqrt(Real(M2_L(2,2),dp)),"# M_(L,22)"
   Write(io_L,104) 33,Sqrt(Real(M2_L(3,3),dp)),"# M_(L,33)"
   Write(io_L,104) 34,Sqrt(Real(M2_E(1,1),dp)),"# M_(E,11)"
   Write(io_L,104) 35,Sqrt(Real(M2_E(2,2),dp)),"# M_(E,22)"
   Write(io_L,104) 36,Sqrt(Real(M2_E(3,3),dp)),"# M_(E,33)"
   Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
   Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
   Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
   Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
   Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
   Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
   Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
   Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
   Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"

 100 Format(a)
 104 Format(i4,2x,1P,e16.8,2x,a)
      
 End Subroutine LH_write_ext

 Subroutine LH_write_msq(io_L, c_sd, id_sd, c_su, id_su)
 Implicit None
  Integer, Intent(in) :: io_L
  Character(len=4), Intent(out) :: c_sd(6), c_su(6)
  Integer, Intent(out) :: id_sd(6), id_su(6)

  Integer :: i1


  If (GenerationMixing) Then

   id_sd(1) = 1000001
   id_sd(2) = 1000003
   id_sd(3) = 1000005
   id_sd(4) = 2000001
   id_sd(5) = 2000003
   id_sd(6) = 2000005
   Do i1=1,6
    c_sd(i1) = "~d_"//Bu(i1)
    Write(io_L,102) id_sd(i1),msdown(i1),"# "//Trim(c_sd(i1))
   End Do

   id_su(1) = 1000002
   id_su(2) = 1000004
   id_su(3) = 1000006
   id_su(4) = 2000002
   id_su(5) = 2000004
   id_su(6) = 2000006
   Do i1=1,6
    c_su(i1) = "~u_"//Bu(i1)
    Write(io_L,102) id_su(i1),msup(i1),"# "//Trim(c_su(i1))
   End Do

  Else ! .not.GenerationMixing

   If (Abs(rsdown(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000001,msdown(1),"# ~d_L"
    Write(io_L,102) 2000001,msdown(2),"# ~d_R"
    id_sd(1) = 1000001
    id_sd(2) = 2000001
    c_sd(1) = "~d_L"
    c_sd(2) = "~d_R"
   Else
    Write(io_L,102) 1000001,msdown(2),"# ~d_L"
    Write(io_L,102) 2000001,msdown(1),"# ~d_R"
    id_sd(2) = 1000001
    id_sd(1) = 2000001
    c_sd(2) = "~d_L"
    c_sd(1) = "~d_R"
   End If
   If (Abs(rsup(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000002,msup(1),"# ~u_L"
    Write(io_L,102) 2000002,msup(2),"# ~u_R"
    id_su(1) = 1000002
    id_su(2) = 2000002
    c_su(1) = "~u_L"
    c_su(2) = "~u_R"
   Else
    Write(io_L,102) 1000002,msup(2),"# ~u_L"
    Write(io_L,102) 2000002,msup(1),"# ~u_R"
    id_su(2) = 1000002
    id_su(1) = 2000002
    c_su(2) = "~u_L"
    c_su(1) = "~u_R"
   End If
   If (Abs(rsdown(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000003,msdown(3),"# ~s_L"
    Write(io_L,102) 2000003,msdown(4),"# ~s_R"
    id_sd(3) = 1000003
    id_sd(4) = 2000003
    c_sd(3) = "~s_L"
    c_sd(4) = "~s_R"
   Else
    Write(io_L,102) 1000003,msdown(4),"# ~s_L"
    Write(io_L,102) 2000003,msdown(3),"# ~s_R"
    id_sd(4) = 1000003
    id_sd(3) = 2000003
    c_sd(4) = "~s_L"
    c_sd(3) = "~s_R"
   End If
   If (Abs(rsup(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000004,msup(3),"# ~c_L"
    Write(io_L,102) 2000004,msup(4),"# ~c_R"
    id_su(3) = 1000004
    id_su(4) = 2000004
    c_su(3) = "~c_L"
    c_su(4) = "~c_R"
   Else
    Write(io_L,102) 1000004,msup(4),"# ~c_L"
    Write(io_L,102) 2000004,msup(3),"# ~c_R"
    id_su(4) = 1000004
    id_su(3) = 2000004
    c_su(4) = "~c_L"
    c_su(3) = "~c_R"
   End If
   Write(io_L,102) 1000005,msdown(5),"# ~b_1"
   Write(io_L,102) 2000005,msdown(6),"# ~b_2"
   id_sd(5) = 1000005
   id_sd(6) = 2000005
   c_sd(5) = "~b_1"
   c_sd(6) = "~b_2"
   Write(io_L,102) 1000006,msup(5),"# ~t_1"
   Write(io_L,102) 2000006,msup(6),"# ~t_2"
   id_su(5) = 1000006
   id_su(6) = 2000006
   c_su(5) = "~t_1"
   c_su(6) = "~t_2"
  End If ! GenerationMixing
  
 102 Format(1x,i9,3x,1P,e16.8,2x,a)

 End Subroutine LH_write_msq

 Subroutine LH_write_rsq(io_L, RUsq_ckm, RDsq_ckm)
 Implicit None
  Integer, Intent(in) :: io_L
  Complex(dp), Intent(in) :: RUsq_ckm(6,6), RDsq_ckm(6,6)
  Integer :: i1, i2

  If (generationmixing) Then
   Write(io_L,100) "Block USQmix  # u-sqark mixing matrix"
   Do i1=1,6
    Do i2=1,6
     Write(io_L,105) i1,i2,Real(RUsq_ckm(i1,i2),dp) &
                   & ,"# R_Su("//bu(i1)//","//bu(i2)//")"
    End Do
   End Do
   If (Maxval(Abs(Aimag((RUsq_ckm)))).Ne.0._dp) Then
    Write(io_L,100) "Block IMUSQmix  # imaginiary parts of u-sqark mixing matrix"
    Do i1=1,6
     Do i2=1,6
      Write(io_L,105) i1,i2,Aimag(RUsq_ckm(i1,i2)) &
                    & ,"# Im[R_Su("//bu(i1)//","//bu(i2)//")]"
     End Do
    End Do
   End If
   Write(io_L,100) "Block DSQmix  # d-squark mixing matrix"
   Do i1=1,6
    Do i2=1,6
     Write(io_L,105) i1,i2,Real(RDsq_ckm(i1,i2),dp) &
                   & ,"# R_Sd("//bu(i1)//","//bu(i2)//")"
    End Do
   End Do
   If (Maxval(Abs(Aimag((RDsq_ckm)))).Ne.0._dp) Then
    Write(io_L,100) "Block IMDSQmix  # imaginiary parts of d-sqark mixing matrix"
    Do i1=1,6
     Do i2=1,6
      Write(io_L,105) i1,i2,Aimag(Rdsq_ckm(i1,i2)) &
                    & ,"# Im[R_Sd("//bu(i1)//","//bu(i2)//")]"
     End Do
    End Do
   End If

  Else ! .not.GenerationMixing

   Write(io_L,100) "Block stopmix  # stop mixing matrix"
   Write(io_L,105) 1,1,Real(RUsq_ckm(5,5),dp),"# R_st(1,1)"
   Write(io_L,105) 1,2,Real(RUsq_ckm(5,6),dp),"# R_st(1,2)"
   Write(io_L,105) 2,1,Real(RUsq_ckm(6,5),dp),"# R_st(2,1)"
   Write(io_L,105) 2,2,Real(RUsq_ckm(6,6),dp),"# R_st(2,2)"
   Write(io_L,100) "Block sbotmix  # sbottom mixing matrix"
   Write(io_L,105) 1,1,Real(RDsq_ckm(5,5),dp),"# R_sb(1,1)"
   Write(io_L,105) 1,2,Real(RDsq_ckm(5,6),dp),"# R_sb(1,2)"
   Write(io_L,105) 2,1,Real(RDsq_ckm(6,5),dp),"# R_sb(2,1)"
   Write(io_L,105) 2,2,Real(RDsq_ckm(6,6),dp),"# R_sb(2,2)"

  End If

 100 Format(a)
 105 Format(1x,2i3,3x,1P,e16.8,3x,a)

 End Subroutine LH_write_rsq

 Subroutine PutUpperCase(name)
 Implicit None
  Character(len=80), Intent(inout) :: name
  Integer :: len=80, i1
  Do i1=1,len
   If (name(i1:i1).Eq."a") name(i1:i1) = "A"
   If (name(i1:i1).Eq."b") name(i1:i1) = "B"
   If (name(i1:i1).Eq."c") name(i1:i1) = "C"
   If (name(i1:i1).Eq."d") name(i1:i1) = "D"
   If (name(i1:i1).Eq."e") name(i1:i1) = "E"
   If (name(i1:i1).Eq."f") name(i1:i1) = "F"
   If (name(i1:i1).Eq."g") name(i1:i1) = "G"
   If (name(i1:i1).Eq."h") name(i1:i1) = "H"
   If (name(i1:i1).Eq."i") name(i1:i1) = "I"
   If (name(i1:i1).Eq."j") name(i1:i1) = "J"
   If (name(i1:i1).Eq."k") name(i1:i1) = "K"
   If (name(i1:i1).Eq."l") name(i1:i1) = "L"
   If (name(i1:i1).Eq."m") name(i1:i1) = "M"
   If (name(i1:i1).Eq."n") name(i1:i1) = "N"
   If (name(i1:i1).Eq."o") name(i1:i1) = "O"
   If (name(i1:i1).Eq."p") name(i1:i1) = "P"
   If (name(i1:i1).Eq."q") name(i1:i1) = "Q"
   If (name(i1:i1).Eq."r") name(i1:i1) = "R"
   If (name(i1:i1).Eq."s") name(i1:i1) = "S"
   If (name(i1:i1).Eq."t") name(i1:i1) = "T"
   If (name(i1:i1).Eq."u") name(i1:i1) = "U"
   If (name(i1:i1).Eq."v") name(i1:i1) = "V"
   If (name(i1:i1).Eq."w") name(i1:i1) = "W"
   If (name(i1:i1).Eq."x") name(i1:i1) = "X"
   If (name(i1:i1).Eq."y") name(i1:i1) = "Y"
   If (name(i1:i1).Eq."z") name(i1:i1) = "Z"
  End Do
 End Subroutine PutUpperCase


  Subroutine ReadMatrixC(io, nmax, mat, ic, mat_name, kont, fill)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax, io, ic
   Integer, Intent(in), Optional :: fill
   Complex(dp), Intent(inout) :: mat(nmax, nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadMatrixC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -308
     Call AddError(308)
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 1
     kont = -308
     Call AddError(308)
     Call TerminateProgram()
    End If

    If (ic.Eq.0) Then
     mat(i1,i2) = Cmplx(0._dp,Aimag(mat(i1,i2)),dp) + wert
     If (Present(fill).And.(i1.Ne.i2)) &
       &  mat(i2,i1) = Cmplx(0._dp, Aimag(mat(i2,i1)), dp) + wert
    Else If (ic.Eq.1) Then
     mat(i1,i2) = Real(mat(i1,i2),dp) + Cmplx(0._dp, wert, dp)
     !-------------------------------------------------------------
     ! if fill==1 -> matrix is hermitian
     ! if fill==2 -> matrix is complex symmetric
     !-------------------------------------------------------------
     If (Present(fill).And.(i1.Ne.i2)) Then
      If (fill.Eq.1) mat(i2,i1) = Real(mat(i2,i1),dp) - Cmplx(0._dp, wert, dp)
      If (fill.Eq.2) mat(i2,i1) = Real(mat(i2,i1),dp) + Cmplx(0._dp, wert, dp)
     End If
    End If

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadMatrixC

  Subroutine ReadMatrixR(io, nmax, mat, mat_name, kont)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax, io
   Real(dp), Intent(inout) :: mat(nmax, nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadMatrixR"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -309
     Call AddError(309) 
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 1
     kont = -309
     Call AddError(309) 
     Call TerminateProgram()
    End If

    mat(i1,i2) = wert

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadMatrixR
  
  Subroutine ReadTensorC(io, nmax, mat, ic, mat_name, kont)
  Implicit None
   Character(len=*) :: mat_name
   Integer, Intent(in) :: nmax, io, ic
   Complex(dp), Intent(inout) :: mat(nmax, nmax, nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1, i2, i3
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadTensorC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, i2, i3, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -312
     Call AddError(312)
     Call TerminateProgram()
    End If
    If ((i2.Lt.1).Or.(i2.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i2=",i2
     Iname = Iname - 1
     kont = -312 
     Call AddError(312)
     Call TerminateProgram()
    End If
    If ((i3.Lt.1).Or.(i3.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//mat_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i3=",i3
     Iname = Iname - 2
     kont = -312
      Call AddError(312)
    Call TerminateProgram()
    End If

    If (ic.Eq.0) mat(i1,i2,i3) = Cmplx(0._dp, Aimag(mat(i1,i2,i3)), dp) + wert
    If (ic.Eq.1) mat(i1,i2,i3) = mat(i1,i2,i3) + Cmplx(0._dp, wert, dp)

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadTensorC
  
  Subroutine ReadVectorC(io, nmax, vec, ic, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io, ic
   Complex(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorC"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -310
     Call AddError(310)
     Call TerminateProgram()
    End If

    If (ic.Eq.0) vec(i1) = Cmplx(0._dp, Aimag(vec(i1)), dp) + wert
    If (ic.Eq.1) vec(i1) = Real(vec(i1),dp) + Cmplx(0._dp, wert, dp)

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadVectorC
  
  Subroutine ReadVectorR(io, nmax, vec, vec_name, kont)
  Implicit None
   Character(len=*) :: vec_name
   Integer, Intent(in) :: nmax, io
   Real(dp), Intent(inout) :: vec(nmax)
   Integer, Intent(out) :: kont

   Character(len=80) :: read_line
   Integer :: i1
   Real(dp) :: wert

   kont = 0

   Iname = Iname + 1
   NameOfUnit(Iname) = "ReadVectorR"
   Do 
    Read(io,*,End=200) read_line
!     Write(*,*) read_line
    If (read_line(1:1).Eq."#") Cycle ! ignore comments
    Backspace(io)                    ! resetting to the beginning of the line
    If ((read_line(1:1).Eq."B").Or.(read_line(1:1).Eq."b") ) Then
     Iname = Iname - 1
     Return ! new block
    End If

    Read(io,*) i1, wert, read_line

    If ((i1.Lt.1).Or.(i1.Gt.nmax)) Then
     Write(ErrCan,*) "Problem while reading "//vec_name//" in routine"// &
        & Trim(NameOfUnit(Iname))//", index i1=",i1
     Iname = Iname - 1
     kont = -311
     Call AddError(311)
     Call TerminateProgram()
    End If

    vec(i1) = wert

   End Do

   200 Iname = Iname - 1
   Return

  End Subroutine ReadVectorR

 Subroutine SetWriteMinBR(wert)
 !-------------------------------------------------------------------
 ! sets the minimal branching ratio (=wert) appearing in the output
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: wert
  BrMin = wert
 End Subroutine SetWriteMinBR


 Subroutine SetWriteMinSig(wert)
 !-------------------------------------------------------------------
 ! sets the minimal cross section (=wert) appearing in the output
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: wert
  SigMin = wert
 End Subroutine SetWriteMinSig


 Subroutine WriteSPhenoOutputLHA1(io_L, m_GUT)
 Implicit None

  Integer, Intent(in) :: io_L
  real(dp) :: m_GUT

  Integer :: kont, n_run, i1,i2,ierr
   Real(dp) :: g, gp, delta
  !------------------------------
  ! masses and mixing angles
  !------------------------------
  Complex(dp), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L, uU_R
  Logical :: WriteOut, check

  Integer, Save  :: f_c=1, k1, k2, k3, k4
  Character(len=11) :: f_n
  Logical :: non_minimal
  Integer :: i3
  Real(dp) :: Yu(3), Yd(3), Q, MaxCont, mat6r(6,6), nr(4,4), mnr(4)
  Complex(dp), Dimension(3,3) :: CKM_Q
  Complex(dp), Dimension(6,6) :: RUsq_ckm, RDsq_ckm
  Character(len=10) :: name, c_m, c_snu(3), c_sle(6)     &
      & , c_su(6), c_sd(6), c_slep(6), c_grav, zeit
  Character(len=6) :: c_lm(3), c_lp(3), c_nu(3), c_gl, c_cp(5) &
     & , c_cm(5), c_c0(8), c_phot, c_sp(7), c_sm(7)
  Character(len=1) :: c_u(3), c_d(3)
  Character(len=4) :: c_S0(6), c_P0(5)
  Character(len=8) :: Datum
  Character(len=13) :: c_sfermion 
  Character(len=80) :: name2
  Integer :: n_n, n_c, n_s0, n_p0, n_spm, n_sl, n_sn, n_sd, n_su, ds 
  Integer :: id_n(8), id_sp(7), id_sm(7), id_S0(6), id_P0(5), id_c(5) &
      & , id_snu(3), i_zaehl, id_sle(6), id_f, id_fp, id_su(6), ii, jj &
      & , id_sd(6), pos
  Integer, Parameter ::id_A0 = 36, id_Hp = 37                            &
      & , id_W = 24, id_Z = 23, id_ph = 22, id_gl = 21                   &
      & , id_l(3) = (/ 11, 13, 15 /), id_nu(3) = (/ 12, 14, 16 /)        &
      & , id_u(3) = (/ 2, 4, 6 /), id_d(3) = (/ 1, 3, 5 /)               &
      & , id_grav = 1000039, id_glu = 1000021
  Logical :: use_flavour_states=.True.
  Real(dp) :: T3, YL, YR, ML, MR, mSf(2), mSf2(2)
  Complex(dp) :: Rsf(2,2), A, Y

  Iname = Iname + 1
  NameOfUnit(Iname) = "Calculate_Omega_h_sq"


  c_phot = "photon"
  c_grav = "~Gravitino"
  c_lm(1) = "e-"
  c_lp(1) = "e+"
  c_lm(2) = "mu-"
  c_lp(2) = "mu+"
  c_lm(3) = "tau-"
  c_lp(3) = "tau+"

  If (GenerationMixing) Then
   c_nu(1) = "nu_1"
   c_nu(2) = "nu_2"
   c_nu(3) = "nu_3"
  Else
   c_nu(1) = "nu_e"
   c_nu(2) = "nu_mu"
   c_nu(3) = "nu_tau"
  End If
  c_d(1) = "d"
  c_d(2) = "s"
  c_d(3) = "b"
  c_u(1) = "u"
  c_u(2) = "c"
  c_u(3) = "t"
  c_gl = "~g"
  c_cp(1) = "chi_1+"
  c_cp(2) = "chi_2+"
  c_cm(1) = "chi_1-"
  c_cm(2) = "chi_2-"
  c_c0(1) = "chi_10"
  c_c0(2) = "chi_20"
  c_c0(3) = "chi_30"
  c_c0(4) = "chi_40"
  id_n(1) = 1000022
  id_n(2) = 1000023
  id_n(3) = 1000025
  id_n(4) = 1000035

  c_sp(1) = "H+"
  c_sm(1) = "H-"
  id_sp = 0
  id_sp(1) = 37
  id_sm = - id_sp  

  n_sd = 6
  n_su = 6
 
  id_c(1) = 1000024
  id_c(2) = 1000037

   c_s0(1) = "h0"
   c_s0(2) = "H0"
   id_S0(1) = 25
   id_S0(2) = 35
   c_P0(1) = "A0"
   id_P0(1) = 36
   n_n = 4
   n_c = 2
   n_s0 = 2
   n_p0 = 1 
   n_spm = 1 
   n_sl = 6
   n_sn = 3

  Q = Sqrt( GetRenormalizationScale() )
  gp = gauge(1)
  g = gauge(2)
  Open(io_L,file="SPheno_1.spc",status="unknown")
  !--------------------------------------------------------
  ! General information
  !--------------------------------------------------------
  ! Les Houches standard
  !-----------------------
   Write(io_L,100) "# SUSY Les Houches Accord 2.beta - MSSM spectrum + Decays"
   Write(io_L,100) "# SPheno 3.xx beta"
   Write(io_L,100) &
     & "# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101"
   Write(io_L,100) "# in case of problems send email to porod@ific.uv.es"
   Write(io_L,100) "# Created: today"
   Write(io_L,100) "Block SPINFO         # Program information"
   Write(io_L,100) "     1   SPheno      # spectrum calculator"
   Write(io_L,100) "     2   3.xx beta   # version number"

   !--------------------------------------
   ! model information
   !--------------------------------------
   If (HighScaleModel.Eq."mSugra") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    Write(io_L,100) "#"

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
   Else If ((HighScaleModel(1:5).Eq."SUGRA").Or.   &
          & (HighScaleModel(1:5).Eq."SEESA")) Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    1    # mSUGRA model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,100) "    3    2    # mSUGRA model + SU(5)"
    Else If (HighScaleModel.Eq."SUGRA_NuR1") Then
     Write(io_L,100) "    3    3    # mSUGRA model + nu_R at a common scale"
     Write(io_L,100) "Block MnuR  # mass scale of the right handed neutrinos"
     Write(io_L,101) 1,MnuR(1),"# m_nu_R      "     
    Else If (HighScaleModel.Eq."SUGRA_NuR") Then
     Y_nu_0 = Transpose(Y_nu_0) ! in the RGEs the transposed Yukawas are used
     Write(io_L,100) "    3    4    # mSUGRA model + three  nu_R, ~nu_R"
     Write(io_L,106) "Block Ynu0 Q=",m_GUT,"# (GUT scale)"
     Do i1=1,3
      Do i2=1,3
       Write(io_L,105) i2,i1,Real(Y_nu_0(i2,i1),dp), &
             & "# Y_(nu,"//bu(i2)//bu(i1)//")"  
      End Do
     End Do
     If (Maxval(Abs(Aimag(Y_nu_0))).Gt.0._dp) Then
      Write(io_L,106) "Block IMYnu0 Q=",m_GUT,"# (GUT scale)"
      Do i1=1,3
       Do i2=1,3
        Write(io_L,105) i2,i1,Aimag(Y_nu_0(i2,i1)) &
             & ,"# Im(Y_(nu,"//bu(i2)//bu(i1)//"))"  
       End Do
      End Do
     End If
    Else If (HighScaleModel.Eq."SEESAW_II") Then
     Write(io_L,100) "    3    5    # mSUGRA model + Higgs triplett"
     Write(io_L,106) "Block YT0 Q=",m_GUT,"# (GUT scale)"
     Do i1=1,3
      Do i2=i1,3
       Write(io_L,105) i2,i1,Real(Y_T_0(i2,i1),dp),"# Y_(T,"//bu(i2)//bu(i1)//")"  
      End Do
     End Do
     If (Maxval(Abs(Aimag(Y_T_0))).Gt.0._dp) Then
      Write(io_L,106) "Block IMYT0 Q=",m_GUT,"# (GUT scale)"
      Do i1=1,3
       Do i2=i1,3
        Write(io_L,105) i2,i1,Aimag(Y_T_0(i2,i1)) &
             & ,"# Im(Y_(T,"//bu(i2)//bu(i1)//"))"  
       End Do
      End Do
     End If
     Write(io_L,106) "Block Higgs3 Q=",m_GUT,"# (GUT scale)"
     Write(io_L,101) 1,M_H3(1),"# m_H3      "
     Write(io_L,101) 2,Real(lam12_0(1),dp),"# Re(lambda_1)"
     If (Aimag(lam12_0(1)).Ne.0._dp) &
         Write(io_L,101) 3,Aimag(lam12_0(1)),"# Im(lambda_1)"
     Write(io_L,101) 4,Real(lam12_0(2),dp),"# Re(lambda_2)"
     If (Aimag(lam12_0(2)).Ne.0._dp) &
         Write(io_L,101) 5,Aimag(lam12_0(2)),"# Im(lambda_2)"
     If (Fifteen_plet) Then
      Write(io_L,100) "    6    1               # using RGEs for 15-plet"
     Else
      Write(io_L,100) "    6    0               # using RGEs for 3-plet"
     End If
      
    End If
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,Sqrt(Real(M2_E_0(1,1),dp)),"# m0      "
    Write(io_L,101) 2,Real(Mi_0(1),dp),"# m12     "
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
    Write(io_L,101) 5,Real(AoY_l_0(1,1),dp),"# A0"
    If (HighScaleModel.Eq."SUGRA_SU5") Then
     Write(io_L,101) 7,Max(M_SO_10,m_GUT),"# SO(10) scale"
     Write(io_L,101) 8,D_SO_10,"# D-terms at SO(10) scale"
    End If
    Write(io_L,100) "#"

    Write(io_L,106) "Block gauge Q=",m_GUT,"# (GUT scale)"
    Write(io_L,104) 1,gauge_0(1),"# g'(Q)^DRbar"
    Write(io_L,104) 2,gauge_0(2),"# g(Q)^DRbar"
    Write(io_L,104) 3,gauge_0(3),"# g3(Q)^DRbar"
    
   Else If (HighScaleModel.Eq."AMSB") Then
    non_minimal = .False.
    Write(io_L,100) "Block MODSEL  # Model selection"
    Write(io_L,100) "    1    3    # mAMSB model"
    If (GenerationMixing) Write(io_L,100) &
      & " 6 1                      # switching on flavour violation"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 1,m0_amsb,"# M_0" 
    Write(io_L,101) 2,m_32,"# m_3/2, gravitino mass"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z   "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"

    
   Else
    non_minimal = .True.
    Write(io_L,100) "# Either the general MSSM or a model has been used"
    Write(io_L,100) &
      & "# which has not yet been implemented in the LesHouches standard"
    Write(io_L,100) "Block MINPAR  # Input parameters"
    Write(io_L,101) 3,tanb_mZ,"# tanb at m_Z    "
    Write(io_L,101) 4, Real(phase_mu,dp),"# Sign(mu)"
   End If 

   !---------------------------------------------------
   ! parameters + masses for SPheno.spc
   !---------------------------------------------------
   Write(io_L,100) "Block SMINPUTS  # SM parameters"
   Write(io_L,102) 1, 1._dp / alpha_MSbar(mZ, mW),"# alpha_em^-1(MZ)^MSbar"
   Write(io_L,102) 2,G_F,"# G_mu [GeV^-2]"
   Write(io_L,102) 3,alphaS_MZ,"# alpha_s(MZ)^MSbar"
   Write(io_L,102) 4,mZ,"# m_Z(pole)"
   Write(io_L,102) 5,mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) 6,mf_u(3),"# m_t(pole)"
   Write(io_L,102) 7,mf_l(3),"# m_tau(pole)"

   !----------------------------------------------------------------
   ! in the case of GenerationMixing all parameters and mixings
   ! are given in the SuperCKM basis for squarks
   !----------------------------------------------------------------
   If (GenerationMixing) Then
    Write(io_L,100) "Block VCKMIN  # CKM matrix, Wolfenstein parameterization"
    Write(io_L,102) 1, lam_wolf,"# lambda"
    Write(io_L,102) 2,A_wolf,"# A"
    Write(io_L,102) 3,rho_wolf,"# rho bar"
    Write(io_L,102) 4,eta_wolf,"# eta bar"
    Call Switch_to_superCKM(Y_d, Y_u, A_d, A_u, M2_D, M2_Q, M2_U         &
              &, Ad_sckm, Au_sckm, M2D_sckm, M2Q_sckm, M2U_sckm, .False. &
              &, RSdown, RSup, Rdsq_ckm, RUsq_ckm, CKM_Q, Yd, Yu )

   Else ! .non.GenerationMixing

    Do i1=1,3
     Yu(i1) = Real(Y_u(i1,i1),dp)
     Yd(i1) = Real(Y_d(i1,i1),dp)
    End Do
    Ad_sckm = A_d
    Au_sckm = A_u

    M2D_SCKM = M2_D
    M2U_SCKM = M2_U
    M2Q_SCKM = M2_Q

    RUsq_ckm = RSup
    RDsq_ckm = RSdown

   End If

   If (non_minimal) Then
    Write(io_L,100) "Block EXTPAR  # "
    Write(io_L,104) 0,Q , "# scale Q where the parameters below are defined"
    Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
    Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
    Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
    Write(io_L,104) 11,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t"
    Write(io_L,104) 12,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b"
    Write(io_L,104) 13,Real(A_l(3,3)/y_l(3,3),dp), "# A_l"
    Write(io_L,104) 23,Real(mu ,dp), "# mu "
    Write(io_L,104) 31,Sqrt(Real(M2_L(1,1),dp)),"# M_(L,11)"
    Write(io_L,104) 32,Sqrt(Real(M2_L(2,2),dp)),"# M_(L,22)"
    Write(io_L,104) 33,Sqrt(Real(M2_L(3,3),dp)),"# M_(L,33)"
    Write(io_L,104) 34,Sqrt(Real(M2_E(1,1),dp)),"# M_(E,11)"
    Write(io_L,104) 35,Sqrt(Real(M2_E(2,2),dp)),"# M_(E,22)"
    Write(io_L,104) 36,Sqrt(Real(M2_E(3,3),dp)),"# M_(E,33)"
    Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
    Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
    Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
    Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
    Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
    Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
    Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
    Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
    Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"
   End If
      
! couplings
  Write(io_L,106) "Block gauge Q=",Q,"# (SUSY scale)"
  Write(io_L,104) 1,gauge(1),"# g'(Q)^DRbar"
  Write(io_L,104) 2,gauge(2),"# g(Q)^DRbar"
  Write(io_L,104) 3,gauge(3),"# g3(Q)^DRbar"

  Write(io_L,106) "Block Yu Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yu(1), "# Y_u(Q)^DRbar"
  Write(io_L,107) 2,2,Yu(2), "# Y_c(Q)^DRbar"
  Write(io_L,107) 3,3,Yu(3), "# Y_t(Q)^DRbar"

  Write(io_L,106) "Block Yd Q=",Q,"# (SUSY scale)"
  Write(io_L,107) 1,1,Yd(1), "# Y_d(Q)^DRbar"
  Write(io_L,107) 2,2,Yd(2), "# Y_s(Q)^DRbar"
  Write(io_L,107) 3,3,Yd(3), "# Y_b(Q)^DRbar"

  If (GenerationMixing) Then 
                             
   Write(io_L,106) "Block VCKM Q=",Q,"# Re(CKM) at the SUSY scale"
   Do i1=1,3
    Do i2=1,3
     Write(io_L,107) i1,i2,Real(CKM_Q(i1,i2),dp),"# Re(V_"//Bu(i1)//Bu(i2)//")"
    End Do
   End Do
   If (Maxval(Abs(Aimag(CKM_Q))).Gt.0._dp) Then
    Write(io_L,106) "Block IMVCKM Q=",Q,"# Im(CKM) at the SUSY scale"
    Do i1=1,3
     Do i2=1,3
      Write(io_L,107) i1,i2,Aimag(CKM_Q(i1,i2)),"# Im(V_"//Bu(i1)//Bu(i2)//")"
     End Do
    End Do
   End If

  End If ! generationmixing

  Write(io_L,106) "Block Ye Q=",Q,"# (SUSY scale)"

  ierr = 0
  If (GenerationMixing) Then
   !-------------------------------------------
   ! check if any off-diagonal term is non-zero
   !-------------------------------------------
   Do i1=1,3
    Do i2=1,3
     If ((i1.Ne.i2).And.(Abs(Y_l(i2,i1)).Ne.0._dp)) ierr = ierr + 1
    End Do
   End Do
  End If

  If (ierr.Ne.0) Then
   Do i1=1,3
    Do i2=1,3
     Write(io_L,105) i2,i1,Real(Y_l(i2,i1),dp),"# Y_(l,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
  Else 
   Write(io_L,107) 1,1,Real(y_l(1,1),dp), "# Y_e(Q)^DRbar"
   Write(io_L,107) 2,2,Real(y_l(2,2),dp), "# Y_mu(Q)^DRbar"
   Write(io_L,107) 3,3,Real(y_l(3,3),dp), "# Y_tau(Q)^DRbar"
  End If
  If (Maxval(Abs(Aimag(Y_l))).Gt.0._dp) Then
   Write(io_L,106) "Block IMYe Q=",Q,"# (SUSY scale)"
   If (GenerationMixing) Then
    Do i1=1,3
      Do i2=1,3
      Write(io_L,105) i2,i1,Aimag(Y_l(i2,i1)),"# Im(Y_(l,"//bu(i2)//bu(i1)//"))"
     End Do
    End Do
   Else 
    Write(io_L,107) 1,1,Aimag(y_l(1,1)), "# Im(Y_e)(Q)^DRbar"
    Write(io_L,107) 2,2,Aimag(y_l(2,2)), "# Im(Y_mu)(Q)^DRbar"
    Write(io_L,107) 3,3,Aimag(y_l(3,3)), "# Im(Y_tau)(Q)^DRbar"
   End If
  End If

   Write(io_L,106) "Block Au Q=",Q,"# (SUSY scale)"
   Write(io_L,107) 1,1,Real(Au_sckm(1,1)/y_u(1,1),dp), "# A_u(Q)^DRbar"
   Write(io_L,107) 2,2,Real(Au_sckm(2,2)/y_u(2,2),dp), "# A_c(Q)^DRbar"
   Write(io_L,107) 3,3,Real(Au_sckm(3,3)/y_u(3,3),dp), "# A_t(Q)^DRbar"

   Write(io_L,106) "Block Ad Q=",Q,"# (SUSY scale)"
   Write(io_L,107) 1,1,Real(Ad_sckm(1,1)/y_d(1,1),dp), "# A_d(Q)^DRbar"
   Write(io_L,107) 2,2,Real(Ad_sckm(2,2)/y_d(2,2),dp), "# A_s(Q)^DRbar"
   Write(io_L,107) 3,3,Real(Ad_sckm(3,3)/y_d(3,3),dp), "# A_b(Q)^DRbar"

   Write(io_L,106) "Block Ae Q=",Q,"# (SUSY scale)"
   Write(io_L,107) 1,1,Real(A_l(1,1)/y_l(1,1),dp), "# A_e(Q)^DRbar"
   Write(io_L,107) 2,2,Real(A_l(2,2)/y_l(2,2),dp), "# A_mu(Q)^DRbar"
   Write(io_L,107) 3,3,Real(A_l(3,3)/y_l(3,3),dp), "# A_tau(Q)^DRbar"


  Write(io_L,106) "Block MSOFT Q=",Q,"# soft SUSY breaking masses at Q"
  Write(io_L,104) 1,Real(Mi(1),dp),"# M_1"
  Write(io_L,104) 2,Real(Mi(2),dp),"# M_2"
  Write(io_L,104) 3,Real(Mi(3),dp),"# M_3"
  Write(io_L,104) 21,M2_H(1),"# M^2_(H,d)"
  Write(io_L,104) 22,M2_H(2),"# M^2_(H,u)"

  Write(io_L,104) 31,Sqrt(Real(M2_L(1,1),dp)),"# M_(L,11)"
  Write(io_L,104) 32,Sqrt(Real(M2_L(2,2),dp)),"# M_(L,22)"
  Write(io_L,104) 33,Sqrt(Real(M2_L(3,3),dp)),"# M_(L,33)"
  Write(io_L,104) 34,Sqrt(Real(M2_E(1,1),dp)),"# M_(E,11)"
  Write(io_L,104) 35,Sqrt(Real(M2_E(2,2),dp)),"# M_(E,22)"
  Write(io_L,104) 36,Sqrt(Real(M2_E(3,3),dp)),"# M_(E,33)"
  Write(io_L,104) 41,Sqrt(Real(M2Q_SCKM(1,1),dp)),"# M_(Q,11)"
  Write(io_L,104) 42,Sqrt(Real(M2Q_SCKM(2,2),dp)),"# M_(Q,22)"
  Write(io_L,104) 43,Sqrt(Real(M2Q_SCKM(3,3),dp)),"# M_(Q,33)"
  Write(io_L,104) 44,Sqrt(Real(M2U_SCKM(1,1),dp)),"# M_(U,11)"
  Write(io_L,104) 45,Sqrt(Real(M2U_SCKM(2,2),dp)),"# M_(U,22)"
  Write(io_L,104) 46,Sqrt(Real(M2U_SCKM(3,3),dp)),"# M_(U,33)"
  Write(io_L,104) 47,Sqrt(Real(M2D_SCKM(1,1),dp)),"# M_(D,11)"
  Write(io_L,104) 48,Sqrt(Real(M2D_SCKM(2,2),dp)),"# M_(D,22)"
  Write(io_L,104) 49,Sqrt(Real(M2D_SCKM(3,3),dp)),"# M_(D,33)"

  If (GenerationMixing) Then
   Write(io_L,106) "Block MSL2 Q=",Q,"# M^2_L soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2_L(i2,i1),dp)  &
                  & ,"# M^2_(L,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2_L))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSL2 Q=",Q  &
                  & ,"# Im(M^2_L) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2_L(i2,i1)) &
                 &   ,"# Im(M^2_(L,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSE2 Q=",Q,"# M^2_E soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2_E(i2,i1),dp)  &
                  & ,"# M^2_(E,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2_E))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSE2 Q=",Q  &
                  & ,"# Im(M^2_E) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2_E(i2,i1)) &
                 &   ,"# Im(M^2_(E,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSQ2 Q=",Q,"# M^2_Q soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2Q_SCKM(i2,i1),dp) &
                   & ,"# M^2_(Q,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2Q_SCKM))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSQ2 Q=",Q  &
                  & ,"# Im(M^2_Q) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2Q_SCKM(i2,i1)) &
                 &   ,"# Im(M^2_(Q,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSU2 Q=",Q,"# M^2_U soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2U_SCKM(i1,i2),dp)  &
                  & ,"# M^2_(U,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2U_SCKM))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSU2 Q=",Q  &
                  & ,"# Im(M^2_U) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2U_SCKM(i1,i2)) &
                 &   ,"# Im(M^2_(U,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

   Write(io_L,106) "Block MSD2 Q=",Q,"# M^2_D soft SUSY breaking masses at Q"
   Do i2=1,3
    Do i1=1,3
     Write(io_L,105) i2,i1,Real(M2D_SCKM(i1,i2),dp)  &
                  & ,"# M^2_(D,"//bu(i2)//bu(i1)//")"  
    End Do
   End Do
   If (Maxval(Abs(Aimag(M2D_SCKM))).Gt.0._dp) Then
    Write(io_L,106) "Block IMMSD2 Q=",Q  &
                  & ,"# Im(M^2_D) soft SUSY breaking masses at Q"
    Do i2=1,3
     Do i1=1,3
      Write(io_L,105) i2,i1,Aimag(M2D_SCKM(i1,i2)) &
                 &   ,"# Im(M^2_(D,"//bu(i2)//bu(i1)//"))"  
     End Do
    End Do
   End If

  End If

  If (HighScaleModel.Eq."NMSSM") Then
   Write(io_L,104) 61,Real(h0,dp),"# lambda"
   Write(io_L,104) 62,0.5_dp*Real(lam,dp),"# kappa"
   Write(io_L,104) 63,Real(ao_h0,dp),"# A_lambda"
   Write(io_L,104) 64,Real(Ao_lam,dp),"# A_kappa"
   Write(io_L,104) 65,Real(oosqrt2*vP*h0,dp),"# mu_eff"

  Else If (HighScaleModel.Eq."RPexplicit") Then
   Write(io_L,106) "Block RVKAPPA Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(eps(1),dp),"# epsilon_1"
   Write(io_L,102) 2,Real(eps(2),dp),"# epsilon_2"
   Write(io_L,102) 3,Real(eps(3),dp),"# epsilon_3"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVKAPPA Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(eps(1)),"# Im(epsilon_1)"
    Write(io_L,102) 2,Aimag(eps(2)),"# Im(epsilon_2)"
    Write(io_L,102) 3,Aimag(eps(3)),"# Im(epsilon_3)"
   End If

   Write(io_L,106) "Block RVD Q=",Q,"# bilinear RP parameters at Q"
   Write(io_L,102) 1,Real(Beps(1),dp),"# Re( B_1 epsilon_1)"
   Write(io_L,102) 2,Real(Beps(2),dp),"# Re( B_2 epsilon_2)"
   Write(io_L,102) 3,Real(Beps(3),dp),"# Re( B_3 epsilon_3)"
   If (Maxval(Abs(Aimag(eps))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVD Q=",Q,"# imaginary parts of bilinear RP parameters at Q"
    Write(io_L,102) 1,Aimag(Beps(1)),"# Im( B_1 epsilon_1)"
    Write(io_L,102) 2,Aimag(Beps(2)),"# Im( B_2 epsilon_2)"
    Write(io_L,102) 3,Aimag(Beps(3)),"# Im( B_3 epsilon_3)"
   End If

   Write(io_L,106) "Block RVSNVEV Q=",Q,"# sneutrino vevs at Q"
   Write(io_L,102) 1,vevL(1),"# v_L_1"
   Write(io_L,102) 2,vevL(2),"# v_L_2"
   Write(io_L,102) 3,vevL(3),"# v_L_3"

   If (Maxval(Abs(Rp_lam)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDA Q=",Q,"# lambda_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lam(i1,i2,i3),dp) &
                     & ,"# lambda_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDA Q=",Q,"# Im(lambda_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lam(i1,i2,i3)) &
            & ,"# Im(lambda_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   If (Maxval(Abs(Rp_lamp)).Gt.0._dp) Then
    Write(io_L,106) "Block RVLAMPBDAP Q=",Q,"# lambda'_ijk at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Real(Rp_lamp(i1,i2,i3),dp) &
                     & ,"# lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)
      End Do
     End Do
    End Do
   End If
   If (Maxval(Abs(Aimag(Rp_lam))).Gt.0._dp) Then
    Write(io_L,106) "Block IMRVLAMPBDAP Q=",Q,"# Im(lambda'_ijk) at Q"
    Do i1=1,3
     Do i2=1,3
      Do i3=1,3
       Write(io_L,109) i1,i2,i3,Aimag(Rp_lamp(i1,i2,i3)) &
            & ,"# Im(lambda'_"//Bu(i1)//Bu(i2)//Bu(i3)//")"
      End Do
     End Do
    End Do
   End If

   Write(io_L,100) "Block SPhenoRP  # additional RP parameters"
   Write(io_L,102) 4,Lam_Ex(1),"# Lambda_1 = v_d epsilon_1 + mu v_L1"
   Write(io_L,102) 5,Lam_Ex(2),"# Lambda_2 = v_d epsilon_2 + mu v_L2"
   Write(io_L,102) 6,Lam_Ex(3),"# Lambda_3 = v_d epsilon_3 + mu v_L3"
   Write(io_L,102) 7,(mN7(3)**2-mN7(2)**2)*1.e18_dp,"# m^2_atm [eV^2]"
   Write(io_L,102) 8,(mN7(2)**2-mN7(1)**2)*1.e18_dp,"# m^2_sol [eV^2]"
   Write(io_L,102) 9,Abs(N7(3,6)/N7(3,7))**2,"# tan^2 theta_atm"
   Write(io_L,102) 10,Abs(N7(2,5)/N7(1,5))**2,"# tan^2 theta_sol"
   Write(io_L,102) 11,Abs(N7(3,5))**2,"# U_e3^2"
   Write(io_L,102) 15,vevSM(1),"# v_d"
   Write(io_L,102) 16,vevSM(2),"# v_u"
  End If


   Write(io_L,100) "Block MASS  # Mass spectrum"
   Write(io_L,100) "#   PDG code      mass          particle"
   Write(io_L,102) id_u(2),mf_u(2),"# m_c(m_c), MSbar"
   Write(io_L,102) id_d(3),mf_d(3),"# m_b(m_b), MSbar"
   Write(io_L,102) id_u(3),mf_u(3),"# m_t(pole)"
   Write(io_L,102) id_Z,mZ,"# m_Z(pole)"
   Write(io_L,102) id_W,mW,"# W+"
   If (HighScaleModel.Ne."RPexplicit") &
        & Write(io_L,102) id_l(3),mf_l(3),"# m_tau(pole)"
   If (HighScaleModel.Eq."NMSSM") Then
    Write(io_L,102) 25,mS03(1),"# leightest neutral scalar" 
    Write(io_L,102) 35,mS03(2),"# second neutral scalar" 
    Write(io_L,102) 45,mS03(3),"# third neutral scalar" 
    Write(io_L,102) 36,mP03(2),"# leighter pseudoscalar" 
    Write(io_L,102) 46,mP03(3),"# heavier pseudoscalar" 
    Write(io_L,102) 37,mSpm(2),"# H+"
   Else If (HighScaleModel.Eq."RPexplicit") Then
!    Write(io_L,102) 12,mN7(1),"# lightest neutrino"
!    Write(io_L,102) 14,mN7(2),"# second lightest neutrino"
!    Write(io_L,102) 16,mN7(3),"# heaviest neutrino"
    Write(io_L,102) id_s0(1),mS05(1),"# leightest neutral scalar" 
    Write(io_L,102) id_s0(2),mS05(2),"# 2nd neutral scalar" 
    Write(io_L,102) id_s0(3),mS05(3),"# 3rd neutral scalar" 
    Write(io_L,102) id_s0(4),mS05(4),"# 4th neutral scalar" 
    Write(io_L,102) id_s0(5),mS05(5),"# 5th neutral scalar" 
    Write(io_L,102) id_P0(1),mP05(2),"# leightest pseudoscalar" 
    Write(io_L,102) id_P0(2),mP05(3),"# 2nd pseudoscalar" 
    Write(io_L,102) id_P0(3),mP05(4),"# 3rd pseudoscalar" 
    Write(io_L,102) id_p0(4),mP05(5),"# 4th pseudoscalar"
    Write(io_L,102) Abs(id_sp(1)),mSpm8(2),"# leightest charged scalar" 
    Write(io_L,102) Abs(id_sp(2)),mSpm8(3),"# 2nd charged scalar" 
    Write(io_L,102) Abs(id_sp(3)),mSpm8(4),"# 3rd charged scalar" 
    Write(io_L,102) Abs(id_sp(4)),mSpm8(5),"# 4th charged scalar" 
    Write(io_L,102) Abs(id_sp(5)),mSpm8(6),"# 5th charged scalar" 
    Write(io_L,102) Abs(id_sp(6)),mSpm8(7),"# 6th charged scalar" 
    Write(io_L,102) Abs(id_sp(7)),mSpm8(8),"# 7th charged scalar" 
    
   Else
    Write(io_L,102) 25,mS0(1),"# h0" 
    Write(io_L,102) 35,mS0(2),"# H0" 
    Write(io_L,102) 36,mP0(2),"# A0" 
    Write(io_L,102) 37,mSpm(2),"# H+"
   End If
! squarks

  If (GenerationMixing) Then
   If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
     i_zaehl = 1
     mat6R = Abs(RDsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_sd(ii) = 1000001
       c_sd(ii) = "~d_L"
      Case(2)
       id_sd(ii) = 1000003
       c_sd(ii) = "~s_L-"
      Case(4)
       id_sd(ii) = 2000001
       c_sd(ii) = "~d_R"
      Case(5)
       id_sd(ii) = 2000003
       c_sd(ii) = "~s_R"
      Case default
       id_sd(ii) = 1000000 * i_zaehl + 5
       c_sd(ii) = "~b_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6
      Write(io_L,102) id_sd(ii),msdown(ii),"# "//Trim(c_sd(ii))
     End Do

     i_zaehl = 1
     mat6R = Abs(RUsq_ckm)
     Do i1=1,6
      MaxCont = Maxval(mat6R)
      Call FindPosition(6, mat6R, MaxCont, ii, jj)
      Select Case(jj)
      Case(1)
       id_su(ii) = 1000002
       c_su(ii) = "~u_L"
      Case(2)
       id_su(ii) = 1000004
       c_su(ii) = "~c_L-"
      Case(4)
       id_su(ii) = 2000002
       c_su(ii) = "~u_R"
      Case(5)
       id_su(ii) = 2000004
       c_su(ii) = "~c_R"
      Case default
       id_su(ii) = 1000000 * i_zaehl + 6
       c_su(ii) = "~t_"//bu(i_zaehl)//"-"
       i_zaehl = I_zaehl + 1
      End Select
      mat6R(ii,:) = 0._dp
      mat6R(:,jj) = 0._dp
     End Do
     Do ii=1,6
      Write(io_L,102) id_su(ii),msup(ii),"# "//Trim(c_su(ii))
     End Do
   
    Else ! use mass ordering

     id_sd(1) = 1000001
     id_sd(2) = 1000003
     id_sd(3) = 1000005
     id_sd(4) = 2000001
     id_sd(5) = 2000003
     id_sd(6) = 2000005
     Do i1=1,6
      c_sd(i1) = "~d_"//Bu(i1)
      Write(io_L,102) id_sd(i1),msdown(i1),"# "//Trim(c_sd(i1))
     End Do

     id_su(1) = 1000002
     id_su(2) = 1000004
     id_su(3) = 1000006
     id_su(4) = 2000002
     id_su(5) = 2000004
     id_su(6) = 2000006
     Do i1=1,6
      c_su(i1) = "~u_"//Bu(i1)
      Write(io_L,102) id_su(i1),msup(i1),"# "//Trim(c_su(i1))
     End Do
    End If
  Else ! .not.GenerationMixing

   If (Abs(rsdown(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000001,msdown(1),"# ~d_L"
    Write(io_L,102) 2000001,msdown(2),"# ~d_R"
    id_sd(1) = 1000001
    id_sd(2) = 2000001
    c_sd(1) = "~d_L"
    c_sd(2) = "~d_R"
   Else
    Write(io_L,102) 1000001,msdown(2),"# ~d_L"
    Write(io_L,102) 2000001,msdown(1),"# ~d_R"
    id_sd(2) = 1000001
    id_sd(1) = 2000001
    c_sd(2) = "~d_L"
    c_sd(1) = "~d_R"
   End If
   If (Abs(rsup(1,1)).Gt.0.5_dp) Then
    Write(io_L,102) 1000002,msup(1),"# ~u_L"
    Write(io_L,102) 2000002,msup(2),"# ~u_R"
    id_su(1) = 1000002
    id_su(2) = 2000002
    c_su(1) = "~u_L"
    c_su(2) = "~u_R"
   Else
    Write(io_L,102) 1000002,msup(2),"# ~u_L"
    Write(io_L,102) 2000002,msup(1),"# ~u_R"
    id_su(2) = 1000002
    id_su(1) = 2000002
    c_su(2) = "~u_L"
    c_su(1) = "~u_R"
   End If
   If (Abs(rsdown(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000003,msdown(3),"# ~s_L"
    Write(io_L,102) 2000003,msdown(4),"# ~s_R"
    id_sd(3) = 1000003
    id_sd(4) = 2000003
    c_sd(3) = "~s_L"
    c_sd(4) = "~s_R"
   Else
    Write(io_L,102) 1000003,msdown(4),"# ~s_L"
    Write(io_L,102) 2000003,msdown(3),"# ~s_R"
    id_sd(4) = 1000003
    id_sd(3) = 2000003
    c_sd(4) = "~s_L"
    c_sd(3) = "~s_R"
   End If
   If (Abs(rsup(3,3)).Gt.0.5_dp) Then
    Write(io_L,102) 1000004,msup(3),"# ~c_L"
    Write(io_L,102) 2000004,msup(4),"# ~c_R"
    id_su(3) = 1000004
    id_su(4) = 2000004
    c_su(3) = "~c_L"
    c_su(4) = "~c_R"
   Else
    Write(io_L,102) 1000004,msup(4),"# ~c_L"
    Write(io_L,102) 2000004,msup(3),"# ~c_R"
    id_su(4) = 1000004
    id_su(3) = 2000004
    c_su(4) = "~c_L"
    c_su(3) = "~c_R"
   End If
   Write(io_L,102) 1000005,msdown(5),"# ~b_1"
   Write(io_L,102) 2000005,msdown(6),"# ~b_2"
   id_sd(5) = 1000005
   id_sd(6) = 2000005
   c_sd(5) = "~b_1"
   c_sd(6) = "~b_2"
   Write(io_L,102) 1000006,msup(5),"# ~t_1"
   Write(io_L,102) 2000006,msup(6),"# ~t_2"
   id_su(5) = 1000006
   id_su(6) = 2000006
   c_su(5) = "~t_1"
   c_su(6) = "~t_2"
  End If ! GenerationMixing
   
   If (HighScaleModel.Ne."RPexplicit") Then

! sleptons

    If (GenerationMixing) Then
     If (Use_Flavour_States) Then ! using flavour ordering, old fashionnd
 
      Do i1=1,3
       MaxCont = Maxval(Abs(Rsneut(i1,:))**2)
       If (Abs(Rsneut(i1,1))**2.Eq.maxCont) Then
        id_snu(i1) = 1000012
        c_snu(i1) = "~nu_eL"
       Else If (Abs(Rsneut(i1,2))**2.Eq.maxCont) Then
        id_snu(i1) = 1000014
        c_snu(i1) = "~nu_muL"
       Else
        id_snu(i1) = 1000016
        c_snu(i1) = "~nu_tauL"
       End If
       Write(io_L,102) id_snu(i1),msneut(i1),"# "//Trim(c_snu(i1))
      End Do

      i_zaehl = 1
      mat6R = Abs(RSlepton)
      Do i1=1,6
       MaxCont = Maxval(mat6R)
       Call FindPosition(6, mat6R, MaxCont, ii, jj)
       Select Case (jj)
       Case(1)
        id_sle(ii) = 1000011
        c_sle(ii) = "~e_L-"
        c_slep(ii) = "~e_L+"
       Case(2)
        id_sle(ii) = 1000013
        c_sle(ii) = "~mu_L-"
        c_slep(ii) = "~mu_L+"
       Case(4)
        id_sle(ii) = 2000011
        c_sle(ii) = "~e_R-"
        c_slep(ii) = "~e_R+"
       Case(5)
        id_sle(ii) = 2000013
        c_sle(ii) = "~mu_R-"
        c_slep(ii) = "~mu_R+"
       Case default
        id_sle(ii) = 1000000 * i_zaehl + 15
        c_sle(ii) = "~tau_"//bu(i_zaehl)//"-"
        c_slep(ii) = "~tau_"//bu(i_zaehl)//"+"
        i_zaehl = I_zaehl + 1
       End Select
       mat6R(ii,:) = 0._dp
       mat6R(:,jj) = 0._dp
      End Do
      Do ii=1,6
       Write(io_L,102) id_sle(ii),mslepton(ii),"# "//Trim(c_sle(ii))
      End Do

     Else ! using mass ordering

      id_snu(1) = 1000012
      id_snu(2) = 1000014
      id_snu(3) = 1000016
      Do i1=1,3
       c_snu(i1) = "~nu_"//Bu(i1)
       Write(io_L,102) id_snu(i1),mSneut(i1),"# "//Trim(c_snu(i1))
      End Do

      id_sle(1) = 1000011
      id_sle(2) = 1000013
      id_sle(3) = 1000015
      id_sle(4) = 2000011
      id_sle(5) = 2000013
      id_sle(6) = 2000015
      Do i1=1,6
       c_sle(i1) = "~l_"//Bu(i1)
       c_slep(i1) = "~l+_"//Bu(i1)
       Write(io_L,102) id_sle(i1),mSlepton(i1),"# "//Trim(c_sle(i1))
      End Do
     End If

    Else ! .not.GenerationMixing

     id_snu = (/ 1000012, 1000014, 1000016 /)
     If (Abs(rslepton(1,1)).Gt.0.5_dp) Then
      Write(io_L,102) 1000011,mslepton(1),"# ~e_L-"
      Write(io_L,102) 2000011,mslepton(2),"# ~e_R-"
      id_sle(1) = 1000011
      id_sle(2) = 2000011
      c_sle(1) = "~e_L-"
      c_sle(2) = "~e_R-"
      c_slep(1) = "~e_L+"
      c_slep(2) = "~e_R+"
     Else
      Write(io_L,102) 1000011,mslepton(2),"# ~e_L-"
      Write(io_L,102) 2000011,mslepton(1),"# ~e_R-"
      id_sle(2) = 1000011
      id_sle(1) = 2000011
      c_sle(2) = "~e_L-"
      c_sle(1) = "~e_R-"
      c_slep(2) = "~e_L+"
      c_slep(1) = "~e_R+"
     End If
     Write(io_L,102) 1000012,msneut(1),"# ~nu_eL"
     c_snu(1) = "~nu_eL"
     If (Abs(rslepton(3,3)).Gt.0.5_dp) Then
      Write(io_L,102) 1000013,mslepton(3),"# ~mu_L-"
      Write(io_L,102) 2000013,mslepton(4),"# ~mu_R-"
      id_sle(3) = 1000013
      id_sle(4) = 2000013
      c_sle(3) = "~mu_L-"
      c_sle(4) = "~mu_R-"
      c_slep(3) = "~mu_L+"
      c_slep(4) = "~mu_R+"
     Else
      Write(io_L,102) 1000013,mslepton(4),"# ~mu_L-"
      Write(io_L,102) 2000013,mslepton(3),"# ~mu_R-"
      id_sle(4) = 1000013
      id_sle(3) = 2000013
      c_sle(4) = "~mu_L-"
      c_sle(3) = "~mu_R-"
      c_slep(4) = "~mu_L+"
      c_slep(3) = "~mu_R+"
     End If
     Write(io_L,102) 1000014,msneut(2),"# ~nu_muL"
     c_snu(2) = "~nu_muL"
     Write(io_L,102) 1000015,mslepton(5),"# ~tau_1-"
     Write(io_L,102) 2000015,mslepton(6),"# ~tau_2-"
     id_sle(5) = 1000015
     id_sle(6) = 2000015
     c_sle(5) = "~tau_1-"
     c_sle(6) = "~tau_2-"
     c_slep(5) = "~tau_1+"
     c_slep(6) = "~tau_2+"
     Write(io_L,102) 1000016,msneut(3),"# ~nu_tauL"
     c_snu(3) = "~nu_tauL"
    End If
   End If ! GenerationMixing

 ! gauginos/higgsinos
   Write(io_L,102) 1000021,mglu,"# ~g"
   ! checking for negative sign

    Do i1=1,4
     If (Sum(Abs(Real(N(i1,:)))).Lt.0.1_dp) Then
      mNr(i1) = - mN(i1)
      nr(i1,1:4) = Aimag(n(i1,:))
     Else   
      mNr(i1) =  mN(i1)
      nr(i1,1:4) = Real(n(i1,:),dp)
     End If
    End Do

    Write(io_L,102) 1000022,mnr(1),"# ~chi_10" 
    Write(io_L,102) 1000023,mnr(2),"# ~chi_20" 
    Write(io_L,102) 1000025,mnr(3),"# ~chi_30" 
    Write(io_L,102) 1000035,mnr(4),"# ~chi_40" 
    Write(io_L,102) 1000024,mc(1),"# ~chi_1+" 
    Write(io_L,102) 1000037,mc(2),"# ~chi_2+"

   If (Maxval(MnuR).Gt.0._dp) Then
    Write(io_L,100) "# masses of right handed neutrinos"
    Write(io_L,100) "Block MnuR"
    Write(io_L,102) 1,mnur(1),"# m_nu_R_1"
    Write(io_L,102) 2,mnur(2),"# m_nu_R_2"
    Write(io_L,102) 3,mnur(3),"# m_nu_R_3"
   End If

! Mixing matrices 
  Write(io_L,100) "# Higgs mixing"
   Write(io_L,100) "Block alpha # Effective Higgs mixing angle"
   Write(io_L,108) -Asin(RS0(1,1)),"# alpha"
   Write(io_L,103) "Block Hmix Q=",Q, "# Higgs mixing parameters"
   Write(io_L,104) 1,Real(mu,dp),"# mu"
   Write(io_L,104) 2,tanb_Q,"# tan[beta](Q)"
   Write(io_L,104) 3,vev_Q,"# v(Q)"
   Write(io_L,104) 4,mA2_Q,"# m^2_A(Q)"

!--------------------
! add correct order
!--------------------


   T3 = 0.5_dp
   Yl = 1._dp / 3._dp
   Yr = -4._dp / 3._dp
   Ml = M2Q_sckm(3,3)
   Mr = M2U_sckm(3,3)
   A = Au_sckm(3,3)
   Y = Yu(3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
   Write(io_L,100) "Block stopmix  # stop mixing matrix"
   Write(io_L,105) 1,1,Real(Rsf(1,1),dp),"# R_st(1,1)"
   Write(io_L,105) 1,2,Real(Rsf(1,2),dp),"# R_st(1,2)"
   Write(io_L,105) 2,1,Real(Rsf(2,1),dp),"# R_st(2,1)"
   Write(io_L,105) 2,2,Real(Rsf(2,2),dp),"# R_st(2,2)"

   T3 = -0.5_dp
   Yl = 1._dp / 3._dp
   Yr = 2._dp / 3._dp
   Ml = M2Q_sckm(3,3)
   Mr = M2d_sckm(3,3)
   A = Ad_sckm(3,3)
   Y = Yd(3)
   Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
   Write(io_L,100) "Block sbotmix  # sbottom mixing matrix"
   Write(io_L,105) 1,1,Real(Rsf(1,1),dp),"# R_sb(1,1)"
   Write(io_L,105) 1,2,Real(Rsf(1,2),dp),"# R_sb(1,2)"
   Write(io_L,105) 2,1,Real(Rsf(2,1),dp),"# R_sb(2,1)"
   Write(io_L,105) 2,2,Real(Rsf(2,2),dp),"# R_sb(2,2)"

   T3 = -0.5_dp
   Yl = -1._dp 
   Yr = 2._dp
    Ml = M2_L(3,3)
    Mr = M2_E(3,3)
    A = A_l(3,3)
    Y = Y_l(3,3)
    Call SfermionMass(Ml, Mr, A, mu, vevSM, Y, T3, Yl, Yr, g, gp, kont &
                    &, msf, msf2, Rsf)
    Write(io_L,100) "Block staumix  # stau mixing matrix"
   Write(io_L,105) 1,1,Real(Rsf(1,1),dp),"# R_sta(1,1)"
   Write(io_L,105) 1,2,Real(Rsf(1,2),dp),"# R_sta(1,2)"
   Write(io_L,105) 2,1,Real(Rsf(2,1),dp),"# R_sta(2,1)"
   Write(io_L,105) 2,2,Real(Rsf(2,2),dp),"# R_sta(2,2)"

   Write(io_L,100) "Block Nmix  # neutralino mixing matrix"
   Do i1=1,n_n
    Do i2=1,n_n
     name = "# N("//Char(48+i1)//","//Char(48+i2)//")"
     Write(io_L,105) i1,i2,nr(i1,i2),name
    End Do
   End Do

   Write(io_L,100) "Block Umix  # chargino U mixing matrix"
   Write(io_L,105) 1,1,Real(U(1,1),dp),"# U(1,1)"
   Write(io_L,105) 1,2,Real(U(1,2),dp),"# U(1,2)"
   Write(io_L,105) 2,1,Real(U(2,1),dp),"# U(2,1)"
   Write(io_L,105) 2,2,Real(U(2,2),dp),"# U(2,2)"
   Write(io_L,100) "Block Vmix  # chargino V mixing matrix"
   Write(io_L,105) 1,1,Real(V(1,1),dp),"# V(1,1)"
   Write(io_L,105) 1,2,Real(V(1,2),dp),"# V(1,2)"
   Write(io_L,105) 2,1,Real(V(2,1),dp),"# V(2,1)"
   Write(io_L,105) 2,2,Real(V(2,2),dp),"# V(2,2)"

 99 Format(1x,i5,3x,a)
100 Format(a)
101 Format(2x,i3,2x,1P,e16.8,2x,a) 
102 Format(1x,i9,3x,1P,e16.8,2x,a)
103 Format(a13,1P,e16.8,2x,a)
104 Format(i4,2x,1P,e16.8,2x,a)
105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)
107 Format(2i3,3x,1P,e16.8,3x,a)
108 Format(9x,1P,E16.8,0P,3x,a)
109 Format(1x,3i3,3x,1P,e16.8,3x,a)
  Close(io_L)
  

!   if (f_c.lt.10) then
!    f_n = "omg000"//bu(f_c)//".out"
!   else if (f_c.lt.100) then
!    k1 = f_c/10
!    f_n = "omg00"//bu(k1)//bu(f_c-k1*10)//".out"
!   else if (f_c.lt.1000) then
!    k1 = f_c/100
!    k2 = (f_c - k1 * 100)/10
!    f_n = "omg0"//bu(k1)//bu(k2)//bu(f_c-k1*100-k2*10)//".out"
!   else
!    k1 = f_c/1000
!    k2 = (f_c - k1 * 1000)/100
!    k3 = (f_c - k1 * 1000 - k2*100)/10
!    f_n = "omg"//bu(k1)//bu(k2)//bu(k3)//bu(f_c-k1*1000-k2*100-k3*10)//".out"
!   end if   

  Iname = Iname - 1
  
 End Subroutine WriteSPhenoOutputLHA1


 Subroutine WriteComplexMatrix(n, matrix, OnlyDiagonal)
 !---------------------------------------------------------------
 ! simplifies the writing of complex matrices in various places
 ! written by Werner Porod
 ! 15.11.01
 !---------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n
  Complex(dp), Intent(in) :: matrix(:,:)
  Logical, Intent(in) :: OnlyDiagonal

  Integer :: i1, i2, dim_matrix
  Complex(dp), Allocatable :: vec(:)

  dim_matrix = Size( matrix, dim=1)

  If (OnlyDiagonal) Then
   Allocate( vec(dim_matrix) )
   Do i1=1,dim_matrix
    vec(i1) = matrix(i1,i1)
   End Do
   If (Maxval( Abs( Aimag(vec) ) ).Eq.0._dp) Then
    Write(n,*) Real(vec)
   Else
    Write(n,*) vec
   End If

   Deallocate( vec )
    
  Else

   If (Maxval( Abs( Aimag(matrix) )) .Lt.1.e-15_dp) Then
    If (dim_matrix.Le.4) Then
     Do i1=1,dim_matrix
      If (dim_matrix.Eq.2) Write(n,902) Real(matrix(i1,:))
      If (dim_matrix.Eq.3) Write(n,903) Real(matrix(i1,:))
      If (dim_matrix.Eq.4) Write(n,904) Real(matrix(i1,:))
     End Do
    Else
     Do i1=1,dim_matrix
      Do i2=1,dim_matrix
       Write(n,901) i1,i2,Real(matrix(i1,i2))
      End Do
     End Do
    End If
   Else
    If (dim_matrix.Eq.2) Then
     Do i1=1,dim_matrix
      Write(n,802) matrix(i1,:)
     End Do
    Else
     Do i1=1,dim_matrix
      Do i2=1,dim_matrix
       Write(n,801) i1,i2,matrix(i1,i2)
      End Do
     End Do
    End If
   End If
  End If

 801 Format(t2,2i4,"  (",f16.5,",",f16.5,")")
 802 Format(t2,"(",f16.5,",",f16.5,")     (",f16.5,",",f16.5,")")


 901 Format(t2,2i4,f16.5)
 902 Format(t2,2f16.5)
 903 Format(t2,3f16.5)
 904 Format(t2,4f16.5)

 End Subroutine WriteComplexMatrix


 Subroutine WriteCrossSection(io, name, sigma, sigmin)
 Implicit None
  Integer, Intent(in) :: io
  Real(dp), Intent(in) :: sigma(:,:), sigmin
  Character(len=*), Intent(in) :: name

  Integer :: len1, len2, i1, i2

  len1= Size(sigma, dim=1)
  len2= Size(sigma, dim=2)

  If (Maxval(sigma).Eq.0._dp) Then
    Write(io,*) " "//name//" : kinematically not possible"
   Else If (Maxval(sigma).Lt.SigMin) Then
    Write(io,*) " "//name//" : cross sections below",SigMin,"fb"
   Else 
    Write(io,*) " "//name
    Do i1=1,len1
     Do i2=i1,len2
      If (Sigma(i1,i2).Gt.SigMin) Write(io,900) i1,i2,Sigma(i1,i2) 
     End Do
    End Do
   End If
  Write(io,*) " "

900 Format(2i4,f16.7," fb")

 End Subroutine WriteCrossSection


 Subroutine WriteCrossSection1(io, name, sigma, sigmin)
 Implicit None
  Integer, Intent(in) :: io
  Real(dp), Intent(in) :: sigma, sigmin
  Character(len=*), Intent(in) :: name

  If (sigma.Eq.0._dp) Then
    Write(io,*) " "//name//" : kinematically not possible"
   Else If (sigma.Lt.SigMin) Then
    Write(io,*) " "//name//" : cross section below",Real(SigMin),"fb"
   Else 
    Write(io,*) " "//name
    Write(io,900) sigma
   End If
  Write(io,*) " "

900 Format(t9,f16.7," fb")

 End Subroutine WriteCrossSection1

 
 Subroutine WriteDecays2(n, titel, names, gP, BR, gT, prez)
 Implicit None
  Integer, Intent(in) :: n
  Character (len=*), Intent(in) :: titel
  Character (len=30), Intent(in) :: names(:)
  Real(dp), Intent(in) :: gP(:), BR(:), gT, prez

  Integer :: i1, dim1

  If (sum(gp).Eq.0._dp) Then! case of stable particles
   Write (n,*) titel//" : stable"
    Write (n,*) ' '
   Return 
  End If

  dim1 = Size( gP )

  If (sum(gP).Lt.prez/100._dp) Then
  ! precision smaller than smallest partial width -> rescaling
   Write (n,*) titel
   Write(n,*) " Attention, the widths are given in eV!"
   Do i1=1,dim1
    If (BR(i1).Ge.prez) Then
     Write (n,100) names(i1), 1.e9_dp*gP(i1), BR(i1)
    End If
   End Do
   If (gT.gt.0._dp) then ! to allow the mixing of 2- and 3-body decays
    Write (n,101) 'Total width :                  ',1.e9_dp*gT
    Write (n,*) ' '
   End if
  Else
   Write (n,*) titel
   Do i1=1,dim1
    If (BR(i1).Ge.prez) Then
     Write (n,100) names(i1), gP(i1), BR(i1)
    End If
   End Do
   If (gT.gt.0._dp) then ! to allow the mixing of 2- and 3-body decays
    Write (n,101) 'Total width :                  ',gT
    Write (n,*) ' '
   End if
  End If

 100 Format(t3,a31,2f14.8)
 101 Format(t3,a31,f14.8)

 End Subroutine WriteDecays2



 Subroutine WriteMassesMSSM(io, WithComments, mGlu, PhaseGlu, mC, U, V, mN, N &
           &, mSneut, Rsneut, mSlepton, RSlepton, mSdown, RSdown, mSup, RSup  &
           &, mP0, mS0, RS0, mSpm, GenerationMixing)
 
 Implicit None
  Integer, Intent(in) :: io
  Logical, Intent(in) :: WithComments, GenerationMixing
  Real(dp), Intent(in) :: mP0(2), mS0(2), RS0(2,2), mSpm(2), mC(2)  &
    & , mN(4), mGlu, mSneut(3), mSlepton(6), mSdown(6), mSup(6)
 Complex(dp), Intent(in) :: PhaseGlu, U(2,2), V(2,2), N(4,4), Rsneut(3,3)     &
    & , RSlepton(6,6), RSdown(6,6), RSup(6,6)

 If (WithComments) Write(io,*) "Masses and mixing matrices"
 If (WithComments) Then
  If (Aimag(phaseglu).Eq.0._dp) Then
   Write(io,*) "Gluino : ",mglu, Real(phaseglu)
  Else
   Write(io,*) "Gluino : ",mglu, phaseglu
  End If
 Else
  If (Aimag(phaseglu).Eq.0._dp) Then
   Write(io,*) mglu, Real(phaseglu)
  Else
   Write(io,*) mglu, phaseglu
  End If
 End If
 If (WithComments) Write(io,*) " "

 If (WithComments) Write(io,*) "Charginos"
 Write(io,*) mC
 If (WithComments) Write(io,*) "  U"
 Call WriteComplexMatrix(io, U, .False.)
 If (WithComments) Write(io,*) "  V"
 Call WriteComplexMatrix(io, V, .False.)
 If (WithComments) Write(io,*) " "

 If (WithComments) Write(io,*) "Neutralinos"
 Write(io,*) mN
 If (WithComments) Write(io,*) "  N"
 Call WriteComplexMatrix(io, N, .False.)
 If (WithComments) Write(io,*) " "

 If (GenerationMixing) Then
  If (WithComments) Write(io,*) "sneutrinos"
  Write(io,*) mSneut
  If (WithComments) Write(io,*) "  R_sn"
  Call WriteComplexMatrix(io, RSneut, .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "sleptons"
  Write(io,*) mSlepton
  If (WithComments) Write(io,*) "  R_sl"
  Call WriteComplexMatrix(io, RSlepton, .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "u-squarks"
  Write(io,*) mSup
  If (WithComments) Write(io,*) "  R_su"
  Call WriteComplexMatrix(io, RSup, .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "d-squarks"
  Write(io,*) mSdown
  If (WithComments) Write(io,*) "  R_sd"
  Call WriteComplexMatrix(io, RSdown, .False.)
  If (WithComments) Write(io,*) " "

 Else
  If (WithComments) Then
   Write(io,*) "e-sneutrino mass   : ",mSneut(1)
   Write(io,*) "mu-sneutrino mass  : ",mSneut(2)
   Write(io,*) "tau-sneutrino mass : ",mSneut(3)
   Write(io,*) " "
  Else
   Write(io,*) mSneut
  End If 

  If (WithComments) Write(io,*) "selectron masses"
  Write(io,*) mSlepton(1:2)
  If (WithComments) Write(io,*) " R_e"
  Call WriteComplexMatrix(io, RSlepton(1:2,1:2), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "smuon masses"
  Write(io,*) mSlepton(3:4)
  If (WithComments) Write(io,*) " R_mu"
  Call WriteComplexMatrix(io, RSlepton(3:4,3:4), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "stau masses"
  Write(io,*) mSlepton(5:6)
  If (WithComments) Write(io,*) " R_tau"
  Call WriteComplexMatrix(io, RSlepton(5:6,5:6), .False.)
  If (WithComments) Write(io,*) " "
  
  If (WithComments) Write(io,*) "u-squark masses"
  Write(io,*) mSup(1:2)
  If (WithComments) Write(io,*) " R_u"
  Call WriteComplexMatrix(io, RSup(1:2,1:2), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "c-squark masses"
  Write(io,*) mSup(3:4)
  If (WithComments) Write(io,*) " R_c"
  Call WriteComplexMatrix(io, RSup(3:4,3:4), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "t-squark masses"
  Write(io,*) mSup(5:6)
  If (WithComments) Write(io,*) " R_t"
  Call WriteComplexMatrix(io, RSup(5:6,5:6), .False.)
  If (WithComments) Write(io,*) " "
  
  If (WithComments) Write(io,*) "d-squark masses"
  Write(io,*) mSdown(1:2)
  If (WithComments) Write(io,*) " R_d"
  Call WriteComplexMatrix(io, RSdown(1:2,1:2), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "s-squark masses"
  Write(io,*) mSdown(3:4)
  If (WithComments) Write(io,*) " R_s"
  Call WriteComplexMatrix(io, RSdown(3:4,3:4), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "b-squark masses"
  Write(io,*) mSdown(5:6)
  If (WithComments) Write(io,*) " R_b"
  Call WriteComplexMatrix(io, RSdown(5:6,5:6), .False.)
  If (WithComments) Write(io,*) " "
  
 End If 

 If (WithComments) Write(io,*) "m_A0, m_H+"
 Write(io,*) mP0(2), mSpm(2)
 If (WithComments) Write(io,*) " "
 If (WithComments) Write(io,*) "m_h0, m_H0"
 Write(io,*) mS0
 If (WithComments) Write(io,*) " R_S0"
 Call WriteRealMatrix(io, RS0, .False.)
 If (WithComments) Write(io,*) " "

 End Subroutine WriteMassesMSSM


 Subroutine WriteMassesNMSSM(io, WithComments, mGlu, PhaseGlu, mC, U, V, mN, N &
           &, mSneut, Rsneut, mSlepton, RSlepton, mSdown, RSdown, mSup, RSup  &
           &, mP0, RP0, mS0, RS0, mSpm, GenerationMixing)
 
 Implicit None
  Integer, Intent(in) :: io
  Logical, Intent(in) :: WithComments, GenerationMixing
  Real(dp), Intent(in) :: mP0(3), RP0(3,3), mS0(3), RS0(3,3), mSpm(2), mC(2)  &
    & , mN(5), mGlu, mSneut(3), mSlepton(6), mSdown(6), mSup(6)
 Complex(dp), Intent(in) :: PhaseGlu, U(2,2), V(2,2), N(5,5), Rsneut(3,3)     &
    & , RSlepton(6,6), RSdown(6,6), RSup(6,6)

 If (WithComments) Write(io,*) "Masses and mixing matrices"
 If (WithComments) Then
  If (Aimag(phaseglu).Eq.0._dp) Then
   Write(io,*) "Gluino : ",mglu, Real(phaseglu)
  Else
   Write(io,*) "Gluino : ",mglu, phaseglu
  End If
 Else
  If (Aimag(phaseglu).Eq.0._dp) Then
   Write(io,*) mglu, Real(phaseglu)
  Else
   Write(io,*) mglu, phaseglu
  End If
 End If
 If (WithComments) Write(io,*) " "

 If (WithComments) Write(io,*) "Charginos"
 Write(io,*) mC
 If (WithComments) Write(io,*) "  U"
 Call WriteComplexMatrix(io, U, .False.)
 If (WithComments) Write(io,*) "  V"
 Call WriteComplexMatrix(io, V, .False.)
 If (WithComments) Write(io,*) " "

 If (WithComments) Write(io,*) "Neutralinos"
 Write(io,*) mN
 If (WithComments) Write(io,*) "  N"
 Call WriteComplexMatrix(io, N, .False.)
 If (WithComments) Write(io,*) " "

 If (GenerationMixing) Then
  If (WithComments) Write(io,*) "sneutrinos"
  Write(io,*) mSneut
  If (WithComments) Write(io,*) "  R_sn"
  Call WriteComplexMatrix(io, RSneut, .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "sleptons"
  Write(io,*) mSlepton
  If (WithComments) Write(io,*) "  R_sl"
  Call WriteComplexMatrix(io, RSlepton, .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "u-squarks"
  Write(io,*) mSup
  If (WithComments) Write(io,*) "  R_su"
  Call WriteComplexMatrix(io, RSup, .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "d-squarks"
  Write(io,*) mSdown
  If (WithComments) Write(io,*) "  R_sd"
  Call WriteComplexMatrix(io, RSdown, .False.)
  If (WithComments) Write(io,*) " "

 Else
  If (WithComments) Then
   Write(io,*) "e-sneutrino mass   : ",mSneut(1)
   Write(io,*) "mu-sneutrino mass  : ",mSneut(2)
   Write(io,*) "tau-sneutrino mass : ",mSneut(3)
   Write(io,*) " "
  Else
   Write(io,*) mSneut
  End If 

  If (WithComments) Write(io,*) "selectron masses"
  Write(io,*) mSlepton(1:2)
  If (WithComments) Write(io,*) " R_e"
  Call WriteComplexMatrix(io, RSlepton(1:2,1:2), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "smuon masses"
  Write(io,*) mSlepton(3:4)
  If (WithComments) Write(io,*) " R_mu"
  Call WriteComplexMatrix(io, RSlepton(3:4,3:4), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "stau masses"
  Write(io,*) mSlepton(5:6)
  If (WithComments) Write(io,*) " R_tau"
  Call WriteComplexMatrix(io, RSlepton(5:6,5:6), .False.)
  If (WithComments) Write(io,*) " "
  
  If (WithComments) Write(io,*) "u-squark masses"
  Write(io,*) mSup(1:2)
  If (WithComments) Write(io,*) " R_u"
  Call WriteComplexMatrix(io, RSup(1:2,1:2), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "c-squark masses"
  Write(io,*) mSup(3:4)
  If (WithComments) Write(io,*) " R_c"
  Call WriteComplexMatrix(io, RSup(3:4,3:4), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "t-squark masses"
  Write(io,*) mSup(5:6)
  If (WithComments) Write(io,*) " R_t"
  Call WriteComplexMatrix(io, RSup(5:6,5:6), .False.)
  If (WithComments) Write(io,*) " "
  
  If (WithComments) Write(io,*) "d-squark masses"
  Write(io,*) mSdown(1:2)
  If (WithComments) Write(io,*) " R_d"
  Call WriteComplexMatrix(io, RSdown(1:2,1:2), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "s-squark masses"
  Write(io,*) mSdown(3:4)
  If (WithComments) Write(io,*) " R_s"
  Call WriteComplexMatrix(io, RSdown(3:4,3:4), .False.)
  If (WithComments) Write(io,*) " "

  If (WithComments) Write(io,*) "b-squark masses"
  Write(io,*) mSdown(5:6)
  If (WithComments) Write(io,*) " R_b"
  Call WriteComplexMatrix(io, RSdown(5:6,5:6), .False.)
  If (WithComments) Write(io,*) " "
  
 End If 

 If (WithComments) then
   Write(io,*) "m_H+", mSpm(2)
 else
   Write(io,*)  mSpm(2)
 end if
 If (WithComments) Write(io,*) " "
 If (WithComments) Write(io,*) "m_G0 m_P0_1, m_P0_2"
 Write(io,*) mP0
 If (WithComments) Write(io,*) " R_P0"
 Call WriteRealMatrix(io, RP03, .False.)
 If (WithComments) Write(io,*) " "
 If (WithComments) Write(io,*) "m_S0_1, m_S0_2, m_S0_3"
 Write(io,*) mS0
 If (WithComments) Write(io,*) " R_S0"
 Call WriteRealMatrix(io, RS03, .False.)
 If (WithComments) Write(io,*) " "

 End Subroutine WriteMassesNMSSM


 Subroutine WriteMatrixBlockC(io, len, mat_in, scale, name, com1, com2, sym, tr)
 Implicit None
  Integer, Intent(in) :: io, len
  Real(dp), Intent(in) :: scale
  Complex(dp), Intent(in) :: mat_in(len,len)
  Character(len=*), Intent(in) :: name, com1, com2
  Logical, Intent(in), Optional :: tr, sym

  Integer :: i1, i2
  Real(dp) :: mat(len,len)
  Logical :: trans, symmetric

  trans = .False.
  If (Present(tr)) trans = tr
  symmetric = .False.
  If (Present(sym)) symmetric = sym

  mat = Real(mat_in,dp)

  Write(io,106) "Block "//Trim(name)//" Q=",scale,"# "//Trim(com1)
  If (trans) mat = Transpose(mat)

  If (symmetric) Then
   Do i2=1,len
    Do i1=i2,len
     Write(io,105) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do

  Else

   Do i2=1,len
    Do i1=1,len
     Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do
  End If

  mat = Aimag(mat_in)
  If (Maxval(Abs(mat)).ne.0._dp) then
   Write(io,106) "Block IM"//Trim(name)//" Q=",scale,"# "//Trim(com1)
   If (trans) mat = Transpose(mat)

   If (symmetric) Then
    Do i2=1,len
     Do i1=i2,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do

   Else

    Do i2=1,len
     Do i1=1,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do
   End If
  End If

105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockC


 Subroutine WriteMatrixBlockC2(io, len, mat_in, name, com1, com2, sym, tr)
 Implicit None
  Integer, Intent(in) :: io, len
  Complex(dp), Intent(in) :: mat_in(len,len)
  Character(len=*), Intent(in) :: name, com1, com2
  Logical, Intent(in), Optional :: tr, sym

  Integer :: i1, i2
  Real(dp) :: mat(len,len)
  Logical :: trans, symmetric

  trans = .False.
  If (Present(tr)) trans = tr
  symmetric = .False.
  If (Present(sym)) symmetric = sym

  mat = Real(mat_in,dp)

  Write(io,106) "Block "//Trim(name)//" # "//Trim(com1)
  If (trans) mat = Transpose(mat)

  If (symmetric) Then
   Do i2=1,len
    Do i1=i2,len
     Write(io,105) i2,i1,mat(i2,i1) &
                 & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do

  Else

   Do i2=1,len
    Do i1=1,len
     Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Re["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
    End Do
   End Do
  End If

  mat = Aimag(mat_in)
  If (Maxval(Abs(mat)).ne.0._dp) then
   Write(io,106) "Block IM"//Trim(name)//" # imaginary parts of "//Trim(com1)
   If (trans) mat = Transpose(mat)

   If (symmetric) Then
    Do i2=1,len
     Do i1=i2,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do

   Else

    Do i2=1,len
     Do i1=1,len
      Write(io,105) i2,i1,mat(i2,i1) &
           & ,"# Im["//Trim(com2)//bu(i2)//","//bu(i1)//")]"
     End Do
    End Do
   End If
  End If

105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockC2


 Subroutine WriteMatrixBlockR(io, len, mat, scale, name, com1, com2, com3, tr)
 Implicit None
  Integer, Intent(in) :: io, len
  Real(dp), Intent(in) :: mat(len,len), scale
  Character(len=*), Intent(in) :: name, com1, com2, com3
  Logical, Intent(in), Optional :: tr

  Integer :: i1, i2
  Logical :: trans

  trans = .False.
  If (Present(tr)) trans = tr

  Write(io,106) "Block "//Trim(name)//" Q=",scale,"# "//Trim(com1)
  Do i2=1,len
   Do i1=1,len
    If (trans) Then
     Write(io,105) i2,i1,mat(i1,i2),"# "//Trim(com2)//bu(i2)//","//bu(i1)//Trim(com3)
    Else
     Write(io,105) i2,i1,mat(i2,i1),"# "//Trim(com2)//bu(i2)//","//bu(i1)//Trim(com3)
    End If

   End Do
  End Do

105 Format(1x,2i3,3x,1P,e16.8,3x,a)
106 Format(a,1P,e16.8,2x,a)

 End Subroutine WriteMatrixBlockR

 Subroutine WriteMSSMParameters(n, CKM, Yd, Yu, WithComments)
 !-----------------------------------------------------------------------
 ! writes out all MSSM parameters
 ! written by Werner Porod
 ! 15.11.01: writting MSSM part,
 !           WriteBilinear and WriteTrilinear should be used later for 
 !           R-parity violation
 !-----------------------------------------------------------------------
 Implicit None
  Complex(dp), Intent(in) :: CKM(3,3)
  Integer, Intent(in) :: n
  Logical, Intent(in) :: WithComments
  Real(dp), Intent(in) :: Yd(3), Yu(3)

  Logical :: OnlyDiagonal

  OnlyDiagonal = .Not.GenerationMixing

103 Format(1P,3e16.8)

  If (WithComments) Write(n,*) "       g'             g             g_3"
  Write(n,103) gauge
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "      Y_e            Y_mu          Y_tau"
  Write(n,903) Real(Y_l(1,1),dp),Real(Y_l(2,2),dp),Real(Y_l(3,3),dp)
!   Call WriteComplexMatrix(n, Y_l, OnlyDiagonal)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "      Y_u            Y_c            Y_t"
  Write(n,903) Yu
 903 Format(t2,3f16.5)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "      Y_d            Y_s            Y_b"
  Write(n,903) Yd
!   Call WriteComplexMatrix(n, Y_d, OnlyDiagonal)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) " CKM matrix"
  If (GenerationMixing) Call WriteComplexMatrix(n, CKM, OnlyDiagonal)
  If (WithComments) Write(n,*) " "


  If (WithComments) Write(n,*) "Gaugino mass parameters"
  If (Maxval(Abs(Aimag(Mi))).Eq.0._dp) Then
   Write(n,103) Real( Mi )
  Else
   Write(n,*) Mi
  End If
  If (WithComments) Write(n,*) " "

  
  If (WithComments) Write(n,*) "tan(beta), mu, B"
  If ((Aimag(mu).Eq.0._dp).And.(Aimag(B).Eq.0._dp)) Then
   Write(n,103) tanb, Real(mu), Real(B)
  Else
   Write(n,*) tanb, mu, B
  End If
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "Slepton mass parameters"
  If (WithComments) Write(n,*) "A_l"
  Call WriteComplexMatrix(n, A_l, OnlyDiagonal)
  If (WithComments) Write(n,*) "M2_E"
  Call WriteComplexMatrix(n, M2_E, OnlyDiagonal)
  If (WithComments) Write(n,*) "M2_L"
  Call WriteComplexMatrix(n, M2_L, OnlyDiagonal)
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "Squark mass parameters"
  If (GenerationMixing) Then
   If (WithComments) Write(n,*) "A_d"
   Call WriteComplexMatrix(n, Ad_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "A_u"
   Call WriteComplexMatrix(n, Au_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_D"
   Call WriteComplexMatrix(n, M2D_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_U"
   Call WriteComplexMatrix(n, M2U_sckm, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_Q"
   Call WriteComplexMatrix(n, M2Q_sckm, OnlyDiagonal)
  Else
   If (WithComments) Write(n,*) "A_d"
   Call WriteComplexMatrix(n, A_d, OnlyDiagonal)
   If (WithComments) Write(n,*) "A_u"
   Call WriteComplexMatrix(n, A_u, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_D"
   Call WriteComplexMatrix(n, M2_D, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_U"
   Call WriteComplexMatrix(n, M2_U, OnlyDiagonal)
   If (WithComments) Write(n,*) "M2_Q"
   Call WriteComplexMatrix(n, M2_Q, OnlyDiagonal)
  End If
  If (WithComments) Write(n,*) " "

  If (WithComments) Write(n,*) "Higgs mass parameters"
  Write(n,102) M2_H
  If (WithComments) Write(n,*) " "

102 Format(1P,2e16.8)

 End Subroutine WriteMSSMParameters


 Subroutine WriteRealMatrix(n, matrix, OnlyDiagonal)
 !---------------------------------------------------------------
 ! simplifies the writing of complex matrices in various places
 ! written by Werner Porod
 ! 15.11.01
 !---------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: n
  Real(dp), Intent(in) :: matrix(:,:)
  Logical, Intent(in) :: OnlyDiagonal

  Integer :: i1, i2, dim_matrix
  Real(dp), Allocatable :: vec(:)

  dim_matrix = Size( matrix, dim=1)

  If (OnlyDiagonal) Then
   Allocate( vec(dim_matrix) )
   Do i1=1,dim_matrix
    vec(i1) = matrix(i1,i1)
   End Do
   Write(n,*) vec
   Deallocate( vec )
    
  Else
    If (dim_matrix.Le.4) Then
     Do i1=1,dim_matrix
      If (dim_matrix.Eq.2) Write(n,902) matrix(i1,:)
      If (dim_matrix.Eq.3) Write(n,903) matrix(i1,:)
      If (dim_matrix.Eq.4) Write(n,904) matrix(i1,:)
     End Do
    Else
     Do i1=1,dim_matrix
      Do i2=1,dim_matrix
       Write(n,901) i1,i2,matrix(i1,i2)
      End Do
     End Do
    End If
  End If

 901 Format(t2,2i4,f16.5)
 902 Format(t2,2f16.5)
 903 Format(t2,3f16.5)
 904 Format(t2,4f16.5)

 End Subroutine WriteRealMatrix


End Module InputOutput

