Module Chargino3Decays

! load modules
Use Control
Use ThreeBodyPhaseSpace
! load modules

Contains


 Subroutine CharginoThreeBodyDecays(n_in, mC, mZ, gZ, L_nu, R_nu          &
    & , mf_l, L_e, R_e, mf_u, L_u, R_u, mf_d, L_d, R_d, cpl_CCZ_L, cpl_CCZ_R &
    & , mN, mW, gW, cpl_NuLW, cpl_UDW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L       &
    & , cpl_CNW_R, mSpm, gSpm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L        &
    & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gS0, cpl_NNS0_L         &
    & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R           &
    & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gP0, cpl_NNP0_L &
    & , cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L, cpl_DDP0_R           &
    & , cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mUSquark, gUSquark   &
    & , cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu, cpl_UGSu_L     &
    & , cpl_UGSu_R, mDSquark, gDSquark, cpl_DNSd_L, cpl_DNSd_R, cpl_CUSd_L   &
    & , cpl_CUSd_R, cpl_DGSd_L, cpl_DGSd_R, mSneutrino, gSneutrino           &
    & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gSlepton &
    & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                     &
    & , GenerationMixing, k_neut, epsI, deltaM, Check_Real_States            &
    & , gMajaron, gCnunu, gCll, gCdd, gCuu, gNln, gNDU, gCgdu &
    & , gCCCC, GCCNN, gT, gP, BR)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a chargino into fermions 
 ! in case of spontaneous R-parity
 ! violation also the decays into majoron are included. The idea is that the
 ! couplings are calculated in another subroutine and that this routine
 ! can be used in a model independent way.
 ! input:
 !    n_in ....... index of the decaying chargino, in case n_in < 0,
 !                 the decay widhts of all charginos will be calculated
 !    mN ......... chargino masses
 !    mC ......... chargino masses
 !
 ! output:
 ! 2-body decay modes
 !  gMajaron(i) ......... chargino_i + majoron J
 !  gPhoton(i) .......... chargino_i + photon
 ! 3-body decay modes
 !  gNll(i,j,k) ......... Neutralino_i + lepton_j + lepton_k (j,k=1 e
 !                                                            j,k=2 mu
 !                                                            j,k=3 tau)
 !  gNqq(i,j) ......... Neutralino_i + quark_j + quark_j (j,k=1 u
 !                                                        j,k=2 d
 !                                                        j,k=3 c
 !                                                        j,k=4 s
 !                                                        j,k=5 t
 !                                                        j,k=6 b)
 !  gCln(i,j) ........... Chargino_i + l_j nu (j=1 -> e nu
 !                                              j=2 -> mu nu
 !                                              j=3 -> tau nu )
 !  gCDU(i,j) .......... Chargino_i + q_j q'_j (j=1 -> u d
 !                                                j=2 -> c s
 !                                                j=3 -> t b )
 !  gP(i) ............... the above described partial widths in one
 !                        array, useful for further processing
 !  gT .................. total width
 !  BR(i) ............... branching ratios, optional
 ! written by Werner Porod, 16.05.2001
 !  - taking the code from the Routine NeutralinoDecays.f as basis 
 !  13.09.03: including the possiblity to check for real intermediates states
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: n_in, k_neut
  Real(dp), Intent(in) :: mN(:), mZ, gZ, L_nu, R_nu, mf_l(3), L_e, R_e      &
     &, mf_u(3), L_u, R_u, mf_d(3), L_d, R_d, mC(:), mW, gW                 &
     &, mSpm(:), gSpm(:), mS0(:), gS0(:), mP0(:), gP0(:), mUSquark(6)       &
     &, gUSquark(6), mGlu, mDSquark(6), gDSquark(6), mSneutrino(:)          &
     &, gSneutrino(:), mSlepton(:), gSlepton(:), epsI, deltaM
  Real(dp), Intent(inout) :: gMajaron(:,:), gT(:), gCnunu(:,:,:,:) &
     & , gCll(:,:,:,:), gCdd(:,:,:,:), gCuu(:,:,:,:), gNln(:,:,:,:)           &
     & , gNDU(:,:,:,:), gCCNN(:,:,:,:), gCCCC(:,:,:,:), gCgdu(:,:,:)
  Real(dp), Intent(inout), Optional :: BR(:,:), gP(:,:)
  Complex(dp), Intent(in) :: cpl_NNZ_L(:,:), cpl_NNZ_R(:,:), cpl_CCZ_L(:,:)  &
     & , cpl_CCZ_R(:,:), cpl_CNW_L(:,:), cpl_CNW_R(:,:), cpl_UDW(:,:)        &
     & , cpl_SmpCN_L(:,:,:), cpl_SmpCN_R(:,:,:), cpl_SmpLNu_L(:,:,:)         &
     & , cpl_SmpLNu_R(:,:,:), cpl_SmpDU_L(:,:,:), cpl_SmpDU_R(:,:,:)         &
     & , cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:), cpl_CCS0_L(:,:,:)             &
     & , cpl_CCS0_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_R(:,:,:)             &
     & , cpl_LLS0_L(:,:,:), cpl_LLS0_R(:,:,:), cpl_UUS0_L(:,:,:)             &
     & , cpl_UUS0_R(:,:,:), cpl_NNP0_L(:,:,:), cpl_NNP0_R(:,:,:)             &
     & , cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:), cpl_DDP0_L(:,:,:)             &
     & , cpl_DDP0_R(:,:,:), cpl_LLP0_L(:,:,:), cpl_LLP0_R(:,:,:)             &
     & , cpl_UUP0_L(:,:,:), cpl_UUP0_R(:,:,:), cpl_UNSu_L(:,:,:)             &
     & , cpl_UNSu_R(:,:,:), cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)             &
     & , cpl_UGSu_L(:,:), cpl_UGSu_R(:,:), cpl_DNSd_L(:,:,:)                 &
     & , cpl_DNSd_R(:,:,:), cpl_CUSd_L(:,:,:), cpl_CUSd_R(:,:,:)             &
     & , cpl_DGSd_L(:,:), cpl_DGSd_R(:,:), cpl_NuNSn_L(:,:,:)                &
     & , cpl_NuNSn_R(:,:,:), cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:)            &
     & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:), cpl_CNuSl_L(:,:,:)            &
     & , cpl_CNuSl_R(:,:,:), cpl_NuLW(:,:)
  Logical, Intent(in) :: GenerationMixing, Check_Real_States

  Integer :: n_neut, n_char, n_S0, n_P0, n_Spm, i_start, i_end, i_run       &
     &, i1, i2, i3, i4, n_length, n_Z4, n_S04, n_Sf4, n_ZS04, n_ZSf8 &
     & , n_S0P04, n_S0Sf8, n_CSf4, n_Sf8, nf, n_W4, n_WSpm, n_WSf8, n_Z8    &
     & , n_ZS08, n_S0P08, n_W8, n_ZW8, n_WSpm8, n_WSpm4, n_Sl, n_Sd, n_Su   &
     & , n_Sn
  Real(dp) :: factor(3), m_nu(3),  gCff(3,3), gNffp(3,3), gCCC, gCNN, diffM &
     & , gW_in, gZ_in
  Real(dp), Allocatable :: IntegralsZ4(:,:), IntegralsS04(:,:)      &
     & , IntegralsSf4(:,:), IntegralsW4(:,:), gS0_in(:), gP0_in(:)  &
     & , gSpm_in(:), g_Su(:), g_Sd(:), g_Sl(:), g_Sn(:)
  Complex(dp), Allocatable :: IntegralsZS04(:,:), IntegralsZSf8(:,:)     &
      & , IntegralsS0P04(:,:), IntegralsS0Sf8(:,:), IntegralsCSf4(:,:)   &
      & , IntegralsSf8(:,:), IntegralsWSpm4(:,:), IntegralsWSf8(:,:)     &
      & , IntegralsZ8(:,:), IntegralsZS08(:,:), IntegralsS0P08(:,:)      &
      & , IntegralsW8(:,:), IntegralsZW8(:,:), IntegralsWSpm8(:,:)
  logical :: check

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CharginoThreeBodyDecays'

  !--------------------
  ! checking for model 
  !--------------------
  n_neut = Size( mN )
  n_char = Size( mC )
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )
  n_Spm = Size( mSpm )
  n_Su = Size( gUSquark )
  n_Sd = Size( gDSquark )
  n_Sl = Size( gSlepton )
  n_Sn = Size( gSneutrino )

  allocate( gS0_in(n_S0) )
  allocate( gP0_in(n_P0) )
  allocate( gSpm_in(n_Spm) )
  allocate( g_Su(n_Su) )
  allocate( g_Sd(n_Sd) )
  allocate( g_Sn(n_Sn) )
  allocate( g_Sl(n_Sl) )

   Allocate( IntegralsZ4(200,8) )
   Allocate( IntegralsZ8(200,12) )
   Allocate( IntegralsW4(200,8) )
   Allocate( IntegralsW8(200,12) )
   Allocate( IntegralsZW8(200,12) )
   Allocate( IntegralsS04(2500,10) )
   Allocate( IntegralsSf4(2500,10) )
   Allocate( IntegralsZS04(2500,12) )
   Allocate( IntegralsZS08(2500,16) )
   Allocate( IntegralsWSpm4(2500,12) )
   Allocate( IntegralsWSpm8(2500,16) )
   Allocate( IntegralsZSf8(2500,16) )
   Allocate( IntegralsWSf8(2500,16) )
   Allocate( IntegralsS0P04(2500,12) )
   Allocate( IntegralsS0P08(2500,16) )
   Allocate( IntegralsS0Sf8(2500,16) )
   Allocate( IntegralsCSf4(2500,12) )
   Allocate( IntegralsSf8(2500,16) )

  If (n_in.Lt.0) Then
   i_start = 1
   i_end = n_char

  Else If ( (n_in.Ge.1).And.(n_in.Le.n_char) ) Then 
   i_start = n_in 
   i_end = n_in

  Else
   If (ErrorLevel.Ge.-1) Then
    Write (ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write (ErrCan,*) 'Value of n_in out of range, (n_in,n_char) = ',n_in,n_char
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
   Iname = Iname - 1
   Return
  End If

  m_nu = 0._dp

   If (Check_Real_States) Then
    gZ_in = 0._dp
    gW_in = 0._dp
    gS0_in = 0._dp
    gP0_in = 0._dp
    gSpm_in = 0._dp
    g_Su = 0._dp
    g_Sd = 0._dp 
    g_Sl = 0._dp 
    g_Sn =  0._dp
   Else
    gZ_in = gZ
    gW_in = gW
    gS0_in = gS0
    gP0_in = gP0
    gSpm_in = gSpm
    g_Su = gUSquark
    g_Sd = gDSquark
    g_Sl = gSlepton
    g_Sn = gSneutrino
   End If
   check = Check_Real_States

  Do i_run = i_start, i_end
   !-----------------
   ! intitialisation
   !-----------------
   gT(i_run) = 0._dp
   gMajaron(i_run,:) = 0._dp
   gCnunu(i_run,:,:,:) = 0._dp
   gCdd(i_run,:,:,:) = 0._dp
   gCuu(i_run,:,:,:) = 0._dp
   gCgdu(i_run,:,:) = 0._dp
   gCll(i_run,:,:,:) = 0._dp
   gCCNN(i_run,:,:,:) = 0._dp
   gCCCC(i_run,:,:,:) = 0._dp
   gNln(i_run,:,:,:) = 0._dp
   gNDU(i_run,:,:,:) = 0._dp

   IntegralsZ4 = 0._dp
   IntegralsW4 = 0._dp
   IntegralsS04 = 0._dp
   IntegralsSf4 = 0._dp
   IntegralsZS04 = 0._dp
   IntegralsWSpm4 = 0._dp
   IntegralsZ8 = 0._dp
   IntegralsZS08 = 0._dp
   IntegralsZSf8 = 0._dp
   IntegralsWSf8 = 0._dp
   IntegralsS0P04 = 0._dp
   IntegralsS0P08 = 0._dp
   IntegralsS0Sf8 = 0._dp
   IntegralsCSf4 = 0._dp
   IntegralsSf8 = 0._dp
   n_Z4 = 0
   n_W4 = 0
   n_S04 = 0
   n_Sf4 = 0
   n_ZS04 = 0
   n_WSpm = 0
   n_WSpm8 = 0
   n_Z8 = 0
   n_W8 = 0
   n_ZW8 = 0
   n_ZSf8 = 0
   n_ZS08 = 0
   n_WSpm4 = 0
   n_WSf8 = 0
   n_S0P04 = 0
   n_S0P08 = 0
   n_S0Sf8 = 0
   n_CSf4 = 0
   n_Sf8 = 0
   factor(1) = oo512pi3 / Abs(mC(i_run))**3   ! for 3-body decays
   factor(2) = 3._dp * factor(1)            ! including color factor
   factor(3) = 4._dp * factor(1)            ! for decays into gluinos
   !--------------------------------------
   ! decays into a chargino + 2 fermions
   !--------------------------------------
   Do i1 = 1, i_run - 1
    gCuu(i_run, i1,:,:) = 0._dp
    gCdd(i_run, i1,:,:) = 0._dp
    gCll(i_run, i1,:,:) = 0._dp
    gCnunu(i_run,i1,:,:) = 0._dp

    Call ChimToChimff(i_run, i1, ' u u ', mC, mZ, gZ_in, Cpl_CCZ_L, Cpl_CCZ_R &
     & , 3, mf_u, L_u, R_u, IntegralsZ4, n_Z4, mS0, gS0_in, cpl_CCS0_L        &
     & , cpl_CCS0_R, cpl_UUS0_L, cpl_UUS0_R, IntegralsS04, n_S04, mP0, gP0_in &
     & , cpl_CCP0_L, cpl_CCP0_R, cpl_UUP0_L, cpl_UUP0_R, mDSquark, g_Sd       &
     & , cpl_CUSd_L, cpl_CUSd_R, IntegralsSf4, n_Sf4, IntegralsZS04, n_ZS04   &
     & , IntegralsZSf8, n_ZSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8       &
     & , n_S0Sf8, IntegralsCSf4, n_CSf4, deltaM, epsI, GenerationMixing       &
     & , check, factor(2), gCff)
    gCuu(i_run, i1,:,:) = gCff
    gT(i_run) = gT(i_run) + Sum( gCff )

    Call ChimToChimff(i_run, i1, ' d d ', mC, mZ, gZ_in, Cpl_CCZ_L, Cpl_CCZ_R &
     & , 3, mf_d, L_d, R_d, IntegralsZ4, n_Z4, mS0, gS0_in, cpl_CCS0_L        &
     & , cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R, IntegralsS04, n_S04, mP0, gP0_in &
     & , cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L, cpl_DDP0_R, mUSquark, g_Su       &
     & , cpl_CDSu_L, cpl_CDSu_R, IntegralsSf4, n_Sf4, IntegralsZS04, n_ZS04   &
     & , IntegralsZSf8, n_ZSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8       &
     & , n_S0Sf8, IntegralsCSf4, n_CSf4, deltaM, epsI, GenerationMixing       &
     & , check, factor(2), gCff)
    gCdd(i_run, i1,:,:) = gCff
    gT(i_run) = gT(i_run) + Sum( gCff )

    If (n_char.Lt.5) Then
     nf = 5 - n_char

     Call ChimToChimff(i_run, i1, ' l l ', mC, mZ, gZ_in, Cpl_CCZ_L           &
      & , Cpl_CCZ_R, nf, mf_l, L_e, R_e, IntegralsZ4, n_Z4, mS0, gS0_in       &
      & , cpl_CCS0_L, cpl_CCS0_R, cpl_LLS0_L, cpl_LLS0_R, IntegralsS04, n_S04 &
      & , mP0, gP0_in, cpl_CCP0_L, cpl_CCP0_R, cpl_LLP0_L, cpl_LLP0_R         &
      & , mSneutrino, g_Sn, cpl_CLSn_L, cpl_CLSn_R, IntegralsSf4, n_Sf4       &
      & , IntegralsZS04, n_ZS04, IntegralsZSf8, n_ZSf8, IntegralsS0P04        &
      & , n_S0P04, IntegralsS0Sf8, n_S0Sf8, IntegralsCSf4, n_CSf4             &
      & , deltaM, epsI, GenerationMixing, check, factor(1), gCff)
     gT(i_run) = gT(i_run) + Sum( gCff )
     gCll(i_run, i1,:,:) = gCff

     Call  ChimToChimNuNu(i_run, i1, mC, mZ, gZ_in, Cpl_CCZ_L, Cpl_CCZ_R, nf  &
      & , L_nu, R_nu, IntegralsZ4, n_Z4, mSlepton, g_Sl, cpl_CNuSl_L          &
      & , cpl_CNuSl_R, IntegralsSf4, n_Sf4, IntegralsZSf8, n_ZSf8             &
      & , IntegralsCSf4, n_CSf4, deltaM, epsI, GenerationMixing, check        &
      & , factor(1), gCff)
     gCnunu(i_run, i1,:,:) = gCff
     gT(i_run) = gT(i_run) +  Sum(gCff)

    End If

   End Do

   !--------------------------------------
   ! decay into neutralinos + 2 SM fermions
   !--------------------------------------
   Do i1=1,n_neut
    gNDU(i_run,i1,:,:) = 0._dp
    gNln(i_run,i1,:,:) = 0._dp
    If (Abs(mC(i_run)).Gt.Abs(mN(i1))) Then
     Call ChimToChi0ffp(i_run, i1, mC, mN, 3, mf_d, mf_u, mW, gW_in          &
      & , Cpl_CNW_L, Cpl_CNW_R, Cpl_UDW, mSpm, gSpm_in, Cpl_SmpCN_L          &
      & , Cpl_SmpCN_R, Cpl_SmpDU_L, Cpl_SmpDU_R, mDSquark, g_Sd, cpl_DNSd_L  &
      & , cpl_DNSD_R, cpl_CUSd_L, cpl_CUSd_R, mUSquark, g_Su, cpl_UNSu_L     &
      & , cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, IntegralsW4, n_W4              &
      & , IntegralsS04, n_S04, IntegralsSf4, n_Sf4, IntegralsWSpm4, n_WSpm   &
      & , IntegralsWSf8, n_WSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8     &
      & , n_S0Sf8, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8, deltaM, epsI  &
      & , GenerationMixing, check, factor(2), gNffp)
     gNDU(i_run,i1,:,:) = gNffp
     gT(i_run) = gT(i_run) + Sum( gNffp )

     nf = 5 - n_char
     If (nf.Gt.0) Then
      Call ChimToChi0ffp(i_run, i1, mC, mN, nf, mf_l, m_nu, mW, gW_in        &
       & , Cpl_CNW_L, Cpl_CNW_R, Cpl_NuLW, mSpm, gSpm_in, Cpl_SmpCN_L        &
       & , Cpl_SmpCN_R, Cpl_SmpLNu_L, Cpl_SmpLNu_R, mSlepton, g_Sl           &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, mSneutrino      &
       & , g_Sn, cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R            &
       & , IntegralsW4, n_W4, IntegralsS04, n_S04, IntegralsSf4, n_Sf4       &
       & , IntegralsWSpm4, n_WSpm, IntegralsWSf8, n_WSf8, IntegralsS0P04     &
       & , n_S0P04, IntegralsS0Sf8, n_S0Sf8, IntegralsCSf4, n_CSf4           &
       & , IntegralsSf8, n_Sf8, deltaM, epsI, GenerationMixing, check        &
       & , factor(1), gNffp)

      gNln(i_run,i1,:,:) = gNffp
      gT(i_run) = gT(i_run) + Sum( gNffp )
     End If
    End If
   End Do

   !--------------------------------------
   ! decay into a gluino + 2 quarks
   !--------------------------------------
   gCGDU(i_run,:,:) = 0._dp
   If (Abs(mC(i_run)).Gt.Abs(mglu)) Then
    Call ChimToGffp(i_run, mC, mGlu, mf_d, mf_u, mDSquark, g_Sd       &
    & , cpl_DGSd_L, cpl_DGSD_R, cpl_CUSd_L, cpl_CUSd_R, mUSquark, g_Su &
    & , cpl_UGSu_L, cpl_UGSu_R, cpl_CDSu_L, cpl_CDSu_R                     &
    & , IntegralsSf4, n_Sf4, IntegralsCSf4, n_CSf4, IntegralsSf8, n_Sf8 &
    & , deltaM, epsI, GenerationMixing, check, factor(3), gNffp)
     gCGDU(i_run,:,:) = gNffp
     gT(i_run) = gT(i_run) + Sum( gNffp )
   End If

   If (k_neut.Ne.0) Then
   !-------------------------------
   ! decay into 3 charginos
   !-------------------------------  
   Do i1=1,i_run-1
    Do i2=1,i1
     Do i3=1,i2
      diffM = Abs(mC(i_run)) - Abs(mC(i1)) - Abs(mC(i2)) - Abs(mC(i3))
      If (diffM.Gt.0._dp) Then
       Call ChimToChimChipChim(i_run, i1, i2, i3, mC, mZ, gZ_in, cpl_CCZ_L   &
        & , cpl_CCZ_R, mS0, gS0_in, cpl_CCS0_L, cpl_CCS0_R, mP0, gP0_in      &
        & , cpl_CCP0_L, cpl_CCP0_R, IntegralsZ4, n_Z4, IntegralsS04, n_S04   &
        & , IntegralsZ8, n_Z8, IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08  &
        & , IntegralsS0P04, n_S0P04, IntegralsS0P08, n_S0P08, deltaM, epsI   &
        & , check, factor(1), gCCC)
       gCCCC(i_run,i1,i2,i3) = gCCC
       gT(i_run) = gT(i_run) + gCCC
      End If
     End Do
    End Do
   End Do
   !-------------------------------------
   ! decay into chargino + 2 neutralinos
   !-------------------------------------
   Do i1=1,i_run-1
    Do i2=1,n_neut
     Do i3=1,i2
      diffM = Abs(mC(i_run)) - Abs(mC(i1)) - Abs(mN(i2)) - Abs(mN(i3))
      If (diffM.Gt.0._dp) Then
       Call ChimToChim2Chi0(i_run, i1, i2, i3, mC, mN, mZ, gZ_in, cpl_CCZ_L  &
        & , cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, mW, gW_in, cpl_CNW_L, cpl_CNW_R &
        & , mSpm, gSpm_in, cpl_SmpCN_L, cpl_SmpCN_R, mS0, gS0_in, cpl_CCS0_L &
        & , cpl_CCS0_R, cpl_NNS0_L, cpl_NNS0_R, mP0, gP0_in, cpl_CCP0_L      &
        & , cpl_CCP0_R, cpl_NNP0_L, cpl_NNP0_R, IntegralsZ4, n_Z4            &
        & , IntegralsW4, n_W4, IntegralsS04, n_S04, IntegralsZW8, n_ZW8      &
        & , IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08, IntegralsW8, n_W8  &
        & , IntegralsWSpm4, n_WSpm4, IntegralsWSpm8, n_WSpm8, IntegralsS0P04 &
        & , n_S0P04, IntegralsS0P08, n_S0P08, deltaM, epsI, check, factor(1) &
        & , gCNN)
       gCCNN(i_run,i1,i2,i3) = gCNN
       gT(i_run) = gT(i_run) + gCNN
      End If
     End Do
    End Do
   End Do
  End If ! k_neut

  End Do ! i_run

  If (Present(gP)) Then
   Do i1=i_start,i_end
    gP(i1,:) = 0._dp
    If (GenerationMixing) Then
     n_length = 1
     Do i2=1,n_Neut
      Do i3=1,3
       Do i4=1,3  
        gP(i1,n_length) = gNDU(i1,i2,i3,i4)
        n_length = n_length + 1
       End Do
      End Do
      Do i3=1,5 - n_char
       gP(i1,n_length) = Sum(gNln(i1,i2,i3,:))
       n_length = n_length + 1
      End Do
     End Do

     Do i3=1,3
      Do i4=1,3  
       gP(i1,n_length) = gCGDU(i1,i3,i4)
       n_length = n_length + 1
      End Do
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       Do i4=1,3  
        gP(i1,n_length) = gCuu(i1,i2,i3,i4)
        n_length = n_length + 1
       End Do
      End Do
      Do i3=1,3
       Do i4=1,3  
        gP(i1,n_length) = gCdd(i1,i2,i3,i4)
        n_length = n_length + 1
       End Do
      End Do
      Do i3=1,5-n_char
       Do i4=1,5-n_char
        gP(i1,n_length) = gCll(i1,i2,i3,i4)
        n_length = n_length + 1
       End Do
      End Do
      if (n_char.lt.5) then
       gP(i1,n_length+9) = sum(gCnunu(i1,i2,:,:))
       n_length = n_length + 1
      end if

     End Do

    Else ! .not.GenerationMixing

     n_length = 1
     Do i2=1,n_Neut
      Do i3=1,3
       gP(i1,n_length) = gNDU(i1,i2,i3,i3)
       n_length = n_length + 1
      End Do

      Do i3=1,5-n_char
       gP(i1,n_length) = gNln(i1,i2,i3,i3)
       n_length = n_length + 1
      End Do
     End Do

     Do i3=1,3
      gP(i1,n_length) = gCGDU(i1,i3,i3)
      n_length = n_length + 1
     End Do

     Do i2=1,i1-1
      Do i3=1,3
       gP(i1,n_length) = gCuu(i1,i2,i3,i3)
       n_length = n_length + 1
      End Do
      Do i3=1,3
       gP(i1,n_length) = gCdd(i1,i2,i3,i3)
       n_length = n_length + 1
      End Do
      Do i3=1,5-n_char
       gP(i1,n_length) = gCll(i1,i2,i3,i3)
       n_length = n_length + 1
      End Do
      if (n_char.lt.5) then
       gP(i1,n_length+9) = sum(gCnunu(i1,i2,:,:))
       n_length = n_length + 1
      end if
     End Do

    End If

    If (Present(BR).And.(gT(i1).Gt.0._dp)) Then
     BR(i1,:) = gP(i1,:) / gT(i1)
    Else If (Present(BR)) Then
     BR(i1,:) = 0._dp
    End If
   End Do
  End If

  Deallocate( IntegralsZ4, IntegralsZ8, IntegralsW4, IntegralsW8      &
    & , IntegralsZW8, IntegralsS04, IntegralsSf4, IntegralsZS04       &
    & , IntegralsZS08, IntegralsWSpm4, IntegralsWSpm8, IntegralsZSf8  &
    & , IntegralsWSf8, IntegralsS0P04, IntegralsS0P08, IntegralsS0Sf8 &
    & , IntegralsCSf4, IntegralsSf8 )
  Deallocate( gS0_in, gP0_in, gSpm_in, g_Su, g_Sd, g_Sn, g_Sl)

  Iname = Iname - 1

 End Subroutine CharginoThreeBodyDecays

 Subroutine ChimToChimff(i_in, i_out, state, mC, mZ, gZ, Cpl_CCZ_L, Cpl_CCZ_R &
    & , n_f, mf, L_f, R_f, IntegralsZ4, n_Z4, mS0, gS0, cpl_CCS0_L, cpl_CCS0_R&
    & , cpl_FFS0_L, cpl_FFS0_R, IntegralsS04, n_S04, mP0, gP0, cpl_CCP0_L     &
    & , cpl_CCP0_R, cpl_FFP0_L, cpl_FFP0_R, mSfp, gSfp, cpl_CFSfp_L           &
    & , cpl_CFSfp_R, IntegralsSf4, n_Sf4, IntegralsZS04, n_ZS04               &
    & , IntegralsZSf8, n_ZSf8, IntegralsS0P04, n_S0P04, IntegralsS0Sf8        &
    & , n_S0Sf8, IntegralsCSf4, n_CSf4, deltaM, epsI                          &
    & , GenerationMixing, check, fac, gCff, WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a neutralino to another Neutralino + fermion pair
 ! Written by Werner Porod, 04.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !--------------------------------------------------------------------------
 Implicit None
  Character(len=5), Intent(in) :: state 
  Integer, Intent(in) :: i_in, i_out, n_f
  Integer, Intent(inout) :: n_Z4, n_S04, n_Sf4, n_ZS04, n_ZSf8, n_S0P04 &
     & , n_S0Sf8, n_CSf4
  Logical, Intent(in) :: check

  Real(dp), Intent(in) :: mC(:), mZ, gZ, mf(:), L_f, R_f, mS0(:), gS0(:) &
      & , mP0(:), gP0(:), mSfp(:), gSfp(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsS04(:,:)         &
      & , IntegralsSf4(:,:)
  Complex(dp), Intent(in) :: Cpl_CCZ_L(:,:), Cpl_CCZ_R(:,:)             &
      & , cpl_CCS0_L(:,:,:), cpl_CCS0_R(:,:,:), cpl_FFS0_L(:,:,:)       &
      & , cpl_FFS0_R(:,:,:), cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:)       &
      & , cpl_FFP0_L(:,:,:), cpl_FFP0_R(:,:,:), cpl_CFSfp_L(:,:,:)       &
      & , cpl_CFSfp_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsZS04(:,:), IntegralsZSf8(:,:)  &
      & , IntegralsS0P04(:,:), IntegralsS0Sf8(:,:), IntegralsCSf4(:,:)
  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(out) :: gCff(:,:)

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution

  Integer :: Isum, i1, i2, i3, i4, n_S0, n_P0, i_gen, n_Sfp
  Real(dp) :: Boson2(2), mass(4), resR, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8)

  Real(dp), Allocatable :: gCffSum(:,:,:)
  Character(len=20), Allocatable :: Contribution(:,:,:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChimToChimff'

  mass(1) = mC(i_in)
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )
  n_Sfp = Size( mSfp )
  Isum = (1 + n_S0 + n_P0 + n_Sfp)**2
  Allocate( gCffSum(n_f, n_f, Isum) )
  Allocate( Contribution(n_f, n_f, Isum) )
   
  gCffSum = 0._dp
  Contribution = ' '

  !-----
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ
  coup1(1) = Cpl_CCZ_L(i_out,i_in)
  coup1(2) = Cpl_CCZ_R(i_out,i_in)
  coup1(3) = L_f
  coup1(4) = R_f
  Do i1=1,n_f
   mass(2) = mC(i_out)
   mass(3) = -mf(i1)
   mass(4) = mf(i1)
   Contribution(i1,i1,1) = 'Z f_'//Bu(i1)//' f_'//Bu(i1)
   Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, gCffSum(i1,i1,Isum), check )
  End Do

  !-------------------
  ! S_0 S_0, diagonal
  !-------------------
  Do i1=1,n_S0
   Isum = Isum + 1
   Boson2(1) = mS0(i1)
   Boson2(2) = gS0(i1)
   coup1(1) = cpl_CCS0_L(i_out,i_in,i1) 
   coup1(2) = cpl_CCS0_R(i_out,i_in,i1)
   Do i2=1,n_f
    coup1(3) = cpl_FFS0_L(i2,i2,i1)
    coup1(4) = cpl_FFS0_R(i2,i2,i1)
    mass(2) = mC(i_out)
    mass(3) = - mf(i2)
    mass(4) = mf(i2)
    Contribution(i2,i2,Isum) = 'S0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsS04, n_S04, gCffSum(i2, i2, Isum), check )
   End Do
  End Do

  !-------------------
  ! P_0 P_0, diagonal
  !-------------------
  Do i1=1,n_P0
   Isum = Isum + 1
   Boson2(1) = mP0(i1)
   Boson2(2) = gP0(i1)
   coup1(1) = cpl_CCP0_L(i_out,i_in,i1) 
   coup1(2) = cpl_CCP0_R(i_out,i_in,i1)
   Do i2=1,n_f
    coup1(3) = cpl_FFP0_L(i2,i2,i1)
    coup1(4) = cpl_FFP0_R(i2,i2,i1)
    mass(2) = mC(i_out)
    mass(3) = - mf(i2)
    mass(4) = mf(i2)
    Contribution(i2,i2,Isum) = 'P0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsS04, n_S04, gCffSum(i2, i2, Isum), check )
   End Do
  End Do

  !-----------------------------
  ! sfermion sfermion, diagonal
  !-----------------------------
  If (GenerationMixing) Then
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    Do i1=1,n_f
     Do i3=1,n_f
      coup1(1) = cpl_CFSfp_L(i_in,i1,i2)
      coup1(2) = cpl_CFSfp_R(i_in,i1,i2)
      coup1(3) = Conjg(cpl_CFSfp_R(i_out,i3,i2))
      coup1(4) = Conjg(cpl_CFSfp_L(i_out,i3,i2))
      mass(2) = mf(i1)
      mass(3) = - mf(i3)
      mass(4) = mC(i_out)        
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gCffSum(i1,i3,Isum) = resR
      Contribution(i1,i3,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i3)
     End Do ! i3 fermions
    End Do  ! i1 fermions
   End Do   ! i2 sfermions

  Else If (n_Sfp.Gt.3) Then ! Sneutrinos are special
   Do i2=1,2*n_f  ! fermion generation
    i1 = (i2+1)/2
     Isum = Isum + 1
     Boson2(1) = mSfp(i2)
     Boson2(2) = gSfp(i2)
     coup1(1) = cpl_CFSfp_L(i_in,i1,i2)
     coup1(2) = cpl_CFSfp_R(i_in,i1,i2)
     coup1(3) = Conjg(cpl_CFSfp_R(i_out,i1,i2))
     coup1(4) = Conjg(cpl_CFSfp_L(i_out,i1,i2))
     mass(2) = mf(i1)
     mass(3) = - mf(i1)
     mass(4) = mC(i_out)        
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSf4, n_Sf4, resR, check)
     gCffSum(i1,i1,Isum) = resR
     Contribution(i1,i1,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
   End Do   ! i1 u-quarks

  Else  ! Sneutrinos are special
   Do i2=1,n_f  ! fermion generation
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    coup1(1) = cpl_CFSfp_L(i_in,i2,i2)
    coup1(2) = cpl_CFSfp_R(i_in,i2,i2)
    coup1(3) = Conjg(cpl_CFSfp_R(i_out,i2,i2))
    coup1(4) = Conjg(cpl_CFSfp_L(i_out,i2,i2))
    mass(2) = mf(i2)
    mass(3) = - mf(i2)
    mass(4) = mC(i_out)        
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsSf4, n_Sf4, resR, check)
    gCffSum(i2,i2,Isum) = resR
    Contribution(i2,i2,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i2)//' f_'//Bu(i2)
   End Do
  End If    ! GenerationMixing

  !--------
  ! Z - S0
  !--------
  Boson4(1) = mZ
  Boson4(2) = gZ
  Do i1=1,n_S0 
   Isum = Isum + 1
   Boson4(3) = mS0(i1)
   Boson4(4) = gS0(i1)
   coup2(1) = Cpl_CCZ_L(i_out,i_in)
   coup2(2) = Cpl_CCZ_R(i_out,i_in)
   coup2(3) = Conjg(cpl_CCS0_R(i_out,i_in,i1))
   coup2(4) = Conjg(cpl_CCS0_L(i_out,i_in,i1))
   coup2(5) = L_f
   coup2(6) = R_f
   Do i2=1,n_f
    mass(2) = mC(i_out)
    mass(3) = -mf(i2)
    mass(4) = mf(i2)
    coup2(7) = Conjg(cpl_FFS0_R(i2,i2,i1))
    coup2(8) = Conjg(cpl_FFS0_L(i2,i2,i1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsZS04, n_ZS04, resC, check)
    gCffSum(i2,i2,Isum) = 2._dp * Real(resC,dp)
    Contribution(i2,i2,Isum) = 'Z S0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
   End Do
  End Do

  !--------
  ! Z - P0
  !--------
  Do i1=1,n_P0 
   Isum = Isum + 1
   Boson4(3) = mP0(i1)
   Boson4(4) = gP0(i1)
   coup2(1) = Cpl_CCZ_L(i_out,i_in)
   coup2(2) = Cpl_CCZ_R(i_out,i_in)
   coup2(3) = Conjg(cpl_CCP0_R(i_out,i_in,i1))
   coup2(4) = Conjg(cpl_CCP0_L(i_out,i_in,i1))
   coup2(5) = L_f
   coup2(6) = R_f
   Do i2=1,n_f
    mass(2) = mC(i_out)
    mass(3) = -mf(i2)
    mass(4) = mf(i2)
    coup2(7) = Conjg(cpl_FFP0_R(i2,i2,i1))
    coup2(8) = Conjg(cpl_FFP0_L(i2,i2,i1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsZS04, n_ZS04, resC, check)
    gCffSum(i2,i2,Isum) = 2._dp * Real(resC,dp)
    Contribution(i2,i2,Isum) = 'Z P0_'//Bu(i1)//' f_'//Bu(i2)//' f_'//Bu(i2)
   End Do
  End Do

  !--------------------------
  ! Z boson - Sfermion_{xyz} 
  !-------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ

  coup2(1) = Cpl_CCZ_L(i_out, i_in)
  coup2(2) = Cpl_CCZ_R(i_out, i_in)
  coup2(5) = L_f
  coup2(6) = R_f

  If (GenerationMixing) Then
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)
    Do i1=1,n_f! generation
     mass(2) = mf(i1)
     mass(3) = -mf(i1)
     mass(4) = mC(i_out)
     coup2(3) = Conjg( cpl_CFSfp_R(i_in, i1, i2) )
     coup2(4) = Conjg( cpl_CFSfp_L(i_in, i1, i2) )
     coup2(7) = cpl_CFSfp_L(i_out, i1, i2)
     coup2(8) = cpl_CFSfp_R(i_out, i1, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     gCffSum(i1,i1,Isum) = - 2._dp *Real(resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
    End Do
   End Do

  Else If (n_Sfp.Gt.3) Then ! Sneutrinos are special
   Do i1=1,n_f! generation
    i_gen = i1*2-1
    Do i2=i_gen,i_gen+1
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)
     mass(2) = mf(i1)
     mass(3) = -mf(i1)
     mass(4) = mC(i_out)
     coup2(3) = Conjg( cpl_CFSfp_R(i_in, i1, i2) )
     coup2(4) = Conjg( cpl_CFSfp_L(i_in, i1, i2) )
     coup2(7) = cpl_CFSfp_L(i_out, i1, i2)
     coup2(8) = cpl_CFSfp_R(i_out, i1, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     gCffSum(i1,i1,Isum) = - 2._dp *Real(resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
    End Do
   End Do

  Else ! Sneutrinos are special
   Do i1=1,n_f! generation
    Isum = Isum + 1
    Boson4(3) = mSfp(i1)
    Boson4(4) = gSfp(i1)
    mass(2) = mf(i1)
    mass(3) = -mf(i1)
    mass(4) = mC(i_out)
    coup2(3) = Conjg( cpl_CFSfp_R(i_in, i1, i1) )
    coup2(4) = Conjg( cpl_CFSfp_L(i_in, i1, i1) )
    coup2(7) = cpl_CFSfp_L(i_out, i1, i1)
    coup2(8) = cpl_CFSfp_R(i_out, i1, i1)
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZSf8, n_ZSf8, resC, check)
    gCffSum(i1,i1,Isum) = - 2._dp *Real(resC,dp)
    Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i1)//' f_'//Bu(i1)//' f_'//Bu(i1)
   End Do

 End If

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  s-channel
  !-------------------------------------------------
  Do i1 = 1,n_S0
   Boson4(1) = mS0(i1)
   Boson4(2) = gS0(i1)
   Do i2 = i1+1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(i2)
    Boson4(4) = gS0(i2)
    coup2(1) = cpl_CCS0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCS0_R(i_out,i_in,i1)
    coup2(3) = Conjg(cpl_CCS0_R(i_out,i_in,i2))
    coup2(4) = Conjg(cpl_CCS0_L(i_out,i_in,i2))
    Do i3=1,n_f
     mass(2) = mC(i_out)
     mass(3) = -mf(i3)
     mass(4) = mf(i3)
     coup2(5) = cpl_FFS0_L(i3,i3,i1)
     coup2(6) = cpl_FFS0_R(i3,i3,i1)
     coup2(7) = Conjg(cpl_FFS0_R(i3,i3,i2))
     coup2(8) = Conjg(cpl_FFS0_L(i3,i3,i2))
     Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsS0P04, n_S0P04, resC, check)
     gCffSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
       &  'S0_'//Bu(i1)//' S0_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
    End Do
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do i1 = 1,n_S0
   Boson4(1) = mS0(i1)
   Boson4(2) = gS0(i1)
   Do i2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(i2)
    Boson4(4) = gP0(i2)
    coup2(1) = cpl_CCS0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCS0_R(i_out,i_in,i1)
    coup2(3) = Conjg(cpl_CCP0_R(i_out,i_in,i2))
    coup2(4) = Conjg(cpl_CCP0_L(i_out,i_in,i2))
    Do i3=1,n_f
    Contribution(i3,i3,Isum) = &
       &  'S0_'//Bu(i1)//' P0_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     mass(2) = mC(i_out)
     mass(3) = -mf(i3)
     mass(4) = mf(i3)
     coup2(5) = cpl_FFS0_L(i3,i3,i1)
     coup2(6) = cpl_FFS0_R(i3,i3,i1)
     coup2(7) = Conjg(cpl_FFP0_R(i3,i3,i2))
     coup2(8) = Conjg(cpl_FFP0_L(i3,i3,i2))
     Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsS0P04, n_S0P04, resC, check)
     gCffSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    End Do
   End Do
  End Do

  !--------------------------------------
  ! S0 boson  s-channel - S-fermion_{xyz} 
  !--------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_S0
    Boson4(1) = mS0(i1)
    Boson4(2) = gS0(i1)
    coup2(1) = cpl_CCS0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCS0_R(i_out,i_in,i1)
    Do i2 = 1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)
     Do i3 =1,n_f
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mC(i_out)
      coup2(3) = Conjg(cpl_CFSfp_R(i_in, i3,i2))
      coup2(4) = Conjg(cpl_CFSfp_L(i_in, i3,i2))
      coup2(5) = cpl_FFS0_L(i3,i3,i1)
      coup2(6) = cpl_FFS0_R(i3,i3,i1)
      coup2(7) = cpl_CFSfp_L(i_out, i3,i2)
      coup2(8) = cpl_CFSfp_R(i_out, i3,i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      gCffSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
      Contribution(i3,i3,Isum) = &
        &  'S0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     End Do
    End Do
   End Do

  Else If (n_Sfp.Gt.3) Then ! Sneutrinos are Special
   Do i1 = 1,n_S0
    Boson4(1) = mS0(i1)
    Boson4(2) = gS0(i1)
    coup2(1) = cpl_CCS0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCS0_R(i_out,i_in,i1)
    Do i3 =1,n_f
     i_gen = 2*i3-1
     Do i2 =i_gen,i_gen+1
      Isum = Isum + 1
      Boson4(3) = mSfp(i2)
      Boson4(4) = gSfp(i2)
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mC(i_out)
      coup2(3) = Conjg(cpl_CFSfp_R(i_in, i3,i2))
      coup2(4) = Conjg(cpl_CFSfp_L(i_in, i3,i2))
      coup2(5) = cpl_FFS0_L(i3,i3,i1)
      coup2(6) = cpl_FFS0_R(i3,i3,i1)
      coup2(7) = cpl_CFSfp_L(i_out, i3,i2)
      coup2(8) = cpl_CFSfp_R(i_out, i3,i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      gCffSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
      Contribution(i3,i3,Isum) = &
        &  'S0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     End Do
    End Do
   End Do

  Else ! Sneutrinos are Special
   Do i1 = 1,n_S0
    Boson4(1) = mS0(i1)
    Boson4(2) = gS0(i1)
    coup2(1) = cpl_CCS0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCS0_R(i_out,i_in,i1)
    Do i3 =1,n_f
     Isum = Isum + 1
     Boson4(3) = mSfp(i1)
     Boson4(4) = gSfp(i1)
     mass(2) = mf(i3)
     mass(3) = -mf(i3)
     mass(4) = mC(i_out)
     coup2(3) = Conjg(cpl_CFSfp_R(i_in, i3,i1))
     coup2(4) = Conjg(cpl_CFSfp_L(i_in, i3,i1))
     coup2(5) = cpl_FFS0_L(i3,i3,i1)
     coup2(6) = cpl_FFS0_R(i3,i3,i1)
     coup2(7) = cpl_CFSfp_L(i_out, i3,i1)
     coup2(8) = cpl_CFSfp_R(i_out, i3,i1)
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsS0Sf8, n_S0Sf8, resC, check)
     gCffSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
       &  'S0_'//Bu(i1)//' Sf_'//Bu(i1)//' f_'//Bu(i3)//' f_'//Bu(i3)
    End Do
   End Do

  End If
  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do i1 = 1,n_P0
   Boson4(1) = mP0(i1)
   Boson4(2) = gP0(i1)
   Do i2 = i1+1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(i2)
    Boson4(4) = gP0(i2)
    coup2(1) = cpl_CCP0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCP0_R(i_out,i_in,i1)
    coup2(3) = Conjg(cpl_CCP0_R(i_out,i_in,i2))
    coup2(4) = Conjg(cpl_CCP0_L(i_out,i_in,i2))
    Do i3=1,n_f
    Contribution(i3,i3,Isum) = &
       &  'P0_'//Bu(i1)//' P0_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     mass(2) = mC(i_out)
     mass(3) = -mf(i3)
     mass(4) = mf(i3)
     coup2(5) = cpl_FFP0_L(i3,i3,i1)
     coup2(6) = cpl_FFP0_R(i3,i3,i1)
     coup2(7) = Conjg(cpl_FFP0_R(i3,i3,i2))
     coup2(8) = Conjg(cpl_FFP0_L(i3,i3,i2))
     Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsS0P04, n_S0P04, resC, check)
     gCffSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    End Do
   End Do
  End Do

  !--------------------------------------
  ! P0 boson  s-channel - S-fermion_{xyz} 
  !--------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_P0
    Boson4(1) = mP0(i1)
    Boson4(2) = gP0(i1)
    coup2(1) = cpl_CCP0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCP0_R(i_out,i_in,i1)
    Do i2 = 1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)
     Do i3 =1,n_f
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mC(i_out)
      coup2(3) = Conjg(cpl_CFSfp_R(i_in, i3,i2))
      coup2(4) = Conjg(cpl_CFSfp_L(i_in, i3,i2))
      coup2(5) = cpl_FFP0_L(i3,i3,i1)
      coup2(6) = cpl_FFP0_R(i3,i3,i1)
      coup2(7) = cpl_CFSfp_L(i_out, i3,i2)
      coup2(8) = cpl_CFSfp_R(i_out, i3,i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      gCffSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
      Contribution(i3,i3,Isum) = &
         &  'P0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
     End Do
    End Do
   End Do

  Else If (n_Sfp.Gt.3) Then ! Sneutrinos are special
   Do i1 = 1,n_P0
    Boson4(1) = mP0(i1)
    Boson4(2) = gP0(i1)
    coup2(1) = cpl_CCP0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCP0_R(i_out,i_in,i1)
    Do i3 =1,n_f
     i_gen = 2*i3-1
     Do i2 =i_gen,i_gen+1
      Isum = Isum + 1
      Contribution(i3,i3,Isum) = &
         &  'P0_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' f_'//Bu(i3)
      Boson4(3) = mSfp(i2)
      Boson4(4) = gSfp(i2)
      mass(2) = mf(i3)
      mass(3) = -mf(i3)
      mass(4) = mC(i_out)
      coup2(3) = Conjg(cpl_CFSfp_R(i_in, i3,i2))
      coup2(4) = Conjg(cpl_CFSfp_L(i_in, i3,i2))
      coup2(5) = cpl_FFP0_L(i3,i3,i1)
      coup2(6) = cpl_FFP0_R(i3,i3,i1)
      coup2(7) = cpl_CFSfp_L(i_out, i3,i2)
      coup2(8) = cpl_CFSfp_R(i_out, i3,i2)
      Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0Sf8, n_S0Sf8, resC, check)
      gCffSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     End Do
    End Do
   End Do

  Else  ! Sneutrinos are special
   Do i1 = 1,n_P0
    Boson4(1) = mP0(i1)
    Boson4(2) = gP0(i1)
    coup2(1) = cpl_CCP0_L(i_out,i_in,i1)
    coup2(2) = cpl_CCP0_R(i_out,i_in,i1)
    Do i3 =1,n_f
     Isum = Isum + 1
     Boson4(3) = mSfp(i3)
     Boson4(4) = gSfp(i3)
     mass(2) = mf(i3)
     mass(3) = -mf(i3)
     mass(4) = mC(i_out)
     coup2(3) = Conjg(cpl_CFSfp_R(i_in, i3,i3))
     coup2(4) = Conjg(cpl_CFSfp_L(i_in, i3,i3))
     coup2(5) = cpl_FFP0_L(i3,i3,i1)
     coup2(6) = cpl_FFP0_R(i3,i3,i1)
     coup2(7) = cpl_CFSfp_L(i_out, i3,i3)
     coup2(8) = cpl_CFSfp_R(i_out, i3,i3)
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsS0Sf8, n_S0Sf8, resC, check)
     gCffSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
        &  'P0_'//Bu(i1)//' Sf_'//Bu(i3)//' f_'//Bu(i3)//' f_'//Bu(i3)
    End Do
   End Do

  End If

  !----------------------------------
  ! Sfermion - Sfermion, non-diagonal
  !----------------------------------
  If (GenerationMixing) Then
   Do i3=1,n_Sfp-1
    Boson4(1) = mSfp(i3)
    Boson4(2) = gSfp(i3)
    Do i4=i3+1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i4)
     Boson4(4) = gSfp(i4)
     Do i1 = 1,n_f ! fermions
      Do i2 = 1,n_f
       mass(2) = mf(i1)
       mass(3) = -mf(i2)
       mass(4) = mC(i_out)
       coup2(1) = cpl_CFSfp_L(i_in,i1,i3)
       coup2(2) = cpl_CFSfp_R(i_in,i1,i3)
       coup2(3) = Conjg(cpl_CFSfp_R(i_in,i1,i4))
       coup2(4) = Conjg(cpl_CFSfp_L(i_in,i1,i4))
       coup2(5) = Conjg(cpl_CFSfp_R(i_out,i2,i3))
       coup2(6) = Conjg(cpl_CFSfp_L(i_out,i2,i3))
       coup2(7) = cpl_CFSfp_L(i_out,i2,i4)
       coup2(8) = cpl_CFSfp_R(i_out,i2,i4)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsCSf4, n_CSf4, resC, check)
       gCffSum(i1,i2,Isum) = 2._dp * Real(resC,dp)
       Contribution(i1,i2,Isum) = &
         &  'tt Sf_'//Bu(i3)//' Sf_'//Bu(i4)//' f_'//Bu(i1)//' f_'//Bu(i2)
      End Do
     End Do
    End Do
   End Do

  Else If (n_Sfp.Gt.3) Then ! .no.GenerationMixing and no sneutrinos

   Do i1 = 1,n_f
    i2 = 2*i1 - 1
    i3 = 2*i1
    Isum = Isum + 1
    Boson4(1) = mSfp(i2)
    Boson4(2) = gSfp(i2)
    Boson4(3) = mSfp(i3)
    Boson4(4) = gSfp(i3)
    mass(2) = mf(i1)
    mass(3) = -mf(i1)
    mass(4) = mC(i_out)
    coup2(1) = cpl_CFSfp_L(i_in,i1,i2)
    coup2(2) = cpl_CFSfp_R(i_in,i1,i2)
    coup2(3) = Conjg(cpl_CFSfp_R(i_in,i1,i3))
    coup2(4) = Conjg(cpl_CFSfp_L(i_in,i1,i3))
    coup2(5) = Conjg(cpl_CFSfp_R(i_out,i1,i2))
    coup2(6) = Conjg(cpl_CFSfp_L(i_out,i1,i2))
    coup2(7) = cpl_CFSfp_L(i_out,i1,i3)
    coup2(8) = cpl_CFSfp_R(i_out,i1,i3)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsCSf4, n_CSf4, resC, check)
    gCffSum(i1,i1,Isum) = 2._dp * Real(resC,dp)
    Contribution(i1,i1,Isum) = &
        &  'tt Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i1)
   End Do

  End If

  !----------
  ! Summing
  !----------
  gCff = 0._dp
  Do i1=1,n_f
   If (GenerationMixing) Then
    Do i2=1,n_f
     gCff(i1,i2) = Sum( gCffSum(i1,i2,1:Isum) )
     If (gCff(i1,i2).Lt.0._dp) Then
      Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//state//') < 0 :' &
      & ,i1,i2,gCff(i1,i2)
      Write (ErrCan,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gCffSum(i1,i2,i3).Ne.0._dp) &
        &      Write (ErrCan,*) Contribution(i1,i2,i3),gCffSum(i1,i2,i3)
      End Do
      gCff(i1,i2) = 0._dp
     End If
    End Do

   Else
    gCff(i1,i1) = Sum( gCffSum(i1,i1,1:Isum) )
    If (gCff(i1,i1).Lt.0._dp) Then
     Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//state//') < 0 :' &
      & ,i1,i1,gCff(i1,i1)
     Write (ErrCan,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gCffSum(i1,i1,i3).Ne.0._dp) &
        & Write (ErrCan,*) Contribution(i1,i1,i3),gCffSum(i1,i1,i3)
     End Do
     gCff(i1,i1) = 0._dp
    End If
   End If
  End Do

  gCff = fac * gCff

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gCffSum = gCffSum * fac

   If (GenerationMixing) Then
    Do i1=1,n_f
     Do i2=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//state//') :' &
      & ,i1,i2,gCff(i1,i2)
     Write (n_out,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gCffSum(i1,i2,i3).Ne.0._dp) &
        &       Write (n_out,*) Contribution(i1,i2,i3),gCffSum(i1,i2,i3)
      End Do
      Write (n_out,*) ' '
     End Do
    End Do

   Else
    Do i1=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//state//') :' &
      & ,i1,i1,gCff(i1,i1)
     Write (n_out,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gCffSum(i1,i1,i3).Ne.0._dp) &
        &      Write (n_out,*) Contribution(i1,i1,i3),gCffSum(i1,i1,i3)
     End Do
     Write (n_out,*) ' '
    End Do

   End If
  End If

  Deallocate( gCffSum, Contribution)

  Iname = Iname - 1

 End Subroutine ChimToChimff


 Subroutine ChimToChimNuNu(i_in, i_out, mC, mZ, gZ, Cpl_CCZ_L, Cpl_CCZ_R, n_f &
    & , L_f, R_f, IntegralsZ4, n_Z4, mSfp, gSfp, cpl_CFSfp_L, cpl_CFSfp_R     &
    & , IntegralsSf4, n_Sf4, IntegralsZSf8, n_ZSf8, IntegralsCSf4, n_CSf4     &
    & , deltaM, epsI, GenerationMixing, check, fac, gCff, WriteContribution   &
    & , n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a neutralino to another Neutralino + fermion pair
 ! Written by Werner Porod, 04.06.2001
 !--------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in, i_out, n_f
  Integer, Intent(inout) :: n_Z4, n_Sf4, n_ZSf8, n_CSf4

  Real(dp), Intent(in) :: mC(:), mZ, gZ, L_f, R_f, mSfp(:), gSfp(:) &
      & , deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsSf4(:,:)
  Complex(dp), Intent(in) :: Cpl_CCZ_L(:,:), Cpl_CCZ_R(:,:)             &
      & , cpl_CFSfp_L(:,:,:), cpl_CFSfp_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsZSf8(:,:), IntegralsCSf4(:,:)
  Logical, Intent(in) :: GenerationMixing, check

  Real(dp), Intent(out) :: gCff(:,:)

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution

  Integer :: Isum, i1, i2, i3, i4, i_gen, n_Sfp
  Real(dp) :: Boson2(2), mass(4), resR, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8)

  Real(dp), Allocatable :: gCffSum(:,:,:)
  Character(len=20), Allocatable :: Contribution(:,:,:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChimToChimNuNu'

  mass(1) = mC(i_in)
  n_Sfp = Size( mSfp )
  Isum = (1 + n_Sfp)**2
  Allocate( gCffSum(n_f, n_f, Isum) )
  Allocate( Contribution(n_f, n_f, Isum) )
   
  gCffSum = 0._dp
  Contribution = ' '

  !-----
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ
  coup1(1) = Cpl_CCZ_L(i_out,i_in)
  coup1(2) = Cpl_CCZ_R(i_out,i_in)
  coup1(3) = L_f
  coup1(4) = R_f
  mass(2) = mC(i_out)
  mass(3) = 0._dp
  mass(4) = 0._dp
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resR, check )
  Do i1=1,n_f
   Contribution(i1,i1,1) = 'Z f_'//Bu(i1)//' f_'//Bu(i1)
   gCffSum(i1,i1,Isum) = resR
  End Do

  !-----------------------------
  ! sfermion sfermion, diagonal
  !-----------------------------
  mass(2) = 0._dp
  mass(3) = 0._dp
  mass(4) = mC(i_out)        
  If (GenerationMixing) Then
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    Do i1=1,n_f
     Do i3=1,n_f
      coup1(1) = cpl_CFSfp_L(i_in,i1,i2)
      coup1(2) = cpl_CFSfp_R(i_in,i1,i2)
      coup1(3) = Conjg(cpl_CFSfp_R(i_out,i3,i2))
      coup1(4) = Conjg(cpl_CFSfp_L(i_out,i3,i2))
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gCffSum(i1,i3,Isum) = resR
      Contribution(i1,i3,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i3)
     End Do ! i3 fermions
    End Do  ! i1 fermions
   End Do   ! i2 sfermions

  Else
   Do i2=1,n_f  ! fermion generation
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    coup1(1) = cpl_CFSfp_L(i_in,i2,i2)
    coup1(2) = cpl_CFSfp_R(i_in,i2,i2)
    coup1(3) = Conjg(cpl_CFSfp_R(i_out,i2,i2))
    coup1(4) = Conjg(cpl_CFSfp_L(i_out,i2,i2))
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsSf4, n_Sf4, resR, check)
    gCffSum(i2,i2,Isum) = resR
    Contribution(i2,i2,Isum) = 'Sf^2_'//Bu(i2)//' f_'//Bu(i2)//' f_'//Bu(i2)
   End Do
  End If    ! GenerationMixing

  !--------------------------
  ! Z boson - Sfermion_{xyz} 
  !-------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ

  mass(2) = 0._dp
  mass(3) = 0._dp
  mass(4) = mC(i_out)
  coup2(1) = Cpl_CCZ_L(i_out, i_in)
  coup2(2) = Cpl_CCZ_R(i_out, i_in)
  coup2(5) = L_f
  coup2(6) = R_f

  If (GenerationMixing) Then
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)
    Do i1=1,n_f! generation
     coup2(3) = Conjg( cpl_CFSfp_R(i_in, i1, i2) )
     coup2(4) = Conjg( cpl_CFSfp_L(i_in, i1, i2) )
     coup2(7) = cpl_CFSfp_L(i_out, i1, i2)
     coup2(8) = cpl_CFSfp_R(i_out, i1, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     gCffSum(i1,i1,Isum) = - 2._dp *Real(resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
    End Do
   End Do

  Else
   Do i1=1,n_f! generation
    i_gen = i1*2-1
    Do i2=i_gen,i_gen+1
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)
     coup2(3) = Conjg( cpl_CFSfp_R(i_in, i1, i2) )
     coup2(4) = Conjg( cpl_CFSfp_L(i_in, i1, i2) )
     coup2(7) = cpl_CFSfp_L(i_out, i1, i2)
     coup2(8) = cpl_CFSfp_R(i_out, i1, i2)
     Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsZSf8, n_ZSf8, resC, check)
     gCffSum(i1,i1,Isum) = - 2._dp *Real(resC,dp)
     Contribution(i1,i1,Isum) = 'Z Sf_'//Bu(i2)//' f_'//Bu(i1)//' f_'//Bu(i1)
    End Do
   End Do

  End If

  !----------------------------------
  ! Sfermion - Sfermion, non-diagonal
  !----------------------------------
  mass(2) = 0._dp
  mass(3) = 0._dp
  mass(4) = mC(i_out)
  If (GenerationMixing) Then
   Do i3=1,2*n_f-1
    Boson4(1) = mSfp(i3)
    Boson4(2) = gSfp(i3)
    Do i4=i3+1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSfp(i4)
     Boson4(4) = gSfp(i4)
     Do i1 = 1,n_f ! fermions
      Do i2 = 1,n_f
       coup2(1) = cpl_CFSfp_L(i_in,i1,i3)
       coup2(2) = cpl_CFSfp_R(i_in,i1,i3)
       coup2(3) = Conjg(cpl_CFSfp_R(i_in,i1,i4))
       coup2(4) = Conjg(cpl_CFSfp_L(i_in,i1,i4))
       coup2(5) = Conjg(cpl_CFSfp_R(i_out,i2,i3))
       coup2(6) = Conjg(cpl_CFSfp_L(i_out,i2,i3))
       coup2(7) = cpl_CFSfp_L(i_out,i2,i4)
       coup2(8) = cpl_CFSfp_R(i_out,i2,i4)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsCSf4, n_CSf4, resC, check)
       gCffSum(i1,i2,Isum) = 2._dp * Real(resC,dp)
       Contribution(i1,i2,Isum) = &
         &  'tt Sf_'//Bu(i3)//' Sf_'//Bu(i4)//' f_'//Bu(i1)//' f_'//Bu(i2)
      End Do
     End Do
    End Do
   End Do

  Else

   Do i1 = 1,n_f
    i2 = 2*i1 - 1
    i3 = 2*i1
    Isum = Isum + 1
    Boson4(1) = mSfp(i2)
    Boson4(2) = gSfp(i2)
    Boson4(3) = mSfp(i3)
    Boson4(4) = gSfp(i3)
    coup2(1) = cpl_CFSfp_L(i_in,i1,i2)
    coup2(2) = cpl_CFSfp_R(i_in,i1,i2)
    coup2(3) = Conjg(cpl_CFSfp_R(i_in,i1,i3))
    coup2(4) = Conjg(cpl_CFSfp_L(i_in,i1,i3))
    coup2(5) = Conjg(cpl_CFSfp_R(i_out,i1,i2))
    coup2(6) = Conjg(cpl_CFSfp_L(i_out,i1,i2))
    coup2(7) = cpl_CFSfp_L(i_out,i1,i3)
    coup2(8) = cpl_CFSfp_R(i_out,i1,i3)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsCSf4, n_CSf4, resC, check)
    gCffSum(i1,i1,Isum) = 2._dp * Real(resC,dp)
    Contribution(i1,i1,Isum) = &
        &  'tt Sf_'//Bu(i2)//' Sf_'//Bu(i3)//' f_'//Bu(i1)//' f_'//Bu(i1)
   End Do

  End If

  !----------
  ! Summing
  !----------
  gCff = 0._dp
  Do i1=1,n_f
   If (GenerationMixing) Then
    Do i2=1,n_f
     gCff(i1,i2) = Sum( gCffSum(i1,i2,1:Isum) )
     If (gCff(i1,i2).Lt.0._dp) Then
      Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//' nu nu) < 0 :' &
      & ,i1,i2,gCff(i1,i2)
      Write (ErrCan,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gCffSum(i1,i2,i3).Ne.0._dp) &
        &      Write (ErrCan,*) Contribution(i1,i2,i3),gCffSum(i1,i2,i3)
      End Do
      gCff(i1,i2) = 0._dp
     End If
    End Do

   Else
    gCff(i1,i1) = Sum( gCffSum(i1,i1,1:Isum) )
    If (gCff(i1,i1).Lt.0._dp) Then
     Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//' nu nu) < 0 :' &
      & ,i1,i1,gCff(i1,i1)
     Write (ErrCan,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gCffSum(i1,i1,i3).Ne.0._dp) &
        & Write (ErrCan,*) Contribution(i1,i1,i3),gCffSum(i1,i1,i3)
     End Do
     gCff(i1,i1) = 0._dp
    End If
   End If
  End Do

  gCff = fac * gCff

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gCffSum = gCffSum * fac

   If (GenerationMixing) Then
    Do i1=1,n_f
     Do i2=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//' nu nu) :' &
      & ,i1,i2,gCff(i1,i2)
     Write (n_out,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gCffSum(i1,i2,i3).Ne.0._dp) &
        &       Write (n_out,*) Contribution(i1,i2,i3),gCffSum(i1,i2,i3)
      End Do
      Write (n_out,*) ' '
     End Do
    End Do

   Else
    Do i1=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi^-_'//Bu(i_out)//' nu nu) :' &
      & ,i1,i1,gCff(i1,i1)
     Write (n_out,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gCffSum(i1,i1,i3).Ne.0._dp) &
        &      Write (n_out,*) Contribution(i1,i1,i3),gCffSum(i1,i1,i3)
     End Do
     Write (n_out,*) ' '
    End Do

   End If
  End If

  Deallocate( gCffSum, Contribution)

  Iname = Iname - 1

 End Subroutine ChimToChimNuNu

 Subroutine ChimToChi0ffp(i_in, i_out, mC, mN, n_f, mf, mfp, mW, gW          &
    & , Cpl_CNW_L, Cpl_CNW_R, Cpl_FpFW, mSpm, gSpm, Cpl_SmpCN_L, Cpl_SmpCN_R &
    & , Cpl_SmpFFp_L, Cpl_SmpFFp_R, mSf, gSf, cpl_FNSf_L, cpl_FNSf_R         &
    & , cpl_CFpSf_L, cpl_CFpSf_R, mSfp, gSfp, cpl_FpNSfp_L, cpl_FpNSfp_R     &
    & , cpl_CFSfp_L, cpl_CFSfp_R                                             &
    & , IntegralsW4, n_W4, IntegralsSpm4, n_Spm4, IntegralsSf4, n_Sf4        &
    & , IntegralsWSpm4, n_WSpm, IntegralsWSf8, n_WSf8, IntegralsSpmC4        &
    & , n_SpmC4, IntegralsSpmSf8, n_SpmSf8, IntegralsSfC4, n_SfC4            &
    & , IntegralsSf8, n_Sf8, deltaM, epsI                                    &
    & , GenerationMixing, check, fac, gNffp, WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a Chargino to a neutralino + fermion pair
 ! Written by Werner Porod, 30.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !--------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in, i_out, n_f
  Integer, Intent(inout) :: n_W4, n_Spm4, n_Sf4, n_WSpm, n_WSf8, n_SpmC4 &
     & , n_SpmSf8, n_SfC4, n_Sf8
  
  Real(dp), Intent(in) :: mN(:), mC(:), mf(:), mfp(:), mW, gW, mSpm(:)     &
      & , gSpm(:), mSf(:), gSf(:), mSfp(:), gSfp(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsW4(:,:), IntegralsSpm4(:,:)         &
      & , IntegralsSf4(:,:)

  Complex(dp), Intent(in) :: Cpl_CNW_L(:,:), Cpl_CNW_R(:,:), Cpl_FpFW(:,:)  &
      & , Cpl_SmpCN_L(:,:,:), Cpl_SmpCN_R(:,:,:), Cpl_SmpFFp_L(:,:,:)       &
      & , Cpl_SmpFFp_R(:,:,:), cpl_FNSf_L(:,:,:), cpl_FNSf_R(:,:,:)         &
      & , cpl_CFpSf_L(:,:,:), cpl_CFpSf_R(:,:,:), cpl_FpNSfp_L(:,:,:)       &
      & , cpl_FpNSfp_R(:,:,:), cpl_CFSfp_L(:,:,:), cpl_CFSfp_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsWSpm4(:,:), IntegralsWSf8(:,:)   &
      & , IntegralsSpmC4(:,:), IntegralsSpmSf8(:,:), IntegralsSfC4(:,:)   &
      & , IntegralsSf8(:,:)

  Logical, Intent(in) :: GenerationMixing, check

  Real(dp), Intent(out) :: gNffp(:,:)

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution

  Integer :: Isum, i1, i2, i3, i4, n_Spm, n_Sfp
  Real(dp) :: Boson2(2), mass(4), resR, Boson4(4)
  Complex(dp) :: coup1(4), resC, coup2(8)

  Real(dp), Allocatable :: gNffpSum(:,:,:)
  Character(len=20), Allocatable :: Contribution(:,:,:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChimToChi0ffp'

  If (n_f.Le.0) Then ! useful for the case of R-parity violation
   gNffp = 0._dp
   Iname = Iname - 1
   Return
  End If

  mass(1) = mC(i_in)

  n_Spm = Size( mSpm )
  n_Sfp = Size( mSfp )

  !-------------------------------------------
  ! to cover the case of sneutrinos correctly
  !-------------------------------------------
  If (n_Sfp.Le.3) Then
   n_Sfp = n_f
  Else
   n_Sfp = 2 * n_f
  End If

  Isum = (1 + n_Spm + 4 * n_f)**2
  Allocate( gNffpSum(n_f, n_f, Isum) )
  Allocate( Contribution(n_f, n_f, Isum) )

  gNffpSum = 0._dp
  Contribution = ' '

  !-----
  ! W W
  !-----
  Isum = 1
  Boson2(1) = mW
  Boson2(2) = gW

  coup1(1) = Conjg( Cpl_CNW_L(i_in,i_out) )
  coup1(2) = Conjg( Cpl_CNW_R(i_in,i_out) )
  coup1(4) = 0._dp

  If (GenerationMixing) Then
   Do i2=1,n_f
    Do i3=1,n_f
     mass(2) = mN(i_out)
     mass(3) = -mfp(i3)
     coup1(3) = Conjg(cpl_FpFW(i3,i2) )
     mass(4) = mf(i2)
     Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                          &, n_W4, resR, check)
     gNffpSum(i2, i3, Isum) = resR
     Contribution(i2,i3,1) = 'W f_'//Bu(i2)//' fp_'//Bu(i3)
    End Do
   End Do

  Else
   Do i2=1,n_f
    mass(2) = mN(i_out)
    mass(3) = -mfp(i2)
    coup1(3) = Conjg(cpl_FpFW(i2,i2)) 
    mass(4) = mf(i2)
    Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                         &, n_W4, resR, check)
    gNffpSum(i2, i2, Isum) = resR
    Contribution(i2,i2,1) = 'W f_'//Bu(i2)//' fp_'//Bu(i2)
   End Do
  End If

  !-------------------
  ! S^- S^-, diagonal
  !-------------------
  Do i1=1,n_Spm
   Isum = Isum + 1

   coup1(1) = Conjg( Cpl_SmpCN_R(i1, i_in,i_out) )
   coup1(2) = Conjg( Cpl_SmpCN_L(i1, i_in,i_out) )
   Boson2(1) = mSpm(i1)
   Boson2(2) = gSpm(i1)

   If (GenerationMixing) Then
    Do i2=1,n_f
     Do i3=1,n_f
      mass(2) = mN(i_out)
      mass(3) = -mfp(i3)
      mass(4) = mf(i2)
      coup1(3) = cpl_SmpFFp_L(i1,i2,i3)
      coup1(4) = cpl_SmpFFp_R(i1,i2,i3)
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSpm4 ,n_Spm4, resR, check)
      gNffpSum(i2,i3,Isum) = resR
      Contribution(i2,i3,Isum) = 'S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i3)
     End Do
    End Do

   Else
    Do i2=1,n_f
     mass(2) = mN(i_out)
     mass(3) = -mfp(i2)
     mass(4) = mf(i2)
     coup1(3) = cpl_SmpFFp_L(i1,i2,i2)
     coup1(4) = cpl_SmpFFp_R(i1,i2,i2)
     Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                          &, IntegralsSpm4 ,n_Spm4, resR, check)
     gNffpSum(i2,i2,Isum) = resR
     Contribution(i2,i2,Isum) = 'S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i2)
    End Do

   End If
  End Do

  !-------------------
  ! Sfp Sfp, diagonal
  !-------------------
  
  If (GenerationMixing) Then
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    Do i3=1,n_f
     coup1(1) = cpl_CFSfp_L(i_in,i3,i2)
     coup1(2) = cpl_CFSfp_R(i_in,i3,i2)
     Do i1=1,n_f
      mass(2) = mf(i3)
      mass(3) = -mfp(i1)
      mass(4) = mN(i_out)
      coup1(3) = Conjg(cpl_FpNSfp_R(i1,i_out,i2))
      coup1(4) = Conjg(cpl_FpNSfp_L(i1,i_out,i2))
       Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gNffpSum(i3,i1,Isum) = resR
      Contribution(i3,i1,Isum) = 'Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i1)
     End Do
    End Do
   End Do

  Else

   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)

    If (n_Sfp.Le.3) Then
     i1 = i2
    Else
     i1 = (i2+1)/2
    End If

    mass(2) = mf(i1)
    mass(3) = -mfp(i1)
    mass(4) = mN(i_out)
    coup1(1) = cpl_CFSfp_L(i_in,i1,i2)
    coup1(2) = cpl_CFSfp_R(i_in,i1,i2)
    coup1(3) = Conjg(cpl_FpNSfp_R(i1,i_out,i2))
    coup1(4) = Conjg(cpl_FpNSfp_L(i1,i_out,i2))
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsSf4, n_Sf4, resR, check)
    gNffpSum(i1,i1,Isum) = resR
    Contribution(i1,i1,Isum) = 'Sfp_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do

  End If

  !-------------------------
  ! Sf Sf, diagonal
  !-------------------------
  If (GenerationMixing) Then
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    Do i3=1,n_f
     coup1(1) = Conjg(cpl_CFpSf_R(i_in,i3,i2))
     coup1(2) = Conjg(cpl_CFpSf_L(i_in,i3,i2))
     Do i1=1,n_f
      mass(2) = mfp(i3)
      mass(3) = - mN(i_out)
      mass(4) = mf(i1)
      coup1(3) = cpl_FNSf_L(i1,i_out,i2)
      coup1(4) = cpl_FNSf_R(i1,i_out,i2)
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gNffpSum(i1,i3,Isum) = resR
      Contribution(i3,i1,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i1)
     End Do
    End Do
   End Do
  Else

   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    i1 = (i2+1)/2

    mass(2) = mfp(i1)
    mass(3) = - mN(i_out)
    mass(4) = mf(i1)
    coup1(1) = Conjg(cpl_CFpSf_R(i_in,i1,i2))
    coup1(2) = Conjg(cpl_CFpSf_L(i_in,i1,i2))
    coup1(3) = cpl_FNSf_L(i1,i_out,i2)
    coup1(4) = cpl_FNSf_R(i1,i_out,i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsSf4, n_Sf4, resR, check)
    gNffpSum(i1,i1,Isum) = resR
    Contribution(i1,i1,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do
  End If

  !-------------------------------
  ! W boson - S+ boson  s-channel
  !-------------------------------
  Boson4(1) = mW
  Boson4(2) = gW
  coup2(1) = Conjg( Cpl_CNW_L(i_in,i_out) )
  coup2(2) = Conjg( Cpl_CNW_R(i_in,i_out) )
  coup2(6) = 0._dp

  If (GenerationMixing) Then
   Do i1=1,n_Spm 
    Isum = Isum + 1
    Boson4(3) = mSpm(i1)
    Boson4(4) = gSpm(i1)

    coup2(3) = Cpl_SmpCN_L(i1,i_in,i_out)
    coup2(4) = Cpl_SmpCN_R(i1,i_in,i_out)
    Do i2=1,n_f
     Do i3=1,n_f
      mass(2) = mN(i_out)
      mass(3) = -mfp(i3)
      mass(4) = mf(i2)
      coup2(5) = Conjg( cpl_FpFW(i3,i2) )
      coup2(7) = Conjg( cpl_SmpFFp_R(i1,i2,i3) )
      coup2(8) = Conjg( cpl_SmpFFp_L(i1,i2,i3) )
      Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsWSpm4, n_WSpm, resC, check)
      gNffpSum(i2,i3,Isum) = 2._dp * Real(resC,dp)
      Contribution(i2,i3,Isum) = &
        & 'W S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i3)
     End Do
    End Do
   End Do

  Else
   Do i1=1,n_Spm 
    Isum = Isum + 1
    Boson4(3) = mSpm(i1)
    Boson4(4) = gSpm(i1)

    coup2(3) = Cpl_SmpCN_L(i1,i_in,i_out)
    coup2(4) = Cpl_SmpCN_R(i1,i_in,i_out)
    Do i2=1,n_f
     mass(2) = mN(i_out)
     mass(3) = -mfp(i2)
     mass(4) = mf(i2)
     coup2(5) = Conjg( cpl_FpFW(i2,i2) )
     coup2(7) = Conjg( cpl_SmpFFp_R(i1,i2,i2) )
     coup2(8) = Conjg( cpl_SmpFFp_L(i1,i2,i2) )
     Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsWSpm4, n_WSpm, resC, check)
     gNffpSum(i2,i2,Isum) = 2._dp * Real(resC,dp)
     Contribution(i2,i2,Isum) = &
        & 'W S^-_'//Bu(i1)//' f_'//Bu(i2)//' fp_'//Bu(i2)
    End Do
   End Do

  End If

  !---------------------------------
  ! W boson - Sfp_{xyz}  t-channel
  !---------------------------------
  coup2(1) = Conjg( Cpl_CNW_L(i_in,i_out) )
  coup2(2) = Conjg( Cpl_CNW_R(i_in,i_out) )
  coup2(6) = 0._dp 

  If (GenerationMixing) Then
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)
    Do i3=1,n_f
      coup2(3) = Conjg( cpl_CFSfp_R(i_in, i3, i2) )
      coup2(4) = Conjg( cpl_CFSfp_L(i_in, i3, i2) )

     Do i1=1,n_f
      mass(2) = mf(i3)
      mass(3) = -mfp(i1)
      mass(4) = mN(i_out)

      coup2(5) = Conjg( cpl_FpFW(i1,i3) )
      coup2(7) = cpl_FpNSfp_L(i1, i_out, i2)
      coup2(8) = cpl_FpNSfp_R(i1, i_out, i2)
      Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                 &, IntegralsWSf8, n_WSf8, resC, check)
      gNffpSum(i3,i1,Isum) = - 2._dp * Real(resC,dp)
      Contribution(i3,i1,Isum) = &
        & 'W Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i1)
     End Do
    End Do
   End Do

  Else
   Do i2=1,n_Sfp
    Isum = Isum + 1
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)

    If (n_Sfp.Le.3) Then
     i1 = i2
    Else
     i1 = (i2+1)/2
    End If

    mass(2) = mf(i1)
    mass(3) = -mfp(i1)
    mass(4) = mN(i_out)

    coup2(3) = Conjg( cpl_CFSfp_R(i_in, i1, i2) )
    coup2(4) = Conjg( cpl_CFSfp_L(i_in, i1, i2) )
    coup2(5) = Conjg( cpl_FpFW(i1,i1) )
    coup2(7) = cpl_FpNSfp_L(i1, i_out, i2)
    coup2(8) = cpl_FpNSfp_R(i1, i_out, i2)
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                               &, IntegralsWSf8, n_WSf8, resC, check)
    gNffpSum(i1,i1,Isum) = - 2._dp * Real(resC,dp)
    Contribution(i1,i1,Isum) = &
      & 'W Sfp_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do

  End If

  !---------------------------------
  ! W boson - Sf_{xyz}  u-channel
  !---------------------------------
  coup2(1) = Conjg( Cpl_CNW_L(i_in,i_out) )
  coup2(2) = Conjg( Cpl_CNW_R(i_in,i_out) )
  coup2(5) = 0._dp 

  If (GenerationMixing) Then
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)

    Do i3=1,n_f
     coup2(3) = cpl_CFpSf_L(i_in, i3, i2)
     coup2(4) = cpl_CFpSf_R(i_in, i3, i2)

     Do i1=1,n_f
      mass(2) = mfp(i3)
      mass(3) = -mf(i1)
      mass(4) = mN(i_out)

      coup2(6) = cpl_FpFW(i3,i1) 
      coup2(7) = Conjg( cpl_FNSf_R(i1, i_out, i2) )
      coup2(8) = Conjg( cpl_FNSf_L(i1, i_out, i2) )
      Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                                &, IntegralsWsf8, n_WSf8, resC, check)
      gNffpSum(i1,i3,Isum) = 2._dp * Real(resC,dp)
      Contribution(i1,i3,Isum) = &
        & 'W Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i3)
     End Do
    End Do
   End Do

  Else
   Do i2=1,2*n_f
    Isum = Isum + 1
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)

    i1 = (i2+1)/2

    mass(2) = mf(i1)
    mass(3) = -mfp(i1)
    mass(4) = mN(i_out)

    coup2(3) = cpl_CFpSf_L(i_in, i1, i2)
    coup2(4) = cpl_CFpSf_R(i_in, i1, i2)
    coup2(6) = cpl_FpFW(i1,i1) 
    coup2(7) = Conjg( cpl_FNSf_R(i1, i_out, i2) )
    coup2(8) = Conjg( cpl_FNSf_L(i1, i_out, i2) )
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsWsf8, n_WSf8, resC, check)
    gNffpSum(i1,i1,Isum) = 2._dp * Real(resC,dp)
    Contribution(i1,i1,Isum) = &
        & 'W Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do

  End If

  !-------------------------------------------
  ! S+ boson  s-channel - S+ boson  s-channel
  !-------------------------------------------
  Do i1 = 1,n_Spm
   Boson4(1) = mSpm(i1)
   Boson4(2) = gSpm(i1)
   coup2(1) = Conjg(Cpl_SmpCN_L(i1,i_in,i_out))
   coup2(2) = Conjg(Cpl_SmpCN_R(i1,i_in,i_out))
   Do i2 = i1+1,n_Spm
    Isum = Isum + 1
    Boson4(3) = mSpm(i2)
    Boson4(4) = gSpm(i2)
    coup2(3) = Cpl_SmpCN_L(i2,i_in,i_out)
    coup2(4) = Cpl_SmpCN_R(i2,i_in,i_out)

    If (GenerationMixing) Then
     Do i3=1,n_f
      Do i4=1,n_f
       mass(2) = mN(i_out)
       mass(3) = -mfp(i4)
       mass(4) = mf(i3)

       coup2(5) = cpl_SmpFFp_L(i1,i3,i4)
       coup2(6) = cpl_SmpFFp_R(i1,i3,i4)
       coup2(7) = Conjg(cpl_SmpFFp_R(i2,i3,i4) )
       coup2(8) = Conjg(cpl_SmpFFp_L(i2,i3,i4) )
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsSpmC4, n_SpmC4, resC, check)
       gNffpSum(i3,i4,Isum) = 2._dp * Real(resC,dp)
       Contribution(i3,i4,Isum) = &
        & 'S^-_'//Bu(i1)//' S^-_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i4)
      End Do
     End Do

    Else
     Do i3=1,n_f
      mass(2) = mN(i_out)
      mass(3) = -mfp(i3)
      mass(4) = mf(i3)
      coup2(5) = cpl_SmpFFp_L(i1,i3,i3)
      coup2(6) = cpl_SmpFFp_R(i1,i3,i3)
      coup2(7) = Conjg(cpl_SmpFFp_R(i2,i3,i3) )
      coup2(8) = Conjg(cpl_SmpFFp_L(i2,i3,i3) )
      Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsSpmC4, n_SpmC4, resC, check)
      gNffpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
      Contribution(i3,i3,Isum) = &
        & 'S^-_'//Bu(i1)//' S^-_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
     End Do

    End If
   End Do
  End Do

  !---------------------------------------------
  ! S+ boson  s-channel - Sfp_{xyz}  t-channel
  !---------------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Conjg(Cpl_SmpCN_R(i1,i_in,i_out))
    coup2(2) = Conjg(Cpl_SmpCN_L(i1,i_in,i_out))
    Do i2 = 1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)
     Do i4 =1,n_f
      coup2(3) = Conjg(cpl_CFSfp_L(i_in,i4,i2))
      coup2(4) = Conjg(cpl_CFSfp_R(i_in,i4,i2))
      Do i3=1,n_f
       mass(2) = mf(i4)
       mass(3) = -mfp(i3)
       mass(4) = mN(i_out)
       coup2(5) = cpl_SmpFFp_L(i1,i4,i3)
       coup2(6) = cpl_SmpFFp_R(i1,i4,i3)
       coup2(7) = cpl_FpNSfp_R(i3,i_out, i2)
       coup2(8) = cpl_FpNSfp_L(i3,i_out, i2)
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSpmSf8, n_SpmSf8, resC, check)
       gNffpSum(i4,i3,Isum) = - 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'S^-_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Conjg(Cpl_SmpCN_R(i1,i_in,i_out))
    coup2(2) = Conjg(Cpl_SmpCN_L(i1,i_in,i_out))
    Do i2 = 1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)

     If (n_Sfp.Le.3) Then
      i3 = i2
     Else
      i3 = (i2+1)/2
     End If

     mass(2) = mf(i3)
     mass(3) = -mfp(i3)
     mass(4) = mN(i_out)

     coup2(3) = Conjg(cpl_CFSfp_L(i_in,i3,i2))
     coup2(4) = Conjg(cpl_CFSfp_R(i_in,i3,i2))
     coup2(5) = cpl_SmpFFp_L(i1,i3,i3)
     coup2(6) = cpl_SmpFFp_R(i1,i3,i3)
     coup2(7) = cpl_FpNSfp_R(i3,i_out, i2)
     coup2(8) = cpl_FpNSfp_L(i3,i_out, i2)
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsSpmSf8, n_SpmSf8, resC, check)
     gNffpSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
        & 'S^-_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
    End Do
   End Do

  End If

  !---------------------------------------------
  ! S+ boson  s-channel - Sf_{xyz}  u-channel
  !---------------------------------------------
  If (GenerationMixing) Then
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Conjg(Cpl_SmpCN_R(i1,i_in,i_out) )
    coup2(2) = Conjg(Cpl_SmpCN_L(i1,i_in,i_out) )
    Do i2 = 1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i4 =1,n_f
      coup2(3) = cpl_CFpSf_L(i_in,i4,i2)
      coup2(4) = cpl_CFpSf_R(i_in,i4,i2)
      Do i3=1,n_f
       mass(2) = mfp(i4)
       mass(3) = -mf(i3)
       mass(4) = mN(i_out)

       coup2(5) = cpl_SmpFFp_L(i1,i3,i4)
       coup2(6) = cpl_SmpFFp_R(i1,i3,i4)
       coup2(7) = Conjg(cpl_FNSf_R(i3, i_out, i2))
       coup2(8) = Conjg(cpl_FNSf_L(i3, i_out, i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                            &, IntegralsSpmSf8, n_SpmSf8, resC, check)
       gNffpSum(i3,i4,Isum) = - 2._dp * Real(resC,dp)      
       Contribution(i3,i4,Isum) = &
        & 'S^-_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i4)
      End Do
     End Do
    End Do
   End Do

  Else
   Do i1 = 1,n_Spm
    Boson4(1) = mSpm(i1)
    Boson4(2) = gSpm(i1)
    coup2(1) = Conjg(Cpl_SmpCN_R(i1,i_in,i_out) )
    coup2(2) = Conjg(Cpl_SmpCN_L(i1,i_in,i_out) )
    Do i2 = 1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)

     i3 = (i2+1)/2

     mass(2) = mfp(i3)
     mass(3) = -mf(i3)
     mass(4) = mN(i_out)

     coup2(3) = cpl_CFpSf_L(i_in,i3,i2)
     coup2(4) = cpl_CFpSf_R(i_in,i3,i2)
     coup2(5) = cpl_SmpFFp_L(i1,i3,i3)
     coup2(6) = cpl_SmpFFp_R(i1,i3,i3)
     coup2(7) = Conjg(cpl_FNSf_R(i3, i_out, i2))
     coup2(8) = Conjg(cpl_FNSf_L(i3, i_out, i2))
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsSpmSf8, n_SpmSf8, resC, check)
     gNffpSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)      
     Contribution(i3,i3,Isum) = &
      & 'S^-_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
    End Do
   End Do

  End If

  !-----------------------------------------------
  ! Sfp_{xyz}  t-channel - Sfp_{xyz}  t-channel
  !-----------------------------------------------
  If (GenerationMixing) Then
   Do i1=1,n_Sfp-1
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2=i1+1,n_Sfp
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)

     Do i4 = 1,n_f
      coup2(1) = cpl_CFSfp_L(i_in,i4,i1)
      coup2(2) = cpl_CFSfp_R(i_in,i4,i1)
      coup2(3) = Conjg(cpl_CFSfp_R(i_in,i4,i2))
      coup2(4) = Conjg(cpl_CFSfp_L(i_in,i4,i2))
      Do i3=1,n_f
       mass(2) = mf(i4)
       mass(3) = -mfp(i3)
       mass(4) = mN(i_out)
       coup2(5) = Conjg(cpl_FpNSfp_R(i3,i_out,i1))
       coup2(6) = Conjg(cpl_FpNSfp_L(i3,i_out,i1))
       coup2(7) = cpl_FpNSfp_L(i3,i_out,i2)
       coup2(8) = cpl_FpNSfp_R(i3,i_out,i2)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsSfC4, n_SfC4, resC, check)
       gNffpSum(i4,i3,Isum) = 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else If (n_Sfp.Gt.3) Then ! no  contribution in case of sneutrinos
   Do i1=1,n_Sfp-1,2
    i2 = i1+1
    Isum = Isum + 1

    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)

    i3 = (i1+1)/2

    mass(2) = mf(i3)
    mass(3) = -mfp(i3)
    mass(4) = mN(i_out)

    coup2(1) = cpl_CFSfp_L(i_in,i3,i1)
    coup2(2) = cpl_CFSfp_R(i_in,i3,i1)
    coup2(3) = Conjg(cpl_CFSfp_R(i_in,i3,i2))
    coup2(4) = Conjg(cpl_CFSfp_L(i_in,i3,i2))
    coup2(5) = Conjg(cpl_FpNSfp_R(i3,i_out,i1))
    coup2(6) = Conjg(cpl_FpNSfp_L(i3,i_out,i1))
    coup2(7) = cpl_FpNSfp_L(i3,i_out,i2)
    coup2(8) = cpl_FpNSfp_R(i3,i_out,i2)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsSfC4, n_SfC4, resC, check)
    gNffpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    Contribution(i3,i3,Isum) = &
     & 'Sfp_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
   End Do
  End If

  !------------------------------------------------
  ! Sfp_{xyz}  t-channel - Sf_{xyz} u-channel
  !------------------------------------------------
  If (GenerationMixing) Then
   Do i1= 1,n_Sfp
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2 =1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i4 = 1,n_f
      coup2(1) = cpl_CFSfp_L(i_in,i4,i1)
      coup2(2) = cpl_CFSfp_R(i_in,i4,i1)
      Do i3 =1, n_f
       mass(2) = mfp(i3)
       mass(3) = -mN(i_out)
       mass(4) = mf(i4)
       coup2(3) = cpl_CFpSf_L(i_in, i3, i2)
       coup2(4) = cpl_CFpSf_R(i_in, i3, i2)
       coup2(5) = Conjg(cpl_FpNSfp_R(i3,i_out,i1))
       coup2(6) = Conjg(cpl_FpNSfp_L(i3,i_out,i1))
       coup2(7) = Conjg(cpl_FNSf_R(i4, i_out, i2))
       coup2(8) = Conjg(cpl_FNSf_L(i4, i_out, i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSf8, n_Sf8, resC, check)
       gNffpSum(i4,i3,Isum) = - 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else 
   Do i1= 1,n_Sfp
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    If (n_Sfp.Le.3) Then
     i3 = i1
    Else
     i3 = (i1+1)/2
    End If
    coup2(1) = cpl_CFSfp_L(i_in,i3,i1)
    coup2(2) = cpl_CFSfp_R(i_in,i3,i1)
    coup2(5) = Conjg(cpl_FpNSfp_R(i3,i_out,i1))
    coup2(6) = Conjg(cpl_FpNSfp_L(i3,i_out,i1))

    Do i2 =2*i3-1,2*i3
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)

     mass(2) = mf(i3)
     mass(3) = -mN(i_out)
     mass(4) = mfp(i3)

     coup2(3) = cpl_CFpSf_L(i_in, i3, i2)
     coup2(4) = cpl_CFpSf_R(i_in, i3, i2)
     coup2(7) = Conjg(cpl_FNSf_R(i3, i_out, i2))
     coup2(8) = Conjg(cpl_FNSf_L(i3, i_out, i2))
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsSf8, n_Sf8, resC, check)
     gNffpSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
    End Do
   End Do

  End If

  !-----------------------------------------------
  ! Sf_{xyz}  u-channel - Sf_{xyz}  u-channel
  !-----------------------------------------------
  If (GenerationMixing) Then

   Do i1=1,2*n_f-1
    Boson4(1) = mSf(i1)
    Boson4(2) = gSf(i1)
    Do i2=i1+1,2*n_f
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i4 = 1,n_f
      coup2(1) = Conjg(cpl_CFpSf_R(i_in, i4, i1))
      coup2(2) = Conjg(cpl_CFpSf_L(i_in, i4, i1))
      coup2(3) = cpl_CFpSf_L(i_in, i4, i2)
      coup2(4) = cpl_CFpSf_R(i_in, i4, i2)
      Do i3=1,n_f
       mass(2) = mfp(i4)
       mass(3) = - mN(i_out)
       mass(4) = mf(i3)
       coup2(5) = cpl_FNSf_L(i3, i_out, i1)
       coup2(6) = cpl_FNSf_R(i3, i_out, i1)
       coup2(7) = Conjg(cpl_FNSf_R(i3, i_out, i2))
       coup2(8) = Conjg(cpl_FNSf_L(i3, i_out, i2))
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsSfC4, n_SfC4, resC, check)
       gNffpSum(i3,i4,Isum) = 2._dp * Real(resC,dp)
       Contribution(i3,i4,Isum) = &
        & 'Sf_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i4)
      End Do
     End Do
    End Do
   End Do

  Else
   Do i1=1,2*n_f-1,2
    Isum = Isum + 1
    i2 = i1+1
    Boson4(1) = mSf(i1)
    Boson4(2) = gSf(i1)
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)

    i3 = (i1+1)/2

    mass(2) = mfp(i3)
    mass(3) = - mN(i_out)
    mass(4) = mf(i3)

    coup2(1) = Conjg(cpl_CFpSf_R(i_in, i3, i1))
    coup2(2) = Conjg(cpl_CFpSf_L(i_in, i3, i1))
    coup2(3) = cpl_CFpSf_L(i_in, i3, i2)
    coup2(4) = cpl_CFpSf_R(i_in, i3, i2)
    coup2(5) = cpl_FNSf_L(i3, i_out, i1)
    coup2(6) = cpl_FNSf_R(i3, i_out, i1)
    coup2(7) = Conjg(cpl_FNSf_R(i3, i_out, i2))
    coup2(8) = Conjg(cpl_FNSf_L(i3, i_out, i2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSfC4, n_SfC4, resC, check)
    gNffpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    Contribution(i3,i3,Isum) = &
      & 'Sf_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
   End Do
  End If


  !----------
  ! Summing
  !----------
  gNffp = 0._dp
  Do i1=1,n_f
   If (GenerationMixing) Then
    Do i2=1,n_f
     gNffp(i1,i2) = Sum( gNffpSum(i1,i2,1:Isum) )
     If (gNffp(i1,i2).Lt.0._dp) Then
      Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi_^-'//Bu(i_in)//' -> Chi_'//Bu(i_out)//') < 0 :' &
      & ,i1,i2,gNffp(i1,i2)
      Write (ErrCan,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gNffpSum(i1,i2,i3).Ne.0._dp) &
        &      Write (ErrCan,*) Contribution(i1,i2,i3),gNffpSum(i1,i2,i3)
      End Do
      gNffp(i1,i2) = 0._dp
     End If
    End Do

   Else
    gNffp(i1,i1) = Sum( gNffpSum(i1,i1,1:Isum) )
    If (gNffp(i1,i1).Lt.0._dp) Then
     Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//') < 0 :' &
      & ,i1,i1,gNffp(i1,i1)
     Write (ErrCan,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gNffpSum(i1,i1,i3).Ne.0._dp) &
        & Write (ErrCan,*) Contribution(i1,i1,i3),gNffpSum(i1,i1,i3)
     End Do
     gNffp(i1,i1) = 0._dp
    End If
   End If
  End Do

  gNffp = fac * gNffp

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gNffpSum = gNffpSum * fac

   If (GenerationMixing) Then
    Do i1=1,n_f
     Do i2=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//') :' &
      & ,i1,i2,gNffp(i1,i2)
     Write (n_out,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gNffpSum(i1,i2,i3).Ne.0._dp) &
        &       Write (n_out,*) Contribution(i1,i2,i3),gNffpSum(i1,i2,i3)
      End Do
      Write (n_out,*) ' '
     End Do
    End Do

   Else
    Do i1=1,n_f
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Chi_'//Bu(i_out)//') :' &
      & ,i1,i1,gNffp(i1,i1)
     Write (n_out,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gNffpSum(i1,i1,i3).Ne.0._dp) &
        &      Write (n_out,*) Contribution(i1,i1,i3),gNffpSum(i1,i1,i3)
     End Do
     Write (n_out,*) ' '
    End Do

   End If
  End If

  Deallocate( gNffpSum, Contribution)

  Iname = Iname - 1

 End Subroutine ChimToChi0ffp
 
 Subroutine ChimToChim2Chi0(i_in, i1, i2, i3, mC, mN, mZ, gZ, cpl_CCZ_L    &
    & , cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, mW, gW, cpl_CNW_L, cpl_CNW_R, mSpm   &
    & , gSpm, cpl_SmpCN_L, cpl_SmpCN_R, mS0, gS0, cpl_CCS0_L, cpl_CCS0_R      &
    & , cpl_NNS0_L, cpl_NNS0_R, mP0, gP0, cpl_CCP0_L, cpl_CCP0_R, cpl_NNP0_L  &
    & , cpl_NNP0_R, IntegralsZ4, n_Z4, IntegralsW4, n_W4, IntegralsS04, n_S04 &
    & , IntegralsZW8, n_ZW8, IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08     &
    & , IntegralsW8, n_W8, IntegralsWSpm4, n_WSpm4, IntegralsWSpm8, n_WSpm8   &
    & , IntegralsS0P04, n_S0P04, IntegralsS0P08, n_S0P08                      &
    & , deltaM, epsI, check, fac, gCNN, WriteContribution, n_out)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a neutralino into another neutralino +
 ! two charginos
 ! input:
 !    i_in ....... index of the decaying neutralino, in case n_in < 0,
 !                 the decay widhts of all neutralinos will be calculated
 !    mN ......... neutralino masses
 !    mC ......... chargino masses
 !
 ! output:
 ! written by Werner Porod, 26.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, i1, i2, i3
  Integer, Intent(inout) :: n_Z4, n_W4, n_S04, n_ZS04, n_S0P04, n_ZW8, n_ZS08 &
     & , n_W8, n_S0P08, n_WSpm4, n_WSpm8

  Real(dp), Intent(in) :: mN(:), mZ, gZ, mC(:), mW, gW, mSpm(:), gSpm(:)   &
     & , mS0(:), gS0(:), mP0(:), gP0(:), epsI, deltaM, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsW4(:,:)   &
     & , IntegralsS04(:,:)

  Complex(dp), Intent(in) :: cpl_NNZ_L(:,:), cpl_NNZ_R(:,:), cpl_CCZ_L(:,:)  &
     & , cpl_CCZ_R(:,:), cpl_CNW_L(:,:), cpl_CNW_R(:,:), cpl_SmpCN_L(:,:,:)  &
     & , cpl_SmpCN_R(:,:,:), cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:)            &
     & , cpl_CCS0_L(:,:,:), cpl_CCS0_R(:,:,:), cpl_NNP0_L(:,:,:)             &
     & , cpl_NNP0_R(:,:,:), cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:)

  Complex(dp), Intent(inout) :: IntegralsZS04(:,:), IntegralsS0P04(:,:) &
      & , IntegralsZW8(:,:), IntegralsZS08(:,:), IntegralsS0P08(:,:)    &
      & , IntegralsW8(:,:), IntegralsWSpm4(:,:), IntegralsWSpm8(:,:)

  Real(dp), Intent(out) :: gCNN

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution, check

  Integer :: Isum, n1, n2, n_S0, n_P0, n_Spm
  Real(dp) :: Boson2(2), Boson4(4), mass(4), resR, resRa
  Complex(dp) :: coup1(4), coup2(8), resC, resCa

  Real(dp), Allocatable :: gCNNSum(:)
  Character(len=12), Allocatable :: Contribution(:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChimToChim2Chi0'

  mass(1) = mC(i_in)
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )
  n_Spm = Size( mSpm )

  Isum = (3 + 2 * n_Spm + n_S0 + n_P0)**2
  Allocate( gCNNSum(Isum) )
  Allocate( Contribution(Isum) )
   
  gCNNSum = 0._dp
  Contribution = ' '

  !-----
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ
  mass(2) = mC(i1)
  mass(3) = -mN(i3)
  mass(4) = mN(i2) 
  coup1(1) = Cpl_CCZ_L(i1,i_in)
  coup1(2) = Cpl_CCZ_R(i1,i_in)
  coup1(3) = Cpl_NNZ_L(i2,i3)
  coup1(4) = Cpl_NNZ_R(i2,i3)
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resR, check)
  gCNNSum(Isum) = resR
  Contribution(Isum) = 'Z'

  !--------------
  ! W W
  !--------------
  Isum = Isum + 1
  Boson2(1) = mW
  Boson2(2) = gW
  !---------------
  ! t-channel
  !---------------
  coup1(1) = Conjg(Cpl_CNW_L(i_in,i3))
  coup1(2) = Conjg(Cpl_CNW_R(i_in,i3))
  coup1(3) = Cpl_CNW_L(i1,i2)
  coup1(4) = Cpl_CNW_R(i1,i2)
  mass(2) = mN(i3)
  mass(3) = -mN(i2)
  mass(4) = mC(i1)                 
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                       &, n_W4, resR, check)
  !---------------
  ! u-channel
  !---------------
  coup1(1) = Conjg(Cpl_CNW_L(i_in,i2))
  coup1(2) = Conjg(Cpl_CNW_R(i_in,i2))
  coup1(3) = Cpl_CNW_L(i1,i3)
  coup1(4) = Cpl_CNW_R(i1,i3)
  mass(2) = mN(i2)
  mass(3) = -mN(i3)
  mass(4) = mC(i1)                 
  Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsW4 &
                       &, n_W4, resRa, check)
  gCNNSum(Isum) = resR+resRa
  Contribution(Isum) = 'W'

  !---------------------------------------------
  ! decays via S^\pm, without interference part 
  !---------------------------------------------
  Do n1=1,n_Spm
   Isum = Isum + 1
   Boson2(1) = mSpm(n1)
   Boson2(2) = gSpm(n1)
   !------------------
   ! t-channel
   !------------------
   mass(2) = mN(i3)
   mass(3) = -mN(i2)
   mass(4) = mC(i1)
   coup1(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i3))
   coup1(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i3))
   coup1(3) = Cpl_SmpCN_L(n1,i1,i2)
   coup1(4) = Cpl_SmpCN_R(n1,i1,i2)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04 ,n_S04, resRa, check)
   !------------------
   ! u-channel
   !------------------
   mass(2) = mN(i2)
   mass(3) = -mN(i3)
   mass(4) = mC(i1)
   coup1(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i2))
   coup1(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i2))
   coup1(3) = Cpl_SmpCN_L(n1,i1,i3)
   coup1(4) = Cpl_SmpCN_R(n1,i1,i3)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04 ,n_S04, resR, check)
   gCNNSum(Isum) = resR+resRa
   Contribution(Isum) = 'S-_'//Bu(n1)
  End Do  ! n1 S^\pm

  !-------------------
  ! Scalar S_0
  !-------------------
  mass(2) = mC(i1)
  mass(3) = -mN(i3)
  mass(4) = mN(i2)
  Do n1=1,n_S0
   Isum = Isum + 1
   Boson2(1) = mS0(n1)
   Boson2(2) = gS0(n1)
   coup1(1) = cpl_CCS0_L(i1,i_in,n1)
   coup1(2) = cpl_CCS0_R(i1,i_in,n1)
   coup1(3) = cpl_NNS0_L(i3,i2,n1)
   coup1(4) = cpl_NNS0_R(i3,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   gCNNSum(Isum) = resR
   Contribution(Isum) = 'S0_'//Bu(n1)
  End Do ! n1 ,  S^0_n1

  !-------------------
  ! Pseudoscalar P_0
  !-------------------
  mass(2) = mC(i1)
  mass(3) = -mN(i3)
  mass(4) = mN(i2)
  Do n1=1,n_P0
   Isum = Isum + 1
   Boson2(1) = mP0(n1)
   Boson2(2) = gP0(n1)
   coup1(1) = cpl_CCP0_L(i1,i_in,n1)
   coup1(2) = cpl_CCP0_R(i1,i_in,n1)
   coup1(3) = cpl_NNP0_L(i2,i3,n1)
   coup1(4) = cpl_NNP0_R(i2,i3,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   gCNNSum(Isum) = resR
   Contribution(Isum) = 'P0_'//Bu(n1)
  End Do ! n1 ,  P^0_n1

  !-------
  ! Z - W
  !-------
  Boson4(1) = mZ
  Boson4(2) = gZ
  Boson4(3) = mW
  Boson4(4) = gW
  Isum = Isum + 1
  !-----------
  ! t-s channel
  !-----------
  mass(2) = mN(i3)
  mass(3) = -mN(i2)
  mass(4) = mC(i1)
  coup2(1) = cpl_CCZ_L(i1, i_in)
  coup2(2) = cpl_CCZ_R(i1, i_in)
  coup2(3) = cpl_CNW_L(i_in, i3)
  coup2(4) = cpl_CNW_R(i_in, i3)
  coup2(5) = cpl_NNZ_R(i2, i3)
  coup2(6) = cpl_NNZ_L(i2, i3)
  coup2(7) = Conjg(cpl_CNW_L(i1, i2))
  coup2(8) = Conjg(cpl_CNW_R(i1, i2))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZW8, n_ZW8, resC, check)
  !-----------
  ! u-s channel
  !-----------
  mass(2) = mN(i2)
  mass(3) = -mN(i3)
  mass(4) = mC(i1)
  coup2(1) = cpl_CCZ_L(i1, i_in)
  coup2(2) = cpl_CCZ_R(i1, i_in)
  coup2(3) = cpl_CNW_L(i_in, i2)
  coup2(4) = cpl_CNW_R(i_in, i2)
  coup2(5) = cpl_NNZ_L(i2, i3)
  coup2(6) = cpl_NNZ_R(i2, i3)
  coup2(7) = Conjg(cpl_CNW_L(i1, i3))
  coup2(8) = Conjg(cpl_CNW_R(i1, i3))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZW8, n_ZW8, resCa, check)
  gCNNSum(Isum) = 2._dp * Real(resC-resCa,dp)
  Contribution(Isum) = 'Z W'

  !-------------------------------
  ! Z boson - S^0 boson  s-channel
  !-------------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ

  mass(2) = mC(i1)
  mass(3) = -mN(i3)
  mass(4) = mN(i2)
  coup2(1) = Cpl_CCZ_L(i1,i_in)
  coup2(2) = Cpl_CCZ_R(i1,i_in)
  coup2(5) = Cpl_NNZ_L(i2,i3)
  coup2(6) = Cpl_NNZ_R(i2,i3)
  Do n1=1,n_S0 
   Isum = Isum + 1
   Boson4(3) = mS0(n1)
   Boson4(4) = gS0(n1)
   coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n1))
   coup2(7) = Conjg(cpl_NNS0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_NNS0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resC, check)
   gCNNSum(Isum) = 2._dp * Real(resC,dp)
   Contribution(Isum) = 'Z S0_'//Bu(n1)
  End Do  ! n1 n_S0

  !-------------------------------
  ! Z boson - P^0 boson  s-channel
  !-------------------------------
  Do n1=1,n_P0 
   Isum = Isum + 1
   Boson4(3) = mP0(n1)
   Boson4(4) = gP0(n1)
   coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n1))
   coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsZS04, n_ZS04, resC, check)
   gCNNSum(Isum) = 2._dp * Real(resC,dp)
   Contribution(Isum) = 'Z P0_'//Bu(n1)
  End Do  ! n1 n_P0

  !---------------------------------
  ! Z boson - S+ boson  t-u channel
  !---------------------------------
  coup2(1) = Cpl_CCZ_L(i1, i_in)
  coup2(2) = Cpl_CCZ_R(i1, i_in)
  Do n1=1,n_Spm  
   Isum = Isum + 1
   Boson4(3) = mSpm(n1)
   Boson4(4) = gSpm(n1)
   !-------------
   ! s-t channel
   !-------------
   mass(2) = mN(i3)
   mass(3) = -mN(i2)
   mass(4) = mC(i1)
   coup2(3) = Cpl_SmpCN_L(n1, i_in, i3)
   coup2(4) = Cpl_SmpCN_R(n1, i_in, i3)
   coup2(5) = cpl_NNZ_R(i2, i3)
   coup2(6) = cpl_NNZ_L(i2, i3)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1, i1, i2))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1, i1, i2))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsZS08, n_ZS08, resC, check)
   !-------------
   ! s-u channel
   !-------------
   mass(2) = mN(i2)
   mass(3) = -mN(i3)
   mass(4) = mC(i1)
   coup2(3) = Cpl_SmpCN_L(n1, i_in, i2)
   coup2(4) = Cpl_SmpCN_R(n1, i_in, i2)
   coup2(5) = cpl_NNZ_L(i2, i3)
   coup2(6) = cpl_NNZ_R(i2, i3)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1, i1, i3))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1, i1, i3))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsZS08, n_ZS08, resCa, check)
   gCNNSum(Isum) = 2._dp * Real(resC-resCa,dp)
   Contribution(Isum) = 'Z S-_'//Bu(n1)
  End Do ! n1 S^+m

  !-------------------------------
  ! W boson - W-boson
  !-------------------------------
  Boson4(1) = mW
  Boson4(2) = gW
  Boson4(3) = mW
  Boson4(4) = gW
  Isum = Isum + 1

  mass(2) = mN(i2)
  mass(3) = -mC(i1)
  mass(4) = mN(i3)
  coup2(1) = Conjg(Cpl_CNW_L(i_in,i3))
  coup2(2) = Conjg(Cpl_CNW_R(i_in,i3))
  coup2(3) = Cpl_CNW_L(i_in,i2)
  coup2(4) = Cpl_CNW_R(i_in,i2)
  coup2(5) = Cpl_CNW_R(i1,i2)
  coup2(6) = Cpl_CNW_L(i1,i2)
  coup2(7) = Conjg(Cpl_CNW_R(i1,i3))
  coup2(8) = Conjg(Cpl_CNW_L(i1,i3))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsW8, n_W8, resC, check)
  gCNNSum(Isum) = -2._dp * Real(resC,dp)
  Contribution(Isum) = 'W W st'

  !---------------------------------
  ! W boson - S+ boson  same-channel
  !---------------------------------
  Do n1=1,n_Spm 
   Isum = Isum + 1
   Boson4(3) = mSpm(n1)
   Boson4(4) = gSpm(n1)
   !-----------------
   ! t-channel
   !-----------------
   mass(2) = mN(i3)
   mass(3) = -mN(i2)
   mass(4) = mC(i1)
   coup2(1) = Conjg(Cpl_CNW_L(i_in,i3))
   coup2(2) = Conjg(Cpl_CNW_R(i_in,i3))
   coup2(3) = Cpl_SmpCN_L(n1,i_in,i3)
   coup2(4) = Cpl_SmpCN_R(n1,i_in,i3)
   coup2(5) = Cpl_CNW_L(i1,i2)
   coup2(6) = Cpl_CNW_R(i1,i2)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1,i1,i2))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1,i1,i2))
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm4, n_WSpm4, resC, check)
   !-----------------
   ! u-channel
   !-----------------
   mass(2) = mN(i2)
   mass(3) = -mN(i3)
   mass(4) = mC(i1)
   coup2(1) = Conjg(Cpl_CNW_L(i_in,i2))
   coup2(2) = Conjg(Cpl_CNW_R(i_in,i2))
   coup2(3) = Cpl_SmpCN_L(n1,i_in,i2)
   coup2(4) = Cpl_SmpCN_R(n1,i_in,i2)
   coup2(5) = Cpl_CNW_L(i1,i3)
   coup2(6) = Cpl_CNW_R(i1,i3)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1,i1,i3))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1,i1,i3))
   Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm4, n_WSpm4, resCa, check)
   gCNNSum(Isum) = 2._dp * Real(resC+resCa,dp)
   Contribution(Isum) = 'W S-_'//Bu(n1)//' t-t'
  End Do  ! n1 Spm(n1)

  !---------------------------------
  ! W boson - S+ boson  t-u channel
  !---------------------------------
  Do n1=1,n_Spm  
   Isum = Isum + 1
   Boson4(3) = mSpm(n1)
   Boson4(4) = gSpm(n1)
   !-------------
   ! t-u channel
   !-------------
   mass(2) = mN(i2)
   mass(3) = -mC(i1)
   mass(4) = mN(i3)
   coup2(1) = Conjg(Cpl_CNW_L(i_in, i3))
   coup2(2) = Conjg(Cpl_CNW_R(i_in, i3))
   coup2(3) = Cpl_SmpCN_L(n1,i_in,i2)
   coup2(4) = Cpl_SmpCN_R(n1,i_in,i2)
   coup2(5) = Cpl_CNW_R(i1,i2)
   coup2(6) = Cpl_CNW_L(i1,i2)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1,i1,i3))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1,i1,i3))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resC, check)
   !-------------
   ! u-t channel
   !-------------
   mass(2) = mN(i3)
   mass(3) = -mC(i1)
   mass(4) = mN(i2)
   coup2(1) = Conjg(Cpl_CNW_L(i_in, i2))
   coup2(2) = Conjg(Cpl_CNW_R(i_in, i2))
   coup2(3) = Cpl_SmpCN_L(n1,i_in,i3)
   coup2(4) = Cpl_SmpCN_R(n1,i_in,i3)
   coup2(5) = Cpl_CNW_R(i1,i3)
   coup2(6) = Cpl_CNW_L(i1,i3)
   coup2(7) = Conjg(Cpl_SmpCN_R(n1,i1,i2))
   coup2(8) = Conjg(Cpl_SmpCN_L(n1,i1,i2))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resCa, check)
   gCNNSum(Isum) = 2._dp * Real(resC+resCa,dp)
   Contribution(Isum) = 'W S-_'//Bu(n1)//' s-t'
  End Do ! n1 S^+m

  !---------------------------------
  ! W boson - S0 boson  t-s channel
  !---------------------------------
  Do n1=1,n_S0  
   Isum = Isum + 1
   Boson4(3) = mS0(n1)
   Boson4(4) = gS0(n1)
   !-------------
   ! t-s channel
   !-------------
   mass(2) = mC(i1)
   mass(3) = -mN(i2)
   mass(4) = mN(i3)
   coup2(1) = Conjg(Cpl_CNW_L(i_in, i3))
   coup2(2) = Conjg(Cpl_CNW_R(i_in, i3))
   coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n1))
   coup2(5) = Cpl_CNW_L(i1,i2)
   coup2(6) = Cpl_CNW_R(i1,i2)
   coup2(7) = Conjg(cpl_NNS0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_NNS0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resC, check)
   !-------------
   ! u-s channel
   !-------------
   mass(2) = mC(i1)
   mass(3) = -mN(i3)
   mass(4) = mN(i2)
   coup2(1) = Conjg(Cpl_CNW_L(i_in, i2))
   coup2(2) = Conjg(Cpl_CNW_R(i_in, i2))
   coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n1))
   coup2(5) = Cpl_CNW_L(i1,i3)
   coup2(6) = Cpl_CNW_R(i1,i3)
   coup2(7) = Conjg(cpl_NNS0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_NNS0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resCa, check)
   gCNNSum(Isum) = -2._dp * Real(resC+resCa,dp)
   Contribution(Isum) = 'W S0_'//Bu(n1)//' s-t'
  End Do ! n1 S^0

  !---------------------------------
  ! W boson - P0 boson  t-s channel
  !---------------------------------
  Do n1=1,n_P0  
   Isum = Isum + 1
   Boson4(3) = mP0(n1)
   Boson4(4) = gP0(n1)
   !-------------
   ! t-s channel
   !-------------
   mass(2) = mC(i1)
   mass(3) = -mN(i2)
   mass(4) = mN(i3)
   coup2(1) = Conjg(Cpl_CNW_L(i_in, i3))
   coup2(2) = Conjg(Cpl_CNW_R(i_in, i3))
   coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n1))
   coup2(5) = Cpl_CNW_L(i1,i2)
   coup2(6) = Cpl_CNW_R(i1,i2)
   coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resC, check)
   !-------------
   ! u-s channel
   !-------------
   mass(2) = mC(i1)
   mass(3) = -mN(i3)
   mass(4) = mN(i2)
   coup2(1) = Conjg(Cpl_CNW_L(i_in, i2))
   coup2(2) = Conjg(Cpl_CNW_R(i_in, i2))
   coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n1))
   coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n1))
   coup2(5) = Cpl_CNW_L(i1,i3)
   coup2(6) = Cpl_CNW_R(i1,i3)
   coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n1))
   coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n1))
   Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                             &, IntegralsWSpm8, n_WSpm8, resCa, check)
   gCNNSum(Isum) = -2._dp * Real(resC+resCa,dp)
   Contribution(Isum) = 'W P0_'//Bu(n1)//' s-t'
  End Do ! n1 P^0

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  s-channel
  !-------------------------------------------------
  mass(2) = mC(i1)
  mass(3) = -mN(i3)
  mass(4) = mN(i2)
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   coup2(1) = cpl_CCS0_L(i1,i_in,n1)
   coup2(2) = cpl_CCS0_R(i1,i_in,n1)
   coup2(5) = cpl_NNS0_L(i2,i3,n1)
   coup2(6) = cpl_NNS0_R(i2,i3,n1)
   Do n2 = n1+1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n2))
    coup2(7) = Conjg(cpl_NNS0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    gCNNSum(Isum) = 2._dp * Real(resC,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//'S0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  mass(2) = mC(i1)
  mass(3) = -mN(i3)
  mass(4) = mN(i2)
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   coup2(1) = cpl_CCS0_L(i1,i_in,n1)
   coup2(2) = cpl_CCS0_R(i1,i_in,n1)
   coup2(5) = cpl_NNS0_L(i2,i3,n1)
   coup2(6) = cpl_NNS0_R(i2,i3,n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n2))
    coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    gCNNSum(Isum) = 2._dp * Real(resC,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//'P0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^+_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 = 1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mC(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i3))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i3))
    coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n2))
    coup2(5) = Cpl_SmpCN_L(n1,i1,i2)
    coup2(6) = Cpl_SmpCN_R(n1,i1,i2)
    coup2(7) = Conjg(cpl_NNS0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNS0_L(i2,i3,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mC(i1)
    mass(3) = -mN(i3)
    mass(4) = mN(i2)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i2))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i2))
    coup2(5) = Cpl_SmpCN_L(n1,i1,i3)
    coup2(6) = Cpl_SmpCN_R(n1,i1,i3)
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    gCNNSum(Isum) = - 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S-_'//Bu(n1)//'S0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  mass(2) = mC(i1)
  mass(3) = -mN(i3)
  mass(4) = mN(i2)
  Do n1 = 1,n_P0
   Boson4(1) = mP0(n1)
   Boson4(2) = gP0(n1)
   coup2(1) = cpl_CCP0_L(i1,i_in,n1)
   coup2(2) = cpl_CCP0_R(i1,i_in,n1)
   coup2(5) = cpl_NNP0_L(i2,i3,n1)
   coup2(6) = cpl_NNP0_R(i2,i3,n1)
   Do n2 = n1+1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n2))
    coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    gCNNSum(Isum) = 2._dp * Real(resC,dp)
    Contribution(Isum) = 'P0_'//Bu(n1)//'P0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------------
  ! P^0_i boson  s-channel - S^+_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mC(i1)
    mass(3) = -mN(i2)
    mass(4) = mN(i3)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i3))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i3))
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n2))
    coup2(5) = Cpl_SmpCN_L(n1,i1,i2)
    coup2(6) = Cpl_SmpCN_R(n1,i1,i2)
    coup2(7) = Conjg(cpl_NNP0_R(i2,i3,n2))
    coup2(8) = Conjg(cpl_NNP0_L(i2,i3,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-------------
    ! s-u channel
    !-------------
    mass(2) = mC(i1)
    mass(3) = -mN(i3)
    mass(4) = mN(i2)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i2))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i2))
    coup2(5) = Cpl_SmpCN_L(n1,i1,i3)
    coup2(6) = Cpl_SmpCN_R(n1,i1,i3)
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    gCNNSum(Isum) = - 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S-_'//Bu(n1)//'P0_'//Bu(n2)
   End Do
  End Do

  !-------------------------------------------
  ! S+ boson  t-channel - S+ boson  t-channel
  !-------------------------------------------
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 = n1+1,n_Spm
    Isum = Isum + 1
    Boson4(3) = mSpm(n2)
    Boson4(4) = gSpm(n2)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mN(i3)
    mass(3) = -mN(i2)
    mass(4) = mC(i1)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i3))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i3))
    coup2(3) = Cpl_SmpCN_L(n2,i_in,i3)
    coup2(4) = Cpl_SmpCN_R(n2,i_in,i3)
    coup2(5) = Cpl_SmpCN_L(n1,i1,i2)
    coup2(6) = Cpl_SmpCN_R(n1,i1,i2)
    coup2(7) = Conjg(Cpl_SmpCN_R(n2,i1,i2))
    coup2(8) = Conjg(Cpl_SmpCN_L(n2,i1,i2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! u-channel
    !-----------
    mass(2) = mN(i2)
    mass(3) = -mN(i3)
    mass(4) = mC(i1)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i2))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i2))
    coup2(3) = Cpl_SmpCN_L(n2,i_in,i2)
    coup2(4) = Cpl_SmpCN_R(n2,i_in,i2)
    coup2(5) = Cpl_SmpCN_L(n1,i1,i3)
    coup2(6) = Cpl_SmpCN_R(n1,i1,i3)
    coup2(7) = Conjg(Cpl_SmpCN_R(n2,i1,i3))
    coup2(8) = Conjg(Cpl_SmpCN_L(n2,i1,i3))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    gCNNSum(Isum) = 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S-_'//Bu(n1)//'S-_'//Bu(n2)//' t-t'
   End Do
  End Do

  !---------------------------------------------
  ! S+ boson  t-channel - S+ boson  u-channel
  !---------------------------------------------
  mass(2) = mN(i2)
  mass(3) = -mC(i1)
  mass(4) = mN(i3)
  Do n1 = 1,n_Spm
   Boson4(1) = mSpm(n1)
   Boson4(2) = gSpm(n1)
   Do n2 =1,n_Spm
    Isum = Isum + 1
    Boson4(3) = mSpm(n2)
    Boson4(4) = gSpm(n2)
    coup2(1) = Conjg(Cpl_SmpCN_R(n1,i_in,i3))
    coup2(2) = Conjg(Cpl_SmpCN_L(n1,i_in,i3))
    coup2(3) = Cpl_SmpCN_L(n2,i_in,i2)
    coup2(4) = Cpl_SmpCN_R(n2,i_in,i2)
    coup2(5) = Cpl_SmpCN_L(n1,i1,i2)
    coup2(6) = Cpl_SmpCN_R(n1,i1,i2)
    coup2(7) = Conjg(Cpl_SmpCN_R(n2,i1,i3))
    coup2(8) = Conjg(Cpl_SmpCN_L(n2,i1,i3))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    gCNNSum(Isum) = - 2._dp * Real(resC,dp) 
    Contribution(Isum) = 'S-_'//Bu(n1)//'S-_'//Bu(n2)//' t-u'
   End Do
  End Do

  !----------
  ! Summing
  !----------
  If (i2.Eq.i3) gCNNSum = 0.5_dp * gCNNSum ! identical particles
  gCNN = Sum( gCNNSum(1:Isum) )
  If (gCNN.Lt.0._dp) Then
   Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write (ErrCan,*) &
     & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//' Chi^-_'//Bu(i2)// &
     & ' Chi^-_'//Bu(i3)//') < 0 :',i_in,i1,i2,i3,gCNN
   Write (ErrCan,*) 'The different contributions are :'
   Do n1=1,Isum
    If (gCNNSum(n1).Ne.0._dp)  Write (ErrCan,*) Contribution(n1),gCNNSum(n1)
   End Do
   gCNN = 0._dp
  End If

  gCNN = fac * gCNN

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gCNNSum = gCNNSum * fac
 
   Write (n_out,*) &
     & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//' Chi^-_'//Bu(i2)// &
     & ' Chi^-_'//Bu(i3)//') :',i_in,i1,i2,i3
   Write (n_out,*) 'The different contributions are :'
   Do n1=1,Isum
    If (gCNNSum(n1).Ne.0._dp)  Write (n_out,*) Contribution(n1),gCNNSum(n1)
   End Do
  End If

  Deallocate( gCNNSum, Contribution )
   
  Iname = Iname - 1

 End Subroutine ChimToChim2Chi0
 
 Subroutine ChimToGffp(i_in, mC, mGlu, mf, mfp, mSf, gSf                 &
    & , cpl_FGSf_L, cpl_FGSf_R, cpl_CFpSf_L, cpl_CFpSf_R, mSfp, gSfp     &
    & , cpl_FpGSfp_L, cpl_FpGSfp_R, cpl_CFSfp_L, cpl_CFSfp_R             &
    & , IntegralsSf4, n_Sf4, IntegralsSfC4, n_SfC4, IntegralsSf8, n_Sf8  &
    & , deltaM, epsI, GenerationMixing, check, fac, gGffp                &
    & , WriteContribution, n_out)
 !--------------------------------------------------------------------------
 ! Calculates the decay of a Chargino to a neutralino + fermion pair
 ! Written by Werner Porod, 30.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !--------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in
  Integer, Intent(inout) :: n_Sf4, n_SfC4, n_Sf8
  
  Real(dp), Intent(in) :: mC(:), mGlu, mf(:), mfp(:), mSf(:), gSf(:) &
      & , mSfp(:), gSfp(:), deltaM, epsI, fac
  Real(dp), Intent(inout) :: IntegralsSf4(:,:)

  Complex(dp), Intent(in) :: cpl_FGSf_L(:,:), cpl_FGSf_R(:,:)         &
      & , cpl_CFpSf_L(:,:,:), cpl_CFpSf_R(:,:,:), cpl_FpGSfp_L(:,:)     &
      & , cpl_FpGSfp_R(:,:), cpl_CFSfp_L(:,:,:), cpl_CFSfp_R(:,:,:)
  Complex(dp), Intent(inout) :: IntegralsSfC4(:,:), IntegralsSf8(:,:)

  Logical, Intent(in) :: GenerationMixing, check

  Real(dp), Intent(out) :: gGffp(:,:)

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution

  Integer :: Isum, i1, i2, i3, i4
  Real(dp) :: Boson2(2), mass(4), resR, Boson4(4), gGffpSum(3,3,100)
  Complex(dp) :: coup1(4), resC, coup2(8)

  Character(len=20) :: Contribution(3,3,100)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChimToGffp'

  mass(1) = mC(i_in)


  gGffpSum = 0._dp
  Contribution = ' '

  !-------------------
  ! Sfp Sfp, diagonal
  !-------------------
  Isum = 0
    
  If (GenerationMixing) Then
   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)
    Do i3=1,3
     coup1(1) = cpl_CFSfp_L(i_in,i3,i2)
     coup1(2) = cpl_CFSfp_R(i_in,i3,i2)
     Do i1=1,3
      mass(2) = mf(i3)
      mass(3) = -mfp(i1)
      mass(4) = mGlu
      coup1(3) = Conjg(cpl_FpGSfp_R(i1,i2))
      coup1(4) = Conjg(cpl_FpGSfp_L(i1,i2))
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gGffpSum(i3,i1,Isum) = resR
      Contribution(i3,i1,Isum) = 'Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i1)
     End Do
    End Do
   End Do

  Else

   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSfp(i2)
    Boson2(2) = gSfp(i2)

    i1 = (i2+1)/2

    mass(2) = mf(i1)
    mass(3) = -mfp(i1)
    mass(4) = mGlu
    coup1(1) = cpl_CFSfp_L(i_in,i1,i2)
    coup1(2) = cpl_CFSfp_R(i_in,i1,i2)
    coup1(3) = Conjg(cpl_FpGSfp_R(i1,i2))
    coup1(4) = Conjg(cpl_FpGSfp_L(i1,i2))
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsSf4, n_Sf4, resR, check)
    gGffpSum(i1,i1,Isum) = resR
    Contribution(i1,i1,Isum) = 'Sfp_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do

  End If

  !-------------------------
  ! Sf Sf, diagonal
  !-------------------------
  If (GenerationMixing) Then
   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    Do i3=1,3
     coup1(1) = Conjg(cpl_CFpSf_R(i_in,i3,i2))
     coup1(2) = Conjg(cpl_CFpSf_L(i_in,i3,i2))
     Do i1=1,3
      mass(2) = mfp(i3)
      mass(3) = - mGlu
      mass(4) = mf(i1)
      coup1(3) = cpl_FGSf_L(i1,i2)
      coup1(4) = cpl_FGSf_R(i1,i2)
      Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                           &, IntegralsSf4, n_Sf4, resR, check)
      gGffpSum(i1,i3,Isum) = resR
      Contribution(i3,i1,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i1)
     End Do
    End Do
   End Do
  Else

   Do i2=1,6
    Isum = Isum + 1
    Boson2(1) = mSf(i2)
    Boson2(2) = gSf(i2)

    i1 = (i2+1)/2

    mass(2) = mfp(i1)
    mass(3) = - mGlu
    mass(4) = mf(i1)
    coup1(1) = Conjg(cpl_CFpSf_R(i_in,i1,i2))
    coup1(2) = Conjg(cpl_CFpSf_L(i_in,i1,i2))
    coup1(3) = cpl_FGSf_L(i1,i2)
    coup1(4) = cpl_FGSf_R(i1,i2)
    Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsSf4, n_Sf4, resR, check)
    gGffpSum(i1,i1,Isum) = resR
    Contribution(i1,i1,Isum) = 'Sf_'//Bu(i2)//' f_'//Bu(i1)//' fp_'//Bu(i1)
   End Do
  End If

  !-----------------------------------------------
  ! Sfp_{xyz}  t-channel - Sfp_{xyz}  t-channel
  !-----------------------------------------------
  If (GenerationMixing) Then
   Do i1=1,5
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2=i1+1,6
     Isum = Isum + 1
     Boson4(3) = mSfp(i2)
     Boson4(4) = gSfp(i2)

     Do i4 = 1,3
      coup2(1) = cpl_CFSfp_L(i_in,i4,i1)
      coup2(2) = cpl_CFSfp_R(i_in,i4,i1)
      coup2(3) = Conjg(cpl_CFSfp_R(i_in,i4,i2))
      coup2(4) = Conjg(cpl_CFSfp_L(i_in,i4,i2))
      Do i3=1,3
       mass(2) = mf(i4)
       mass(3) = -mfp(i3)
       mass(4) = mGlu
       coup2(5) = Conjg(cpl_FpGSfp_R(i3,i1))
       coup2(6) = Conjg(cpl_FpGSfp_L(i3,i1))
       coup2(7) = cpl_FpGSfp_L(i3,i2)
       coup2(8) = cpl_FpGSfp_R(i3,i2)
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsSfC4, n_SfC4, resC, check)
       gGffpSum(i4,i3,Isum) = 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else 
   Do i1=1,5,2
    i2 = i1+1
    Isum = Isum + 1

    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Boson4(3) = mSfp(i2)
    Boson4(4) = gSfp(i2)

    i3 = (i1+1)/2

    mass(2) = mf(i3)
    mass(3) = -mfp(i3)
    mass(4) = mGlu

    coup2(1) = cpl_CFSfp_L(i_in,i3,i1)
    coup2(2) = cpl_CFSfp_R(i_in,i3,i1)
    coup2(3) = Conjg(cpl_CFSfp_R(i_in,i3,i2))
    coup2(4) = Conjg(cpl_CFSfp_L(i_in,i3,i2))
    coup2(5) = Conjg(cpl_FpGSfp_R(i3,i1))
    coup2(6) = Conjg(cpl_FpGSfp_L(i3,i1))
    coup2(7) = cpl_FpGSfp_L(i3,i2)
    coup2(8) = cpl_FpGSfp_R(i3,i2)
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                          &, IntegralsSfC4, n_SfC4, resC, check)
    gGffpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    Contribution(i3,i3,Isum) = &
     & 'Sfp_'//Bu(i1)//' Sfp_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
   End Do
  End If

  !------------------------------------------------
  ! Sfp_{xyz}  t-channel - Sf_{xyz} u-channel
  !------------------------------------------------
  If (GenerationMixing) Then
   Do i1= 1,6
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    Do i2 =1,6
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i4 = 1,3
      coup2(1) = cpl_CFSfp_L(i_in,i4,i1)
      coup2(2) = cpl_CFSfp_R(i_in,i4,i1)
      Do i3 =1, 3
       mass(2) = mfp(i3)
       mass(3) = -mGlu
       mass(4) = mf(i4)
       coup2(3) = cpl_CFpSf_L(i_in, i3, i2)
       coup2(4) = cpl_CFpSf_R(i_in, i3, i2)
       coup2(5) = Conjg(cpl_FpGSfp_R(i3,i1))
       coup2(6) = Conjg(cpl_FpGSfp_L(i3,i1))
       coup2(7) = Conjg(cpl_FGSf_R(i4,  i2))
       coup2(8) = Conjg(cpl_FGSf_L(i4,  i2))
       Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSf8, n_Sf8, resC, check)
       gGffpSum(i4,i3,Isum) = - 2._dp * Real(resC,dp)
       Contribution(i4,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i4)//' fp_'//Bu(i3)
      End Do
     End Do
    End Do
   End Do

  Else 
   Do i1= 1,6
    Boson4(1) = mSfp(i1)
    Boson4(2) = gSfp(i1)
    i3 = (i1+1)/2

    coup2(1) = cpl_CFSfp_L(i_in,i3,i1)
    coup2(2) = cpl_CFSfp_R(i_in,i3,i1)
    coup2(5) = Conjg(cpl_FpGSfp_R(i3,i1))
    coup2(6) = Conjg(cpl_FpGSfp_L(i3,i1))

    Do i2 =2*i3-1,i3
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)

     mass(2) = mf(i3)
     mass(3) = -mGlu
     mass(4) = mfp(i3)

     coup2(3) = cpl_CFpSf_L(i_in, i3, i2)
     coup2(4) = cpl_CFpSf_R(i_in, i3, i2)
     coup2(7) = Conjg(cpl_FGSf_R(i3,  i2))
     coup2(8) = Conjg(cpl_FGSf_L(i3,  i2))
     Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsSf8, n_Sf8, resC, check)
     gGffpSum(i3,i3,Isum) = - 2._dp * Real(resC,dp)
     Contribution(i3,i3,Isum) = &
        & 'Sfp_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
    End Do
   End Do

  End If

  !-----------------------------------------------
  ! Sf_{xyz}  u-channel - Sf_{xyz}  u-channel
  !-----------------------------------------------
  If (GenerationMixing) Then

   Do i1=1,5
    Boson4(1) = mSf(i1)
    Boson4(2) = gSf(i1)
    Do i2=i1+1,6
     Isum = Isum + 1
     Boson4(3) = mSf(i2)
     Boson4(4) = gSf(i2)
     Do i4 = 1,3
      coup2(1) = Conjg(cpl_CFpSf_R(i_in, i4, i1))
      coup2(2) = Conjg(cpl_CFpSf_L(i_in, i4, i1))
      coup2(3) = cpl_CFpSf_L(i_in, i4, i2)
      coup2(4) = cpl_CFpSf_R(i_in, i4, i2)
      Do i3=1,3
       mass(2) = mfp(i4)
       mass(3) = - mGlu
       mass(4) = mf(i3)
       coup2(5) = cpl_FGSf_L(i3,  i1)
       coup2(6) = cpl_FGSf_R(i3,  i1)
       coup2(7) = Conjg(cpl_FGSf_R(i3,  i2))
       coup2(8) = Conjg(cpl_FGSf_L(i3,  i2))
       Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsSfC4, n_SfC4, resC, check)
       gGffpSum(i3,i4,Isum) = 2._dp * Real(resC,dp)
       Contribution(i3,i4,Isum) = &
        & 'Sf_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i4)
      End Do
     End Do
    End Do
   End Do

  Else
   Do i1=1,5,2
    Isum = Isum + 1
    i2 = i1+1
    Boson4(1) = mSf(i1)
    Boson4(2) = gSf(i1)
    Boson4(3) = mSf(i2)
    Boson4(4) = gSf(i2)

    i3 = (i1+1)/2

    mass(2) = mfp(i3)
    mass(3) = - mGlu
    mass(4) = mf(i3)

    coup2(1) = Conjg(cpl_CFpSf_R(i_in, i3, i1))
    coup2(2) = Conjg(cpl_CFpSf_L(i_in, i3, i1))
    coup2(3) = cpl_CFpSf_L(i_in, i3, i2)
    coup2(4) = cpl_CFpSf_R(i_in, i3, i2)
    coup2(5) = cpl_FGSf_L(i3,  i1)
    coup2(6) = cpl_FGSf_R(i3,  i1)
    coup2(7) = Conjg(cpl_FGSf_R(i3,  i2))
    coup2(8) = Conjg(cpl_FGSf_L(i3,  i2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsSfC4, n_SfC4, resC, check)
    gGffpSum(i3,i3,Isum) = 2._dp * Real(resC,dp)
    Contribution(i3,i3,Isum) = &
      & 'Sf_'//Bu(i1)//' Sf_'//Bu(i2)//' f_'//Bu(i3)//' fp_'//Bu(i3)
   End Do
  End If


  !----------
  ! Summing
  !----------
  gGffp = 0._dp
  Do i1=1,3
   If (GenerationMixing) Then
    Do i2=1,3
     gGffp(i1,i2) = Sum( gGffpSum(i1,i2,1:Isum) )
     If (gGffp(i1,i2).Lt.0._dp) Then
      Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi_^-'//Bu(i_in)//' -> Gluino) < 0 :' &
      & ,i1,i2,gGffp(i1,i2)
      Write (ErrCan,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gGffpSum(i1,i2,i3).Ne.0._dp) &
        &      Write (ErrCan,*) Contribution(i1,i2,i3),gGffpSum(i1,i2,i3)
      End Do
      gGffp(i1,i2) = 0._dp
     End If
    End Do

   Else
    gGffp(i1,i1) = Sum( gGffpSum(i1,i1,1:Isum) )
    If (gGffp(i1,i1).Lt.0._dp) Then
     Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
     Write (ErrCan,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Gluino) < 0 :' &
      & ,i1,i1,gGffp(i1,i1)
     Write (ErrCan,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gGffpSum(i1,i1,i3).Ne.0._dp) &
        & Write (ErrCan,*) Contribution(i1,i1,i3),gGffpSum(i1,i1,i3)
     End Do
     gGffp(i1,i1) = 0._dp
    End If
   End If
  End Do

  gGffp = fac * gGffp

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gGffpSum = gGffpSum * fac

   If (GenerationMixing) Then
    Do i1=1,3
     Do i2=1,3
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Gluino) :' &
      & ,i1,i2,gGffp(i1,i2)
     Write (n_out,*) 'The different contributions are :'
      Do i3=1,Isum
      If (gGffpSum(i1,i2,i3).Ne.0._dp) &
        &       Write (n_out,*) Contribution(i1,i2,i3),gGffpSum(i1,i2,i3)
      End Do
      Write (n_out,*) ' '
     End Do
    End Do

   Else
    Do i1=1,3
     Write (n_out,*) &
      & 'Gamma(Chi^-_'//Bu(i_in)//' -> Gluino) :' &
      & ,i1,i1,gGffp(i1,i1)
     Write (n_out,*) 'The different contributions are :'
     Do i3=1,Isum
      If (gGffpSum(i1,i1,i3).Ne.0._dp) &
        &      Write (n_out,*) Contribution(i1,i1,i3),gGffpSum(i1,i1,i3)
     End Do
     Write (n_out,*) ' '
    End Do

   End If
  End If

  Iname = Iname - 1

 End Subroutine ChimToGffp

 Subroutine ChimToChimChipChim(i_in, i1, i2, i3, mC, mZ, gZ, cpl_CCZ_L        &
    & , cpl_CCZ_R, mS0, gS0, cpl_CCS0_L, cpl_CCS0_R, mP0, gP0, cpl_CCP0_L     &
    & , cpl_CCP0_R, IntegralsZ4, n_Z4, IntegralsS04, n_S04, IntegralsZ8, n_Z8 &
    & , IntegralsZS04, n_ZS04, IntegralsZS08, n_ZS08, IntegralsS0P04, n_S0P04 &
    & , IntegralsS0P08, n_S0P08                                               &
    & , deltaM, epsI, check, fac, gCCC, WriteContribution, n_out)
 !------------------------------------------------------------------
 ! calculates all 3-body decays of a chargino into 3 charginos
 ! input:
 !    i_in ....... index of the decaying neutralino
 !                 the decay widhts of all neutralinos will be calculated
 !    mN ......... neutralino masses
 !
 ! output:
 ! written by Werner Porod, 26.06.2001
 ! 13.09.03: introducing new variable check: if .True. then only contributions
 !           are calculated if the intermediated states are off-shell
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, i1, i2, i3
  Integer, Intent(inout) :: n_Z4, n_S04, n_ZS04, n_S0P04, n_Z8, n_ZS08, n_S0P08

  Real(dp), Intent(in) :: mC(:), mZ, gZ, mS0(:), gS0(:), mP0(:), gP0(:) &
     & , epsI, deltaM, fac
  Real(dp), Intent(inout) :: IntegralsZ4(:,:), IntegralsS04(:,:)

  Complex(dp), Intent(in) :: Cpl_CCZ_L(:,:), Cpl_CCZ_R(:,:)             &
      & , cpl_CCS0_L(:,:,:), cpl_CCS0_R(:,:,:), cpl_CCP0_L(:,:,:)       &
      & , cpl_CCP0_R(:,:,:)

  Complex(dp), Intent(inout) :: IntegralsZS04(:,:), IntegralsS0P04(:,:) &
      & , IntegralsZ8(:,:), IntegralsZS08(:,:), IntegralsS0P08(:,:)

  Real(dp), Intent(out) :: gCCC

  Integer, Optional :: n_out
  Logical, Optional :: WriteContribution, check

  Integer :: Isum, n1, n2, n_S0, n_P0
  Real(dp) :: Boson2(2), Boson4(4), mass(4), resR, resRa
  Complex(dp) :: coup1(4), coup2(8), resC, resCa

  Real(dp), Allocatable :: gCCCSum(:)
  Character(len=15), Allocatable :: Contribution(:)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChimToChimChipChim'

  mass(1) = mC(i_in)
  n_S0 = Size( mS0 )
  n_P0 = Size( mP0 )

  Isum = 4 * (1 + n_S0 + n_P0)**2
  Allocate( gCCCSum(Isum) )
  Allocate( Contribution(Isum) )
   
  gCCCSum = 0._dp
  Contribution = ' '

  !-----
  ! Z Z
  !-----
  Isum = 1
  Boson2(1) = mZ
  Boson2(2) = gZ

   !--------------
   ! s-channel
   !--------------
   coup1(1) = Cpl_CCZ_L(i1,i_in)
   coup1(2) = Cpl_CCZ_R(i1,i_in)
   coup1(3) = Cpl_CCZ_L(i3,i2)
   coup1(4) = Cpl_CCZ_R(i3,i2)
   mass(2) = mC(i1)
   mass(3) = -mC(i2)
   mass(4) = mC(i3)       
   Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resR, check)
   !--------------
   ! t-channel
   !--------------
   coup1(1) = Cpl_CCZ_L(i3,i_in)
   coup1(2) = Cpl_CCZ_R(i3,i_in)
   coup1(3) = Cpl_CCZ_L(i1,i2)
   coup1(4) = Cpl_CCZ_R(i1,i2)
   mass(2) = mC(i3)
   mass(3) = -mC(i2)
   mass(4) = mC(i1)       
   Call IntegrateGaugeSS(Boson2, mass, coup1, deltaM, epsI, IntegralsZ4 &
                       &, n_Z4, resRa, check)
   gCCCSum(Isum) = resR + resRa
   Contribution(Isum) = 'Z'

  !-------------------
  ! Scalar S_0
  !-------------------
  Do n1=1,n_S0
   Isum = Isum + 1
   Boson2(1) = mS0(n1)
   Boson2(2) = gS0(n1)

   !--------------
   ! s-channel
   !--------------
   mass(2) = mC(i1)
   mass(3) = -mC(i2)
   mass(4) = mC(i3)
   coup1(1) = cpl_CCS0_L(i1,i_in,n1)
   coup1(2) = cpl_CCS0_R(i1,i_in,n1)
   coup1(3) = cpl_CCS0_L(i3,i2,n1)
   coup1(4) = cpl_CCS0_R(i3,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   !--------------
   ! t-channel
   !--------------
   mass(2) = mC(i3)
   mass(3) = -mC(i2)
   mass(4) = mC(i1)
   coup1(1) = cpl_CCS0_L(i3,i_in,n1)
   coup1(2) = cpl_CCS0_R(i3,i_in,n1)
   coup1(3) = cpl_CCS0_L(i1,i2,n1)
   coup1(4) = cpl_CCS0_R(i1,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                         &, IntegralsS04, n_S04, resRa, check)
   gCCCSum(Isum) = resR+resRa
   Contribution(Isum) = 'S0_'//Bu(n1)
  End Do ! n1 ,  S^0_n1

  !-------------------
  ! Pseudoscalar P_0
  !-------------------
  Do n1=1,n_P0
   Isum = Isum + 1
   Boson2(1) = mP0(n1)
   Boson2(2) = gP0(n1)
    
   !--------------
   ! s-channel
   !--------------
   mass(2) = mC(i1)
   mass(3) = -mC(i2)
   mass(4) = mC(i3)

   coup1(1) = cpl_CCP0_L(i1,i_in,n1)
   coup1(2) = cpl_CCP0_R(i1,i_in,n1)
   coup1(3) = cpl_CCP0_L(i3,i2,n1)
   coup1(4) = cpl_CCP0_R(i3,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resR, check)
   !--------------
   ! u-channel
   !--------------
   mass(2) = mC(i3)
   mass(3) = -mC(i2)
   mass(4) = mC(i1)

   coup1(1) = cpl_CCP0_L(i3,i_in,n1)
   coup1(2) = cpl_CCP0_R(i3,i_in,n1)
   coup1(3) = cpl_CCP0_L(i1,i2,n1)
   coup1(4) = cpl_CCP0_R(i1,i2,n1)
   Call IntegrateScalarSS(Boson2, mass, coup1, deltaM, epsI &
                        &, IntegralsS04, n_S04, resRa, check)
   gCCCSum(Isum) = resR+resRa
   Contribution(Isum) = 'P0_'//Bu(n1)
  End Do ! n1 ,  P^0_n1

  !-------------------------------
  ! Z boson - Z-boson
  !-------------------------------
  Boson4(1) = mZ
  Boson4(2) = gZ
  Boson4(3) = mZ
  Boson4(4) = gZ
  Isum = Isum + 1

  mass(2) = mC(i3)
  mass(3) = -mC(i2)
  mass(4) = mC(i1)
  coup2(1) = Cpl_CCZ_L(i1,i_in)
  coup2(2) = Cpl_CCZ_R(i1,i_in)
  coup2(3) = Conjg(Cpl_CCZ_L(i3,i_in))
  coup2(4) = Conjg(Cpl_CCZ_R(i3,i_in))
  coup2(5) = Cpl_CCZ_R(i3,i2)
  coup2(6) = Cpl_CCZ_L(i3,i2)
  coup2(7) = Conjg(Cpl_CCZ_R(i1,i2))
  coup2(8) = Conjg(Cpl_CCZ_L(i1,i2))
  Call IntegrategaugeST(Boson4, Mass, coup2, deltaM, epsI   &
                      &, IntegralsZ8, n_Z8, resC, check)
  gCCCSum(Isum) = -2._dp * Real(resC,dp)
  Contribution(Isum) = 'Z Z s-t'

  !-------------------------------
  ! Z boson - S^0 boson  s-channel
  !-------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_S0 
    Isum = Isum + 1
    Boson4(3) = mS0(n1)
    Boson4(4) = gS0(n1)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = Cpl_CCZ_L(i1,i_in)
    coup2(2) = Cpl_CCZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n1))
    coup2(5) = Cpl_CCZ_L(i3,i2)
    coup2(6) = Cpl_CCZ_R(i3,i2)
    coup2(7) = Conjg(cpl_CCS0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_CCS0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = Cpl_CCZ_L(i3,i_in)
    coup2(2) = Cpl_CCZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_CCS0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_CCS0_L(i3,i_in,n1))
    coup2(5) = Cpl_CCZ_L(i1,i2)
    coup2(6) = Cpl_CCZ_R(i1,i2)
    coup2(7) = Conjg(cpl_CCS0_R(i1,i2,n1))
    coup2(8) = Conjg(cpl_CCS0_L(i1,i2,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resCa, check)
    gCCCSum(Isum) = 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'Z S0_'//Bu(n1)//' s-s'
   End Do  ! n1 n_S0

   !----------------------------------
   ! Z boson - S^0 boson  s-t channel
   !----------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_S0 
    Isum = Isum + 1
    Boson4(3) = mS0(n1)
    Boson4(4) = gS0(n1)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = Cpl_CCZ_L(i1,i_in)
    coup2(2) = Cpl_CCZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_CCS0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_CCS0_L(i3,i_in,n1))
    coup2(5) = Cpl_CCZ_L(i3,i2)
    coup2(6) = Cpl_CCZ_R(i3,i2)
    coup2(7) = Conjg(cpl_CCS0_R(i1,i2,n1))
    coup2(8) = Conjg(cpl_CCS0_L(i1,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resC, check)
    !-------------
    ! t-s channel
    !-------------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = Cpl_CCZ_L(i3,i_in)
    coup2(2) = Cpl_CCZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n1))
    coup2(5) = Cpl_CCZ_R(i1,i2)
    coup2(6) = Cpl_CCZ_L(i1,i2)
    coup2(7) = Conjg(cpl_CCS0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_CCS0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCa, check)
    gCCCSum(Isum) = 2._dp * Real(resCa-resC)
    Contribution(Isum) = 'Z S0_'//Bu(n1)//' s-t'
   End Do  ! n1 n_S0

  !-------------------------------
  ! Z boson - P^0 boson  s-channel
  !-------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_P0 
    Isum = Isum + 1
    Boson4(3) = mP0(n1)
    Boson4(4) = gP0(n1)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = Cpl_CCZ_L(i1,i_in)
    coup2(2) = Cpl_CCZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n1))
    coup2(5) = Cpl_CCZ_L(i3,i2)
    coup2(6) = Cpl_CCZ_R(i3,i2)
    coup2(7) = Conjg(cpl_CCP0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_CCP0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = Cpl_CCZ_L(i3,i_in)
    coup2(2) = Cpl_CCZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_CCP0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_CCP0_L(i3,i_in,n1))
    coup2(5) = Cpl_CCZ_L(i1,i2)
    coup2(6) = Cpl_CCZ_R(i1,i2)
    coup2(7) = Conjg(cpl_CCP0_R(i1,i2,n1))
    coup2(8) = Conjg(cpl_CCP0_L(i1,i2,n1))
    Call IntegrateGaugeSscalarS(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS04, n_ZS04, resCa, check)
    gCCCSum(Isum) = 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'Z P0_'//Bu(n1)//' s-s'
   End Do  ! n1 n_P0

   !----------------------------------
   ! Z boson - P^0 boson  s-t channel
   !----------------------------------
   Boson4(1) = mZ
   Boson4(2) = gZ

   Do n1=1,n_P0 
    Isum = Isum + 1
    Boson4(3) = mP0(n1)
    Boson4(4) = gP0(n1)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = Cpl_CCZ_L(i1,i_in)
    coup2(2) = Cpl_CCZ_R(i1,i_in)
    coup2(3) = Conjg(cpl_CCP0_R(i3,i_in,n1))
    coup2(4) = Conjg(cpl_CCP0_L(i3,i_in,n1))
    coup2(5) = Cpl_CCZ_L(i3,i2)
    coup2(6) = Cpl_CCZ_R(i3,i2)
    coup2(7) = Conjg(cpl_CCP0_R(i1,i2,n1))
    coup2(8) = Conjg(cpl_CCP0_L(i1,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resC, check)
    !-------------
    ! t-s channel
    !-------------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = Cpl_CCZ_L(i3,i_in)
    coup2(2) = Cpl_CCZ_R(i3,i_in)
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n1))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n1))
    coup2(5) = Cpl_CCZ_R(i1,i2)
    coup2(6) = Cpl_CCZ_L(i1,i2)
    coup2(7) = Conjg(cpl_CCP0_R(i3,i2,n1))
    coup2(8) = Conjg(cpl_CCP0_L(i3,i2,n1))
    Call IntegrateGaugeSscalarT(Boson4, Mass, coup2, deltaM, epsI &
                              &, IntegralsZS08, n_ZS08, resCa, check)
    gCCCSum(Isum) = 2._dp * Real(resCa-resC)
    Contribution(Isum) = 'Z P0_'//Bu(n1)//' s-t'
   End Do  ! n1 n_P0

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  s-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = n1+1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = cpl_CCS0_L(i1,i_in,n1)
    coup2(2) = cpl_CCS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_CCS0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCS0_L(i1,i_in,n2))
    coup2(5) = cpl_CCS0_L(i3,i2,n1)
    coup2(6) = cpl_CCS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_CCS0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_CCS0_L(i3,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = cpl_CCS0_L(i3,i_in,n1)
    coup2(2) = cpl_CCS0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_CCS0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_CCS0_L(i3,i_in,n2))
    coup2(5) = cpl_CCS0_L(i1,i2,n1)
    coup2(6) = cpl_CCS0_R(i1,i2,n1)
    coup2(7) = Conjg(cpl_CCS0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_CCS0_L(i1,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    gCCCSum(Isum) = 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' S0_'//Bu(n2)//' s-s'
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - S^0_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = 1,n_S0
    Isum = Isum + 1
    Boson4(3) = mS0(n2)
    Boson4(4) = gS0(n2)
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = cpl_CCS0_L(i1,i_in,n1)
    coup2(2) = cpl_CCS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_CCS0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_CCS0_L(i3,i_in,n2))
    coup2(5) = cpl_CCS0_L(i3,i2,n1)
    coup2(6) = cpl_CCS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_CCS0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_CCS0_L(i1,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    gCCCSum(Isum) = - 2._dp * Real(resC,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' S0_'//Bu(n2)//' s-t'
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = cpl_CCS0_L(i1,i_in,n1)
    coup2(2) = cpl_CCS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n2))
    coup2(5) = cpl_CCS0_L(i3,i2,n1)
    coup2(6) = cpl_CCS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i3,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = cpl_CCS0_L(i3,i_in,n1)
    coup2(2) = cpl_CCS0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_CCP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i3,i_in,n2))
    coup2(5) = cpl_CCS0_L(i1,i2,n1)
    coup2(6) = cpl_CCS0_R(i1,i2,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i1,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    gCCCSum(Isum) = 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' P0_'//Bu(n2)//' s-s'
   End Do
  End Do

  !-------------------------------------------------
  ! S^0_i boson  s-channel - P^0_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_S0
   Boson4(1) = mS0(n1)
   Boson4(2) = gS0(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = cpl_CCS0_L(i1,i_in,n1)
    coup2(2) = cpl_CCS0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_CCP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i3,i_in,n2))
    coup2(5) = cpl_CCS0_L(i3,i2,n1)
    coup2(6) = cpl_CCS0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i1,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    !-----------
    ! t-s channel
    !-----------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = cpl_CCS0_L(i3,i_in,n1)
    coup2(2) = cpl_CCS0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n2))
    coup2(5) = cpl_CCS0_L(i1,i2,n1)
    coup2(6) = cpl_CCS0_R(i1,i2,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i3,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resCa, check)
    gCCCSum(Isum) = - 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'S0_'//Bu(n1)//' P0_'//Bu(n2)//' s-t'
   End Do
  End Do

  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  s-channel
  !-------------------------------------------------
  Do n1 = 1,n_P0
   Boson4(1) = mP0(n1)
   Boson4(2) = gP0(n1)
   Do n2 = n1+1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-----------
    ! s-channel
    !-----------
    mass(2) = mC(i1)
    mass(3) = -mC(i2)
    mass(4) = mC(i3)
    coup2(1) = cpl_CCP0_L(i1,i_in,n1)
    coup2(2) = cpl_CCP0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_CCP0_R(i1,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i1,i_in,n2))
    coup2(5) = cpl_CCP0_L(i3,i2,n1)
    coup2(6) = cpl_CCP0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i3,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i3,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resC, check)
    !-----------
    ! t-channel
    !-----------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = cpl_CCP0_L(i3,i_in,n1)
    coup2(2) = cpl_CCP0_R(i3,i_in,n1)
    coup2(3) = Conjg(cpl_CCP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i3,i_in,n2))
    coup2(5) = cpl_CCP0_L(i1,i2,n1)
    coup2(6) = cpl_CCP0_R(i1,i2,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i1,i2,n2))
    Call IntegrateScalarS1S2(Boson4, Mass, coup2, deltaM, epsI &
                           &, IntegralsS0P04, n_S0P04, resCa, check)
    gCCCSum(Isum) = 2._dp * Real(resC+resCa,dp)
    Contribution(Isum) = 'P0_'//Bu(n1)//' P0_'//Bu(n2)//' s-s'
   End Do
  End Do

  !-------------------------------------------------
  ! P^0_i boson  s-channel - P^0_j boson  t-channel
  !-------------------------------------------------
  Do n1 = 1,n_P0
   Boson4(1) = mP0(n1)
   Boson4(2) = gP0(n1)
   Do n2 = 1,n_P0
    Isum = Isum + 1
    Boson4(3) = mP0(n2)
    Boson4(4) = gP0(n2)
    !-------------
    ! s-t channel
    !-------------
    mass(2) = mC(i3)
    mass(3) = -mC(i2)
    mass(4) = mC(i1)
    coup2(1) = cpl_CCP0_L(i1,i_in,n1)
    coup2(2) = cpl_CCP0_R(i1,i_in,n1)
    coup2(3) = Conjg(cpl_CCP0_R(i3,i_in,n2))
    coup2(4) = Conjg(cpl_CCP0_L(i3,i_in,n2))
    coup2(5) = cpl_CCP0_L(i3,i2,n1)
    coup2(6) = cpl_CCP0_R(i3,i2,n1)
    coup2(7) = Conjg(cpl_CCP0_R(i1,i2,n2))
    coup2(8) = Conjg(cpl_CCP0_L(i1,i2,n2))
    Call IntegrateScalarST(Boson4, Mass, coup2, deltaM, epsI &
                         &, IntegralsS0P08, n_S0P08, resC, check)
    gCCCSum(Isum) = - 2._dp * Real(resC,dp)
    Contribution(Isum) = 'P0_'//Bu(n1)//' P0_'//Bu(n2)//' s-t'
   End Do
  End Do

  !----------
  ! Summing
  !----------
  gCCC = Sum( gCCCSum(1:Isum) )
  If (gCCC.Lt.0._dp) Then
   Write (ErrCan,*) 'Error in Subroutine '//NameOfUnit(Iname)
   Write (ErrCan,*) &
     & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//Bu(i2)//Bu(i3)//') < 0 :' &
     & ,i_in,i1,i2,i3, gCCC
   Write (ErrCan,*) 'The different contributions are :'
   Do n1=1,Isum
    If (gCCCSum(n1).Ne.0._dp)  Write (ErrCan,*) Contribution(n1),gCCCSum(n1)
   End Do
   gCCC = 0._dp
  End If

  !----------------------------------------------
  ! symmetry 
  !----------------------------------------------
  If (i1.Eq.i3) gCCC = 0.5_dp * gCCC 

  gCCC = fac * gCCC

  !---------------------------
  ! for detailed information
  !---------------------------
  If (Present(WriteContribution).And.Present(n_out)) Then

   gCCCSum = gCCCSum * fac
 
   Write (ErrCan,*) &
     & 'Gamma(Chi_'//Bu(i_in)//' -> Chi_'//Bu(i1)//Bu(i2)//Bu(i3)//') :' &
     & ,i_in,i1,i2,i3
   Write (ErrCan,*) 'The different contributions are :'
   Do n1=1,Isum
    If (gCCCSum(n1).Ne.0._dp)  Write (ErrCan,*) Contribution(n1),gCCCSum(n1)
   End Do
  End If

  Deallocate( gCCCSum, Contribution )
   
  Iname = Iname - 1

 End Subroutine ChimToChimChipChim

End Module Chargino3Decays

