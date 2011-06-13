Module BranchingRatios 

! load modules
Use Chargino3Decays
Use Control
Use Couplings
Use Gluino3Decays
Use LoopCouplings
Use Model_Data
Use Neut3Decays
Use StandardModel
Use Slepton3BodyDecays
Use Stop3BodyDecays
Use SusyDecays
! load modules

! interfaces
Interface CalculateBR
 Module Procedure CalculateBR_MSSM, CalculateBR_NMSSM, CalculateBR_RPeps
End Interface
! interfaces

! private variables
!Real(dp), Private :: erg, factor(3), invfac(2)
!Integer, private :: i2
Contains


 Subroutine CalculateBR_MSSM(gauge, mGlu, PhaseGlu, mC, U, V, mN, N, mSneut  &
     & , RSneut, mSlepton, RSlepton, mSup, RSup, mSdown, RSdown, uL_L, uL_R  &
     & , uD_L, uD_R, uU_L, uU_R, mS0, RS0                                    &
     & , mP0, RP0, mSpm, RSpm, epsI, deltaM, CTBD, kont, fac3                &
     & , Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, vevSM, Fgmsb, m32, grav_fac       &
     & , gP_Sl, gT_Sl, BR_Sl, gP_Sn, gT_Sn, BR_Sn, gP_Sd, gT_Sd, BR_Sd       &
     & , gP_Su, gT_Su, BR_Su, gP_C, gT_C, BR_C                               &
     & , gT_N, gP_N2, BR_N2, gP_N3, BR_N3, gP_Glu, gT_Glu, BR_Glu            &
     & , gP_P0, gT_P0, BR_P0       &
     & , gP_S0, gT_S0, BR_S0        &
     & , gP_Spm, gT_Spm, BR_Spm) 
 !------------------------------------------------------------------
 ! Calculates the branching of SUSY particles within the MSSM
 ! it is assumed that the SUSY couplings as well as the parameters
 ! are stored in the public variables in the module SusyParameters.
 ! Input: - gauge(i) .... the gauge couplings
 !        - epsI ........ precision to which the integrals for the
 !                        3-body decays are evolved
 !        - deltaM ...... maximal ratio of mass over phasespace in 
 !                        3-body decays where the masses are set 0
 !                        in the calculation of the phase space.
 !        - kont ........ control variable, is 0 if everything is o.k.
 !        - CTBD ........ logical variable, it .true. then all 3-body
 !          decays are calculated, even if 2-body decays are avaiable.
 !          Otherwisethe 3-body decays are only calculated if 2-body decays
 !          are kinematically forbidden
 !        - fac3 ........ if a total two-body decay width devided by the
 !          the correcponding mass is smaller than fac3, the width
 !          will be recalculated using the routines for 3-body decays
 !  the variable GenerationMixing is taken from the Module InputOutput
 ! The exact form of the output depends partly on the variable
 ! GenerationMixing. 
 ! output:
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  
 ! written by Werner Porod, 25.04.02
 ! 16.09.02: instead of using globally defined couplings, now locally
 !           defined couplings are used.
 ! 29.11.02: adding three-body decays of the lighter stop
 ! 09.12.02: adding decays:  neutralino_1 -> gravitino + photon
 !                           sleptons -> lepton + photon
 ! 14.08.03: adding decays:  neutralino_i -> neutralino_j + photon
 !                           gluino -> neutralino_i + gluon
 ! 21.08.03: adding new possibility concerning 2-body decay modes versus
 !           3-body decay modes: if the new variable calc_2and3_body is
 !           set true, then both possibilities will be calculated.
 !           However, in the 3-body modes on-shell particles will be
 !           negelected. This requires also a change in the output
 !           format. Moreover, also the charge conjugated final states
 !           will be printed in future. 
 !           Starting with the gluino
 ! 10.09.03: changing ordering for neutralino 2-body decays such that
 !           charge conjugated final states are included
 ! 11.09.03: merging branching ratio array for two- and three-body
 !           decay modes of neutralinos
 ! 19.01.09: m32 is given in eV, defining new variable m_grav to get
 !           value in GeV
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: epsI, deltaM, gauge(3), fac3
  Real(dp), Intent(in) :: mGlu, mC(2), mN(4), mSneut(3), mSlepton(6)         &
         & , mSdown(6), mSup(6), mP0(2), RP0(2,2), mS0(2), RS0(2,2), mSpm(2)
  Complex(dp), Intent(in) :: PhaseGlu, RSpm(2,2), U(2,2), V(2,2), N(4,4) &
         & , RSlepton(6,6), Rsneut(3,3), RSup(6,6), RSdown(6,6), uL_L(3,3) &
         & , uL_R(3,3), uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)
  Complex(dp), Intent(in) :: A_d(3,3), A_l(3,3), A_u(3,3), Y_d(3,3), Y_l(3,3) &
         & , Y_u(3,3), mu
  Real(dp), Intent(in) :: vevSM(2), Fgmsb, m32, grav_fac
  Logical, Intent(in) :: CTBD
  Real(dp), Intent(inout) ::  gP_Sl(6,45), gT_Sl(6), BR_Sl(6,45)             &
         & , gP_Sn(3,30), gT_Sn(3), BR_Sn(3,30)                              &
         & , gP_Sd(6,54), gT_Sd(6), BR_Sd(6,54)                              &
         & , gP_Su(6,66), gT_Su(6), BR_Su(6,66)                              &
         & , gP_C(2,267), gT_C(2), BR_C(2,267)                               &
         & , gT_N(4), gP_N2(4,200), BR_N2(4,200), gP_N3(4,400), BR_N3(4,400) &
         & , gP_Glu(230), gT_Glu, BR_Glu(230)                                &
         & , gP_P0(2,200), gT_P0(2), BR_P0(2,200)                            &
         & , gP_S0(2,200), gT_S0(2), BR_S0(2,200)                            &
         & , gP_Spm(2,200), gT_Spm(2), BR_Spm(2,200)

  Integer :: i1, i2, k_neut
  Real(dp) :: m_grav, gNNnunu(4,3,3,3)   &
    & , gNNll(4,3,3,3), gNNdd(4,3,3,3), gNNuu(4,3,3,3), gNCln(4,2,3,3)       &
    & , gNCDU(4,2,3,3), gNgdd(4,3,3), gNguu(4,3,3), gT_N3(4), gGNuu(4,3,3)   &
    & , gGNdd(4,3,3), gGgluon(4), gGCdu(2,3,3), gCuu(2,1,3,3), gCdd(2,1,3,3) &
    & , gCll(2,1,3,3), gCnunu(2,1,3,3), gCNln(2,4,3,3), gCNDU(2,4,3,3)       &
    & , gCgdu(2,3,3), gStWB(6,3), sinW, gCCMajaron(2,1), gCCCC(2,1,1,1)      &
    & , gCCNN(2,1,4,4), gNNCC(4,3,2,2), gNNNN(4,3,3,3), gNNPhoton(4,3)

  Complex(dp) :: cpl_SmpSlSn(2,6,3), cpl_SmpSdSu(2,6,6)   &
      & , cpl_SmpSnSl(2,3,6), cpl_SmpSuSd(2,6,6), cpl_SmpP03(2,2,2)    &
      & , cpl_SmpP0W(2,2), cpl_SmpS03(2,2,2), cpl_SmpS0W(2,2)          &
      & , cpl_SmpLNu_L(2,3,3), cpl_SmpLNu_R(2,3,3), cpl_SmpDU_L(2,3,3) &
      & , cpl_SmpDU_R(8,3,3), cpl_SmpZ(2,2), cpl_DUW(3,3)
  Real(dp) :: cpl_LLZ_L, cpl_LLZ_R, cpl_DDZ_L, cpl_DDZ_R  &
      & , cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R
  Complex(dp) :: cpl_CCZ_L(2,2), cpl_CCZ_R(2,2)           &
      & , cpl_NNZ_L(4,4), cpl_NNZ_R(4,4), cpl_NNS0_L(4,4,2)            &
      & , cpl_NNS0_R(4,4,2), cpl_NNP0_L(4,4,2), cpl_NNP0_R(4,4,2) 
  Complex(dp) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,4,6), cpl_UNSu_R(3,4,6)         &
      & , cpl_LNSl_L(3,4,6), cpl_LNSl_R(3,4,6), cpl_NuNSn_L(3,4,3)      & 
      & , cpl_NuNSn_R(3,4,3), cpl_DDP0_L(3,3,2), cpl_LLP0_L(3,3,2)      &
      & , cpl_UUP0_L(3,3,2), cpl_DDP0_R(3,3,2), cpl_LLP0_R(3,3,2)       &
      & , cpl_UUP0_R(3,3,2), cpl_DDS0_L(3,3,2), cpl_LLS0_L(3,3,2)       &
      & , cpl_UUS0_L(3,3,2), cpl_DDS0_R(3,3,2), cpl_LLS0_R(3,3,2)       &
      & , cpl_UUS0_R(3,3,2)
  Complex(dp) :: cpl_CUSd_L(2,3,6), cpl_CUSd_R(2,3,6)      &
      & , cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), cpl_CLSn_L(2,3,3)       &
      & , cpl_CLSn_R(2,3,3), cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6)
  Complex(dp) :: cpl_GlGlS0(2), cpl_GGS0(2)
  Complex(dp) :: cpl_P0SdSd(2,6,6), cpl_P0SuSu(2,6,6), cpl_P0SlSl(2,6,6) &
      & , cpl_P0SnSn(2,3,3), cpl_P0S0Z(2,2) 
  Real(dp) :: cpl_P0S03(2,2,2), cpl_S0WWvirt(2), cpl_S0ZZvirt(2), vev
  Complex(dp) :: cpl_S0SdSd(2,6,6), cpl_S0SuSu(2,6,6) &
      & , cpl_S0SlSl(2,6,6), cpl_S0SnSn(2,3,3), cpl_LNuW(3,3)
  Real(dp) :: cpl_S03(2,2,2), cpl_S0WW(2), cpl_S0ZZ(2), cpl_FFpW
  Complex(dp) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6), cpl_SlSnW(6,3) &
      & , cpl_SnSlW(3,6), cpl_SdSdZ(6,6), cpl_SlSlZ(6,6), cpl_SnSnZ(3,3)     &
      & , cpl_SuSuZ(6,6)
  Complex(dp) :: cpl_CCP0_L(2,2,2), cpl_CCP0_R(2,2,2)    &
      & , cpl_CCS0_L(2,2,2), cpl_CCS0_R(2,2,2), cpl_CNW_L(2,4)        &
      & , cpl_CNW_R(2,4), cpl_SmpCN_L(2,2,4), cpl_SmpCN_R(2,2,4), cpl_NGP &
      & , cpl_NGZ, cpl_NGH

  Real(dp) :: mSup2(6), mSdown2(6), mSneut2(3), mSlepton2(6), mC2(2), mN2(4) &
     & , mP02(2), gStCNeu(2), gStWBNeu, gStBSnL(3,3,3), gStBSlN(3,6,3), tanb &
     & , sinW2, cosW
  Complex(dp) :: coup, g_u(3), g_d(3), g_l(3), g_c(2), g_sl(6), g_sd(6), g_su(6)
  Real(dp) :: g_W, g_Hp
  Real(dp), Parameter :: mf_nu(3)=0._dp
  Logical :: OnlySM

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CalculateBR_MSSM'

  !----------------------------------------------------------
  ! first all couplings are calculated, see module Couplings
  !----------------------------------------------------------
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)

  m_grav = 1.e-9_dp * m32

  Call AllCouplings(gauge, Y_l, uL_L, uL_R, Y_d, uD_L, uD_R, Y_u, uU_L, uU_R   &
    & , vevSM, RSpm, RP0, RS0, U, V, N, mu, PhaseGlu, RSlepton, A_l, Rsneut    &
    & , RSup, A_u, RSdown, A_d                                                 &
    & , cpl_SmpSlSn, cpl_SmpSdSu, cpl_SmpSnSl, cpl_SmpSuSd, cpl_SmpP03         &
    & , cpl_SmpP0W, cpl_SmpS03, cpl_SmpS0W, cpl_SmpLNu_L, cpl_SmpLNu_R         &
    & , cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_LLZ_L, cpl_LLZ_R      &
    & , cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R   &
    & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R     &
    & , cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L, cpl_GDSd_R, cpl_DNSd_L             &
    & , cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R, cpl_UNSu_L, cpl_UNSu_R             &
    & , cpl_LNSl_L, cpl_LNSl_R, cpl_NuNSn_L, cpl_NuNSn_R, cpl_DDP0_L           &
    & , cpl_LLP0_L, cpl_UUP0_L, cpl_DDP0_R, cpl_LLP0_R, cpl_UUP0_R             &
    & , cpl_DDS0_L, cpl_LLS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_LLS0_R             &
    & , cpl_UUS0_R, cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R             &
    & , cpl_CLSn_L, cpl_CLSn_R, cpl_CNuSl_L, cpl_CNuSl_R, cpl_GlGlS0           &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0SlSl, cpl_P0SnSn, cpl_P0S0Z, cpl_P0S03   &
    & , cpl_S0SdSd, cpl_S0SuSu, cpl_S0SlSl, cpl_S0SnSn, cpl_S03, cpl_S0WW      &
    & , cpl_S0ZZ, cpl_FFpW, cpl_LNuW, cpl_SdSuW, cpl_SuSdW, cpl_SlSnW          &
    & , cpl_SnSlW, cpl_SdSdZ, cpl_SlSlZ, cpl_SnSnZ, cpl_SuSuZ, cpl_CCP0_L      &
    & , cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L, cpl_CNW_R               &
    & , cpl_SmpCN_L, cpl_SmpCN_R, GenerationMixing)


   
  gP_Sn = 0._dp
  gT_Sn = 0._dp
  BR_Sn = 0._dp
  Call SfermionTwoBodyDecays(-1, mSneut, mf_nu, mf_l                       &
          &, mN, cpl_NuNSn_L, cpl_NuNSn_R, mC, cpl_CLSn_L, cpl_CLSn_R       &
          &, mSlepton, mW, cpl_SnSlW, mZ, cpl_SnSnZ, mSpm, cpl_SmpSnSl      &
          &, mP0, cpl_P0SnSn, mS0, cpl_S0SnSn                               &
          &, 1, GenerationMixing                                            &
          &, gP_Sn, gT_Sn, BR_Sn)

  gP_Sl = 0._dp
  gT_Sl = 0._dp
  BR_Sl = 0._dp
  Call SfermionTwoBodyDecays(-1, mSlepton, mf_l, mf_nu                      &
          &, mN, cpl_LNSl_L, cpl_LNSl_R, mC, cpl_CNuSl_L, cpl_CNuSl_R       &
          &, mSneut, mW, cpl_SlSnW, mZ, cpl_SlSlZ, mSpm, cpl_SmpSlSn        &
          &, mP0, cpl_P0SlSl, mS0, cpl_S0SlSl                               &
          &, 2, GenerationMixing                                            &
          &, gP_Sl, gT_Sl, BR_Sl)

  !-------------------------------------
  ! calculation of decay into gravitino
  ! relevant for GMSB
  !-------------------------------------
   If ((GenerationMixing).And.(mslepton(6).Gt.m_grav)) Then
    Write(ErrCan,*) "Problem in routine "//NameOfUnit(Iname)
    Write(ErrCan,*) "Decay of a slepton into gravitino and lepton is not"
    Write(ErrCan,*) "included in case of generation mixing"
   Else ! .not.GenerationMixing
    
    Do i1=1,6
     If (mSlepton(i1).Gt.m_grav) Then
      gP_Sl(i1,13) = oo16Pi * mSlepton(i1) * (mSlepton(i1)**2-m_grav**2)**2 &
                 &   / (grav_fac*Fgmsb)**2
     Else
      gP_Sl(i1,13) = 0._dp
     End If
     If (gP_Sl(i1,13).Gt.0._dp) Then
      BR_Sl(i1,1:12) =  BR_Sl(i1,1:12) * gT_Sl(i1) / (gT_Sl(i1) + gP_Sl(i1,13))
      gT_Sl(i1) = gT_Sl(i1) + gP_Sl(i1,13)
      BR_Sl(i1,13) = gP_Sl(i1,13) / gT_Sl(i1)
     End If
    End Do
   End If ! GenerationMixing

  gP_Su = 0._dp
  gT_Su = 0._dp
  BR_Su = 0._dp
  Call SfermionTwoBodyDecays(-1, mSup, mf_u, mf_d                           &
          &, mN, cpl_UNSu_L, cpl_UNSu_R, mC, cpl_CDSu_L, cpl_CDSu_R         &
          &, mSdown, mW, cpl_SuSdW, mZ, cpl_SuSuZ, mSpm, cpl_SmpSuSd        &
          &, mP0, cpl_P0SuSu, ms0, cpl_S0SuSu                               &
          &, 0, GenerationMixing                                            &
          &, gP_Su, gT_Su, BR_Su, mGlu, cpl_GUSu_L, cpl_GUSu_R)

  gP_Sd = 0._dp
  gT_Sd = 0._dp
  BR_Sd = 0._dp
  Call SfermionTwoBodyDecays(-1, mSdown, mf_d, mf_u                         &
          &, mN, cpl_DNSd_L, cpl_DNSd_R, mC, cpl_CUSd_L, cpl_CUSd_R         &
          &, mSup, mW, cpl_SdSuW, mZ, cpl_SdSdZ, mSpm, cpl_SmpSdSu          &
          &, mP0, cpl_P0SdSd, mS0, cpl_S0SdSd                               &
          &, 0, GenerationMixing                                            &
          &, gP_Sd, gT_Sd, BR_Sd, mGlu, cpl_GDSd_L, cpl_GDSd_R)


  gP_Glu = 0._dp
  gT_Glu = 0._dp
  BR_Glu = 0._dp

  If (.Not.CTBD) Then
   Call GluinoTwoBodyDecays(mGlu, mSdown, cpl_GDSd_L, cpl_GDSd_R          &  
        & , mf_d, mSup, cpl_GUSu_L, cpl_GUSu_R, mf_u, 0, GenerationMixing &  
        & , gP_Glu(1:72), gT_Glu, BR_Glu(1:72) )

   If (.Not.GenerationMixing) Then
    gP_Glu(25) = GluinoToStopC( mglu, mSup**2, Rsup, mSdown(5:6)**2, mP0**2 &
            &    , vevSM(2)/vevSM(1), y_d(3,3), mu, A_d(3,3), cpl_GUSu_R, kont)
    gP_Glu(26) = gP_Glu(25)
   End If

   If (gT_Glu.Lt.fac3*mglu) Then
    Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .False. , gGNdd, gGNuu           &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )
    If (GenerationMixing) Then
     gP_Glu(1:72) = 0._dp
    Else
     gP_Glu(1:24) = 0._dp
     gP_Glu(27:72) = 0._dp
    End If

   Else ! calculate only 3-body modes via virtual particles

    Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .True. , gGNdd, gGNuu            &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )

   End If

  Else ! calculation of 3-body decay modes is enforced

   If (.Not.GenerationMixing) Then
    gP_Glu(25) = GluinoToStopC( mglu, mSup**2, Rsup, mSdown(5:6)**2, mP0**2 &
           &    , vevSM(2)/vevSM(1), y_d(3,3), mu, A_d(3,3), cpl_GUSu_R, kont)
    gP_Glu(26) = gP_Glu(25)
   End If
   Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .False. , gGNdd, gGNuu           &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )

  End If
  gT_Glu = Sum(gP_Glu) 
  If (gT_Glu.Gt.0._dp) BR_Glu = gP_Glu / gT_Glu 
  BR_Glu2 = 0._dp
  gP_Glu2 = 0._dp
  gP_Glu2(1:76) = gP_Glu(1:76) 
  BR_Glu2(1:76) = BR_Glu(1:76)
  BR_Glu3 = 0._dp
  gP_Glu3 = 0._dp
  gP_Glu3(1:150) = gP_Glu(77:226) 
  BR_Glu3(1:150) = BR_Glu(77:226) 

  gP_S0 = 0._dp
  gT_S0 = 0._dp
  BR_S0 = 0._dp
  vev = Sqrt(vevSM(1)**2+vevSM(2)**2)
  !------------------------------------------------------------------------
  ! couplings to a real and a virtual gauge boson
  !------------------------------------------------------------------------
  cpl_S0WWvirt = cpl_S0WW / vev
  cpl_S0ZZvirt = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * cpl_S0ZZ / vev
  !------------------------------------------------------------------------
  ! loop induced couplings to photons and gluons
  !------------------------------------------------------------------------
  cpl_GGS0 = 0._dp
  mC2 = mC**2
  mSdown2 = mSdown**2
  mSup2 = mSup**2
  mSlepton2 = mSlepton**2
  Do i1=1,2
   g_W = (vevSM(1) * RS0(i1,1) + vevSM(2) * RS0(i1,2) ) / vev
   g_u = RS0(i1,2) * vev / vevSM(2)
   g_d = RS0(i1,1) * vev / vevSM(1)
   g_l = RS0(i1,1) * vev / vevSM(1)
   Forall(i2=1:2) g_c(i2) = cpl_CCS0_L(i2,i2,i1) * mW /( gauge(2) * mC(i2))
   g_Hp = cpl_SmpS03(2,2,i1) * mW /( gauge(2) * mSpm2(2))
   Forall(i2=1:6) g_sl(i2) = - 0.5_dp * cpl_S0SlSl(i1,i2,i2) * vevSM(1) &
                           &          / mSlepton2(i2)
   Forall(i2=1:6) g_sd(i2) = - 0.5_dp * cpl_S0SdSd(i1,i2,i2) * vevSM(1) &
                           &          / mSdown2(i2)
   Forall(i2=1:6) g_su(i2) = - 0.5_dp * cpl_S0SuSu(i1,i2,i2) * vevSM(2) &
                           &          / mSup2(i2)
   Call CoupScalarPhoton(mS02(i1), mW2, g_W, mf_u2, g_u, mf_d2, g_d, mf_l2   &
    & , g_l, mC2, g_c, mSpm2(2), g_Hp, mSup2, g_su, mSdown2, g_sd, mSlepton2 &
    & , g_sl, cpl_GGS0(i1), coup )
   cpl_GGS0(i1) = cpl_GGS0(i1) * oo4pi * gauge(2)**2 * sinW2
   Call CoupScalarGluon(mS02(i1), mf_u2, g_u, mf_d2, g_d, mSup2, g_su &
                      & , mSdown2, g_sd, cpl_GlGlS0(i1), coup )
   cpl_GlGlS0(i1) = cpl_GlGlS0(i1) * oo4pi * gauge(3)**2
  end do

  Call ScalarTwoBodyDecays(-1, mS0, cpl_S03, cpl_GlGlS0, cpl_GGS0            &
          &, mf_l, cpl_LLS0_L, cpl_LLS0_R, mf_d, cpl_DDS0_L, cpl_DDS0_R      & 
          &, mf_u, cpl_UUS0_L, cpl_UUS0_R, mSlepton, cpl_S0SlSl              &
          &, mSneut, cpl_S0SnSn, mSdown, cpl_S0SdSd, mSup, cpl_S0SuSu        & 
          &, mN, cpl_NNS0_L, cpl_NNS0_R, mC, cpl_CCS0_L, cpl_CCS0_R          &
          &, mW, cpl_S0WW, cpl_S0WWvirt, mZ, cpl_S0ZZ, cpl_S0ZZvirt          &
          &, mSpm, cpl_SmpS03, mP0, cpl_P0S03, cpl_P0S0Z, cpl_SmpS0W, mglu   &
          &, GenerationMixing, gP_S0, gT_S0, BR_S0)

  gP_P0 = 0._dp
  gT_P0 = 0._dp
  BR_P0 = 0._dp
  Call PseudoscalarTwoBodyDecays(2, mP0                                      &
          &, mf_l, cpl_LLP0_L, cpl_LLP0_R, mf_d, cpl_DDP0_L, cpl_DDP0_R      & 
          &, mf_u, cpl_UUP0_L, cpl_UUP0_R, mSlepton, cpl_P0SlSl              &
          &, mSdown, cpl_P0SdSd, mSup, cpl_P0SuSu                            & 
          &, mN, cpl_NNP0_L, cpl_NNP0_R, mC, cpl_CCP0_L, cpl_CCP0_R          &
          &, mSpm, cpl_SmpP03, mS0, cpl_P0S03, mZ, cpl_P0S0Z, mW, cpl_SmpP0W &
          &, mglu, GenerationMixing, gP_P0, gT_P0, BR_P0)
  gT_P0(1) = gamZ ! needed for 3-body decays

  gP_Spm = 0._dp
  gT_Spm = 0._dp
  BR_Spm = 0._dp
  Call  ChargedscalarTwoBodyDecays(2, mSpm                                   &
          &, mf_l, cpl_SmpLNu_L, cpl_SmpLNu_R                                &
          &, mf_d, mf_u, cpl_SmpDU_L, cpl_SmpDU_R                            &
          &, mSlepton, mSneut, cpl_SmpSlSn, mSdown, mSup, cpl_SmpSdSu        & 
          &, mN, mC, cpl_SmpCN_L, cpl_SmpCN_R, mW, mZ, cpl_SmpZ              &
          &, mP0, cpl_SmpP03, cpl_SmpP0W, mS0, cpl_SmpS03, cpl_SmpS0W        &
          &, 1, GenerationMixing, gP_Spm, gT_Spm, BR_Spm)
  gT_Spm(1) = gamW ! needed for 3-body decays


  gP_C = 0._dp
  gT_C = 0._dp
  BR_C = 0._dp
  Do i1=1,2
   If (.Not.CTBD)  Then
    Call CharginoTwoBodyDecays(i1, mC, mSlepton, cpl_CNuSl_L, cpl_CNuSl_R    &
        & , mSneut, cpl_CLSn_L, cpl_CLSn_R, mf_l, mSdown, cpl_CUSd_L         &
        & , cpl_CUSd_R, mf_u, mSup, cpl_CDSu_L, cpl_CDSu_R, mf_d             &
        & , mN, mW, cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R     &
        & , mZ, cpl_CCZ_L, cpl_CCZ_R, mP0, cpl_CCP0_L, cpl_CCP0_R            &
        & , mS0, cpl_CCS0_L, cpl_CCS0_R                                      &
        & , 1, GenerationMixing, gP_C(:,1:63), gT_C, BR_C(:,1:63) )

    If (gT_C(i1).Lt. fac3*Abs(mC(i1))) Then
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .False., gCCMajaron        &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gT_C  &
      & , gP_C(:,64:267), BR_C(:,64:267))
     gP_C(i1,1:63) = 0._dp

    Else ! calculate 3-body decays with virtual interemdiate states only
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .True., gCCMajaron         &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gT_C  &
      & , gP_C(:,64:267), BR_C(:,64:267))
    End If

   Else ! enforce calculation of 3-body final states 
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .False., gCCMajaron        &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gT_C  &
      & , gP_C(:,64:267), BR_C(:,64:267))

   End If

   gT_C(i1) = Sum(gP_C(i1,:))
   If (gT_C(i1).Gt.0._dp) BR_C(i1,:) = gP_C(i1,:) / gT_C(i1)
  End Do
  BR_C2 = 0._dp
  BR_C2(:,1:63) = BR_C(:,1:63)
  gP_C2 = 0._dp
  gP_C2(:,1:63) = gP_C(:,1:63)
  BR_C3 = 0._dp
  BR_C3(:,1:204) = BR_C(:,64:267)
  gP_C3 = 0._dp
  gP_C3(:,1:204) = gP_C(:,64:267)
  
  gP_N2 = 0._dp
  BR_N2 = 0._dp
  gP_N3 = 0._dp
  BR_N3 = 0._dp
  OnlySM = .True.

  Do i1=1,4
   If (.Not.CTBD) Then
    sinW = Sqrt(sinW2)
    cosW = Sqrt(1._dp - sinW2)    
    cpl_NGP = (N(i1,1)*cosW+N(i1,2)*sinW) / (grav_fac*Fgmsb)
    cpl_NGZ = (- N(i1,1)*sinW+ N(i1,2)*cosW) / (grav_fac*Fgmsb)
    cpl_NGH = (N(i1,3) * RS0(1,1) - N(i1,4)*RS0(1,2)) / (grav_fac*Fgmsb)
    Call NeutralinoTwoBodyDecays(i1, mN, mSlepton, cpl_LNSl_L, cpl_LNSl_R   &
       & , mf_l, mSneut, cpl_NuNSn_L, cpl_NuNSn_R, mSdown, cpl_DNSd_L       &
       & , cpl_DNSd_R, mf_d, mSup, cpl_UNSu_L, cpl_UNSu_R, mf_u, mC, mW     &
       & , cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R, mZ         &
       & , cpl_NNZ_L, cpl_NNZ_R, mP0, cpl_NNP0_L, cpl_NNP0_R, mS0           &
       & , cpl_NNS0_L, cpl_NNS0_R, m_grav, cpl_NGP, cpl_NGZ, cpl_NGH, 1        &
       & , GenerationMixing, gP_N2, gT_N, BR_N2 )

    If (gT_N(i1).Lt. fac3*Abs(mN(i1))) Then
     Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L            &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .False., 0._dp            &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 2500, gT_N3, gP_N3, BR_N3 )
     gP_N2 = 0._dp

    Else ! calculate only 3-body via virtual particles
     Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L            &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .True., 0.0_dp            &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 2500, gT_N3, gP_N3, BR_N3 )
    End If

   Else ! calculation of 3-body decay modes is enforced
    OnlySM = .True.
    If (Size(mN).Gt.4) OnlySM = .False.
      Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L           &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .False., 0.0_dp           &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 2500, gT_N3, gP_N3, BR_N3 )

    End If ! CTBD

   !------------------------------------------
   ! decay into photons are 2-body decays
   !------------------------------------------
   If (i1.Gt.1) Then
    gP_N2(i1,151:153) =  gP_N3(i1,1:3)
    gP_N3(i1,1:197) = gP_N3(i1,4:200)
   End If
   gT_N(i1) = Sum(gP_N2(i1,:)) + Sum(gP_N3(i1,:))
   If (gT_N(i1).Gt.0._dp) Then
    BR_N2(i1,:) = gP_N2(i1,:) / gT_N(i1)
    BR_N3(i1,:) = gP_N3(i1,:) / gT_N(i1)
   End If
  End Do

  !-------------------------------------------------------
  ! in the case of the lighter stop it is possible that
  ! all two-body decay modes are kinematically forbidden 
  !-------------------------------------------------------
  mSup2 = mSup**2
  mSdown2 = mSdown**2
  mSneut2 = mSneut**2
  mSlepton2 = mSlepton**2
  mC2 = mC**2
  mN2 = mN**2
  mP02 = mP0**2
  tanb = vevSM(2) / vevSM(1)
  If ((gT_Su(5).Lt.fac3*mSup(5)).Or.(CTBD.and.(.not.GenerationMixing))) Then ! calculation including widths
   Call StopDecays3(mSup, mSup2, RSup, Y_d(3,3), tanb, mSdown2, gT_Sd        &
        & , A_d(3,3), mu, mP02, mN, mN2, gauge(2), cpl_UNSu_L, cpl_UNSu_R    &
        & , cpl_DNSd_L, cpl_DNSd_R, cpl_SdSuW, mSneut2, mC, mC2, gT_C        &
        & , cpl_CLSn_L, cpl_CLSn_R, cpl_CDSu_L, cpl_CDSu_R, cpl_CNW_L        &
        & , cpl_CNW_R, mSlepton2, cpl_CNuSl_R, epsI, gstCneu, gStWBNeu       &
        & , gStBSnL, gStBSlN, kont, .False.)
   gT_Su(5) =  gstWbneu + Sum(gstCneu) + Sum(gStBSnL) + Sum(gStBSlN)
   gP_Su(5,:) = 0._dp
   gP_Su(5,55:56) = gstCneu
   gP_Su(5,57) = gstWbneu
   gP_Su(5,58) = gStBSnL(3,1,1)
   gP_Su(5,59) = gStBSnL(3,2,2)
   gP_Su(5,60) = gStBSnL(3,3,3)
   gP_Su(5,61) = gStBSlN(3,1,1)
   gP_Su(5,62) = gStBSlN(3,2,1)
   gP_Su(5,63) = gStBSlN(3,3,2)
   gP_Su(5,64) = gStBSlN(3,4,2)
   gP_Su(5,65) = gStBSlN(3,5,3)
   gP_Su(5,66) = gStBSlN(3,6,3)
   BR_Su(5,:) = gP_Su(5,:) / gT_Su(5)

  Else  ! calculation excluding widths
   Call StopDecays3(mSup, mSup2, RSup, Y_d(3,3), tanb, mSdown2, gT_Sd        &
        & , A_d(3,3), mu, mP02, mN, mN2, gauge(2), cpl_UNSu_L, cpl_UNSu_R    &
        & , cpl_DNSd_L, cpl_DNSd_R, cpl_SdSuW, mSneut2, mC, mC2, gT_C        &
        & , cpl_CLSn_L, cpl_CLSn_R, cpl_CDSu_L, cpl_CDSu_R, cpl_CNW_L        &
        & , cpl_CNW_R, mSlepton2, cpl_CNuSl_R, epsI, gstCneu, gStWBNeu       &
        & , gStBSnL, gStBSlN, kont, .True.)

   gP_Su(5,55:56) = gstCneu
   gP_Su(5,57) = gstWbneu
   gP_Su(5,58) = gStBSnL(3,1,1)
   gP_Su(5,59) = gStBSnL(3,2,2)
   gP_Su(5,60) = gStBSnL(3,3,3)
   gP_Su(5,61) = gStBSlN(3,1,1)
   gP_Su(5,62) = gStBSlN(3,2,1)
   gP_Su(5,63) = gStBSlN(3,3,2)
   gP_Su(5,64) = gStBSlN(3,4,2)
   gP_Su(5,65) = gStBSlN(3,5,3)
   gP_Su(5,66) = gStBSlN(3,6,3)
   gT_Su(5) = Sum(gP_Su(5,:))
   If (gT_Su(5).Gt.0._dp) BR_Su(5,:) = gP_Su(5,:) / gT_Su(5)
   gP_Su3 = 0._dp
   BR_Su3 = 0._dp
   gP_Su3(5,1:10) = gP_Su(5,57:66) 
   BR_Su3(5,1:10) = BR_Su(5,57:66) 

  End If

  Iname = Iname - 1

 End Subroutine CalculateBR_MSSM

 Subroutine CalculateBR_NMSSM(gauge, mGlu, PhaseGlu, mC, U, V, mN, N, mSneut  &
     & , RSneut, mSlepton, RSlepton, mSup, RSup, mSdown, RSdown, uL_L, uL_R   &
     & , uD_L, uD_R, uU_L, uU_R, mS0, RS0, mP0, RP0, mSpm, RSpm               &
     & , epsI, deltaM, CTBD, kont, fac3, Y_d, A_d, Y_l, A_l, Y_u, A_u         &
     & , mu, h0, Ah0, lam, Alam, vevSM, vP, Fgmsb, m32, grav_fac              &
     & , gP_Sl, gT_Sl, BR_Sl, gP_Sn, gT_Sn, BR_Sn, gP_Sd, gT_Sd, BR_Sd        &
     & , gP_Su, gT_Su, BR_Su, gP_C, gT_C, BR_C                                &
     & , gP_N, gT_N, BR_N, gP_Glu, gT_Glu, BR_Glu, gP_P0, gT_P0, BR_P0        &
     & , gP_S0, gT_S0, BR_S0, gP_Spm, gT_Spm, BR_Spm) 
 !------------------------------------------------------------------
 ! Calculates the branching of SUSY particles within the MSSM
 ! it is assumed that the SUSY couplings as well as the parameters
 ! are stored in the public variables in the module SusyParameters.
 ! Input: - gauge(i) .... the gauge couplings
 !        - epsI ........ precision to which the integrals for the
 !                        3-body decays are evolved
 !        - deltaM ...... maximal ratio of mass over phasespace in 
 !                        3-body decays where the masses are set 0
 !                        in the calculation of the phase space.
 !        - kont ........ control variable, is 0 if everything is o.k.
 !        - CTBD ........ logical variable, it .true. then all 3-body
 !          decays are calculated, even if 2-body decays are avaiable.
 !          Otherwisethe 3-body decays are only calculated if 2-body decays
 !          are kinematically forbidden
 !        - fac3 ........ if a total two-body decay width devided by the
 !          the correcponding mass is smaller than fac3, the width
 !          will be recalculated using the routines for 3-body decays
 !  the variable GenerationMixing is taken from the Module InputOutput
 ! The exact form of the output depends partly on the variable
 ! GenerationMixing. 
 ! output:
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  
 ! written by Werner Porod, 25.04.02
 ! 16.09.02: instead of using globally defined couplings, now locally
 !           defined couplings are used.
 ! 29.11.02: adding three-body decays of the lighter stop
 ! 09.12.02: adding decays:  neutralino_1 -> gravitino + photon
 !                           sleptons -> lepton + photon
 ! 14.08.03: adding decays:  neutralino_i -> neutralino_j + photon
 !                           gluino -> neutralino_i + gluon
 ! 21.08.03: adding new possibility concerning 2-body decay modes versus
 !           3-body decay modes: if the new variable calc_2and3_body is
 !           set true, then both possibilities will be calculated.
 !           However, in the 3-body modes on-shell particles will be
 !           negelected. This requires also a change in the output
 !           format. Moreover, also the charge conjugated final states
 !           will be printed in future. 
 !           Starting with the gluino
 ! 10.09.03: changing ordering for neutralino 2-body decays such that
 !           charge conjugated final states are included
 ! 11.09.03: merging branching ratio array for two- and three-body
 !           decay modes of neutralinos
 ! 19.01.09: m32 is given in eV, defining new variable m_grav to get
 !           value in GeV
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: epsI, deltaM, gauge(3), fac3
  Real(dp), Intent(in) :: mGlu, mC(2), mN(5), mSneut(3), mSlepton(6)         &
         & , mSdown(6), mSup(6), mP0(3), RP0(3,3), mS0(3), RS0(3,3), mSpm(2)
  Complex(dp), Intent(in) :: PhaseGlu, RSpm(2,2), U(2,2), V(2,2), N(5,5) &
         & , RSlepton(6,6), Rsneut(3,3), RSup(6,6), RSdown(6,6), uL_L(3,3) &
         & , uL_R(3,3), uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)
  Complex(dp), Intent(in) :: A_d(3,3), A_l(3,3), A_u(3,3), Y_d(3,3), Y_l(3,3) &
         & , Y_u(3,3), mu, h0, Ah0, lam, Alam
  Real(dp), Intent(in) :: vevSM(2), vP, Fgmsb, m32, grav_fac
  Logical, Intent(in) :: CTBD
  Real(dp), Intent(inout) ::  gP_Sl(6,45), gT_Sl(6), BR_Sl(6,45)             &
         & , gP_Sn(3,30), gT_Sn(3), BR_Sn(3,30)                              &
         & , gP_Sd(6,54), gT_Sd(6), BR_Sd(6,54)                              &
         & , gP_Su(6,66), gT_Su(6), BR_Su(6,66)                              &
         & , gP_C(2,267), gT_C(2), BR_C(2,267)                               &
         & , gP_N(5,350), gT_N(5), BR_N(5,350)                               &
         & , gP_Glu(230), gT_Glu, BR_Glu(230)                                &
         & , gP_P0(3,200), gT_P0(3), BR_P0(3,200)                            &
         & , gP_S0(3,200), gT_S0(3), BR_S0(3,200)                            &
         & , gP_Spm(2,200), gT_Spm(2), BR_Spm(2,200)

  Integer :: i1, k_neut
  Real(dp) :: m_grav, gNNnunu(5,4,3,3)   &
    & , gNNll(5,4,3,3), gNNdd(5,4,3,3), gNNuu(5,4,3,3), gNCln(5,2,3,3)       &
    & , gNCDU(5,2,3,3), gNgdd(5,3,3), gNguu(5,3,3), gT_N3(5), gGNuu(5,3,3)   &
    & , gGNdd(5,3,3), gGgluon(5), gGCdu(2,3,3), gCuu(2,1,3,3), gCdd(2,1,3,3) &
    & , gCll(2,1,3,3), gCnunu(2,1,3,3), gCNln(2,5,3,3), gCNDU(2,5,3,3)       &
    & , gCgdu(2,3,3), gStWB(6,3), sinW, gCCMajaron(2,1), gCCCC(2,1,1,1)      &
    & , gCCNN(2,1,5,5), gNNCC(5,4,2,2), gNNNN(5,4,4,4), gNNPhoton(5,4)

  Complex(dp) :: cpl_SmpSlSn(2,6,3), cpl_SmpSdSu(2,6,6)   &
      & , cpl_SmpSnSl(2,3,6), cpl_SmpSuSd(2,6,6), cpl_SmpP03(2,2,3)    &
      & , cpl_SmpP0W(2,3), cpl_SmpS03(2,2,3), cpl_SmpS0W(2,3)          &
      & , cpl_SmpLNu_L(2,3,3), cpl_SmpLNu_R(2,3,3), cpl_SmpDU_L(2,3,3) &
      & , cpl_SmpDU_R(8,3,3), cpl_SmpZ(2,2), cpl_DUW(3,3)
  Real(dp) :: cpl_LLZ_L, cpl_LLZ_R, cpl_DDZ_L, cpl_DDZ_R  &
      & , cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R
  Complex(dp) :: cpl_CCZ_L(2,2), cpl_CCZ_R(2,2)           &
      & , cpl_NNZ_L(5,5), cpl_NNZ_R(5,5), cpl_NNS0_L(5,5,3)            &
      & , cpl_NNS0_R(5,5,3), cpl_NNP0_L(5,5,3), cpl_NNP0_R(5,5,3) 
  Complex(dp) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,5,6), cpl_DNSd_R(3,5,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,5,6), cpl_UNSu_R(3,5,6)         &
      & , cpl_LNSl_L(3,5,6), cpl_LNSl_R(3,5,6), cpl_NuNSn_L(3,5,3)      & 
      & , cpl_NuNSn_R(3,5,3), cpl_DDP0_L(3,3,3), cpl_LLP0_L(3,3,3)      &
      & , cpl_UUP0_L(3,3,3), cpl_DDP0_R(3,3,3), cpl_LLP0_R(3,3,3)       &
      & , cpl_UUP0_R(3,3,3), cpl_DDS0_L(3,3,3), cpl_LLS0_L(3,3,3)       &
      & , cpl_UUS0_L(3,3,3), cpl_DDS0_R(3,3,3), cpl_LLS0_R(3,3,3)       &
      & , cpl_UUS0_R(3,3,3)
  Complex(dp) :: cpl_CUSd_L(2,3,6), cpl_CUSd_R(2,3,6)      &
      & , cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), cpl_CLSn_L(2,3,3)       &
      & , cpl_CLSn_R(2,3,3), cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6)
  Complex(dp) :: cpl_GlGlS0(3), cpl_GGS0(3)
  Complex(dp) :: cpl_P0SdSd(3,6,6), cpl_P0SuSu(3,6,6), cpl_P0SlSl(3,6,6) &
      & , cpl_P0SnSn(3,3,3), cpl_P0S0Z(3,3) 
  Real(dp) :: cpl_P0S03(3,3,3), cpl_S0WWvirt(3), cpl_S0ZZvirt(3), vev
  Complex(dp) :: cpl_S0SdSd(3,6,6), cpl_S0SuSu(3,6,6) &
      & , cpl_S0SlSl(3,6,6), cpl_S0SnSn(3,3,3), cpl_LNuW(3,3)
  Real(dp) :: cpl_S03(3,3,3), cpl_S0WW(3), cpl_S0ZZ(3), cpl_FFpW
  Complex(dp) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6), cpl_SlSnW(6,3) &
      & , cpl_SnSlW(3,6), cpl_SdSdZ(6,6), cpl_SlSlZ(6,6), cpl_SnSnZ(3,3)     &
      & , cpl_SuSuZ(6,6)
  Complex(dp) :: cpl_CCP0_L(2,2,3), cpl_CCP0_R(2,2,3)                     &
      & , cpl_CCS0_L(2,2,3), cpl_CCS0_R(2,2,3), cpl_CNW_L(2,5)            &
      & , cpl_CNW_R(2,5), cpl_SmpCN_L(2,2,5), cpl_SmpCN_R(2,2,5), cpl_NGP &
      & , cpl_NGZ, cpl_NGH

  Real(dp) :: mSup2(6), mSdown2(6), mSneut2(3), mSlepton2(6), mC2(2), mN2(5) &
     & , mP02(3), gStCNeu(2), gStWBNeu, gStBSnL(3,3,3), gStBSlN(3,6,3), tanb &
     & , sinW2, cosW
  Real(dp), Parameter :: mf_nu(3)=0._dp
  Logical :: OnlySM

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CalculateBR_NMSSM'

  !----------------------------------------------------------
  ! first all couplings are calculated, see module Couplings
  !----------------------------------------------------------
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  m_grav = 1.e-9_dp * m32
  
  Call AllCouplings(gauge, Y_l, uL_L, uL_R, Y_d, uD_L, uD_R, Y_u, uU_L, uU_R &
    & , vevSM, h0, Ah0, lam, Alam, vP, RSpm, RP0, RS0, U, V, N, mu, PhaseGlu &
    & , RSlepton, A_l, Rsneut, RSup, A_u, RSdown, A_d                        &
    & , cpl_SmpSlSn, cpl_SmpSdSu, cpl_SmpSnSl, cpl_SmpSuSd, cpl_SmpP03       &
    & , cpl_SmpP0W, cpl_SmpS03, cpl_SmpS0W, cpl_SmpLNu_L, cpl_SmpLNu_R       &
    & , cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_LLZ_L, cpl_LLZ_R    &
    & , cpl_DDZ_L, cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_NuNuZ_L, cpl_NuNuZ_R &
    & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R   &
    & , cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L, cpl_GDSd_R, cpl_DNSd_L           &
    & , cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R, cpl_UNSu_L, cpl_UNSu_R           &
    & , cpl_LNSl_L, cpl_LNSl_R, cpl_NuNSn_L, cpl_NuNSn_R, cpl_DDP0_L         &
    & , cpl_LLP0_L, cpl_UUP0_L, cpl_DDP0_R, cpl_LLP0_R, cpl_UUP0_R           &
    & , cpl_DDS0_L, cpl_LLS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_LLS0_R           &
    & , cpl_UUS0_R, cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R           &
    & , cpl_CLSn_L, cpl_CLSn_R, cpl_CNuSl_L, cpl_CNuSl_R, cpl_GlGlS0         &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0SlSl, cpl_P0SnSn, cpl_P0S0Z, cpl_P0S03 &
    & , cpl_S0SdSd, cpl_S0SuSu, cpl_S0SlSl, cpl_S0SnSn, cpl_S03, cpl_S0WW    &
    & , cpl_S0ZZ, cpl_FFpW, cpl_LNuW, cpl_SdSuW, cpl_SuSdW, cpl_SlSnW        &
    & , cpl_SnSlW, cpl_SdSdZ, cpl_SlSlZ, cpl_SnSnZ, cpl_SuSuZ, cpl_CCP0_L    &
    & , cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L, cpl_CNW_R             &
    & , cpl_SmpCN_L, cpl_SmpCN_R, GenerationMixing)


   
  gP_Sn = 0._dp
  gT_Sn = 0._dp
  BR_Sn = 0._dp
  Call SfermionTwoBodyDecays(-1, mSneut, mf_nu, mf_l                       &
          &, mN, cpl_NuNSn_L, cpl_NuNSn_R, mC, cpl_CLSn_L, cpl_CLSn_R       &
          &, mSlepton, mW, cpl_SnSlW, mZ, cpl_SnSnZ, mSpm, cpl_SmpSnSl      &
          &, mP0, cpl_P0SnSn, mS0, cpl_S0SnSn                               &
          &, 1, GenerationMixing                                            &
          &, gP_Sn, gT_Sn, BR_Sn)

  gP_Sl = 0._dp
  gT_Sl = 0._dp
  BR_Sl = 0._dp
  Call SfermionTwoBodyDecays(-1, mSlepton, mf_l, mf_nu                      &
          &, mN, cpl_LNSl_L, cpl_LNSl_R, mC, cpl_CNuSl_L, cpl_CNuSl_R       &
          &, mSneut, mW, cpl_SlSnW, mZ, cpl_SlSlZ, mSpm, cpl_SmpSlSn        &
          &, mP0, cpl_P0SlSl, mS0, cpl_S0SlSl                               &
          &, 2, GenerationMixing                                            &
          &, gP_Sl, gT_Sl, BR_Sl)

  !-------------------------------------
  ! calculation of decay into gravitino
  ! relevant for GMSB
  !-------------------------------------
   If ((GenerationMixing).and.(mslepton(6).gt.m_grav)) Then
    Write(ErrCan,*) "Problem in routine "//NameOfUnit(Iname)
    Write(ErrCan,*) "Decay of a slepton into gravitino and lepton is not"
    Write(ErrCan,*) "included in case of generation mixing"
   Else ! .not.GenerationMixing
    
    Do i1=1,6
     If (mSlepton(i1).Gt.m_grav) Then
      gP_Sl(i1,13) = oo16Pi * mSlepton(i1) * (mSlepton(i1)**2-m_grav**2)**2 &
                 &   / (grav_fac * Fgmsb)**2
     Else
      gP_Sl(i1,13) = 0._dp
     End If
     If (gP_Sl(i1,13).Gt.0._dp) Then
      BR_Sl(i1,1:12) =  BR_Sl(i1,1:12) * gT_Sl(i1) / (gT_Sl(i1) + gP_Sl(i1,13))
      gT_Sl(i1) = gT_Sl(i1) + gP_Sl(i1,13)
      BR_Sl(i1,13) = gP_Sl(i1,13) / gT_Sl(i1)
     End If
    End Do
   End If ! GenerationMixing

  gP_Su = 0._dp
  gT_Su = 0._dp
  BR_Su = 0._dp
  Call SfermionTwoBodyDecays(-1, mSup, mf_u, mf_d                           &
          &, mN, cpl_UNSu_L, cpl_UNSu_R, mC, cpl_CDSu_L, cpl_CDSu_R         &
          &, mSdown, mW, cpl_SuSdW, mZ, cpl_SuSuZ, mSpm, cpl_SmpSuSd        &
          &, mP0, cpl_P0SuSu, ms0, cpl_S0SuSu                               &
          &, 0, GenerationMixing                                            &
          &, gP_Su, gT_Su, BR_Su, mGlu, cpl_GUSu_L, cpl_GUSu_R)

  gP_Sd = 0._dp
  gT_Sd = 0._dp
  BR_Sd = 0._dp
  Call SfermionTwoBodyDecays(-1, mSdown, mf_d, mf_u                         &
          &, mN, cpl_DNSd_L, cpl_DNSd_R, mC, cpl_CUSd_L, cpl_CUSd_R         &
          &, mSup, mW, cpl_SdSuW, mZ, cpl_SdSdZ, mSpm, cpl_SmpSdSu          &
          &, mP0, cpl_P0SdSd, mS0, cpl_S0SdSd                               &
          &, 0, GenerationMixing                                            &
          &, gP_Sd, gT_Sd, BR_Sd, mGlu, cpl_GDSd_L, cpl_GDSd_R)


  gP_Glu = 0._dp
  gT_Glu = 0._dp
  BR_Glu = 0._dp

  If (.Not.CTBD) Then
   Call GluinoTwoBodyDecays(mGlu, mSdown, cpl_GDSd_L, cpl_GDSd_R          &  
        & , mf_d, mSup, cpl_GUSu_L, cpl_GUSu_R, mf_u, 0, GenerationMixing &  
        & , gP_Glu(1:72), gT_Glu, BR_Glu(1:72) )

   If (.Not.GenerationMixing) Then
    gP_Glu(25) = GluinoToStopC( mglu, mSup**2, Rsup, mSdown(5:6)**2, mP0**2 &
            &    , vevSM(2)/vevSM(1), y_d(3,3), mu, A_d(3,3), cpl_GUSu_R, kont)
    gP_Glu(26) = gP_Glu(25)
   End If

   If (gT_Glu.Lt.fac3*mglu) Then
    Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .False. , gGNdd, gGNuu           &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )
    If (GenerationMixing) Then
     gP_Glu(1:72) = 0._dp
    Else
     gP_Glu(1:24) = 0._dp
     gP_Glu(27:72) = 0._dp
    End If

   Else ! calculate only 3-body modes via virtual particles

    Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .True. , gGNdd, gGNuu            &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )

   End If

  Else ! calculation of 3-body decay modes is enforced

   If (.Not.GenerationMixing) Then
    gP_Glu(25) = GluinoToStopC( mglu, mSup**2, Rsup, mSdown(5:6)**2, mP0**2 &
           &    , vevSM(2)/vevSM(1), y_d(3,3), mu, A_d(3,3), cpl_GUSu_R, kont)
    gP_Glu(26) = gP_Glu(25)
   End If
   Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .False. , gGNdd, gGNuu           &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )

  End If
  gT_Glu = Sum(gP_Glu) 
  If (gT_Glu.gt.0._dp) BR_Glu = gP_Glu / gT_Glu 
  BR_Glu2 = 0._dp
  gP_Glu2 = 0._dp
  gP_Glu2(1:77) = gP_Glu(1:77) 
  BR_Glu2(1:77) = BR_Glu(1:77) 
  BR_Glu3 = 0._dp
  gP_Glu3 = 0._dp
  gP_Glu3(1:150) = gP_Glu(78:227) 
  BR_Glu3(1:150) = BR_Glu(78:227) 

  gP_S0 = 0._dp
  gT_S0 = 0._dp
  BR_S0 = 0._dp
  vev = Sqrt(vevSM(1)**2+vevSM(2)**2)
  cpl_S0WWvirt = cpl_S0WW / vev
  cpl_S0ZZvirt = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * cpl_S0ZZ / vev
  cpl_GGS0 = 0._dp
  Call ScalarTwoBodyDecays(-1, mS0, cpl_S03, cpl_GlGlS0, cpl_GGS0            &
          &, mf_l, cpl_LLS0_L, cpl_LLS0_R, mf_d, cpl_DDS0_L, cpl_DDS0_R      & 
          &, mf_u, cpl_UUS0_L, cpl_UUS0_R, mSlepton, cpl_S0SlSl              &
          &, mSneut, cpl_S0SnSn, mSdown, cpl_S0SdSd, mSup, cpl_S0SuSu        & 
          &, mN, cpl_NNS0_L, cpl_NNS0_R, mC, cpl_CCS0_L, cpl_CCS0_R          &
          &, mW, cpl_S0WW, cpl_S0WWvirt, mZ, cpl_S0ZZ, cpl_S0ZZvirt          &
          &, mSpm, cpl_SmpS03, mP0, cpl_P0S03, cpl_P0S0Z, cpl_SmpS0W, mglu   &
          &, GenerationMixing, gP_S0, gT_S0, BR_S0)

  gP_P0 = 0._dp
  gT_P0 = 0._dp
  BR_P0 = 0._dp
  Call PseudoscalarTwoBodyDecays(2, mP0                                      &
          &, mf_l, cpl_LLP0_L, cpl_LLP0_R, mf_d, cpl_DDP0_L, cpl_DDP0_R      & 
          &, mf_u, cpl_UUP0_L, cpl_UUP0_R, mSlepton, cpl_P0SlSl              &
          &, mSdown, cpl_P0SdSd, mSup, cpl_P0SuSu                            & 
          &, mN, cpl_NNP0_L, cpl_NNP0_R, mC, cpl_CCP0_L, cpl_CCP0_R          &
          &, mSpm, cpl_SmpP03, mS0, cpl_P0S03, mZ, cpl_P0S0Z, mW, cpl_SmpP0W &
          &, mglu, GenerationMixing, gP_P0, gT_P0, BR_P0)
  gT_P0(1) = gamZ ! needed for 3-body decays

  gP_Spm = 0._dp
  gT_Spm = 0._dp
  BR_Spm = 0._dp
  Call  ChargedscalarTwoBodyDecays(2, mSpm                                   &
          &, mf_l, cpl_SmpLNu_L, cpl_SmpLNu_R                                &
          &, mf_d, mf_u, cpl_SmpDU_L, cpl_SmpDU_R                            &
          &, mSlepton, mSneut, cpl_SmpSlSn, mSdown, mSup, cpl_SmpSdSu        & 
          &, mN, mC, cpl_SmpCN_L, cpl_SmpCN_R, mW, mZ, cpl_SmpZ              &
          &, mP0, cpl_SmpP03, cpl_SmpP0W, mS0, cpl_SmpS03, cpl_SmpS0W        &
          &, 1, GenerationMixing, gP_Spm, gT_Spm, BR_Spm)
  gT_Spm(1) = gamW ! needed for 3-body decays


  gP_C = 0._dp
  gT_C = 0._dp
  BR_C = 0._dp
  Do i1=1,2
   If (.Not.CTBD)  Then
    Call CharginoTwoBodyDecays(i1, mC, mSlepton, cpl_CNuSl_L, cpl_CNuSl_R    &
        & , mSneut, cpl_CLSn_L, cpl_CLSn_R, mf_l, mSdown, cpl_CUSd_L         &
        & , cpl_CUSd_R, mf_u, mSup, cpl_CDSu_L, cpl_CDSu_R, mf_d             &
        & , mN, mW, cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R     &
        & , mZ, cpl_CCZ_L, cpl_CCZ_R, mP0, cpl_CCP0_L, cpl_CCP0_R            &
        & , mS0, cpl_CCS0_L, cpl_CCS0_R                                      &
        & , 1, GenerationMixing, gP_C(:,1:63), gT_C, BR_C(:,1:63) )
    
    If (gT_C(i1).Lt. fac3*Abs(mC(i1))) Then
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .False., gCCMajaron        &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gT_C  &
      & , gP_C(:,64:267), BR_C(:,64:267))
     gP_C(i1,1:63) = 0._dp

    Else ! calculate 3-body decays with virtual interemdiate states only
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .True., gCCMajaron         &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gT_C  &
      & , gP_C(:,64:267), BR_C(:,64:267))
    End If

   Else ! enforce calculation of 3-body final states 
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .False., gCCMajaron        &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gT_C  &
      & , gP_C(:,64:267), BR_C(:,64:267))

   End If

   gT_C(i1) = Sum(gP_C(i1,:))
   If (gT_C(i1).Gt.0._dp) BR_C(i1,:) = gP_C(i1,:) / gT_C(i1)
  End Do
  BR_C2 = 0._dp
  BR_C2(:,1:63) = BR_C(:,1:63)
  gP_C2 = 0._dp
  gP_C2(:,1:63) = gP_C(:,1:63)
  BR_C3 = 0._dp
  BR_C3(:,1:204) = BR_C(:,64:267)
  gP_C3 = 0._dp
  gP_C3(:,1:204) = gP_C(:,64:267)

  gP_N = 0._dp
  gT_N = 0._dp
  BR_N = 0._dp
  OnlySM = .True.
  If (Size(mN).Gt.5) OnlySM = .False.
  Do i1=1,5
   If (.Not.CTBD) Then
    sinW = Sqrt(sinW2)
    cosW = Sqrt(1._dp - sinW2)    
    cpl_NGP = (N(i1,1)*cosW+N(i1,2)*sinW) / (grav_fac * Fgmsb)
    cpl_NGZ = (- N(i1,1)*sinW+ N(i1,2)*cosW) / (grav_fac * Fgmsb)
    cpl_NGH = (N(i1,3) * RS0(1,1) - N(i1,4)*RS0(1,2)) / (grav_fac * Fgmsb)
    Call NeutralinoTwoBodyDecays(i1, mN, mSlepton, cpl_LNSl_L, cpl_LNSl_R   &
       & , mf_l, mSneut, cpl_NuNSn_L, cpl_NuNSn_R, mSdown, cpl_DNSd_L       &
       & , cpl_DNSd_R, mf_d, mSup, cpl_UNSu_L, cpl_UNSu_R, mf_u, mC, mW     &
       & , cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R, mZ         &
       & , cpl_NNZ_L, cpl_NNZ_R, mP0, cpl_NNP0_L, cpl_NNP0_R, mS0           &
       & , cpl_NNS0_L, cpl_NNS0_R, m_grav, cpl_NGP, cpl_NGZ, cpl_NGH, 1        &
       & , GenerationMixing, gP_N, gT_N, BR_N )


    If (gT_N(i1).Lt. fac3*Abs(mN(i1))) Then
     Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L            &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .False., 0.0_dp           &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 2500, gT_N3, gP_N(:,150:350)            &
       & , BR_N(:,150:350), .True. )
     gP_N(i1,1:149) = 0._dp

    Else ! calculate only 3-body via virtual particles
     Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L            &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .True., 0.0_dp            &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 2500, gT_N3, gP_N(:,150:350)            &
       & , BR_N(:,150:350), .True. )
    End If

   Else ! calculation of 3-body decay modes is enforced
    OnlySM = .True.
    If (Size(mN).Gt.4) OnlySM = .False.
      Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L           &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .False., 0.0_dp           &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 2500, gT_N3, gP_N(:,150:350)            &
       & , BR_N(:,150:350), .True. )

    End If ! CTBD


   gT_N(i1) = Sum(gP_N(i1,:))
   If (gT_N(i1).Gt.0._dp) BR_N(i1,:) = gP_N(i1,:) / gT_N(i1)

  End Do
  gP_N2 = 0._dp
  BR_N2 = 0._dp
  gP_N2(:,1:153) =  gP_N(:,1:153)
  BR_N2(:,1:153) =  BR_N(:,1:153)
  gP_N3 = 0._dp
  BR_N3 = 0._dp
  gP_N3(:,1:197) =  gP_N(:,154:350)
  BR_N3(:,1:197) =  BR_N(:,154:350)
  !-------------------------------------------------------
  ! in the case of the lighter stop it is possible that
  ! all two-body decay modes are kinematically forbidden 
  !-------------------------------------------------------
  mSup2 = mSup**2
  mSdown2 = mSdown**2
  mSneut2 = mSneut**2
  mSlepton2 = mSlepton**2
  mC2 = mC**2
  mN2 = mN**2
  mP02 = mP0**2
  tanb = vevSM(2) / vevSM(1)
  If ((gT_Su(5).Lt.fac3*mSup(5)).Or.(CTBD)) Then ! calculation including widths
   Call StopDecays3(mSup, mSup2, RSup, Y_d(3,3), tanb, mSdown2, gT_Sd        &
        & , A_d(3,3), mu, mP02, mN, mN2, gauge(2), cpl_UNSu_L, cpl_UNSu_R    &
        & , cpl_DNSd_L, cpl_DNSd_R, cpl_SdSuW, mSneut2, mC, mC2, gT_C        &
        & , cpl_CLSn_L, cpl_CLSn_R, cpl_CDSu_L, cpl_CDSu_R, cpl_CNW_L        &
        & , cpl_CNW_R, mSlepton2, cpl_CNuSl_R, epsI, gstCneu, gStWBNeu       &
        & , gStBSnL, gStBSlN, kont, .False.)
   gT_Su(5) =  gstWbneu + Sum(gstCneu) + Sum(gStBSnL) + Sum(gStBSlN)
   gP_Su(5,:) = 0._dp
   gP_Su(5,55:56) = gstCneu
   gP_Su(5,57) = gstWbneu
   gP_Su(5,58) = gStBSnL(3,1,1)
   gP_Su(5,59) = gStBSnL(3,2,2)
   gP_Su(5,60) = gStBSnL(3,3,3)
   gP_Su(5,61) = gStBSlN(3,1,1)
   gP_Su(5,62) = gStBSlN(3,2,1)
   gP_Su(5,63) = gStBSlN(3,3,2)
   gP_Su(5,64) = gStBSlN(3,4,2)
   gP_Su(5,65) = gStBSlN(3,5,3)
   gP_Su(5,66) = gStBSlN(3,6,3)
   BR_Su(5,:) = gP_Su(5,:) / gT_Su(5)

  Else  ! calculation excluding widths
   Call StopDecays3(mSup, mSup2, RSup, Y_d(3,3), tanb, mSdown2, gT_Sd        &
        & , A_d(3,3), mu, mP02, mN, mN2, gauge(2), cpl_UNSu_L, cpl_UNSu_R    &
        & , cpl_DNSd_L, cpl_DNSd_R, cpl_SdSuW, mSneut2, mC, mC2, gT_C        &
        & , cpl_CLSn_L, cpl_CLSn_R, cpl_CDSu_L, cpl_CDSu_R, cpl_CNW_L        &
        & , cpl_CNW_R, mSlepton2, cpl_CNuSl_R, epsI, gstCneu, gStWBNeu       &
        & , gStBSnL, gStBSlN, kont, .True.)

   gP_Su(5,55:56) = gstCneu
   gP_Su(5,57) = gstWbneu
   gP_Su(5,58) = gStBSnL(3,1,1)
   gP_Su(5,59) = gStBSnL(3,2,2)
   gP_Su(5,60) = gStBSnL(3,3,3)
   gP_Su(5,61) = gStBSlN(3,1,1)
   gP_Su(5,62) = gStBSlN(3,2,1)
   gP_Su(5,63) = gStBSlN(3,3,2)
   gP_Su(5,64) = gStBSlN(3,4,2)
   gP_Su(5,65) = gStBSlN(3,5,3)
   gP_Su(5,66) = gStBSlN(3,6,3)
   gT_Su(5) = Sum(gP_Su(5,:))
   If (gT_Su(5).Gt.0._dp) BR_Su(5,:) = gP_Su(5,:) / gT_Su(5)
   gP_Su3 = 0._dp
   BR_Su3 = 0._dp
   gP_Su3(5,1:10) = gP_Su(5,57:66) 
   BR_Su3(5,1:10) = BR_Su(5,57:66) 
  End If

  Iname = Iname - 1

 End Subroutine CalculateBR_NMSSM

 Subroutine CalculateBR_RPeps(gauge, mGlu, PhaseGlu, mC, U, V, mN, N, mSup   &
     & , RSup, mSdown, RSdown, uD_L, uD_R, uU_L, uU_R, mS0, RS0  &
     & , mP0, RP0, mSpm, RSpm, epsI, deltaM, CTBD, kont, fac3                &
     & , Y_d, A_d, Y_l, A_l, Y_u, A_u, mu, eps, vevSM, vL, Fgmsb, m32        &
     & , grav_fac, gP_Sd, gT_Sd, BR_Sd, gP_Su, gT_Su, BR_Su                  &
     & , gT_C, gP_C2, BR_C2, gP_C3, BR_C3          &
     & , gT_N, gP_N2, BR_N2, gP_N3, BR_N3, gP_Glu, gT_Glu, BR_Glu            &
     & , gP_P0, gT_P0, BR_P0, gP_S0, gT_S0, BR_S0, gP_Spm, gT_Spm, BR_Spm) 
 !------------------------------------------------------------------
 ! Calculates the branching of SUSY particles within the MSSM
 ! it is assumed that the SUSY couplings as well as the parameters
 ! are stored in the public variables in the module SusyParameters.
 ! Input: - gauge(i) .... the gauge couplings
 !        - epsI ........ precision to which the integrals for the
 !                        3-body decays are evolved
 !        - deltaM ...... maximal ratio of mass over phasespace in 
 !                        3-body decays where the masses are set 0
 !                        in the calculation of the phase space.
 !        - kont ........ control variable, is 0 if everything is o.k.
 !        - CTBD ........ logical variable, it .true. then all 3-body
 !          decays are calculated, even if 2-body decays are avaiable.
 !          Otherwisethe 3-body decays are only calculated if 2-body decays
 !          are kinematically forbidden
 !        - fac3 ........ if a total two-body decay width devided by the
 !          the correcponding mass is smaller than fac3, the width
 !          will be recalculated using the routines for 3-body decays
 !  the variable GenerationMixing is taken from the Module InputOutput
 ! The exact form of the output depends partly on the variable
 ! GenerationMixing. 
 ! output:
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  - gT_Sn(i) .... total width of sneutrino_i, i=1-3
 !  - BR_Sn(i,j) .. branching ratios of sneutrino_i
 !    GenerationMixing=.false.
 !     j=1-4 ... neutralino_j
 !     j=5,6 ... chargino 1,2
 !     j=7,8 ... W + slepton_(2*i-1,2*i)
 !     j=9,10 .. H+ + slepton_(2*i-1,2*i)
 !    GenerationMixing=.true.
 !     j=1-4 ............. neutralino_j
 !     j=5,10 ............ e+chargino_1,2 , myon+chargino_1,2, tau+chargino_1,2
 !     j=11,16 ........... W + slepton_1-6
 !     j=17,22 ........... H+ + slepton_1-6
 !     j=22+j  ........... Z + sneutrino j=1,i-1 (formally, actually =0)
 !     j=22+j,22+j*2  .... A0 + sneutrino j=1,i-1
 !     j=22+j,22+j*3  .... h0 + sneutrino j=1,i-1
 !     j=22+j,22+j*4  .... H0 + sneutrino j=1,i-1
 !  
 ! written by Werner Porod, 25.04.02
 ! 16.09.02: instead of using globally defined couplings, now locally
 !           defined couplings are used.
 ! 29.11.02: adding three-body decays of the lighter stop
 ! 09.12.02: adding decays:  neutralino_1 -> gravitino + photon
 !                           sleptons -> lepton + photon
 ! 14.08.03: adding decays:  neutralino_i -> neutralino_j + photon
 !                           gluino -> neutralino_i + gluon
 ! 21.08.03: adding new possibility concerning 2-body decay modes versus
 !           3-body decay modes: if the new variable calc_2and3_body is
 !           set true, then both possibilities will be calculated.
 !           However, in the 3-body modes on-shell particles will be
 !           negelected. This requires also a change in the output
 !           format. Moreover, also the charge conjugated final states
 !           will be printed in future. 
 !           Starting with the gluino
 ! 10.09.03: changing ordering for neutralino 2-body decays such that
 !           charge conjugated final states are included
 ! 11.09.03: merging branching ratio array for two- and three-body
 !           decay modes of neutralinos
 ! 19.01.09: m32 is given in eV, defining new variable m_grav to get
 !           value in GeV
 !------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Real(dp), Intent(in) :: epsI, deltaM, gauge(3), fac3
  Real(dp), Intent(in) :: mGlu, mC(5), mN(7)         &
         & , mSdown(6), mSup(6), mP0(5), RP0(5,5), mS0(5), RS0(5,5), mSpm(8)
  Complex(dp), Intent(in) :: PhaseGlu, RSpm(8,8), U(5,5), V(5,5), N(7,7) &
         & , RSup(6,6), RSdown(6,6), uD_L(3,3), uD_R(3,3), uU_L(3,3), uU_R(3,3)
  Complex(dp), Intent(in) :: A_d(3,3), A_l(3,3), A_u(3,3), Y_d(3,3), Y_l(3,3) &
         & , Y_u(3,3), mu, eps(3)
  Real(dp), Intent(in) :: vevSM(2), Fgmsb, m32, vL(3), grav_fac
  Logical, Intent(in) :: CTBD
  Real(dp), Intent(inout) :: gP_Sd(6,54), gT_Sd(6), BR_Sd(6,54)              &
         & , gP_Su(6,66), gT_Su(6), BR_Su(6,66)                              &
         & , gT_C(2), gP_C2(2,120), BR_C2(2,120), gP_C3(2,300), BR_C3(2,300) &
         & , gT_N(4), gP_N2(4,200), BR_N2(4,200), gP_N3(4,400), BR_N3(4,400) &
         & , gP_Glu(230), gT_Glu, BR_Glu(230)                                &
         & , gP_P0(5,200), gT_P0(5), BR_P0(5,200)                            &
         & , gP_S0(5,200), gT_S0(5), BR_S0(5,200)                            &
         & , gP_Spm(8,200), gT_Spm(8), BR_Spm(8,200)

  Integer :: i1, i2, i3, z1, z2, k_neut
  Real(dp) :: m_grav, gNNnunu(7,6,3,3)   &
    & , gNNll(7,6,3,3), gNNdd(7,6,3,3), gNNuu(7,6,3,3), gNCln(7,5,3,3)       &
    & , gNCDU(7,5,3,3), gNgdd(7,6,3), gNguu(7,6,3), gGNuu(7,3,3)             &
    & , gGNdd(7,3,3), gGgluon(7), gGCdu(5,3,3), gCuu(5,4,3,3), gCdd(5,4,3,3) &
    & , gCll(5,4,3,3), gCnunu(5,4,3,3), gCNln(5,7,3,3), gCNDU(5,7,3,3)       &
    & , gCgdu(5,3,3), gStWB(6,3), sinW, gCCMajaron(5,4), gCCCC(5,4,4,4)      &
    & , gCCNN(5,4,7,7), gNNCC(7,6,5,5), gNNNN(7,6,6,6), gNNPhoton(7,6)       &
    & , gP_C(5,400), BR_C(5,400), gTa_C(5), gP_N(7,485), BR_N(7,485), gTa_N(7)

 !----------------------------
 ! dummy couplings
 !----------------------------
 Complex(dp) :: cpl_P0SlSl(5,6,6)=0._dp, cpl_CNuSl_R(5,3,6)=0._dp &
      & , cpl_CNuSl_L(5,3,6)=0._dp, cpl_CLSn_R(5,3,3)=0._dp       &
      & , cpl_CLSn_L(5,3,3)=0._dp, cpl_LLP0_R(3,3,5)=0._dp        &
      & , cpl_LLP0_L(3,3,5)=0._dp, cpl_LLS0_R(3,3,5)=0._dp        &
      & , cpl_LLS0_L(3,3,5)=0._dp, cpl_NuNSn_R(3,7,3)=0._dp       &
      & , cpl_NuNSn_L(3,7,3)=0._dp, cpl_LNSl_R(3,7,6)=0._dp       &
      & , cpl_LNSl_L(3,7,6)=0._dp, cpl_SmpLNu_R(8,3,3)=0._dp      &
      & , cpl_SmpLNu_L(8,3,3)=0._dp, cpl_SmpSlSn(8,6,3)=0._dp     &
      & , cpl_S0SlSl(5,6,6)=0._dp, cpl_S0SnSn(5,3,3)=0._dp        &
      & , cpl_LNuW(3,3)=0._dp
 Real(dp) :: cpl_LLZ_L = 0._dp, cpl_LLZ_R = 0._dp, cpl_NuNuZ_L = 0._dp &
      & , cpl_NuNuZ_R = 0._dp
 !----------------------------
 ! couplings
 !----------------------------
  Complex(dp) :: cpl_SmpSdSu(8,6,6), cpl_SmpSuSd(8,6,6)       &
      & , cpl_SmpP03(8,8,5), cpl_SmpP0W(8,5), cpl_SmpS03(8,8,5)            &
      & , cpl_SmpS0W(8,5), cpl_SmpDU_L(8,3,3), cpl_SmpDU_R(8,3,3)          &
      & , cpl_SmpZ(8,8), cpl_DUW(3,3)
  Real(dp) :: cpl_DDZ_L = 0._dp, cpl_DDZ_R = 0._dp, cpl_UUZ_L = 0._dp      &
      & , cpl_UUZ_R = 0._dp
  Complex(dp) :: cpl_CCZ_L(5,5), cpl_CCZ_R(5,5)           &
      & , cpl_NNZ_L(7,7), cpl_NNZ_R(7,7), cpl_NNS0_L(7,7,5)            &
      & , cpl_NNS0_R(7,7,5), cpl_NNP0_L(7,7,5), cpl_NNP0_R(7,7,5) 
  Complex(dp) :: cpl_GDSd_L(3,6), cpl_GDSd_R(3,6)          &
      & , cpl_DNSd_L(3,7,6), cpl_DNSd_R(3,7,6), cpl_GUSu_L(3,6)         &
      & , cpl_GUSu_R(3,6), cpl_UNSu_L(3,7,6), cpl_UNSu_R(3,7,6)         &
      & , cpl_DDP0_L(3,3,5), cpl_UUP0_L(3,3,5), cpl_DDP0_R(3,3,5)       &
      & , cpl_UUP0_R(3,3,5), cpl_DDS0_L(3,3,5), cpl_UUS0_L(3,3,5)       &
      & , cpl_DDS0_R(3,3,5), cpl_UUS0_R(3,3,5)
  Complex(dp) :: cpl_CUSd_L(5,3,6), cpl_CUSd_R(5,3,6)      &
      & , cpl_CDSu_L(5,3,6), cpl_CDSu_R(5,3,6)
  Complex(dp) :: cpl_GlGlS0(5), cpl_GGS0(5)
  Complex(dp) :: cpl_P0SdSd(5,6,6), cpl_P0SuSu(5,6,6), cpl_P0S0Z(5,5) 
  Real(dp) :: cpl_P0S03(5,5,5)
  Complex(dp) :: cpl_S0SdSd(5,6,6), cpl_S0SuSu(5,6,6)
  Real(dp) :: cpl_S03(5,5,5), cpl_S0WW(5), cpl_S0ZZ(5), cpl_S0WWvirt(5) &
      & , cpl_S0ZZvirt(5)
  Complex(dp) :: cpl_SdSuW(6,6), cpl_SuSdW(6,6) &
      & , cpl_SdSdZ(6,6), cpl_SuSuZ(6,6), cpl_NGP, cpl_NGZ, cpl_NGH
  Complex(dp) :: cpl_CCP0_L(5,5,5), cpl_CCP0_R(5,5,5)    &
      & , cpl_CCS0_L(5,5,5), cpl_CCS0_R(5,5,5), cpl_CNW_L(5,7)        &
      & , cpl_CNW_R(5,7), cpl_SmpCN_L(8,5,7), cpl_SmpCN_R(8,5,7)

  Real(dp) :: mSup2(6), mSdown2(6), mSneut2(3), mSlepton2(6), mC2(5), mN2(7) &
     & , mP02(5), gStCNeu(2), gStWBNeu, gStBSnL(3,3,3), gStBSlN(3,6,3), tanb &
     & , sinW2, cosW, vev, gam(2)
  Real(dp) :: mSlepton(6) = 1.e18_dp, mSneut(3) =  1.e18_dp
  Complex(dp) :: bi(4)
  Real(dp), Parameter :: mf_nu(3)=0._dp
  Logical :: OnlySM

  Iname = Iname + 1
  NameOfUnit(Iname) = 'CalculateBR_RPeps'

  !----------------------------------------------------------
  ! first all couplings are calculated, see module Couplings
  !----------------------------------------------------------
  sinW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  bi(1) = mu
  bi(2:4) = eps
  m_grav = 1.e-9_dp * m32

  Call AllCouplings(gauge, Y_l, Y_d, uD_L, uD_R, Y_u, uU_L, uU_R, vevSM  &
    &  , vL, RSpm, RP0, RS0, U, V, N, bi, PhaseGlu, A_l, RSup, A_u    &
    &  , RSdown, A_d, GenerationMixing                                  &
    & , cpl_SmpSdSu, cpl_SmpSuSd, cpl_SmpP03, cpl_SmpP0W, cpl_SmpS03          &
    & , cpl_SmpS0W, cpl_SmpDU_L, cpl_SmpDU_R, cpl_SmpZ, cpl_DUW, cpl_DDZ_L    &
    & , cpl_DDZ_R, cpl_UUZ_L, cpl_UUZ_R, cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L      &
    & , cpl_NNZ_R, cpl_NNS0_L, cpl_NNS0_R, cpl_NNP0_L, cpl_NNP0_R, cpl_GDSd_L &
    & , cpl_GDSd_R, cpl_DNSd_L, cpl_DNSd_R, cpl_GUSu_L, cpl_GUSu_R            &
    & , cpl_UNSu_L, cpl_UNSu_R, cpl_DDP0_L, cpl_UUP0_L, cpl_DDP0_R            &
    & , cpl_UUP0_R, cpl_DDS0_L, cpl_UUS0_L, cpl_DDS0_R, cpl_UUS0_R            &
    & , cpl_CUSd_L, cpl_CUSd_R, cpl_CDSu_L, cpl_CDSu_R, cpl_GlGlS0            &
    & , cpl_P0SdSd, cpl_P0SuSu, cpl_P0S0Z, cpl_P0S03, cpl_S0SdSd, cpl_S0SuSu  &
    & , cpl_S03, cpl_S0WW, cpl_S0ZZ, cpl_SdSuW, cpl_SuSdW, cpl_SdSdZ          &
    & , cpl_SuSuZ, cpl_CCP0_L, cpl_CCP0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_CNW_L  &
    & , cpl_CNW_R, cpl_SmpCN_L, cpl_SmpCN_R)


  gP_Su = 0._dp
  gT_Su = 0._dp
  BR_Su = 0._dp
  Call SfermionTwoBodyDecays(-1, mSup, mf_u, mf_d                           &
          &, mN, cpl_UNSu_L, cpl_UNSu_R, mC, cpl_CDSu_L, cpl_CDSu_R         &
          &, mSdown, mW, cpl_SuSdW, mZ, cpl_SuSuZ, mSpm, cpl_SmpSuSd        &
          &, mP0, cpl_P0SuSu, ms0, cpl_S0SuSu                               &
          &, 0, GenerationMixing                                            &
          &, gP_Su, gT_Su, BR_Su, mGlu, cpl_GUSu_L, cpl_GUSu_R)

  gP_Sd = 0._dp
  gT_Sd = 0._dp
  BR_Sd = 0._dp
  Call SfermionTwoBodyDecays(-1, mSdown, mf_d, mf_u                         &
          &, mN, cpl_DNSd_L, cpl_DNSd_R, mC, cpl_CUSd_L, cpl_CUSd_R         &
          &, mSup, mW, cpl_SdSuW, mZ, cpl_SdSdZ, mSpm, cpl_SmpSdSu          &
          &, mP0, cpl_P0SdSd, mS0, cpl_S0SdSd                               &
          &, 0, GenerationMixing                                            &
          &, gP_Sd, gT_Sd, BR_Sd, mGlu, cpl_GDSd_L, cpl_GDSd_R)

  gP_Glu = 0._dp
  gT_Glu = 0._dp
  BR_Glu = 0._dp

  If (.Not.CTBD) Then
   Call GluinoTwoBodyDecays(mGlu, mSdown, cpl_GDSd_L, cpl_GDSd_R          &  
        & , mf_d, mSup, cpl_GUSu_L, cpl_GUSu_R, mf_u, 0, GenerationMixing &  
        & , gP_Glu(1:72), gT_Glu, BR_Glu(1:72) )

   If (.Not.GenerationMixing) Then
    gP_Glu(25) = GluinoToStopC( mglu, mSup**2, Rsup, mSdown(5:6)**2, mP0**2 &
            &    , vevSM(2)/vevSM(1), y_d(3,3), mu, A_d(3,3), cpl_GUSu_R, kont)
    gP_Glu(26) = gP_Glu(25)
   End If

   If (gT_Glu.Lt.fac3*mglu) Then
    Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .False. , gGNdd, gGNuu           &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )
    If (GenerationMixing) Then
     gP_Glu(1:72) = 0._dp
    Else
     gP_Glu(1:24) = 0._dp
     gP_Glu(27:72) = 0._dp
    End If

   Else ! calculate only 3-body modes via virtual particles

    Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .True. , gGNdd, gGNuu            &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )

   End If

  Else ! calculation of 3-body decay modes is enforced

   If (.Not.GenerationMixing) Then
    gP_Glu(25) = GluinoToStopC( mglu, mSup**2, Rsup, mSdown(5:6)**2, mP0**2 &
           &    , vevSM(2)/vevSM(1), y_d(3,3), mu, A_d(3,3), cpl_GUSu_R, kont)
    gP_Glu(26) = gP_Glu(25)
   End If
   Call GluinoThreeBodyDecays(mglu, mN, mC, mf_u, g_T, mf_d, mSup , gT_Su  &
       & , gauge(3), cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R         &
       & , cpl_GUSu_L, cpl_GUSu_R, cpl_SdSuW, mSdown, gT_sd, cpl_DNSd_L     &
       & , cpl_DNSd_R , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R      &
       & , GenerationMixing, epsI, deltaM, .False. , gGNdd, gGNuu           &
       & , gGCdu, gGgluon, gStWB, gT_Glu, gP_Glu(73:230), BR_Glu(73:230) )

  End If
  gT_Glu = Sum(gP_Glu) 
  If (gT_Glu.Gt.0._dp) BR_Glu = gP_Glu / gT_Glu 
  BR_Glu2 = 0._dp
  gP_Glu2 = 0._dp
  BR_Glu3 = 0._dp
  gP_Glu3 = 0._dp
  If (GenerationMixing) Then
   gP_Glu2(1:79) = gP_Glu(1:79) 
   BR_Glu2(1:79) = BR_Glu(1:79)
   gP_Glu3(1:148) = gP_Glu(79:226) 
   BR_Glu3(1:148) = BR_Glu(79:226) 
  Else
   gP_Glu2(1:26) = gP_Glu(1:26) 
   BR_Glu2(1:26) = BR_Glu(1:26)
   gP_Glu2(27:33) = gP_Glu(73:79) 
   BR_Glu2(27:33) = BR_Glu(73:79)
   gP_Glu3(1:151) = gP_Glu(80:230) 
   BR_Glu3(1:151) = BR_Glu(80:230) 
  End If

  gP_S0 = 0._dp
  gT_S0 = 0._dp
  BR_S0 = 0._dp
  vev = Sqrt(vevSM(1)**2+vevSM(2)**2)
  cpl_S0WWvirt = cpl_S0WW / vev
  cpl_S0ZZvirt = Sqrt(7._dp/12._dp-10._dp/9._dp*sinW2+40._dp/27._dp*sinW2**2) &
             & * cpl_S0ZZ / vev
  cpl_GGS0 = 0._dp
  Call ScalarTwoBodyDecays(-1, mS0, cpl_S03, cpl_GlGlS0, cpl_GGS0            &
          &, mf_l, cpl_LLS0_L, cpl_LLS0_R, mf_d, cpl_DDS0_L, cpl_DDS0_R      & 
          &, mf_u, cpl_UUS0_L, cpl_UUS0_R, mSlepton, cpl_S0SlSl              &
          &, mSneut, cpl_S0SnSn, mSdown, cpl_S0SdSd, mSup, cpl_S0SuSu        & 
          &, mN, cpl_NNS0_L, cpl_NNS0_R, mC, cpl_CCS0_L, cpl_CCS0_R          &
          &, mW, cpl_S0WW, cpl_S0WWvirt, mZ, cpl_S0ZZ, cpl_S0ZZvirt          &
          &, mSpm, cpl_SmpS03, mP0, cpl_P0S03, cpl_P0S0Z, cpl_SmpS0W, mglu   &
          &, GenerationMixing, gP_S0, gT_S0, BR_S0)

  gP_P0 = 0._dp
  gT_P0 = 0._dp
  BR_P0 = 0._dp
  Do i1=2,5
    Call PseudoscalarTwoBodyDecays(i1, mP0                                   &
          &, mf_l, cpl_LLP0_L, cpl_LLP0_R, mf_d, cpl_DDP0_L, cpl_DDP0_R      & 
          &, mf_u, cpl_UUP0_L, cpl_UUP0_R, mSlepton, cpl_P0SlSl              &
          &, mSdown, cpl_P0SdSd, mSup, cpl_P0SuSu                            & 
          &, mN, cpl_NNP0_L, cpl_NNP0_R, mC, cpl_CCP0_L, cpl_CCP0_R          &
          &, mSpm, cpl_SmpP03, mS0, cpl_P0S03, mZ, cpl_P0S0Z, mW, cpl_SmpP0W &
          &, mglu, GenerationMixing, gP_P0, gT_P0, BR_P0)
  End Do
  gT_P0(1) = gamZ ! needed for 3-body decays

  gP_Spm = 0._dp
  gT_Spm = 0._dp
  BR_Spm = 0._dp

  Do i1=2,8
    Call  ChargedscalarTwoBodyDecays(i1, mSpm                                &
          &, mf_l, cpl_SmpLNu_L, cpl_SmpLNu_R                                &
          &, mf_d, mf_u, cpl_SmpDU_L, cpl_SmpDU_R                            &
          &, mSlepton, mSneut, cpl_SmpSlSn, mSdown, mSup, cpl_SmpSdSu        & 
          &, mN, mC, cpl_SmpCN_L, cpl_SmpCN_R, mW, mZ, cpl_SmpZ              &
          &, mP0, cpl_SmpP03, cpl_SmpP0W, mS0, cpl_SmpS03, cpl_SmpS0W        &
          &, 1, GenerationMixing, gP_Spm, gT_Spm, BR_Spm)
    If (mSpm(i1).Gt.(m_grav+mf_l(3))) Then ! m_3/2 in eV given
      Do i2=1,3
       gP_Spm(i1,177+i2) = oo16Pi * mSpm(i1)**5  &
              & * (Abs(RSpm(i1,2+i2))**2+Abs(RSpm(i1,5+i2))**2) &
              & / (grav_fac * Fgmsb)**2
      end do
      gT_Spm(i1) = Sum(gP_Spm(i1,:))
      BR_Spm(i1,:) = gP_Spm(i1,:) / gT_Spm(i1)
    end if
     if ((abs(mN(4)).gt.mSpm(2)).and.(i1.le.4)) then
      z1 = 181 
      z2 = 191
      Do i2=1,3
       Do i3=i2,3
        call Slepton_to_Slepton_ll(i1, 2, i2, i3, mSpm, mN(4:7) &
                & , cpl_SmpCN_L(:,:,4:7) , cpl_SmpCN_R(:,:,4:7) , 1.e-5_dp, gam)
        gP_Spm(i1,z1) = gam(1)
        z1 = z1 + 1
        if (i2.ne.i3) gP_Spm(i1,z1) = gam(1)
        if (i2.ne.i3) z1 = z1 + 1
        gP_Spm(i1,z2) = gam(2)
        z2 = z2 + 1
       end do
      end do
      gT_Spm(i1) = Sum(gP_Spm(i1,:))
      BR_Spm(i1,:) = gP_Spm(i1,:) / gT_Spm(i1)
     end if

  End Do
  gT_Spm(1) = gamW ! needed for 3-body decays

  gT_C = 0._dp
  gTa_C = 0._dp
  gP_C = 0._dp
  BR_C = 0._dp
  gP_C2 = 0._dp
  BR_C2 = 0._dp
  gP_C3 = 0._dp
  BR_C3 = 0._dp
  Do i1=4,5
   If (.Not.CTBD)  Then
    Call CharginoTwoBodyDecays(i1, mC, mSlepton, cpl_CNuSl_L, cpl_CNuSl_R    &
        & , mSneut, cpl_CLSn_L, cpl_CLSn_R, mf_l, mSdown, cpl_CUSd_L         &
        & , cpl_CUSd_R, mf_u, mSup, cpl_CDSu_L, cpl_CDSu_R, mf_d             &
        & , mN, mW, cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R     &
        & , mZ, cpl_CCZ_L, cpl_CCZ_R, mP0, cpl_CCP0_L, cpl_CCP0_R            &
        & , mS0, cpl_CCS0_L, cpl_CCS0_R                                      &
        & , 1, GenerationMixing, gP_C(:,1:108), gTa_C, BR_C(:,1:108) )

    If (gTa_C(i1).Lt. fac3*Abs(mC(i1))) Then
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .False., gCCMajaron        &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gTa_C &
      & , gP_C(:,109:400), BR_C(:,109:400))
     gP_C(i1,1:108) = 0._dp

    Else ! calculate 3-body decays with virtual interemdiate states only
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .True., gCCMajaron         &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gTa_C &
      & , gP_C(:,109:400), BR_C(:,109:400))
    End If

   Else ! enforce calculation of 3-body final states 
     Call CharginoThreeBodyDecays(i1, mC, mZ, gamZ, cpl_NuNuZ_L, cpl_NuNuZ_R &
      & , mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L, cpl_UUZ_R, mf_d       &
      & , cpl_DDZ_L, cpl_DDZ_R, cpl_CCZ_L, cpl_CCZ_R                         &
      & , mN, mW, gamW, cpl_LNuW, cpl_DUW, cpl_NNZ_L, cpl_NNZ_R, cpl_CNW_L   &
      & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L    &
      & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L     &
      & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R         &
      & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0         &
      & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L         &
      & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup   &
      & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu        &
      & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R      &
      & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn      &
      & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl  &
      & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R                   &
      & , GenerationMixing, k_neut, epsI, deltaM, .False., gCCMajaron        &
      & , gCnunu, gCll, gCdd, gCuu, gCNln, gCNDU, gCgdu, gCCCC, gCCNN, gTa_C &
      & , gP_C(:,109:400), BR_C(:,109:400))
   End If

   gT_C(i1-3) = Sum(gP_C(i1,:))
   If (gT_C(i1-3).Gt.0._dp) BR_C(i1,:) = gP_C(i1,:) / gT_C(i1-3)
  End Do
  BR_C2 = 0._dp
  BR_C2(:,1:108) = BR_C(4:5,1:108)
  gP_C2 = 0._dp
  gP_C2(:,1:108) = gP_C(4:5,1:108)
  BR_C3 = 0._dp
  BR_C3(:,1:292) = BR_C(4:5,109:400)
  gP_C3 = 0._dp
  gP_C3(:,1:292) = gP_C(4:5,109:400)
  
  gP_N2 = 0._dp
  BR_N2 = 0._dp
  gP_N3 = 0._dp
  BR_N3 = 0._dp
  OnlySM = .True.

  Do i1=4,7
   If (.Not.CTBD) Then
    sinW = Sqrt(sinW2)
    cosW = Sqrt(1._dp - sinW2)    
    cpl_NGP = (N(i1,1)*cosW+N(i1,2)*sinW) / (grav_fac * Fgmsb)
    cpl_NGZ = (- N(i1,1)*sinW+ N(i1,2)*cosW) / (grav_fac * Fgmsb)
    cpl_NGH = (N(i1,3) * RS0(1,1) - N(i1,4)*RS0(1,2)) / (grav_fac * Fgmsb)
    Call NeutralinoTwoBodyDecays(i1, mN, mSlepton, cpl_LNSl_L, cpl_LNSl_R   &
       & , mf_l, mSneut, cpl_NuNSn_L, cpl_NuNSn_R, mSdown, cpl_DNSd_L       &
       & , cpl_DNSd_R, mf_d, mSup, cpl_UNSu_L, cpl_UNSu_R, mf_u, mC, mW     &
       & , cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R, mZ         &
       & , cpl_NNZ_L, cpl_NNZ_R, mP0, cpl_NNP0_L, cpl_NNP0_R, mS0           &
       & , cpl_NNS0_L, cpl_NNS0_R, m_grav, cpl_NGP, cpl_NGZ, cpl_NGH, 1        &
       & , GenerationMixing, gP_N(:,1:185), gTa_N, BR_N(:,1:185) )

    If (gTa_N(i1).Lt. fac3*Abs(mN(i1))) Then
     Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L            &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .False., 0.0_dp           &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 4000, gTa_N, gP_N(:,186:485)            &
       & , BR_N(:,186:485) )
     gP_N(:,1:185) = 0._dp

    Else ! calculate only 3-body via virtual particles
     Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L            &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .True., 0.0_dp            &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 4000, gTa_N, gP_N(:,186:485)            &
       & , BR_N(:,186:485) )
    End If

   Else ! calculation of 3-body decay modes is enforced
    OnlySM = .True.
    If (Size(mN).Gt.4) OnlySM = .False.
      Call NeutralinoThreeBodyDecays(i1, mN, mZ, gamZ, cpl_NuNuZ_L           &
       & , cpl_NuNuZ_R, mf_l, cpl_LLZ_L, cpl_LLZ_R, mf_u, cpl_UUZ_L          &
       & , cpl_UUZ_R, mf_d, cpl_DDZ_L, cpl_DDZ_R, cpl_NNZ_L, cpl_NNZ_R       &
       & , mC, mW, gamW, cpl_LNuW, cpl_DUW, cpl_CCZ_L, cpl_CCZ_R, cpl_CNW_L  &
       & , cpl_CNW_R, mSpm, gT_Spm, cpl_SmpCN_L, cpl_SmpCN_R, cpl_SmpLNu_L   &
       & , cpl_SmpLNu_R, cpl_SmpDU_L, cpl_SmpDU_R, mS0, gT_S0, cpl_NNS0_L    &
       & , cpl_NNS0_R, cpl_CCS0_L, cpl_CCS0_R, cpl_DDS0_L, cpl_DDS0_R        &
       & , cpl_LLS0_L, cpl_LLS0_R, cpl_UUS0_L, cpl_UUS0_R, mP0, gT_P0        &
       & , cpl_NNP0_L, cpl_NNP0_R, cpl_CCP0_L, cpl_CCP0_R, cpl_DDP0_L        &
       & , cpl_DDP0_R, cpl_LLP0_L, cpl_LLP0_R, cpl_UUP0_L, cpl_UUP0_R, mSup  &
       & , gT_Su, cpl_UNSu_L, cpl_UNSu_R, cpl_CDSu_L, cpl_CDSu_R, mGlu       &
       & , cpl_GUSu_L, cpl_GUSu_R, mSdown, gT_Sd, cpl_DNSd_L, cpl_DNSd_R     &
       & , cpl_CUSd_L, cpl_CUSd_R, cpl_GDSd_L, cpl_GDSd_R, mSneut, gT_Sn     &
       & , cpl_NuNSn_L, cpl_NuNSn_R, cpl_CLSn_L, cpl_CLSn_R, mSlepton, gT_Sl &
       & , cpl_LNSl_L, cpl_LNSl_R, cpl_CNuSl_L, cpl_CNuSl_R, gauge(2), sinW2 &
       & , GenerationMixing, OnlySM, epsI, deltaM, .False., 0.0_dp           &
       & , gNNPhoton, gNNll, gNNnunu, gNNdd, gNNuu, gNCln, gNCDU, gNgdd      &
       & , gNguu, gNNNN, GNNCC, 200, 4000, gTa_N, gP_N(:,186:485)            &
       & , BR_N(:,186:485) )

    End If ! CTBD

   !------------------------------------------
   ! decay into photons are 2-body decays
   !------------------------------------------
   gP_N2(i1-3,1:191) =  gP_N(i1,1:191)
   gP_N3(i1-3,1:294) = gP_N(i1,192:485)
   gT_N(i1-3) = Sum(gP_N2(i1-3,:)) + Sum(gP_N3(i1-3,:))

   If (gT_N(i1-3).Gt.0._dp) Then
    BR_N2(i1-3,:) = gP_N2(i1-3,:) / gT_N(i1-3)
    BR_N3(i1-3,:) = gP_N3(i1-3,:) / gT_N(i1-3)
   End If
  End Do

  !-------------------------------------------------------
  ! in the case of the lighter stop it is possible that
  ! all two-body decay modes are kinematically forbidden 
  !-------------------------------------------------------
  mSup2 = mSup**2
  mSdown2 = mSdown**2
  mSneut2 = mSneut**2
  mSlepton2 = mSlepton**2
  mC2 = mC**2
  mN2 = mN**2
  mP02 = mP0**2
  tanb = vevSM(2) / vevSM(1)
  If ((gT_Su(5).Lt.fac3*mSup(5)).Or.(CTBD)) Then ! calculation including widths
   Call StopDecays3(mSup, mSup2, RSup, Y_d(3,3), tanb, mSdown2, gT_Sd        &
        & , A_d(3,3), mu, mP02, mN, mN2, gauge(2), cpl_UNSu_L, cpl_UNSu_R    &
        & , cpl_DNSd_L, cpl_DNSd_R, cpl_SdSuW, mSneut2, mC, mC2, gT_C        &
        & , cpl_CLSn_L, cpl_CLSn_R, cpl_CDSu_L, cpl_CDSu_R, cpl_CNW_L        &
        & , cpl_CNW_R, mSlepton2, cpl_CNuSl_R, epsI, gstCneu, gStWBNeu       &
        & , gStBSnL, gStBSlN, kont, .False.)
   gT_Su(5) =  gstWbneu + Sum(gstCneu) + Sum(gStBSnL) + Sum(gStBSlN)
   gP_Su(5,:) = 0._dp
   gP_Su(5,55:56) = gstCneu
   gP_Su(5,57) = gstWbneu
   gP_Su(5,58) = gStBSnL(3,1,1)
   gP_Su(5,59) = gStBSnL(3,2,2)
   gP_Su(5,60) = gStBSnL(3,3,3)
   gP_Su(5,61) = gStBSlN(3,1,1)
   gP_Su(5,62) = gStBSlN(3,2,1)
   gP_Su(5,63) = gStBSlN(3,3,2)
   gP_Su(5,64) = gStBSlN(3,4,2)
   gP_Su(5,65) = gStBSlN(3,5,3)
   gP_Su(5,66) = gStBSlN(3,6,3)
   BR_Su(5,:) = gP_Su(5,:) / gT_Su(5)

  Else  ! calculation excluding widths
   Call StopDecays3(mSup, mSup2, RSup, Y_d(3,3), tanb, mSdown2, gT_Sd        &
        & , A_d(3,3), mu, mP02, mN, mN2, gauge(2), cpl_UNSu_L, cpl_UNSu_R    &
        & , cpl_DNSd_L, cpl_DNSd_R, cpl_SdSuW, mSneut2, mC, mC2, gT_C        &
        & , cpl_CLSn_L, cpl_CLSn_R, cpl_CDSu_L, cpl_CDSu_R, cpl_CNW_L        &
        & , cpl_CNW_R, mSlepton2, cpl_CNuSl_R, epsI, gstCneu, gStWBNeu       &
        & , gStBSnL, gStBSlN, kont, .True.)

   gP_Su(5,55:56) = gstCneu
   gP_Su(5,57) = gstWbneu
   gP_Su(5,58) = gStBSnL(3,1,1)
   gP_Su(5,59) = gStBSnL(3,2,2)
   gP_Su(5,60) = gStBSnL(3,3,3)
   gP_Su(5,61) = gStBSlN(3,1,1)
   gP_Su(5,62) = gStBSlN(3,2,1)
   gP_Su(5,63) = gStBSlN(3,3,2)
   gP_Su(5,64) = gStBSlN(3,4,2)
   gP_Su(5,65) = gStBSlN(3,5,3)
   gP_Su(5,66) = gStBSlN(3,6,3)
   gT_Su(5) = Sum(gP_Su(5,:))
   If (gT_Su(5).Gt.0._dp) BR_Su(5,:) = gP_Su(5,:) / gT_Su(5)
   gP_Su3 = 0._dp
   BR_Su3 = 0._dp
   gP_Su3(5,1:10) = gP_Su(5,57:66) 
   BR_Su3(5,1:10) = BR_Su(5,57:66) 

  End If

  Iname = Iname - 1

 End Subroutine CalculateBR_RPeps

End Module BranchingRatios

