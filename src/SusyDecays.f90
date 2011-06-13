Module SusyDecays
! comments
! In this module the routines for the decays of SUSY particles are
! stored. 

! load modules
Use Control
Use DecayFunctions
Use LoopCouplings
Use Mathematics, Only:  Li2
! load modules

Contains


 Subroutine ChargedscalarTwoBodyDecays(i_in, mSpm                            &
          &, mf_l, cpl_SmpLNu_L, cpl_SmpLNu_R                                &
          &, mf_d, mf_u, cpl_SmpDU_L, cpl_SmpDU_R                            &
          &, mSlepton, mSneutrino, cpl_SmpSlSn, mSdown, mSup, cpl_SmpSdSu    & 
          &, mN, mC, cpl_SmpCN_L, cpl_SmpCN_R, mW, mZ, cpl_SmpZ              &
          &, mP0, cpl_SmpP03, cpl_SmpP0W, mS0, cpl_SmpS03, cpl_SmpS0W        &
          &, k_neut, GenerationMixing, gP, gT, BR)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of pseudoscalars
 ! input:
 !  i_in ................. specifies the decaying pseudoscalar. The decays of
 !                         all pseudoscalars are calculated if i_in < 0.
 !  mSpm(i) .............. masses of charged scalars
 !  mf_l(i) .............. lepton masses
 !  cpl_SmpLNu_L(i,j,k) .. left coupling charged scalar - lepton - neutrino
 !  cpl_SmpLNu_R(i,j,k) .. right coupling charged scalar - lepton - neutrino
 !  mf_d(i) ............. d-quark masses
 !  mf_u(i) ............. u-quark masses
 !  cpl_SmpDU_L(i,j,k) .. left coupling charged scalar d-quark u-quark
 !  cpl_SmpDU_R(i,j,k) .. right coupling charged scalar d-quark u-quark
 !  mSlepton(i) ......... slepton masses
 !  mSneutrino(i) ....... sneutrino masses
 !  cpl_SmpSlSn(i,j,k) .. coupling charged scalar slepton sneutrino
 !  mSdown(i) ........... d-squark masses
 !  mSup(i) ............. u-squark masses
 !  cpl_SmpSdSu(i,j,k) .. coupling charged scalar d-squark u-squark
 !  mN(i) ............... neutralino masses
 !  mC(i) ............... chargino masses
 !  cpl_SmpCN_L(i,j,k) .. left charged scalar - chargino - neutralino coupling
 !  cpl_SmpCN_R(i,j,k) .. right charged scalar - chargino - neutralino coupling
 !  mW .................. mass of the W-boson
 !  mZ .................. mass of the Z-boson
 !  cpl_SmpZ(i,j) ....... coupling charged scalar - scalar - Z
 !  mP0(i) .............. pseudoscalar masses
 !  cpl_SmpP03(i,j,k) ... charged scalar - charged scalar - pseudo scalar 
 !  cpl_SmpP0W(i) ....... charged scalar - pseudo scalar - W coupling
 !  mS0(i) .............. scalar masses
 !  cpl_SmpS03(i,j,k) ... charged scalar - charged scalar - scalar coupling
 !  cpl_SmpS0W(i,j) ..... charged scalar - scalar - W coupling
 !  k_neut .............. summing over neutrinos if =1; summing over all
 !                        SM-fermions if=2
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output: 
 !  depends on the values of k_neut and GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 30.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut
  Real(dp), Intent(in) :: mSpm(:), mf_l(3), mf_d(3), mf_u(3), mSlepton(6) &
          & , mSneutrino(3), mSdown(6), mSup(6), mN(:), mC(:), mW, mZ     &
          & , mP0(:), mS0(:)
  Complex(dp), Intent(in) ::  cpl_SmpLNu_L(:,:,:), cpl_SmpLNu_R(:,:,:)    &
          & , cpl_SmpDU_L(:,:,:), cpl_SmpDU_R(:,:,:), cpl_SmpSlSn(:,:,:)  &
          & , cpl_SmpSdSu(:,:,:), cpl_SmpCN_L(:,:,:), cpl_SmpCN_R(:,:,:)  &
          & , cpl_SmpP03(:,:,:), cpl_SmpP0W(:,:), cpl_SmpZ(:,:)           &
          & , cpl_SmpS03(:,:,:), cpl_SmpS0W(:,:)
  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(inout) :: gP(:,:), gT(:)
  Real(dp), Optional, Intent(inout) :: BR(:,:)

  Integer :: i1, i2, n_neut, n_char, i_start, i_end, i_count &
         & , n_Spm, i3, n_P0, n_S0, i4
  Real(dp) :: gam, m_in, m1, m2
  Complex(dp) :: coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ChargedscalarTwoBodyDecays'


  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  If (i_in.Lt.0) Then
   i_start = 2
   i_end = n_Spm
   gT = 0._dp
   gP = 0._dp

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_Spm) ) Then 
   i_start = i_in 
   i_end = i_in
   gT(i_in) = 0._dp
   gP(i_in,:) = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_Spm) = ',i_in,n_Spm
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
   Iname = Iname - 1
   Return
  End If

  Do i1 = i_start, i_end
   m_in = mSpm(i1)
   If (m_in.Eq.0._dp) Cycle ! massless particle

   i_count = 1
   If (GenerationMixing) Then
    !---------------
    ! into leptons
    !---------------
    Do i2 = 1, 5-n_char
     Do i3 = 1, 5-n_char
      Call ScalarToTwoFermions(m_in, mf_l(i2), 0._dp, cpl_SmpLNu_L(i1,i2,i3) &
                             &, cpl_SmpLNu_R(i1,i2,i3), gam )
      If ((k_neut.Eq.1).Or.(k_neut.Eq.2) ) Then
       gP(i1, i_count) = gP(i1, i_count) + gam
      Else
       gP(i1, i_count) = gam
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If (k_neut.Eq.1) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    If (k_neut.Eq.2) Then
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End If
    !--------------
    ! into quarks
    !--------------
    Do i2 = 1, 3
     Do i3 = 1,3
      Call ScalarToTwoFermions(m_in, mf_d(i2), mf_u(i3), cpl_SmpDU_L(i1,i2,i3)&
                             &, cpl_SmpDU_R(i1,i2,i3), gam )
      If (k_neut.Eq.2) Then
       gP(i1, i_count) = gP(i1, i_count) + 3._dp * gam ! color
      Else
       gP(i1, i_count) = 3._dp * gam ! color
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    If (k_neut.Eq.2) Then
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End If
    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,2*(5-n_char)
     Do i3 = 1,5-n_char
      Call ScalarToTwoScalars(m_in, mSlepton(i2), mSneutrino(i3)   &
                             &, cpl_SmpSlSn(i1,i2,i3), gP(i1, i_count) )
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End Do
    End Do
    !------------------
    ! into squarks
    !------------------
    Do i2 = 1,6
     Do i3 = 1,6
      Call ScalarToTwoScalars(m_in, mSdown(i2), mSup(i3)   &
                             &, cpl_SmpSdSu(i1,i2,i3), gP(i1, i_count) )
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! color
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End Do
    End Do

   Else  ! GenerationMixing = .FALSE.
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1, 5 - n_char
     Call ScalarToTwoFermions(m_in, mf_l(i2), 0._dp, cpl_SmpLNu_L(i1,i2,i2) &
                            &, cpl_SmpLNu_R(i1,i2,i2), gam )
     If (k_neut.Eq.2) Then
      gP(i1, i_count) = gP(i1, i_count) + gam
     Else
      gP(i1, i_count) = gam
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    If (k_neut.Eq.2) Then
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End If
    !------------------
    ! into quarks
    !------------------
    Do i2 = 1, 3
     Call ScalarToTwoFermions(m_in, mf_d(i2), mf_u(i2), cpl_SmpDU_L(i1,i2,i2) &
                            &, cpl_SmpDU_R(i1,i2,i2), gam )
     If (k_neut.Eq.2) Then
      gP(i1, i_count) = gP(i1, i_count) + 3._dp * gam ! colour
     Else
      gP(i1, i_count) = 3._dp * gam ! colour
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    If (k_neut.Eq.2) Then
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End If
    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,(5-n_char)
     Do i3 = 1,2
      coupC = cpl_SmpSlSn(i1,(i2-1)*2+i3, i2)
      m1 = mSlepton((i2-1)*2+i3)
      m2 = mSneutrino(i2)
      Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End Do
    End Do
    !------------------
    ! into squarks
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      Do i4 = 1,2
       coupC = cpl_SmpSdSu(i1,(i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSdown((i2-1)*2+i3)
       m2 = mSup((i2-1)*2+i4)
       Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
       gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End Do
     End Do
    End Do

   End If    ! GenerationMixing

   !------------------------
   ! Charginos Neutralinos
   !------------------------
   Do i2 = 1, n_char
    Do i3 = 1, n_neut
     Call ScalarToTwoFermions(m_in, mC(i2), mN(i3), cpl_SmpCN_L(i1,i2,i3) &
                             &, cpl_SmpCN_R(i1,i2,i3), gP(i1, i_count) )
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do
   
   !-----------------------
   ! pseudoscalar + W
   !-----------------------
   Do i2 = 2,n_P0
    coupC = cpl_SmpP0W(i1, i2)
    Call ScalarToScalarVectorBoson(m_in, mP0(i2), mW, coupC, gP(i1, i_count) )
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   !-------------------------------
   ! charged scalar + pseudoscalar
   !-------------------------------
   Do i2 = 2,i1-1
    Do i3 = 2,n_P0
     coupC = cpl_SmpP03(i1, i2, i3)
     Call ScalarToTwoScalars(m_in, mSpm(i2), mP0(i3), coupC, gP(i1, i_count) )
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   !--------------
   ! scalar + W
   !--------------
   Do i2 = 1,n_S0
    coupC = cpl_SmpS0W(i1, i2)
    Call ScalarToScalarVectorBoson(m_in, mS0(i2), mW, coupC, gP(i1, i_count) )
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   !-------------------------
   ! charged scalar + scalar
   !-------------------------
   Do i2 = 2,i1-1
    Do i3 = 1,n_S0
     coupC = cpl_SmpS03(i1, i2, i3)
     Call ScalarToTwoScalars(m_in, mSpm(i2), mS0(i3), coupC, gP(i1, i_count) )
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   !----------------------
   ! charged scalar + Z
   !---------------------
   Do i2 = 2,i1-1
    coupC = cpl_SmpZ(i1, i2)
    Call ScalarToScalarVectorBoson(m_in, mSpm(i2), mZ, coupC, gP(i1, i_count) )
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   If ((Present(BR)).And.(gT(i1).Eq.0)) Then
    BR(i1,:) = 0._dp
   Else If (Present(BR)) Then
    BR(i1,:) = gP(i1,:) / gT(i1)
   End If


  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine ChargedscalarTwoBodyDecays

 Subroutine CharginoTwoBodyDecays(i_in, mC                                  &
          &, mSlepton, cpl_CNuSl_L, cpl_CNuSl_R                             &
          &, mSneutrino, cpl_CLSn_L, cpl_CLSn_R, mf_l                       &  
          &, mSdown, cpl_CUSd_L, cpl_CUSd_R, mf_u                           &  
          &, mSup, cpl_CDSu_L, cpl_CDSu_R, mf_d                             &  
          &, mN, mW, cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R   &
          &, mZ, cpl_CCZ_L, cpl_CCZ_R, mP0, cpl_CCP0_L, cpl_CCP0_R          &
          &, mS0, cpl_CCS0_L, cpl_CCS0_R                                    &
          &, k_neut, GenerationMixing                                       &
          &, gP, gT, BR )
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of charginos:
 ! input:
 !  i_in ................. specifies the decaying chargino. The decays of
 !                all charginos are calculated if n_in < 0. In the case of
 !                generation diagonal models, the charginos are ordered
 !                according to their generation
 !  mC(i) ............... chargino masses
 !  mSlepton(i) ......... slepton masses
 !  cpl_CNuSl_L(i,j,k) .. left chargino neutrino slepton coupling
 !  cpl_CNuSl_R(i,j,k) .. right chargino neutrino slepton coupling
 !  mSneutrino(i) ....... sneutrino masses
 !  cpl_CLSn_L(i,j,k) ... left chargino lepton sneutrino coupling
 !  cpl_CLSn_R(i,j,k) ... right chargino lepton sneutrino coupling
 !  mf_l(i) ............. lepton masses
 !  mSdown(i) ........... d-squark masses
 !  cpl_CUSd_L(i,j,k) ... left chargino u-quark d-squark coupling
 !  cpl_CUSd_R(i,j,k) ... right chargino u-quark d-squark coupling
 !  mf_u(i) ............. u-quark masses
 !  mSup(i) ............. u-squark masses
 !  cpl_CDSu_L(i,j,k) ... left chargino d-quark u-squark coupling
 !  cpl_CDSu_R(i,j,k) ... right chargino d-quark u-squark coupling
 !  mf_d(i) ............. d-quark masses
 !  mN(i) ............... neutralino masses
 !  mW .................. mass of the W-boson
 !  cpl_CNW_L(i,j) ...... left chargino neutralino W coupling
 !  cpl_CNW_R(i,j) ...... right chargino neutralino W coupling
 !  mSpm(i) ............. masses of charged scalars
 !  cpl_SmpCNW_L(i,j,k) . left charged scalar - chargino -neutralino coupling
 !  cpl_SmpCNW_R(i,j,k) . right charged scalar - chargino - neutralino coupling
 !  mZ .................. mass of the Z-boson
 !  cpl_CCZ_L(i,j) ...... left chargino Z coupling
 !  cpl_CCZ_R(i,j) ...... right chargino Z coupling
 !  mP0(i) .............. pseudoscalar masses
 !  cpl_CCP0_L(i,j,k) ... left chargino pseudoscalar coupling
 !  cpl_CCP0_R(i,j,k) ... right chargino pseudoscalar coupling
 !  mS0(i) .............. scalar masses
 !  cpl_CCS0_L(i,j,k) ... left chargino scalar coupling
 !  cpl_CCS0_R(i,j,k) ... right chargino scalar coupling
 !  k_neut ................ if =1 .... summing over neutrinos 
 !                          if =2 .... summing over all SM-fermions 
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output:
 !  depends on the values of k_neut and GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 26.04.2001
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut
  Real(dp), Intent(in) ::  mC(:), mSlepton(6), mSneutrino(3), mf_l(3) &
         & , mSdown(6), mf_d(3), mSup(6), mf_u(3), mN(:), mW, mSpm(:), mZ  &
         & , mP0(:), mS0(:)
  Complex(dp), Intent(in) :: cpl_CNuSl_L(:,:,:), cpl_CNuSl_R(:,:,:)        &
         & , cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:), cpl_CUSd_L(:,:,:)       &
         & , cpl_CUSd_R(:,:,:), cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)       &
         & , cpl_CNW_L(:,:), cpl_CNW_R(:,:), cpl_SmpCN_L(:,:,:)            &
         & , cpl_SmpCN_R(:,:,:), cpl_CCZ_L(:,:), cpl_CCZ_R(:,:)            &
         & , cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:), cpl_CCS0_L(:,:,:)       &
         & , cpl_CCS0_R(:,:,:)
  Real(dp), Intent(inout) :: gP(:,:), gT(:)
  Real(dp), Optional, Intent(inout) :: BR(:,:)
  Logical, Intent(in) :: GenerationMixing
 
 
  Integer :: i1, i2, n_neut, n_char, i_gen, i_start, i_end, i_count &
         & , n_Spm, i3, n_P0, n_S0
  Real(dp) :: gam, m_in 
  Complex(dp) :: coupLC, coupRC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CharginoTwoBodyDecays'


  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_char
   gT = 0._dp
   gP = 0._dp

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_char) ) Then 
   i_start = i_in 
   i_end = i_in
   gT(i_in) = 0._dp
   gP(i_in,:) = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_char) = ',i_in,n_char
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
   Iname = Iname - 1
   Return
  End If

  Do i1 = i_start, i_end
   m_in = mC(i1)
   If (Abs(m_in).Eq.0._dp) Cycle ! massless particle

   i_count = 1

   If (GenerationMixing) Then
    !----------------------------
    ! Slepton neutrino, summing over neutrinos if k_neut=1 or 2
    !----------------------------
    Do i2 = 1,2*(5-n_char) ! in case of R-parity violation
     Do i3 = 1,3
      coupLC = cpl_CNuSl_L(i1,i3,i2)
      coupRC = cpl_CNuSl_R(i1,i3,i2)
      Call  FermionToFermionScalar(m_in, 0._dp, mSlepton(i2), coupLC, coupRC &
                                 &, gam)
      If ((k_neut.Eq.1).Or.(k_neut.Eq.2)) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
      Else
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If ((k_neut.Eq.1).Or.(k_neut.Eq.2)) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    !----------------------------------------------------
    ! Sneutrino lepton, summing over leptons if k_neut=2
    !----------------------------------------------------
    Do i2 = 1,5 - n_char ! in case of R-parity violation
     Do i3 = 1,3
      coupLC = cpl_CLSn_L(i1,i3,i2)
      coupRC = cpl_CLSn_R(i1,i3,i2)
      Call  FermionToFermionScalar(m_in, mf_l(i3), mSneutrino(i2), coupLC &
                                 &, coupRC, gam)
      If (k_neut.Eq.2) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
      Else
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If (k_neut.Eq.2) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    !----------------------------------------------------
    ! u-Squark d-quark, summing over quarks if k_neut=2
    !----------------------------------------------------
    Do i2 = 1,6
     Do i3 = 1,3
      coupLC = cpl_CDSu_L(i1,i3,i2)
      coupRC = cpl_CDSu_R(i1,i3,i2)
      Call  FermionToFermionScalar(m_in, mf_d(i3), mSup(i2), coupLC, coupRC &
                                 &, gam)
      If (k_neut.Eq.2) Then
       gP(i1, i_count) = gP(i1, i_count) + 3._dp * gam  ! colour 
      Else
       gP(i1, i_count) = 3._dp * gam  ! colour 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If (k_neut.Eq.2) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    !----------------------------------------------------
    ! d-Squark u-quark, summing over quarks if k_neut=2
    !----------------------------------------------------
    Do i2 = 1,6
     Do i3 = 1,3
      coupLC = cpl_CUSd_L(i1,i3,i2)
      coupRC = cpl_CUSd_R(i1,i3,i2)
      Call  FermionToFermionScalar(m_in, mf_u(i3), mSdown(i2), coupLC, coupRC &
                                 &, gam)
      If (k_neut.Eq.2) Then
       gP(i1, i_count) = gP(i1, i_count) + 3._dp * gam  ! colour 
      Else
       gP(i1, i_count) = 3._dp * gam  ! colour 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If (k_neut.Eq.2) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
      
   Else ! GenerationMixing = .FALSE.
    !----------------------------
    ! Slepton neutrino
    !----------------------------
    Do i2 = 1,5 - n_char ! in case of R-parity violation
     Do i3 = 1,2
      i_gen = (i2-1)*2 + i3
      coupLC = cpl_CNuSl_L(i1,i2,i_gen)
      coupRC = cpl_CNuSl_R(i1,i2,i_gen)
      Call  FermionToFermionScalar(m_in, 0._dp, mSlepton(i_gen), coupLC &
                                 &, coupRC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End Do
    End Do
    !------------------
    ! Sneutrino lepton
    !------------------
    Do i2 = 1,5 - n_char ! in case of R-parity violation
     coupLC = cpl_CLSn_L(i1,i2,i2)
     coupRC = cpl_CLSn_R(i1,i2,i2)
     Call  FermionToFermionScalar(m_in, mf_l(i2), mSneutrino(i2), coupLC &
                                 &, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! d-Squark u-quark
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      i_gen = (i2-1)*2 + i3
      coupLC = cpl_CUSd_L(i1,i2,i_gen)
      coupRC = cpl_CUSd_R(i1,i2,i_gen)
      Call  FermionToFermionScalar(m_in, mf_u(i2), mSdown(i_gen), coupLC &
                                 &, coupRC, gam)
      gP(i1, i_count) = 3._dp * gam ! colour 
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End Do
    End Do
    !------------------
    ! u-Squark d-quark
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      i_gen = (i2-1)*2 + i3
      coupLC = cpl_CDSu_L(i1,i2,i_gen)
      coupRC = cpl_CDSu_R(i1,i2,i_gen)
      Call  FermionToFermionScalar(m_in, mf_d(i2), mSup(i_gen), coupLC &
                                 &, coupRC, gam)
      gP(i1, i_count) = 3._dp * gam  ! colour 
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End Do
    End Do

   End If ! GenerationMixing

   !------------------
   ! neutralino W
   !------------------
   Do i2 =1, n_neut
    coupLC = cpl_CNW_L(i1,i2)
    coupRC = cpl_CNW_R(i1,i2)
    Call FermionToFermionVectorBoson(m_in, mN(i2), mW, coupLC, coupRC, gam)
    gP(i1, i_count) = gam 
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   !---------------------------
   ! charged scalar neutralino  
   !--------------------------
   Do i2 = 2, n_Spm
    Do i3 = 1, n_neut
     coupLC = cpl_SmpCN_L(i2,i1,i3)
     coupRC = cpl_SmpCN_R(i2,i1,i3)
     Call FermionToFermionScalar(m_in, mN(i3), mSpm(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   !------------------
   ! chargino Z
   !------------------
   Do i2 =1, i1-1
    coupLC = cpl_CCZ_L(i1,i2)
    coupRC = cpl_CcZ_R(i1,i2)
    Call FermionToFermionVectorBoson(m_in, mC(i2), mZ, coupLC, coupRC, gam)
    gP(i1, i_count) = gam 
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   !------------------------
   ! pseudoscalar chargino
   !------------------------
   Do i2 = 2, n_P0
    Do i3 = 1, i1-1
     coupLC = cpl_CCP0_L(i1,i3,i2)
     coupRC = cpl_CCP0_R(i1,i3,i2)
     Call FermionToFermionScalar(m_in, mC(i3), mP0(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   !-----------------
   ! scalar chargino
   !-----------------
   Do i2 = 1, n_S0
    Do i3 = 1, i1-1
     coupLC = cpl_CCS0_L(i1,i3,i2)
     coupRC = cpl_CCS0_R(i1,i3,i2)
     Call FermionToFermionScalar(m_in, mC(i3), mS0(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   If ((Present(BR)).And.(gT(i1).Eq.0)) Then
    BR(i1,:) = 0._dp
   Else If (Present(BR)) Then
    BR(i1,:) = gP(i1,:) / gT(i1)
   End If

  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine CharginoTwoBodyDecays

 Subroutine GluinoTwoBodyDecays(mGlu, mSdown, cpl_DGSd_L, cpl_DGSd_R, mf_d  &  
                              &, mSup, cpl_UGSu_L, cpl_UGSu_R, mf_u         &  
                              &, k_neut, GenerationMixing                   &
                              &, gP, gT, BR )
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of gluinos at tree-level:
 ! input:
 !  mGlu ................ gluino masses
 !  mSdown(i) ........... d-squark masses
 !  cpl_DGSd_L(i,j,k) ... left d-quark gluino d-squark coupling
 !  cpl_DGSd_R(i,j,k) ... right d-quark gluino d-squark coupling
 !  mf_d(i) ............. d-quark masses
 !  mSup(i) ............. u-squark masses
 !  cpl_UGSu_L(i,j,k) ... left u-quark gluino u-squark coupling
 !  cpl_UGSu_R(i,j,k) ... right u-quark gluino u-squark coupling
 !  mf_u(i) ............. u-quark masses
 !  k_neut .............. if =1 .... summing over all quarks
 !  GenerationMixing .... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output:
 !  depends on the values of k_neut and GenerationMixing.
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 26.04.2001
 ! 24.08.03: up to now there has been a sum over charged conjuagted states
 !           due to the need for the Les Houches Interface this will be
 !           removed
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: k_neut
  Real(dp), Intent(in) ::  mGlu, mSdown(6), mf_d(3), mSup(6), mf_u(3)
  Complex(dp), Intent(in) :: cpl_DGSd_L(:,:), cpl_DGSd_R(:,:)   &
                          &, cpl_UGSu_L(:,:), cpl_UGSu_R(:,:)
  Real(dp), Intent(inout) :: gP(:), gT
  Real(dp), Optional, Intent(inout) :: BR(:)
  Logical, Intent(in) :: GenerationMixing
 
  Integer :: i1, i2, i_gen, i_count 
  Real(dp) :: gam
  Complex(dp) :: coupLC, coupRC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GluinoTwoBodyDecays'

  gT = 0._dp
  gP = 0._dp

  If (mglu.Eq.0._dp) Then ! massless particle
   Iname = Iname - 1
   Return
  End If

  i_count = 1

  If (GenerationMixing) Then
   !----------------------------------------------------
   ! u-Squark u-quark, summing over quarks if k_neut=1
   !----------------------------------------------------
   Do i1 = 1,6
    Do i2 = 1,3
     coupLC = cpl_UGSu_L(i2,i1)
     coupRC = cpl_UGSu_R(i2,i1)
     Call  FermionToFermionScalar(mglu, mf_u(i2), mSup(i1), coupLC, coupRC &
                                &, gam)
     If (k_neut.Eq.1) Then
      gP(i_count) = gP(i_count) +  gam
      gP(i_count+1) = gP(i_count+1) +  gam
     Else
      gP(i_count) = gam 
      gP(i_count+1) = gam 
      i_count = i_count + 2
     End If
    End Do
    If (k_neut.Eq.1)  i_count = i_count + 2
   End Do
   !----------------------------------------------------
   ! d-Squark d-quark, summing over quarks if k_neut=2
   !----------------------------------------------------
   Do i1 = 1,6
    Do i2 = 1,3
     coupLC = cpl_DGSd_L(i2,i1)
     coupRC = cpl_DGSd_R(i2,i1)
     Call  FermionToFermionScalar(mglu, mf_d(i2), mSdown(i1), coupLC, coupRC &
                                &, gam)
     If (k_neut.Eq.1) Then
      gP(i_count) = gP(i_count) + gam
      gP(i_count+1) = gP(i_count+1) +  gam
     Else
      gP(i_count) = gam
      gP(i_count+1) = gam 
      i_count = i_count + 2
     End If
    End Do
    If (k_neut.Eq.1) i_count = i_count + 2
   End Do
     
  Else ! GenerationMixing = .FALSE.
   !------------------
   ! u-Squark u-quark
   !------------------
   Do i1 = 1,3
    Do i2 = 1,2
     i_gen = (i1-1)*2 + i2
     coupLC = cpl_UGSu_L(i1,i_gen)
     coupRC = cpl_UGSu_R(i1,i_gen)
     Call  FermionToFermionScalar(mglu, mf_u(i1), mSup(i_gen), coupLC &
                                &, coupRC, gam)
     gP(i_count) = gam
     gP(i_count+1) = gam
     i_count = i_count + 2
    End Do
   End Do
   !------------------
   ! d-Squark d-quark
   !------------------
   Do i1 = 1,3
    Do i2 = 1,2
     i_gen = (i1-1)*2 + i2
     coupLC = cpl_DGSd_L(i1,i_gen)
     coupRC = cpl_DGSd_R(i1,i_gen)
     Call  FermionToFermionScalar(mglu, mf_d(i1), mSdown(i_gen), coupLC &
                                &, coupRC, gam)
     gP(i_count) = gam
     gP(i_count+1) = gam
     i_count = i_count + 2
    End Do
   End Do
  End If ! GenerationMixing

  !--------------------------------------------------------------
  ! summation over colour gives a factor 2 
  !-------------------------------------------------------------
  gP = 2._dp * gP
  gT = Sum(gP)

  If ((Present(BR)).And.(gT.Eq.0)) Then
   BR = 0._dp
  Else If (Present(BR)) Then
   BR = gP / gT
  End If

  Iname = Iname - 1

 End Subroutine GluinoTwoBodyDecays

 Subroutine NeutralinoTwoBodyDecays(i_in, mN                                &
          &, mSlepton, cpl_LNSl_L, cpl_LNSl_R, mf_l                         &
          &, mSneutrino, cpl_NuNSn_L, cpl_NuNSn_R                           &  
          &, mSdown, cpl_DNSd_L, cpl_DNSd_R, mf_d                           &  
          &, mSup, cpl_UNSu_L, cpl_UNSu_R, mf_u                             &  
          &, mC, mW, cpl_CNW_L, cpl_CNW_R, mSpm, cpl_SmpCN_L, cpl_SmpCN_R   &
          &, mZ, cpl_NNZ_L, cpl_NNZ_R, mP0, cpl_NNP0_L, cpl_NNP0_R          &
          &, mS0, cpl_NNS0_L, cpl_NNS0_R, m32, cpl_NGP, cpl_NGZ, cpl_NGH    &
          &, k_neut, GenerationMixing                                       &
          &, gP, gT, BR )
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of neutralinos:
 ! input:
 !  i_in ................. specifies the decaying neutralino. The decays of
 !                all neutralinos are calculated if n_in < 0. In the case of
 !                generation diagonal models, the neutralinos are ordered
 !                according to their generation
 !  mN(i) ............... neutralino masses
 !  mSlepton(i) ......... slepton masses
 !  cpl_LNSl_L(i,j,k) ... left lepton neutralino slepton coupling
 !  cpl_LNSl_R(i,j,k) ... right lepton neutralino slepton coupling
 !  mf_l(i) ............. lepton masses
 !  mSneutrino(i) ....... sneutrino masses
 !  cpl_NuNSn_L(i,j,k) .. left neutrino neutralino sneutrino coupling
 !  cpl_NuNSn_R(i,j,k) .. right neutrino neutralino sneutrino coupling
 !  mSdown(i) ........... d-squark masses
 !  cpl_DNSd_L(i,j,k) ... left d-quark neutralino d-squark coupling
 !  cpl_DNSd_R(i,j,k) ... right d-quark neutralino d-squark coupling
 !  mf_d(i) ............. d-quark masses
 !  mSup(i) ............. u-squark masses
 !  cpl_UNSu_L(i,j,k) ... left u-quark neutralino u-squark coupling
 !  cpl_UNSu_R(i,j,k) ... right u-quark neutralino u-squark coupling
 !  mf_u(i) ............. u-quark masses
 !  mC(i) ............... chargino masses
 !  mW .................. mass of the W-boson
 !  cpl_CNW_L(i,j) ...... left chargino neutralino W coupling
 !  cpl_CNW_R(i,j) ...... right chargino neutralino W coupling
 !  mSpm(i) ............. masses of charged scalars
 !  cpl_SmpCNW_L(i,j,k) . left charged scalar - chargino -neutralino coupling
 !  cpl_SmpCNW_R(i,j,k) . right charged scalar - chargino - neutralino coupling
 !  mZ .................. mass of the Z-boson
 !  cpl_NNZ_L(i,j) ...... left neutralino Z coupling
 !  cpl_NNZ_R(i,j) ...... right neutralino Z coupling
 !  mP0(i) .............. pseudoscalar masses
 !  cpl_NNP0_L(i,j,k) ... left neutralino pseudoscalar coupling
 !  cpl_NNP0_R(i,j,k) ... right neutralino pseudoscalar coupling
 !  mS0(i) .............. scalar masses
 !  cpl_NNS0_L(i,j,k) ... left neutralino scalar coupling
 !  cpl_NNS0_R(i,j,k) ... right neutralino scalar coupling
 !  , m32, cpl_NGP, cpl_NGZ, cpl_NGH
 !  k_neut ................ if =1 .... summing over neutrinos 
 !                          if =2 .... summing over all SM-fermions 
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output:
 !  depends on the values of k_neut and GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 26.04.2001
 ! 10.09.03: give now explicitly branching ratios of charge conjugated states
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut
  Real(dp), Intent(in) ::  mN(:), mSlepton(6), mSneutrino(3), mf_l(3) &
         & , mSdown(6), mf_d(3), mSup(6), mf_u(3), mC(:), mW, mSpm(:), mZ  &
         & , mP0(:), mS0(:), m32
  Complex(dp), Intent(in) :: cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)          &
         & , cpl_NuNSn_L(:,:,:), cpl_NuNSn_R(:,:,:), cpl_DNSd_L(:,:,:)     &
         & , cpl_DNSd_R(:,:,:), cpl_UNSu_L(:,:,:), cpl_UNSu_R(:,:,:)       &
         & , cpl_CNW_L(:,:), cpl_CNW_R(:,:), cpl_SmpCN_L(:,:,:)            &
         & , cpl_SmpCN_R(:,:,:), cpl_NNZ_L(:,:), cpl_NNZ_R(:,:)            &
         & , cpl_NNP0_L(:,:,:), cpl_NNP0_R(:,:,:), cpl_NNS0_L(:,:,:)       &
         & , cpl_NNS0_R(:,:,:), cpl_NGP, cpl_NGZ, cpl_NGH
  Real(dp), Intent(inout) :: gP(:,:), gT(:)
  Real(dp), Optional, Intent(inout) :: BR(:,:)
  Logical, Intent(in) :: GenerationMixing
 
 
  Integer :: i1, i2, n_neut, n_char, i_gen, i_start, i_end, i_count &
         & , n_Spm, i3, n_P0, n_S0
  Real(dp) :: gam, m_in 
  Complex(dp) :: coupLC, coupRC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutralinoTwoBodyDecays'

  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_neut
   gT = 0._dp
   gP = 0._dp

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_neut) ) Then 
   i_start = i_in 
   i_end = i_in
   gT(i_in) = 0._dp
   gP(i_in,:) = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_neut) = ',i_in,n_neut
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
    Iname = Iname - 1
    Return
   End If

  Do i1 = i_start, i_end
   m_in = mN(i1)
   If (Abs(m_in).Eq.0._dp) Cycle ! massless particle
   i_count = 1

   If (GenerationMixing) Then
    !----------------------------
    ! Slepton lepton, summing over neutrinos if k_neut=2
    !----------------------------
    Do i2 = 1,2*(5-n_char)
     Do i3 = 1,5 - n_char 
      coupLC = cpl_LNSl_L(i3,i1,i2)
      coupRC = cpl_LNSl_R(i3,i1,i2)
      Call  FermionToFermionScalar(m_in, mf_l(i3), mSlepton(i2), coupLC &
                                 &, coupRC, gam)
      If (k_neut.Eq.2) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
       gP(i1, i_count+1) = gP(i1, i_count+1) + gam 
      Else
       gP(i1, i_count) = gam 
       gP(i1, i_count+1) = gam 
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      End If
     End Do
     If (k_neut.Eq.2) Then
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     End If
    End Do
    !-----------------------------------------------------------
    ! Sneutrino neutrino, summing over neutrinos if k_neut=1 or 2
    !-----------------------------------------------------------
    Do i2 = 1,5 - n_char
     Do i3 = 1,5 - n_char  
      coupLC = cpl_NuNSn_L(i3,i1,i2)
      coupRC = cpl_NuNSn_R(i3,i1,i2)
      Call  FermionToFermionScalar(m_in, 0._dp, mSneutrino(i2), coupLC &
                                 &, coupRC, gam)
      If ((k_neut.Eq.1).Or.(k_neut.Eq.2)) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
       gP(i1, i_count+1) = gP(i1, i_count+1) + gam 
      Else
       gP(i1, i_count) = gam 
       gP(i1, i_count+1) = gam 
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      End If
     End Do
     If ((k_neut.Eq.1).Or.(k_neut.Eq.2)) Then
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     End If
    End Do
    !----------------------------------------------------
    ! u-Squark u-quark, summing over quarks if k_neut=2
    !----------------------------------------------------
    Do i2 = 1,6
     Do i3 = 1,3
      coupLC = cpl_UNSu_L(i3,i1,i2)
      coupRC = cpl_UNSu_R(i3,i1,i2)
      Call  FermionToFermionScalar(m_in, mf_u(i3), mSup(i2), coupLC, coupRC &
                                 &, gam)
      If (k_neut.Eq.2) Then
       gP(i1, i_count) = gP(i1, i_count) + 3._dp * gam  ! colour 
       gP(i1, i_count+1) = gP(i1, i_count+1) + 3._dp * gam  ! colour 
      Else
       gP(i1, i_count) = 3._dp * gam  ! colour 
       gP(i1, i_count+1) = 3._dp * gam  ! colour 
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      End If
     End Do
     If (k_neut.Eq.2) Then
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     End If
    End Do
    !----------------------------------------------------
    ! d-Squark d-quark, summing over quarks if k_neut=2
    !----------------------------------------------------
    Do i2 = 1,6
     Do i3 = 1,3
      coupLC = cpl_DNSd_L(i3,i1,i2)
      coupRC = cpl_DNSd_R(i3,i1,i2)
      Call  FermionToFermionScalar(m_in, mf_d(i3), mSdown(i2), coupLC, coupRC &
                                 &, gam)
      If (k_neut.Eq.2) Then
       gP(i1, i_count) = gP(i1, i_count) + 3._dp * gam  ! colour 
       gP(i1, i_count+1) = gP(i1, i_count+1) + 3._dp * gam  ! colour 
      Else
       gP(i1, i_count) = 3._dp * gam  ! colour 
       gP(i1, i_count+1) = 3._dp * gam  ! colour 
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      End If
     End Do
     If (k_neut.Eq.2) Then
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     End If
    End Do
      
   Else ! GenerationMixing = .FALSE.
    !----------------------------
    ! Slepton lepton
    !----------------------------
    Do i2 = 1,5 - n_char
     Do i3 = 1,2
      i_gen = (i2-1)*2 + i3
      coupLC = cpl_LNSl_L(i2,i1,i_gen)
      coupRC = cpl_LNSl_R(i2,i1,i_gen)
      Call  FermionToFermionScalar(m_in, mf_l(i3), mSlepton(i_gen), coupLC &
                                 &, coupRC, gam)
      gP(i1, i_count) = gam 
      gP(i1, i_count+1) = gam 
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     End Do
    End Do
    !--------------------
    ! Sneutrino neutrino
    !--------------------
    Do i2 = 1,5 - n_char
     coupLC = cpl_NuNSn_L(i2,i1,i2)
     coupRC = cpl_NuNSn_R(i2,i1,i2)
     Call  FermionToFermionScalar(m_in, 0._dp, mSneutrino(i2), coupLC &
                                 &, coupRC, gam)
      gP(i1, i_count) = gam 
      gP(i1, i_count+1) = gam 
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
    End Do
    !------------------
    ! u-Squark u-quark
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      i_gen = (i2-1)*2 + i3
      coupLC = cpl_UNSu_L(i2,i1,i_gen)
      coupRC = cpl_UNSu_R(i2,i1,i_gen)
      Call  FermionToFermionScalar(m_in, mf_u(i2), mSup(i_gen), coupLC &
                                 &, coupRC, gam)
      gP(i1, i_count) = 3._dp * gam  ! colour 
      gP(i1, i_count+1) = gP(i1, i_count)
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     End Do
    End Do
    !------------------
    ! d-Squark d-quark
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      i_gen = (i2-1)*2 + i3
      coupLC = cpl_DNSd_L(i2,i1,i_gen)
      coupRC = cpl_DNSd_R(i2,i1,i_gen)
      Call  FermionToFermionScalar(m_in, mf_d(i2), mSdown(i_gen), coupLC &
                                 &, coupRC, gam)
      gP(i1, i_count) = 3._dp * gam  ! colour 
      gP(i1, i_count+1) = gP(i1, i_count)
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     End Do
    End Do

   End If ! GenerationMixing

   !------------------
   ! chargino W
   !------------------
   Do i2 =1, n_char
    coupLC = cpl_CNW_L(i2,i1)
    coupRC = cpl_CNW_R(i2,i1)
    Call FermionToFermionVectorBoson(m_in, mC(i2), mW, coupLC, coupRC, gam)
    gP(i1, i_count) = gam 
    gP(i1, i_count+1) = gam 
    gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
    i_count = i_count + 2
   End Do

   !---------------------------
   ! charged scalar chargino
   !--------------------------
   Do i2 = 2, n_Spm
    Do i3 = 1, n_char
     coupLC = cpl_SmpCN_L(i2,i3,i1)
     coupRC = cpl_SmpCN_R(i2,i3,i1)
     Call FermionToFermionScalar(m_in, mC(i3), mSpm(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gP(i1, i_count+1) = gam 
     gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
     i_count = i_count + 2
    End Do
   End Do
   !------------------
   ! neutralino Z
   !------------------
   Do i2 =1, i1-1
    coupLC = cpl_NNZ_L(i1,i2)
    coupRC = cpl_NNZ_R(i1,i2)
    Call FermionToFermionVectorBoson(m_in, mN(i2), mZ, coupLC, coupRC, gam)
    gP(i1, i_count) = gam 
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   !------------------------
   ! pseudoscalar neutralino
   !------------------------
   Do i2 = 2, n_P0
    Do i3 = 1, i1-1
     coupLC = cpl_NNP0_L(i1,i3,i2)
     coupRC = cpl_NNP0_R(i1,i3,i2)
     Call FermionToFermionScalar(m_in, mN(i3), mP0(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   !-----------------
   ! scalar neutralino
   !-----------------
   Do i2 = 1, n_S0
    Do i3 = 1, i1-1
     coupLC = cpl_NNS0_L(i1,i3,i2)
     coupRC = cpl_NNS0_R(i1,i3,i2)
     Call FermionToFermionScalar(m_in, mN(i3), mS0(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do
   If (Abs(mN(i1)).Gt.m32) Then
    !-----------------------------------------
    ! gravitino photon
    !-----------------------------------------
    gP(i1, i_count) = oo16pi * Abs(cpl_NGP)**2 * Abs(mN(i1))**5
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
    !-----------------------------------------
    ! gravitino Z
    !-----------------------------------------
    If (Abs(mN(i1)).Gt.(m32+mZ)) Then
     gP(i1, i_count) = oo16pi * Abs(cpl_NGZ)**2 * Abs(mN(i1))  &
                     &                          * (mN(i1)**4 - mZ**4)
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End If
   !-----------------------------------------
   ! gravitino h0
   !-----------------------------------------
    If (Abs(mN(i1)).Gt.(m32+mS0(1))) Then ! h0
     gP(i1, i_count) = oo16pi * Abs(cpl_NGH)**2 * Abs(mN(i1)) &
                    &        * (mN(i1)**4 - mS0(1)**4)
     gT(i1) = gT(i1) + gP(i1, i_count)
    End If
   End If

   If ((Present(BR)).And.(gT(i1).Eq.0)) Then
    BR(i1,:) = 0._dp
   Else If (Present(BR)) Then
    BR(i1,:) = gP(i1,:) / gT(i1)
   End If

  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine NeutralinoTwoBodyDecays

 Subroutine PseudoscalarTwoBodyDecays(i_in, mP0                              &
          &, mf_l, cpl_LLP0_L, cpl_LLP0_R, mf_d, cpl_DDP0_L, cpl_DDP0_R      & 
          &, mf_u, cpl_UUP0_L, cpl_UUP0_R, mSlepton, cpl_P0SlSl              &
          &, mSdown, cpl_P0SdSd, mSup, cpl_P0SuSu                            & 
          &, mN, cpl_NNP0_L, cpl_NNP0_R, mC, cpl_CCP0_L, cpl_CCP0_R          &
          &, mSpm, cpl_SmpP03, mS0, cpl_P0S03, mZ, cpl_P0S0Z, mW, cpl_SmpP0W &
          &, mglu, GenerationMixing, gP, gT, BR)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of pseudoscalars
 ! input:
 !  i_in ................. specifies the decaying pseudoscalar. The decays of
 !                         all pseudoscalars are calculated if i_in < 0.
 !  mP0(i) .............. scalar masses
 !  mf_l(i) ............. lepton masses
 !  cpl_LLP0_L(i,j,k) ... left coupling lepton-lepton-pseudoscalar
 !  cpl_LLP0_R(i,j,k) ... right coupling lepton-lepton-pseudoscalar
 !  mf_d(i) ............. d-quark masses
 !  cpl_DDP0_L(i,j,k) ... left coupling d-quark d-quark pseudoscalar
 !  cpl_DDP0_R(i,j,k) ... right coupling d-quark d-quark pseudoscalar
 !  mf_u(i) ............. u-quark masses
 !  cpl_UUP0_L(i,j,k) ... left coupling u-quark u-quark pseudoscalar
 !  cpl_UUP0_R(i,j,k) ... right coupling u-quark u-quark pseudoscalar
 !  mSlepton(i) ......... slepton masses
 !  cpl_P0SlSl(i,j,k) ... coupling pseudoscalar slepton slepton
 !  mSdown(i) ........... d-squark masses
 !  cpl_P0SdSd(i,j,k) ... coupling pseudoscalar d-squark d-squark
 !  mSup(i) ............. u-squark masses
 !  cpl_P0SuSu(i,j,k) ... coupling pseudoscalar u-squark u-squark
 !  mN(i) ............... neutralino masses
 !  cpl_NNP0_L(i,j,k) ... left neutralino-neutralino-pseudoscalar coupling
 !  cpl_NNP0_R(i,j,k) ... right neutralino-neutralino-pseudoscalar coupling
 !  mC(i) ............... chargino masses
 !  cpl_CCP0_L(i,j,k) ... left chargino-chargino-pseudoscalar coupling
 !  cpl_CCP0_R(i,j,k) ... right chargino-chargino-pseudoscalar coupling
 !  mW .................. mass of the W-boson
 !  cpl_S0WW(i) ......... scalar-W-W coupling
 !  mZ ................... mass of the Z-boson
 !  cpl_S0ZZ(i) ......... scalar-Z-Z coupling
 !  mSpm(i) .............. masses of charged scalars
 !  cpl_SmpS03(i,j,k) .... charged scalar - charged scalar - scalar coupling
 !  mP0(i) ............... pseudoscalar masses
 !  cpl_P0S03(i,j,k) ..... pseudoscalar - pseudoscalar - scalar coupling
 !  cpl_P0S0Z(i,j) ....... pseudoscalar-scalar-Z coupling
 !  cpl_SmpS0W(i,j) ...... charged scalar - scalar - W coupling
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output: 
 !  depends on the value of  GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 30.04.2001
 ! 15.11.02: adding QCD corrections for decays into fermions
 ! 14.09.03: adding charge conjugated states to output
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mP0(:), mf_l(3), mf_d(3), mf_u(3), mSlepton(6) &
          & , mSdown(6), mSup(6), mN(:), mC(:), mW, mZ, mSpm(:), mS0(:), mglu
  Real(dp), Intent(in) ::  cpl_P0S03(:,:,:)
  Complex(dp), Intent(in) ::  cpl_LLP0_L(:,:,:)            &
          & , cpl_LLP0_R(:,:,:), cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)    &
          & , cpl_UUP0_L(:,:,:), cpl_UUP0_R(:,:,:), cpl_P0SlSl(:,:,:)    &
          & , cpl_P0SdSd(:,:,:), cpl_P0SuSu(:,:,:)    &
          & , cpl_NNP0_L(:,:,:), cpl_NNP0_R(:,:,:), cpl_CCP0_L(:,:,:)    &
          & , cpl_CCP0_R(:,:,:), cpl_SmpP03(:,:,:), cpl_SmpP0W(:,:)      &
          & , cpl_P0S0Z(:,:)
  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(inout) :: gP(:,:), gT(:)
  Real(dp), Optional, Intent(inout) :: BR(:,:)

  Integer :: i1, i2, n_neut, n_char, i_start, i_end, i_count &
         & , n_Spm, i3, n_P0, n_S0, i4
  Real(dp) :: m_in, m1, m2, alpha_3
  Complex(dp) :: coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'PseudoscalarTwoBodyDecays'

  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  If (i_in.Lt.0) Then
   i_start = 2
   i_end = n_P0
   gT = 0._dp
   gP = 0._dp

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_P0) ) Then 
   i_start = i_in 
   i_end = i_in
   gT(i_in) = 0._dp
   gP(i_in,:) = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_P0) = ',i_in,n_P0
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
   Iname = Iname - 1
   Return
  End If

  Do i1 = i_start, i_end
   m_in = mP0(i1)
   If (m_in.Eq.0._dp) Cycle ! massless particle

   alpha_3 = AlphaSDR(m_in, mGlu, mSup, mSdown) ! needed for QCD corrections
   i_count = 1
   If (GenerationMixing) Then
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1, 5 - n_char
     Do i3 = i2, 5 - n_char
      Call ScalarToTwoFermions(m_in, mf_l(i2), mf_l(i3), cpl_LLP0_L(i2,i3,i1) &
                             &, cpl_LLP0_R(i2,i3,i1), gP(i1, i_count) )
      If (i2.Ne.i3) then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      end if
     End Do
    End Do
    !------------------
    ! into d-quarks
    !------------------
    Do i2 = 1, 3
     Do i3 = i2,3
      Call ScalarToTwoFermions(m_in, mf_d(i2), mf_d(i3), cpl_DDP0_L(i2,i3,i1) &
                             &, cpl_DDP0_R(i2,i3,i1), gP(i1, i_count) )
      gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_d(i2), m_in, alpha_3)
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      end if
     End Do
    End Do
    !------------------
    ! into u-quarks
    !------------------
    Do i2 = 1, 3
     Do i3 = i2,3
      Call ScalarToTwoFermions(m_in, mf_u(i2), mf_u(i3), cpl_UUP0_L(i2,i3,i1) &
                             &, cpl_UUP0_R(i2,i3,i1), gP(i1, i_count) )
      gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_u(i2), m_in, alpha_3)
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      end if
     End Do
    End Do
    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,2*(5 - n_char)
     Do i3 = i2,2*(5 - n_char)
      Call ScalarToTwoScalars(m_in, mSlepton(i2), mSlepton(i3)   &
                             &, cpl_P0SlSl(i1,i2,i3), gP(i1, i_count) )
      If (i2.Ne.i3) then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      end if
     End Do
    End Do
    !------------------
    ! into d-squarks
    !------------------
    Do i2 = 1,6
     Do i3 = i2,6
      Call ScalarToTwoScalars(m_in, mSdown(i2), mSdown(i3)   &
                             &, cpl_P0SdSd(i1,i2,i3), gP(i1, i_count) )
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      end if
     End Do
    End Do
    !------------------
    ! into u-squarks
    !------------------
    Do i2 = 1,6
     Do i3 = i2,6
      Call ScalarToTwoScalars(m_in, mSup(i2), mSup(i3)   &
                             &, cpl_P0SuSu(i1,i2,i3), gP(i1, i_count) )
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      end if
     End Do
    End Do

   Else  ! GenerationMixing = .FALSE.
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1, 5 - n_char
     Call ScalarToTwoFermions(m_in, mf_l(i2), mf_l(i2), cpl_LLP0_L(i2,i2,i1) &
                            &, cpl_LLP0_R(i2,i2,i1), gP(i1, i_count) )
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into d-quarks
    !------------------
    Do i2 = 1, 3
     Call ScalarToTwoFermions(m_in, mf_d(i2), mf_d(i2), cpl_DDP0_L(i2,i2,i1) &
                            &, cpl_DDP0_R(i2,i2,i1), gP(i1, i_count) )
     gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_d(i2), m_in, alpha_3)
     gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into u-quarks
    !------------------
    Do i2 = 1, 3
     Call ScalarToTwoFermions(m_in, mf_u(i2), mf_u(i2), cpl_UUP0_L(i2,i2,i1) &
                            &, cpl_UUP0_R(i2,i2,i1), gP(i1, i_count) )
     gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_u(i2), m_in, alpha_3)
     gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,5 - n_char
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = cpl_P0SlSl(i1,(i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSlepton((i2-1)*2+i3)
       m2 = mSlepton((i2-1)*2+i4)
       Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
       If (i3.Eq.i4) Then
        gT(i1) = gT(i1) + gP(i1, i_count)
        i_count = i_count + 1
       Else
        gP(i1, i_count+1) = gP(i1, i_count) ! charge
        gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do
    !------------------
    ! into d-squarks
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = cpl_P0SdSd(i1,(i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSdown((i2-1)*2+i3)
       m2 = mSdown((i2-1)*2+i4)
       Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
       gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
       If (i3.Eq.i4) Then
        gT(i1) = gT(i1) + gP(i1, i_count)
        i_count = i_count + 1
       Else
        gP(i1, i_count+1) = gP(i1, i_count) ! colour + charge
        gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do
    !------------------
    ! into u-squarks
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = cpl_P0SuSu(i1,(i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSup((i2-1)*2+i3)
       m2 = mSup((i2-1)*2+i4)
       Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
       gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
       If (i3.Eq.i4) Then
        gT(i1) = gT(i1) + gP(i1, i_count)
        i_count = i_count + 1
       Else
        gP(i1, i_count+1) = gP(i1, i_count) ! colour + charge
        gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do

   End If    ! GenerationMixing

   !-------------
   ! Neutralinos
   !-------------
   Do i2 = 1, n_neut
    Do i3 = i2, n_neut
     Call ScalarToTwoFermions(m_in, mN(i2), mN(i3), cpl_NNP0_L(i2,i3, i1) &
                             &, cpl_NNP0_R(i2,i3, i1), gP(i1, i_count) )
     If (i2.Eq.i3) gP(i1, i_count) = 0.5_dp * gP(i1, i_count) ! Majorana
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do
   
   !-------------
   ! Charginos
   !-------------
   Do i2 = 1, n_char
    Do i3 = i2, n_char
     Call ScalarToTwoFermions(m_in, mC(i2), mC(i3), cpl_CCP0_L(i2,i3, i1) &
                             &, cpl_CCP0_R(i2,i3, i1), gP(i1, i_count) )
     If (i2.Ne.i3) then
      gP(i1, i_count+1) = gP(i1, i_count) ! charge
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     else
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     end if
    End Do
   End Do
   
   !--------------------
   ! charged scalar + W
   !--------------------
   Do i2 = 2,n_Spm
    coupC = cpl_SmpP0W(i2, i1)
    Call ScalarToScalarVectorBoson(m_in, mSpm(i2), mW, coupC, gP(i1, i_count) )
    gP(i1, i_count+1) = gP(i1, i_count)
    gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
    i_count = i_count + 2
   End Do

   !-------------------
   ! 2 charged scalars
   !-------------------
   Do i2 = 2,n_Spm
    Do i3 = i2,n_Spm
     coupC = cpl_SmpP03(i2, i3, i1)
     Call ScalarToTwoScalars(m_in, mSpm(i2), mSpm(i3), coupC, gP(i1, i_count) )
     If (i2.Ne.i3) then
      gP(i1, i_count+1) = gP(i1, i_count) ! charge
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     else
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     end if
    End Do
   End Do

   !-------------
   ! scalars + Z
   !-------------
   Do i2 = 1,n_S0
    coupC = cpl_P0S0Z(i1, i2)
    Call ScalarToScalarVectorBoson(m_in, mS0(i2), mZ, coupC, gP(i1, i_count) )
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   !------------------------
   ! pseudoscalar + scalar 
   !------------------------
   Do i2 = 2,i1-1
    Do i3 = 1,n_S0
     coupC = cpl_P0S03(i1, i2, i3)
     Call ScalarToTwoScalars(m_in, mP0(i2), mS0(i3), coupC, gP(i1, i_count) )
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   If ((Present(BR)).And.(gT(i1).Eq.0)) Then
    BR(i1,:) = 0._dp
   Else If (Present(BR)) Then
    BR(i1,:) = gP(i1,:) / gT(i1)
   End If


  End Do ! i1
 
  Iname = Iname - 1
 contains

  Real(dp) Function FFqcd(mf, mA, alpha_s)
  implicit none
   Real(dp) , Intent(in) :: mf, mA, alpha_s
   Real(dp) :: fac, beta, beta2, ratio, R_beta_1, Ln_R_beta_1, Ln_beta

   FFqcd = 0._dp
   ratio = mf / mA
   if (ratio.ge.0.5_dp) return ! decay is kinematically forbitten

   if (ratio.ge.0.495_dp) return ! Coloumb singularity

    beta2 = 1._dp - 4._dp * ratio**2
    beta = Sqrt(beta2)
    
    R_beta_1 = (1. - beta) / (1._dp + beta)
    Ln_beta = Log(beta)
    Ln_R_beta_1 = Log(R_beta_1)

    fac = (19._dp + 2._dp * beta2 + 3._dp * beta**4) / (16._dp * beta)      &
      &     * (-Ln_R_beta_1)                                                &
      & + 0.375_dp * (7._dp - beta2) - 3._dp * Log(4._dp/(1._dp - beta**2)) &
      & - 4._dp * Ln_beta                                                   &
      & + (1._dp + beta**2)                                                 &
      &       * ( 4._dp * Li2(R_beta_1) + 2._dp * Li2(- R_beta_1)           &
      &         + Ln_R_beta_1 * ( 3._dp * Log(2._dp/(1._dp + beta))         &
      &                         + 2._dp * Ln_beta )  ) / beta
    fac =  fac - 3._dp * Log(ratio)  ! absorb large logarithms in mass

    FFqcd = 1._dp + 5._dp * alpha_s * fac * oo3pi 

  end  Function FFqcd

 End Subroutine PseudoscalarTwoBodyDecays

 Subroutine ScalarTwoBodyDecays(i_in, mS0, cpl_S03, cpl_GlGlS0, cpl_GGS0     &
          &, mf_l, cpl_LLS0_L, cpl_LLS0_R, mf_d, cpl_DDS0_L, cpl_DDS0_R      & 
          &, mf_u, cpl_UUS0_L, cpl_UUS0_R, mSlepton, cpl_S0SlSl              &
          &, mSneutrino, cpl_S0SnSn, mSdown, cpl_S0SdSd, mSup, cpl_S0SuSu    & 
          &, mN, cpl_NNS0_L, cpl_NNS0_R, mC, cpl_CCS0_L, cpl_CCS0_R          &
          &, mW, cpl_S0WW, cpl_S0WWvirt, mZ, cpl_S0ZZ, cpl_S0ZZvirt          &
          &, mSpm, cpl_SmpS03, mP0, cpl_P0S03, cpl_P0S0Z, cpl_SmpS0W, mglu   &
          &, GenerationMixing, gP, gT, BR)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of scalars
 ! input:
 !  i_in ................. specifies the decaying scalar. The decays of
 !                         all scalarss are calculated if i_in < 0.
 !  mS0(i) .............. scalar masses
 !  cpl_S03(i,j,k) ...... trilinear self interaction of scalars
 !  cpl_GlGlS0(i1) ...... gluon-gluon-scalar coupling
 !  cpl_GGS0(i1) ........ photon-photon-scalar coupling
 !  mf_l(i) ............. lepton masses
 !  cpl_LLS0_L(i,j,k) ... left coupling lepton-lepton-scalar
 !  cpl_LLS0_R(i,j,k) ... right coupling lepton-lepton-scalar
 !  mf_d(i) ............. d-quark masses
 !  cpl_DDS0_L(i,j,k) ... left coupling d-quark d-quark scalar
 !  cpl_DDS0_R(i,j,k) ... right coupling d-quark d-quark scalar
 !  mf_u(i) ............. u-quark masses
 !  cpl_UUS0_L(i,j,k) ... left coupling u-quark u-quark scalar
 !  cpl_UUS0_R(i,j,k) ... right coupling u-quark u-quark scalar
 !  mSlepton(i) ......... slepton masses
 !  cpl_S0SlSl(i,j,k) ... coupling scalar slepton slepton
 !  mSneutrino(i) ....... sneutrino masses
 !  cpl_S0SnSn(i,j,k) ... coupling scalar sneutrino sneutrino
 !  mSdown(i) ........... d-squark masses
 !  cpl_S0SdSd(i,j,k) ... coupling scalar d-squark d-squark
 !  mSup(i) ............. u-squark masses
 !  cpl_S0SuSu(i,j,k) ... coupling scalar u-squark u-squark
 !  mN(i) ............... neutralino masses
 !  cpl_NNS0_L(i,j,k) ... left neutralino-neutralino-scalar coupling
 !  cpl_NNS0_R(i,j,k) ... right neutralino-neutralino-scalar coupling
 !  mC(i) ............... chargino masses
 !  cpl_CCS0_L(i,j,k) ... left chargino-chargino-scalar coupling
 !  cpl_CCS0_R(i,j,k) ... right chargino-chargino-scalar coupling
 !  mW .................. mass of the W-boson
 !  cpl_S0WW(i) ......... scalar-W-W coupling
 !  mZ ................... mass of the Z-boson
 !  cpl_S0ZZ(i) ......... scalar-Z-Z coupling
 !  mSpm(i) .............. masses of charged scalars
 !  cpl_SmpS03(i,j,k) .... charged scalar - charged scalar - scalar coupling
 !  mP0(i) ............... pseudoscalar masses
 !  cpl_P0S03(i,j,k) ..... pseudoscalar - pseudoscalar - scalar coupling
 !  cpl_P0S0Z(i,j) ....... pseudoscalar-scalar-Z coupling
 !  cpl_SmpS0W(i,j) ...... charged scalar - scalar - W coupling
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 ! output: 
 !  depends on the value of  GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  gP(:,:) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 30.04.2001
 ! 15.11.02: adding QCD corrections for decays into fermions
 ! 14.09.03: adding charge conjugated states to output
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mS0(:), mf_l(3), mf_d(3), mf_u(3), mSlepton(6) &
          & , mSneutrino(3), mSdown(6), mSup(6), mN(:), mC(:), mW, mZ    &
          & , mSpm(:), mP0(:), mglu
  Real(dp), Intent(in) :: cpl_S03(:,:,:), cpl_S0WW(:), cpl_S0ZZ(:)       &
          & , cpl_P0S03(:,:,:), cpl_S0WWvirt(:), cpl_S0ZZvirt(:)
  Complex(dp), Intent(in) :: cpl_GlGlS0(:), cpl_GGS0(:), cpl_LLS0_L(:,:,:) &
          & , cpl_LLS0_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_R(:,:,:)      &
          & , cpl_UUS0_L(:,:,:), cpl_UUS0_R(:,:,:), cpl_S0SlSl(:,:,:)      &
          & , cpl_S0SnSn(:,:,:), cpl_S0SdSd(:,:,:), cpl_S0SuSu(:,:,:)      &
          & , cpl_NNS0_L(:,:,:), cpl_NNS0_R(:,:,:), cpl_CCS0_L(:,:,:)      &
          & , cpl_CCS0_R(:,:,:), cpl_SmpS03(:,:,:), cpl_SmpS0W(:,:)        &
          & , cpl_P0S0Z(:,:)
  Logical, Intent(in) :: GenerationMixing

  Real(dp), Intent(inout) :: gP(:,:), gT(:)
  Real(dp), Optional, Intent(inout) :: BR(:,:)

  Integer :: i1, i2, n_neut, n_char, i_start, i_end, i_count &
         & , n_Spm, i3, n_P0, n_S0, i4
  Real(dp) :: m_in, m1, m2, alpha_3
  Complex(dp) :: coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ScalarTwoBodyDecays'

  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_S0
   gT = 0._dp
   gP = 0._dp

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_S0) ) Then 
   i_start = i_in 
   i_end = i_in
   gT(i_in) = 0._dp
   gP(i_in,:) = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_S0) = ',i_in,n_S0
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
   Iname = Iname - 1
   Return
  End If

  Do i1 = i_start, i_end
   m_in = mS0(i1)
   If (m_in.Eq.0._dp) Cycle ! massless particle

   alpha_3 = AlphaSDR(m_in, mGlu, mSup, mSdown) ! needed for QCD corrections
   i_count = 1
   If (GenerationMixing) Then
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1,5-n_char 
     Do i3 = i2,5-n_char
      Call ScalarToTwoFermions(m_in, mf_l(i2), mf_l(i3), cpl_LLS0_L(i2,i3,i1) &
                             &, cpl_LLS0_R(i2,i3,i1), gP(i1, i_count) )
      If (i2.Ne.i3) Then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      Else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into d-quarks
    !------------------
    Do i2 = 1, 3
     Do i3 = i2,3
      Call ScalarToTwoFermions(m_in, mf_d(i2), mf_d(i3), cpl_DDS0_L(i2,i3,i1) &
                             &, cpl_DDS0_R(i2,i3,i1), gP(i1, i_count) )
      gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_d(i2), m_in, alpha_3)
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) Then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      Else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into u-quarks
    !------------------
    Do i2 = 1, 3
     Do i3 = i2,3
      Call ScalarToTwoFermions(m_in, mf_u(i2), mf_u(i3), cpl_UUS0_L(i2,i3,i1) &
                             &, cpl_UUS0_R(i2,i3,i1), gP(i1, i_count) )
      gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_u(i2), m_in, alpha_3)
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) Then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      Else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,2*(5-n_char)
     Do i3 = i2,2*(5-n_char)
      Call ScalarToTwoScalars(m_in, mSlepton(i2), mSlepton(i3)   &
                             &, cpl_S0SlSl(i1,i2,i3), gP(i1, i_count) )
      If (i2.Ne.i3) Then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      Else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into sneutrinos
    !------------------
    Do i2 = 1, 5-n_char
     Do i3 = i2,5-n_char
      Call ScalarToTwoScalars(m_in, mSneutrino(i2), mSneutrino(i3)   &
                             &, cpl_S0SnSn(i1,i2,i3), gP(i1, i_count) )
      If (i2.Ne.i3) Then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      Else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into d-squarks
    !------------------
    Do i2 = 1,6
     Do i3 = i2,6
      Call ScalarToTwoScalars(m_in, mSdown(i2), mSdown(i3)   &
                             &, cpl_S0SdSd(i1,i2,i3), gP(i1, i_count) )
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) Then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      Else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do
    !------------------
    ! into u-squarks
    !------------------
    Do i2 = 1,6
     Do i3 = i2,6
      Call ScalarToTwoScalars(m_in, mSup(i2), mSup(i3)   &
                             &, cpl_S0SuSu(i1,i2,i3), gP(i1, i_count) )
      gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
      If (i2.Ne.i3) Then
       gP(i1, i_count+1) = gP(i1, i_count) ! charge
       gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
       i_count = i_count + 2
      Else
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
    End Do

   Else  ! GenerationMixing = .FALSE.
    !------------------
    ! into leptons
    !------------------
    Do i2 = 1, 5-n_char
     Call ScalarToTwoFermions(m_in, mf_l(i2), mf_l(i2), cpl_LLS0_L(i2,i2,i1) &
                            &, cpl_LLS0_R(i2,i2,i1), gP(i1, i_count) )
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into d-quarks
    !------------------
    Do i2 = 1, 3
     Call ScalarToTwoFermions(m_in, mf_d(i2), mf_d(i2), cpl_DDS0_L(i2,i2,i1) &
                            &, cpl_DDS0_R(i2,i2,i1), gP(i1, i_count) )
     gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_d(i2), m_in, alpha_3)
     gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into u-quarks
    !------------------
    Do i2 = 1, 3
     Call ScalarToTwoFermions(m_in, mf_u(i2), mf_u(i2), cpl_UUS0_L(i2,i2,i1) &
                            &, cpl_UUS0_R(i2,i2,i1), gP(i1, i_count) )
     gP(i1, i_count) = gP(i1, i_count) * FFqcd(mf_u(i2), m_in, alpha_3)
     gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into sleptons
    !------------------
    Do i2 = 1,5-n_char
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = cpl_S0SlSl(i1,(i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSlepton((i2-1)*2+i3)
       m2 = mSlepton((i2-1)*2+i4)
       Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
       If (i3.Eq.i4) Then
        gT(i1) = gT(i1) + gP(i1, i_count)
        i_count = i_count + 1
       Else
        gP(i1, i_count+1) = gP(i1, i_count) ! charge
        gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do
    !------------------
    ! into sneutrinos
    !------------------
    Do i2 = 1, 5-n_char
     Call ScalarToTwoScalars(m_in, mSneutrino(i2), mSneutrino(i2)   &
                            &, cpl_S0SnSn(i1,i2,i2), gP(i1, i_count) )
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
    !------------------
    ! into d-squarks
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = cpl_S0SdSd(i1,(i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSdown((i2-1)*2+i3)
       m2 = mSdown((i2-1)*2+i4)
       Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
       gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
       If (i3.Eq.i4) Then
        gT(i1) = gT(i1) + gP(i1, i_count)
        i_count = i_count + 1
       Else
        gP(i1, i_count+1) = gP(i1, i_count) ! colour + charge
        gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do
    !------------------
    ! into u-squarks
    !------------------
    Do i2 = 1,3
     Do i3 = 1,2
      Do i4 = i3,2
       coupC = cpl_S0SuSu(i1,(i2-1)*2+i3, (i2-1)*2+i4)
       m1 = mSup((i2-1)*2+i3)
       m2 = mSup((i2-1)*2+i4)
       Call ScalarToTwoScalars(m_in, m1, m2, coupC, gP(i1, i_count) )
       gP(i1, i_count) = 3._dp * gP(i1, i_count) ! colour
       If (i3.Eq.i4) Then
        gT(i1) = gT(i1) + gP(i1, i_count)
        i_count = i_count + 1
       Else
        gP(i1, i_count+1) = gP(i1, i_count) ! colour + charge
        gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
        i_count = i_count + 2
       End If
      End Do
     End Do
    End Do

   End If    ! GenerationMixing

   !-------------
   ! Neutralinos
   !-------------
   Do i2 = 1, n_neut
    Do i3 = i2, n_neut
     Call ScalarToTwoFermions(m_in, mN(i2), mN(i3), cpl_NNS0_L(i2,i3, i1) &
                             &, cpl_NNS0_R(i2,i3, i1), gP(i1, i_count) )
     If (i2.Eq.i3) gP(i1, i_count) = 0.5_dp * gP(i1, i_count) ! Majorana
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do
   
   !-------------
   ! Charginos
   !-------------
   Do i2 = 1, n_char
    Do i3 = i2, n_char
     Call ScalarToTwoFermions(m_in, mC(i2), mC(i3), cpl_CCS0_L(i2,i3, i1) &
                             &, cpl_CCS0_R(i2,i3, i1), gP(i1, i_count) )
     If (i2.Ne.i3) Then
      gP(i1, i_count+1) = gP(i1, i_count) ! charge
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     Else
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
   End Do
   
   !------
   ! Z Z
   !------
   coupC = cpl_S0ZZ(i1)
   Call ScalarToTwoVectorBosons(m_in, mZ, coupC, gP(i1, i_count) )
   gP(i1, i_count) = 0.5_dp * gP(i1, i_count) ! identical particles 
   gT(i1) = gT(i1) + gP(i1, i_count)
   i_count = i_count + 1

   !------
   ! W W
   !------
   coupC = cpl_S0WW(i1)
   Call ScalarToTwoVectorBosons(m_in, mW, coupC, gP(i1, i_count) )
   gT(i1) = gT(i1) + gP(i1, i_count)
   i_count = i_count + 1

   !-------------------
   ! pseudoscalars + Z
   !-------------------
   Do i2 = 2,n_P0
    coupC = cpl_P0S0Z(i2, i1)
    Call ScalarToScalarVectorBoson(m_in, mP0(i2), mZ, coupC, gP(i1, i_count) )
    gT(i1) = gT(i1) + gP(i1, i_count)
    i_count = i_count + 1
   End Do

   !-------------------
   ! 2 pseudoscalars 
   !-------------------
   Do i2 = 2,n_P0
    Do i3 = i2,n_P0
     coupC = cpl_P0S03(i2, i3, i1)
     Call ScalarToTwoScalars(m_in, mP0(i2), mP0(i3), coupC, gP(i1, i_count) )
     If (i2.Eq.i3) gP(i1,i_count)= 0.5_dp*gP(i1,i_count) ! identical particles 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   !-------------------
   ! 2 scalars 
   !-------------------
   Do i2 = 1,i1-1
    Do i3 = i2,i1-1
     coupC = cpl_S03(i1, i2, i3)
     Call ScalarToTwoScalars(m_in, mS0(i2), mS0(i3), coupC, gP(i1, i_count) )
     If (i2.Eq.i3) gP(i1,i_count)= 0.5_dp*gP(i1,i_count) ! identical particles 
     gT(i1) = gT(i1) + gP(i1, i_count)
     i_count = i_count + 1
    End Do
   End Do

   !--------------------
   ! charged scalar + W
   !--------------------
   Do i2 = 2,n_Spm
    coupC = cpl_SmpS0W(i2, i1)
    Call ScalarToScalarVectorBoson(m_in, mSpm(i2), mW, coupC, gP(i1, i_count) )
    gP(i1, i_count+1) = gP(i1, i_count)
    gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
    i_count = i_count + 2
   End Do

   !-------------------
   ! 2 charged scalars
   !-------------------
   Do i2 = 2,n_Spm
    Do i3 = i2,n_Spm
     coupC = cpl_SmpS03(i2, i3, i1)
     Call ScalarToTwoScalars(m_in, mSpm(i2), mSpm(i3), coupC, gP(i1, i_count) )
     If (i2.Ne.i3) Then
      gP(i1, i_count+1) = gP(i1, i_count) ! charge
      gT(i1) = gT(i1) + 2._dp * gP(i1, i_count)
      i_count = i_count + 2
     Else
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
   End Do

   !--------------
   ! two gluons
   !--------------
   gP(i1, i_count) = G_F * m_in**3 * oosqrt2 * oo36pi3 * Abs(cpl_GlGlS0(i1))**2
   gT(i1) = gT(i1) + gP(i1, i_count)
   i_count = i_count + 1
   !--------------
   ! two photons
   !--------------
   gP(i1, i_count) = G_F * m_in**3 * oosqrt2 * oo128pi3 * Abs(cpl_GGS0(i1))**2
   gT(i1) = gT(i1) + gP(i1, i_count)
   i_count = i_count + 1

   !-----------------------
   ! W W^*
   !-----------------------
   Call  ScalarToVectorbosonsVR(m_in, mW, cpl_S0WWvirt(i1), gP(i1, i_count) )
   gT(i1) = gT(i1) + gP(i1, i_count)
   gP(i1, i_count) = 0.5_dp * gP(i1, i_count) ! formula is for sum over charges
   gP(i1, i_count + 1) = gP(i1, i_count)
   i_count = i_count + 2

   !-----------------------
   ! Z Z^*
   !-----------------------
   Call  ScalarToVectorbosonsVR(m_in, mZ, cpl_S0ZZvirt(i1), gP(i1, i_count) )
   gT(i1) = gT(i1) + gP(i1, i_count)
   i_count = i_count + 1


   If ((Present(BR)).And.(gT(i1).Eq.0)) Then
    BR(i1,:) = 0._dp
   Else If (Present(BR)) Then
    BR(i1,:) = gP(i1,:) / gT(i1)
   End If


  End Do ! i1
 
  Iname = Iname - 1

 Contains

  Real(dp) Function FFqcd(mf, mA, alpha_s)
  Implicit None
   Real(dp) , Intent(in) :: mf, mA, alpha_s
   Real(dp) :: fac, beta, beta2, ratio, R_beta_1, Ln_R_beta_1, Ln_beta

   FFqcd = 0._dp
   ratio = mf / mA
   If (ratio.Ge.0.5_dp) Return ! decay is kinematically forbitten

   If (ratio.Ge.0.495_dp) Return ! Coloumb singularity

    beta2 = 1._dp - 4._dp * ratio**2
    beta = Sqrt(beta2)
    
    R_beta_1 = (1. - beta) / (1._dp + beta)
    Ln_beta = Log(beta)
    Ln_R_beta_1 = Log(R_beta_1)

    fac = (3._dp + 34._dp * beta2 - 13._dp * beta**4) / (16._dp * beta**3)  &
      &     * (-Ln_R_beta_1)                                                &
      & + 0.375_dp * (7._dp - beta2) - 3._dp * Log(4._dp/(1._dp - beta**2)) &
      & - 4._dp * Ln_beta                                                   &
      & + (1._dp + beta**2)                                                 &
      &       * ( 4._dp * Li2(R_beta_1) + 2._dp * Li2(- R_beta_1)           &
      &         + Ln_R_beta_1 * ( 3._dp * Log(2._dp/(1._dp + beta))         &
      &                         + 2._dp * Ln_beta )  ) / beta

    fac = fac - 3._dp * Log(ratio)  ! absorb large logarithms in mass

    FFqcd = 1._dp + 5._dp * alpha_s * fac * oo3pi 

  End  Function FFqcd

 End Subroutine ScalarTwoBodyDecays


 Subroutine SfermionTwoBodyDecays(i_in, mSf, mf, mfp                    &
          &, mN, cpl_FNSf_L, cpl_FNSf_R, mC, cpl_CFpSf_L, cpl_CFpSf_R   &
          &, mSfp, mW, cpl_SfSfpW, mZ, cpl_SfSfZ, mSpm, cpl_SmpSfSfp    &
          &, mP0, cpl_P0SfSf, mS0, cpl_S0SfSf                           &
          &, k_neut, GenerationMixing                                   &
          &, gP, gT, BR                                                 &
          &, mG, cpl_GQSq_L, cpl_GQSq_R)
 !-----------------------------------------------------------------------
 ! Calculates the 2-body decays of sfermions:
 ! input:
 !  i_in ................. specifies the decaying sfermion. The decays of
 !                all sfermions are calculated if n_in < 0. In the case of
 !                generation diagonal models, the sfermions are ordered
 !                according to their generation
 !  mSf(i) .............. sfermion masses
 !  mf(i) ............... the corresponding fermion masses
 !  mfp(i) .............. the corresponding fermion' masses
 !  mN(i) ............... neutralino masses
 !  cpl_FNSf_L(i,j,k) ... left fermion-neutralino-sfermion coupling
 !  cpl_FNSf_R(i,j,k) ... right fermion-neutralino-sfermion coupling
 !  mC(i) ............... chargino masses
 !  cpl_CFpSf_L(i,j,k) ... left chargino fermion' sfermion coupling
 !  cpl_CFpSf_R(i,j,k) ... right chargino fermion' sfermion coupling
 !  mSfp(i) .............. the corresponding sfermion' masses
 !  mW ................... mass of the W-boson
 !  cpl_SfSfpW(i,j) ...... coupling sfermion-sfermion'-W
 !  mZ ................... mass of the Z-boson
 !  cpl_SfSfZ(i,j) ....... coupling sfermion-sfermion-Z
 !  mSpm(i) .............. masses of charged scalars
 !  cpl_SmpSfSfp(i,j,k) .. charged scalar - sfermion - sfermion' coupling
 !  mP0(i) ............... pseudoscalar masses
 !  cpl_P0SfSf(i,j,k) .... pseudoscalar - sfermion - sfermion coupling
 !  mS0(i) ............... scalar masses
 !  cpl_S0SfSf(i,j,k) .... scalar - sfermion - sfermion coupling
 !  k_neut ................ if =1 .... summing over fermions in the neutralino
 !                                     final states
 !                          if =2 .... summing over fermions in the chargino
 !                                     final states
 !                          if =3 .... summing over fermions in the chargino
 !                                     and neutralino final states
 !  GenerationMixing ..... mixing between the generations is taken into 
 !                         account if =.TRUE. 
 !  mG ................... Gluino mass, optional
 !  cpl_GQSq_L(i,j) ... left gluino quark squark coupling, optional
 !  cpl_GQSq_R(i,j) ... right gluino quark squark coupling, optional
 ! output: 
 !  depends on the values of k_neut and GenerationMixing and also on the
 !  lengths of mN, mC, mSpm, mP0, and mS0 which are measured by n_neut,
 !  n_char, n_Spm, n_P0 and n_S0, respectively (inside the subroutine).
 !  In addition the variables n_sfer (n_sferp) give the lengths of the
 !  sfermion (sfermions') depending on the sfermion:
 !    n_sfer, n_sferp = 3 for sneutrinos and 6 otherwise 
 !  gP(i,j) ...... partial widths
 !  gamT(:) ...... total width
 !  BR(:,:) the corresponding branching ratios, optional
 ! written by Werner Porod, 16.04.2001
 !  19.04.2001: adding interface for decay into gluinos
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in, k_neut
  Real(dp), Intent(in) :: mSf(:), mf(3), mfp(3), mN(:), mC(:), mSfp(:) &
               & , mW, mZ, mSpm(:), mP0(:), mS0(:)
  Real(dp), Intent(in), Optional :: mG
  Complex(dp), Intent(in) :: cpl_FNSf_L(:,:,:), cpl_FNSf_R(:,:,:)    &
                              & , cpl_CFpSf_L(:,:,:), cpl_CFpSf_R(:,:,:)
  Complex(dp), Intent(in) :: cpl_SfSfpW(:,:), cpl_SmpSfSfp(:,:,:)    &
             & , cpl_SfSfZ(:,:), cpl_P0SfSf(:,:,:), cpl_S0SfSf(:,:,:)
  Complex(dp), Intent(in), Optional :: cpl_GQSq_L(:,:), cpl_GQSq_R(:,:)
  Real(dp), Intent(inout) :: gP(:,:), gT(:)
  Real(dp), Optional, Intent(inout) :: BR(:,:)
  Logical, Intent(in) :: GenerationMixing

  Integer :: i1, i2, n_sfer, n_neut, n_char, i_gen, i_start, i_end, i_count &
         & , n_sferp, n_Spm, i3, n_P0, n_S0
  Real(dp) :: gam, m_in 
  Complex(dp) :: coupLC, coupRC, coupC
  !-----------------
  ! Initialization
  !-----------------
  Iname = Iname + 1
  NameOfUnit(Iname) = 'SfermionTwoBodyDecays'

  n_sfer = Size(mSf)
  n_sferp = Size(mSfp)
  n_char = Size(mC)
  n_neut = Size(mN)
  n_Spm = Size(mSpm)
  n_P0 = Size(mP0)
  n_S0 = Size(mS0)

  If (i_in.Lt.0) Then
   i_start = 1 
   i_end = n_sfer
   gT = 0._dp
   gP = 0._dp

  Else If ( (i_in.Ge.1).And.(i_in.Le.n_sfer) ) Then 
   i_start = i_in 
   i_end = i_in
   gT(i_in) = 0._dp
   gP(i_in,:) = 0._dp

  Else
   If (ErrorLevel.Ge.-1) Then
    Write(ErrCan,*) 'Problem in subroutine '//NameOfUnit(Iname)
    Write(ErrCan,*) 'Value of i_in out of range, (i_in,n_sfer) = ',i_in,n_sfer
   End If

   If (ErrorLevel.Gt.0) Call TerminateProgram

   If (Present(BR)) BR = 0._dp
   Iname = Iname - 1
   Return
  End If

  Do i1 = i_start, i_end
   m_in = msf(i1)
   i_count = 1
   If (GenerationMixing) Then
    !---------------------------------------------------------------------
    ! into neutralinos, if k_neut=1 or k_neut=3 summing over all fermions
    !---------------------------------------------------------------------
    Do i2 = 1, n_neut
     Do i_gen = 1,3
      coupLC = cpl_FNSf_L(i_gen,i2,i1)
      coupRC = cpl_FNSf_R(i_gen,i2,i1)
      Call ScalarToTwoFermions(m_in, mf(i_gen), mN(i2), coupLC, coupRC, gam)
      If ((k_neut.Eq.1).Or.(k_neut.Eq.3)) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
      Else
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If ((k_neut.Eq.1).Or.(k_neut.Eq.3)) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    !--------------------------------------------------------------------
    ! into charginos, if k_neut=2 or k_neut=3 summing over all fermions
    !--------------------------------------------------------------------
    Do i2 = 1, n_char
     Do i_gen = 1,3
      coupLC = cpl_CFpSf_L(i2, i_gen, i1)
      coupRC = cpl_CFpSf_R(i2, i_gen, i1)

      Call ScalarToTwoFermions(m_in, mfp(i_gen), mC(i2), coupLC, coupRC, gam)
      If ((k_neut.Eq.2).Or.(k_neut.Eq.3)) Then
       gP(i1, i_count) = gP(i1, i_count) + gam 
      Else
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)
       i_count = i_count + 1
      End If
     End Do
     If ((k_neut.Eq.2).Or.(k_neut.Eq.3)) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End Do
    !----------------------------------------------------
    ! into gluino, if k_neut=3 summing over all quarks
    !----------------------------------------------------
    If (Present(mG).And.Present(cpl_GQSq_L).And.Present(cpl_GQSq_R) ) Then
     Do i_gen = 1,3
      coupLC = cpl_GQSq_L(i_gen, i1)
      coupRC = cpl_GQSq_R(i_gen, i1)
      Call ScalarToTwoFermions(m_in, mf(i_gen), mG, coupLC, coupRC, gam)
      If (k_neut.Eq.3) Then
       gP(i1, i_count) = gP(i1, i_count) + 16._dp * gam / 3._dp ! Colour factor
      Else
       gP(i1, i_count) = 16._dp * gam / 3._dp ! Colour factor 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End If
     End Do
     If (k_neut.Eq.3) Then
      gT(i1) = gT(i1) + gP(i1, i_count)
      i_count = i_count + 1
     End If
    End If
    !-----------------
    ! W-boson
    !-----------------
    Do i2=1,n_sferp
     coupC = cpl_SfSfpW(i1, i2)
     Call ScalarToScalarVectorBoson(m_in,mSfp(i2),mW,coupC,gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !-----------------
    ! charged scalar
    !-----------------
    Do i2=2,n_Spm
     Do i3=1,n_sferp
      coupC = cpl_SmpSfSfp(i2, i1, i3)
      Call ScalarToTwoScalars(m_in, mSfp(i3), mSpm(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End Do
    !-----------------
    ! Z-boson
    !-----------------
    Do i2=1,i1-1
     coupC = cpl_SfSfZ(i1, i2)
     Call ScalarToScalarVectorBoson(m_in, mSf(i2), mZ, coupC, gam)
     gP(i1, i_count) = gam 

     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !-----------------
    ! pseudoscalar
    !-----------------
    Do i2=2,n_P0
     Do i3=1,i1-1
      coupC = cpl_P0SfSf(i2, i1, i3)
      Call ScalarToTwoScalars(m_in, mSf(i3), mP0(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End Do
    !-----------------
    ! scalar
    !-----------------
    Do i2=1,n_S0
     Do i3=1,i1-1
      coupC = cpl_S0SfSf(i2, i1, i3)
      Call ScalarToTwoScalars(m_in, mSf(i3), mS0(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End Do

   Else  ! GenerationMixing = .FALSE.
    If (n_sfer.Eq.6) Then
     i_gen = (i1+1) / 2
    Else
     i_gen = i1
    Endif
    !--------------------------
    ! into neutralinos
    !--------------------------
    Do i2 = 1,n_neut
     coupLC = cpl_FNSf_L(i_gen,i2,i1)
     coupRC = cpl_FNSf_R(i_gen,i2,i1)
     Call ScalarToTwoFermions(m_in, mf(i_gen), mN(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !--------------------------
    ! into charginos
    !--------------------------
    Do i2 = 1, n_char
     coupLC = cpl_CFpSf_L(i2, i_gen, i1)
     coupRC = cpl_CFpSf_R(i2, i_gen, i1)
     Call ScalarToTwoFermions(m_in, mfp(i_gen), mC(i2), coupLC, coupRC, gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End Do
    !--------------------------
    ! into gluino
    !--------------------------
    If (Present(mG).And.Present(cpl_GQSq_L).And.Present(cpl_GQSq_R) ) Then
     coupLC = cpl_GQSq_L(i_gen, i1)
     coupRC = cpl_GQSq_R(i_gen, i1)
     Call ScalarToTwoFermions(m_in, mf(i_gen), mG, coupLC, coupRC, gam)
     gP(i1, i_count) = 16._dp * gam / 3._dp ! Colour factor 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    End If

    !-----------------
    ! W-boson
    !-----------------
    If (n_sferp.Eq.3) Then
     coupC = cpl_SfSfpW(i1, i_gen)
     Call ScalarToScalarVectorBoson(m_in,mSfp(i_gen),mW,coupC,gam)
     gP(i1, i_count) = gam 
     gT(i1) = gT(i1) + gP(i1, i_count)       
     i_count = i_count + 1
    Else
     Do i2 =1,2
      coupC = cpl_SfSfpW(i1, (i_gen-1)*2 + i2)
      Call ScalarToScalarVectorBoson(m_in,mSfp((i_gen-1)*2 + i2),mW,coupC,gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    End If   

    !-----------------
    ! charged scalar
    !-----------------
    If (n_sferp.Eq.3) Then
     Do i2 = 2, n_Spm
      coupC = cpl_SmpSfSfp(i2, i1, i_gen)
      Call ScalarToTwoScalars(m_in, mSfp(i_gen), mSpm(i2), coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End Do
    Else
     Do i2 =2,n_Spm
      Do i3 =1,2
       coupC = cpl_SmpSfSfp(i2, i1, (i_gen-1)*2 + i3)
       Call ScalarToTwoScalars(m_in, mSfp((i_gen-1)*2+i3), mSpm(i2), coupC &
                             &, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End Do
     End Do
    End If   

    !-----------------
    ! Z-boson
    !-----------------
    If (n_sfer.Ne.3) Then ! no Sneutrino
     If (i1.Eq.2) Then
      coupC = cpl_SfSfZ(1, 2)
      Call ScalarToScalarVectorBoson(m_in, mSf(1), mZ, coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     Else If (i1.Eq.4) Then
      coupC = cpl_SfSfZ(3, 4)
      Call ScalarToScalarVectorBoson(m_in, mSf(3), mZ, coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     Else If (i1.Eq.6) Then
      coupC = cpl_SfSfZ(5, 6)
      Call ScalarToScalarVectorBoson(m_in, mSf(5), mZ, coupC, gam)
      gP(i1, i_count) = gam 
      gT(i1) = gT(i1) + gP(i1, i_count)       
      i_count = i_count + 1
     End If   
    End If   

    !-----------------
    ! pseudoscalar
    !-----------------
    If (n_sfer.Ne.3) Then ! no Sneutrino
     Do i2=2,n_P0
      If (i1.Eq.2) Then
       coupC = cpl_P0SfSf(i2, 1, 2)
       Call ScalarToTwoScalars(m_in, mSf(1), mP0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.4) Then
       coupC = cpl_P0SfSf(i2, 3, 4)
       Call ScalarToTwoScalars(m_in, mSf(3), mP0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.6) Then
       coupC = cpl_P0SfSf(i2, 5, 6)
       Call ScalarToTwoScalars(m_in, mSf(5), mP0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End If   
     End Do
    End If   

    !-----------------
    ! scalar
    !-----------------
    If (n_sfer.Ne.3) Then ! no Sneutrino
     Do i2=1,n_S0
      If (i1.Eq.2) Then
       coupC = cpl_S0SfSf(i2, 1, 2)
       Call ScalarToTwoScalars(m_in, mSf(1), mS0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.4) Then
       coupC = cpl_S0SfSf(i2, 3, 4)
       Call ScalarToTwoScalars(m_in, mSf(3), mS0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      Else If (i1.Eq.6) Then
       coupC = cpl_S0SfSf(i2, 5, 6)
       Call ScalarToTwoScalars(m_in, mSf(5), mS0(i2), coupC, gam)
       gP(i1, i_count) = gam 
       gT(i1) = gT(i1) + gP(i1, i_count)       
       i_count = i_count + 1
      End If   
     End Do
    End If   

   End If    ! GenerationMixing


   If ((Present(BR)).And.(gT(i1).Eq.0)) Then
    BR(i1,:) = 0._dp
   Else If (Present(BR)) Then
    BR(i1,:) = gP(i1,:) / gT(i1)
   End If

  End Do ! i1
 
  Iname = Iname - 1

 End Subroutine SfermionTwoBodyDecays

End Module  SusyDecays

