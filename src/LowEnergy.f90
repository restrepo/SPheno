Module LowEnergy

! load modules
Use Control
Use Couplings
Use LoopCouplings
Use LoopFunctions
Use LoopMasses
Use StandardModel
! load modules

! interfaces

Interface C_10_10p
 Module Procedure C_QdQdLL_AAp
End Interface

Interface C_11_11p
 Module Procedure C_QdQdNuNu_LR
End Interface

Interface Gminus2
 Module Procedure Gminus2a, Gminus2b
End Interface

Interface LtoLpGamma
 Module Procedure  LtoLpGammaMSSM, LtoLpGammaEPS3
End Interface
! interfaces 

Real(dp) :: L_edm_contri(5,2)
! private variables
Real(dp), Private, Parameter :: sa(7) = (/ -1._dp / 3._dp, 1._dp / 4._dp &
    & ,  -1._dp / 5._dp,  1._dp / 6._dp,  -1._dp / 7._dp,   1._dp / 8._dp     &
    & ,  -1._dp / 9._dp /)                                                    &
    & , sb(7) = (/  1._dp / 6._dp, -1._dp /12._dp,  1._dp /20._dp             &
    &            , -1._dp /30._dp,  1._dp /42._dp, -1._dp /56._dp             &
    &            ,  1._dp /72._dp  /)                                         &
    & , sd(7) = (/ -1._dp / 20._dp,  1._dp / 15._dp, -1._dp / 14._dp          &
    &            ,  1._dp / 14._dp, -5._dp / 72._dp,  1._dp / 15._dp          &
    &            , -7._dp /110._dp /)                                         &
    & , se(7) = (/  1._dp / 30._dp, -1._dp / 30._dp,  1._dp / 35._dp          &
    &            , -1._dp / 42._dp,  5._dp /252._dp, -1._dp / 60._dp          &
    &            ,  7._dp /495._dp /)
!-----------------------------------------------------------
! constants for edms, can be changed with routine EDMinit 
!-----------------------------------------------------------
Real(dp), Private :: ecmfactor = 1.973e-14_dp, elimit = 4.3e-27_dp            &
    & , nlimit = 1.1e-25_dp, qcdCorrection = 1.53_dp, CqcdCorrection = 3.4_dp &
    & , chiralMass = 1.19_dp, deltaUp = 0.746_dp, deltaDown = -0.508_dp       &
    & , deltaStrange = -0.226_dp
!-----------------------------------------------------------------
! constants for B-physics, can be changed with routine BmesonInit 
!-----------------------------------------------------------------
Real(dp), Private :: MBq(2) = (/ 5.2794_dp, 5.3696_dp /), etaB = 0.55_dp  &
    & , FBhatB(2) = (/ 0.216_dp, 1.21_dp*0.216_dp /)                      &
!    & , FBhatB(2) = (/ 0.244_dp, 1.21_dp*0.244_dp /)                      &
    & , TauB(2) = (/ 1.537_dp, 1.461_dp/) ! in pico seconds
Real(dp), Private :: mBm(2) = (/ 5.279_dp, 6.286_dp /)  &
    & , TauBm(2) = (/ 1.638e-12_dp, 0.46e-12_dp /) ! in seconds

Real(dp), Parameter, Private :: e_u = 2._dp/3._dp, e_d=-1._dp/3._dp
Logical, Private :: WriteLLpGamma=.False.
!Private :: C_7, C_7p, C_8, C_8p, F1, F2, F3, F4, F3Gamma
Private :: F1, F2, F3, F4, F3Gamma
! private variables

Contains 


 Subroutine BR_lj_to_3li(j, i, gp, g, Y_l, uL_L, uL_R, mSlep2, RSl, mN, N &
                      & , mSnu2, RSn, mC, U, V, mS02, RS0, mP02, RP0, A_l &
                      & , mu, vevSM, Br)
 !-------------------------------------------------------------------------
 ! calculates the lepton number violating decay of lepton into 3 lighter
 ! leptons
 ! Formulas by E. Arganda and M.J. Herrero, Phys. Rev. D73, 055003 (2006).
 ! Based on routines by E. Arganda and M.J. Herrero
 ! 20.07.2009: adding lepton mixing matrices, thanks to António Figueiredo
 !-------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: j, i ! index of decaying lepton and the index of
                              ! daughters
  Real(dp), Intent(in) :: gp, g             ! U(1) and SU(2) gauge couplings
  Complex(dp), Intent(in) :: Y_L(3,3)       ! lepton Yukawa couplings
  Complex(dp), Intent(in) :: uL_L(3,3), uL_R(3,3) ! the corresponding mixing matrices
  Real(dp), Intent(in) :: mSlep2(6)         ! slepton masses squared
  Complex(dp), Intent(in) :: RSl(6,6)       ! slepton mixing matrix
  Real(dp), Intent(in) :: mN(4)             ! neutralino masses
  Complex(dp), Intent(in) :: N(4,4)         ! neutralino mixing matrix 
  Real(dp), Intent(in) :: mSnu2(3)          ! sneutrino masses squared
  Complex(dp), Intent(in) :: RSn(3,3)       ! sneutrino mixing matrix
  Real(dp), Intent(in) :: mC(2)             ! chargino masses
  Complex(dp), Intent(in) :: U(2,2), V(2,2) ! chargino mixing matrices
  Real(dp), Intent(in) :: mS02(2)           ! neutral Higgs boson masses squared
  Real(dp), Intent(in) :: RS0(2,2)          ! neutral Higgs mixing matrix
  Real(dp), Intent(in) :: mP02(2)           ! pseudoscalar Higgs masses squared
  Real(dp), Intent(in) :: RP0(2,2)          ! pseudoscalar Higgs mixing matrix
  Complex(dp), Intent(in) :: A_l(3,3)     ! trilinear Higgs - slepton couplings
  Complex(dp), Intent(in) :: mu             ! mu parameter of the superpotential
  Real(dp), Intent(in) ::  vevSM(2)         ! Higgs vevs (v_d, v_u)
  Real(dp), Intent(out) :: BR
  
  Integer :: i1, i2, i3, i4
  Real(dp) :: mN2(4), mC2(2), fun1, fun2, X_ax, cosW, cosW2, sinW2 &
      & , c_Zee_L, c_Zee_R, fun3, fun4
  Complex(dp) :: c_LNSl_L(3,4,6), c_LNSl_R(3,4,6), QSlepton(6,6), QSnu(3,3) &
            &  , c_CLSn_L(2,3,3), c_CLSn_R(2,3,3), ELn(4,4), ERn(4,4)       &
            &  , ELc(2,2), ERc(2,2), DR(3,4,4), DL(3,4,4)                   &
            &  , WR(3,2,2), WL(3,2,2), GSlepton(3,6,6), bi(1)               &
            &  , c_HLL_L(3,3), c_HLL_R(3,3), GSnu(3,3,3)
  Complex(dp) :: A1L, A1R, A2L, A2R, FL, FR, FLL, FLR, FRL, FRR &
      & , B1L, B1R, B2L, B2R, B3L, B3R, B4L, B4R                &
      & , H2L, H2R, H3L, H3R, wertL, wertR, fun1c, fun2c

  Real(dp), Parameter :: oo576Pi2 = 1._dp / (576._dp * Pi2)
  Complex(dp), Parameter :: mat0(3,3) = 0._dp

  Iname = Iname + 1
  NameOfUnit(Iname) = "BR_lj_to_3li"

  BR = 0._dp

  !------------
  ! couplings
  !------------
  c_LNSl_L = 0._dp
  c_LNSl_R = 0._dp
  Do i1=1,6
   Do i2=1,4
    Call CoupNeutralinoSlepton3(i, i2, i1, gp, g, RSl, uL_L, uL_R, Y_l, N &
                                  &, c_LNSl_L(i,i2,i1), c_LNSl_R(i,i2,i1) )
    Call CoupNeutralinoSlepton3(j, i2, i1, gp, g, RSl, uL_L, uL_R, Y_l, N &
                                  &, c_LNSl_L(j,i2,i1), c_LNSl_R(j,i2,i1) )
   End Do
  End Do

  c_CLSn_L = 0._dp
  c_CLSn_R = 0._dp
  Do i1=1,3
   Do i2=1,2
    Call CoupCharginoSfermion3(i2, j, i1, g, -0.5_dp, RSn, Y_l, mat0, uL_L &
                  & , uL_R, U, V, c_CLSn_L(i2, j, i1), c_CLSn_R(i2, j, i1) )
    Call CoupCharginoSfermion3(i2, i, i1, g, -0.5_dp, RSn, Y_l, mat0, uL_L &
                  & , uL_R, U, V, c_CLSn_L(i2, i, i1), c_CLSn_R(i2, i, i1) )
   End Do
  End Do

  cosW2 = mW2 / mZ2 !g**2 / (gp**2 + g**2)
  cosW = Sqrt(cosW2)
  Do i1=1,4
   Do i2=1,4
    Call CoupNeutralinoZ(i1, i2, N, g, cosW, ELn(i1,i2), ERn(i1,i2) )
   End Do
  End Do

  Do i1=1,2
   Do i2=1,2
    Call CoupCharginoZ(i1, i2, U, V, g, cosW, ELc(i1,i2), ERc(i1,i2) )
   End Do
  End Do

  sinW2 = 1._dp - cosW2
  Do i1=1,6
   Do i2=1,6
    Call CoupSleptonZ(i1, i2, g, sinW2, RSl, QSlepton(i1,i2) )
   End Do
  End Do
  Do i1=1,3
   Do i2=1,3
    Call CoupSneutrinoZ3(i1, i2, g, sinW2, RSn, QSnu(i1,i2) )
   End Do
  End Do
  Qslepton = - Qslepton ! sign difference in paper by Arganda/Herrero
  Qsnu = - Qsnu ! sign difference in paper by  Arganda/Herrero

  Call CoupFermionZ(-0.5_dp, -1._dp, g, sinW2, c_Zee_L,c_Zee_R)

  DR = 0._dp
  DL = 0._dp
  Do i2=1,4
   Do i3=1,4
    Do i1=1,2
     Call CoupNeutralinoScalar(i2,i3,i1, N, RS0, gp, g &
                             &, DL(i1,i2,i3), DR(i1,i2,i3) )
    End Do
     Call CoupNeutralinoPseudoScalar(i2,i3,2, N, RP0, gp, g &
                             &, DL(3,i2,i3), DR(3,i2,i3) )
   End Do
  End Do

  WR = 0._dp
  WL = 0._dp
  Do i1=1,2
   Do i2=1,2
    Do i3=1,2
     Call CoupCharginoScalar(i1,i2,i3, U, V, RS0, g, WL(i3,i2,i1), WR(i3,i2,i1))
    End Do
    Call CoupCharginoPseudoScalar(i1,i2,2, U, V, RP0, g &
                                & , WL(3,i2,i1), WR(3,i2,i1))
   End Do
  End Do

  GSlepton = 0._dp
  bi = mu
  Do i1=1,6
   Do i2=1,6
    Do i3=1,2
     Call CoupScalarSfermion3MSSM_3(i3, i1, i2, RS0, -0.5_dp, -1._dp, Y_L &
                 & , Rsl, A_l, mu, vevSM, gp, g, GSlepton(i3, i1, i2) )
    End Do
    Call CoupPseudoScalarSfermion3_3(2, i1, i2, RP0, -0.5_dp, Y_l, Rsl &
                 & , A_l, bi, GSlepton(3, i1, i2) )
   End Do
  End Do

  GSnu = 0._dp
  bi = mu
  Do i1=1,3
   Do i2=1,3
    Do i3=1,2
     Call CoupScalarSfermion3MSSM_3(i3, i1, i2, RS0, 0.5_dp, 0._dp, mat0 &
                 & , Rsn, mat0, mu, vevSM, gp, g, GSnu(i3, i1, i2) )
    End Do
   End Do
  End Do

  c_HLL_L = 0._dp
  c_HLL_R = 0._dp
  Do i1=1,3
   Do i2=1,2
    Call CoupFermionScalar3(i1, i1, i2, -0.5_dp, Y_l, uL_L, uL_R, RS0 &
                          & , c_HLL_L(i2,i1), c_HLL_R(i2,i1) )
   End Do
   Call CoupFermionPseudoScalar3(i1, i1, 2, -0.5_dp, Y_l, uL_L, uL_R, RP0 &
                          & , c_HLL_L(3,i1), c_HLL_R(3,i1) )
  End Do

  !---------------------------------------------
  ! Penguin contributions
  !---------------------------------------------
  ! neutralinos
  !--------------
  A1L = 0._dp
  A1R = 0._dp
  A2L = 0._dp
  A2R = 0._dp
  mN2 = mN**2

  Do i1=1,6
   Do i2=1,4
    X_ax = mN2(i2)/mSlep2(i1)
    fun1 = ( 2._dp - 9._dp*X_ax + 18._dp*X_ax**2 - 11._dp*X_ax**3 & 
        & + 6._dp*X_ax**3 * Log(X_ax) ) /  (1._dp-X_ax)**4
    A1L = A1L &
      & + c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1)) * fun1 / mSlep2(i1)
    A1R = A1R &
      & + c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1)) * fun1 / mSlep2(i1)

    fun1 = 2._dp * F2(X_ax)
    fun2 = 2._dp * F4(X_ax)
    A2L = A2L + ( ( c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1))      & 
        &         + c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1))      &
        &           *mf_l(i)/mf_l(j) ) * fun1                         &
        &       + c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1))        &
        &         * mN(i2)/mf_l(j)*fun2 ) / mSlep2(i1)
    A2R = A2R + ( ( c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1))      & 
        &         + c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1))      &
        &           *mf_l(i)/mf_l(j) ) * fun1                         &
        &       + c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1))        &
        &         * mN(i2)/mf_l(j)*fun2 ) / mSlep2(i1)
   End Do
  End Do

  !--------------
  ! charginos
  !--------------
  mC2 = mC**2
  Do i1=1,3
   Do i2=1,2
    X_ax=mC2(i2)/mSnu2(i1)
    fun1 = ( 16._dp - 45._dp*X_ax + 36._dp * X_ax**2 - 7._dp*X_ax**3 &
         & + 6._dp * (2._dp-3._dp*X_ax) * Log(X_ax) ) / (1._dp-X_ax)**4
    A1L = A1L - c_CLSn_R(i2,i,i1) * Conjg(c_CLSn_R(i2,j,i1)) * fun1 / mSnu2(i1)
    A1R = A1R - c_CLSn_L(i2,i,i1) * Conjg(c_CLSn_L(i2,j,i1)) * fun1 / mSnu2(i1)
    fun1 = 2._dp * F1(X_ax)
    fun2 = 2._dp * F3(X_ax)
    A2L = A2L - ( ( c_CLSn_L(i2,i,i1)*Conjg(c_CLSn_L(i2,j,i1))                 & 
        &         + c_CLSn_R(i2,i,i1)*Conjg(c_CLSn_R(i2,j,i1))     &
        &           *mf_l(i)/mf_l(j) ) *fun1                       &
        &       + c_CLSn_L(i2,i,i1) * Conjg(c_CLSn_R(i2,j,i1))     &
        &          * mC(i2) / mf_l(j) * fun2 ) / mSnu2(i1)
    A2R = A2R - ( ( c_CLSn_R(i2,i,i1)*Conjg(c_CLSn_R(i2,j,i1))                 & 
        &         + c_CLSn_L(i2,i,i1)*Conjg(c_CLSn_L(i2,j,i1))     &
        &           *mf_l(i)/mf_l(j) ) *fun1                       &
        &       + c_CLSn_R(i2,i,i1) * Conjg(c_CLSn_L(i2,j,i1))     &
        &          * mC(i2) / mf_l(j) * fun2 ) / mSnu2(i1)
   End Do
  End Do

  A1L = oo576Pi2 * A1L 
  A1R = oo576Pi2 * A1R 
  A2L = oo32Pi2 * A2L 
  A2R = oo32Pi2 * A2R 

  !-----------
  ! Z-penguin
  !-----------
  FL = 0._dp
  FR = 0._DP

  Do i1=1,6
   Do i2=1,4
    Do i3=1,4
     fun1 = 0.5_dp * vertexC0tilde(mSlep2(i1),mN2(i2),mN2(i3))
     fun2 = mN(i2) * mN(i3) * C0_3m(mSlep2(i1),mN2(i2),mN2(i3))
     FL = FL - c_LNSl_R(i,i3,i1)*Conjg(c_LNSl_R(j,i2,i1))             & 
      &   * ( ERn(i3,i2) * fun1 - ELn(i3,i2) * fun2 )
     FR = FR - c_LNSl_L(i,i3,i1)*Conjg(c_LNSl_L(j,i2,i1))             & 
      &   * ( ELn(i3,i2) * fun1 - ERn(i3,i2) * fun2 )
    End Do
   End Do
  End Do

   Do i1=1,4
    Do i2=1,6
     Do i3=1,6
      fun1 = 0.5_dp * Qslepton(i2,i3)  &
           &        * vertexC0tilde(mN2(i1),mSlep2(i2),mSlep2(i3)) 
      FL = FL - c_LNSl_R(i,i1,i2) * Conjg(c_LNSl_R(j,i1,i3)) * fun1
      FR = FR - c_LNSl_L(i,i1,i2) * Conjg(c_LNSl_L(j,i1,i3)) * fun1
     End Do
    End Do
   End Do

   Do i1=1,6
    Do i2=1,4
     fun1 = B1(0._dp, mN2(i2), mSlep2(i1) )
     FL = FL - c_LNSl_R(i,i2,i1) * Conjg(c_LNSl_R(j,i2,i1)) * c_Zee_L * fun1
     FR = FR - c_LNSl_L(i,i2,i1) * Conjg(c_LNSl_L(j,i2,i1)) * c_Zee_R * fun1
    End Do
   End Do

   Do i1=1,3
    Do i2=1,2
     Do i3=1,2
      fun1 = 0.5_dp * vertexC0tilde(mSnu2(i1), mC2(i2), mC2(i3) )
      fun2 = C0_3m(mSnu2(i1), mC2(i2), mC2(i3) )
      FL = FL - c_CLSn_R(i3,i,i1)*Conjg(c_CLSn_R(i2,j,i1) )             & 
         &  * ( ERc(i3,i2) * fun1 - ELc(i3,i2) * mC(i2) * mC(i3) * fun2 )
      FR = FR - c_CLSn_L(i3,i,i1)*Conjg(c_CLSn_L(i2,j,i1) )             & 
         &  * ( ELc(i3,i2) * fun1 - ERc(i3,i2) * mC(i2) * mC(i3) * fun2 )
     End Do
    End Do
   End Do

   Do i1=1,3
    Do i2=1,3
     Do i3=1,2
      fun1 = 0.5_dp * vertexC0tilde(mC2(i3), mSnu2(i1), mSnu2(i2) )
      FL = FL - c_CLSn_R(i3,i,i1) * Conjg(c_CLSn_R(i3,j,i2)) &
         &       * Qsnu(i1,i2) * fun1
      FR = FR - c_CLSn_L(i3,i,i1) * Conjg(c_CLSn_L(i3,j,i2)) &
         &       * Qsnu(i1,i2) * fun1
     End Do
    End Do
   End Do

   Do i1=1,3
    Do i2=1,2
     fun1 = B1(0._dp, mC2(i2), mSnu2(i1) )
     FL = FL - c_CLSn_R(i2,i,i1)*Conjg(c_CLSn_R(i2,j,i1)) * c_Zee_L * fun1
     FR = FR - c_CLSn_L(i2,i,i1)*Conjg(c_CLSn_L(i2,j,i1)) * c_Zee_R * fun1
    End Do
   End Do

  FL = oo16pi2 * FL
  FR = oo16pi2 * FR

  fun1 = 1._dp / (g**2 * sinW2 * mZ2)
  FLL = - fun1 * c_Zee_L * FL 
  FLR = - fun1 * c_Zee_R * FL 
  FRL = - fun1 * c_Zee_L * FR 
  FRR = - fun1 * c_Zee_R * FR

  !-------------------
  ! box contributions
  !-------------------
  B1L = 0._dp
  B1R = 0._dp
  B2L = 0._dp
  B2R = 0._dp
  B3L = 0._dp
  B3R = 0._dp
  B4L = 0._dp
  B4R = 0._dp


  Do i1=1,6
   Do i2=1,6
    Do i3=1,4
     Do i4=1,4
      fun1 = D27_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSlep2(i2) )
      fun2 = mN(i3) * mN(i4)  &
         &   * D0_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSlep2(i2) )

      B1L = B1L + 2._dp * fun1 * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1)) &
          &                    * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2)) & 
          &      + fun2 * c_LNSl_R(i,i4,i2) * c_LNSl_R(i,i4,i1)               &
          &             * Conjg(c_LNSl_R(j,i3,i1)) * Conjg(c_LNSl_R(i,i3,i2))

      B1R = B1R + 2._dp * fun1 * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1)) &
          &                    * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2)) & 
          &      + fun2 * c_LNSl_L(i,i4,i2) * c_LNSl_L(i,i4,i1)               &
          &             * Conjg(c_LNSl_L(j,i3,i1)) * Conjg(c_LNSl_L(i,i3,i2))

      B2L = B2L + fun1 * ( c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1))        &
          &                * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2))      &
          &              - c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_R(j,i3,i1))        &
          &                * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_L(i,i3,i2))      &
          &              + c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_L(i,i3,i2))        &
          &                *Conjg(c_LNSl_R(j,i3,i1)) * c_LNSl_L(i,i4,i1) )     &
          &     - 0.5_dp * fun2 * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1)) &
          &                     * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2))

      B2R = B2R + fun1 * ( c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1))        &
          &                * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2))      &
          &              - c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_L(j,i3,i1))        &
          &                * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_R(i,i3,i2))      &
          &              + c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_R(i,i3,i2))        &
          &                *Conjg(c_LNSl_L(j,i3,i1)) * c_LNSl_R(i,i4,i1) )     &
          &     - 0.5_dp * fun2 * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1)) &
          &                     * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2))

      B3L = B3L + fun2 * ( c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i1))        &
          &              * c_LNSl_L(i,i4,i1) * Conjg(c_LNSl_R(i,i4,i2))        &
          &              + 0.5_dp * c_LNSl_L(i,i4,i2) * c_LNSl_L(i,i4,i1)      &
          &               * Conjg(c_LNSl_R(j,i3,i1))* Conjg(c_LNSl_R(i,i3,i2)) )

      B3R = B3R + fun2 * ( c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i1))        &
          &              * c_LNSl_R(i,i4,i1) * Conjg(c_LNSl_L(i,i4,i2))        &
          &              + 0.5_dp * c_LNSl_R(i,i4,i2) * c_LNSl_R(i,i4,i1)      &
          &               * Conjg(c_LNSl_L(j,i3,i1))* Conjg(c_LNSl_L(i,i3,i2)) )

      B4L = B4L + fun2 * Conjg(c_LNSl_R(j,i3,i1)) * Conjg(c_LNSl_R(i,i3,i2))  &
         &            * c_LNSl_L(i,i4,i2) * c_LNSl_L(i,i4,i1) / 8._dp

      B4R = B4R + fun2 * Conjg(c_LNSl_L(j,i3,i1)) * Conjg(c_LNSl_L(i,i3,i2))  &
          &            * c_LNSl_R(i,i4,i2) * c_LNSl_R(i,i4,i1) / 8._dp

     End Do
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i2=1,3
    Do i3=1,2
     Do i4=1,2
      fun1 = D27_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSnu2(i2) )
      fun2 = mC(i3) * mC(i4) * D0_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSnu2(i2))

      B1L = B1L + 2._dp * fun1 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1)) &
          &                    * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))

      B2L = B2L + fun1 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1))          &
          &            * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))          & 
          &     - 0.5_dp * fun2 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1)) &
          &                     * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))

      B3L = B3L + fun2 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i1))    &
          &            * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))

      B1R = B1R + 2._dp * fun1 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1)) &
          &                    * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))

      B2R = B2R + fun1 * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1))          &
          &            * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))          & 
          &     - 0.5_dp * fun2 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1)) &
          &                     * c_CLSn_L(i4,i,i1) * Conjg(c_CLSn_R(i4,i,i2))

      B3R = B3R + fun2 * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i1))    &
          &            * c_CLSn_R(i4,i,i1) * Conjg(c_CLSn_L(i4,i,i2))

     End Do
    End Do
   End Do
  End Do

  fun1 = oo16pi2 / (g**2*sinW2)
  B1L = fun1 * B1L
  B1R = fun1 * B1R
  B2L = fun1 * B2L
  B2R = fun1 * B2R
  B3L = fun1 * B3L
  B3R = fun1 * B3R
  B4L = fun1 * B4L
  B4R = fun1 * B4R
  !---------------------------
  ! Higgs contributions
  !---------------------------
  H2L = 0._dp
  H2R = 0._dp
  H3L = 0._dp
  H3R = 0._dp

  Do i1=1,3
   wertL = 0._dp
   wertR = 0._dp
   Do i2=1,4
    Do i3=1,4
     fun1 = B0(0._dp, mN2(i2), mN2(i3))
     Do i4=1,6
      fun2 = C0_3m(mSlep2(i4), mN2(i2), mN2(i3) )
      fun3 = vertexC12(mSlep2(i4), mN2(i2), mN2(i3) )
      fun4 = vertexC11(mSlep2(i4), mN2(i2), mN2(i3) )

      wertL = wertL &
        &   + ( fun1 + mSlep2(i4) * fun2 +mf_l2(j) * fun3                    &
        &       + mf_l2(i)*( fun4 - fun3 ) )                                 &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(i) * mf_l(j) *( fun4 + fun2 )                             &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(i) * mN(i3) *( fun4 - fun3 + fun2 )                       &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(j) * mN(i3) * fun3                                        &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(i) * mN(i2) * ( fun4 - fun3 )                             &
        &      * c_LNSl_R(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(j)* mN(i2) * ( fun3 + fun2 )                              &
        &      * c_LNSl_L(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mN(i2) * mN(i3) * fun2                                         &
        &      * c_LNSl_L(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4))

      wertR = wertR &
        &   + ( fun1 + mSlep2(i4) * fun2 +mf_l2(j) * fun3                    &
        &       + mf_l2(i)*( fun4 - fun3 ) )                                 &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(i) * mf_l(j) *( fun4 + fun2 )                             &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(i) * mN(i3) *( fun4 - fun3 + fun2 )                       &
        &      * c_LNSl_L(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(j) * mN(i3) * fun3                                        &
        &      * c_LNSl_R(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mf_l(i) * mN(i2) * ( fun4 - fun3 )                             &
        &      * c_LNSl_L(i,i2,i4) * DL(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4)) &
        &   + mf_l(j)* mN(i2) * ( fun3 + fun2 )                              &
        &      * c_LNSl_R(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_R(j,i3,i4)) &
        &   + mN(i2) * mN(i3) * fun2                                         &
        &      * c_LNSl_R(i,i2,i4) * DR(i1,i2,i3) * Conjg(c_LNSl_L(j,i3,i4))
     End Do
    End Do
   End Do

   Do i2=1,6
    Do i3=1,6
     Do i4=1,4
      fun1 = vertexC11(mN2(i4), mSlep2(i2), mSlep2(i3))
      fun2 = vertexC12(mN2(i4), mSlep2(i2), mSlep2(i3))
      fun3 = C0_3m(mN2(i4), mSlep2(i2), mSlep2(i3))

      wertL = wertL                                                            &
          &  - ( mf_l(i) * (fun1-fun2)                                         &
          &              *c_LNSl_R(i,i4,i2) *Conjg(c_LNSl_R(j,i4,i3))          &
          &    + mf_l(j) * fun2 * c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_L(j,i4,i3)) &
          &    - mN(i4) * fun3 * c_LNSl_L(i,i4,i2) * Conjg(c_LNSl_R(j,i4,i3))  &
          &    ) * GSlepton(i1,i2,i3)
      wertR = wertR                                                            &
          &  - ( mf_l(i) * (fun1-fun2)                                         &
          &              *c_LNSl_L(i,i4,i2) *Conjg(c_LNSl_L(j,i4,i3))          &
          &    + mf_l(j) * fun2 * c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_R(j,i4,i3)) &
          &    - mN(i4) * fun3 * c_LNSl_R(i,i4,i2) * Conjg(c_LNSl_L(j,i4,i3))  &
          &    ) * GSlepton(i1,i2,i3)
     End Do
    End Do
   End Do

   Do i2=1,6
    Do i3=1,4
     fun1 = B1(0._dp, mN2(i3), mSlep2(i2))
     fun2 = B0(0._dp, mN2(i3), mSlep2(i2))

     wertL = wertL                                                     &
         & + (-mf_l2(i) * fun1                                         &
         &              * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2)) &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  & 
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             *c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))   & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             *c_LNSl_L(i,i3,i2)*Conjg(c_LNSl_R(j,i3,i2))     &
         &   ) * c_HLL_L(i1,j) / (mf_l2(i) - mf_l2(j) )                &
         & + (-mf_l2(j) * fun1                                         &
         &              * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2)) & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  &
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  &
         &   ) * c_HLL_L(i1,i) / (mf_l2(j)-mf_l2(i))
     wertR = wertR                                                     &
         & + (-mf_l2(i) * fun1                                         &
         &              * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2)) &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  & 
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             *c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))   & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             *c_LNSl_R(i,i3,i2)*Conjg(c_LNSl_L(j,i3,i2))     &
         &   ) * c_HLL_R(i1,j) / (mf_l2(i) - mf_l2(j) )                &
         & + (-mf_l2(j) * fun1                                         &
         &              * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2)) & 
         &   + mf_l(j) * mN(i3) * fun2                                 &
         &             * c_LNSl_L(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  &
         &   - mf_l(i) * mf_l(j) * fun1                                &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_R(j,i3,i2))  &
         &   + mf_l(i) * mN(i3) * fun2                                 &
         &             * c_LNSl_R(i,i3,i2) * Conjg(c_LNSl_L(j,i3,i2))  &
         &   ) * c_HLL_R(i1,i) / (mf_l2(j)-mf_l2(i))
      End Do
   End Do

   Do i2=1,2
    Do i3=1,2
     fun1 = B0(0._dp, mC2(i2), mC2(i3))
     Do i4=1,3
      fun2 = C0_3m(mSnu2(i4), mC2(i2), mC2(i3))
      fun3 = vertexC12(mSnu2(i4), mC2(i2), mC2(i3))
      fun4 = vertexC11(mSnu2(i4), mC2(i2), mC2(i3))

      wertL = wertL                                                         &
        & + ( fun1 + mSnu2(i4) * fun2 + mf_l2(j) * fun3                     &
        &   + mf_l2(i) * (fun4 -fun3) )                                     &
        &     * c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(i) * mf_l(j) * (fun4  +fun2)                             &
        &     *c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4))  &
        &   + mf_l(i) * mC(i3) * (fun4 - fun3 + fun2)                       &
        &     * c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(j) * mC(i3) * fun3                                       &
        &     * c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(i) * mC(i2) *(fun4 -fun3)                                &
        &     * c_CLSn_R(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(j) * mC(i2) * (fun3 + fun2)                              &
        &     * c_CLSn_L(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mC(i2) * mC(i3) * fun2                                        &
        &     * c_CLSn_L(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4))

      wertR = wertR                                                         &
        & + ( fun1 + mSnu2(i4) * fun2 + mf_l2(j) * fun3                     &
        &   + mf_l2(i) * (fun4 -fun3) )                                     &
        &     * c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(i) * mf_l(j) * (fun4  +fun2)                             &
        &     *c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4))  &
        &   + mf_l(i) * mC(i3) * (fun4 - fun3 + fun2)                       &
        &     * c_CLSn_L(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(j) * mC(i3) * fun3                                       &
        &     * c_CLSn_R(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mf_l(i) * mC(i2) *(fun4 -fun3)                                &
        &     * c_CLSn_L(i2,i,i4) * WL(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4)) &
        &   + mf_l(j) * mC(i2) * (fun3 + fun2)                              &
        &     * c_CLSn_R(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_R(i3,j,i4)) &
        &   + mC(i2) * mC(i3) * fun2                                        &
        &     * c_CLSn_R(i2,i,i4) * WR(i1,i2,i3) * Conjg(c_CLSn_L(i3,j,i4))
     End Do
    End Do
   End Do

   If (i1.Lt.3) Then ! pseudoscalar does not contribute
    Do i2=1,3
     Do i3=1,3
      Do i4=1,2
       fun1 = vertexC11(mC2(i4), mSnu2(i2), mSnu2(i3))
       fun2 = vertexC12(mC2(i4), mSnu2(i2), mSnu2(i3))
       fun3 = C0_3m(mC2(i4), mSnu2(i2), mSnu2(i3))

       wertL = wertL                                                        &
        & - ( mf_l(i) * (fun1 -fun2)                                        &
        &             * c_CLSn_R(i4,i,i2) *Conjg(c_CLSn_R(i4,j,i3))         & 
        &   + mf_l(j) * fun2 * c_CLSn_L(i4,i,i2) * Conjg(c_CLSn_L(i4,j,i3)) & 
        &   - mC(i4) * fun3 * c_CLSn_L(i4,i,i2) * Conjg(c_CLSn_R(i4,j,i3))  &
        &   ) * GSnu(i1,i2,i3)

       wertR = wertR                                                        &
        & - ( mf_l(i) * (fun1 -fun2)                                        &
        &             * c_CLSn_L(i4,i,i2) *Conjg(c_CLSn_L(i4,j,i3))         & 
        &   + mf_l(j) * fun2 * c_CLSn_R(i4,i,i2) * Conjg(c_CLSn_R(i4,j,i3)) & 
        &   - mC(i4) * fun3 * c_CLSn_R(i4,i,i2) * Conjg(c_CLSn_L(i4,j,i3))  &
        &   ) * GSnu(i1,i2,i3)
      End Do
     End Do
    End Do 
   End If

   Do i2=1,3
    Do i3=1,2
     fun1 = B1(0._dp, mC2(i3), mSnu2(i2))
     fun2 = B0(0._dp, mC2(i3), mSnu2(i2))

     wertL = wertL                                                   &
     & + ( - mf_l2(i) * fun1                                         &
     &       * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))          & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   + mf_l(j) * mC(i3) * fun2                                   &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   )*c_HLL_L(i1,j) / (mf_l2(i)-mf_l2(j))                       &
     & + ( - mf_l2(j) * fun1                                         &
     &                * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2)) & 
     &   + mf_l(j) * mC(i3) * fun2                                   & 
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  & 
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   ) * c_HLL_L(i1,i) / (mf_l2(j)-mf_l2(i))

     wertR = wertR                                                   &
     & + ( - mf_l2(i) * fun1                                         &
     &       * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))          & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  &
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   + mf_l(j) * mC(i3) * fun2                                   &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   )*c_HLL_R(i1,j) / (mf_l2(i)-mf_l2(j))                       &
     & + ( - mf_l2(j) * fun1                                         &
     &                * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2)) & 
     &   + mf_l(j) * mC(i3) * fun2                                   & 
     &             * c_CLSn_L(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   - mf_l(i) * mf_l(j) * fun1                                  & 
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_R(i3,j,i2))    & 
     &   + mf_l(i) * mC(i3) * fun2                                   &
     &             * c_CLSn_R(i3,i,i2) * Conjg(c_CLSn_L(i3,j,i2))    & 
     &   ) * c_HLL_R(i1,i) / (mf_l2(j)-mf_l2(i))
    End Do
   End Do

   If (i1.Lt.3) Then
    fun1c = c_HLL_R(i1,i) / mS02(i1)
    fun2c = c_HLL_L(i1,i) / mS02(i1)
   Else
    fun1c = c_HLL_R(i1,i) / mP02(2)
    fun2c = c_HLL_L(i1,i) / mP02(2)
   End If
   H2L = H2L + 0.5_dp * fun1c * wertL 
   H2R = H2R + 0.5_dp * fun2c * wertR 
   H3L = H3L - fun2c * wertL 
   H3R = H3R - fun1c * wertR 
  End Do
  !------------------------------------------
  ! adding Bbox and Higgs contributions
  ! as the contribute to the same operator
  !------------------------------------------
  fun1 = oo16pi2 / (g**2*sinW2)

  B2L = B2L + fun1 * H2L
  B2R = B2R + fun1 * H2R
  B3L = B3L + fun1 * H3L
  B3R = B3R + fun1 * H3R


  !-------------------------
  ! now the branching ratio
  !-------------------------
  BR = Abs(A1L)**2 + Abs(A1R)**2                                         &
   & - 4._dp* Real(A1L*Conjg(A2R)+A2L*Conjg(A1R),dp)                     &
   & + (Abs(A2L)**2 + Abs(A2R)**2)                                       &
   &     * (16._dp * Log(mf_l(j)/mf_l(i)) - 22._dp ) /3._dp              &
   & + ( Abs(B1L)**2 + Abs(B1R)**2 ) / 6._dp                             &
   & + ( Abs(B2L)**2 + Abs(B2R)**2 ) / 3._dp                             &
   & + ( Abs(B3L)**2 + Abs(B3R)**2 ) / 24._dp                            &
   & + 6._dp * ( Abs(B4L)**2 + Abs(B4R)**2 )                             &
   & - Real( B3L * Conjg(B4L) + B3R * Conjg(B4R),dp)                     &
   & + 2._dp * Real( A1L * Conjg(B1L) + A1R * Conjg(B1R)                 &
   &               + A1L * Conjg(B2L) + A1R * Conjg(B2R), dp) / 3._dp    &
   & - 4._dp * Real( A2R * Conjg(B1L) + A2L * Conjg(B1R)                 &
   &               + A2L * Conjg(B2R) + A2R * Conjg(B2L), dp) / 3._dp    &
   & + ( 2._dp * (Abs(FLL)**2 + Abs(FRR)**2) + Abs(FLR)**2 + Abs(FRL)**2 &
   &   + 2._dp * Real( B1L * Conjg(FLL) + B1R * Conjg(FRR)               &
   &                 + B2L * Conjg(FLR) + B2R * Conjg(FRL), dp )         &
   &   + 4._dp * Real( A1L * Conjg(FLL) + A1R * Conjg(FRR), dp )         &
   &   + 2._dp * Real( A1L * Conjg(FLR) + A1R * Conjg(FRL), dp )         &
   &   - 8._dp * Real( A2R * Conjg(FLL) + A2L * Conjg(FRR), dp )         &
   &   - 4._dp * Real( A2L * Conjg(FRL) + A2R * Conjg(FLR), dp ) ) / 3._dp

  !----------------------------------------------------------------------
  ! taking alpha(Q=0) instead of alpha(m_Z) as this contains most of the
  ! running of the Wilson coefficients
  !----------------------------------------------------------------------
  BR = oo32pi * Alpha**2 * mf_l(j)**5 * BR

  If (j.Eq.2) Then
   BR = BR / GammaMu
  Else If  (j.Eq.3) Then
   BR = BR / GammaTau
  End If

  Iname = Iname - 1

 Contains

 Complex(dp) Function vertexC0tilde(masa1c,masa2c,masa3c)
   Real(dp), Intent(in) :: masa1c,masa2c,masa3c

   vertexC0tilde = 0._dp

   If (masa2c.Eq.masa3c) Then

   vertexC0tilde=(masa1c*(masa1c-masa2c)-masa1c**2*Log(masa1c)           & 
        & +(2._dp*masa1c-masa2c)*masa2c*Log(masa2c))/(masa1c-masa2c)**2   

   Else

   vertexC0tilde=1._dp-1._dp/(masa2c-masa3c)*((masa1c**2*Log(masa1c)      & 
        & -masa2c**2*Log(masa2c))/(masa1c-masa2c)-(masa1c**2*Log(masa1c)  & 
        & -masa3c**2*Log(masa3c))/(masa1c-masa3c))

   End If

 End Function vertexC0tilde

 Complex(dp) Function vertexC11(masa1c,masa2c,masa3c)
   Real(dp), Intent(in) :: masa1c,masa2c,masa3c

   vertexC11 = 0._dp

   If (masa2c.Eq.masa3c) Then

   vertexC11=-(3._dp*masa1c**2-4._dp*masa1c*masa2c+masa2c**2 & 
        & -2._dp*masa1c**2*Log(masa1c/masa2c)) & 
        & /(2._dp*(masa1c-masa2c)**3)   

   Else

   vertexC11=masa1c*(masa1c**2/masa2c**2-masa1c/masa2c & 
        & -masa1c**2/masa3c**2+masa1c**3/(masa2c*masa3c**2) & 
        & +masa1c/masa3c-masa1c**3/(masa2c**2*masa3c) & 
        & +masa1c**2/masa3c**2*Log(masa1c/masa2c) & 
        & -2._dp*masa1c**3/(masa2c*masa3c**2)*Log(masa1c/masa2c) & 
        & -masa1c**2/masa2c**2*Log(masa1c/masa3c) & 
        & +2._dp*masa1c**3/(masa2c**2*masa3c)*Log(masa1c/masa3c) & 
        & -Log(masa2c/masa3c)+2._dp*masa1c/masa2c*Log(masa2c/masa3c) & 
        & +2._dp*masa1c/masa3c*Log(masa2c/masa3c) & 
        & -4._dp*masa1c**2/(masa2c*masa3c)*Log(masa2c/masa3c)) & 
        & /(2._dp*(1._dp-masa1c/masa2c)**2*masa2c  & 
        & *(1._dp-masa1c/masa3c)**2*(masa1c/masa2c-masa1c/masa3c) & 
        & *masa3c)

   End If

 End Function vertexC11

 Complex(dp) Function vertexC12(masa1c,masa2c,masa3c)
   Real(dp), Intent(in) :: masa1c,masa2c,masa3c

   vertexC12 = 0._dp

   If (masa2c.Eq.masa3c) Then

   vertexC12=-(3._dp*masa1c**2-4._dp*masa1c*masa2c+masa2c**2 & 
        & -2._dp*masa1c**2*Log(masa1c/masa2c)) & 
        & /(4._dp*(masa1c-masa2c)**3)   

   Else

   vertexC12=((masa1c-masa2c)*(masa1c-masa3c)*(masa2c-masa3c)*masa3c      & 
        & +masa1c**2*(masa2c**2*Log(masa1c/masa2c)+masa3c                 & 
        & *(-2._dp*masa2c+masa3c)*Log(masa1c/masa3c))                     & 
        & +masa2c**2*(2._dp*masa1c-masa3c)*masa3c*Log(masa2c/masa3c))     &
        & /(2._dp*(masa1c-masa2c)*(masa1c-masa3c)**2*(masa2c-masa3c)**2)

   End If

 End Function vertexC12

 End Subroutine BR_lj_to_3li

 Real(dp) Function Bm_to_l_nu(i, j, mH2, tanb, RSpm, Yd, RuL, RdR, Yl, vevSM)
 !-------------------------------------------------------------------
 ! calculates the branching ratio of the rare decay B- -> tau nu
 ! based on the formula by W.Hou, PRD48 (1993) 2342 
 !-------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i           ! index of the lepton generation
  Integer, Intent(in) :: j           ! j=1 -> B_u, j=2-> B_c
  Real(dp), Intent(in) :: mH2, tanb  ! m^2_H+, tan(beta)
  Complex(dp), Intent(in), Optional :: RSpm(2,2) & ! H+ mixing matrix
                                   & , Yd(3,3)   & ! d-quark Yukawas
                                   & , RuL(3,3)  & ! L u-quark mixing matrix
                                   & , RdR(3,3)  & ! R d-quark mixing matrix
                                   & , Yl(3,3)     ! lepton Yukawa coupling
  Real(dp), Intent(in), Optional :: vevSM(2)       !

  real(dp) :: r_fac
  Complex(dp) :: cL, cR
  Complex(dp), Parameter :: Zero3(3,3) = (0._dp,0._dp)

  If (    Present(RSpm).And.Present(Yd).And.Present(RuL).And.Present(RdR) &
     .and.Present(Yl).and.Present(vevSM)) then
   call CoupChargedScalarFermion3(2, 3, j, RSpm, Yd, id3C, RdR, Zero3 &
                                    &, RuL, id3C, cL, cR)
   r_fac = 0.25_dp * Abs(cL * Yl(i,i) / CKM(j,3))**2 * vevSM(1)**4 &
         &         / (mf_l2(i) * mf_d_mZ(3)**2 ) 
  else
   r_fac = 1
  end if

   Bm_to_l_nu = oo8pi * G_F**2 * mf_l2(i) * mBm(j)                   &
            &   * (1._dp -  mf_l2(i) / mBm(j)**2)**2                 &
            &   * Abs(CKM(j,3))**2 * FBhatB(j)**2 *  tauBm(j) / hbar    &
            &   * (1._dp - r_fac * (tanb * mBm(j))**2 / mH2)**2

 End Function Bm_to_l_nu


 Subroutine BtoQGamma(I_f, mf_d, gauge, mf_u, mW, Y_d, RdL, RdR, Y_u, RuL, RuR &
   & , mSpm2, RSpm, mC, U, V, mSup2, RSup, A_u, mSdown2, RSdown, A_d, mglu     &
   & , phi_g, mN, N, mu, mS02, RS0, mP02, RP0, vevSM, Bratio, c7_o, c7p_o      &
   & , c8_o, c8p_o, A_CP, i_scheme)
 !-------------------------------------------
 ! gives 10^4 Br(b->s gamma)
 ! input:
 !  I_f .................. if I_f==1 -> Q=d, I_f==1 -> Q=s
 !  mf_d(i) .............. MSbar d-quark masses
 !  gauge(i) ............. U(1), SU(2) and SU(3) MSbar gauge couplings
 !  mf_u(i) .............. MSbar u-quark masses
 !  mW ............... W-boson mass
 !  Y_d(i,j) ............. MSbar d-type Yukawa couplings
 !  RdL(j,j), RdR(i,j) ... mixing matrix of left- and right d-quarks 
 !  Y_u(i,j) ............. MSbar u-type Yukawa couplings
 !  RuL(j,j), RuR(i,j) ... mixing matrix of left- and right u-quarks
 !  mSpm2(i) ............. mass of G+ and H+ squared
 !  RSpm(i,j) ............ mixing matrix of G+ and H+
 !  mC(i) ................ chargino masses
 !  U(i,j), V(i,j) ....... chargino mixing matrices
 !  mSup2(i) ............. u-squark masses squared
 !  RSup(i,j) ............ u-squark mixing matrix
 !  A_u(i,j) ............. trilinear soft SUSY u-squark Higgs coupling
 !  mSdown2(i) ........... d-squark masses squared
 !  RSdown(i,j) .......... d-squark mixing matrix
 !  A_d(i,j) ............. trilinear soft SUSY d-squark Higgs coupling
 !  mglu ................. gluino mass
 !  phi_g ................ phase of M_3
 !  mN(i) ................ neutralino masses
 !  N(i,j) ............... neutralino mixing matrix
 !  mu ................... mu parameter
 !  mS02(i) .............. neutral Higgs masses squared
 !  RS0(i,j) ............. neutral Higgs mixing matrix
 !  mP02(i) .............. neutral pseudoscalar Higgs masses squared
 !  RP0(i,j) ............. neutral pseudoscalar Higgs mixing matrix
 !  vevSM(i) ............. MSSM vevs (v_d, v_u)
 !  i_scheme (optional) ..  1 ... E_0=1.6 GeV, m_c/m_b=0.23 NLO
 !                          2 ... E_0=1.6 GeV, m_c/m_b=0.29 NLO
 !                          3 ... E_0=m_b/20 GeV, m_c/m_b=0.23 NLO
 !                          4 ... E_0=m_b/20 GeV, m_c/m_b=0.29 NLO
 !                          5 ... LO
 ! output:
 !  Bratio ............... 10^4 BR( b -> Q gamma)
 !           i_scheme=i_t:  1 ... E_0=1.6 GeV, m_c/m_b=0.23 NLO
 !                          2 ... E_0=1.6 GeV, m_c/m_b=0.29 NLO
 !                          3 ... E_0=m_b/20 GeV, m_c/m_b=0.23 NLO
 !                          4 ... E_0=m_b/20 GeV, m_c/m_b=0.29 NLO
 !                          5 ... LO
 !  Wilson coefficients at m_W (including the various contributions)
 !                              all are optional
 !  c7_o(i) .............. C_7: 1 ... total
 !                              2 ... SM contribution  
 !                              3 ... H+ contribution  
 !                              4 ... chargino contribution  
 !                              5 ... gluino contribution  
 !                              6 ... neutralino contribution  
 !                              7 ... neutral Higgs contribution  
 !  c7p_o(i) ............. C'_7: 1 ... total
 !                               2 ... H+ contribution  
 !                               3 ... chargino contribution  
 !                               4 ... gluino contribution  
 !                               5 ... neutralino contribution  
 !                               6 ... neutral Higgs contribution  
 !  c8_o(i) .............. C_8: 1 ... total
 !                              2 ... SM contribution  
 !                              3 ... H+ contribution  
 !                              4 ... chargino contribution  
 !                              5 ... gluino contribution  
 !                              6 ... neutralino contribution  
 !                              7 ... neutral Higgs contribution  
 !  c8p_o(i) ............. C'_8: 1 ... total
 !                               2 ... H+ contribution  
 !                               3 ... chargino contribution  
 !                               4 ... gluino contribution  
 !                               5 ... neutralino contribution  
 !                               6 ... neutral Higgs contribution  
 ! written by Werner Porod, 11.11.02
 ! 13.03.02: -adding the possiblity to switch between d and s by the switch if
 !           -adding gluino and neutralino contributions
 !           -creating interface for neutral Higgs boson, are set to 0 for
 !            the moment being
 ! 02.10.03: changing numbers in front of R_7,8 as given in coll. with 
 !           Enrico and Tobias, adding i_scheme for the various possiblities
 !           of the a^x_y - arrays, default = 4
 !           i_scheme=i_t:  1 ... E_0=1.6 GeV, m_c/m_b=0.23 NLO
 !                          2 ... E_0=1.6 GeV, m_c/m_b=0.29 NLO
 !                          3 ... E_0=m_b/20 GeV, m_c/m_b=0.23 NLO
 !                          4 ... E_0=m_b/20 GeV, m_c/m_b=0.29 NLO
 !                          5 ... LO
 !-------------------------------------------
 Implicit None
  Integer, Intent(in) :: I_f
  Real(dp), Intent(in) :: mf_d(3), gauge(3), mf_u(3), mW, mSpm2(2), mC(2)   &
    & , mSup2(6), mSdown2(6), mglu, mN(4), mS02(2), RS0(2,2), mP02(2)       &
    & , RP0(2,2), vevSM(2)
  Complex(dp), Intent(in), Dimension(3,3) :: RdL, RdR, RuL, RuR, Y_d, Y_u
  Complex(dp), Intent(in) :: U(2,2), V(2,2), RSpm(2,2), RSup(6,6)    &
    & , RSdown(6,6), phi_g, N(4,4), mu, A_u(3,3), A_d(3,3)
  Real(dp), Intent(out) :: Bratio
  Complex(dp), Intent(out), Optional :: c7_o(7), c7p_o(6), c8_o(7), c8p_o(6)
  Real(dp), Intent(out), Optional :: A_CP
  Integer, Intent(in), Optional :: i_scheme

  Integer :: i1, i2, i3, i_t
  Real(dp) :: g, gp, gs, mSup2_in(6)
  Complex(dp) :: c7(7), c7p(6), c8(7), c8p(6), r7, r7p, r8, r8p, cpl_uWd(3,3) &
               & , cpl_CSQQp_L(2,3,3) , cpl_CSQQp_R(2,3,3), cpl_CDSu_L(2,3,6) &
               & , cpl_CDSu_R(2,3,6), RSup_in(6,6), mix(6,6), cpl_DGSd_L(3,6) &
               & , cpl_DGSd_R(3,6), cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6)      &
               & , cpl_DDS0_L(3,3,2), cpl_DDS0_R(3,3,2), cpl_DDP0_L(3,3,2)    &
               & , cpl_DDP0_R(3,3,2), norm, CKM(3,3), epsq, epsqC
  Complex(dp) :: delta_C7_0, delta_C7p_0, delta_C8_0, delta_C8p_0, delta_C7_1 &
               & , delta_C7p_1, delta_C8_1, delta_C8p_1
  Real(dp), Parameter :: T3=-0.5_dp, ed = -1._dp / 3._dp
  Real(dp), Dimension(5), Parameter ::                                         &
    &   a = (/  7.8221_dp, 6.9120_dp, 8.1819_dp, 7.1714_dp, 7.9699_dp /)       &
    & , a_77 = (/ 0.8161_dp, 0.8161_dp, 0.8283_dp, 0.8283_dp, 0.9338_dp /)     &
    & , ar_7 = (/ 4.8802_dp, 4.5689_dp, 4.9228_dp, 4.6035_dp, 5.3314_dp /)     &
    & , ai_7 = (/ 0.3546_dp, 0.2167_dp, 0.3322_dp, 0.2029_dp, 0.0_dp    /)     &
    & , a_88 = (/ 0.0197_dp, 0.0197_dp, 0.0986_dp, 0.0986_dp, 0.0066_dp /)     &
    & , ar_8 = (/ 0.5680_dp, 0.5463_dp, 0.7810_dp, 0.7600_dp, 0.4498_dp /)     &
    & , ai_8 = (/-0.0987_dp,-0.1105_dp,-0.0963_dp,-0.1091_dp, 0.0_dp /)        &
    & , a_eps = (/ 0.4384_dp, 0.3787_dp, 0.8598_dp, 0.7097_dp, 0.0_dp /)       &
    & , ar_eps = (/-1.6981_dp,-2.6679_dp,-1.3329_dp,-2.4935_dp, 0.0_dp /)      &
    & , ai_eps = (/ 2.4997_dp, 2.8956_dp, 2.5274_dp, 2.9127_dp, 0.0_dp /)      &
    & , ar_87 = (/ 0.1923_dp, 0.1923_dp, 0.2025_dp, 0.2025_dp, 0.1576_dp /)    &
    & , ar_7eps = (/-0.7827_dp,-1.0940_dp,-0.8092_dp,-1.1285_dp, 0.0_dp /)     &
    & , ar_8eps = (/ -0.0601_dp,-0.0819_dp,-0.0573_dp,-0.0783_dp, 0.0_dp /)    &
    & , ai_87 = (/-0.0487_dp,-0.0487_dp,-0.0487_dp,-0.0487_dp, 0.0_dp /)       &
    & , ai_7eps = (/-0.9067_dp,-1.0447_dp,-0.9291_dp,-1.0585_dp, 0.0_dp /)     &
    & , ai_8eps = (/-0.0661_dp,-0.0779_dp,-0.0637_dp,-0.0765_dp, 0.0_dp /)

  Iname = Iname + 1
  NameOfUnit(Iname) = "BtoQGamma"

  !---------------
  ! couplings
  !---------------
  gp = gauge(1)
  g = gauge(2)
  gs = gauge(3)

  CKM =  Matmul(RuL, Transpose(Conjg(RdL)) )

  cpl_uWd = g * oosqrt2 * CKM

  norm = 0.5_dp * cpl_uWd(3,3) * Conjg(cpl_uWd(3,i_f)) / mW**2

  cpl_CSQQp_L = 0._dp
  cpl_CSQQp_R = 0._dp
  Do i1=1,3
   Do i2=1,3
    Call CoupChargedScalarFermion3(2, i1, i2, RSpm, Y_d, RdL, RdR, Y_U &
                     &, RuL, RuR, cpl_CSQQp_L(2,i1,i2), cpl_CSQQp_R(2,i1,i2) )
!    Write(*,*) "Cs",i1,i2, cpl_CSQQp_L(2,i1,i2), cpl_CSQQp_R(2,i1,i2)
   End Do
  End Do
  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp
  If (GenerationMixing) Then
   RSup_in = RSup
   mSup2_in = mSup2
  Else ! need to reorder things to be in agreement with ordering of scalar ups
   RSup_in = 0._dp
   RSup_in(1,1) = RSup(1,1)
   RSup_in(1,4) = RSup(1,2)
   RSup_in(4,1) = RSup(2,1)
   RSup_in(4,4) = RSup(2,2)
   RSup_in(2,2) = RSup(3,3)
   RSup_in(2,5) = RSup(3,4)
   RSup_in(5,2) = RSup(4,3)
   RSup_in(5,5) = RSup(4,4)
   RSup_in(3,3) = RSup(5,5)
   RSup_in(3,6) = RSup(5,6)
   RSup_in(6,3) = RSup(6,5)
   RSup_in(6,6) = RSup(6,6)
   mix = 0._dp
   mix(1:3,1:3) = rUL
   mix(4:6,4:6) = rUr
   RSup_in = Matmul(RSup_in, mix)
   mSup2_in(1) = mSup2(1)
   mSup2_in(2) = mSup2(3)
   mSup2_in(3) = mSup2(5)
   mSup2_in(4) = mSup2(2)
   mSup2_in(5) = mSup2(4)
   mSup2_in(6) = mSup2(6)
  End If   
   Do i2=1,2
    Do i3=1,6
     Call CoupCharginoSfermion3(i2, 3, i3, g, -0.5_dp, RSup_in, Y_D, Y_U &
           &, RdL, RdR, U, V, cpl_CDSu_L(i2, 3, i3), cpl_CDSu_R(i2, 3, i3) )
     Call CoupCharginoSfermion3(i2, I_f, i3, g, -0.5_dp, RSup_in, Y_D, Y_U &
           &, RdL, RdR, U, V, cpl_CDSu_L(i2, I_f, i3), cpl_CDSu_R(i2, I_f, i3) )
    End Do
   End Do

  cpl_DGSd_L = 0._dp
  cpl_DGSd_R = 0._dp
  cpl_DNSd_L = 0._dp
  cpl_DNSd_R = 0._dp
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  If (GenerationMixing) Then
   Do i1=1,6
    Call CoupGluinoSquark3(gs, phi_g, 3, i1, Rsdown, RdL, RdR &
                         & , cpl_DGSd_L(3,i1), cpl_DGSd_R(3,i1) )
    Call CoupGluinoSquark3(gs, phi_g, I_f, i1, Rsdown, RdL, RdR &
                         & , cpl_DGSd_L(I_f,i1), cpl_DGSd_R(I_f,i1) )
    Do i2=1,4
     Call CoupNeutralinoSfermion3(3, i2, i1, gp, g, T3, ed, RSdown, RdL, RdR &
                        & , Y_d, N, cpl_DNSd_L(3,i2,i1), cpl_DNSd_R(3,i2,i1))
     Call CoupNeutralinoSfermion3(I_f, i2, i1, gp, g, T3, ed, RSdown, RdL, RdR &
                       & , Y_d, N, cpl_DNSd_L(I_f,i2,i1), cpl_DNSd_R(I_f,i2,i1))
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     Do i3=1,2
    Call  CoupFermionScalar31L_eff(i1, i2, i3, T3, ed, gauge, Y_d, RdL, RdR     &
     & , RS0, vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2_in     &
     & , RSup, Y_u, A_u, mC, U, V, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
    Call CoupFermionPseudoScalar31L_eff(i1, i2, i3, T3, ed, gauge, Y_d, RdL     &
    & , RdR, RP0, vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2_in &
    & , RSup, Y_u, A_u, mC, U, V, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
  End If
  !---------------------
  ! wilson coefficients
  !---------------------
  Call C_7(3, I_f, mW, mf_u, cpl_uWd, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R  &
        & , mC, mSup2_in, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L  &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L        &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c7 )
  Call C_7p(3, I_f, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R               &
         & , mC, mSup2_in, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L  &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L         &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c7p  )
  Call C_8(3, I_f, mW, mf_u, cpl_uWd  &
        & , mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R                             &
         & , mC, mSup2_in, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L  &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L         &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c8)
  Call C_8p(3, I_f, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R               &
         & , mC, mSup2_in, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L  &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L         &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c8p )

  c7 = c7 / norm
  c7p = c7p / norm
  c8 = c8 / norm
  c8p = c8p / norm

  r7 = c7(1) / c7(2)
  r7p = c7p(1) / c7(2)
  r8 = c8(1) / c8(2)
  r8p = c8p(1) / c8(2)

  epsq = Conjg(CKM(1,I_f)) * CKM(1,3) / (Conjg(CKM(3,I_f)) * CKM(3,3))
  epsqC = Conjg(epsq)

  If (Present(i_scheme)) Then
   i_t = i_scheme
  Else
   i_t = 3
  End If

  if (i_t.ne.0) then
   Bratio = a(i_t) + a_77(i_t) * ( Abs(r7)**2 + Abs(r7p)**2)                    &
       & + ar_7(i_t) * Real(r7,dp) + ai_7(i_t) * Aimag(r7)                     &
       & + a_88(i_t) * ( Abs(r8)**2 + Abs(r8p)**2)                             &
       & + ar_8(i_t) * Real(r8,dp) + ai_8(i_t) * Aimag(r8)                     &
       & + a_eps(i_t) * Abs(epsq)**2 + ar_eps(i_t) * Real(epsq,dp)             &
       & + ai_eps(i_t) * Aimag(epsq)                                           &
       & + ar_87(i_t) * Real(r8*Conjg(r7)+r8p*Conjg(r7p),dp)                   &
       & + ai_87(i_t) * Aimag(r8*Conjg(r7)+r8p*Conjg(r7p))                     &
       & + ar_7eps(i_t) * Real(r7*epsqC,dp) + ar_8eps(i_t) * Real(r8*epsqC,dp) &
       & + ai_7eps(i_t) * Aimag(r7*epsqC) + ai_8eps(i_t) * Aimag(r8*epsqC) 

  !-----------------------
  ! is 10^4 BR
  !-----------------------
   Bratio = 2.567e-1_dp * Abs(Conjg(CKM(3,I_f))*CKM(3,3)/CKM(2,3))**2 * Bratio

  Else ! E.Lunghi, J.Matias, hep-ph/0612166, eq.18 still missing 1-loop part
   delta_C7_0 = sum(C7(3:7))
   delta_C7p_0 = C7p(1)
   delta_C8_0 = sum(C8(3:7))
   delta_C8p_0 = C8p(1)
   delta_C7_1 = 0._dp
   delta_C7p_1 = 0._dp
   delta_C8_1 = 0._dp
   delta_C8p_1 = 0._dp
   Bratio = 2.98_dp + 4.743_dp * (Abs(delta_C7_0)**2 + Abs(delta_C7p_0)**2 ) &
        & + 0.789_dp *  (Abs(delta_C8_0)**2 + Abs(delta_C8p_0)**2 )          &
        & + Real( (-7.184_dp,0.612_dp) * delta_C7_0                          &
        &       + (-2.225_dp,-0.557_dp) * delta_C8_0                         &
        &       + (2.454_dp,-0.884_dp) * ( delta_C7_0 * conjg(delta_C8_0)    &
        &                                +delta_C7p_0 * conjg(delta_C8p_0) ) &
        &       , dp)  
  end if

  If (Present(A_CP)) Then
   A_CP = Aimag( ai_7(i_t)*(r7) + ai_8(i_t)*(r8)                       &
      &        + ai_87(i_t) * (r8*Conjg(r7)+r8p*Conjg(r7p))            &
      &        + ai_7eps(i_t) * r7*epsqC + ai_8eps(i_t) * r8*epsqC )   &
      & / ( a(i_t) + a_77(i_t) * ( Abs(r7)**2 + Abs(r7p)**2)           &
      &   + ar_7(i_t) * Real(r7,dp) + ar_8(i_t) * Real(r8,dp)          &
      &   + a_88(i_t) * ( Abs(r8)**2 + Abs(r8p)**2)                    &
      &   + ar_87(i_t) * Real(r8*Conjg(r7)+r8p*Conjg(r7p),dp) )
  End If

  If (Present(c7_o)) c7_o = c7
  If (Present(c7p_o)) c7p_o = c7p
  If (Present(c8_o)) c8_o = c8
  If (Present(c8p_o)) c8p_o = c8p

  Iname = Iname - 1

 End Subroutine BtoQGamma

 Subroutine BToSLL(gauge, mf_d, mf_u, mW, Y_d, RdL, RdR, Y_u, RuL, RuR     &
   & , Y_l, RlL, RlR, mSneut2, Rsneut, mSlepton2, Rslepton                 &
   & , mSpm2, RSpm, mC, U, V, mSup2, RSup, A_u, mSdown2, RSdown, A_d, mglu &
   & , phi_g, mN, N, mu, mS02, RS0, mP02, RP0, vevSM, BtoSEE, BtoSMuMu     &
   & , c7_o, c7p_o, c8_o, c8p_o, c9_o, c9p_o, c10_o, c10p_o)
 !-------------------------------------------
 ! gives BR(b->s e+ e-) and BR(b->s e+ e-) using formulas 12/13 of
 ! T. Huber et al., hep-ph/0512066
 ! input: all running quantities are given at Q=M_W
 !  gauge(i) ............. U(1), SU(2) and SU(3) MSbar gauge couplings
 !  mf_d(i) .............. MSbar d-quark masses
 !  mf_u(i) .............. MSbar u-quark masses
 !  mW ............... W-boson mass
 !  Y_d(i,j) ............. MSbar d-type Yukawa couplings
 !  RdL(j,j), RdR(i,j) ... mixing matrix of left- and right d-quarks 
 !  Y_u(i,j) ............. MSbar u-type Yukawa couplings
 !  RuL(j,j), RuR(i,j) ... mixing matrix of left- and right u-quarks
 !  Y_l(i,j) ............. MSbar lepton Yukawa couplings
 !  RlL(j,j), RlR(i,j) ... mixing matrix of left- and right leptons
 !  mSneut2(i) ........... sneutrino masses squared
 !  RSneut(i,j) .......... sneutrino mixing matrix
 !  mSlepton2(i) ......... slepton masses squared
 !  RSlepton(i,j) ........ slepton mixing matrix
 !  mSpm2(i) ............. mass of G+ and H+ squared
 !  RSpm(i,j) ............ mixing matrix of G+ and H+
 !  mC(i) ................ chargino masses
 !  U(i,j), V(i,j) ....... chargino mixing matrices
 !  mSup2(i) ............. u-squark masses squared
 !  RSup(i,j) ............ u-squark mixing matrix
 !  A_u(i,j) ............. trilinear soft SUSY u-squark Higgs coupling
 !  mSdown2(i) ........... d-squark masses squared
 !  RSdown(i,j) .......... d-squark mixing matrix
 !  A_d(i,j) ............. trilinear soft SUSY d-squark Higgs coupling
 !  mglu ................. gluino mass
 !  phi_g ................ phase of M_3
 !  mN(i) ................ neutralino masses
 !  N(i,j) ............... neutralino mixing matrix
 !  mu ................... mu parameter
 !  mS02(i) .............. neutral Higgs masses squared
 !  RS0(i,j) ............. neutral Higgs mixing matrix
 !  mP02(i) .............. neutral pseudoscalar Higgs masses squared
 !  RP0(i,j) ............. neutral pseudoscalar Higgs mixing matrix
 !  vevSM(i) ............. MSSM vevs (v_d, v_u)
 ! output:
 !  BtoSEE ............... BR(b -> s e+ e-) 
 !  BtoSMuMu ............. BR(b -> s mu+ mu-)
 !  Wilson coefficients at m_W (including the various contributions)
 !                              all are optional
 !  c7_o(i) .............. C_7: 1 ... total
 !                              2 ... SM contribution  
 !                              3 ... H+ contribution  
 !                              4 ... chargino contribution  
 !                              5 ... gluino contribution  
 !                              6 ... neutralino contribution  
 !                              7 ... neutral Higgs contribution  
 !  c7p_o(i) ............. C'_7: 1 ... total
 !                               2 ... H+ contribution  
 !                               3 ... chargino contribution  
 !                               4 ... gluino contribution  
 !                               5 ... neutralino contribution  
 !                               6 ... neutral Higgs contribution  
 !  c8_o(i) .............. C_8: 1 ... total
 !                              2 ... SM contribution  
 !                              3 ... H+ contribution  
 !                              4 ... chargino contribution  
 !                              5 ... gluino contribution  
 !                              6 ... neutralino contribution  
 !                              7 ... neutral Higgs contribution  
 !  c8p_o(i) ............. C'_8: 1 ... total
 !                               2 ... H+ contribution  
 !                               3 ... chargino contribution  
 !                               4 ... gluino contribution  
 !                               5 ... neutralino contribution  
 !                               6 ... neutral Higgs contribution  
 !  c9_o(i) .............. C_9: 1 ... total
 !                              2 ... SM contribution  
 !                              3 ... H+ contribution  
 !                              4 ... chargino contribution  
 !                              5 ... neutralino contribution  
 !                              6 ... gluino contribution  
 !  c9p_o(i) ............. C'_9: 1 ... total
 !                               2 ... SM contribution  
 !                               3 ... H+ contribution  
 !                               4 ... chargino contribution  
 !                               5 ... neutralino contribution  
 !                               6 ... gluino contribution  
 ! written by Werner Porod, 07.12.05
 !-------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_d(3), gauge(3), mf_u(3), mW, mSpm2(2), mC(2)   &
    & , mSup2(6), mSdown2(6), mglu, mN(4), mS02(2), RS0(2,2), mP02(2)       &
    & , RP0(2,2), vevSM(2), mSneut2(3), mSlepton2(6)
  Complex(dp), Intent(in), Dimension(3,3) :: RdL, RdR, RuL, RuR, Y_d, Y_u   &
    & , RlL, RlR, Y_l
  Complex(dp), Intent(in) :: U(2,2), V(2,2), RSpm(2,2), RSup(6,6), Rsneut(3,3) &
    & , RSdown(6,6), phi_g, N(4,4), mu, A_u(3,3), A_d(3,3), Rslepton(6,6)
  Real(dp), Intent(out) :: BtoSEE, BToSMuMu
  Complex(dp), Intent(out), Optional :: c7_o(7), c7p_o(6), c8_o(7), c8p_o(6) &
    & , c9_o(2,6), c9p_o(2,6), c10_o(2,6), c10p_o(2,6)

  Integer :: i1, i2, i3, i4
  Real(dp) :: g, gp, gs
  Complex(dp) :: c7(7), c7p(6), c8(7), c8p(6), r7, r7p, r8, r8p, cpl_uWd(3,3) &
               & , cpl_CSQQp_L(2,3,3) , cpl_CSQQp_R(2,3,3), cpl_CDSu_L(2,3,6) &
               & , cpl_CDSu_R(2,3,6), RSup_in(6,6), mix(6,6), cpl_DGSd_L(3,6) &
               & , cpl_DGSd_R(3,6), cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6)      &
               & , cpl_DDS0_L(3,3,2), cpl_DDS0_R(3,3,2), cpl_DDP0_L(3,3,2)    &
               & , cpl_DDP0_R(3,3,2), norm, CKM(3,3)
  Complex(dp) :: r9(2), r9p(2), r10(2), r10p(2), c9(2), c9p(2)         &
    & , c10(2), c10p(2), c9_c(2,6), c9p_c(2,6), c10_c(2,6), c10p_c(2,6)
  Complex(dp) :: kappa_q, cpl_CLSn_L(2,3,3), cpl_CLSn_R(2,3,3)    &
             & , cpl_LNSl_L(3,4,6), cpl_LNSl_R(3,4,6)
  Real(dp) :: mC2(2), mN2(4), a_s, sW2, tanb, e2, mG2, mW2
  Real(dp), Parameter :: T3=-0.5_dp, ed = -1._dp / 3._dp
  Complex(dp), Parameter :: Null3(3,3) = 0._dp

  Iname = Iname + 1
  NameOfUnit(Iname) = "BToSLL"

  !---------------
  ! couplings
  !---------------
  gp = gauge(1)
  g = gauge(2)
  gs = gauge(3)

  CKM =  Matmul(RuL, Transpose(Conjg(RdL)) )

  cpl_uWd = g * oosqrt2 * CKM
  mW2 = mW**2

  norm = 0.5_dp * cpl_uWd(3,3) * Conjg(cpl_uWd(3,2)) / mW2

  cpl_CSQQp_L = 0._dp
  cpl_CSQQp_R = 0._dp
  Do i1=1,3
   Do i2=1,3
    Call CoupChargedScalarFermion3(2, i1, i2, RSpm, Y_d, RdL, RdR, Y_U &
                     &, RuL, RuR, cpl_CSQQp_L(2,i1,i2), cpl_CSQQp_R(2,i1,i2) )
   End Do
  End Do

  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp
  If (GenerationMixing) Then
   RSup_in = RSup
  Else ! need to reorder things to be in agreement with ordering of scalar ups
   RSup_in = 0._dp
   RSup_in(1,1) = RSup(1,1)
   RSup_in(1,4) = RSup(1,2)
   RSup_in(4,1) = RSup(2,1)
   RSup_in(4,4) = RSup(2,2)
   RSup_in(2,2) = RSup(3,3)
   RSup_in(2,5) = RSup(3,4)
   RSup_in(5,2) = RSup(4,3)
   RSup_in(5,5) = RSup(4,4)
   RSup_in(3,3) = RSup(5,5)
   RSup_in(3,6) = RSup(5,6)
   RSup_in(6,3) = RSup(6,5)
   RSup_in(6,6) = RSup(6,6)
   mix = 0._dp
   mix(1:3,1:3) = rUL
   mix(4:6,4:6) = rUr
   RSup_in = Matmul(RSup_in, mix)
  End If   
  Do i1=2,3
   Do i2=1,2
    Do i3=1,6
     Call CoupCharginoSfermion3(i2, i1, i3, g, -0.5_dp, RSup_in, Y_D, Y_U &
           &, RdL, RdR, U, V, cpl_CDSu_L(i2, i1, i3), cpl_CDSu_R(i2, i1, i3) )
    End Do
   End Do
  End Do

  cpl_DGSd_L = 0._dp
  cpl_DGSd_R = 0._dp
  cpl_DNSd_L = 0._dp
  cpl_DNSd_R = 0._dp
  cpl_DDS0_L = 0._dp
  cpl_DDS0_R = 0._dp
  cpl_DDP0_L = 0._dp
  cpl_DDP0_R = 0._dp
  cpl_LNSl_L = 0._dp
  cpl_LNSl_R = 0._dp
  cpl_CLSn_L = 0._dp
  cpl_CLSn_R = 0._dp
  If (GenerationMixing) Then
   Do i1=1,6
    Do i3=2,3
     Call CoupGluinoSquark3(gs, phi_g, i3, i1, Rsdown, RdL, RdR     &
                          & , cpl_DGSd_L(i3,i1), cpl_DGSd_R(i3,i1) )
     Do i2=1,4
      Call CoupNeutralinoSfermion3(i3, i2, i1, gp, g, T3, ed, RSdown, RdL, RdR &
                        & , Y_d, N, cpl_DNSd_L(i3,i2,i1), cpl_DNSd_R(i3,i2,i1))
     End Do
    End Do
   End Do
   Do i1=2,3
    Do i2=2,3
     Do i3=1,2
    Call  CoupFermionScalar31L_eff(i1, i2, i3, T3, ed, gauge, Y_d, RdL, RdR    &
     & , RS0, vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup &
     & , Y_u, A_u, mC, U, V, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
    Call CoupFermionPseudoScalar31L_eff(i1, i2, i3, T3, ed, gauge, Y_d, RdL    &
     & , RdR, RP0, vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2  &
     & , RSup, Y_u, A_u, mC, U, V, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
   Do i1=1,3
    Do i2=1,4
     Do i3=1,6  
      Call CoupNeutralinoSlepton(i1, i2, i3, gauge(1), gauge(2), RSlepton &
          & , RlL, RlR, Y_l, N, cpl_LNSl_L(i1,i2,i3), cpl_LNSl_R(i1,i2,i3) )
     End Do
    End Do
    Do i2=1,2
     Do i3=1,3
     Call CoupCharginoSfermion3(i2, i1, i3, gauge(2), -0.5_dp, RSneut, Y_l &
        & , Null3, id3C, id3C, U, V, cpl_CLSn_L(i2,i1,i3), cpl_CLSn_R(i2,i1,i3))
     End Do
    End Do
   End Do

  Else ! no generation mixing

   Do i1=2,3
    Do i2=2,3
     Do i3=1,2
    Call  CoupFermionScalar31L_eff(i1, i2, i3, T3, ed, gauge, Y_d, RdL, RdR    &
     & , RS0, vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup &
     & , Y_u, A_u, mC, U, V, cpl_DDS0_L(i1,i2,i3), cpl_DDS0_R(i1,i2,i3))
    Call CoupFermionPseudoScalar31L_eff(i1, i2, i3, T3, ed, gauge, Y_d, RdL    &
     & , RdR, RP0, vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2  &
     & , RSup, Y_u, A_u, mC, U, V, cpl_DDP0_L(i1,i2,i3), cpl_DDP0_R(i1,i2,i3))
     End Do
    End Do
   End Do
   Do i1=1,3
    Do i2=1,4
     Do i3=1,2
      i4 = (i1-1)*2 + i3  
      Call CoupNeutralinoSlepton(i1, i2, i4, gauge(1), gauge(2), RSlepton &
          & , RlL, RlR, Y_l, N, cpl_LNSl_L(i1,i2,i4), cpl_LNSl_R(i1,i2,i4) )
     End Do
    End Do
    Do i2=1,2
     Call CoupCharginoSfermion3(i2, i1, i1, gauge(2), -0.5_dp, RSneut, Y_l &
        & , Null3, id3C, id3C, U, V, cpl_CLSn_L(i2,i1,i1), cpl_CLSn_R(i2,i1,i1))
    End Do
   End Do

  End If
  !---------------------
  ! wilson coefficients
  !---------------------
  Call C_7(3, 2, mW, mf_u, cpl_uWd, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R  &
        & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L   &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L      &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c7 )
  Call C_7p(3, 2, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R              &
         & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L  &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L      &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c7p  )
  Call C_8(3, 2, mW, mf_u, cpl_uWd                                         &
        & , mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R                          &
         & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L  &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L      &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c8)
  Call C_8p(3, 2, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R              &
         & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L  &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L      &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, c8p )

  mC2 = mC**2
  mN2 = mN**2
  mG2 = mGlu**2

  a_s = gauge(3)**2 * oo4pi
  sW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  tanb = vevSM(2) / vevSM(1)
  e2 = 4._dp * pi * gauge(2)**2 * sW2
  kappa_q = 8._dp * Sqrt2 * G_F * e2 * CKM(3,3) * Conjg( CKM(3,2) )
  kappa_q = 1._dp / kappa_q
  Call C_9(3, 2, 1, mW, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2, mN2        &
       & , mSneut2, mSlepton2, mSup2, mSdown2, tanb, U, V, N, RSup, RSdown    &
       & , cpl_CDSu_L, cpl_CDSu_R, cpl_CLSn_L, cpl_CLSn_R, cpl_LNSl_L         &
       & , cpl_LNSl_R, cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R         &
       & , a_s, kappa_q, e2, sW2, .False.                                     &
       & , c9(1), c9p(1), c9_c(1,:), c9p_c(1,:) )
  Call C_9(3, 2, 2, mW, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2, mN2        &
       & , mSneut2, mSlepton2, mSup2, mSdown2, tanb, U, V, N, RSup, RSdown    &
       & , cpl_CDSu_L, cpl_CDSu_R, cpl_CLSn_L, cpl_CLSn_R, cpl_LNSl_L         &
       & , cpl_LNSl_R, cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R         &
       & , a_s, kappa_q, e2, sW2, .False.                                     &
       & , c9(2), c9p(2), c9_c(2,:), c9p_c(2,:) )
  Call C_10_10p(3, 2, 1, mW, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2, mN2   &
       & , mSneut2, mSlepton2, mSup2, mSdown2, tanb, U, V, N, RSup, RSdown    &
       & , cpl_CDSu_L, cpl_CDSu_R, cpl_CLSn_L, cpl_CLSn_R, cpl_LNSl_L         &
       & , cpl_LNSl_R, cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R         &
       & , a_s, kappa_q, e2, sW2, .False., c10(1), c10p(1), c10_c(1,:)        &
       & , c10p_c(1,:))
  Call C_10_10p(3, 2, 2, mW, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2, mN2   &
       & , mSneut2, mSlepton2, mSup2, mSdown2, tanb, U, V, N, RSup, RSdown    &
       & , cpl_CDSu_L, cpl_CDSu_R, cpl_CLSn_L, cpl_CLSn_R, cpl_LNSl_L         &
       & , cpl_LNSl_R, cpl_DGSd_L, cpl_DGSd_R, cpl_DNSd_L, cpl_DNSd_R         &
       & , a_s, kappa_q, e2, sW2, .False., c10(2), c10p(2), c10_c(2,:)        &
       & , c10p_c(2,:))

  c7 = c7 / norm
  c7p = c7p / norm
  c8 = c8 / norm
  c8p = c8p / norm

  r7 = c7(1) / c7(2)
  r7p = c7p(1) / c7(2)
  r8 = c8(1) / c8(2)
  r8p = c8p(1) / c8(2)

  r9 = c9 / c9_c(:,2)
  r9p = c9p / c9_c(:,2)
  r10 = c10 / c10_c(:,2)
  r10p = c10p / c10_c(:,2)

  BtoSEE = (2.3278_dp - 1.655e-3_dp * Aimag(R10(1))                            &
       & + 5.e-4_dp * Aimag(r10(1)*Conjg(r8) + r10p(1)*Conjg(r8p) )            &
       & + 5.24e-2_dp * Aimag(r7) + 5.04e-3_dp * Aimag(r8)                     &
       & + 2.266e-2_dp * Aimag(r7 * Conjg(r8) + r7p * Conjg(r8p) )             &
       & + 4.96e-3_dp * Aimag(r7 * Conjg(r9(1)) + r7p * Conjg(r9p(1)) )        &
       & + 2.61e-2_dp * Aimag(r8 * Conjg(r9(1)) + r8p * Conjg(r9p(1)) )        &
       & - 6.51e-3_dp * Aimag(r9(1)) - 0.5426_dp * Real( r10(1), dp)           &
       & - 2.578e-2_dp * Real(r7,dp) - 1.28e-2_dp * Real(r8,dp)                &
       & + 1.53e-2_dp * Real(r7*Conjg(r10(1)) + r7p*Conjg(r10p(1)), dp )       &
       & + 6.74e-2_dp * Real(r7 * Conjg(r8) + r7p * Conjg(r8p), dp )           &
       & - 0.86996_dp * Real(r7*Conjg(r9(1)) + r7p*Conjg(r9p(1)), dp )         &
       & + 1.85e-3_dp * Real(r8*Conjg(r10(1)) + r8p*Conjg(r10p(1)), dp )       &
       & - 9.926e-2_dp * Real(r8*Conjg(r9(1)) + r8p*Conjg(r9p(1)), dp )        &
       & + 2.841_dp * Real(r9(1),dp) + 0.2813_dp * (Abs(r7)**2 + Abs(r7p)**2 ) &
       & - 0.10705_dp * Real( r9(1) * Conjg(r10(1))                            &
       &                    + r9p(1) * Conjg(r10p(1)), dp)                     &
       & + 11.0367_dp * (Abs(r10(1))**2 + Abs(r10p(1))**2 )                    &
       & + 1.528_dp * (Abs(r9(1))**2 + Abs(r9p(1))**2 )                        &
       & + 3.765e-3_dp * (Abs(r8)**2 + Abs(r8p)**2 ) ) * 1.e-7_dp

  BtoSMuMu = (2.1913_dp - 1.655e-3_dp * Aimag(R10(2))                          &
      & + 5.e-4_dp * Aimag(r10(2)*Conjg(r8) + r10p(2)*Conjg(r8p) )             &
      & + 5.35e-2_dp * Aimag(r7) + 5.13e-3_dp * Aimag(r8)                      &
      & + 2.266e-2_dp * Aimag(r7 * Conjg(r8) + r7p * Conjg(r8p) )              &
      & + 4.96e-3_dp * Aimag(r7 * Conjg(r9(2)) + r7p * Conjg(r9p(2)) )         &
      & + 2.61e-2_dp * Aimag(r8 * Conjg(r9(2)) + r8p * Conjg(r9p(2)) )         &
      & - 1.18e-2_dp * Aimag(r9(2)) - 0.5426_dp * Real( r10(2), dp)            &
      & + 2.81e-2_dp * Real(r7,dp) - 8.66e-3_dp * Real(r8,dp)                  &
      & + 1.53e-2_dp * Real(r7*Conjg(r10(2)) + r7p*Conjg(r10p(2)), dp )        &
      & + 6.859e-2_dp * Real(r7 * Conjg(r8) + r7p * Conjg(r8p), dp )           &
      & - 0.8554_dp * Real(r7*Conjg(r9(2)) + r7p*Conjg(r9p(2)), dp )           &
      & + 1.85e-3_dp * Real(r8*Conjg(r10(2)) + r8p*Conjg(r10p(2)), dp )        &
      & - 9.81e-2_dp * Real(r8*Conjg(r9(2)) + r8p*Conjg(r9p(2)), dp )          &
      & + 2.7008_dp * Real(r9(2),dp) + 0.2889_dp * (Abs(r7)**2 + Abs(r7p)**2 ) &
      & - 0.10705_dp * Real( r9(2) * Conjg(r10(2))                             &
      &                    + r9p(2) * Conjg(r10p(2)), dp)                      &
      & + 10.7687_dp * (Abs(r10(2))**2 + Abs(r10p(2))**2 )                     &
      & + 1.4892_dp * (Abs(r9(2))**2 + Abs(r9p(2))**2 )                        &
      & + 3.81e-3_dp * (Abs(r8)**2 + Abs(r8p)**2 ) ) * 1.e-7_dp

  If (Present(c7_o)) c7_o = c7
  If (Present(c7p_o)) c7p_o = c7p
  If (Present(c8_o)) c8_o = c8
  If (Present(c8p_o)) c8p_o = c8p
  If (Present(c9_o)) c9_o = c9_c
  If (Present(c9p_o)) c9p_o = c9p_c
  If (Present(c10_o)) c10_o = c10_c
  If (Present(c10p_o)) c10p_o = c10p_c

  Iname = Iname - 1

 End Subroutine BToSLL

  Subroutine B_To_SNuNu(gauge, mf_d, mf_u, mW, Y_d, RdL, RdR, Y_u, RuL, RuR &
   & , Y_l, mSneut2, Rsneut, mSlepton2, Rslepton, mSpm2, RSpm, mC, U, V     &
   & , mSup2, RSup, mSdown2, RSdown, mglu, phi_g, mN, N, vevSM, l_QCD       &
   & , BtoSNuNu, c11_o, c11p_o)
 !-------------------------------------------
 ! gives BR(b->s nu nu) using formulas 5.10 of
 ! Bobeth et al., NPB630 (2002) 87
 ! input: all running quantities are given at Q=M_W
 !  gauge(i) ............. U(1), SU(2) and SU(3) MSbar gauge couplings
 !  mf_d(i) .............. MSbar d-quark masses
 !  mf_u(i) .............. MSbar u-quark masses
 !  mW ............... W-boson mass
 !  Y_d(i,j) ............. MSbar d-type Yukawa couplings
 !  RdL(j,j), RdR(i,j) ... mixing matrix of left- and right d-quarks 
 !  Y_u(i,j) ............. MSbar u-type Yukawa couplings
 !  RuL(j,j), RuR(i,j) ... mixing matrix of left- and right u-quarks
 !  Y_l(i,j) ............. MSbar lepton Yukawa couplings
 !  mSneut2(i) ........... sneutrino masses squared
 !  RSneut(i,j) .......... sneutrino mixing matrix
 !  mSlepton2(i) ......... slepton masses squared
 !  RSlepton(i,j) ........ slepton mixing matrix
 !  mSpm2(i) ............. mass of G+ and H+ squared
 !  RSpm(i,j) ............ mixing matrix of G+ and H+
 !  mC(i) ................ chargino masses
 !  U(i,j), V(i,j) ....... chargino mixing matrices
 !  mSup2(i) ............. u-squark masses squared
 !  RSup(i,j) ............ u-squark mixing matrix
 !  mSdown2(i) ........... d-squark masses squared
 !  RSdown(i,j) .......... d-squark mixing matrix
 !  mglu ................. gluino mass
 !  phi_g ................ phase of M_3
 !  mN(i) ................ neutralino masses
 !  N(i,j) ............... neutralino mixing matrix
 !  vevSM(i) ............. MSSM vevs (v_d, v_u)
 !  l_QCD ................ if true, then QCD corrections will be included
 !                         this is not yet completed
 ! output:
 !  BtoSNuNu ............. BR(b -> s nu nu) 
 !  Wilson coefficients at m_W (including the various contributions)
 !                              all are optional
 !  c11_o(i) .............. C_11: 1 ... total
 !                                2 ... SM contribution  
 !                                3 ... H+ contribution  
 !                                4 ... chargino contribution  
 !                                5 ... neutralino contribution  
 !                                6 ... gluino contribution  
 !  c11p_o(i) ............. C'_11: 1 ... total
 !                                 2 ... SM contribution  
 !                                 3 ... H+ contribution  
 !                                 4 ... chargino contribution  
 !                                 5 ... neutralino contribution  
 !                                 6 ... gluino contribution  
 ! written by Werner Porod, 10.10.07
 !-------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_d(3), gauge(3), mf_u(3), mW, mSpm2(2), mC(2)   &
    & , mSup2(6), mSdown2(6), mglu, mN(4), vevSM(2), mSneut2(3), mSlepton2(6)
  Complex(dp), Intent(in), Dimension(3,3) :: RdL, RdR, RuL, RuR, Y_d, Y_u, Y_l
  Complex(dp), Intent(in) :: U(2,2), V(2,2), RSpm(2,2), RSup(6,6)    &
    & , Rsneut(3,3), RSdown(6,6), phi_g, N(4,4), Rslepton(6,6)
  Logical, Intent(in) :: l_QCD
  Real(dp), Intent(out) :: BToSNuNu
  Complex(dp), Intent(out), Optional :: c11_o(3,6), c11p_o(3,6)

  Integer :: i1, i2, i3, i4
  Real(dp) :: g, gp, gs
  Complex(dp) :: cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), RSup_in(6,6), mix(6,6)  &
    & , cpl_DGSd_L(3,6), cpl_DGSd_R(3,6), cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6)
  Complex(dp) :: c11(3), c11p(3), c11_c(6), c11p_c(6)
  Complex(dp) :: kappa_q, cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6)    &
             & , cpl_NuNSn_R(3,4,3), CKM(3,3)
  Real(dp) :: mC2(2), mN2(4), a_s, sW2, tanb, e2, mG2, mW2, a_ew
  Real(dp), Parameter :: T3=-0.5_dp, ed = -1._dp / 3._dp
  Complex(dp), Parameter :: Null3(3,3) = 0._dp

 
  Iname = Iname + 1
  NameOfUnit(Iname) = "B_To_SNuNu"


  !---------------
  ! couplings
  !---------------
  gp = gauge(1)
  g = gauge(2)
  gs = gauge(3)

  CKM =  Matmul(RuL, Transpose(Conjg(RdL)) )

  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp
  If (GenerationMixing) Then
   RSup_in = RSup
  Else ! need to reorder things to be in agreement with ordering of scalar ups
   RSup_in = 0._dp
   RSup_in(1,1) = RSup(1,1)
   RSup_in(1,4) = RSup(1,2)
   RSup_in(4,1) = RSup(2,1)
   RSup_in(4,4) = RSup(2,2)
   RSup_in(2,2) = RSup(3,3)
   RSup_in(2,5) = RSup(3,4)
   RSup_in(5,2) = RSup(4,3)
   RSup_in(5,5) = RSup(4,4)
   RSup_in(3,3) = RSup(5,5)
   RSup_in(3,6) = RSup(5,6)
   RSup_in(6,3) = RSup(6,5)
   RSup_in(6,6) = RSup(6,6)
   mix = 0._dp
   mix(1:3,1:3) = rUL
   mix(4:6,4:6) = rUr
   RSup_in = Matmul(RSup_in, mix)
  End If   
  Do i1=2,3
   Do i2=1,2
    Do i3=1,6
     Call CoupCharginoSfermion3(i2, i1, i3, g, -0.5_dp, RSup_in, Y_D, Y_U &
           &, RdL, RdR, U, V, cpl_CDSu_L(i2, i1, i3), cpl_CDSu_R(i2, i1, i3) )
    End Do
   End Do
  End Do

  cpl_DGSd_L = 0._dp
  cpl_DGSd_R = 0._dp
  cpl_DNSd_L = 0._dp
  cpl_DNSd_R = 0._dp
  cpl_NuNSn_R = 0._dp
  cpl_CNuSl_L = 0._dp
  cpl_CNuSl_R = 0._dp
  If (GenerationMixing) Then
   Do i1=1,6
    Do i3=2,3
     Call CoupGluinoSquark3(gs, phi_g, i3, i1, Rsdown, RdL, RdR     &
                          & , cpl_DGSd_L(i3,i1), cpl_DGSd_R(i3,i1) )
     Do i2=1,4
      Call CoupNeutralinoSfermion3(i3, i2, i1, gp, g, T3, ed, RSdown, RdL, RdR&
                        & , Y_d, N, cpl_DNSd_L(i3,i2,i1), cpl_DNSd_R(i3,i2,i1))
     End Do
    End Do
   End Do
   Do i1=1,3
    Do i2=1,4
     Do i3=1,3  
      Call  CoupNeutralinoSneutrino3(i1, i2, i3, gp, g, N, Rsneut, id3C &
        & , cpl_NuNSn_R(i1,i2,i3))
     End Do
    End Do
    Do i2=1,2
     Do i3=1,6
     Call CoupCharginoSfermion3(i2, i1, i3, g, 0.5_dp, RSlepton, Y_l &
        & , Null3, id3C, id3C, U, V, cpl_CNuSl_L(i2,i1,i3)                  &
        & , cpl_CNuSl_R(i2,i1,i3))
     End Do
    End Do
   End Do

  Else ! no generation mixing

!    Do i1=1,4
!     Call CoupNeutralinoSneutrino(i1, gp, g, N,cpl_NuNSn_R(1,1,i1))
!    end do
!    cpl_NuNSn_R(2,2,:) = pl_NuNSn_R(1,1,:)
!    cpl_NuNSn_R(3,3,:) = pl_NuNSn_R(1,1,:)
! 
!    Do i1=1,3
!     Do i2=1,4
!      Do i3=1,2
!       i4 = (i1-1)*2 + i3  
!       Call Call CoupNeutralinoSneutrino(i1, gU1, gSU2, Nmat, coupRC)
!      End Do
!     End Do
!     Do i2=1,2
!      Do i3=1,2
!       i4 = (i1-1)*2 + i3  
!       Call CoupCharginoSfermion3(i2, i1, i1, g, 0.5_dp, RSlepton, Y_l, Null3 &
!           & , id3C, id3C, U, V, cpl_CNuSl_L(i2,i1,i4), cpl_CNuSl_R(i2,i1,i4))
!      End Do
!     End Do
!    End Do

  End If

  mW2 = mW**2
  mG2 = mglu**2
  mC2 = mC**2
  mN2 = mN**2
  tanb = vevsM(2) / vevSM(1)

  a_s = gauge(3)**2 * oo4pi
  sW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  tanb = vevSM(2) / vevSM(1)
  e2 = 4._dp * pi * gauge(2)**2 * sW2
  kappa_q = 8._dp * Sqrt2 * G_F * e2 * CKM(3,3) * Conjg( CKM(3,2) )
  kappa_q = 1._dp / kappa_q

  c11 = 0._dp
  c11p = 0._dp
  Do i1=1,3
   Call C_11_11p(3, 2, i1, mW, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2      &
     & , mN2, mSneut2, mSlepton2, mSup2, mSdown2, tanb, U, V, N, RSup, RSdown &
     & , cpl_CDSu_L, cpl_CDSu_R, cpl_CNuSl_R, cpl_NuNSn_R, cpl_DNSd_L         &
     & , cpl_DNSd_R, cpl_DGSd_L, cpl_DGSd_R, a_s, kappa_q, e2, sW2, l_QCD     &
     & , c11(i1), c11p(i1), c11_c, c11p_c)
   If (Present(c11_o)) c11_o(i1,:) = c11_c
   If (Present(c11p_o)) c11p_o(i1,:) = c11p_c
  End Do

  BToSNuNu = 0._dp
  Do i1=1,3
   BToSNuNu = BToSNuNu + Abs(C11(i1))**2 + Abs(C11p(i1))**2 &
            &          - 0.08_dp * Real(C11(i1)*Conjg(C11p(i1)),dp)
  End Do
  BToSNuNu = 5.39e-6_dp * Abs(CKM(3,2)/CKM(2,3))**2 * BToSNuNu

  Iname = Iname - 1

 End Subroutine B_To_SNuNu

 Subroutine C1_to_C6(i_loop, mu, as, mt, mW, sW2, Ci)
 !----------------------------------------------------------------------------
 ! calculates C1 to C6 in the SM up to 2-loop order
 ! Formulas are taken from C.Bobeth et al., NPB574 (2000) 291,  hep-ph/9904413
 ! input:
 !   i_loop ......... loop order, 0 ... tree, n ... n-loops
 !   mu ............. scale at which the Wilson coefficients are given
 !   as ............. alpha_s(mu)
 !   mt ............. top mass in the MSbar scheme at Q=mu 
 !   mW ............. W-boson mass
 !   sW2 ............ sin(theta_W) squared
 ! output
 !   Ci(i) .......... Wilson coefficients C_i
 ! written by Werner Porod, 02.12.2005
 !----------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_loop
  Real(dp), Intent(in) :: as, mu, mt, mW, sW2
  Complex(dp), Intent(out) :: Ci(6)

  Real(dp) :: asp, L, x, Et0_x

  Iname = Iname + 1
  NameOfUnit(Iname) = "C1_to_C6"
 
  Ci = 0._dp
  !-------------------------------------------
  ! tree level
  !-------------------------------------------
  Ci(2) = -1._dp
  If (i_loop.Eq.0) Then
   Iname = Iname - 1
   Return
  End If
  !-------------------------------------------
  ! 1-loop level
  !-------------------------------------------
  asp = oo4pi * as
  L = Log(mu**2 / mW**2)
  x = (mt/mW)**2
  Et0_x = E_t_0(x)

  Ci(1) = asp * (-15._dp - 6._dp * L)
  Ci(4) = asp * ( (7._dp - 6._dp * L) / 9._dp +  Et0_x)

  If (i_loop.Eq.1) Then
   Iname = Iname - 1
   Return
  End If

  !-------------------------------------------
  ! warning if higher loop level is required
  !-------------------------------------------
  Write(ErrCan,*) "Warning from Routine C1_to_C6"
  Write(ErrCan,*) "required loop order",i_loop,"not yet implemented"
  If (ErrorLevel.Eq.2) Call TerminateProgram()

 End Subroutine C1_to_C6


 Subroutine C_7(i, j, mW, mf_u, cpl_uWd  &
        & , mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R         &
        & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res   )
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  e m_i / (16 Pi^2) \bar{d}_L,j,a \sigma^\mu\nu d_R,i,a F_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mW, mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_uWd(3,3), cpl_CSQQp_L(:,:,:)         &
     & , cpl_CSQQp_R(:,:,:), cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)     &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_R(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(7)
  Integer :: i1, i2, n_i
  Real(dp) :: xt
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_7"

  res = 0._dp
  !----------------
  ! standard model
  !----------------
  Do i1=1,3
   xt = (mf_u(i1)/mW)**2
   res(2) = res(2) - Conjg( cpl_uWd(i1, j ) ) * cpl_uWd(i1, i)   &
       &             * xt * (2._dp * F1(xt) + 3._dp * F2(xt) )
  End Do
  res(2) = 0.25_dp * res(2) / mW**2

  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(3) = res(3) &
      & - ( Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2)  &
      &     * (e_u * F1(xt)   + F2(xt) )                          &
      &   + Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2)  &
      &     * mf_u(i2) * (e_u * F3(xt)  + F4(xt) ) / mf_d(i) )   / mSpm2(i1)
   End Do
  End Do
  res(3) = 0.25_dp * res(3) 
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(4) = res(4) &
      & + ( Conjg(cpl_CDSu_R(i1,i,i2)) * cpl_CDSu_R(i1,j,i2)    &
      &    * (e_u * F2(xt)  + F1(xt) )                          &
      &   + Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_R(i1,j,i2)  &
      &      * mC(i1) * (e_u * F4(xt)  + F3(xt) ) / mf_d(i)  )  / mSup2(i2)
   End Do
  End Do
  res(4) = 0.25_dp * res(4) 

  If (GenerationMixing) Then
   !----------
   ! gluino
   !----------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(5) = res(5) + ( Conjg(cpl_DGSd_R(i,i2)) * cpl_DGSd_R(j,i2) * F2(xt) &
      &               + Conjg(cpl_DGSd_L(i,i2)) * cpl_DGSd_R(j,i2)          &
      &                 * mGlu * F4(xt) / mf_d(i) )   / mSdown2(i2)
   End Do
   res(5) = e_d * 4._dp * res(5) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(6) = res(6) &
      & + ( Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2) * F2(xt)   &
      &   + Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2)            &
      &         * mN(i1) * F4(xt) / mf_d(i) )   / mSdown2(i2)
    End Do
   End Do
   res(6) = 0.25_dp * e_d * res(6) 
   !----------
   ! Higgs
   !----------
   Do i1=1,3
    Do i2=1,2
     xt = mf_d(i1)**2/mS02(i2)
     res(7) = res(7) &
       & + ( Conjg(cpl_DDS0_R(i,i1,i2)) * cpl_DDS0_R(i1,j,i2) * F2(xt) &
       &   + Conjg( cpl_DDS0_L(i,i1,i2) ) * cpl_DDS0_R(i1,j,i2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mS02(i2)
    End Do
     xt = mf_d(i1)**2/mP02(2)
     res(7) = res(7) &
       & + ( Conjg(cpl_DDP0_R(i1,i,2)) * cpl_DDP0_R(i1,j,2) * F2(xt) &
       &   + Conjg( cpl_DDP0_L(i1,i,2) ) * cpl_DDP0_R(i1,j,2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mP02(2)
   End Do
   res(7) = 0.25_dp * e_d * res(7) 
  End If

  res(1) = Sum( res(2:7) )
  
  Iname = Iname - 1

 End Subroutine C_7

 Subroutine C_7p(i, j, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R &
         & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res  )
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  e m_i / (16 Pi^2) \bar{d}_R,j,a \sigma^\mu\nu d_L,i,a F_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_CSQQp_L(:,:,:), cpl_CSQQp_R(:,:,:)   &
      & , cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)                        &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_r(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(6)

  Integer :: i1, i2, n_i
  Real(dp) :: xt

  Iname = Iname + 1
  NameOfUnit(Iname) = "C_7p"

  !----------------
  ! standard model
  !----------------
  res = 0._dp
  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(2) = res(2) &
      & - ( Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2)    &
      &          * (e_u * F1(xt)  + F2(xt) )                        &
      &   + Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2)    &
      &          * mf_u(i2) * (e_u * F3(xt)  + F4(xt) ) / mf_d(i) ) &
      &  / mSpm2(i1)
   End Do
  End Do
  res(2) = 0.25_dp * res(2)
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(3) = res(3) &
      & + ( Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2)    &
      &                * (e_u * F2(xt)  + F1(xt) )                &
      &   + Conjg( cpl_CDSu_R(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2)    &
      &         * mC(i1) * (e_u * F4(xt)  + F3(xt) ) / mf_d(i))   &
      &  / mSup2(i2)
   End Do
  End Do
  res(3) = 0.25_dp * res(3)

  If (GenerationMixing) Then
   !----------
   ! gluino
   !----------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(4) = res(4) &
      & + ( Conjg( cpl_DGSd_L(i,i2) ) * cpl_DGSd_L(j,i2) * F2(xt)  &
      &   + Conjg( cpl_DGSd_R(i,i2) ) * cpl_DGSd_L(j,i2)           &
      &     * mGlu * F4(xt) / mf_d(i) )    / mSdown2(i2)
   End Do
   res(4) = e_d * 4._dp * res(4) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(5) = res(5) &
      & + ( Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2) * F2(xt)  &
      &   + Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2)           &
      &         * mN(i1) * F4(xt) / mf_d(i)  )     / mSdown2(i2)
    End Do
   End Do
   res(5) = 0.25_dp * e_d * res(5) 
   !----------
   ! Higgs
   !----------
   Do i1=1,3
    Do i2=1,2
     xt = mf_d(i1)**2/mS02(i2)
     res(6) = res(6) &
       & + ( Conjg(cpl_DDS0_L(i1,i,i2)) * cpl_DDS0_L(i1,j,i2) * F2(xt) &
       &   + Conjg( cpl_DDS0_R(i1,i,i2) ) * cpl_DDS0_L(i1,j,i2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mS02(i2)
    End Do
     xt = mf_d(i1)**2/mP02(2)
     res(6) = res(6) &
       & + ( Conjg(cpl_DDP0_L(i,i1,2)) * cpl_DDP0_L(i1,j,2) * F2(xt) &
       &   + Conjg( cpl_DDP0_R(i,i1,2) ) * cpl_DDP0_L(i1,j,2)  &
       &      * mf_d(i1) *  F4(xt) / mf_d(i)  )  / mP02(2)
   End Do
   res(6) = 0.25_dp * e_d * res(6) 
  End If

  res(1) = Sum( res(2:6) )
  
  Iname = Iname - 1

 End Subroutine C_7p

 Subroutine C_8(i, j, mW, mf_u, cpl_uWd  &
        & , mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R   &
        & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res)
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  g_s m_i / (16 Pi^2) \bar{d}_L,j,a \sigma^\mu\nu d_R,i,b t_c^ab G^c_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mW, mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_uWd(3,3), cpl_CSQQp_L(:,:,:)         &
     & , cpl_CSQQp_R(:,:,:), cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)     &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_r(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(7)

  Integer :: i1, i2, n_i
  Real(dp) :: xt
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_8"

  res = 0._dp
  !----------------------------
  ! standard model
  !----------------------------
  Do i1=1,3
   xt = (mf_u(i1)/mW)**2
   res(2) = res(2) - Conjg( cpl_uWd(i1, j ) ) * cpl_uWd(i1, i) * xt * F1(xt)
  End Do
  res(2) = 0.75_dp * res(2) / mW**2
  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(3) = res(3) &
      & - ( Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2) * F1(xt)  &
      &   + Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_R(i1,j,i2)  &
      &                *  mf_u(i2) * F3(xt) / mf_d(i) )      / mSpm2(i1)
   End Do
  End Do
  res(3) = 0.25_dp * res(3)
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(4) = res(4) &
      & + ( Conjg( cpl_CDSu_R(i1,i,i2) ) * cpl_CDSu_R(i1,j,i2) *  F2(xt)    &
      &   + Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_R(i1,j,i2)              &
      &                          * mC(i1) * F4(xt) / mf_d(i))    / mSup2(i2)
   End Do
  End Do
  res(4) = 0.25_dp * res(4)

  If (GenerationMixing) Then
   !----------
   ! gluino
   !----------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(5) = res(5) &
      & + ( Conjg( cpl_DGSd_R(i,i2) ) * cpl_DGSd_R(j,i2)       &
      &                * (4._dp * F1(xt)  + F1(xt) )           &
      &      + Conjg( cpl_DGSd_L(i,i2) ) * cpl_DGSd_R(j,i2)    &
      &         * mGlu * (4._dp * F3(xt) + F3(xt) ) / mf_d(i)) &
      &  / mSdown2(i2)
   End Do
   res(5) = - 0.25_dp * res(5) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(6) = res(6) &
      & + ( Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2) * F2(xt)  &
      &   + Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_R(j,i1,i2)           &
      &         * mN(i1) * F4(xt) / mf_d(i) )   / mSdown2(i2)
    End Do
   End Do
   res(6) = 0.25_dp * res(6)
   !----------
   ! Higgs
   !----------
   Do i1=1,3
    Do i2=1,2
     xt = mf_d(i1)**2/mS02(i2)
     res(7) = res(7) &
       & - ( Conjg( cpl_DDS0_R(i,i1,i2) ) * cpl_DDS0_R(i1,j,i2) * F1(xt)  &
       &   + Conjg( cpl_DDS0_L(i,i1,i2) ) * cpl_DDS0_R(i1,j,i2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mS02(i2)
    End Do
     xt = mf_d(i1)**2/mP02(2)
     res(7) = res(7) &
       & - ( Conjg( cpl_DDP0_R(i,i1,2) ) * cpl_DDP0_R(i1,j,2) * F1(xt)  &
       &   + Conjg( cpl_DDP0_L(i,i1,2) ) * cpl_DDP0_R(i1,j,2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mP02(2)
   End Do
   res(7) = 0.25_dp * res(7)
  End If

  res(1) = Sum(res(2:7))

  Iname = Iname - 1

 End Subroutine C_8

 Subroutine C_8p(i, j, mf_u, mf_d, mSpm2, cpl_CSQQp_L, cpl_CSQQp_R &
         & , mC, mSup2, cpl_CDSu_L, cpl_CDSu_R, mGlu, mSdown2, cpl_DGSd_L    &
        & , cpl_DGSd_R, mN,  cpl_DNSd_L, cpl_DNSd_R, mS02, cpl_DDS0_L       &
        & , cpl_DDS0_r, mP02, cpl_DDP0_L, cpl_DDP0_R, res )
 !-------------------------------------------------------------------------
 ! coefficient of the operator:
 !  g_s m_i / (16 Pi^2) \bar{d}_R,j,a \sigma^\mu\nu d_L,i,b t_c^ab G^c_\mu\nu
 ! a,b, are SU(3) indices
 ! written by Werner Porod, 11.11.02
 !-------------------------------------------------------------------------
 
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mf_u(3), mSpm2(:), mC(:), mSup2(:), mf_d(3)  &
     & , mSdown2(:), mGlu, mN(:), mS02(:), mP02(:)
  Complex(dp), Intent(in) :: cpl_CSQQp_L(:,:,:), cpl_CSQQp_R(:,:,:)   &
      & , cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)                        &
     & , cpl_DGSd_L(:,:),  cpl_DGSd_R(:,:), cpl_DNSd_L(:,:,:)         &
     & , cpl_DNSd_R(:,:,:), cpl_DDS0_L(:,:,:), cpl_DDS0_r(:,:,:)      &
     & , cpl_DDP0_L(:,:,:), cpl_DDP0_R(:,:,:)
  Complex(dp), Intent(out) :: res(6)

  Integer :: i1, i2, n_i
  Real(dp) :: xt

  Iname = Iname + 1
  NameOfUnit(Iname) = "C_8p"

  !----------------------------
  ! standard model, tree level
  !----------------------------
  res = 0._dp
  !---------------------
  ! charged Higgs boson
  !---------------------
  n_i = Size( mSpm2 )
  Do i1=2,n_i
   Do i2=1,3
    xt = mf_u(i2)**2/mSpm2(i1)
    res(2) = res(2) &
      & - ( Conjg( cpl_CSQQp_L(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2) * F1(xt)   &
      &   + Conjg( cpl_CSQQp_R(i1,i,i2) ) * cpl_CSQQp_L(i1,j,i2)            &
      &                * mf_u(i2) * F3(xt) / mf_d(i) ) / mSpm2(i1)
   End Do
  End Do
  res(2) = 0.25_dp * res(2)
  !----------
  ! chargino
  !----------
  n_i = Size( mC )
  Do i1=1,n_i
   Do i2=1,6
    xt = mC(i1)**2/mSup2(i2)
    res(3) = res(3) &
      & + ( Conjg( cpl_CDSu_L(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2) * F2(xt)    &
      &   + Conjg( cpl_CDSu_R(i1,i,i2) ) * cpl_CDSu_L(i1,j,i2)             &
      &                   * mC(i1) * F4(xt) / mf_d(i)    )     / mSup2(i2)
   End Do
  End Do
  res(3) = 0.25_dp * res(3)

  If (GenerationMixing) Then
   !--------
   ! gluino
   !--------
   Do i2=1,6
    xt = mGlu**2/mSdown2(i2)
    res(4) = res(4) &
      & + ( Conjg( cpl_DGSd_L(i,i2) ) * cpl_DGSd_L(j,i2)         &
      &                * (4._dp * F1(xt)  + F1(xt) )             &
      &   + Conjg( cpl_DGSd_R(i,i2) ) * cpl_DGSd_L(j,i2)         &
      &         * mGlu * (4._dp * F3(xt) + F3(xt) ) / mf_d(i))   &
      &  / mSdown2(i2)
   End Do
   res(4) = - 0.25_dp * res(4) / 3._dp
   !----------
   ! neutralino
   !----------
   n_i = Size( mN )
   Do i1=1,n_i
    Do i2=1,6
     xt = mN(i1)**2/mSdown2(i2)
     res(5) = res(5) &
      & + ( Conjg( cpl_DNSd_L(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2) * F2(xt)  &
      &   + Conjg( cpl_DNSd_R(i,i1,i2) ) * cpl_DNSd_L(j,i1,i2)           &
      &                      * mN(i1) * F4(xt) / mf_d(i)     )   / mSdown2(i2)
    End Do
   End Do
   res(5) = 0.25_dp * res(5)
   !----------
   ! Higgs
   !----------
   Do i1=1,3
    Do i2=1,2
     xt = mf_d(i1)**2/mS02(i2)
     res(6) = res(6) &
       & - ( Conjg( cpl_DDS0_L(i,i1,i2) ) * cpl_DDS0_L(i1,j,i2) * F1(xt)  &
       &   + Conjg( cpl_DDS0_R(i,i1,i2) ) * cpl_DDS0_L(i1,j,i2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mS02(i2)
    End Do
     xt = mf_d(i1)**2/mP02(2)
     res(6) = res(6) &
       & - ( Conjg( cpl_DDP0_L(i,i1,2) ) * cpl_DDP0_L(i1,j,2) * F1(xt)  &
       &   + Conjg( cpl_DDP0_R(i,i1,2) ) * cpl_DDP0_L(i1,j,2)  &
       &                *  mf_d(i1) * F3(xt) / mf_d(i) )      / mP02(2)
   End Do
   res(6) = 0.25_dp * res(6)
  End If

  res(1) = Sum( res(2:6) )
  
  Iname = Iname - 1

 End Subroutine C_8p

 Subroutine C_9(i, j, k, mu, mf_d, mf_u, mf_l, mW2, mHp2, mGlu2, mC2        &
       & , mN2, mSnu2, mSlept2, mSqu2, mSqd2, tanb, U, V, N, RSup, RSdown   &
       & , c_CDSu_L, c_CDSu_R, c_CLSn_L, c_CLSn_R, c_LNSl_L, c_LNSl_R       &
       & , c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R                           &
       & , a_s, kappa_q, e2, sW2, l_QCD, c9, c9p, c9_c, c9p_c)
 !---------------------------------------------------------------------------
 ! in this routine the coefficients c_A and c_A' of the effective operator
 !
 !             2 G_F    alpha
 !    H_eff = -------------------- V_ti V_tj^* (c_9 O_A + c_9' O_A')
 !             Sqrt(2) 2 Pi s^2_W
 !
 ! with O_A = (\bar(q) \gamma^\mu P_L b) (\bar(l) \gamma_\mu l)
 !      O_A' = (\bar(q) \gamma^\mu P_R b) (\bar(l) \gamma_\mu l)
 ! for the transitions b -> q l+ l- with q=s,d;  s -> d mu+ mu-;
 ! B_q -> l+ l- and K_L -> mu+ mu- are calculated.
 ! The formulas include QCD corrections and are based on
 ! C.bobeth, A.J.Buras, F.Krueger, J.Urban, NPB 630, (2002) 87
 ! written by Werner Porod, 14.12.2004
 !---------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j              ! indices of the d-type quarks, i>j
  Integer, Intent(in) :: k                 ! generation index of the lepton l
  Real(dp), Intent(in) :: mu               ! renormalization scale
  Real(dp), Intent(in) :: mf_d(3), mf_u(3) ! running d- and u-quark masses at mu
  Real(dp), Intent(in) :: mf_l(3)          ! lepton masses
  Real(dp), Intent(in) :: mW2              ! W-boson mass squared
  Real(dp), Intent(in) :: mHp2             ! mass of charged Higgs boson squared
  Real(dp), Intent(in) :: mglu2            ! gluino mass squared
  Real(dp), Intent(in) :: mC2(2)           ! chargino masses squared
  Real(dp), Intent(in) :: mN2(4)           ! neutralino masses squared
  Real(dp), Intent(in) :: mSnu2(3)         ! sneutrino masses squared
  Real(dp), Intent(in) :: mSlept2(6)       ! slepton masses squared
  Real(dp), Intent(in) :: mSqu2(6)         ! u-squark masses squared
  Real(dp), Intent(in) :: mSqd2(6)         ! d-squark masses squared
  Real(dp), Intent(in) :: a_s              ! alpha_s(mu)
  Real(dp), Intent(in) :: tanb             ! tan(beta) = v_2/v_1
  ! coupling chargino - d-quark - u-squark
  Complex(dp), Intent(in) ::  c_CDSu_L(2,3,6), c_CDSu_R(2,3,6)
  ! coupling neutralino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DNSd_L(3,4,6), c_DNSd_R(3,4,6)
  ! coupling gluino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DGSd_L(3,6), c_DGSd_R(3,6)
  ! coupling chargino - lepton - sneutrino
  Complex(dp), Intent(in) ::  c_CLSn_L(2,3,3), c_CLSn_R(2,3,3)
  ! coupling neutralino - lepton - slepton
  Complex(dp), Intent(in) ::  c_LNSl_L(3,4,6), c_LNSl_R(3,4,6)
  Complex(dp), Intent(in) :: kappa_Q    ! = (8 sqrt(2) G_F e^2 V_ti V^*_tj)^-1
  Real(dp), Intent(in) :: e2            ! = e^2 = g^2 sin^2(theta_W)
  Real(dp), Intent(in) :: sW2           ! sin^2(theta_W)
  Logical, Intent(in) :: l_QCd          ! if .true. QCD corrections are included
  ! chargino mixing matrices
  Complex(dp), Intent(in) :: U(2,2), V(2,2)
  ! neutralino mixing matrix
  Complex(dp), Intent(in) :: N(4,4)
  ! squark mixing matrices
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6)

  Complex(dp), Intent(out) :: C9, C9p ! Wilson coefficients
  Complex(dp), Intent(out), Optional :: C9_c(6), C9p_c(6) ! Wilson coefficients
                                                 ! + individual contributions

  Integer :: i1, i2, i3, i4
  Real(dp) :: x, y, z, x_ij(2,2), y_ai(6,2), v_fi(3,2), r_a(6), s_ai(6,4) &
      & , n_bi(6,4), n_ij(4,4), Lt, Lua(6), a_sp, x_in
  Complex(dp) :: ZNN(4,4), B_L(4), B_R(4), C_L(5), C_R(5), D_L(5), D_R(5) &
      &  , fact, fac(2)

  !----------------------------------
  ! Initialisation
  !----------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_9"

  c9 = 0._dp
  c9p = 0._dp
  B_L = 0._dp
  B_R = 0._dp
  C_L = 0._dp
  C_R = 0._dp
  D_L = 0._dp
  D_R = 0._dp

  x = mf_u(3)**2 / mW2
  y = mf_u(3)**2 / mHp2
  z = mHp2 / mW2
  x_ij(1,1) = 1._dp
  x_ij(2,2) = 1._dp
  x_ij(1,2) = mC2(1) / mC2(2)
  x_ij(2,1) = 1._dp / x_ij(1,2)
  r_a = mSqd2/mGlu2
  Do i1=1,2
   Do i2=1,6
    y_ai(i2,i1) = mSqu2(i2) / mC2(i1)
   End Do
   Do i2=1,3
    v_fi(i2,i1) = mSnu2(i2) / mC2(i1)
   End Do
  End Do
  Do i1=1,4
   n_ij(i1,i1) = 1._dp
   ZNN(i1,i1) = N(i1,4) * Conjg( N(i1,4) ) - N(i1,3) * Conjg( N(i1,3) )
   Do i2=i1+1,4
    n_ij(i1,i2) = mN2(i1)/mN2(i2)
    n_ij(i2,i1) = 1._dp / n_ij(i1,i2)
    ZNN(i1,i2) = N(i1,4) * Conjg( N(i2,4) ) - N(i1,3) * Conjg( N(i2,3) )
    ZNN(i2,i1) = Conjg( ZNN(i1,i2) ) 
   End Do
   Do i2=1,6
    s_ai(i2,i1) = mSqd2(i2) / mN2(i1)
    n_bi(i2,i1) = mSlept2(i2) / mN2(i1)
   End Do
  End Do

  Lt = Log(mu**2 / mf_u(3)**2)
  Lua = Log(mu**2 / mSqu2)

  a_sp = oo4Pi * a_s
  !-------------------------------------------
  ! SM contribution
  !-------------------------------------------
  If (l_QCD) Then
   C_L(1) = 0.25_dp * (f_1_0(x) + a_sp * (f_1_1(x)+ 8._dp* x * Df_1_0(x) * Lt)) 
   B_L(1) = 0.25_dp * (f_2_0(x) + a_sp * (f_10_1(x) + 8._dp* x *Df_2_0(x) * Lt))
  Else
   C_L(1) = 0.25_dp * f_1_0(x)
   B_L(1) = 0.25_dp * f_2_0(x)
  End If 
  D_L(1) = h_0_1(x)

   C9 = ((1._dp - 4._dp * sW2) * C_L(1) - B_L(1)) / sW2 - D_L(1)
   C9p = ((1._dp - 4._dp * sW2) * C_R(1) - B_R(1)) / sW2 - D_R(1)
  !------------------------------------------
  ! H+ contribution
  !------------------------------------------
  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_2_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  Else
   fact = f_2_0(y) 
  End If

  C_L(2) = - 0.125_dp * x * fact / tanb**2 ! m^2_H/m^2_W * y = x
  C_R(2) = 0.125_dp * mf_d(i) * mf_d(j) * tanb**2 * fact / mW2

  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_7_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  End If

  B_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * fact * mf_l(k)**2 &
         &       / (16._dp * mHp2 * mW2)

  D_L(2) = y * (2._dp/3._dp * h_0_5(y) - h_0_6(y)) / tanb**2
  D_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * mf_l(k)**2  / (16._dp * mHp2**2) &
         & * y * (2._dp/3._dp * h_0_5(y) - h_0_6(y))

  C9 = (1._dp - 4._dp * sW2) * C_L(2) / sW2 - D_L(2)
  C9p = ((1._dp - 4._dp * sW2) * C_R(2) - B_R(2)) / sW2 - D_R(2)

  !------------------------------------------------------------------
  ! Chargino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2))                                      &
         &   * ( f_3_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &              + y_ai(i3,i2) * Dyf_3_0(x_ij(i1,i2),y_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &            + y_ai(i3,i2) * Dyf_4_0(x_ij(i1,i2),y_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2)) * f_3_0(x_ij(i1,i2),y_ai(i3,i2))
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2))
     End If
     C_L(3) = C_L(3) - c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i2,i,i3))     &
         &             * ( fac(1) * U(i1,1) * Conjg(U(i2,1))            &
         &               - fac(2) * V(i2,1) * Conjg(V(i1,1)) )
     C_R(3) = C_R(3) - c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i2,i,i3))     &
          &             * ( fac(1) * V(i2,1) * Conjg(V(i1,1))           &
          &               - fac(2) * U(i1,1) * Conjg(U(i2,1)) )
     Do i4=1,3
      If (l_QCD) Then
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      Else
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      End If
      B_L(3) = B_L(3) + Conjg(c_CDSu_R(i2,i,i3))* c_CDSu_R(i1,j,i3) &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_R(i1,k,i4))          &
          &                    * c_CLSn_R(i2,k,i4)                          &
          &           - fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_L(i2,k,i4)       &
          &                    * Conjg(c_CLSn_L(i1,k,i4))    ) / mC2(i2)
      B_R(3) = B_R(3) - Conjg(c_CDSu_L(i2,i,i3))*c_CDSu_L(i1,j,i3) &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_L(i1,k,i4))           &
          &                    * c_CLSn_L(i2,k,i4)                           &
          &           - fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_R(i2,k,i4)        &
          &                    * Conjg(c_CLSn_R(i1,k,i4))    ) / mC2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    x_in = 1._dp / y_ai(i2,i1)
    fac(1) = (2._dp/3._dp * h_0_6(x_in) - h_0_5(x_in)) / mSqu2(i2)
    D_L(3) = D_L(3) + fac(1) * Conjg(c_CDSu_R(i1,i,i2))* c_CDSu_R(i1,j,i2)
    D_R(3) = D_R(3) + fac(1) * Conjg(c_CDSu_L(i1,i,i2))* c_CDSu_L(i1,j,i2)
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(y_ai(i2,i1),y_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(y_ai(i2,i1),y_ai(i3,i1))                         &
           &          + 4._dp * (y_ai(i2,i1)*DxF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    +y_ai(i3,i1)*DyF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( y_ai(i2,i1), y_ai(i3,i1))
     End If
     C_L(3) = C_L(3) - fac(1) * ( RSup(i2,1)*Conjg(RSup(i3,1))      &
         &                      + RSup(i2,2)*Conjg(RSup(i3,2))      &
         &                      + RSup(i2,3)*Conjg(RSup(i3,3))  )   &
         &               * c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i1,i,i2))
     C_R(3) = C_R(3) + fac(1) * ( RSup(i3,4)*Conjg(RSup(i2,4))      &
         &                      + RSup(i3,5)*Conjg(RSup(i2,5))      &
         &                      + RSup(i3,6)*Conjg(RSup(i2,6))  )   &
         &                    * c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i1,i,i2))
    End Do
   End Do
  End Do
  B_L(3) = B_L(3) * Kappa_Q * sW2
  B_R(3) = B_R(3) * Kappa_Q * sW2
  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(3) = C_L(3) * fact
  C_R(3) = C_R(3) * fact
  fact = 16._dp * fact * mW2
  D_L(3) = D_L(3) * fact
  D_R(3) = D_R(3) * fact

   C9 = ((1._dp - 4._dp * sW2) * C_L(3) - B_L(3)) / sW2 - D_L(3)
   C9p = ((1._dp - 4._dp * sW2) * C_R(3) - B_R(3)) / sW2 - D_R(3)

  !------------------------------------------------------------------
  ! neutralino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,4
   Do i2=1,4
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2))                                      &
         &   * ( f_3_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &              + s_ai(i3,i2) * Dyf_3_0(n_ij(i1,i2),s_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &            + s_ai(i3,i2) * Dyf_4_0(n_ij(i1,i2),s_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2)) * f_3_0(n_ij(i1,i2),s_ai(i3,i2))
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2))
     End If

     C_L(4) = C_L(4) - (fac(1) * ZNN(i1,i2) + fac(2) * ZNN(i2,i1) )&
         &            * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i2,i3))
     C_R(4) = C_R(4) - (fac(1) * ZNN(i1,i2) + fac(2) * ZNN(i2,i1) )&
         &            * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i2,i3)) 

     Do i4=1,6
      If (l_QCD) Then
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      Else
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      End If
      B_L(4) = B_L(4) - Conjg(c_DNSd_R(i,i2,i3))* c_DNSd_R(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * (- Conjg(c_LNSl_R(k,i1,i4)) * c_LNSl_R(k,i2,i4)   &
          &                 + Conjg(c_LNSl_L(k,i2,i4)) * c_LNSl_L(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_L(k,i2,i4) * Conjg(c_LNSl_L(k,i1,i4))   &
          &                 - c_LNSl_R(k,i1,i4) * Conjg(c_LNSl_R(k,i2,i4)) ) &
          &           ) / mN2(i2)
      B_R(4) = B_R(4) + Conjg(c_DNSd_L(i,i2,i3))*c_DNSd_L(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * ( Conjg(c_LNSl_L(k,i1,i4)) * c_LNSl_L(k,i2,i4)   &
          &                 + Conjg(c_LNSl_R(k,i2,i4)) * c_LNSl_R(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_R(k,i2,i4) * Conjg(c_LNSl_R(k,i1,i4))   &
          &                 - c_LNSl_L(k,i1,i4) * Conjg(c_LNSl_L(k,i2,i4)) ) &
          &           ) / mN2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    x_in = 1._dp/s_ai(i2,i1)
    fac(1) = h_0_6(x_in) / mSqd2(i2)
    D_L(4) = D_L(4) + fac(1) * Conjg(c_DNSd_R(i,i1,i2))* c_DNSd_R(j,i1,i2)
    D_R(4) = D_R(4) + fac(1) * Conjg(c_DNSd_L(i,i1,i2))* c_DNSd_L(j,i1,i2)
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(s_ai(i2,i1),s_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(s_ai(i2,i1),s_ai(i3,i1))                         &
           &          + 4._dp * (s_ai(i2,i1)*DxF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    +s_ai(i3,i1)*DyF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( s_ai(i2,i1), s_ai(i3,i1))
     End If
     C_L(4) = C_L(4) - fac(1) * ( RSdown(i3,4)*Conjg(RSdown(i2,4))      &
         &                      + RSdown(i3,5)*Conjg(RSdown(i2,5))      &
         &                      + RSdown(i3,6)*Conjg(RSdown(i2,6))  )   &
         &                    * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i1,i2))
     C_R(4) = C_R(4) + fac(1) * ( RSdown(i3,1)*Conjg(RSdown(i2,1))      &
         &                      + RSdown(i3,2)*Conjg(RSdown(i2,2))      &
         &                      + RSdown(i3,3)*Conjg(RSdown(i2,3))  )   &
         &                    * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i1,i2))
    End Do
   End Do
    
  End Do
  B_L(4) = B_L(4) * Kappa_Q * sW2
  B_R(4) = B_R(4) * Kappa_Q * sW2
  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(4) = C_L(4) * fact
  C_R(4) = C_R(4) * fact
  fact = - 4._dp * kappa_Q * e2 / 3._dp
  D_L(4) = D_L(4) * fact
  D_R(4) = D_R(4) * fact

  !------------------------------------------------------------------
  ! gluino contributions
  !------------------------------------------------------------------
   Do i2=1,6
    x_in = 1._dp/r_a(i2)
    fac(1) = h_0_6(x_in) / mSqd2(i2)
    D_L(5) = D_L(5) + fac(1) * Conjg(c_DGSd_R(i,i2))* c_DGSd_R(j,i2)
    D_R(5) = D_R(5) + fac(1) * Conjg(c_DGSd_L(i,i2))* c_DGSd_L(j,i2)
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     Else
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     End If
     C_L(5) = C_L(5) - fac(1) * ( RSdown(i3,4)*Conjg(RSdown(i2,4))      &
         &                      + RSdown(i3,5)*Conjg(RSdown(i2,5))      &
         &                      + RSdown(i3,6)*Conjg(RSdown(i2,6))  )   &
         &                    * c_DGSd_R(j,i3) * Conjg(c_DGSd_R(i,i2))
     C_R(5) = C_R(5) + fac(1) * ( RSdown(i3,1)*Conjg(RSdown(i2,1))      &
         &                      + RSdown(i3,2)*Conjg(RSdown(i2,2))      &
         &                      + RSdown(i3,3)*Conjg(RSdown(i2,3))  )   &
         &                    * c_DGSd_L(j,i3) * Conjg(c_DGSd_L(i,i2))
    End Do
   End Do
  fact = 4._dp * kappa_Q * e2  / (3._dp * mW2)
  C_L(5) = C_L(5) * fact
  C_R(5) = C_R(5) * fact
  fact = - 64._dp * kappa_Q * e2 / 9._dp
  D_L(5) = D_L(5) * fact
  D_R(5) = D_R(5) * fact

  C9 = ((1._dp - 4._dp * sW2) *Sum(C_L) - Sum(B_L)) / sW2 - Sum(D_L)
  C9p = ((1._dp - 4._dp * sW2) * Sum(C_R) - Sum(B_R)) / sW2 - Sum(D_R)
  !---------------------
  ! c-quark contribution
  !---------------------

  C9 = c9 + 38._dp/27._dp - 4._dp * Log(mu**2/mW2)/ 9._dp

  If (Present(c9_c)) Then
   c9_c(1) = c9
   c9_c(2:6) = ((1._dp - 4._dp * sW2) * C_L ) / sW2 - D_L
   c9_c(2) = c9_c(2) + 38._dp/27._dp - 4._dp * Log(mu**2/mW2)/ 9._dp
   c9_c(2:5) =c9_c(2:5) - B_L / sW2
  End If
  If (Present(c9p_c)) Then
   c9p_c(1) = c9p
   c9p_c(2:6) = ((1._dp - 4._dp * sW2) * C_R ) / sW2 - D_R
   c9p_c(2:5) =c9p_c(2:5) - B_R / sW2
  End If

  Iname = Iname - 1

 End Subroutine C_9


 Subroutine C_QdQdLL_AAp(i, j, k, mu, mf_d, mf_u, mf_l, mW2, mHp2, mGlu2, mC2 &
       & , mN2, mSnu2, mSlept2, mSqu2, mSqd2, tanb, U, V, N, RSup, RSdown     &
       & , c_CDSu_L, c_CDSu_R, c_CLSn_L, c_CLSn_R, c_LNSl_L, c_LNSl_R  &
       & , c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R  &
       & , a_s, kappa_q, e2, sW2, l_QCD, c_A, c_Ap, c10_c, c10p_c)
 !---------------------------------------------------------------------------
 ! in this routine the coefficients c_A and c_A' of the effective operator
 !
 !             2 G_F    alpha
 !    H_eff = -------------------- V_ti V_tj^* (c_A O_A + c_A' O_A')
 !             Sqrt(2) 2 Pi s^2_W
 !
 ! with O_A = (\bar(q) \gamma^\mu P_L b) (\bar(l) \gamma_\mu gamma_5 l)
 !      O_A' = (\bar(q) \gamma^\mu P_R b) (\bar(l) \gamma_\mu gamma_5 l)
 ! for the transitions b -> q l+ l- with q=s,d;  s -> d mu+ mu-;
 ! B_q -> l+ l- and K_L -> mu+ mu- are calculated.
 ! The formulas include QCD corrections and are based on
 ! C.Bobeth, A.J.Buras, F.Krueger, J.Urban, NPB 630, (2002) 87
 ! written by Werner Porod, 14.12.2004
 !---------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j              ! indices of the d-type quarks, i>j
  Integer, Intent(in) :: k                 ! generation index of the lepton l
  Real(dp), Intent(in) :: mu               ! renormalization scale
  Real(dp), Intent(in) :: mf_d(3), mf_u(3) ! running d- and u-quark masses at mu
  Real(dp), Intent(in) :: mf_l(3)          ! lepton masses
  Real(dp), Intent(in) :: mW2              ! W-boson mass squared
  Real(dp), Intent(in) :: mHp2             ! mass of charged Higgs boson squared
  Real(dp), Intent(in) :: mglu2            ! gluino mass squared
  Real(dp), Intent(in) :: mC2(2)           ! chargino masses squared
  Real(dp), Intent(in) :: mN2(4)           ! neutralino masses squared
  Real(dp), Intent(in) :: mSnu2(3)         ! sneutrino masses squared
  Real(dp), Intent(in) :: mSlept2(6)       ! slepton masses squared
  Real(dp), Intent(in) :: mSqu2(6)         ! u-squark masses squared
  Real(dp), Intent(in) :: mSqd2(6)         ! d-squark masses squared
  Real(dp), Intent(in) :: a_s              ! alpha_s(mu)
  Real(dp), Intent(in) :: tanb             ! tan(beta) = v_2/v_1
  ! coupling chargino - d-quark - u-squark
  Complex(dp), Intent(in) ::  c_CDSu_L(2,3,6), c_CDSu_R(2,3,6)
  ! coupling neutralino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DNSd_L(3,4,6), c_DNSd_R(3,4,6)
  ! coupling gluino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DGSd_L(3,6), c_DGSd_R(3,6)
  ! coupling chargino - lepton - sneutrino
  Complex(dp), Intent(in) ::  c_CLSn_L(2,3,3), c_CLSn_R(2,3,3)
  ! coupling neutralino - lepton - slepton
  Complex(dp), Intent(in) ::  c_LNSl_L(3,4,6), c_LNSl_R(3,4,6)
  Complex(dp), Intent(in) :: kappa_Q    ! = (8 sqrt(2) G_F e^2 V_ti V^*_tj)^-1
  Real(dp), Intent(in) :: e2            ! = e^2 = g^2 sin^2(theta_W)
  Real(dp), Intent(in) :: sW2           ! sin^2(theta_W)
  Logical, Intent(in) :: l_QCd          ! if .true. QCD corrections are included
  ! chargino mixing matrices
  Complex(dp), Intent(in) :: U(2,2), V(2,2)
  ! neutralino mixing matrix
  Complex(dp), Intent(in) :: N(4,4)
  ! squark mixing matrices
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6)
!  Real(dp), intent(in) ::   !
  Complex(dp), Intent(out) :: C_A, C_Ap ! Wilson coefficients
  Complex(dp), Intent(out), Optional :: C10_c(6),C10p_c(6) ! Wilson coefficients
                                                 ! + individual contributions

  Integer :: i1, i2, i3, i4
  Real(dp) :: x, y, z, x_ij(2,2), y_ai(6,2), v_fi(3,2), r_a(6), s_ai(6,4) &
      & , n_bi(6,4), n_ij(4,4), Lt, Lua(6), a_sp
  Complex(dp) :: ZNN(4,4), fact, fac(2), c_L(5), c_R(5), B_L(4), B_R(4)

  !----------------------------------
  ! Initialisation
  !----------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_QdQdLL_AAp"

  C_A = 0._dp
  C_Ap = 0._dp
  B_L = 0._dp
  B_R = 0._dp
  C_L = 0._dp
  C_R = 0._dp

  x = mf_u(3)**2 / mW2
  y = mf_u(3)**2 / mHp2
  z = mHp2 / mW2
  x_ij(1,1) = 1._dp
  x_ij(2,2) = 1._dp
  x_ij(1,2) = mC2(1) / mC2(2)
  x_ij(2,1) = 1._dp / x_ij(1,2)
  r_a = mSqd2/mGlu2
  Do i1=1,2
   Do i2=1,6
    y_ai(i2,i1) = mSqu2(i2) / mC2(i1)
   End Do
   Do i2=1,3
    v_fi(i2,i1) = mSnu2(i2) / mC2(i1)
   End Do
  End Do
  Do i1=1,4
   n_ij(i1,i1) = 1._dp
   ZNN(i1,i1) = N(i1,4) * Conjg( N(i1,4) ) - N(i1,3) * Conjg( N(i1,3) )
   Do i2=i1+1,4
    n_ij(i1,i2) = mN2(i1)/mN2(i2)
    n_ij(i2,i1) = 1._dp / n_ij(i1,i2)
    ZNN(i1,i2) = N(i1,4) * Conjg( N(i2,4) ) - N(i1,3) * Conjg( N(i2,3) )
    ZNN(i2,i1) = Conjg( ZNN(i1,i2) ) 
   End Do
   Do i2=1,6
    s_ai(i2,i1) = mSqd2(i2) / mN2(i1)
    n_bi(i2,i1) = mSlept2(i2) / mN2(i1)
   End Do
  End Do

  Lt = Log(mu**2 / mf_u(3)**2)
  Lua = Log(mu**2 / mSqu2)

  a_sp = oo4Pi * a_s
  !-------------------------------------------
  ! SM contribution
  !-------------------------------------------
  If (l_QCD) Then
   C_L(1) = - 0.25_dp *(f_1_0(x) + a_sp * (f_1_1(x)+ 8._dp* x * Df_1_0(x) * Lt))
   B_L(1) = 0.25_dp *(f_2_0(x) + a_sp * (f_10_1(x) + 8._dp * x * Df_2_0(x) * Lt))
  Else
   C_L(1) = - 0.25_dp * f_1_0(x)
   B_L(1) = 0.25_dp * f_2_0(x)
  End If
                         
  !------------------------------------------
  ! H+ contribution
  !------------------------------------------
  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_2_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  Else
   fact = f_2_0(y) 
  End If

  C_L(2) = 0.125_dp * x * fact / tanb**2 ! m^2_H/m^2_W * y = x
  C_R(2) = - 0.125_dp * mf_d(i) * mf_d(j) * tanb**2 * fact / mW2

  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_7_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  End If

  B_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * fact * mf_l(k)**2 &
      &              / (16._dp * mHp2 * mW2)

  !------------------------------------------------------------------
  ! Chargino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2))                                      &
         &   * ( f_3_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &              + y_ai(i3,i2) * Dyf_3_0(x_ij(i1,i2),y_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &            + y_ai(i3,i2) * Dyf_4_0(x_ij(i1,i2),y_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2)) * f_3_0(x_ij(i1,i2),y_ai(i3,i2))
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2))
     End If
     C_L(3) = C_L(3) + c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i2,i,i3))     &
         &             * ( fac(1) * U(i1,1) * Conjg(U(i2,1))            &
         &               - fac(2) * V(i2,1) * Conjg(V(i1,1)) )
     C_R(3) = C_R(3) + c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i2,i,i3))     &
          &             * ( fac(1) * V(i2,1) * Conjg(V(i1,1))           &
          &               - fac(2) * U(i1,1) * Conjg(U(i2,1)) )
     Do i4=1,3
      If (l_QCD) Then
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      Else
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), v_fi(i4,i2) )
      End If
      B_L(3) = B_L(3) + Conjg(c_CDSu_R(i2,i,i3))* c_CDSu_R(i1,j,i3)      &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_R(i1,k,i4))       &
          &                    * c_CLSn_R(i2,k,i4)                       &
          &           + fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_L(i2,k,i4)    &
          &                    * Conjg(c_CLSn_L(i1,k,i4))    ) / mC2(i2)
      B_R(3) = B_R(3) - Conjg(c_CDSu_L(i2,i,i3))*c_CDSu_L(i1,j,i3)       &
          &         * ( 0.5_dp * fac(1) * Conjg(c_CLSn_L(i1,k,i4))       &
          &                    * c_CLSn_L(i2,k,i4)                       &
          &           + fac(2) * Sqrt(x_ij(i1,i2))* c_CLSn_R(i2,k,i4)    &
          &                    * Conjg(c_CLSn_R(i1,k,i4))    ) / mC2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(y_ai(i2,i1),y_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(y_ai(i2,i1),y_ai(i3,i1))                         &
           &          + 4._dp * (y_ai(i2,i1)*DxF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    +y_ai(i3,i1)*DyF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( y_ai(i2,i1), y_ai(i3,i1))
     End If
     C_L(3) = C_L(3) + fac(1) * ( RSup(i2,1)*Conjg(RSup(i3,1))      &
         &                      + RSup(i2,2)*Conjg(RSup(i3,2))      &
         &                      + RSup(i2,3)*Conjg(RSup(i3,3))  )   &
         &               * c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i1,i,i2))
     C_R(3) = C_R(3) - fac(1) * ( RSup(i3,4)*Conjg(RSup(i2,4))      &
         &                      + RSup(i3,5)*Conjg(RSup(i2,5))      &
         &                      + RSup(i3,6)*Conjg(RSup(i2,6))  )   &
         &                    * c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i1,i,i2))
    End Do
   End Do
    
  End Do
  B_L(3) = kappa_Q *sW2 * B_L(3)
  B_R(3) = kappa_Q *sW2 * B_R(3)
  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(3) = fact * C_L(3)
  C_R(3) = fact * C_R(3)

  !------------------------------------------------------------------
  ! neutralino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,4
   Do i2=1,4
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2))                                      &
         &   * ( f_3_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &              + s_ai(i3,i2) * Dyf_3_0(n_ij(i1,i2),s_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &            + s_ai(i3,i2) * Dyf_4_0(n_ij(i1,i2),s_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2)) * f_3_0(n_ij(i1,i2),s_ai(i3,i2))
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2))
     End If
     fac(1) = fac(1) * ZNN(i1,i2) + fac(2) * ZNN(i2,i1)
     C_L(4) = C_L(4) + fac(1) * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i2,i3))
     C_R(4) = C_R(4) + fac(1) * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i2,i3)) 
     Do i4=1,6
      If (l_QCD) Then
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      Else
       fac(1) = f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      End If
      B_L(4) = B_L(4) + Conjg(c_DNSd_R(i,i2,i3))* c_DNSd_R(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * ( Conjg(c_LNSl_R(k,i1,i4)) * c_LNSl_R(k,i2,i4)   &
          &                 + Conjg(c_LNSl_L(k,i2,i4)) * c_LNSl_L(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_L(k,i2,i4) * Conjg(c_LNSl_L(k,i1,i4))   &
          &                 + c_LNSl_R(k,i1,i4) * Conjg(c_LNSl_R(k,i2,i4)) ) &
          &           ) / mN2(i2)
      B_R(4) = B_R(4) - Conjg(c_DNSd_L(i,i2,i3))*c_DNSd_L(j,i1,i3) &
          &         * ( 0.5_dp * fac(1)                                      &
          &               * ( Conjg(c_LNSl_L(k,i1,i4)) * c_LNSl_L(k,i2,i4)   &
          &                 + Conjg(c_LNSl_R(k,i2,i4)) * c_LNSl_R(k,i1,i4) ) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                           &
          &               * ( c_LNSl_R(k,i2,i4) * Conjg(c_LNSl_R(k,i1,i4))   &
          &                 + c_LNSl_L(k,i1,i4) * Conjg(c_LNSl_L(k,i2,i4)) ) &
          &           ) / mN2(i2)

     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(s_ai(i2,i1),s_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(s_ai(i2,i1),s_ai(i3,i1))                         &
           &          + 4._dp * (s_ai(i2,i1)*DxF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    +s_ai(i3,i1)*DyF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( s_ai(i2,i1), s_ai(i3,i1))
     End If
     C_L(4) = C_L(4) + fac(1) * ( RSdown(i3,4)*Conjg(RSdown(i2,4))      &
         &                      + RSdown(i3,5)*Conjg(RSdown(i2,5))      &
         &                      + RSdown(i3,6)*Conjg(RSdown(i2,6))  )   &
         &                    * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i1,i2))
     C_R(4) = C_R(4) - fac(1) * ( RSdown(i3,1)*Conjg(RSdown(i2,1))      &
         &                      + RSdown(i3,2)*Conjg(RSdown(i2,2))      &
         &                      + RSdown(i3,3)*Conjg(RSdown(i2,3))  )   &
         &                    * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i1,i2))
    End Do
   End Do
    
  End Do

  B_L(4) = kappa_Q *sW2 * B_L(4)
  B_R(4) = kappa_Q *sW2 * B_R(4)

  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(4) = fact * C_L(4)
  C_R(4) = fact * C_R(4)
  !------------------------------------------------------------------
  ! gluino contributions
  !------------------------------------------------------------------
  fact = 4._dp * kappa_Q * e2  / (3._dp * mW2)

   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     Else
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     End If
     C_L(5) = C_L(5) + fac(1) * ( RSdown(i3,4)*Conjg(RSdown(i2,4))      &
         &                      + RSdown(i3,5)*Conjg(RSdown(i2,5))      &
         &                      + RSdown(i3,6)*Conjg(RSdown(i2,6))  )   &
         &                    * c_DGSd_R(j,i3) * Conjg(c_DGSd_R(i,i2))
     C_R(5) = C_R(5) - fac(1) * ( RSdown(i3,1)*Conjg(RSdown(i2,1))      &
         &                      + RSdown(i3,2)*Conjg(RSdown(i2,2))      &
         &                      + RSdown(i3,3)*Conjg(RSdown(i2,3))  )   &
         &                    * c_DGSd_L(j,i3) * Conjg(c_DGSd_L(i,i2))
    End Do
   End Do
   c_L(5) = c_L(5) * fact
   c_R(5) = c_R(5) * fact

  c_a = (Sum(C_L) + Sum(B_L)) / sW2
  c_ap= (Sum(C_R) + sum(B_R)) / sW2

  If (Present(c10_c)) Then
   c10_c(1) = c_a
   c10_c(2:6) = C_L / sW2
   c10_c(2:5) =c10_c(2:5) + B_L / sW2
  End If
  If (Present(c10p_c)) Then
   c10p_c(1) = c_ap
   c10p_c(2:6) = C_R / sW2
   c10p_c(2:5) =c10p_c(2:5) + B_R / sW2
  End If

  Iname = Iname - 1

 End Subroutine C_QdQdLL_AAp


 Subroutine C_QdQdNuNu_LR(i, j, k, mu, mf_d, mf_u, mf_l, mW2, mHp2, mGlu2, mC2 &
            & , mN2, mSnu2, mSlept2, mSqu2, mSqd2, tanb, U, V, N, RSup, RSdown &
            & , c_CDSu_L, c_CDSu_R, c_CNuSl_R, c_NuNSn_R, c_DNSd_L, c_DNSd_R   &
            & , c_DGSd_L, c_DGSd_R, a_s, kappa_q, e2, sW2, l_QCD               &
            & , c11, c11p, C11_c, C11p_c)
 !---------------------------------------------------------------------------
 ! in this routine the coefficients c_L and c_R of the effective operator
 !
 !             4 G_F    alpha
 !    H_eff = -------------------- V_ti V_tj^* (c11 O_L + c11p O_R)
 !             Sqrt(2) 2 Pi s^2_W
 !
 ! for the transitions b -> q nu nu with q=s,d and s -> d nu nu are
 ! calculated. 
 ! The formulas include QCD corrections and are based on
 ! C.Bobeth, A.J.Buras, F.Krueger, J.Urban, NPB 630, (2002) 87
 ! written by Werner Porod, 10.11.2004
 !---------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j              ! indices of the d-type quarks, i>j
  Integer, Intent(in) :: k                 ! generation index of neutrinos
  Real(dp), Intent(in) :: mu               ! renormalization scale
  Real(dp), Intent(in) :: mf_d(3), mf_u(3) ! running d- and u-quark masses at mu
  Real(dp), Intent(in) :: mf_l(3)          ! lepton masses
  Real(dp), Intent(in) :: mW2              ! W-boson mass squared
  Real(dp), Intent(in) :: mHp2             ! mass of charged Higgs boson squared
  Real(dp), Intent(in) :: mglu2            ! gluino mass squared
  Real(dp), Intent(in) :: mC2(2)           ! chargino masses squared
  Real(dp), Intent(in) :: mN2(4)           ! neutralino masses squared
  Real(dp), Intent(in) :: mSnu2(3)         ! sneutrino masses squared
  Real(dp), Intent(in) :: mSlept2(6)       ! slepton masses squared
  Real(dp), Intent(in) :: mSqu2(6)         ! u-squark masses squared
  Real(dp), Intent(in) :: mSqd2(6)         ! d-squark masses squared
  Real(dp), Intent(in) :: a_s              ! alpha_s(mu)
  Real(dp), Intent(in) :: tanb             ! tan(beta) = v_2/v_1
  Real(dp), Intent(in) :: e2               ! = e^2 = g^2 sin^2(theta_W)
  Real(dp), Intent(in) :: sW2              ! sin^2(theta_W)
  ! chargino mixing matrices
  Complex(dp), Intent(in) :: U(2,2), V(2,2)
  ! neutralino mixing matrix
  Complex(dp), Intent(in) :: N(4,4)
  ! squark mixing matrices
  Complex(dp), Intent(in) :: RSup(6,6), RSdown(6,6)
  ! coupling chargino - d-quark - u-squark
  Complex(dp), Intent(in) ::  c_CDSu_L(2,3,6), c_CDSu_R(2,3,6)
  ! coupling neutralino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DNSd_L(3,4,6), c_DNSd_R(3,4,6)
  ! coupling gluino - d-quark - d-squark
  Complex(dp), Intent(in) ::  c_DGSd_L(3,6), c_DGSd_R(3,6)
  ! coupling chargino - neutrino - slepton
  Complex(dp), Intent(in) ::  c_CNuSl_R(2,3,6)
  ! coupling neutralino - neutrino - sneutrino
  Complex(dp), Intent(in) ::  c_NuNSn_R(3,4,3)
  Complex(dp), Intent(in) :: kappa_Q       ! = (8 sqrt(2) G_F e^2 V_ti V^*_tj)^-1
  Logical, Intent(in) :: l_QCd          ! if .true. QCD corrections are included
  
  Complex(dp), Intent(out) :: c11, c11p ! Wilson coefficients
  Complex(dp), Intent(out), Optional :: C11_c(6),C11p_c(6) ! Wilson coefficients
                                                 ! + individual contributions

  Integer :: i1, i2, i3, i4
  Real(dp) :: x, y, z, x_ij(2,2), y_ai(6,2), z_bi(6,2), r_a(6), s_ai(6,4) &
      & , n_bi(3,4), n_ij(4,4), Lt, Lua(6), a_sp
  Complex(dp) :: ZNN(4,4), fact, c_L(5), c_R(5), b_L(5), b_R(5), fac(2)
  !----------------------------------
  ! Initialisation
  !----------------------------------
  Iname = Iname + 1
  NameOfUnit(Iname) = "C_QdQdNuNu_LR"

  c_L = 0._dp
  c_R = 0._dp
  b_L = 0._dp
  b_R = 0._dp

  x = mf_u(3)**2 / mW2
  y = mf_u(3)**2 / mHp2
  z = mHp2 / mW2
  x_ij(1,1) = 1._dp
  x_ij(2,2) = 1._dp
  x_ij(1,2) = mC2(1) / mC2(2)
  x_ij(2,1) = 1._dp / x_ij(1,2)
  r_a = mSqd2/mGlu2
  Do i1=1,2
   Do i2=1,6
    y_ai(i2,i1) = mSqu2(i2) / mC2(i1)
    z_bi(i2,i1) = mSlept2(i2)**2 / mC2(i1)
   End Do
  End Do
  Do i1=1,4
   n_ij(i1,i1) = 1._dp
   ZNN(i1,i1) = N(i1,4) * Conjg( N(i1,4) ) - N(i1,3) * Conjg( N(i1,3) )
   Do i2=i1+1,4
    n_ij(i1,i2) = mN2(i1)/mN2(i2)
    n_ij(i2,i1) = 1._dp / n_ij(i1,i2)
    ZNN(i1,i2) = N(i1,4) * Conjg( N(i2,4) ) - N(i1,3) * Conjg( N(i2,3) )
    ZNN(i2,i1) = Conjg( ZNN(i1,i2) ) 
   End Do
   Do i2=1,3
    s_ai(i2,i1) = mSqd2(i2) / mN2(i1)
    s_ai(i2+3,i1) = mSqd2(i2+3) / mN2(i1)
    n_bi(i2,i1) = mSnu2(i2) / mN2(i1)
   End Do
  End Do
  Lt = Log(mu**2 / mf_u(3)**2)
  Lua = Log(mu**2 / mSqu2)

  a_sp = oo4Pi * a_s
  !-------------------------------------------
  ! SM contribution
  !-------------------------------------------
  If (l_QCD) Then
   C_L(1) = 0.25_dp * (f_1_0(x) + a_sp* (f_1_1(x) + 8._dp* x* Df_1_0(x) * Lt) )
   B_L(1) = - f_2_0(x) - a_sp * (f_6_1(x) + 8._dp * x * Df_2_0(x) * Lt)
  Else
   C_L(1) = 0.25_dp * f_1_0(x)
   B_L(1) = - f_2_0(x)
  End If
  !------------------------------------------
  ! H+ contribution
  !------------------------------------------
  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_2_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  Else
   fact = f_2_0(y) 
  End If

  C_L(2) = - 0.125_dp * x * fact / tanb**2 ! m^2_H/m^2_W * y = x
  C_R(2) = 0.125_dp * mf_d(i) * mf_d(j) * tanb**2 * fact / mW2

  If (l_QCD) Then
   fact = f_2_0(y) * (1._dp + 8._dp * a_sp * Lt) &
      & + a_sp * (f_7_1(y) + 8._dp * y * Df_2_0(y) * Lt)
  End If

  B_R(2) = - mf_d(i) * mf_d(j) * tanb**4 * fact * mf_l(k)**2 &
         &           / (16._dp * mHp2 * mW2)

  !------------------------------------------------------------------
  ! Chargino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,2
   Do i2=1,2
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2))                                      &
         &   * ( f_3_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &              + y_ai(i3,i2) * Dyf_3_0(x_ij(i1,i2),y_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(x_ij(i1,i2),y_ai(i3,i2))                      &
         &            + y_ai(i3,i2) * Dyf_4_0(x_ij(i1,i2),y_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(x_ij(i1,i2)) * f_3_0(x_ij(i1,i2),y_ai(i3,i2))
      fac(2) = f_4_0(x_ij(i1,i2),y_ai(i3,i2))
     End If
     ! Z contribution
     C_L(3) = C_L(3) - c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i2,i,i3))   &
         &            * ( fac(1) * U(i1,1) * Conjg(U(i2,1))                  &
         &              - fac(2) * V(i2,1) * Conjg(V(i1,1)) )
     C_R(3) = C_R(3) - c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i2,i,i3))   &
         &            * ( fac(1) * V(i2,1) * Conjg(V(i1,1))                  &
         &              - fac(2) * U(i1,1) * Conjg(U(i2,1)) )
     Do i4=1,6
      If (l_QCD) Then
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
      Else
       fac(1) = f_5_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
       fac(2) = f_6_0(x_ij(i1,i2), y_ai(i3,i2), z_bi(i4,i2) )
      End If
      B_L(3) = B_L(3) + fac(1)                    &
          &       * Conjg(c_CNuSl_R(i1,k,i4) * c_CDSu_R(i2,i,i3) )   &
          &       * c_CNuSl_R(i2,k,i4) * c_CDSu_R(i1,j,i3) / mC2(i2)
      B_R(3) = B_R(3) - fac(2) * Sqrt(x_ij(i1,i2))   &
          &       * Conjg(c_CNuSl_R(i1,k,i4) * c_CDSu_L(i2,i,i3) )     &
          &       * c_CNuSl_R(i2,k,i4) * c_CDSu_L(i1,j,i3) / mC2(i2)
     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(y_ai(i2,i1),y_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(y_ai(i2,i1),y_ai(i3,i1))                         &
           &          + 4._dp * (y_ai(i2,i1)*DxF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    +y_ai(i3,i1)*DyF_4_0(y_ai(i2,i1),y_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( y_ai(i2,i1), y_ai(i3,i1))
     End If
     C_L(3) = C_L(3) - fac(1) * ( RSup(i2,1)*Conjg(RSup(i3,1))      &
         &                      + RSup(i2,2)*Conjg(RSup(i3,2))      &
         &                      + RSup(i2,3)*Conjg(RSup(i3,3))  )   &
         &                    * c_CDSu_R(i1,j,i3) * Conjg(c_CDSu_R(i1,i,i2))
     C_R(3) = C_R(3) + fac(1) * ( RSup(i3,4)*Conjg(RSup(i2,4))      &
         &                      + RSup(i3,5)*Conjg(RSup(i2,5))      &
         &                      + RSup(i3,6)*Conjg(RSup(i2,6))  )   &
         &                    * c_CDSu_L(i1,j,i3) * Conjg(c_CDSu_L(i1,i,i2))
    End Do
   End Do
       
  End Do

  B_L(3) = 0.5_dp * kappa_Q * sW2 * B_L(3)
  B_R(3) = kappa_Q * sW2 * B_R(3)

  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(3) = fact * C_L(3)
  C_R(3) = fact * C_R(3)

  !------------------------------------------------------------------
  ! neutralino contributions
  ! please note that the couplings to sfermions are complex conjugated
  ! and left-right exchanged compared to my notation
  !------------------------------------------------------------------
  Do i1=1,4
   Do i2=1,4
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2))                                      &
         &   * ( f_3_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &     + a_sp * ( f_3_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &              + s_ai(i3,i2) * Dyf_3_0(n_ij(i1,i2),s_ai(i3,i2))      &
         &                * 4._dp * Lua(i3) ) )
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2)) * (1._dp+4._dp*a_sp*Lua(i3))  &
         &   + a_sp * ( f_4_1(n_ij(i1,i2),s_ai(i3,i2))                      &
         &            + s_ai(i3,i2) * Dyf_4_0(n_ij(i1,i2),s_ai(i3,i2))      &
                        * 4._dp * Lua(i3))
     Else
      fac(1) = 2._dp * Sqrt(n_ij(i1,i2)) * f_3_0(n_ij(i1,i2),s_ai(i3,i2))
      fac(2) = f_4_0(n_ij(i1,i2),s_ai(i3,i2))
     End If
     fac(1) = - fac(1) * ZNN(i1,i2) - fac(2) * ZNN(i2,i1)
     C_L(4) = C_L(4) + fac(1) * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i2,i3))
     C_R(4) = C_R(4) + fac(1) * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i2,i3))
     Do i4=1,3
      If (l_QCD) Then
       fac(1) = 0.5_dp * f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      Else
       fac(1) = 0.5_dp * f_5_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
       fac(2) = f_6_0(n_ij(i1,i2), s_ai(i3,i2), n_bi(i4,i2) )
      End If
      B_L(4) = B_L(4) + Conjg(c_DNSd_R(i,i2,i3))* c_DNSd_R(j,i1,i3) &
          &         * ( fac(1) * Conjg(c_NuNSn_R(k,i1,i4)) * c_NuNSn_R(k,i2,i4) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                              &
          &               * c_NuNSn_R(k,i1,i4) * Conjg(c_NuNSn_R(k,i2,i4))      &
          &           ) / mN2(i2)
      B_R(4) = B_R(4) - Conjg(c_DNSd_L(i,i2,i3))*c_DNSd_L(j,i1,i3) &
          &         * ( fac(1) * Conjg(c_NuNSn_R(k,i2,i4)) * c_NuNSn_R(k,i1,i4) &
          &           + fac(2) * Sqrt(n_ij(i1,i2))                              &
          &               * c_NuNSn_R(k,i2,i4) * Conjg(c_NuNSn_R(k,i1,i4))      &
          &           ) / mN2(i2)
     End Do
    End Do
   End Do
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0(s_ai(i2,i1),s_ai(i3,i1)) * (1._dp+4._dp*a_sp*Lua(i3))    &
           & + a_sp * (F_5_1(s_ai(i2,i1),s_ai(i3,i1))                         &
           &          + 4._dp * (s_ai(i2,i1)*DxF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    +s_ai(i3,i1)*DyF_4_0(s_ai(i2,i1),s_ai(i3,i1)) &
           &                    ) *Lua(i2) ) 
     Else
      fac(1) = F_4_0( s_ai(i2,i1), s_ai(i3,i1))
     End If
     C_L(4) = C_L(4) - fac(1) * ( RSdown(i3,4)*Conjg(RSdown(i2,4))      &
         &                      + RSdown(i3,5)*Conjg(RSdown(i2,5))      &
         &                      + RSdown(i3,6)*Conjg(RSdown(i2,6))  )   &
         &                    * c_DNSd_R(j,i1,i3) * Conjg(c_DNSd_R(i,i1,i2))
     C_R(4) = C_R(4) + fac(1) * ( RSdown(i3,1)*Conjg(RSdown(i2,1))      &
         &                      + RSdown(i3,2)*Conjg(RSdown(i2,2))      &
         &                      + RSdown(i3,3)*Conjg(RSdown(i2,3))  )   &
         &                    * c_DNSd_L(j,i1,i3) * Conjg(c_DNSd_L(i,i1,i2))
    End Do
   End Do
       
  End Do

  B_L(4) = kappa_Q *sW2 * B_L(4)
  B_R(4) = kappa_Q *sW2 * B_R(4)

  fact = kappa_Q * e2 / (4._dp * mW2)
  C_L(4) = fact * C_L(4)
  C_R(4) = fact * C_R(4)

  !------------------------------------------------------------------
  ! gluino contributions
  !------------------------------------------------------------------
   Do i2=1,6
    Do i3=1,6
     If (l_QCD) Then
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     Else
      fac(1) = F_4_0( r_a(i2), r_a(i3))
     End If
     C_L(5) = C_L(5) - fac(1) * ( RSdown(i3,4)*Conjg(RSdown(i2,4))      &
         &                      + RSdown(i3,5)*Conjg(RSdown(i2,5))      &
         &                      + RSdown(i3,6)*Conjg(RSdown(i2,6))  )   &
         &                    * c_DGSd_R(j,i3) * Conjg(c_DGSd_R(i,i2))
     C_R(5) = C_R(5) + fac(1) * ( RSdown(i3,1)*Conjg(RSdown(i2,1))      &
         &                      + RSdown(i3,2)*Conjg(RSdown(i2,2))      &
         &                      + RSdown(i3,3)*Conjg(RSdown(i2,3))  )   &
         &                    * c_DGSd_L(j,i3) * Conjg(c_DGSd_L(i,i2))
    End Do
   End Do
  fact = 4._dp * kappa_Q * e2  / (3._dp * mW2)
  c_L(5) = c_L(5) * fact
  c_R(5) = c_R(5) * fact

  C11 = Sum(C_L) + sum(B_L)
  C11p = Sum(C_R) + sum(B_R)

  If (Present(c11_c)) then
   C11_c(1) = C11
   C11_c(2:6) = C_L + B_L
  end if
  If (Present(c11p_c)) then
   C11p_c(1) = C11p
   C11p_c(2:6) = C_R + B_R
  end if

  Iname = Iname - 1

 End Subroutine C_QdQdNuNu_LR


 Subroutine Delta_F2_Boxes(i, j, T3, g, Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R    &
     & , mf_u, mf_d, mC, U, V, mN, N, mGlu, phi_g &
     & , mSpm2, RSpm, mSup2, RSup, mSdown2, RSdown            &
     & , B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2, B_SRR1, B_SRR2  )
 !---------------------------------------------------------------------------
 ! gives the box contribution to Delta_F =2 box contributions:
 !     B_VLL  -> bar(q)_Li \gamma_mu q_Lj  * \bar(q)_Li \gamma^mu q_Lj
 !     B_VRR  -> bar(q)_Ri \gamma_mu q_Rj  * \bar(q)_Ri \gamma^mu q_Rj
 !     B_LR1  -> bar(q)_Li \gamma_mu q_Lj  * \bar(q)_Ri \gamma^mu q_Rj
 !     B_LR2  -> bar(q)_Ri q_Lj  * \bar(q)_Li q_Rj
 !     B_SLL1 -> bar(q)_Ri q_Lj  * \bar(q)_Ri q_Lj
 !     B_SLL2 -> bar(q)_Ri \sigma^\mu\nu q_Lj  * \bar(q)_Ri \sigma_\mu\nu q_Lj 
 !     B_SRR1 -> bar(q)_Li q_Rj  * \bar(q)_Li q_Rj
 !     B_SRR2 -> bar(q)_Li \sigma^\mu\nu q_Rj  * \bar(q)_Li \sigma_\mu\nu q_Rj 
 ! using formulas of Buras et al., NPB659 (2003) 3, for charged Higgs and
 ! charginos, Appendix A4
 ! for neutralinos/gluinos the formulas of S.Baek et al. are used,
 ! PRD64 (2001) 095001
 ! Input: i, j .... indices of the outside quarks
 !        T3 ...... isospin of the outside left quarks
 !         mC, mN, mGlu, mSpm, mS0, mP0, RS0, RP0, mSup2, mSdown2
 !         U, V, N, phi_g, RSpm, RSup, RSdown
 !        
 ! Output: G^2_F m^2_W * (V^tb* V^tq)**2 *
 !               (C_VLL, C_VRR, C_SLL1, C_SLL2, C_SRR1, C_SRR2, C_LR1, C_LR2)
 !         
 ! written by Werner Porod, 08 Jul 2003
 !   08.07.03: charged Higgs, for down-quarks
 !---------------------------------------------------------------------------
 implicit none
  Real(dp), Intent(in) ::  mC(2), mN(4), mGlu, mSpm2(2), mSup2(6)        &
     & , mSdown2(6), T3, g(3), mf_u(3), mf_d(3)
  Complex(dp), Intent(in) ::  U(2,2), V(2,2), N(4,4), phi_g, RSpm(2,2)      &
     & , RSup(6,6), RSdown(6,6)
  Complex(dp), Intent(in), Dimension(3,3) :: Y_u, Y_d, Ru_L, Ru_R, Rd_L, Rd_R
  Integer, Intent(in) :: i, j

  Complex(dp), Intent(out) ::  B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2   &
     & , B_SRR1, B_SRR2

  Real(dp) :: D0m2, D2m2, mq2(3), mq(3), mqij, mC2(2), mSf2(6), mSf2p(6)    &
     & , mf(3), mf2(3), mg2, mN2(4)
  Complex(dp) :: c_Hqqp_L(2,3,3), c_Hqqp_R(2,3,3), coupLC, coupRC           &
     & , c_Wqqp_L(3,3), c_CQSq_L(2,3,6), c_CQSq_R(2,3,6), c_QGSq_L(3,6)     &
     & , c_QGSq_R(3,6), c_QNSq_L(3,4,6), c_QNSq_R(3,4,6) 
  Integer :: i1, i2, i3, i4
!  Logical ::
!  Character ::

  B_VLL = 0._dp
  B_VRR = 0._dp
  B_LR1 = 0._dp
  B_LR2 = 0._dp
  B_SLL1 = 0._dp
  B_SLL2 = 0._dp
  B_SRR1 = 0._dp
  B_SRR2 = 0._dp

  if (T3.le.0._dp) then
   c_Wqqp_L = - oosqrt2 * g(2) * CKM
  else
   c_Wqqp_L = - oosqrt2 * g(2) * Transpose(Conjg(CKM))
  end if
  Do i1=1,3
   Do i2=1,3
    Do i3=1,2
     call CoupChargedScalarFermion(i3, i1, i2, RSpm, Y_d, Rd_L, Rd_R, Y_u &
                     &, Ru_L, Ru_R, coupLC, coupRC)
     if (T3.lt.0._dp) then
      c_Hqqp_L(i3,i1,i2) = Conjg(coupRC)
      c_Hqqp_R(i3,i1,i2) = Conjg(coupLC)
     else
      c_Hqqp_L(i3,i2,i1) = coupLC
      c_Hqqp_R(i3,i2,i1) = coupRC
     end if
    end do
   end do
  end do

  if (T3.lt.0._dp) then
   mq = mf_u
   mq2 = mq**2
   mf = mf_d
   mf2 = mf**2
   mSf2 = mSdown2
   mSf2p = mSup2
   Do i1=1,2
    Do i2=1,3
     Do i3=1,6
      call CoupCharginoSfermion3(i1, i2, i3, g(2), T3, RSup, Y_d, Y_u, Rd_L &
                         & , Rd_R, U, V, coupLC, coupRC)
      !------------------------------------
      ! using the notation of Buras et al.
      !------------------------------------
      c_CQSq_L(i1,i2,i3) = Conjg(coupRC)
      c_CQSq_R(i1,i2,i3) = Conjg(coupLC)
     end do
    end do
   end do
   Do i1=1,3
    Do i2=1,6
     call CoupGluinoSquark(g(3), phi_g, i1, i2, RSdown, Rd_L, Rd_R  &
                         & , coupLC, coupRC)
     c_QGSq_L(i1,i2) = coupLC
     c_QGSq_R(i1,i2) = coupRC
     Do i3=1,4
      call CoupNeutralinoSdown(i1, i3, i2, g(1), g(2), RSdown, Rd_L, Rd_R &
                         & , Y_d, N, coupLC, coupRC)
      c_QNSq_L(i1,i3,i2) = coupLC
      c_QNSq_R(i1,i3,i2) = coupRC
     end do
    end do
   end do
  else ! T3 >0

   mq = mf_d
   mq2 = mq**2
   mf = mf_u
   mf2 = mf**2

   mSf2p = mSdown2
   mSf2 = mSup2
   Do i1=1,2
    Do i2=1,3
     Do i3=1,6
      call CoupCharginoSfermion3(i1, i2, i3, g(2), T3, RSdown, Y_d, Y_u, Ru_L &
                         & , Ru_R, U, V, coupLC, coupRC)
      !------------------------------------
      ! using the notation of Buras et al.
      !------------------------------------
      c_CQSq_L(i1,i2,i3) = Conjg(coupRC)
      c_CQSq_R(i1,i2,i3) = Conjg(coupLC)
     end do
    end do
   end do
   Do i1=1,3
    Do i2=1,6
     call CoupGluinoSquark(g(3), phi_g, i1, i2, RSup, Ru_L, Ru_R  &
                         & , coupLC, coupRC)
     c_QGSq_L(i1,i2) = coupLC
     c_QGSq_R(i1,i2) = coupRC
     Do i3=1,4
      call CoupNeutralinoSup(i1, i3, i2, g(1), g(2), RSup, Ru_L, Ru_R &
                         & , Y_u, N, coupLC, coupRC)
      c_QNSq_L(i1,i3,i2) = coupLC
      c_QNSq_R(i1,i3,i2) = coupRC
     end do
    end do
   end do
  end if   
!if (l_mbq(9))  return
  !------------------------------
  ! W+ H+ + H+ H+ contributions
  !------------------------------

  Do i1=1,3
   Do i2=1,3
    mqij = mq(i1) * mq(i2)
    Do i3=1,2
     Do i4=1,2
      D0m2 = - mqij * D0_Bagger(mSpm2(i3), mSpm2(i4), mq2(i1), mq2(i2) )
      D2m2 = 4._dp * D27_Bagger(mSpm2(i3), mSpm2(i4), mq2(i1), mq2(i2) )
      If ((i3.Eq.1).and.(i4.eq.2)) then
       B_VLL = B_VLL - Conjg(c_Wqqp_L(i1,i)) * c_Wqqp_L(i2,j)                &
            &           * Conjg(c_Hqqp_L(2,i,i2)) * c_Hqqp_L(2,j,i1) * D0m2  &
            &        + 0.25_dp *  Conjg(c_Hqqp_L(1,i,i1)) * c_Hqqp_L(1,j,i2) &
            &           * Conjg(c_Hqqp_L(2,i,i2)) * c_Hqqp_L(2,j,i1) * D2m2  
      end if
      If ((i3.Eq.2).and.(i4.eq.2)) then
       B_VLL = B_VLL + 0.125_dp * Conjg(c_Hqqp_L(2,i,i1)) * c_Hqqp_L(2,j,i2) &
            &           * Conjg(c_Hqqp_L(2,i,i2)) * c_Hqqp_L(2,j,i1) * D2m2  
       end if
      B_VRR = B_VRR + 0.125_dp * Conjg(c_Hqqp_R(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &          * Conjg(c_Hqqp_R(i3,i,i2)) * c_Hqqp_R(i3,j,i1) * D2m2
      B_SLL1 = B_SLL1 + 0.5_dp * Conjg(c_Hqqp_R(i4,i,i1)) * c_Hqqp_L(i4,j,i2) &
            &      * Conjg(c_Hqqp_R(i3,i,i2)) * c_Hqqp_L(i3,j,i1) * D0m2
      B_SRR1 = B_SRR1 + 0.5_dp * Conjg(c_Hqqp_L(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &      * Conjg(c_Hqqp_L(i3,i,i2)) * c_Hqqp_R(i3,j,i1) * D0m2
      B_LR1 = B_LR1 + 0.25_dp * Conjg(c_Hqqp_L(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &         * Conjg(c_Hqqp_R(i3,i,i2)) * c_Hqqp_L(i3,j,i1) * D2m2
      B_LR2 = B_LR2 + Conjg(c_Hqqp_R(i4,i,i1)) * c_Hqqp_R(i4,j,i2) &
            &     * Conjg(c_Hqqp_L(i3,i,i2)) * c_Hqqp_L(i3,j,i1) * D0m2
      If (i3.Eq.1) then
       B_LR2 = B_LR2 - Conjg(c_Wqqp_L(i1,i)) * c_Wqqp_L(i2,j) &
            &            * Conjg(c_Hqqp_R(i4,i,i2)) * c_Hqqp_R(i4,j,i1) * D2m2
      end if
     end do
    end do
   end do
  end do
  !------------------------------
  ! Chargino contributions
  !------------------------------
  mC2 = mC**2
  Do i1=1,2
   Do i2=1,2
    mqij = mC(i1) * mC(i2)
    Do i3=1,6
     Do i4=1,6
      D0m2 = - mqij * D0_Bagger(mC2(i1), mC2(i2), mSf2p(i3), mSf2p(i4) )
      D2m2 = 4._dp * D27_Bagger(mC2(i1), mC2(i2), mSf2p(i3), mSf2p(i4) )
      B_VLL = B_VLL + 0.125_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_L(i1,i,i4)) &
            &                  * c_CQSq_L(i1,j,i3) * c_CQSq_L(i2,j,i4) * D2m2
      B_VRR = B_VRR + 0.125_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_R(i1,i,i4)) &
            &                  * c_CQSq_R(i1,j,i3) * c_CQSq_R(i2,j,i4) * D2m2
      B_SLL1 = B_SLL1 - 0.25_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_R(i1,i,i4)) &
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_L(i2,j,i4) * D0m2
      B_SLL2 = B_SLL2 + 0.0625_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_R(i1,i,i4))&
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_L(i2,j,i4) * D0m2
      B_SRR1 = B_SRR1 - 0.25_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_L(i1,i,i4)) &
            &           * c_CQSq_R(i1,j,i3) * c_CQSq_R(i2,j,i4) * D0m2
      B_SRR2 = B_SRR2 + 0.0625_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_L(i1,i,i4))&
            &           * c_CQSq_R(i1,j,i3) * c_CQSq_R(i2,j,i4) * D0m2
      B_LR1 = B_LR1 - 0.5_dp * Conjg(c_CQSq_L(i2,i,i3)*c_CQSq_R(i1,i,i4)) &
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_R(i2,j,i4) * D0m2
      B_LR2 = B_LR2 - 0.5_dp * Conjg(c_CQSq_R(i2,i,i3)*c_CQSq_L(i1,i,i4)) &
            &           * c_CQSq_L(i1,j,i3) * c_CQSq_R(i2,j,i4) * D2m2
     end do
    end do
   end do
  end do

  !-----------------------------------------------------------
  ! gluino contributions
  !-----------------------------------------------------------
  mg2 = mglu**2
  Do i1=1,6
   Do i2=1,6
    D0m2 = mG2 * D0_Bagger(mg2, mg2, mSf2(i1), mSf2(i2) )
    D2m2 = D27_Bagger(mg2, mg2, mSf2(i1), mSf2(i2) )
    B_VLL = B_VLL - 4._dp *(D0m2 + 11*D2m2)                      &
          &               * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2)) &
          &               * c_QGSq_R(j,i1)*c_QGSq_R(j,i2) / 9._dp
    B_VRR = B_VRR - 4._dp *(D0m2 + 11*D2m2)                      &
          &               * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2)) &
          &               * c_QGSq_L(j,i1)*c_QGSq_L(j,i2) / 9._dp
    B_SLL1 = B_SLL1 - 37._dp * D0m2 * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2)) &
          &               * c_QGSq_R(j,i1)*c_QGSq_R(j,i2) / 9._dp
    B_SLL2 = B_SLL2 - D0m2 * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2)) &
          &                * c_QGSq_R(j,i1)*c_QGSq_R(j,i2) / 12._dp
    B_SRR1 = B_SRR1 - 37._dp * D0m2 * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2)) &
          &               * c_QGSq_L(j,i1)*c_QGSq_L(j,i2) / 9._dp
    B_SRR2 = B_SRR2 - D0m2 * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2)) &
          &                * c_QGSq_L(j,i1)*c_QGSq_L(j,i2) / 12._dp
    B_LR1 = B_LR1 + 2._dp * (20._dp * D2m2                             &
          &                    * Conjg(c_QGSq_L(i,i2)*c_QGSq_R(i,i1))  &
          &                    * c_QGSq_L(j,i2)*c_QGSq_R(j,i1)         &
          &                 +  Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2))    &
          &                    *  (-60._dp * D2m2                      &
          &                            *c_QGSq_L(j,i2)*c_QGSq_R(j,i1)  & 
          &                       + D0m2                               &
          &                           *c_QGSq_L(j,i1)*c_QGSq_R(j,i2) ) &
          &                  ) / 9._dp
    B_LR2 = B_LR2 + ( 48._dp * D2m2                                &
          &              * Conjg(c_QGSq_L(i,i2)*c_QGSq_R(i,i1))    &
          &              * c_QGSq_L(j,i2)*c_QGSq_R(j,i1)           &
          &         + 4._dp * Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2)) &
          &              * (22._dp * D2m2                          &
          &                        * c_QGSq_L(j,i2)*c_QGSq_R(j,i1) & 
          &                -21._dp * D0m2                          &
          &                        * c_QGSq_L(j,i1)*c_QGSq_R(j,i2) &
          &                 ) ) / 9._dp
   end do 
  end do
! if (l_mbq(12))  return

  !-------------------------------------------------------------------
  ! gluino/neutralino contributions, the sign in front B27_Bagger is 
  ! the relativ sign to S.Baek et al.
  !-------------------------------------------------------------------
   mN2 = mN**2
  Do i1=1,6
   Do i2=1,6
    Do i3=1,4
     D0m2 = mN(i3) * mGlu * D0_Bagger(mN2(i3), mg2, mSf2(i1), mSf2(i2) )
     D2m2 = D27_Bagger(mN2(i3), mg2, mSf2(i1), mSf2(i2) )
     B_VLL = B_VLL - ( D0m2 * Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2)) &
           &                * c_QGSq_R(j,i1)*c_QGSq_R(j,i2)              & 
           &         + (4._dp * D2m2 * c_QGSq_R(j,i1)                    &
           &                  * Conjg(c_QGSq_R(i,i2)*c_QNSq_R(i,i3,i1))  &
           &           + D0m2 * c_QNSq_R(j,i3,i1)                        &
           &                  * Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2))     &
           &           ) * c_QNSq_R(j,i3,i2)                             &
           &         ) / 3._dp
     B_VRR = B_VRR - ( D0m2 * Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2)) &
           &                * c_QGSq_L(j,i1)*c_QGSq_L(j,i2)              & 
           &         + (4._dp * D2m2 * c_QGSq_L(j,i1)                    &
           &                  * Conjg(c_QGSq_L(i,i2)*c_QNSq_L(i,i3,i1))  &
           &           + D0m2 * c_QNSq_L(j,i3,i1)                        &
           &                  * Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2))     &
           &           ) * c_QNSq_L(j,i3,i2)                             &
           &         ) / 3._dp
     B_SRR1 = B_SRR1 + D0m2 * ( Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2))      &
           &                      *c_QGSq_L(j,i1)*c_QGSq_L(j,i2)                &
           &                  + (-7._dp*Conjg(c_QGSq_R(i,i2)*c_QNSq_R(i,i3,i1)) &
           &                       *c_QGSq_L(j,i1)                              &
           &                    + Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2))          &
           &                      *c_QNSq_L(j,i3,i1) ) *c_QNSq_L(j,i3,i2)       &
           &                  ) / 3._dp
     B_SRR2 = B_SRR2 - D0m2 * ( Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2))  &
           &                    * c_QGSq_L(j,i1)*c_QGSq_L(j,i2)             &
           &                  + ( Conjg(c_QGSq_R(i,i2)*c_QNSq_R(i,i3,i1))   &
           &                     * c_QGSq_L(j,i1)                           & 
           &                    + Conjg(c_QGSq_R(i,i1)*c_QGSq_R(i,i2))      &
           &                      *c_QNSq_L(j,i3,i1) ) * c_QNSq_L(j,i3,i2)  &
            &                 ) / 12._dp
     B_SLL1 = B_SLL1 + D0m2 * ( Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2))      &
           &                      *c_QGSq_R(j,i1)*c_QGSq_R(j,i2)                &
           &                  + (-7._dp*Conjg(c_QGSq_L(i,i2)*c_QNSq_L(i,i3,i1)) &
           &                       *c_QGSq_R(j,i1)                              &
           &                    + Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2))          &
           &                      *c_QNSq_R(j,i3,i1) ) *c_QNSq_R(j,i3,i2)       &
           &                  ) / 3._dp
     B_SLL2 = B_SLL2 - D0m2 * ( Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2))  &
           &                    * c_QGSq_R(j,i1)*c_QGSq_R(j,i2)             &
           &                  + ( Conjg(c_QGSq_L(i,i2)*c_QNSq_L(i,i3,i1))   &
           &                     * c_QGSq_R(j,i1)                           & 
           &                    + Conjg(c_QGSq_L(i,i1)*c_QGSq_L(i,i2))      &
           &                      *c_QNSq_R(j,i3,i1) ) * c_QNSq_R(j,i3,i2)  &
            &                 ) / 12._dp
     B_LR1 = B_LR1 - 2._dp * ( ( D0m2                                          &
           &                    * ( Conjg(c_QGSq_L(i,i2)*c_QNSq_R(i,i3,i1))    &
           &                          * c_QGSq_R(j,i1)*c_QNSq_L(j,i3,i2)       &
           &                      + Conjg(c_QGSq_R(i,i1)*c_QNSq_L(i,i3,i2))    &
           &                          *c_QGSq_L(j,i2)* c_QNSq_R(j,i3,i1)       &
           &                      ) ) / 6._dp                                  &
           &                  + D2m2                                           &
           &                    * (( Conjg(c_QNSq_R(i,i3,i2))*c_QGSq_L(j,i2)   &
           &                       + Conjg(c_QGSq_R(i,i2))*c_QNSq_L(j,i3,i2))  &
           &                      * (Conjg(c_QNSq_L(i,i3,i1))*c_QGSq_R(j,i1)   &
           &                        + Conjg(c_QGSq_L(i,i1))*c_QNSq_R(j,i3,i1)) &
           &                      + (Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_R(i,i3,i2))&
           &                          *c_QGSq_L(j,i1)* c_QGSq_R(j,i2)          &
           &                        + Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2))     &
           &                             *c_QNSq_L(j,i3,i1)*c_QNSq_R(j,i3,i2) )&
           &                        / 3._dp )                                   &
           &                      ) 
     B_LR2 = B_LR2 - 4._dp * ( D0m2                                            &
           &                   * ( Conjg(c_QGSq_L(i,i2)*c_QNSq_R(i,i3,i1))     &
           &                       * c_QGSq_R(j,i1)*c_QNSq_L(j,i3,i2)          &
           &                     + Conjg(c_QGSq_R(i,i1)*c_QNSq_L(i,i3,i2))     &
           &                         *c_QGSq_L(j,i2)* c_QNSq_R(j,i3,i1)        &
           &                     )  / 2._dp                                    &
           &                 + D2m2                                            &
           &                   * ( Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_R(i,i3,i2))  &
           &                        * c_QGSq_L(j,i1)*c_QGSq_R(j,i2)            &
           &                     + (( Conjg(c_QNSq_R(i,i3,i2))*c_QGSq_L(j,i2)  &
           &                        + Conjg(c_QGSq_R(i,i2))*c_QNSq_L(j,i3,i2)) &
           &                      * (Conjg(c_QNSq_L(i,i3,i1))*c_QGSq_R(j,i1)   &
           &                        + Conjg(c_QGSq_L(i,i1))*c_QNSq_R(j,i3,i1)) &
           &                       ) /3._dp                                    &
           &                     + Conjg(c_QGSq_L(i,i1)*c_QGSq_R(i,i2))        &
           &                        *c_QNSq_L(j,i3,i1)* c_QNSq_R(j,i3,i2)      &
           &                       )  )
       
    end do 
   end do 
  end do

  !-------------------------------------------------------------------
  ! neutralino contributions, the sign in front B27_Bagger is 
  ! the relativ sign to S.Baek et al.
  !-------------------------------------------------------------------
  Do i1=1,6
   Do i2=1,6
    Do i3=1,4
     Do i4=1,4
      D0m2 = mN(i3) * mN(i4) * D0_Bagger(mN2(i3), mN2(i4), mSf2(i1), mSf2(i2) )
      D2m2 = D27_Bagger(mN2(i3), mN2(i4), mSf2(i1), mSf2(i2) )
      B_VLL = B_VLL - Conjg(c_QNSq_R(i,i3,i2))*c_QNSq_R(j,i4,i2)              &
            &    * ( 2._dp * D2m2 *Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_R(j,i3,i1) & 
            &      + D0m2*Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_R(j,i4,i1) ) / 4._dp
      B_VRR = B_VRR - Conjg(c_QNSq_L(i,i3,i2))*c_QNSq_L(j,i4,i2)              &
            &    * ( 2._dp * D2m2 *Conjg(c_QNSq_L(i,i4,i1))*c_QNSq_L(j,i3,i1) & 
            &      + D0m2*Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_L(j,i4,i1) ) / 4._dp
      B_SLL1 = B_SLL1 + 0.25_dp * D0m2 * c_QNSq_R(j,i4,i2)                     &
             &           * (2._dp *Conjg(c_QNSq_L(i,i3,i1)*c_QNSq_L(i,i3,i2))  &
             &                    *c_QNSq_R(j,i4,i1)                           & 
             &             - Conjg(c_QNSq_L(i,i3,i2))                          &
             &                * ( Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_R(j,i4,i1)   & 
             &                  - Conjg(c_QNSq_L(i,i4,i1))*c_QNSq_R(j,i3,i1) ) ) 
      B_SLL2 = B_SLL2 - D0m2 * Conjg(c_QNSq_L(i,i3,i2))*c_QNSq_R(j,i4,i2)    &
             &            * ( Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_R(j,i4,i1)     & 
             &              - Conjg(c_QNSq_L(i,i4,i1))*c_QNSq_R(j,i3,i1) )   &
             &            / 16._dp
      B_SRR1 = B_SRR1 + 0.25_dp * D0m2 * c_QNSq_L(j,i4,i2)                     &
             &           * (2._dp *Conjg(c_QNSq_R(i,i3,i1)*c_QNSq_R(i,i3,i2))  &
             &                    *c_QNSq_L(j,i4,i1)                           & 
             &             - Conjg(c_QNSq_R(i,i3,i2))                          &
             &                * ( Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_L(j,i4,i1)   & 
             &                  - Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_L(j,i3,i1) ) ) 
      B_SRR2 = B_SRR2 - D0m2 * Conjg(c_QNSq_R(i,i3,i2))*c_QNSq_L(j,i4,i2)    &
             &            * ( Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_L(j,i4,i1)     & 
             &              - Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_L(j,i3,i1) )   &
             &            / 16._dp
      B_LR1 = B_LR1 + 0.5_dp * Conjg(c_QNSq_R(i,i3,i2))*c_QNSq_R(j,i4,i2)     &
            &    * ( 2._dp * D2m2*Conjg(c_QNSq_L(i,i3,i1))*c_QNSq_L(j,i4,i1) & 
            &      + D0m2 * Conjg(c_QNSq_L(i,i4,i1)) * c_QNSq_L(j,i3,i1) )
      B_LR2 = B_LR2 + 2._dp*D2m2*Conjg(c_QNSq_L(i,i3,i2))*c_QNSq_R(j,i4,i2)  &
            &           * ( Conjg(c_QNSq_R(i,i3,i1))*c_QNSq_L(j,i4,i1)       &
            &             + Conjg(c_QNSq_R(i,i4,i1))*c_QNSq_L(j,i3,i1) ) 
     end do 
    end do 
   end do 
  end do
  
 End Subroutine Delta_F2_Boxes


  Subroutine Delta_MK(mf_u, g, Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R, mC, U, V &
      & , mN, N, mGlu, phi_g, mS02, RS0, mP02, RP0, mSpm2, RSpm, mSup2, RSup &
      & , A_u, mu, mSdown2, RSdown, A_d, vevSM  &
      & , res )
 !---------------------------------------------------------------------------
 ! Input: mf_u, mC, mN, mGlu, mS0, mP0, mSpm
 !        U, V, N, C
 !        
 ! Output: 
 !         res
 !         
 ! written by Werner Porod, 01 Jul 2003
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_u(3), mC(:), mN(:), mGlu, mS02(:), mP02(:) &
     & , mSpm2(:), g(3), mSdown2(6), mSup2(6), RS0(:,:), RP0(:,:), vevSM(2)
  Complex(dp), Intent(in) :: U(:,:), V(:,:), N(:,:), RSpm(:,:), RSup(6,6) &
     & , RSdown(6,6), phi_g, A_u(3,3), A_d(3,3), mu
  Complex(dp), Intent(in), Dimension(3,3) :: Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R
  Complex(dp), Intent(out) :: res

  Real(dp) :: xt
  Complex(dp) :: K_VLL, K_VRR, K_LR1, K_LR2, K_SLL1, K_SLL2, K_SRR1, K_SRR2 &
    & , coupLC, coupRC
  Real(dp), Parameter :: oo4r = 0.25_dp / 0.985_dp, P1_LR = -0.71_dp   &
    & , P2_LR = 0.9_dp, P1_SLL = -0.37_dp, P2_SLL = -0.72_dp           &
    & , T3= -0.5_dp, e_d =-1._dp/3._dp
  Integer :: i1

  xt = mf_u(3)**2 / mW**2

  Call Delta_F2_Boxes(2, 1, T3, g, Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R, mf_u, mf_d_mZ &
     & , mC, U, V, mN, N, mGlu, phi_g, mSpm2, RSpm, mSup2, RSup, mSdown2   &
     & , RSdown, K_VLL, K_VRR, K_LR1, K_LR2, K_SLL1, K_SLL2, K_SRR1, K_SRR2  )
  !------------------------------------------
  ! double penguins
  !------------------------------------------
  Do i1=1,2
   Call  CoupFermionScalar31L_eff(2, 1, i1, T3, e_d, g, Y_d, Rd_L, Rd_R, RS0  &
    & , vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup, Y_u &
    & , A_u, mC, U, V, coupLC, coupRC)
   K_LR2 = K_LR2 - 16._dp * Pi2 * coupLC * coupRC / mS02(i1) 
   K_SLL1 = K_SLL1 - 8._dp * Pi2 * coupLC * coupLC / mS02(i1) 
   K_SRR1 = K_SRR1 - 8._dp * Pi2 * coupRC * coupRC / mS02(i1)
  End Do

  Call CoupFermionPseudoScalar31L_eff(2, 1, 2, T3, e_d, g, Y_d, Rd_L, Rd_R, RP0&
    & , vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup, Y_u &
    & , A_u, mC, U, V, coupLC, coupRC)
  K_LR2 = K_LR2 - 16._dp * Pi2 * coupLC * coupRC / mP02(2) 
  K_SLL1 = K_SLL1 - 8._dp * Pi2 * coupLC * coupLC / mP02(2) 
  K_SRR1 = K_SRR1 - 8._dp * Pi2 * coupRC * coupRC / mP02(2) 

  res = oo6pi2 * MBq(1) * etaB * FBhatB(1)**2                               &
      &        * ( (G_F*mW)**2 * CKM(2,1)**2 * Conjg(CKM(2,2))**2 * S0(xt)  &
      &          + oo4r * (K_VLL + K_VRR) + P1_LR * K_LR1  + P2_LR * K_LR2  &
      &          + P1_SLL * (K_SLL1 + K_SRR1) + P2_SLL * (K_SLL2 + K_SRR2)  )

 End Subroutine Delta_MK


  Subroutine Delta_MB(i, mf_u, g, Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R, mC, U, V &
      & , mN, N, mGlu, phi_g, mS02, RS0, mP02, RP0, mSpm2, RSpm, mSup2, RSup &
      & , A_u, mu, mSdown2, RSdown, A_d, vevSM  &
      & , res )
 !---------------------------------------------------------------------------
 ! Input: mf_u, mC, mN, mGlu, mS0, mP0, mSpm
 !        U, V, N, C
 !        
 ! Output: 
 !         res
 !         
 ! written by Werner Porod, 01 Jul 2003
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_u(3), mC(:), mN(:), mGlu, mS02(:), mP02(:) &
     & , mSpm2(:), g(3), mSdown2(6), mSup2(6), RS0(:,:), RP0(:,:), vevSM(2)
  Complex(dp), Intent(in) :: U(:,:), V(:,:), N(:,:), RSpm(:,:), RSup(6,6) &
     & , RSdown(6,6), phi_g, A_u(3,3), A_d(3,3), mu
  Complex(dp), Intent(in), Dimension(3,3) :: Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R
  Integer, Intent(in) :: i

  Complex(dp), Intent(out) :: res

  Real(dp) :: xt
  Complex(dp) :: B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2, B_SRR1, B_SRR2 &
    & , coupLC, coupRC
  Real(dp), Parameter :: oo4r = 0.25_dp / 0.985_dp, P1_LR = -0.71_dp   &
    & , P2_LR = 0.9_dp, P1_SLL = -0.37_dp, P2_SLL = -0.72_dp           &
    & , T3= -0.5_dp, e_d =-1._dp/3._dp
  Integer :: i1

  xt = mf_u(3)**2 / mW**2

  Call Delta_F2_Boxes(3, i, T3, g, Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R, mf_u, mf_d_mZ &
     & , mC, U, V, mN, N, mGlu, phi_g, mSpm2, RSpm, mSup2, RSup, mSdown2   &
     & , RSdown, B_VLL, B_VRR, B_LR1, B_LR2, B_SLL1, B_SLL2, B_SRR1, B_SRR2  )

  !------------------------------------------
  ! double penguins
  !------------------------------------------
  Do i1=1,2
   Call  CoupFermionScalar31L_eff(3, i, i1, T3, e_d, g, Y_d, Rd_L, Rd_R, RS0  &
    & , vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup, Y_u &
    & , A_u, mC, U, V, coupLC, coupRC)
   B_LR2 = B_LR2 - 16._dp * Pi2 * coupLC * coupRC / mS02(i1) 
   B_SLL1 = B_SLL1 - 8._dp * Pi2 * coupLC * coupLC / mS02(i1) 
   B_SRR1 = B_SRR1 - 8._dp * Pi2 * coupRC * coupRC / mS02(i1)
  End Do

  Call CoupFermionPseudoScalar31L_eff(3, i, 2, T3, e_d, g, Y_d, Rd_L, Rd_R, RP0&
    & , vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup, Y_u &
    & , A_u, mC, U, V, coupLC, coupRC)
  B_LR2 = B_LR2 - 16._dp * Pi2 * coupLC * coupRC / mP02(2) 
  B_SLL1 = B_SLL1 - 8._dp * Pi2 * coupLC * coupLC / mP02(2) 
  B_SRR1 = B_SRR1 - 8._dp * Pi2 * coupRC * coupRC / mP02(2) 

  res = oo6pi2 * MBq(i) * etaB * FBhatB(i)**2                               &
      &        * ( (G_F*mW)**2 * CKM(3,i)**2 * Conjg(CKM(3,3))**2 * S0(xt)  &
      &          + oo4r * (B_VLL + B_VRR) + P1_LR * B_LR1  + P2_LR * B_LR2  &
      &          + P1_SLL * (B_SLL1 + B_SRR1) + P2_SLL * (B_SLL2 + B_SRR2)  )

 End Subroutine Delta_MB


 Real(dp) Function DeltaRho(mZ2, mW2, mP02, RP0, mSneutrino2, RSneutrino     &
           & , mSlepton2, RSlepton, mUsquark2, RUsquark, mDSquark2, RDSquark &
           & , mC, U, V, mN, N)
 !---------------------------------------------------------------------------
 ! calculates the rho-parameter taking into account the complete SM + SUSY
 ! contributions at 1-loop, 2-loop QCD corrections due to t-quark
 ! using tree-level Higgs masses because otherwise there is a problem 
 ! with gauge invariance.
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mZ2, mW2, mP02(2), RP0(2,2), mSneutrino2(3)    &
     & , mSlepton2(6), mUsquark2(6), mDSquark2(6), mC(2), mN(4)
  Complex(dp), Intent(in) :: RSneutrino(3,3), RSlepton(6,6)  &
    & , RUsquark(6,6), RDSquark(6,6), U(2,2), V(2,2), N(4,4)

  Real(dp) :: mC2(2), mN2(4), vev2, v2(2), mHp2(2), mh2(2) &
          & , gSU2, sinW2, t2a, alpha, Rh0(2,2), mu_old, Drho_top
  Complex(dp) :: RSpm(2,2), dmZ2, dmW2

  mC2 = mC**2
  mN2 = mN**2

  gSU2 = Sqrt(4._dp * Sqrt2 * G_F * mW2)

  sinW2 = 1._dp - mW2/mZ2 

  vev2 = 2._dp * mW/gSU2
  v2(1) = -vev2 * rp0(1,1)
  v2(2) = vev2  * rp0(1,2)

  mHp2(1) = mW2
  mHp2(2) = mP02(2) + mW2
  RSpm = RP0

  mh2(1) = 0.5_dp * (mP02(2)+mZ2 - Sqrt( (mP02(2)+mZ2)**2 &
         &                   - 4._dp *mZ2*mP02(2)*(rp0(1,1)**2-rp0(1,2)**2)**2))
  mh2(2) = 0.5_dp * (mP02(2)+mZ2 + Sqrt( (mP02(2)+mZ2)**2 &
         &                   - 4._dp *mZ2*mP02(2)*(rp0(1,1)**2-rp0(1,2)**2)**2))
  t2a = (mP02(2)+mZ2)/(mP02(2)-mZ2) &
      &  *2._dp*rp0(1,1)*rp0(1,2)/(rp0(1,1)**2-rp0(1,2)**2)
  alpha = Atan(t2a)*0.5_dp

  rh0(1,1) = -Sin(alpha)
  rh0(1,2) = Cos(alpha)
  rh0(2,1) = Rh0(1,2)
  rh0(2,2) = -rh0(1,1)
  !--------------------------------------------------------------
  ! precision observable defined at the Z-pole, performing DRbar
  ! renormalization
  ! 1-loop contributions
  !--------------------------------------------------------------
  mu_old = SetRenormalizationScale( mZ2 )
  Call SetLoopMassModel(2,4,2,2,2,6,3)
  Call PiZZT1(0._dp, gSU2, sinW2, v2, mZ2, mW2, mh2, Rh0, mP02, RP0        &
   & , mHp2, RSpm, mSneutrino2, RSneutrino, mSlepton2, RSlepton, mUsquark2 &
   & , RUsquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2, mC, mC2, U, V   &
   & , mN, mN2, N, dmZ2)

  Call PiWWT1(0._dp, gSU2, sinW2, mh2, Rh0, mHp2, RSpm , v2                &
         & , mP02, RP0, mSneutrino2, RSneutrino, mSlepton2, RSlepton       &
         & , mUSquark2, RUSquark, mDSquark2, RDSquark, mf_l2, mf_u2, mf_d2 &
         & , CKM, mN, mN2, N, mC, mC2, U, V, mZ2, mW2, dmW2)

  !-----------------------------------------------------------------
  ! 1-loop top contributions, I am only interested in the SUSY part
  !-----------------------------------------------------------------
  Drho_top = 3 * G_F * mf_u2(3) * oosqrt2 * oo8pi2
  !--------------------------------------------------------------
  ! resetting renormalization scale
  !--------------------------------------------------------------
  mu_old = SetRenormalizationScale( mu_old )

  DeltaRho = dmZ2 / mZ2 - dmW2 / mW2 - Drho_top

 End Function DeltaRho


 Subroutine Gminus2a(gp,g, mf, Y_l, mSn2, mSl2, Rsl, mC, U, V, mN, N, a_mu)
 !-----------------------------------------------------------------------
 ! Calculates g-2 of a lepton using the formulas given in
 ! T.Ibrahi and P.Nath, PRD61, 095008 (2000)
 ! written by Werner Porod, 21.03.2001
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gp, g, mf, mSn2, mSl2(:), mC(:), mN(:)
  Real(dp), Intent(out) :: a_mu
  Complex(dp), Intent(in) :: Y_l, Rsl(:,:), U(:,:), V(:,:), N(:,:)

  Integer :: i1, i2, n_C, n_N
  Real(dp) :: ratio
  Complex(dp) :: coupL, coupR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Gminus2a'
       
  n_C = 2
  n_N = 4

  a_mu = 0._dp
  !-----------------------
  ! chargino contribution
  !-----------------------
  Do i1=1,n_C
   Call CoupCharginoSfermion(i1, 1, g, -0.5_dp, id2C, Y_l, ZeroC, U, V, &
                             coupL, coupR)
   ratio = mSn2 / mC(i1)**2
   a_mu = a_mu                                                        &
     &  - Real( coupL * Conjg(coupR),dp ) * F3gamma(ratio) / mC(i1)      &
     &  + 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F2(ratio) / mC(i1)**2 
  End Do
  !-------------------------
  ! neutralino contribution
  !-------------------------
  Do i1=1,n_N
   Do i2 = 1,2
    Call CoupNeutralinoSlepton(i1, i2, gp, g, RSl, Y_l, N, coupL, coupR)
    ratio = mSl2(i2) / mN(i1)**2
    a_mu = a_mu                                                      &
      &  - 2._dp * Real( coupL * Conjg(coupR),dp ) * F4(ratio) / mN(i1)    &
      &  - 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F1(ratio) / mN(i1)**2
   End Do
  End Do

  a_mu = a_mu * mf * oo16pi2

  Iname = Iname - 1

 End Subroutine Gminus2a

 Subroutine Gminus2b(i_in, mSn2, mC, cpl_CLSn_L, cpl_CLSn_R, mSl2, mN &
                   &, cpl_LNSl_L, cpl_LNSl_R, a_mu, GenerationMixing)
 !-----------------------------------------------------------------------
 ! Calculates g-2 of a lepton using the formulas given in
 ! T.Ibrahim and P.Nath, PRD61, 095008 (2000)
 ! input:
 !  i_in ................ index of the lepton: 1 -> electron, 2 -> muon,
 !                                             3 -> tau
 !  mSn2(i) ............. sneutrino masses squared
 !  mC(i) ............... chargino masses
 !  cpl_CLSn_L(i,j,k) ... left chargino - lepton -sneutrino coupling
 !  cpl_CLSn_R(i,j,k) ... right chargino - lepton -sneutrino coupling
 !  mSl2(i) ............. slepton masses squared
 !  mN(i) ............... neutralino masses squared
 !  cpl_LNSl_L(i,j,k) ... left lepton - neutralino - slepton coupling
 !  cpl_LNSl_R(i,j,k) ... right lepton - neutralino - slepton coupling
 ! output:
 !  a_mu ................ SUSY contribution to g-2
 ! written by Werner Porod, 07.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mSn2(3), mSl2(6), mC(:), mN(:)
  Real(dp), Intent(out) :: a_mu
  Complex(dp), Intent(in) :: cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:) &
                         & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)
  Logical, Intent(in) :: GenerationMixing

  Integer :: i1, i2, n_C, n_N, n_sf
  Real(dp) :: mf, ratio, mC2, mN2
  Complex(dp) :: coupL, coupR

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Gminus2b'
       
  n_C = Size(mC)
  n_N = Size(mN)

  a_mu = 0._dp
  mf = mf_l(i_in)
  !-----------------------
  ! chargino contribution
  !-----------------------
  Do i1=1,n_C
   mC2 =  mC(i1)**2
   If (GenerationMixing) Then
    Do i2 =1,3
     coupL = cpl_CLSn_L(i1,i_in,i2)
     coupR = cpl_CLSn_R(i1,i_in,i2)
     ratio = mSn2(i2) / mC2
     a_mu = a_mu                                                        &
       &  - Real( coupL * Conjg(coupR),dp ) * F3gamma(ratio) / mC(i1)      &
       &  + 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F2(ratio) / mC2 
    End Do
   Else
    coupL = cpl_CLSn_L(i1,i_in,i_in)
    coupR = cpl_CLSn_R(i1,i_in,i_in)
    ratio = mSn2(i_in) / mC2
    a_mu = a_mu                                                        &
      &  - Real( coupL * Conjg(coupR),dp ) * F3gamma(ratio) / mC(i1)      &
      &  + 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F2(ratio) / mC2 
   End If
  End Do
  !-------------------------
  ! neutralino contribution
  !-------------------------
  Do i1=1,n_N
   mN2 = mN(i1)**2
   If (GenerationMixing) Then
    Do i2 = 1,6
     coupL = cpl_LNSl_L(i_in,i1,i2)
     coupR = cpl_LNSl_R(i_in,i1,i2)
     ratio = mSl2(i2) / mN2
     a_mu = a_mu                                                      &
       &  - 2._dp * Real( coupL * Conjg(coupR),dp ) * F4(ratio) / mN(i1)    &
       &  - 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F1(ratio) / mN2
    End Do
   Else
    Do i2 = 1,2
     n_sf = (i_in-1)*2 + i2
     coupL = cpl_LNSl_L(i_in,i1,n_sf)
     coupR = cpl_LNSl_R(i_in,i1,n_sf)
     ratio = mSl2(n_sf) / mN2
     a_mu = a_mu                                                      &
       &  - 2._dp * Real( coupL * Conjg(coupR),dp ) * F4(ratio) / mN(i1)    &
       &  - 2._dp * mf * (Abs(coupL)**2 + Abs(coupR)**2) * F1(ratio) / mN2
    End Do
   End If
  End Do

  a_mu = a_mu * mf * oo16pi2

  Iname = Iname - 1

 End Subroutine Gminus2b
 Subroutine EDMinit
 !----------------------------------------------------------------------
 ! input for the program:
 ! - calculates several numerical constants, data are stored in 
 !   common /EDMnumerical/
 ! written by Thomas Gajdosik, 22.02.2000
 ! - 31.07.01: portation to f90 by Werner Porod 
 !----------------------------------------------------------------------
 Implicit None

 !-------------
 ! EDM data:
 !-------------
  Open(11,file='EDMs.in',status='old',err=200)

  Read (11,800) ecmfactor
  Read (11,800) elimit
  Read (11,800) nlimit
  Read (11,800) qcdCorrection
  Read (11,800) CqcdCorrection
  Read (11,800) chiralMass
  Read (11,800) deltaUp
  Read (11,800) deltaDown
  Read (11,800) deltaStrange
  
  Close(11)
200  Return

 800  Format(f16.7)

 End Subroutine EDMinit

 Subroutine Lepton_EDM(i_in, mN, mSl2, cpl_LNSl_L, cpl_LNSl_R    &
                     & , mC, mSnu2, cpl_CLSn_L, cpl_CLSn_R, EDM  )
 !----------------------------------------------------------------------
 ! calculates the EDM-contributions to leptons in the MSSM
 ! input:
 ! output: EDM(i) .... i=1 Neutralino contribution
 !                     i=2 Chargino contribution
 !                     i=3 total EDM
 ! written by Werner Porod, 31.07.01
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mN(4), mSl2(6), mC(2), mSnu2(3)
  Complex(dp), Intent(in) :: cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:) &
       & , cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:)
  Real(dp), Intent(out) :: EDM(3)

  Integer :: i1, i2, i_sl
  Real(dp) :: mF2, r

  EDM = 0._dp
  !------------------------------------
  ! neutralino contribution
  !------------------------------------
  Do i1=1,4
   mF2 = mN(i1)**2
   Do i2=1,2
    i_sl = 2*(i_in-1)+i2
    r = mF2 / mSl2(i_sl)  
    EDM(1) = EDM(1)                                                   &
   &     - Aimag(cpl_LNSl_L(i_in,i1,i_sl) * Conjg( cpl_LNSl_R(i_in,i1,i_sl)) )&
   &     * FeynFunctionB(r) * mN(i1) / mSl2(i_sl)  
    L_edm_contri(i1,i2) =   - Aimag(cpl_LNSl_L(i_in,i1,i_sl) * Conjg( cpl_LNSl_R(i_in,i1,i_sl)) )&
   &     * FeynFunctionB(r) * mN(i1) / mSl2(i_sl)  
   Enddo
  Enddo

  !------------------------------------
  ! chargino contribution
  !------------------------------------
  Do i1=1,2
   mF2 = mC(i1)**2
   r = mF2 / mSnu2(i_in) 
   EDM(2) = EDM(2)                  &
   &     + Aimag(cpl_CLSn_R(i1,i_in,i_in) * Conjg( cpl_CLSn_L(i1,i_in,i_in)) )&
   &     * FeynFunctionA(r) * mC(i1) / mSnu2(i_in)
   l_edm_contri(5,i1) = Aimag(cpl_CLSn_R(i1,i_in,i_in) * Conjg( cpl_CLSn_L(i1,i_in,i_in)) )&
   &     * FeynFunctionA(r) * mC(i1) / mSnu2(i_in)
  Enddo

  EDM(3) = EDM(1) + EDM(2)
  EDM = oo16pi2 * ecmfactor * EDM 
  l_edm_contri =  oo16pi2 * ecmfactor *l_edm_contri 

 End Subroutine Lepton_EDM

 Subroutine Lepton_EDM3(i_in, mN, mSl2, cpl_LNSl_L, cpl_LNSl_R    &
                     & , mC, mSnu2, cpl_CLSn_L, cpl_CLSn_R, EDM  )
 !----------------------------------------------------------------------
 ! calculates the EDM-contributions to leptons in the MSSM
 ! input:
 ! output: EDM(i) .... i=1 Neutralino contribution
 !                     i=2 Chargino contribution
 !                     i=3 total EDM
 ! written by Werner Porod, 31.07.01
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: i_in
  Real(dp), Intent(in) :: mN(4), mSl2(6), mC(2), mSnu2(3)
  Complex(dp), Intent(in) :: cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:) &
       & , cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:)
  Real(dp), Intent(out) :: EDM(3)

  Integer :: i1, i2
  Real(dp) :: mF2, r

  EDM = 0._dp
  !------------------------------------
  ! neutralino contribution
  !------------------------------------
  Do i1=1,4
   mF2 = mN(i1)**2
   Do i2=1,6
    r = mF2 / mSl2(i2)  
    EDM(1) = EDM(1)                                                   &
   &     - Aimag(cpl_LNSl_L(i_in,i1,i2) * Conjg( cpl_LNSl_R(i_in,i1,i2)) )&
   &     * FeynFunctionB(r) * mN(i1) / mSl2(i2)  
   Enddo
  Enddo

  !------------------------------------
  ! chargino contribution
  !------------------------------------
  Do i1=1,2
   mF2 = mC(i1)**2
   Do i2=1,3
    r = mF2 / mSnu2(i2) 
    EDM(2) = EDM(2)                  &
    &     + Aimag(cpl_CLSn_R(i1,i_in,i2) * Conjg( cpl_CLSn_L(i1,i_in,i2)) )&
    &     * FeynFunctionA(r) * mC(i1) / mSnu2(i2)
   End Do
  Enddo

  EDM(3) = EDM(1) + EDM(2)
  EDM = oo16pi2 * ecmfactor * EDM 

 End Subroutine Lepton_EDM3

 Subroutine Neutron_EDM(mN, mC, mglu, mSup2, c_UNSu_L, c_UNSu_R              &
            & , c_GUSu_L, c_GUSu_R, c_CDSu_L, c_CDSu_R, mSdown2              &
            & , c_DNSd_L, c_DNSd_R, c_GDSd_L, c_GDSd_R, c_CUSd_L, c_CUSd_R   &
            & , GenerationMixing, L_cq, L_rqp, EDMs, EDMS_parts, edm_parts2)
 !-------------------------------------------------------------------
 ! calculates the neutron EDM
 ! written by Werner Porod
 ! 06.07.03: adding charged Higgs boson
 !           adding optional argument EDMS_parts(2,i=1,..,4)
 !           for the contributions: i=1...neutralino
 !                                  i=2...chargino
 !                                  i=3...gluino
 !                                  i=4...charged Higgs
 !-------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mN(4), mC(2), mglu, mSup2(6), mSdown2(6)
  Complex(dp), Intent(in) :: c_UNSu_L(3,4,6), c_UNSu_R(3,4,6), c_GUSu_L(3,6) &
       & , c_GUSu_R(3,6), c_CDSu_L(2,3,6), c_CDSu_R(2,3,6), c_DNSd_L(3,4,6)  &
       & , c_DNSd_R(3,4,6), c_GDSd_L(3,6), c_GDSd_R(3,6), c_CUSd_L(2,3,6)    &
       & , c_CUSd_R(2,3,6)
  Logical, Intent(in) :: GenerationMixing, L_cq, L_rqp
  Real(dp), Intent(out) :: EDMs(2)
  Real(dp), Optional, Intent(out) :: EDMS_parts(2,3), edm_parts2(2,6)

  Integer :: i1, i2
  Real(dp) :: r, B_r, A_r, mF2, e_d, e_u, g_s, d_d_n, dh_d_n, d_u_n, dh_u_n &
    & , d_s_n, d_d_c, dh_d_c, d_u_c, dh_u_c, d_s_c, d_d_g, dh_d_g, d_u_g    &
    & , dh_u_g, d_s_g

  !----------------------------
  ! Initialization
  !----------------------------
  EDMs = 0._dp
  d_d_n = 0._dp  ! neutralino contribution to d-quark EDM
  dh_d_n = 0._dp ! neutralino contribution to d-quark chromoelectric EDM
  d_u_n = 0._dp  ! neutralino contribution to u-quark EDM
  dh_u_n = 0._dp ! neutralino contribution to u-quark chromoelectric EDM
  d_s_n = 0._dp  ! neutralino contribution to s-quark EDM

  d_d_c = 0._dp  ! chargino contribution to d-quark EDM
  dh_d_c = 0._dp ! chargino contribution to d-quark chromoelectric EDM
  d_u_c = 0._dp  ! chargino contribution to u-quark EDM
  dh_u_c = 0._dp ! chargino contribution to u-quark chromoelectric EDM
  d_s_c = 0._dp  ! chargino contribution to s-quark EDM

  d_d_g = 0._dp  ! gluino contribution to d-quark EDM
  dh_d_g = 0._dp ! gluino contribution to d-quark chromoelectric EDM
  d_u_g = 0._dp  ! gluino contribution to u-quark EDM
  dh_u_g = 0._dp ! gluino contribution to u-quark chromoelectric EDM
  d_s_g = 0._dp  ! gluino contribution to s-quark EDM

  e_d = -1._dp / 3._dp
  e_u = 2._dp / 3._dp
  g_s = Sqrt(4._dp * Pi * AlphaS_mZ)

  !------------------------------------
  ! neutralino contribution
  !------------------------------------
  If (GenerationMixing) Then
   Do i1=1,4
    mF2 = mN(i1)**2
    Do i2=1,6
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSdown2(i2)
     d_d_n = d_d_n - Aimag(c_DNSd_L(1,i1,i2) * Conjg(c_DNSd_R(1,i1,i2)) ) * B_r
     If (L_rqp) d_s_n = d_s_n   &
          & - Aimag(c_DNSd_L(2,i1,i2) * Conjg(c_DNSd_R(2,i1,i2)) ) * B_r

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSup2(i2)
     d_u_n = d_u_n - Aimag(c_UNSu_L(1,i1,i2) * Conjg(c_UNSu_R(1,i1,i2)) ) * B_r
    Enddo
   Enddo
 
  Else ! .not.GenerationMixing
   Do i1=1,4
    mF2 = mN(i1)**2
    Do i2=1,2
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSdown2(i2)
     d_d_n = d_d_n - Aimag(c_DNSd_L(1,i1,i2) * Conjg(c_DNSd_R(1,i1,i2)) ) * B_r

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mN(i1) / mSup2(i2)
     d_u_n = d_u_n - Aimag(c_UNSu_L(1,i1,i2) * Conjg(c_UNSu_R(1,i1,i2)) ) * B_r
     If (L_rqp) Then
      r = mF2 / mSdown2(i2+2)
      B_r = FeynFunctionB(r) * mN(i1) / mSdown2(i2+2)
      d_s_n = d_s_n &
          & - Aimag(c_DNSd_L(2,i1,i2+2) * Conjg(c_DNSd_R(2,i1,i2+2)) ) * B_r
     End If
    Enddo
   Enddo
  End If

  If (L_cq) Then
   dh_d_n = g_s * d_d_n
   dh_u_n = g_s * d_u_n
  Endif
  d_d_n = e_d * d_d_n
  d_s_n = e_d * d_s_n
  d_u_n = e_u * d_u_n
  !------------------------------------
  ! chargino contribution
  !------------------------------------
  If (GenerationMixing) Then
   Do i1=1,2
    mF2 = mC(i1)**2
    Do i2=1,6
     r = mF2 / mSup2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSup2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSup2(i2)
     d_d_c = d_d_c - Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) &
           &         * (e_u * B_r + (e_d - e_u) * A_r)
     If (L_cq) dh_d_c = dh_d_c                        &
          & + g_s * Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) * B_r
     If (L_rqp) d_s_c = d_s_c   &
          & - Aimag(c_CDSu_L(i1,2,i2) * Conjg(c_CDSu_R(i1,2,i2)) ) &
          &         * (e_u * B_r + (e_d - e_u) * A_r)

     r = mF2 / mSdown2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSdown2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSdown2(i2)
     d_u_c = d_u_c - Aimag(c_CUSd_L(i1,1,i2) * Conjg(c_CUSd_R(i1,1,i2)) ) &
          &         * (e_d * B_r + (e_u - e_d) * A_r)
     If (L_cq) dh_u_c = dh_u_c                        &
         & + g_s * Aimag(c_CUSd_L(i1,1,i2) * Conjg(c_CUSd_R(i1,1,i2)) ) * B_r
    Enddo
   Enddo
 
  Else ! .not.GenerationMixing
   Do i1=1,2
    mF2 = mC(i1)**2
    Do i2=1,2
     r = mF2 / mSup2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSup2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSup2(i2)
     d_d_c = d_d_c - Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) &
           &         * (e_u * B_r + (e_d - e_u) * A_r)
     If (L_cq) dh_d_c = dh_d_c                        &
          & + g_s * Aimag(c_CDSu_L(i1,1,i2) * Conjg(c_CDSu_R(i1,1,i2)) ) * B_r

     r = mF2 / mSdown2(i2)
     A_r = FeynFunctionA(r) * mC(i1) / mSdown2(i2)
     B_r = FeynFunctionB(r) * mC(i1) / mSdown2(i2)
     d_u_c = d_u_c - Aimag(c_CUSd_L(1,i1,i2) * Conjg(c_CUSd_R(1,i1,i2)) ) &
          &         * (e_d * B_r + (e_u - e_d) * A_r)
     If (L_cq) dh_u_c = dh_u_c                        &
         & + g_s * Aimag(c_CUSd_L(i1,1,i2) * Conjg(c_CUSd_R(i1,1,i2)) ) * B_r

     If (L_rqp) Then
      r = mF2 / mSup2(i2+2)
      A_r = FeynFunctionA(r) * mC(i1) / mSup2(i2+2)
      B_r = FeynFunctionB(r) * mC(i1) / mSup2(i2+2)
      d_s_c = d_s_c &
          & - Aimag(c_CDSu_L(i1,2,i2+2) * Conjg(c_CDSu_R(i1,2,i2+2)) ) &
          &         * (e_u * B_r + (e_d - e_u) * A_r)
     End If
    Enddo
   Enddo
  End If
  !------------------------------------
  ! gluino contribution
  !------------------------------------
  mF2 = mGlu**2
  If (GenerationMixing) Then
    Do i2=1,6
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSdown2(i2)
     d_d_g = d_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSdown2(i2) - B_r/3._dp
      dh_d_g = dh_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * A_r
     end if
     If (L_rqp) d_s_g = d_s_g   &
          & - Aimag(c_GDSd_L(2,i2) * Conjg(c_GDSd_R(2,i2)) ) * B_r

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSup2(i2)
     d_u_g = d_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSup2(i2) - B_r/3._dp
      dh_u_g = dh_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * A_r
     end if
    Enddo
 
  Else ! .not.GenerationMixing
    Do i2=1,2
     r = mF2 / mSdown2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSdown2(i2)
     d_d_g = d_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSdown2(i2) - B_r/3._dp
      dh_d_g = dh_d_g - Aimag(c_GDSd_L(1,i2) * Conjg(c_GDSd_R(1,i2)) ) * A_r
     end if

     r = mF2 / mSup2(i2)
     B_r = FeynFunctionB(r) * mGlu / mSup2(i2)
     d_u_g = d_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * B_r
     if (L_cq) then
      A_r = 3._dp * FeynFunctionA(r) * mGlu / mSup2(i2) - B_r/3._dp
      dh_u_g = dh_u_g - Aimag(c_GUSu_L(1,i2) * Conjg(c_GUSu_R(1,i2)) ) * A_r
     end if
     If (L_rqp) Then
      r = mF2 / mSdown2(i2+2)
      B_r = FeynFunctionB(r) * mGlu / mSdown2(i2+2)
      d_s_g = d_s_g &
          & - Aimag(c_GDSd_L(2,i2+2) * Conjg(c_GDSd_R(2,i2+2)) ) * B_r
     End If
    Enddo
  End If

  d_d_g = 16._dp * e_d * d_d_g / 3._dp
  d_s_g = 16._dp * e_d * d_s_g / 3._dp
  d_u_g = 16._dp * e_u * d_u_g / 3._dp

  !---------------------------------------
  ! neutron EDM in the chiral quark model
  !---------------------------------------
  If (L_cq) Then
   EDMs(1) = qcdCorrection * ( 4._dp * (d_d_n + d_d_c + d_d_g)   &
         &                   - d_u_n - d_u_c -d_u_g) / 3._dp     &
         & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)             &
         &              * ( 4._dp * (dh_d_n + dh_d_c + dh_d_g)   &
         &                - dh_u_n - dh_u_c -dh_u_g) / 3._dp 
if (present(edm_parts2)) then
   edm_parts2(1,1) = d_d_n
   edm_parts2(1,2) = dh_d_n
   edm_parts2(1,3) = d_d_c
   edm_parts2(1,4) = dh_d_c
   edm_parts2(1,5) = d_d_g
   edm_parts2(1,6) = dh_d_g
   edm_parts2(2,1) = d_u_n
   edm_parts2(2,2) = dh_u_n
   edm_parts2(2,3) = d_u_c
   edm_parts2(2,4) = dh_u_c
   edm_parts2(2,5) = d_u_g
   edm_parts2(2,6) = dh_u_g
end if
  End If
  !----------------------------------------------------
  ! neutron EDM in the relativistic quark-parton model
  ! Deltas refer to proton, not to neutron
  !----------------------------------------------------
  If (L_rqp) Then
   EDMs(2) = qcdCorrection * ( DeltaUp      * (d_d_n + d_d_c + d_d_g) &
         &                   + DeltaDown    * (d_u_n + d_u_c + d_u_g) &
         &                   + DeltaStrange * (d_s_n + d_s_c + d_s_g) )
  End If

  EDMs = oo16pi2 * ecmfactor * EDMs

  !-----------------------------------
  ! individual contributions
  !-----------------------------------
  If (present(EDMS_parts)) then
   !  chiral quark model
   EDMs_parts(1,1) = qcdCorrection * ( 4._dp * d_d_n  - d_u_n ) / 3._dp     &
                 & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)                &
         &              * ( 4._dp * dh_d_n - dh_u_n ) / 3._dp
   EDMs_parts(1,2) = qcdCorrection * ( 4._dp * d_d_c  - d_u_c ) / 3._dp     &
                 & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)                &
         &              * ( 4._dp * dh_d_c - dh_u_c ) / 3._dp
   EDMs_parts(1,3) = qcdCorrection * ( 4._dp * d_d_g  - d_u_g ) / 3._dp     &
                 & + CqcdCorrection * Sqrt(Alpha_mZ * oo4pi)                &
         &              * ( 4._dp * dh_d_g - dh_u_g ) / 3._dp
 
   ! relativistic quark-parton model
   EDMs_parts(2,1) = qcdCorrection * ( DeltaDown * d_d_n + DeltaUp * d_u_n  &
                   &                 + DeltaStrange * d_s_n  )
   EDMs_parts(2,2) = qcdCorrection * ( DeltaDown * d_d_c + DeltaUp * d_u_c  &
                   &                 + DeltaStrange * d_s_c  )
   EDMs_parts(2,3) = qcdCorrection * ( DeltaDown * d_d_g + DeltaUp * d_u_g  &
                   &                 + DeltaStrange * d_s_g  )
   
   EDMs_parts = oo16pi2 * ecmfactor * EDMs_parts
  end if

 End Subroutine Neutron_EDM


! Interface LtoLpGamma
 Subroutine LtoLpGammaMSSM(i, j, mSn2, mC, cpl_CLSn_L, cpl_CLSn_R, mSl2, mN &
                   &, cpl_LNSl_L, cpl_LNSl_R, width, Br)
 !-----------------------------------------------------------------------
 ! Calculates the width l -> l' gamma using the formula of
 ! J.Hisano et al., PRD 53, 2442 (1996)
 ! input:
 !  i ................... index of decaying lepton: 2 -> muon, 3 -> tau
 !  j ................... index of final state lepton: 1 -> electron, 2 -> muon
 !  mSn2(i) ............. sneutrino masses squared
 !  mC(i) ............... chargino masses
 !  cpl_CLSn_L(i,j,k) ... left chargino - lepton -sneutrino coupling
 !  cpl_CLSn_R(i,j,k) ... right chargino - lepton -sneutrino coupling
 !  mSl2(i) ............. slepton masses squared
 !  mN(i) ............... neutralino masses squared
 !  cpl_LNSl_L(i,j,k) ... left lepton - neutralino - slepton coupling
 !  cpl_LNSl_R(i,j,k) ... right lepton - neutralino - slepton coupling
 ! output:
 !  width ............... decay width in GeV
 !  Br .................. branching ratio, optional
 ! written by Werner Porod, 09.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mSn2(3), mSl2(6), mC(:), mN(:)
  Real(dp), Intent(out) :: width
  Real(dp), Intent(out), Optional :: Br
  Complex(dp), Intent(in) :: cpl_CLSn_L(:,:,:), cpl_CLSn_R(:,:,:) &
                         & , cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)

  Integer :: i1, i2, n_C, n_N
  Real(dp) :: ratio, mN2, mC2, part(2)
  Complex(dp) :: AL, AR ,aln ,arn, alI, arI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'LtoLpGammaMSSM'
       
  n_C = Size(mC)
  n_N = Size(mN)

  AL = ZeroC
  AR = ZeroC
  !-------------------------------------
  ! neutralino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_N
   mN2 = mN(i1)**2
   Do i2 = 1, 6
    ratio = mSl2(i2) / mN2
    part(1) = 2._dp * mf_l(i) * F1(ratio) / mN2
    part(2) = 2._dp * F4(ratio) / mN(i1)
    ALI = cpl_LNSl_L(j,i1,i2) * ( Conjg( cpl_LNSl_L(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNSl_R(i,i1,i2) ) * part(2) )
    ARI = cpl_LNSl_R(j,i1,i2) * ( Conjg( cpl_LNSl_R(i,i1,i2) ) * part(1) &
        &                       + Conjg( cpl_LNSl_L(i,i1,i2) ) * part(2) )
    AL = AL + ALI
    AR = AR + ARI
    
   End Do
  End Do
  aln = oo32pi2 * al
  arn = oo32pi2 * ar
  !-------------------------------------
  ! chargino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_C
   mC2 = mC(i1)**2
   Do i2 = 1, 3
    ratio = mSn2(i2) / mC2
    part(1) = - 2._dp * mf_l(i) * F2(ratio) / mC2
    part(2) = F3gamma(ratio) / mC(i1)
    ALI = cpl_CLSn_L(i1,j,i2) * ( Conjg( cpl_CLSn_L(i1,i,i2) ) * part(1) &
        &                       + Conjg( cpl_CLSn_R(i1,i,i2) ) * part(2) )
    ARI = cpl_CLSn_R(i1,j,i2) * ( Conjg( cpl_CLSn_R(i1,i,i2) ) * part(1) &
        &                       + Conjg( cpl_CLSn_L(i1,i,i2) ) * part(2) )
    AL = AL + ALI
    AR = AR + ARI
   End Do
  End Do
  AL = oo32pi2 * AL
  AR = oo32pi2 * AR

  width = 0.25_dp * Alpha * mf_l(i)**3 * (Abs(AL)**2 + Abs(AR)**2)

  If (Present(Br)) Then
   If (i.Eq.2) Br = width / GammaMu
   If (i.Eq.3) Br = width / GammaTau
  End If

  Iname = Iname - 1

 End Subroutine LtoLpGammaMSSM

 Subroutine LtoLpGammaEps3(i, j, mC, mP02, cpl_CCP0_L, cpl_CCP0_R   &
                         &, mS02, cpl_CCS0_L, cpl_CCS0_R            &
                         &, mN, mSpm2, cpl_SpmCN_L, cpl_SpmCN_R     &
                         &, mUsquark2, cpl_CDSu_L, cpl_CDSu_R       &
                         &, mDsquark2, cpl_CUSd_L, cpl_CUSd_R       &
                         &, width, Br, GenerationMixing, Amplitude)
 !-----------------------------------------------------------------------
 ! Calculates the width l -> l' gamma based on the formula of
 ! J.Hisano et al., PRD 53, 2442 (1996), this is an extension to the
 ! 3-generation R-parity model
 ! input:
 !  i ................... index of decaying lepton: 2 -> muon, 3 -> tau
 !  j ................... index of final state lepton: 1 -> electron, 2 -> muon
 !  mC(i) ............... chargino masses
 !  mP02(i) ............. pseudoscalar masses squared
 !  cpl_CCP0_L(i,j,k) ... left chargino - chargino - pseudoscalar coupling
 !  cpl_CCP0_R(i,j,k) ... right chargino - chargino - pseudoscalar coupling
 !  mS02(i) ............. scalar masses squared
 !  cpl_CCS0_L(i,j,k) ... left chargino - chargino - scalar coupling
 !  cpl_CCS0_R(i,j,k) ... right chargino - chargino - scalar coupling
 !  mN(i) ............... neutralino masses squared
 !  mSpm2(i) ............. charged scalar masses squared
 !  cpl_SpmCN_L(i,j,k) .. left charged scalar - chargino - neutralino  coupling
 ! cpl_SpmCN_R(i,j,k) .. right charged scalar - chargino - neutralino  coupling
 !  GenerationMixing ... generation mixing between quarks/squarks is taken into
 !                       account if .TRUE., optional
 !  Amplitude(i,j) ..... contains the neutralino the chargino,, the u-quark
 !                       and the d-quark contribution to the amplitude if
 !                       present, optional, the second index specifies L,R
 ! output:
 !  width ............... decay width in GeV
 !  Br .................. branching ratio, optional
 ! written by Werner Porod, 14.05.2001
 !-----------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Real(dp), Intent(in) :: mP02(:), mS02(:), mSpm2(:), mC(:), mN(:)   &
         &, mUsquark2(6), mDsquark2(6)
  Real(dp), Intent(out) :: width
  Real(dp), Intent(out), Optional :: Br
  Complex(dp), Intent(in) :: cpl_CCP0_L(:,:,:), cpl_CCP0_R(:,:,:)    &
                         & , cpl_CCS0_L(:,:,:), cpl_CCS0_R(:,:,:)    &
                         & , cpl_SpmCN_L(:,:,:), cpl_SpmCN_R(:,:,:)  &
                         & , cpl_CDSu_L(:,:,:), cpl_CDSu_R(:,:,:)    &
                         & , cpl_CUSd_L(:,:,:), cpl_CUSd_R(:,:,:)
  Complex(dp), Intent(out), Optional :: Amplitude(4,2)
  Logical, Intent(in), Optional :: GenerationMixing

  Integer :: i1, i2, n_C, n_N, n_P0, n_S0, n_Spm, i_gen
  Real(dp) :: ratio, mN2, mC2, part(2), mF2
  Complex(dp) :: AL, AR, amp(4,2)
  Logical :: generation

  Iname = Iname + 1
  NameOfUnit(Iname) = 'LtoLpGammaEps3'
       
  n_C = Size(mC)
  n_N = Size(mN)
  n_P0 = Size(mP02)
  n_S0 = Size(mS02)
  n_Spm = Size(mSpm2)

  AL = ZeroC
  AR = ZeroC
  amp = ZeroC

  If (Present(GenerationMixing) ) Then
   generation = GenerationMixing
  Else
   generation = .False.
  End If

  !-------------------------------------
  ! neutralino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_N
   mN2 = mN(i1)**2
   Do i2 = 1, n_Spm
    ratio = mSpm2(i2) / mN2
    part(1) = 2._dp * mf_l(i) * F1(ratio) / mN2
    part(2) = 2._dp * F4(ratio) / mN(i1)
    AL = AL + cpl_SpmCN_L(i2,j,i1) * ( Conjg(cpl_SpmCN_L(i2,i,i1) ) * part(1) &
       &                            + Conjg( cpl_SpmCN_R(i2,i,i1) ) * part(2) )
    AR = AR + cpl_SpmCN_R(i2,j,i1) * ( Conjg(cpl_SpmCN_R(i2,i,i1) ) * part(1) &
       &                            + Conjg( cpl_SpmCN_L(i2,i,i1) ) * part(2) )
   End Do
  End Do
  amp(1,1) = AL
  amp(1,2) = AR
  AL = ZeroC
  AR = ZeroC
  !-------------------------------------
  ! chargino contribution to A_L, A_R
  !-------------------------------------
  Do i1 = 1, n_C
   mC2 = mC(i1)**2
   Do i2 = 1, n_P0
    ratio = mP02(i2) / mC2
    part(1) = - 2._dp * mf_l(i) * F2(ratio) / mC2
    part(2) = F3gamma(ratio) / mC(i1)
    AL = AL + cpl_CCP0_L(i1,j,i2) * ( Conjg( cpl_CCP0_L(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCP0_R(i1,i,i2) ) * part(2) )
    AR = AR + cpl_CCP0_R(i1,j,i2) * ( Conjg( cpl_CCP0_R(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCP0_L(i1,i,i2) ) * part(2) )
   End Do
   Do i2 = 1, n_S0
    ratio = mS02(i2) / mC2
    part(1) = - 2._dp * mf_l(i) * F2(ratio) / mC2
    part(2) = F3gamma(ratio) / mC(i1)
    AL = AL + cpl_CCS0_L(i1,j,i2) * ( Conjg( cpl_CCS0_L(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCS0_R(i1,i,i2) ) * part(2) )
    AR = AR + cpl_CCS0_R(i1,j,i2) * ( Conjg( cpl_CCS0_R(i1,i,i2) ) * part(1) &
       &                            + Conjg( cpl_CCS0_L(i1,i,i2) ) * part(2) )
   End Do
  End Do
  amp(2,1) = AL !- amp(1,1)
  amp(2,2) = AR !- amp(1,2)
  AL = ZeroC
  AR = ZeroC
  !----------------------------------------------------------------------
  ! d-quark contribution to A_L, A_R
  ! note that the colour factor 3 cancels the 1/3 contribution of the
  ! quark charges
  !----------------------------------------------------------------------
  If (generation) Then
   Do i1 = 1, 3
    mF2 = mf_d(i1)**2
    Do i2 = 1,6
     ratio = mUsquark2(i2) / mF2
     part(1) = mf_l(i) * ( -2._dp * F2(ratio) - 4._dp * F1(ratio) ) / mF2
     part(2) = (F3gamma(ratio) + 4._dp * F4(ratio) ) / mf_d(i1)
     AL = AL + cpl_CDSu_L(j,i1,i2) * ( Conjg( cpl_CDSu_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CDSu_R(j,i1,i2) * ( Conjg( cpl_CDSu_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  Else
   Do i1 = 1, 3
    mF2 = mf_d(i1)**2
    i_gen = 2*i1-1
    Do i2 = i_gen,i_gen+1
     ratio = mUsquark2(i2) / mF2
     part(1) = mf_l(i) * (- 2._dp * F2(ratio) - 4._dp * F1(ratio) ) / mF2
     part(2) = (F3gamma(ratio) + 4._dp * F4(ratio) ) / mf_d(i1)
     AL = AL + cpl_CDSu_L(j,i1,i2) * ( Conjg( cpl_CDSu_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CDSu_R(j,i1,i2) * ( Conjg( cpl_CDSu_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CDSu_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  End If
  amp(3,1) = AL !- amp(1,1)
  amp(3,2) = AR !- amp(1,2)
  AL = ZeroC
  AR = ZeroC
  !----------------------------------------------------------------------
  ! u-quark contribution to A_L, A_R
  ! note that the colour factor 3 cancels the 1/3 contribution of the
  ! quark charges
  !----------------------------------------------------------------------
  If (generation) Then
   Do i1 = 1, 3
    mF2 = mf_u(i1)**2
    Do i2 = 1,6
     ratio = mDsquark2(i2) / mF2
     part(1) = mf_l(i) * (4._dp * F2(ratio) + 2._dp * F1(ratio) ) / mF2
     part(2) = (- 2._dp * F3gamma(ratio) - 2._dp * F4(ratio) ) / mf_u(i1)
     AL = AL + cpl_CUSd_L(j,i1,i2) * ( Conjg( cpl_CUSd_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CUSd_R(j,i1,i2) * ( Conjg( cpl_CUSd_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  Else
   Do i1 = 1, 3
    mF2 = mf_u(i1)**2
    i_gen = 2*i1-1
    Do i2 = i_gen,i_gen+1
     ratio = mDsquark2(i2) / mF2
     part(1) = mf_l(i) * (4._dp * F2(ratio) + 2._dp * F4(ratio) ) / mF2
     part(2) = (- 2._dp * F3gamma(ratio) - 2._dp * F4(ratio) ) / mf_u(i1)
     AL = AL + cpl_CUSd_L(j,i1,i2) * ( Conjg( cpl_CUSd_L(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_R(i,i1,i2) ) * part(2) )
     AR = AR + cpl_CUSd_R(j,i1,i2) * ( Conjg( cpl_CUSd_R(i,i1,i2) ) * part(1) &
        &                            + Conjg( cpl_CUSd_L(i,i1,i2) ) * part(2) )
    End Do
   End Do
  End If
  amp(4,1) = AL 
  amp(4,2) = AR 

  AL = oo32pi2 * Sum(amp(:,1))
  AR = oo32pi2 * Sum(amp(:,2))

  width = 0.25_dp * Alpha * mf_l(i)**3 * (Abs(AL)**2 + Abs(AR)**2)

  If (Present(Br)) Then
   If (i.Eq.2) Br = width / GammaMu
   If (i.Eq.3) Br = width / GammaTau
  End If

  If (Present(Amplitude))  Amplitude = amp

  Iname = Iname - 1
 
 End Subroutine LtoLpGammaEps3

 Subroutine MuEconversion(gp, g, Y_l, mSlep2, RSl, mN, N, mSnu2, RSn &
                      & , mSup2, RSu, Y_u, uU_L, uU_R, mSdown2, RSd  &
                      & , Y_d, uD_L, uD_R, mC, U, V, Znuc, Nnuc, Br)
 !-------------------------------------------------------------------------
 ! calculates mu-e conversion on heavy nuclei
 ! Formulas by J.Hisano et al., PRD 55 91996) 2442 and
 !             R.Kitano et al., PRD 66 (2002) 096002
 ! written by Werner Porod, 31.08.2006
 !-------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gp, g             ! U(1) and SU(2) gauge couplings
  Complex(dp), Intent(in) :: Y_L(3,3)       ! lepton Yukawa couplings
  Real(dp), Intent(in) :: mSlep2(6)         ! slepton masses squared
  Complex(dp), Intent(in) :: RSl(6,6)       ! slepton mixing matrix
  Real(dp), Intent(in) :: mN(4)             ! neutralino masses
  Complex(dp), Intent(in) :: N(4,4)         ! neutralino mixing matrix 
  Real(dp), Intent(in) :: mSnu2(3)          ! sneutrino masses squared
  Complex(dp), Intent(in) :: RSn(3,3)       ! sneutrino mixing matrix
  Real(dp), Intent(in) :: mC(2)             ! chargino masses
  Complex(dp), Intent(in) :: U(2,2), V(2,2) ! chargino mixing matrices
  Real(dp), Intent(in) :: mSup2(6)          ! u-squark masses squared
  Complex(dp), Intent(in) :: RSu(6,6)       ! u-squark mixing matrix
  Complex(dp), Intent(in) :: Y_u(3,3)       ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: uU_L(3,3), uU_R(3,3) ! u-quark mixing matrices
  Real(dp), Intent(in) :: mSdown2(6)        ! d-squark masses squared
  Complex(dp), Intent(in) :: RSd(6,6)       ! d-squark mixing matrix
  Complex(dp), Intent(in) :: Y_d(3,3)       ! d-quark Yukawa couplings
  Complex(dp), Intent(in) :: uD_L(3,3), uD_R(3,3) ! d-quark mixing matrices
  Integer, Intent(in) :: Znuc, Nnuc   ! proton and neutron number of the nuclei
  Real(dp), Intent(out) :: BR
  
  Integer :: i1, i2, i3, i4
  Real(dp) :: mN2(4), mC2(2), fun1, fun2, X_ax, cosW, cosW2, sinW2 &
      & , c_Zee_L, c_Zee_R, fun3, fun4
  Complex(dp) :: c_LNSl_L(3,4,6), c_LNSl_R(3,4,6), QSlepton(6,6), QSnu(3,3) &
            &  , c_CLSn_L(2,3,3), c_CLSn_R(2,3,3), ELn(4,4), ERn(4,4)       &
            &  , ELc(2,2), ERc(2,2), c_DNSd_L(4,6), c_DNSd_R(4,6)           &
            &  , c_CDSu_L(2,6), c_CDSu_R(2,6), c_UNSu_L(4,6), c_UNSu_R(4,6) &
            &  , c_CUSd_L(2,6), c_CUSd_R(2,6)
  Complex(dp) :: A1L, A1R, A2L, A2R, FL, FR, FLL, FLR, FRL, FRR &
      & , B1L, B1R, B2L, B2R, B3L, B3R, B4L, B4R

  Real(dp), Parameter :: oo576Pi2 = 1._dp / (576._dp * Pi2)
  Complex(dp), Parameter :: mat0(3,3) = 0._dp

  Iname = Iname + 1
  NameOfUnit(Iname) = "MuEconversion"

  BR = 0._dp

  !------------
  ! couplings
  !------------
  c_LNSl_L = 0._dp
  c_LNSl_R = 0._dp
  c_DNSd_L = 0._dp
  c_DNSd_R = 0._dp
  c_UNSu_L = 0._dp
  c_UNSu_R = 0._dp
  Do i1=1,6
   Do i2=1,4
    Call CoupNeutralinoSlepton3(1, i2, i1, gp, g, RSl, id3c, id3c, Y_l, N &
                                  &, c_LNSl_L(1,i2,i1), c_LNSl_R(1,i2,i1) )
    Call CoupNeutralinoSlepton3(2, i2, i1, gp, g, RSl, id3c, id3c, Y_l, N &
                                  &, c_LNSl_L(2,i2,i1), c_LNSl_R(2,i2,i1) )
    Call CoupNeutralinoSdown(1, i2, i1, gp, g, RSd, uD_L, uD_R, Y_d, N &
                                  &, c_DNSd_L(i2,i1), c_DNSd_R(i2,i1) )
    Call CoupNeutralinoSup(1, i2, i1, gp, g, RSu, uU_L, uU_R, Y_u, N &
                                  &, c_UNSu_L(i2,i1), c_UNSu_R(i2,i1) )
   End Do
  End Do

  c_CLSn_L = 0._dp
  c_CLSn_R = 0._dp
  c_CDSu_L = 0._dp
  c_CDSu_R = 0._dp
  c_CUSd_L = 0._dp
  c_CUSd_R = 0._dp
  Do i1=1,3
   Do i2=1,2
    Call CoupCharginoSfermion3(i2, 2, i1, g, -0.5_dp, RSn, Y_l, mat0, id3c &
                  & , id3c, U, V, c_CLSn_L(i2, 2, i1), c_CLSn_R(i2, 2, i1) )
    Call CoupCharginoSfermion3(i2, 1, i1, g, -0.5_dp, RSn, Y_l, mat0, id3c &
                  & , id3c, U, V, c_CLSn_L(i2, 1, i1), c_CLSn_R(i2, 1, i1) )
   End Do
  End Do
  Do i3=1,6
   Do i2=1,2
    Call CoupCharginoSfermion3(i2, 1, i3, g, -0.5_dp, RSd, Y_d, Y_u, uD_L &
                   & , uD_R, U, V, c_CDSu_L(i2, i3), c_CDSu_R(i2, i3) )
    Call CoupCharginoSfermion3(i2, 1, i3, g, 0.5_dp, RSu, Y_d, Y_u, uU_L &
                   & , uU_R, U, V, c_CUSd_L(i2, i3), c_CUSd_R(i2, i3) )
   End Do
  End Do

  cosW2 = mW2 / mZ2 !g**2 / (gp**2 + g**2)
  cosW = Sqrt(cosW2)
  Do i1=1,4
   Do i2=1,4
    Call CoupNeutralinoZ(i1, i2, N, g, cosW, ELn(i1,i2), ERn(i1,i2) )
   End Do
  End Do

  Do i1=1,2
   Do i2=1,2
    Call CoupCharginoZ(i1, i2, U, V, g, cosW, ELc(i1,i2), ERc(i1,i2) )
   End Do
  End Do

  sinW2 = 1._dp - cosW2
  Do i1=1,6
   Do i2=1,6
    Call CoupSleptonZ(i1, i2, g, sinW2, RSl, QSlepton(i1,i2) )
   End Do
  End Do
  Do i1=1,3
   Do i2=1,3
    Call CoupSneutrinoZ3(i1, i2, g, sinW2, RSn, QSnu(i1,i2) )
   End Do
  End Do
  Qslepton = - Qslepton ! sign difference in paper by 
  Qsnu = - Qsnu ! sign difference in paper by 

  Call CoupFermionZ(-0.5_dp, -1._dp, g, sinW2, c_Zee_L,c_Zee_R)

  !---------------------------------------------
  ! Penguin contributions
  !---------------------------------------------
  ! neutralinos
  !--------------
  A1L = 0._dp
  A1R = 0._dp
  A2L = 0._dp
  A2R = 0._dp
  mN2 = mN**2

  Do i1=1,6
   Do i2=1,4
    X_ax = mN2(i2)/mSlep2(i1)
    fun1 = ( 2._dp - 9._dp*X_ax + 18._dp*X_ax**2 - 11._dp*X_ax**3 & 
        & + 6._dp*X_ax**3 * Log(X_ax) ) /  (1._dp-X_ax)**4
    A1L = A1L &
      & + c_LNSl_R(1,i2,i1) * Conjg(c_LNSl_R(2,i2,i1)) * fun1 / mSlep2(i1)
    A1R = A1R &
      & + c_LNSl_L(1,i2,i1) * Conjg(c_LNSl_L(2,i2,i1)) * fun1 / mSlep2(i1)

    fun1 = 2._dp * F2(X_ax)
    fun2 = 2._dp * F4(X_ax)
    A2L = A2L + ( ( c_LNSl_L(1,i2,i1) * Conjg(c_LNSl_L(2,i2,i1))      & 
        &         + c_LNSl_R(1,i2,i1) * Conjg(c_LNSl_R(2,i2,i1))      &
        &           *mf_l(1)/mf_l(2) ) * fun1                         &
        &       + c_LNSl_L(1,i2,i1) * Conjg(c_LNSl_R(2,i2,i1))        &
        &         * mN(i2)/mf_l(2)*fun2 ) / mSlep2(i1)
    A2R = A2R + ( ( c_LNSl_R(1,i2,i1) * Conjg(c_LNSl_R(2,i2,i1))      & 
        &         + c_LNSl_L(1,i2,i1) * Conjg(c_LNSl_L(2,i2,i1))      &
        &           *mf_l(1)/mf_l(2) ) * fun1                         &
        &       + c_LNSl_R(1,i2,i1) * Conjg(c_LNSl_L(2,i2,i1))        &
        &         * mN(i2)/mf_l(2)*fun2 ) / mSlep2(i1)
   End Do
  End Do

  !--------------
  ! charginos
  !--------------
  mC2 = mC**2
  Do i1=1,3
   Do i2=1,2
    X_ax=mC2(i2)/mSnu2(i1)
    fun1 = ( 16._dp - 45._dp*X_ax + 36._dp * X_ax**2 - 7._dp*X_ax**3 &
         & + 6._dp * (2._dp-3._dp*X_ax) * Log(X_ax) ) / (1._dp-X_ax)**4
    A1L = A1L - c_CLSn_R(i2,1,i1) * Conjg(c_CLSn_R(i2,2,i1)) * fun1 / mSnu2(i1)
    A1R = A1R - c_CLSn_L(i2,1,i1) * Conjg(c_CLSn_L(i2,2,i1)) * fun1 / mSnu2(i1)
    fun1 = 2._dp * F1(X_ax)
    fun2 = 2._dp * F3(X_ax)
    A2L = A2L - ( ( c_CLSn_L(i2,1,i1)*Conjg(c_CLSn_L(i2,2,i1))                 & 
        &         + c_CLSn_R(i2,1,i1)*Conjg(c_CLSn_R(i2,2,i1))     &
        &           *mf_l(1)/mf_l(2) ) *fun1                       &
        &       + c_CLSn_L(i2,1,i1) * Conjg(c_CLSn_R(i2,2,i1))     &
        &          * mC(i2) / mf_l(2) * fun2 ) / mSnu2(i1)
    A2R = A2R - ( ( c_CLSn_R(i2,1,i1)*Conjg(c_CLSn_R(i2,2,i1))                 & 
        &         + c_CLSn_L(i2,1,i1)*Conjg(c_CLSn_L(i2,2,i1))     &
        &           *mf_l(1)/mf_l(2) ) *fun1                       &
        &       + c_CLSn_R(i2,1,i1) * Conjg(c_CLSn_L(i2,2,i1))     &
        &          * mC(i2) / mf_l(2) * fun2 ) / mSnu2(i1)
   End Do
  End Do

  A1L = oo576Pi2 * A1L 
  A1R = oo576Pi2 * A1R 
  A2L = oo32Pi2 * A2L 
  A2R = oo32Pi2 * A2R 

  !-----------
  ! Z-penguin
  !-----------
  FL = 0._dp
  FR = 0._DP

  Do i1=1,6
   Do i2=1,4
    Do i3=1,4
     fun1 = 0.5_dp * vertexC0tilde(mSlep2(i1),mN2(i2),mN2(i3))
     fun2 = mN(i2) * mN(i3) * C0_3m(mSlep2(i1),mN2(i2),mN2(i3))
     FL = FL - c_LNSl_R(1,i3,i1)*Conjg(c_LNSl_R(2,i2,i1))             & 
      &   * ( ERn(i3,i2) * fun1 - ELn(i3,i2) * fun2 )
     FR = FR - c_LNSl_L(1,i3,i1)*Conjg(c_LNSl_L(2,i2,i1))             & 
      &   * ( ELn(i3,i2) * fun1 - ERn(i3,i2) * fun2 )
    End Do
   End Do
  End Do

   Do i1=1,4
    Do i2=1,6
     Do i3=1,6
      fun1 = 0.5_dp * Qslepton(i2,i3)  &
           &        * vertexC0tilde(mN2(i1),mSlep2(i2),mSlep2(i3)) 
      FL = FL - c_LNSl_R(1,i1,i2) * Conjg(c_LNSl_R(2,i1,i3)) * fun1
      FR = FR - c_LNSl_L(1,i1,i2) * Conjg(c_LNSl_L(2,i1,i3)) * fun1
     End Do
    End Do
   End Do

   Do i1=1,6
    Do i2=1,4
     fun1 = B1(0._dp, mN2(i2), mSlep2(i1) )
     FL = FL - c_LNSl_R(1,i2,i1) * Conjg(c_LNSl_R(2,i2,i1)) * c_Zee_L * fun1
     FR = FR - c_LNSl_L(1,i2,i1) * Conjg(c_LNSl_L(2,i2,i1)) * c_Zee_R * fun1
    End Do
   End Do

   Do i1=1,3
    Do i2=1,2
     Do i3=1,2
      fun1 = 0.5_dp * vertexC0tilde(mSnu2(i1), mC2(i2), mC2(i3) )
      fun2 = C0_3m(mSnu2(i1), mC2(i2), mC2(i3) )
      FL = FL - c_CLSn_R(i3,1,i1)*Conjg(c_CLSn_R(i2,2,i1) )             & 
         &  * ( ERc(i3,i2) * fun1 - ELc(i3,i2) * mC(i2) * mC(i3) * fun2 )
      FR = FR - c_CLSn_L(i3,1,i1)*Conjg(c_CLSn_L(i2,2,i1) )             & 
         &  * ( ELc(i3,i2) * fun1 - ERc(i3,i2) * mC(i2) * mC(i3) * fun2 )
     End Do
    End Do
   End Do

   Do i1=1,3
    Do i2=1,3
     Do i3=1,2
      fun1 = 0.5_dp * vertexC0tilde(mC2(i3), mSnu2(i1), mSnu2(i2) )
      FL = FL - c_CLSn_R(i3,1,i1) * Conjg(c_CLSn_R(i3,2,i2)) &
         &       * Qsnu(i1,i2) * fun1
      FR = FR - c_CLSn_L(i3,1,i1) * Conjg(c_CLSn_L(i3,2,i2)) &
         &       * Qsnu(i1,i2) * fun1
     End Do
    End Do
   End Do

   Do i1=1,3
    Do i2=1,2
     fun1 = B1(0._dp, mC2(i2), mSnu2(i1) )
     FL = FL - c_CLSn_R(i2,1,i1)*Conjg(c_CLSn_R(i2,2,i1)) * c_Zee_L * fun1
     FR = FR - c_CLSn_L(i2,1,i1)*Conjg(c_CLSn_L(i2,2,i1)) * c_Zee_R * fun1
    End Do
   End Do

  FL = oo16pi2 * FL
  FR = oo16pi2 * FR

  fun1 = 1._dp / (g**2 * sinW2 * mZ2)
  FLL = - fun1 * c_Zee_L * FL 
  FLR = - fun1 * c_Zee_R * FL 
  FRL = - fun1 * c_Zee_L * FR 
  FRR = - fun1 * c_Zee_R * FR

  !-------------------
  ! box contributions
  !-------------------
  B1L = 0._dp   ! u-neutralino
  B1R = 0._dp
  B2L = 0._dp   ! d-neutralino
  B2R = 0._dp
  B3L = 0._dp   ! u-chargino
  B3R = 0._dp
  B4L = 0._dp   ! d-chargino
  B4R = 0._dp


  Do i1=1,6
   Do i2=1,6
    Do i3=1,4
     Do i4=1,4
      fun1 = -0.5_dp * D27_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSup2(i2) ) &
         & + 0.25_dp * mN(i3) * mN(i4)                                      &
         &   * D0_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSup2(i2) )

      B1L = B1L + fun1 * Conjg(c_LNSl_R(2,i3,i1)) * c_LNSl_R(1,i4,i1)  &
          &            * ( c_UNSu_R(i3,i2) * Conjg(c_UNSu_R(i4,i2))    &
          &              - c_UNSu_L(i3,i2) * Conjg(c_UNSu_L(i4,i2)) )

      B1R = B1R + fun1 * Conjg(c_LNSl_L(2,i3,i1)) * c_LNSl_L(1,i4,i1)  &
          &            * ( c_UNSu_L(i3,i2) * Conjg(c_UNSu_L(i4,i2))    &
          &              - c_UNSu_R(i3,i2) * Conjg(c_UNSu_R(i4,i2)) ) 

      fun1 = -0.5_dp * D27_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSdown2(i2) ) &
         & + 0.25_dp * mN(i3) * mN(i4)                                        &
         &   * D0_Bagger(mN2(i3), mN2(i4), mSlep2(i1), mSdown2(i2) )

      B2L = B1L + fun1 * Conjg(c_LNSl_R(2,i3,i1)) * c_LNSl_R(1,i4,i1)  &
          &            * ( c_DNSd_R(i3,i2) * Conjg(c_DNSd_R(i4,i2))    &
          &              - c_DNSd_L(i3,i2) * Conjg(c_DNSd_L(i4,i2)) )

      B2R = B1R + fun1 * Conjg(c_LNSl_L(2,i3,i1)) * c_LNSl_L(1,i4,i1)  &
          &            * ( c_DNSd_L(i3,i2) * Conjg(c_DNSd_L(i4,i2))    &
          &              - c_DNSd_R(i3,i2) * Conjg(c_DNSd_R(i4,i2)) )

     End Do
    End Do
   End Do
  End Do

  Do i1=1,3
   Do i2=1,3
    Do i3=1,2
     Do i4=1,2
      fun1 = 0.5_dp * D27_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSup2(i2) )
      fun2 = -0.25_dp * mC(i3) * mC(i4) &
           &         * D0_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSup2(i2))

      B3L = B3L + fun1 * c_CUSd_R(i3,i2) * Conjg(c_CLSn_R(i3,2,i1)) &
          &            * c_CLSn_R(i4,1,i1) * Conjg(c_CUSd_R(i4,i2)) &
          &     + fun2 * c_CUSd_L(i3,i2) * Conjg(c_CLSn_R(i3,2,i1)) &
          &            * c_CLSn_R(i4,1,i1) * Conjg(c_CUSd_L(i4,i2))

      B3R = B3R + fun1 * c_CUSd_L(i3,i2) * Conjg(c_CLSn_L(i3,2,i1)) &
          &            * c_CLSn_L(i4,1,i1) * Conjg(c_CUSd_L(i4,i2)) &
          &     + fun2 * c_CUSd_R(i3,i2) * Conjg(c_CLSn_L(i3,2,i1)) &
          &            * c_CLSn_L(i4,1,i1) * Conjg(c_CUSd_R(i4,i2))

      fun1 = -0.5_dp * D27_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSdown2(i2) )
      fun2 = 0.25_dp * mC(i3) * mC(i4) &
           &         * D0_Bagger(mC2(i3), mC2(i4), mSnu2(i1), mSdown2(i2))

      B4L = B4L + fun1 * c_CDSu_R(i3,i2) * Conjg(c_CLSn_R(i3,2,i1)) &
          &            * c_CLSn_R(i4,1,i1) * Conjg(c_CDSu_R(i4,i2)) &
          &     + fun2 * c_CDSu_L(i3,i2) * Conjg(c_CLSn_R(i3,2,i1)) &
          &            * c_CLSn_R(i4,1,i1) * Conjg(c_CDSu_L(i4,i2))

      B4R = B4R + fun1 * c_CDSu_L(i3,i2) * Conjg(c_CLSn_L(i3,2,i1)) &
          &            * c_CLSn_L(i4,1,i1) * Conjg(c_CDSu_L(i4,i2)) &
          &     + fun2 * c_CDSu_R(i3,i2) * Conjg(c_CLSn_L(i3,2,i1)) &
          &            * c_CLSn_L(i4,1,i1) * Conjg(c_CDSu_R(i4,i2))

     End Do
    End Do
   End Do
  End Do

  fun1 = oo16pi2 / (g**2*sinW2)
  B1L = fun1 * B1L
  B1R = fun1 * B1R
  B2L = fun1 * B2L
  B2R = fun1 * B2R
  B3L = fun1 * B3L
  B3R = fun1 * B3R
  B4L = fun1 * B4L
  B4R = fun1 * B4R

  !-------------------------
  ! now the branching ratio
  !-------------------------
  BR = Abs(A1L)**2 + Abs(A1R)**2                                         &
   & - 4._dp* Real(A1L*Conjg(A2R)+A2L*Conjg(A1R),dp)                     &
   & + (Abs(A2L)**2 + Abs(A2R)**2)                                       &
   &     * (16._dp * Log(mf_l(2)/mf_l(1)) - 22._dp ) /3._dp              &
   & + ( Abs(B1L)**2 + Abs(B1R)**2 ) / 6._dp                             &
   & + ( Abs(B2L)**2 + Abs(B2R)**2 ) / 3._dp                             &
   & + ( Abs(B3L)**2 + Abs(B3R)**2 ) / 24._dp                            &
   & + 6._dp * ( Abs(B4L)**2 + Abs(B4R)**2 )                             &
   & - Real( B3L * Conjg(B4L) + B3R * Conjg(B4R),dp)                     &
   & + 2._dp * Real( A1L * Conjg(B1L) + A1R * Conjg(B1R)                 &
   &               + A1L * Conjg(B2L) + A1R * Conjg(B2R), dp) / 3._dp    &
   & - 4._dp * Real( A2R * Conjg(B1L) + A2L * Conjg(B1R)                 &
   &               + A2L * Conjg(B2R) + A2R * Conjg(B2L), dp) / 3._dp    &
   & + ( 2._dp * (Abs(FLL)**2 + Abs(FRR)**2) + Abs(FLR)**2 + Abs(FRL)**2 &
   &   + 2._dp * Real( B1L * Conjg(FLL) + B1R * Conjg(FRR)               &
   &                 + B2L * Conjg(FLR) + B2R * Conjg(FRL), dp )         &
   &   + 4._dp * Real( A1L * Conjg(FLL) + A1R * Conjg(FRR), dp )         &
   &   + 2._dp * Real( A1L * Conjg(FLR) + A1R * Conjg(FRL), dp )         &
   &   - 8._dp * Real( A2R * Conjg(FLL) + A2L * Conjg(FRR), dp )         &
   &   - 4._dp * Real( A2L * Conjg(FRL) + A2R * Conjg(FLR), dp ) ) / 3._dp

  !----------------------------------------------------------------------
  ! taking alpha(Q=0) instead of alpha(m_Z) as this contains most of the
  ! running of the Wilson coefficients
  !----------------------------------------------------------------------
  ! Formula according to Hisano et al, needs to be improved by the one
  ! of Kitano
  BR = 4._dp * Alpha**5 * 17.6_dp**4 * 0.54_dp**2 / Znuc * mf_l(2)**5 &
     &       * ( ( Znuc * (A1L - A2R) - 0._dp )**2                    &
               +  ( Znuc * (A1R - A2L) - 0._dp )**2  )
                  
  Iname = Iname - 1

 Contains

  Complex(dp) Function vertexC0tilde(masa1c,masa2c,masa3c)
   Real(dp), Intent(in) :: masa1c,masa2c,masa3c

   vertexC0tilde = 0._dp

   If (masa2c.Eq.masa3c) Then

   vertexC0tilde=(masa1c*(masa1c-masa2c)-masa1c**2*Log(masa1c)           & 
        & +(2._dp*masa1c-masa2c)*masa2c*Log(masa2c))/(masa1c-masa2c)**2   

   Else

   vertexC0tilde=1._dp-1._dp/(masa2c-masa3c)*((masa1c**2*Log(masa1c)      & 
        & -masa2c**2*Log(masa2c))/(masa1c-masa2c)-(masa1c**2*Log(masa1c)  & 
        & -masa3c**2*Log(masa3c))/(masa1c-masa3c))

   End If

  End Function vertexC0tilde

 End Subroutine MuEconversion
! Interface LtoLpGamma


 Subroutine ZtoLiLj(i, j, UseSavedLoopFunctions  &
       & , mC, mC2, cpl_CLSn_L, cpl_CLSn_R, mSnu2, cpl_SnSnZ &
       & , mN, mN2, cpl_LNSl_L, cpl_LNSl_R, mSle2, cpl_SlSlZ &
       & , cpl_CCZ_L, cpl_CCZ_R, cpl_NNZ_L, cpl_NNZ_R        &
       & , BR) 
 !-------------------------------------------------------------------------
 ! in this subroutine the decay branching ratio for the decay
 ! Z -> l_i  l_j is calculated.
 ! The formulas of X.Bi et al, PRD63, 096008, (2001) are used.
 ! Input: i,j ....................... indices of the leptons
 !         UseSavedLoopFunctions .... the Loopfunctions calculated in a
 !                                    previous run are used if set .true.
 ! written by Werner Porod, 17.10.01 
 !-------------------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: i, j
  Logical, Intent(in) :: UseSavedLoopFunctions
  Real(dp), Intent(in) :: mC(2), mC2(2), mSnu2(3), mN(4), mN2(4), mSle2(6)
  Complex(dp), Intent(in) ::  cpl_SnSnZ(:,:), cpl_CLSn_L(:,:,:)            &
     & , cpl_CLSn_R(:,:,:), cpl_LNSl_L(:,:,:), cpl_LNSl_R(:,:,:)           &
     & , cpl_SlSlZ(:,:), cpl_CCZ_L(:,:), cpl_CCZ_R(:,:), cpl_NNZ_L(:,:)    &
     & , cpl_NNZ_R(:,:)
  Real(dp), Intent(out) :: BR

  !--------------------------
  ! to save calculation time
  !--------------------------
  Complex(dp), Save :: C00_C_i_Snu_jk(2,3), C00_N_i_Sle_jk(4,6,6)   &
     & , C00_Snu_i_C_jk(3,2,2), C00_Sle_i_N_jk(6,4,4)               &
     & , C0_C_i_Snu_jk(2,3), C0_N_i_Sle_jk(4,6,6)                 &
     & , C1_C_i_Snu_jk(2,3), C1_N_i_Sle_jk(4,6,6)                 &
     & , C2_C_i_Snu_jk(2,3), C2_N_i_Sle_jk(4,6,6)                 &
     & , C0_Snu_i_C_jk(3,2,2), C0_Sle_i_N_jk(6,4,4)                 &
     & , C1_Snu_i_C_jk(3,2,2), C1_Sle_i_N_jk(6,4,4)                 &
     & , C2_Snu_i_C_jk(3,2,2), C2_Sle_i_N_jk(6,4,4)                 &
     & , B0_Z_CC(2,2), B0_0_CSn(2,3), B0_Z_NN(4,4), B0_0_NSl(4,6)   &
     & , B1_Z_CC(2,2), B1_0_CSn(2,3), B1_Z_NN(4,4), B1_0_Nsl(4,6)!     &

  Integer :: i1, i2, i3
  Real(dp) :: mZ2a
  Complex(dp) :: Ai(4)

  Iname = Iname + 1
  NameOfUnit(Iname) = "ZtoLiLj"

  !------------------------------------------------------
  ! Calculate first the loop functions if necessary
  !------------------------------------------------------
  If (.Not.UseSavedLoopFunctions) Then
   mZ2a = 0._dp ! the leptons are on-shell
   Do i1=1,2
    B0_Z_CC(i1,i1) = B0(mZ2, mC2(i1), mC2(i1) )
    B1_Z_CC(i1,i1) = B1(mZ2, mC2(i1), mC2(i1) )
    Do i2=i1+1,2
     B0_Z_CC(i1,i2) = B0(mZ2, mC2(i1), mC2(i2) )
     B1_Z_CC(i1,i2) = B1(mZ2, mC2(i1), mC2(i2) )
     B0_Z_CC(i2,i1) = B0(mZ2, mC2(i2), mC2(i1) )
     B1_Z_CC(i2,i1) = B1(mZ2, mC2(i2), mC2(i1) )
    End Do
    Do i2=1,3
     B0_0_CSn(i1,i2) = B0(mZ2a, mC2(i1), mSnu2(i2) )
     B1_0_CSn(i1,i2) = B1(mZ2a, mC2(i1), mSnu2(i2) )
     C00_C_i_Snu_jk(i1,i2) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C0_C_i_Snu_jk(i1,i2) = &
          & C0( 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C1_C_i_Snu_jk(i1,i2) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C2_C_i_Snu_jk(i1,i2) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mC2(i1), mSnu2(i2), mSnu2(i2) )
     C0_Snu_i_C_jk(i2,i1,i1)  = &
          & C0( 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     C1_Snu_i_C_jk(i2,i1,i1)  = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     C2_Snu_i_C_jk(i2,i1,i1)  = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     C00_Snu_i_C_jk(i2,i1,i1)  = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1) )
     If (i1.Eq.1) Then 
      C0_Snu_i_C_jk(i2,i1,i1+1)  = &
            & C0( 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C0_Snu_i_C_jk(i2,i1+1,i1) = C0_Snu_i_C_jk(i2,i1,i1+1)
      C1_Snu_i_C_jk(i2,i1,i1+1)  = &
            & Cget("C1  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C1_Snu_i_C_jk(i2,i1+1,i1) = C1_Snu_i_C_jk(i2,i1,i1+1)
      C2_Snu_i_C_jk(i2,i1,i1+1)  = &
            & Cget("C2  ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C2_Snu_i_C_jk(i2,i1+1,i1) = C2_Snu_i_C_jk(i2,i1,i1+1)
      C00_Snu_i_C_jk(i2,i1,i1+1)  = &
            & Cget("C00 ", 0._dp, mZ2, 0._dp, mSnu2(i2), mC2(i1), mC2(i1+1) )
      C00_Snu_i_C_jk(i2,i1+1,i1) = C00_Snu_i_C_jk(i2,i1,i1+1)
     End If
    End Do
   End Do
   Do i1=1,4
    B0_Z_NN(i1,i1) = B0(mZ2, mN2(i1), mN2(i1) )
    B1_Z_NN(i1,i1) = B1(mZ2, mN2(i1), mN2(i1) )
    Do i2=i1+1,4
     B0_Z_NN(i1,i2) = B0(mZ2, mN2(i1), mN2(i2) )
     B1_Z_NN(i1,i2) = B1(mZ2, mN2(i1), mN2(i2) )
     B0_Z_NN(i2,i1) = B0(mZ2, mN2(i2), mN2(i1) )
     B1_Z_NN(i2,i1) = B1(mZ2, mN2(i2), mN2(i1) )
    End Do
    Do i2=1,6
     B0_0_NSl(i1,i2) = B0(mZ2a, mN2(i1), mSle2(i2) )
     B1_0_NSl(i1,i2) = B1(mZ2a, mN2(i1), mSle2(i2) )
     C00_N_i_Sle_jk(i1,i2,i2) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     C0_N_i_Sle_jk(i1,i2,i2) = &
          & C0( 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     C1_N_i_Sle_jk(i1,i2,i2) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     C2_N_i_Sle_jk(i1,i2,i2) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i2) )
     Do i3=i2+1,6
      C00_N_i_Sle_jk(i1,i2,i3) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C00_N_i_Sle_jk(i1,i3,i2) = C00_N_i_Sle_jk(i1,i2,i3)
      C0_N_i_Sle_jk(i1,i2,i3) = &
          & C0( 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C0_N_i_Sle_jk(i1,i3,i2) = C0_N_i_Sle_jk(i1,i2,i3)
      C1_N_i_Sle_jk(i1,i2,i3) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C1_N_i_Sle_jk(i1,i3,i2) = C1_N_i_Sle_jk(i1,i2,i3)
      C2_N_i_Sle_jk(i1,i2,i3) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mN2(i1), mSle2(i2), mSle2(i3) )
      C2_N_i_Sle_jk(i1,i3,i2) = C2_N_i_Sle_jk(i1,i2,i3)
     End Do
     C0_Sle_i_N_jk(i2,i1,i1) = &
          & C0( 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     C1_Sle_i_N_jk(i2,i1,i1) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     C2_Sle_i_N_jk(i2,i1,i1) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     C00_Sle_i_N_jk(i2,i1,i1) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i1) )
     Do i3=i1+1,4
      C0_Sle_i_N_jk(i2,i1,i3) = &
          & C0( 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C0_Sle_i_N_jk(i2,i3,i1) = C0_Sle_i_N_jk(i2,i1,i3)
      C1_Sle_i_N_jk(i2,i1,i3) = &
          & Cget("C1  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C1_Sle_i_N_jk(i2,i3,i1) = C1_Sle_i_N_jk(i2,i1,i3)
      C2_Sle_i_N_jk(i2,i1,i3) = &
          & Cget("C2  ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C2_Sle_i_N_jk(i2,i3,i1) = C2_Sle_i_N_jk(i2,i1,i3)
      C00_Sle_i_N_jk(i2,i1,i3) = &
          & Cget("C00 ", 0._dp, mZ2, 0._dp, mSle2(i2), mN2(i1), mN2(i3) )
      C00_Sle_i_N_jk(i2,i3,i1) = C00_Sle_i_N_jk(i2,i1,i3)
     End Do
    End Do
   End Do
  End If

  !------------------------------------------------------
  ! now the effective couplings
  !------------------------------------------------------
  Ai = 0._dp
  Do i1=1,2
   Do i2=1,3
    Ai(1) = Ai(1)                                                           &
        & + 2._dp * cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_R(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * C00_C_i_Snu_jk(i1,i2)
    Ai(2) = Ai(2)                                                           &
        & + 2._dp * cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_L(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * C00_C_i_Snu_jk(i1,i2)
    Ai(3) = Ai(3)                                                           &
        & + 2._dp * cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_L(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * mC(i1) * ( C0_C_i_Snu_jk(i1,i2)      &
        &                                       + C1_C_i_Snu_jk(i1,i2)      &
        &                                       + C2_C_i_Snu_jk(i1,i2) )
    Ai(4) = Ai(4)                                                           &
        & + 2._dp * cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_R(i1, j, i2) )  &
        &         * cpl_SnSnZ(i2,i2) * mC(i1) * ( C0_C_i_Snu_jk(i1,i2)      &
        &                                       + C1_C_i_Snu_jk(i1,i2)      &
        &                                       + C2_C_i_Snu_jk(i1,i2) )
    Do i3=1,2
     Ai(1) = Ai(1)                                                           &
      &  - cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_R(i3, j, i2) )            &
      &    * ( cpl_CCZ_R(i1,i3) * ( mSnu2(i2) * C0_Snu_i_C_jk(i2,i1,i3)      &
      &                           + B0_Z_CC(i1,i3)                           &
      &                           - C00_Snu_i_C_jk(i2,i1,i3) )               & 
      &      - mC(i1) *mC(i3) *cpl_CCZ_L(i1,i3) * C0_Snu_i_C_jk(i2,i1,i3) )
     Ai(2) = Ai(2)                                                           &
      &  - cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_L(i3, j, i2) )            &
      &    * ( cpl_CCZ_L(i1,i3) * ( mSnu2(i2) * C0_Snu_i_C_jk(i2,i1,i3)      &
      &                           + B0_Z_CC(i1,i3)                           &
      &                           - C00_Snu_i_C_jk(i2,i1,i3) )               & 
      &      - mC(i1) *mC(i3) *cpl_CCZ_R(i1,i3) * C0_Snu_i_C_jk(i2,i1,i3) )
     Ai(3) = Ai(3)                                                           &
      &  + cpl_CLSn_R(i1, i, i2) * Conjg( cpl_CLSn_L(i3, j, i2) )            &
      &    * ( cpl_CCZ_L(i1,i3) * mC(i3) *  C1_Snu_i_C_jk(i2,i1,i3)          &
      &      + cpl_CCZ_R(i1,i3) * mC(i1) *  C1_Snu_i_C_jk(i2,i1,i3)  )
     Ai(4) = Ai(4)                                                           &
      &  + cpl_CLSn_L(i1, i, i2) * Conjg( cpl_CLSn_R(i3, j, i2) )            &
      &    * ( cpl_CCZ_R(i1,i3) * mC(i3) *  C1_Snu_i_C_jk(i2,i1,i3)          &
      &      + cpl_CCZ_L(i1,i3) * mC(i1) *  C1_Snu_i_C_jk(i2,i1,i3)  )
    End Do
   End Do
  End Do
  Do i1=1,4
   Do i2=1,6
    Ai(1) = Ai(1)                                                           &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * C00_N_i_Sle_jk(i1,i2,i2)
    Ai(2) = Ai(2)                                                           &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * C00_N_i_Sle_jk(i1,i2,i2)
    Ai(3) = Ai(3)                                                           &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i2) )
    Ai(4) = Ai(4)                                                           &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i2,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i2) )
    Do i3=i2+1,6
     Ai(1) = Ai(1)                                                          &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * C00_N_i_Sle_jk(i1,i2,i3)             &
        & + 2._dp * cpl_LNSl_R(i, i1, i3) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * C00_N_i_Sle_jk(i1,i3,i2)
     Ai(2) = Ai(2)                                                          &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * C00_N_i_Sle_jk(i1,i2,i3)             &
        & + 2._dp * cpl_LNSl_L(i, i1, i3) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * C00_N_i_Sle_jk(i1,i3,i2)
     Ai(3) = Ai(3)                                                          &
        & + 2._dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_L(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i3) ) &
        & + 2._dp * cpl_LNSl_R(i, i1, i3) * Conjg( cpl_LNSl_L(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i3,i2) )
     Ai(4) = Ai(4)                                                          &
        & + 2._dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_R(j, i1, i3) )  &
        &         * cpl_SlSlZ(i2,i3) * mN(i1) * ( C0_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C1_N_i_Sle_jk(i1,i2,i3)   &
        &                                       + C2_N_i_Sle_jk(i1,i2,i3) ) &
        & + 2._dp * cpl_LNSl_L(i, i1, i3) * Conjg( cpl_LNSl_R(j, i1, i2) )  &
        &         * cpl_SlSlZ(i3,i2) * mN(i1) * ( C0_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C1_N_i_Sle_jk(i1,i3,i2)   &
        &                                       + C2_N_i_Sle_jk(i1,i3,i2) )
    End Do
    Do i3=1,4
     Ai(1) = Ai(1)                                                           &
      &  - 0.5_dp * cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_R(j, i3, i2) )   &
      &    * ( cpl_NNZ_R(i1,i3) * ( mSle2(i2) * C0_Sle_i_N_jk(i2,i1,i3)      &
      &                           + B0_Z_NN(i1,i3)                           &
      &                           - C00_Sle_i_N_jk(i2,i1,i3) )               & 
      &      - mN(i1) *mN(i3) *cpl_NNZ_L(i1,i3) * C0_Sle_i_N_jk(i2,i1,i3) )
     Ai(2) = Ai(2)                                                           &
      &  - 0.5_dp * cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_L(j, i3, i2) )   &
      &    * ( cpl_NNZ_L(i1,i3) * ( mSle2(i2) * C0_Sle_i_N_jk(i2,i1,i3)      &
      &                           + B0_Z_NN(i1,i3)                           &
      &                           - C00_Sle_i_N_jk(i2,i1,i3) )               & 
      &      - mN(i1) *mN(i3) *cpl_NNZ_R(i1,i3) * C0_Sle_i_N_jk(i2,i1,i3) )
     Ai(3) = Ai(3)                                                           &
      &  + cpl_LNSl_R(i, i1, i2) * Conjg( cpl_LNSl_L(j, i3, i2) )            &
      &    * ( cpl_NNZ_L(i1,i3) * mN(i3) *  C1_Sle_i_N_jk(i2,i1,i3)          &
      &      + cpl_NNZ_R(i1,i3) * mN(i1) *  C1_Sle_i_N_jk(i2,i1,i3)  )
     Ai(4) = Ai(4)                                                           &
      &  + cpl_LNSl_L(i, i1, i2) * Conjg( cpl_LNSl_R(j, i3, i2) )            &
      &    * ( cpl_NNZ_R(i1,i3) * mN(i3) *  C1_Sle_i_N_jk(i2,i1,i3)          &
      &      + cpl_NNZ_L(i1,i3) * mN(i1) *  C1_Sle_i_N_jk(i2,i1,i3)  )
    End Do
   End Do
  End Do

  Ai = oo16pi2 * Ai
  BR = oo48Pi * ( 2._dp * (Abs(Ai(1))**2 + Abs(Ai(2))**2 ) * mZ         &
     &          + 0.25_dp * (Abs(Ai(3))**2 + Abs(Ai(4))**2 ) * mZ * mZ2 ) / gamZ

  Iname = Iname - 1
 
 End Subroutine ZtoLiLj

  Subroutine Bs_to_MuMu(i, mf_u, g, Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R        &
      & , mC, U, V, mN, N, mGlu, phi_g, mS02, RS0, mP02, RP0, mSpm2, RSpm   &
      & , mSup2, RSup, A_u, mu, mSdown2, RSdown, A_d, vevSM, mSnu2, mSlept2 &
      & , c_CDSu_L, c_CDSu_R, c_CLSn_L, c_CLSn_R, c_LNSl_L, c_LNSl_R        &
      & , c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R &
      & , res )
! Subroutine Delta_MB(i, mf_u, res )
 !---------------------------------------------------------------------------
 ! Input: mf_u, mC, mN, mGlu, mS0, mP0, mSpm
 !        U, V, N, C
 !        
 ! Output: 
 !         res
 !         
 ! written by Werner Porod, 01 Jul 2003
 ! the eqn numbers ref. to Buras et al., NPB659, 3 (2003), hep-ph/0210145
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_u(3), mC(:), mN(:), mGlu, mS02(:), mP02(:)     &
     & , mSpm2(:), g(3), mSdown2(6), mSup2(6), RS0(:,:), RP0(:,:), vevSM(2) &
     & , mSnu2(3), mSlept2(6)
  Complex(dp), Intent(in) :: U(:,:), V(:,:), N(:,:), RSpm(:,:), RSup(6,6)   &
     & , RSdown(6,6), phi_g, A_u(3,3), A_d(3,3), mu, c_CDSu_L(2,3,6)        &
     & , c_CDSu_R(2,3,6), c_CLSn_L(2,3,3), c_CLSn_R(2,3,3), c_LNSl_L(3,4,6) &
     & , c_LNSl_R(3,4,6), c_DGSd_L(3,6), c_DGSd_R(3,6), c_DNSd_L(3,4,6)     &
     & , c_DNSd_R(3,4,6) 
  Complex(dp), Intent(in), Dimension(3,3) :: Y_u, Ru_L, Ru_R, Y_d, Rd_L, Rd_R
  Integer, Intent(in) :: i

  Real(dp), Intent(out) :: res

  Real(dp) :: xt, a_s, e2, sW2, tanb
  complex(dp) :: kappa_q
  Complex(dp) :: coupLC, coupRC, cs, cp, csp, cpp, fac, ca, cap
  Real(dp), Parameter ::  T3= -0.5_dp, e_d =-1._dp/3._dp
  Integer :: i1

  !-------------------------
  ! SM
  !-------------------------
  xt = mf_u(3)**2 / mW**2

  a_s = oo4pi * g(3)**2
  e2 = 4*pi * alpha
  sW2 = g(1)**2 / (g(1)**2 + g(2)**2)
  kappa_q = 8._dp * Sqrt2 * G_F * e2 * CKM(3,3) * Conjg( CKM(3,i) )
  kappa_q = 1._dp / kappa_q
  tanb = vevSM(2) / vevSM(1)
  Call  C_QdQdLL_AAp(3, i, 2, mW, mf_d_mZ, mf_u, mf_l, mW2, mSpm2(2), mGlu**2     &
       & , mC**2, mN**2, mSnu2, mSlept2, mSup2, mSdown2, tanb, U, V, N, RSup   &
       & , RSdown, c_CDSu_L, c_CDSu_R, c_CLSn_L, c_CLSn_R, c_LNSl_L, c_LNSl_R  &
       & , c_DGSd_L, c_DGSd_R, c_DNSd_L, c_DNSd_R  &
       & , a_s, kappa_q, e2, sW2, .false., cA, cAp)
! I am calculating C_10 and not C_A anymore
  cA = sW2 * cA
  cAp = sW2 * cAp

  !------------------------------------------
  ! double penguins
  !------------------------------------------
  cs = 0._dp
  csp = 0._dp
  Do i1=1,2
    Call  CoupFermionScalar31L_eff(3, i, i1, T3, e_d, g, Y_d, Rd_L, Rd_R, RS0 &
     & , vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup, Y_u &
     & , A_u, mC, U, V, coupLC, coupRC)
   cs = cs - coupLC * RS0(i1,1) / RP0(1,1) / mS02(i1)
   csp = csp - coupRC * RS0(i1,1) / RP0(1,1) / mS02(i1)
  End Do

  Call CoupFermionPseudoScalar31L_eff(3, i, 2, T3, e_d, g, Y_d, Rd_L, Rd_R, RP0&
    & , vevSM, mSdown2, RSdown, A_d, phi_g, mglu, mu, mN, N, mSup2, RSup, Y_u &
    & , A_u, mC, U, V, coupLC, coupRC)
  cp = - coupLC * vevSM(2) / (vevSM(1) * mP02(2))
  cpp = 0._dp
  cpp = coupRC * vevSM(2) / (vevSM(1) * mP02(2))


  !-------------------------------------------------
  ! Eq +6.29+6.33+6.38
  !-------------------------------------------------
!  fac = - 4._dp * Pi2 * MBq(i) * mf_l_mZ(2) / ( mf_d_mZ(3) * mW2 )
  fac = - 4._dp * Pi2 * MBq(i) * mf_l_mZ(2) / ( (mf_d_mZ(3) + mf_d_mZ(i)) * mW2 )
  fac = fac / (G_F**1.5_dp * 2._dp**1.75_dp * CKM(3,i) * Conjg(CKM(3,3)) )

  cs = fac * (cs - csp * mf_d(i)/mf_d(3))
  cp = fac * (cp - cpp * mf_d(i)/mf_d(3))
 !-------------------------------
 ! Eq.6.37
 !-------------------------------
  res = 2.32e-6_dp * (TauB(i)/1.5_dp) * (FBhatB(i)/0.23_dp)**2   &
      &            * (Abs(CKM(3,i)) / 0.04_dp)**2                &
      &            * ( Abs(cs)**2 + Abs(cp + 0.04_dp * (ca-cap))**2  )

 End Subroutine Bs_to_MuMu


 Real(dp) Function F1_0(x)
 !-----------------------------------------------------------
 ! loop function as given by Bobeth et al., NPB 630 (2002) 87
 !-----------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r
  Real(dp), Parameter :: c1(7) = (/ 5._dp/6._dp, -1._dp/24._dp, 0._dp            &
          &                       , 1._dp/120._dp, -1._dp/105._dp, 1._dp/112._dp &
          &                       , -1._dp/126._dp /)

  If (x.Eq.0._dp) Then
   F1_0 = 0._dp

  Else If (x.Eq.1._dp) Then
   F1_0 = 0.75_dp

  Else If (Abs(x-1._dp).Lt.1.e-4_dp) Then
   r = x - 1._dp
   F1_0 = c1(7) * r
   Do i1=1,3
    F1_0 = (c1(7-i1) + F1_0)*r
   End Do
   Do i1=5,6
    F1_0 = (c1(7-i1) + F1_0)*r
   End Do
   F1_0 = 0.75_dp + F1_0

  Else
   F1_0 = - 0.5_dp * x * (6 - x)/(x - 1) &
      &  + 0.5_dp * x * (2 + 3 * x) * Log(x) / (x - 1)**2
  End If

 End Function F1_0


 Real(dp) Function F2_0(x)
 !-----------------------------------------------------------
 ! loop function as given by Bobeth et al., NPB 630 (2002) 87
 !-----------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r
  Real(dp), Parameter :: c1(7) = (/ -1._dp/6._dp, 1._dp/12._dp, -0.05_dp         &
          &                       , 1._dp/30._dp, -1._dp/42._dp, 1._dp/56._dp &
          &                       , -1._dp/72._dp /)

  If (x.Eq.0._dp) Then
   F2_0 = 0._dp

  Else If (x.Eq.1._dp) Then
   F2_0 = - 0.5_dp

  Else If (Abs(x-1._dp).Lt.1.e-4_dp) Then
   r = x - 1._dp
   F2_0 = c1(7) * r
   Do i1=1,6
    F2_0 = (c1(7-i1) + F2_0)*r
   End Do
   F2_0 = -0.5_dp + F2_0

  Else
   F2_0 = - x / (x - 1) + x * Log(x) / (x - 1)**2
  End If

 End Function F2_0


 Real(dp) Function F1(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp), Parameter :: c1(7) = (/ 1._dp / 24._dp, -0.025_dp, 1._dp/60._dp  &
          &                       , -1._dp / 84._dp , 1._dp/112._dp         &
          &                       , -1._dp / 144._dp , 1._dp / 180._dp /)
  Real(dp) :: r
  
  Iname = Iname + 1
  NameOfUnit(Iname) = "F1"

  If (x.Lt.0) Then
   F1 = 0
   If (ErrorLevel.Gt.-2) Then
    Write(ErrCan,*) "Error in Function F1, x=",x
    If (ErrorLevel.Ge.0) Call TerminateProgram
   End If
  End If

  If (x.Eq.0._dp) Then
   F1 = 1._dp / 6._dp
  Else If (x.Eq.1._dp) Then
   F1 = 1._dp / 24._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   F1 = c1(1)
   r = x - 1._dp
   Do i1=1,6
    F1 = F1 + c1(i1+1) * r**i1
   End Do 
  Else
   F1 = ( 2._dp + 3._dp * x - 6._dp * x**2 + x**3 &
      & + 6._dp * x * Log(x) ) / (12._dp * (1._dp - x)**4 )
  End If

  Iname = Iname - 1
 End Function F1

 Real(dp) Function F2(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp), Parameter :: c1(7) = (/ 1._dp/24._dp, -1._dp/ 60._dp         &
          &                       , 1._dp/120._dp, -1._dp / 210._dp      &
          &                       , 1._dp/336._dp, -1._dp / 504._dp       &
          &                       , 1._dp/720._dp  /)
  Real(dp) :: r

  If (x.Eq.0._dp) Then
   F2 = 1._dp / 12._dp
  Else If (x.Eq.1._dp) Then
   F2 = 1._dp / 24._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   F2 = c1(1)
   r = x - 1._dp
   Do i1=1,6
    F2 = F2 + c1(i1+1) * r**i1
   End Do 
  Else
   F2 = ( 1._dp - 6._dp * x + 3._dp * x**2 + 2._dp * x**3 &
      &  - 6._dp * x**2 * Log(x) ) / (12._dp * (1._dp - x)**4 )
  End If

 End Function F2


 Real(dp) Function F3(x)
 Implicit None
  Real(dp), Intent(in) :: x 

  Integer :: i1
  Real(dp) :: r

  If (x.Eq.0._dp) Then
   F3 = Huge(x) ! goes to infinity
  Else If (x.Eq.1._dp) Then
   F3 = 1._dp / 3._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   F3 = 1._dp / 3._dp
   r = x - 1._dp
   Do i1=1,10
    F3 = F3 + r**i1 * (-1)**i1 / (i1 + 3._dp)
   End Do 
  Else
   F3 = (-3._dp + 4._dp * x - x**2 - 2._dp * Log(x) ) &
    & / (2._dp * (1._dp - x)**3 ) 
  End If

 End Function F3

 Real(dp) Function F3gamma(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: x2, r
  Real(dp), Parameter :: c1(10) = (/ -2._dp/3._dp, 1._dp/ 6._dp        &
          &                       , -1._dp/15._dp, 1._dp / 30._dp      &
          &                       , -2._dp/105._dp, 1._dp / 84._dp     &
          &                       , -1._dp/1260._dp, 1._dp/180._dp     &
          &                       , -2._dp/495._dp, 1._dp / 330._dp  /)
      
  If (x.Eq.1._dp) Then
   F3gamma = - 2._dp / 3._dp
  Else If (x.Eq.0._dp) Then
   F3gamma = - 1._dp
  Else If (Abs(1._dp-x).Lt.1.e-2_dp) Then
   F3gamma = c1(1)
   r = x - 1._dp
   Do i1=1,9
    F3gamma = F3gamma + c1(i1+1) * r**i1
   End Do 

  Else
   x2 = x*x
   F3gamma = (1 - 4*x + 3 * x2 - 2 * x2 * Log(x) ) / (x-1)**3
  End If

 End Function F3gamma

 Real(dp) Function F4(x)
 Implicit None
  Real(dp), Intent(in) :: x

  Integer :: i1
  Real(dp) :: r

  If (x.Eq.0._dp) Then
   F4 = 0.5_dp
  Else If (x.Eq.1._dp) Then
   F4 = 1._dp / 6._dp
  Else If (Abs(1._dp - x).Lt.1.e-2_dp) Then
   F4 = 1._dp / 6._dp
   r = x - 1._dp
   Do i1=1,10
    F4 = F4 +  (-r)**i1 / (6._dp + 5._dp * i1 + i1**2)
   End Do 
  Else
   F4 = 0.5_dp * (1+x) / (1-x)**2 + x * Log(x) /(1-x)**3
  End If

 End Function F4


!\section{Function FeynFunctionA}
!\begin{verbatim}
 Real(dp) Function FeynFunctionA(r)
 !----------------------------------------------------------------------
 ! calculation of the Function A(r) ... 0-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.1.01
 !----------------------------------------------------------------------
  Implicit None
   Real(dp), Intent(in) :: r

   Real(dp) :: x, x2
   Integer :: i

   If (r.Eq.1._dp) Then
    FeynFunctionA = - 1._dp / 3._dp
  
   Else

    x = r - 1._dp 
    x2 = x * x
    If (x2 .Lt. 1.e-3_dp) Then
     FeynFunctionA = sa(1)
     Do i=1,6
      FeynFunctionA = FeynFunctionA + sa(i+1)*x**i
     End Do
    Else
     FeynFunctionA = ( 1._dp - 0.5_dp * x - Log(r) / x ) / x2
    End If

   End If

 End Function FeynFunctionA
!\end{verbatim}

!\section{Function FeynFunctionB}
!\begin{verbatim}
 Real(dp) Function FeynFunctionB(r)
 !----------------------------------------------------------------------
 ! calculation of the Function B(r) ... 0-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.01.01
 !----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: r

  Real(dp) :: x, x2
  Integer :: i

  If (r.Eq.1._dp) Then
   FeynFunctionB = 1._dp / 6._dp

  Else
   x = r - 1._dp
   x2 = x * x
   If (x2 .Lt. 1.e-3_dp) Then
    FeynFunctionB  = sb(1)
    Do i=1,6
      FeynFunctionB  = FeynFunctionB  + sb(i+1)*x**i
    End Do
   Else
    FeynFunctionB = ( 1._dp + 0.5_dp * x - r * Log(r) / x ) / x2
   End If

  End If

 End Function FeynFunctionB
!\end{verbatim}

 Real(dp) Function FeynFunctionC(r)
 Implicit None
  Real(dp), Intent(in) :: r

  If (r.Eq.1._dp) Then
   FeynFunctionC = - 19._dp / 18._dp
  Else
   FeynFunctionC = 3._dp * FeynFunctionA(r) - FeynFunctionB(r) / 3._dp
  End If
 End Function FeynFunctionC

!\section{Function FeynFunctionD}
!\begin{verbatim}
 Real(dp) Function FeynFunctionD(r)
 !----------------------------------------------------------------------
 ! calculation of the Function D(r) ... 1-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.01.01
 !---------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: r

  Real(dp) :: x, x2
  Integer :: i

  If (r.Eq.1._dp) Then
   FeynFunctionD = - 0.05_dp
  
  Else
   x = r - 1._dp
   x2 = x*x
   If (x2 .Lt. 1.0d-3) Then
    FeynFunctionD = sd(1)
    Do i=1,6
     FeynFunctionD = FeynFunctionD + sd(i+1)*x**i
    End Do
   Else
    FeynFunctionD = ( 4._dp + x - x2/6._dp - (4._dp + 3._dp*x)*Log(r)/x ) &
                  & / (x2*x2)
   End If
  End If

 End Function FeynFunctionD
!\end{verbatim}

!\section{Function FeynFunctionB}
!\begin{verbatim}
 Real(dp) Function FeynFunctionE(r)
 !--------------------------------------------------------------------- 
 ! calculation of the Function E(r) ... 0-order in the mass
 ! input:
 !   x ............ m_fermion**2 / m_boson**2
 ! written by Werner Porod, 31.01.01
 !----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: r

  Real(dp) :: x, x2
  Integer :: i

  If (r.Eq.1._dp) Then
   FeynFunctionE = 1._dp / 30._dp
  
  Else
   x = r - 1._dp 
   x2 = x*x
   If (x2 .Lt. 1.e-3_dp) Then
    FeynFunctionE = se(1)
    Do i=1,6
     FeynFunctionE = FeynFunctionE + se(i+1)*x**i
    End Do
   Else
    FeynFunctionE = ( 4._dp * r + x2 / 3._dp   &
                  & - 2._dp * r * (1._dp + r) * Log(r) / x ) / (x2*x2)
   End If
  End If

 End Function FeynFunctionE

  Subroutine K_To_PiNuNu(gauge, mf_d, mf_u, mW, mZ,Y_d, RdL, RdR, Y_u      &
   & ,RuL, RuR, Y_l, mSneut2, Rsneut, mSlepton2, Rslepton, mSpm2, RSpm     & 
   & ,mC, U, V, mSup2, RSup, mSdown2, RSdown, mglu, phi_g, mN, N, vevSM    &
   & ,l_QCD, KtoPiNuNu, c11_o, c11p_o)
 !-------------------------------------------
 ! gives BR(K->pi nu nu) using formulas 5.10 of
 ! Bobeth et al., NPB630 (2002) 87
 ! input: all running quantities are given at Q=M_W
 !  gauge(i) ............. U(1), SU(2) and SU(3) MSbar gauge couplings
 !  mf_d(i) .............. MSbar d-quark masses
 !  mf_u(i) .............. MSbar u-quark masses
 !  mW ............... W-boson mass
 !  Y_d(i,j) ............. MSbar d-type Yukawa couplings
 !  RdL(j,j), RdR(i,j) ... mixing matrix of left- and right d-quarks 
 !  Y_u(i,j) ............. MSbar u-type Yukawa couplings
 !  RuL(j,j), RuR(i,j) ... mixing matrix of left- and right u-quarks
 !  Y_l(i,j) ............. MSbar lepton Yukawa couplings
 !  mSneut2(i) ........... sneutrino masses squared
 !  RSneut(i,j) .......... sneutrino mixing matrix
 !  mSlepton2(i) ......... slepton masses squared
 !  RSlepton(i,j) ........ slepton mixing matrix
 !  mSpm2(i) ............. mass of G+ and H+ squared
 !  RSpm(i,j) ............ mixing matrix of G+ and H+
 !  mC(i) ................ chargino masses
 !  U(i,j), V(i,j) ....... chargino mixing matrices
 !  mSup2(i) ............. u-squark masses squared
 !  RSup(i,j) ............ u-squark mixing matrix
 !  mSdown2(i) ........... d-squark masses squared
 !  RSdown(i,j) .......... d-squark mixing matrix
 !  mglu ................. gluino mass
 !  phi_g ................ phase of M_3
 !  mN(i) ................ neutralino masses
 !  N(i,j) ............... neutralino mixing matrix
 !  vevSM(i) ............. MSSM vevs (v_d, v_u)
 !  l_QCD ................ if true, then QCD corrections will be included
 !                         this is not yet completed
 ! output:
 !  BtoSNuNu ............. BR(b -> s nu nu) 
 !  Wilson coefficients at m_W (including the various contributions)
 !                              all are optional
 !  c11_o(i) .............. C_11: 1 ... total
 !                                2 ... SM contribution  
 !                                3 ... H+ contribution  
 !                                4 ... chargino contribution  
 !                                5 ... neutralino contribution  
 !                                6 ... gluino contribution  
 !  c11p_o(i) ............. C'_11: 1 ... total
 !                                 2 ... SM contribution  
 !                                 3 ... H+ contribution  
 !                                 4 ... chargino contribution  
 !                                 5 ... neutralino contribution  
 !                                 6 ... gluino contribution  
 ! written by Joel Jones, 10.05.09
 !-------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: mf_d(3), gauge(3), mf_u(3), mW, mSpm2(2), mC(2)   &
    & , mSup2(6), mSdown2(6), mglu, mN(4), vevSM(2), mSneut2(3), mSlepton2(6) &
    & , mZ
  Complex(dp), Intent(in), Dimension(3,3) :: RdL, RdR, RuL, RuR, Y_d, Y_u, Y_l
  Complex(dp), Intent(in) :: U(2,2), V(2,2), RSpm(2,2), RSup(6,6)    &
    & , Rsneut(3,3), RSdown(6,6), phi_g, N(4,4), Rslepton(6,6)
  Logical, Intent(in) :: l_QCD
  Real(dp), Intent(out) :: KToPiNuNu
  Complex(dp), Intent(out), Optional :: c11_o(3,6), c11p_o(3,6)

  Integer :: i1, i2, i3, i4
  Real(dp) :: g, gp, gs
  Complex(dp) :: cpl_CDSu_L(2,3,6), cpl_CDSu_R(2,3,6), RSup_in(6,6), mix(6,6)  &
    & , cpl_DGSd_L(3,6), cpl_DGSd_R(3,6), cpl_DNSd_L(3,4,6), cpl_DNSd_R(3,4,6)
  Complex(dp) :: c11(3), c11p(3), c11_c(6), c11p_c(6)
  Complex(dp) :: kappa_q, cpl_CNuSl_L(2,3,6), cpl_CNuSl_R(2,3,6)    &
             & , cpl_NuNSn_R(3,4,3)
  Real(dp) :: mC2(2), mN2(4), a_s, sW2, tanb, e2, mG2, mW2, a_ew
  Real(dp), Parameter :: T3=-0.5_dp, ed = -1._dp / 3._dp
  Complex(dp), Parameter :: Null3(3,3) = 0._dp
  !JJP
  Real(dp) ::mu,a_mu,mf_up(3),mf_dp(3),beta0,n_flav,mSup2p(6),mSdown2p(6)
 
  Iname = Iname + 1
  NameOfUnit(Iname) = "K_To_PiNuNu"


  !---------------
  ! couplings
  !---------------
  gp = gauge(1)
  g = gauge(2)
  gs = gauge(3)

  a_s = gauge(3)**2 * oo4pi

  CKM =  Matmul(RuL, Transpose(Conjg(RdL)) )

  cpl_CDSu_L = 0._dp
  cpl_CDSu_R = 0._dp
  If (GenerationMixing) Then
   RSup_in = RSup
  Else ! need to reorder things to be in agreement with ordering of scalar ups
   RSup_in = 0._dp
   RSup_in(1,1) = RSup(1,1)
   RSup_in(1,4) = RSup(1,2)
   RSup_in(4,1) = RSup(2,1)
   RSup_in(4,4) = RSup(2,2)
   RSup_in(2,2) = RSup(3,3)
   RSup_in(2,5) = RSup(3,4)
   RSup_in(5,2) = RSup(4,3)
   RSup_in(5,5) = RSup(4,4)
   RSup_in(3,3) = RSup(5,5)
   RSup_in(3,6) = RSup(5,6)
   RSup_in(6,3) = RSup(6,5)
   RSup_in(6,6) = RSup(6,6)
   mix = 0._dp
   mix(1:3,1:3) = rUL
   mix(4:6,4:6) = rUr
   RSup_in = Matmul(RSup_in, mix)
  End If   

  Do i1=1,3 !Maybe 1,2?
   Do i2=1,2
    Do i3=1,6
     Call CoupCharginoSfermion3(i2, i1, i3, g, -0.5_dp, RSup_in, Y_D, Y_U &
           &, RdL, RdR, U, V, cpl_CDSu_L(i2, i1, i3), cpl_CDSu_R(i2, i1, i3) )
    End Do
   End Do
  End Do

  cpl_DGSd_L = 0._dp
  cpl_DGSd_R = 0._dp
  cpl_DNSd_L = 0._dp
  cpl_DNSd_R = 0._dp
  cpl_NuNSn_R = 0._dp
  cpl_CNuSl_L = 0._dp
  cpl_CNuSl_R = 0._dp

   Do i1=1,6
    Do i3=1,3 !Maybe 1,2?
     Call CoupGluinoSquark3(gs, phi_g, i3, i1, Rsdown, RdL, RdR     &
                          & , cpl_DGSd_L(i3,i1), cpl_DGSd_R(i3,i1) )
     Do i2=1,4
      Call CoupNeutralinoSfermion3(i3, i2, i1, gp, g, T3, ed, RSdown, RdL, RdR&
                        & , Y_d, N, cpl_DNSd_L(i3,i2,i1), cpl_DNSd_R(i3,i2,i1))
     End Do
    End Do
   End Do
   Do i1=1,3
    Do i2=1,4
     Do i3=1,3  
      Call  CoupNeutralinoSneutrino3(i1, i2, i3, gp, g, N, Rsneut, id3C &
        & , cpl_NuNSn_R(i1,i2,i3))
     End Do
    End Do
    Do i2=1,2
     Do i3=1,6
     Call CoupCharginoSfermion3(i2, i1, i3, g, 0.5_dp, RSlepton, Y_l &
        & , Null3, id3C, id3C, U, V, cpl_CNuSl_L(i2,i1,i3)                  &
        & , cpl_CNuSl_R(i2,i1,i3))
     End Do
    End Do
   End Do

  mW2 = mW**2
  mG2 = mglu**2
  mC2 = mC**2
  mN2 = mN**2
  tanb = vevSM(2) / vevSM(1)

  a_s = gauge(3)**2 * oo4pi
  sW2 = gauge(1)**2 / (gauge(1)**2 + gauge(2)**2)
  e2 = 4._dp * pi * gauge(2)**2 * sW2
  kappa_q = 8._dp * Sqrt2 * G_F * e2 * CKM(3,1) * Conjg( CKM(3,2) )
  kappa_q = 1._dp / kappa_q

  c11 = 0._dp
  c11p = 0._dp

  Do i1=1,3
   Call C_11_11p(2, 1, i1, mZ, mf_d, mf_u, mf_l, mW2, mSpm2(2), mG2, mC2    &
     & ,mN2, mSneut2, mSlepton2, mSup2, mSdown2, tanb, U, V, N, RSup, RSdown &
     & , cpl_CDSu_L, cpl_CDSu_R, cpl_CNuSl_R, cpl_NuNSn_R, cpl_DNSd_L         &
     & , cpl_DNSd_R, cpl_DGSd_L, cpl_DGSd_R, a_s, kappa_q, e2, sW2, l_QCD    &
     & , c11(i1), c11p(i1), c11_c, c11p_c)
   If (Present(c11_o)) c11_o(i1,:) = c11_c
   If (Present(c11p_o)) c11p_o(i1,:) = c11p_c
  End Do

  KToPiNuNu = 0._dp
  Do i1=1,3
   KToPiNuNu = KToPiNuNu + Aimag( Conjg(CKM(3,2)) * CKM(3,1) *      &
        & ( C11(i1)+C11p(i1) ))**2
  End Do
  KToPiNuNu = 2._dp * 1.17e-4_dp * 0.944_dp * KToPiNuNu

  Iname = Iname - 1

 End Subroutine K_To_PiNuNu


  Real(dp) Function S0(x)
  Implicit None
   Real(dp), Intent(in) :: x

   S0 = 1._dp - 2.75_dp * x + 0.25_dp * x**2 &
    & - 1.5_dp * x**2 * Log(x) / (1-x)
   S0 = x *  S0 / (1 -x)**2

  End  Function S0

End Module LowEnergy
