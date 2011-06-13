Module RGEs

! load modules
 Use Control
 Use Mathematics
! load modules

! global variables
  !-------------------------------------------------------------------
  ! coefficients for gauge contribution in RGEs for gauge couplings
  !-------------------------------------------------------------------
  Real(dp), Parameter :: b_1(3) = (/ 6.6_dp, 1._dp, -3._dp /)  &
    & , b_2(3,3) = Reshape( Source = (/ 199._dp / 25._dp, 1.8_dp, 2.2_dp    &
    &                                , 5.4_dp,           25._dp,  9._dp     &
    &                                , 17.6_dp,          24._dp, 14._dp /), &
    &                       Shape = (/3, 3/) )
  Real(dp) :: Delta_b_1(3) = 0._dp, Delta_b_2(3,3) = 0._dp
  !-------------------------------------------------------------------
  ! coefficients for Yukawa contribution in RGEs for gauge couplings
  ! the matrix is revised compared to usual notation: tau,bottom,top
  !-------------------------------------------------------------------
  Real(dp), Parameter :: a_2(3,3) =     &
    &       Reshape( Source = (/  3.6_dp,  2._dp,  0._dp       &
    &                          ,  2.8_dp,  6._dp,  4._dp       &
    &                          ,  5.2_dp,  6._dp,  4._dp  /),  &
    &                Shape = (/3, 3/) )
  ! including right handed neutrinos: tau,neutrino,bottom,top
  Real(dp), Parameter :: a_2a(3,4) =     &
    &       Reshape( Source = (/  3.6_dp,  2._dp,  0._dp       &
    &                          ,  1.2_dp,  2._dp,  0._dp       &
    &                          ,  2.8_dp,  6._dp,  4._dp       &
    &                          ,  5.2_dp,  6._dp,  4._dp  /),  &
    &                Shape = (/3, 4/) )
  ! including complete 15-plet: tau, T, bottom,top, S, Z, lam_1, lam_2
  Real(dp), Parameter :: a_2b(3,8) =     &
    &       Reshape( Source = (/  3.6_dp,  2._dp,  0._dp       &
    &                          ,  5.4_dp,  7._dp,  0._dp       &
    &                          ,  2.8_dp,  6._dp,  4._dp       &
    &                          ,  5.2_dp,  6._dp,  4._dp        &
    &                          ,  4.8_dp,  0._dp,  9._dp        &
    &                          ,  2.8_dp,  6._dp,  4._dp        &
    &                          ,  5.4_dp,  7._dp,  0._dp        &
    &                          ,  5.4_dp,  7._dp,  0._dp   /),  &
    &                Shape = (/3, 8/) )
  !-------------------------------------------------------------------
  ! coefficients for contribution in RGEs for Yukawa couplings
  ! the matrix is revised compared to usual notation: tau,bottom,top
  !-------------------------------------------------------------------
  ! gauge contributions
  !-------------------------
  Real(dp), Parameter :: c1_1(3,3) =     &
    &       Reshape( Source = (/-1.8_dp,  -7._dp /15._dp ,-13._dp /15._dp  &
    &                        , -3._dp  ,  -3._dp,           -3._dp         &
    &                        ,  0._dp  , -16._dp/3._dp,   -16._dp/3._dp /) &
    &              , Shape = (/3, 3/) )
  !-------------------------
  ! Yukawa contributions
  !-------------------------
  Real(dp), Parameter :: c2_1(3,3) =     &
    &       Reshape( Source = (/  4._dp,  1._dp,  0._dp    &
    &                          ,  3._dp,  6._dp,  1._dp    &
    &                          ,  0._dp,  1._dp,  6._dp /) &
    &              , Shape = (/3, 3/) )

 Logical, Save :: TwoLoopRGE=.True.
 Real(dp) :: M2S_GUT
 Complex(dp) :: Alam_GUT, Alamp_GUT
! global variables

! private variables
! simplifies matrix multiplication in case of diagonal entries only
Logical, Private, Save :: OnlyDiagonal
! fix to completely decouple heavy states
Logical, Private, Save :: decoupling_heavy_states = .False.
! private variables

!------------------------------------------------
! needed for programs by Florian Staub
!------------------------------------------------
# ifdef SEESAWIII
! Real(dp), parameter :: sqrt0d3 = Sqrt(0.3_dp), sqrt1d2 = Sqrt(1.2_dp) &
!        &  , sqrt7d5 = Sqrt(7.5_dp), sqrt30 = Sqrt(30._dp)             &
!        &  , sqrt4o3 = Sqrt(4._dp/3._dp)
 Real(dp), parameter :: & 
         &   sqrt0d3 = 0.54772255750516611345696978280080213395274469499798_dp &
         & , sqrt1d2 = 1.0954451150103322269139395656016042679054893899960_dp  &
         & , sqrt7d5 = 2.7386127875258305672848489140040106697637234749899_dp  &
         & , sqrt30  = 5.4772255750516611345696978280080213395274469499798_dp  &
         & , sqrt4o3 = 1.1547005383792515290182975610039149112952035025403_dp

 Real(dp), Save :: NGHb3,NGHg3, NGHw3, NGHx3,NGHxb3
 Integer, save :: ThresholdCrossed = 1
 Real(dp),Parameter::id3R(3,3)= Reshape(Source=(/ 1,0,0,&
                             &                    0,1,0,&
                             &                    0,0,1 &
                             & /),shape=(/3,3/))

# endif SEESAWIII
Contains


 Subroutine CouplingsToG(gauge, yuk_l, yuk_d, yuk_u, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)    ! u-quark Yukawa couplings
  Real(dp), Intent(out) :: g1(57)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_d(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_d(i1,i2))
    g1(SumI+32) = Real(yuk_u(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_u(i1,i2))
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CouplingsToG

 Subroutine CouplingsToG2(gauge, yuk_l, yuk_nu, yuk_d, yuk_u, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_nu for MSSM
 ! written by Werner Porod, 04.03.2001
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3) & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)    ! u-quark Yukawa couplings
  Real(dp), Intent(out) :: g1(75)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG2'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CouplingsToG2

 Subroutine CouplingsToG3(gauge, yuk_l, yuk_nu, yuk_d, yuk_u, Mnu, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_nu for MSSM
 ! written by Werner Porod
 ! 30.08.02: taking CouplingsToG2 as basis, adding left neutrino mass
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3) & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)  & ! u-quark Yukawa couplings
                         & ,  Mnu(3,3)     ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: g1(93)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG3'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+68) = Real(Mnu(i1,i2),dp)
    g1(SumI+69) = Aimag(Mnu(i1,i2))
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine CouplingsToG3

 Subroutine CouplingsToG4(gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_T for MSSM
 ! written by Werner Porod, 04.03.2001
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)  & ! triplet Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)  & ! u-quark Yukawa couplings
                         & , lam1        & ! coupling triplet H_d
                         & , lam2          ! coupling triplet H_u
  Real(dp), Intent(out) :: g1(79)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG4'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
   End Do
  End Do
  g1(76) = Real(lam1,dp)
  g1(77) = Aimag(lam1)
  g1(78) = Real(lam2,dp)
  g1(79) = Aimag(lam2)

  Iname = Iname - 1

 End Subroutine CouplingsToG4

 Subroutine CouplingsToG5(gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S &
             & , lam1, lam2, M15, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings to a vector, splitting complex parameters into
 ! real and imaginary part, including Y_T for MSSM
 ! written by Werner Porod, 04.03.2001
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)   ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)  & ! triplet Yukawa couplings
                         & , yuk_d(3,3)  & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)  & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Yukawa couplings 
                         & , yuk_S(3,3)   & ! triplet Yukawa couplings
                         & , lam1        & ! coupling triplet H_d
                         & , lam2          ! coupling triplet H_u
  Real(dp), Intent(in) :: M15(3)            ! M_T, M_Z, M_S 15-plet masses
  Real(dp), Intent(out) :: g1(118)          ! vector containing these couplings

  Integer :: i1, i2, sumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'CouplingsToG5'

  g1(1:3) = gauge
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+68) = Real(yuk_Z(i1,i2),dp)
    g1(SumI+69) = Aimag(yuk_Z(i1,i2))
    g1(SumI+86) = Real(yuk_S(i1,i2),dp)
    g1(SumI+87) = Aimag(yuk_S(i1,i2))
   End Do
  End Do
  g1(112) = Real(lam1,dp)
  g1(113) = Aimag(lam1)
  g1(114) = Real(lam2,dp)
  g1(115) = Aimag(lam2)
  g1(116:118) = M15
  Iname = Iname - 1

 End Subroutine CouplingsToG5

 Function Frge3(t,gy) Result(f)
 !--------------------------------------------------------
 ! 2-loop rge's for g_3, e and m_b in the SM, valid in the
 ! range m_b<Q<m_Z
 ! written by Werner Porod, 3.1.01
 ! 25.09.01: Portation to f90
 !--------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: t  &  ! scale
                     &  , gy(3) !g_s, e, m_b
  Real(dp) :: f(3)

  Real(dp) :: g32, g34, g36, e2, e4, g32e2, beta1, beta2, beta3, Q

  Iname = Iname + 1
  NameOfUnit(Iname) = 'Frge3'

  q = t
  g32 = gy(1)**2
  g34 = gy(1)**4
  g36 = gy(1)**6
  e2 = gy(2)**2
  e4 = gy(2)**4
  g32e2 = g32 * e2 
 !--------
 ! g_3
 !--------
  beta1 = - 23._dp * g32 / 3._dp
  beta2 = 22._dp * g32e2 / 9._dp - 116._dp * g34 / 3._dp
  beta3 = - 9769._dp * g36 / 54._dp
  f(1) = oo16pi2 * gy(1) * ( beta1 + oo16pi2 * ( beta2 + oo16pi2 * beta3) )
 !--------
 ! e
 !--------
  beta1 = 80._dp * e2 / 9._dp
  beta2 = 464._dp * e4 / 27._dp + 176._dp * g32e2 / 9._dp
  f(2) = oo16pi2 * gy(2) * ( beta1 + oo16pi2 * beta2)
 !--------
 ! m_b
 !--------
  beta1 = - 8._dp * g32 - 2._dp * e2 / 3._dp
  beta2 = - 1012._dp * g34 / 9._dp - 8._dp * g32e2 / 9._dp   &
    &    + 397._dp * e4 / 81._dp
  beta3 = - 949.742_dp * g36
  f(3) = oo16pi2 * gy(3) * ( beta1 + oo16pi2 * ( beta2 + oo16pi2 * beta3) )

  Iname = Iname - 1

 End Function Frge3

 Subroutine GToCouplings(g1, gauge, yuk_l, yuk_d, yuk_u)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(57)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine GToCouplings

 Subroutine GToCouplings2(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(75)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3) &  ! neutrino Yukawa couplings
                         & , yuk_d(3,3)  &  ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine GToCouplings2

 Subroutine GToCouplings3(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u, Mnu)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 ! 30.08.02: adding left neutrino mass matrix, taking GToCouplings2 as basis
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(93)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & ,  Mnu(3,3)      ! dim-5 neutrino mass operator
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings3'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Mnu(i1,i2)  = Cmplx( g1(SumI+68), g1(SumI+69),dp )
   End Do
  End Do

  Iname = Iname - 1

 End Subroutine GToCouplings3

 Subroutine GToCouplings4(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(79)            ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings4'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
   End Do
  End Do
  lam1  = Cmplx( g1(76), g1(77),dp )
  lam2  = Cmplx( g1(78), g1(79),dp )

  Iname = Iname - 1

 End Subroutine GToCouplings4

 Subroutine GToCouplings5(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S &
             & , lam1, lam2, M15)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings, putting real and imaginary part
 ! together
 ! written by Werner Porod, 17.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(118)           ! vector containing these couplings
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Yukawa couplings 
                         & , yuk_S(3,3)   & ! triplet Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Real(dp), Intent(out) :: M15(3)           ! M_T, M_Z, M_S 15-plet masses
 
  Integer :: i1, i2, SumI

  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToCouplings5'

  gauge = g1(1:3)
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2)  = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2)  = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    yuk_Z(i1,i2)  = Cmplx( g1(SumI+68), g1(SumI+69),dp )
    yuk_S(i1,i2)  = Cmplx( g1(SumI+86), g1(SumI+87),dp )
   End Do
  End Do
  lam1  = Cmplx( g1(112), g1(113),dp )
  lam2  = Cmplx( g1(114), g1(115),dp )
  M15 = g1(116:118)

  Iname = Iname - 1

 End Subroutine GToCouplings5


 Subroutine GToParameters(g1, gauge, yuk_l, yuk_d, yuk_u  &
                         &, Mhlf, Ae, Ad, Au, Me, Ml, Md, Mq, Mu, Mh, mue, B)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(213)            ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters'

  gauge = g1(1:3)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(56 + 2*i1), g1(57 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+56), g1(SumI+57),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+74), g1(SumI+75),dp )
    Au(i1,i2) = Cmplx( g1(SumI+92), g1(SumI+93),dp )
    Me(i1,i2) = Cmplx( g1(SumI+110), g1(SumI+111),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+128), g1(SumI+129),dp )
    Md(i1,i2) = Cmplx( g1(SumI+146), g1(SumI+147),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+164), g1(SumI+165),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+182), g1(SumI+183),dp )
   End Do
  End Do
  mH = g1(208:209)
  mue = Cmplx( g1(210), g1(211),dp )
  B = Cmplx( g1(212), g1(213),dp )

  Iname = Iname - 1

 End Subroutine GToParameters

 Subroutine GToParameters2(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u  &
        &, Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(267)            ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters2'

  gauge = g1(1:3)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(74 + 2*i1), g1(75 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+74), g1(SumI+75),dp )
    Anu(i1,i2) = Cmplx( g1(SumI+92), g1(SumI+93),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+110), g1(SumI+111),dp )
    Au(i1,i2) = Cmplx( g1(SumI+128), g1(SumI+129),dp )
    Me(i1,i2) = Cmplx( g1(SumI+146), g1(SumI+147),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+164), g1(SumI+165),dp )
    Mr(i1,i2) = Cmplx( g1(SumI+182), g1(SumI+183),dp )
    Md(i1,i2) = Cmplx( g1(SumI+200), g1(SumI+201),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+218), g1(SumI+219),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+236), g1(SumI+237),dp )
   End Do
  End Do
  mH = g1(262:263)
  mue = Cmplx( g1(264), g1(265),dp )
  B = Cmplx( g1(266), g1(267),dp )

  Iname = Iname - 1

 End Subroutine GToParameters2

 Subroutine GToParameters3(g1, gauge, yuk_l, yuk_nu, yuk_d, yuk_u  &
        & , Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue &
        & , B, MnuL)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(285)            ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters3'

  gauge = g1(1:3)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(74 + 2*i1), g1(75 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_nu(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+74), g1(SumI+75),dp )
    Anu(i1,i2) = Cmplx( g1(SumI+92), g1(SumI+93),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+110), g1(SumI+111),dp )
    Au(i1,i2) = Cmplx( g1(SumI+128), g1(SumI+129),dp )
    Me(i1,i2) = Cmplx( g1(SumI+146), g1(SumI+147),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+164), g1(SumI+165),dp )
    Mr(i1,i2) = Cmplx( g1(SumI+182), g1(SumI+183),dp )
    Md(i1,i2) = Cmplx( g1(SumI+200), g1(SumI+201),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+218), g1(SumI+219),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+236), g1(SumI+237),dp )
    MnuL(i1,i2) = Cmplx( g1(SumI+260), g1(SumI+261),dp )
   End Do
  End Do
  mH = g1(262:263)
  mue = Cmplx( g1(264), g1(265),dp )
  B = Cmplx( g1(266), g1(267),dp )

  Iname = Iname - 1

 End Subroutine GToParameters3

 Subroutine GToParameters4(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2  &
         & ,  Mhlf, Ae, AT, Ad, Au, Alam1, Alam2, Me, Ml, Md, Mq, Mu, Mh  &
         & , MT, mue, B, MnuL)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(277)           ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Real(dp), Intent(out) :: MT(2)        ! soft Higgs triplet masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters4'

  gauge = g1(1:3)
  lam1 = Cmplx(g1(76),g1(77),dp)
  lam2 = Cmplx(g1(78),g1(79),dp)
  Alam1 = Cmplx(g1(158),g1(159),dp)
  Alam2 = Cmplx(g1(160),g1(161),dp)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(78 + 2*i1), g1(79 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    Ae(i1,i2) = Cmplx( g1(SumI+78), g1(SumI+79),dp )
    AT(i1,i2) = Cmplx( g1(SumI+96), g1(SumI+97),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+114), g1(SumI+115),dp )
    Au(i1,i2) = Cmplx( g1(SumI+132), g1(SumI+133),dp )
    Me(i1,i2) = Cmplx( g1(SumI+154), g1(SumI+155),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+172), g1(SumI+173),dp )
    Md(i1,i2) = Cmplx( g1(SumI+190), g1(SumI+191),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+208), g1(SumI+209),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+226), g1(SumI+227),dp )
    MnuL(i1,i2) = Cmplx( g1(SumI+252), g1(SumI+253),dp )
   End Do
  End Do
  mH = g1(252:253)
  mT = g1(254:255)
  mue = Cmplx( g1(256), g1(257),dp )
  B = Cmplx( g1(258), g1(259),dp )

  Iname = Iname - 1

 End Subroutine GToParameters4

 Subroutine GToParameters5(g1, gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S &
         & , lam1, lam2, Mhlf, Ae, AT, Ad, Au, AZ, AS, Alam1, Alam2, Me, Ml    &
         & , Md, Mq, Mu, Mh, MT, mZ, mS, MT15, MZ15, MS15, mue, B, MnuL)
 !-----------------------------------------------------------------------
 ! transform a vector to the couplings ans susy mass parameters, 
 ! putting real and imaginary part together
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 24.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g1(356)           ! vector containing the parameters
  Real(dp), Intent(out) :: gauge(3)         ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(out) :: yuk_l(3,3)  & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Z Yukawa couplings
                         & , yuk_S(3,3)   & ! triplet S Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(out) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(out) :: Ae(3,3)  & ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , AZ(3,3)   & ! triplet Z A parameters
                         & , AS(3,3)   & ! triplet S A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(out) :: Me(3,3) & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(out) :: Mh(2)        ! soft Higgs masses squared
  Real(dp), Intent(out) :: MT(2)        ! soft Higgs triplet masses squared
  Real(dp), Intent(out) :: MZ(2)        ! soft Higgs Z triplet masses squared
  Real(dp), Intent(out) :: MS(2)        ! soft Higgs S triplet masses squared
  Real(dp), Intent(out) :: MT15         ! Higgs triplet masses squared
  Real(dp), Intent(out) :: MZ15         ! Higgs Z triplet masses squared
  Real(dp), Intent(out) :: MS15         ! Higgs S triplet masses squared
  Complex(dp), Intent(out) :: mue, B    ! mu, mu*B

  Integer :: i1, i2, SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'GToParameters5'

  gauge = g1(1:3)
  lam1 = Cmplx(g1(112),g1(113),dp)
  lam2 = Cmplx(g1(114),g1(115),dp)
  Alam1 = Cmplx(g1(230),g1(231),dp)
  Alam2 = Cmplx(g1(232),g1(233),dp)
  Do i1=1,3
   Mhlf(i1) = Cmplx( g1(114 + 2*i1), g1(115 + 2*i1),dp )
   Do i2=1,3
    SumI = 6*i1+2*i2
    yuk_l(i1,i2) = Cmplx( g1(SumI-4), g1(SumI-3),dp )
    yuk_T(i1,i2) = Cmplx( g1(SumI+14), g1(SumI+15),dp )
    yuk_d(i1,i2) = Cmplx( g1(SumI+32), g1(SumI+33),dp )
    yuk_u(i1,i2) = Cmplx( g1(SumI+50), g1(SumI+51),dp )
    yuk_Z(i1,i2)  = Cmplx( g1(SumI+68), g1(SumI+69),dp )
    yuk_S(i1,i2)  = Cmplx( g1(SumI+86), g1(SumI+87),dp )

    Ae(i1,i2) = Cmplx( g1(SumI+114), g1(SumI+115),dp )
    AT(i1,i2) = Cmplx( g1(SumI+132), g1(SumI+133),dp )
    Ad(i1,i2) = Cmplx( g1(SumI+150), g1(SumI+151),dp )
    Au(i1,i2) = Cmplx( g1(SumI+168), g1(SumI+169),dp )
    AZ(i1,i2) = Cmplx( g1(SumI+186), g1(SumI+187),dp )
    AS(i1,i2) = Cmplx( g1(SumI+204), g1(SumI+205),dp )

    Me(i1,i2) = Cmplx( g1(SumI+226), g1(SumI+227),dp )
    Ml(i1,i2) = Cmplx( g1(SumI+244), g1(SumI+245),dp )
    Md(i1,i2) = Cmplx( g1(SumI+262), g1(SumI+263),dp )
    Mq(i1,i2) = Cmplx( g1(SumI+280), g1(SumI+281),dp )
    Mu(i1,i2) = Cmplx( g1(SumI+298), g1(SumI+299),dp )

    MnuL(i1,i2) = Cmplx( g1(SumI+331), g1(SumI+332),dp )
   End Do
  End Do

  mH = g1(324:325)
  mT = g1(326:327)
  mZ = g1(328:329)
  mS = g1(330:331)

  MT15 = g1(332)
  MZ15 = g1(333)
  MS15 = g1(334)
  
  mue = Cmplx( g1(335), g1(336),dp )
  B = Cmplx( g1(337), g1(338),dp )

  Iname = Iname - 1

 End Subroutine GToParameters5

 Subroutine ParametersToG(gauge, yuk_l, yuk_d, yuk_u                 &
               &,  Mhlf, Ae, Ad, Au, Me, Ml, Md, Mq, Mu, Mh, mue, B, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)   & ! lepton A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(in) :: Me(3,3)  & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(in) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(in) :: mue, B    ! mu, mu*B
  Real(dp), Intent(out) :: g1(213)            ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG'

  g1(1:3) = gauge
  Do i1=1,3
   g1(56 + 2*i1) = Real(Mhlf(i1),dp)
   g1(57 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_d(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_d(i1,i2))
    g1(SumI+32) = Real(yuk_u(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_u(i1,i2))
    g1(SumI+56) = Real(Ae(i1,i2),dp)
    g1(SumI+57) = Aimag(Ae(i1,i2))
    g1(SumI+74) = Real(Ad(i1,i2),dp)
    g1(SumI+75) = Aimag(Ad(i1,i2))
    g1(SumI+92) = Real(Au(i1,i2),dp)
    g1(SumI+93) = Aimag(Au(i1,i2))
    g1(SumI+110) = Real(Me(i1,i2),dp)
    g1(SumI+111) = Aimag(Me(i1,i2))
    g1(SumI+128) = Real(Ml(i1,i2),dp)
    g1(SumI+129) = Aimag(Ml(i1,i2))
    g1(SumI+146) = Real(Md(i1,i2),dp)
    g1(SumI+147) = Aimag(Md(i1,i2))
    g1(SumI+164) = Real(Mq(i1,i2),dp)
    g1(SumI+165) = Aimag(Mq(i1,i2))
    g1(SumI+182) = Real(Mu(i1,i2),dp)
    g1(SumI+183) = Aimag(Mu(i1,i2))
   End Do
  End Do
  g1(208) = mH(1)
  g1(209) = mH(2)
  g1(210) = Real(mue,dp)
  g1(211) = Aimag(mue)
  g1(212) = Real(B,dp)
  g1(213) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG

 Subroutine ParametersToG2(gauge, yuk_l, yuk_nu, yuk_d, yuk_u             &
         &,  Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)   & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(in) :: Me(3,3)  & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3)   ! R u-squark mass parameters squared
  Real(dp), Intent(in) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(in) :: mue, B    ! mu, mu*B
  Real(dp), Intent(out) :: g1(267)            ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG2'

  g1(1:3) = gauge
  Do i1=1,3
   g1(74 + 2*i1) = Real(Mhlf(i1),dp)
   g1(75 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+74) = Real(Ae(i1,i2),dp)
    g1(SumI+75) = Aimag(Ae(i1,i2))
    g1(SumI+92) = Real(Anu(i1,i2),dp)
    g1(SumI+93) = Aimag(Anu(i1,i2))
    g1(SumI+110) = Real(Ad(i1,i2),dp)
    g1(SumI+111) = Aimag(Ad(i1,i2))
    g1(SumI+128) = Real(Au(i1,i2),dp)
    g1(SumI+129) = Aimag(Au(i1,i2))
    g1(SumI+146) = Real(Me(i1,i2),dp)
    g1(SumI+147) = Aimag(Me(i1,i2))
    g1(SumI+164) = Real(Ml(i1,i2),dp)
    g1(SumI+165) = Aimag(Ml(i1,i2))
    g1(SumI+182) = Real(Mr(i1,i2),dp)
    g1(SumI+183) = Aimag(Mr(i1,i2))
    g1(SumI+200) = Real(Md(i1,i2),dp)
    g1(SumI+201) = Aimag(Md(i1,i2))
    g1(SumI+218) = Real(Mq(i1,i2),dp)
    g1(SumI+219) = Aimag(Mq(i1,i2))
    g1(SumI+236) = Real(Mu(i1,i2),dp)
    g1(SumI+237) = Aimag(Mu(i1,i2))
   End Do
  End Do
  g1(262) = mH(1)
  g1(263) = mH(2)
  g1(264) = Real(mue,dp)
  g1(265) = Aimag(mue)
  g1(266) = Real(B,dp)
  g1(267) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG2

 Subroutine ParametersToG3(gauge, yuk_l, yuk_nu, yuk_d, yuk_u             &
         & ,  Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B   &
         & , MnuL5, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 08.01.07: including dim5 operator for neutrinos
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_nu(3,3)  & ! neutrino Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)     ! u-quark Yukawa couplings
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)   & ! lepton A parameters
                         & , Anu(3,3)  & ! neutrino A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)     ! u-quark A parameters
  Complex(dp), Intent(in) :: Me(3,3)  & ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Mr(3,3) & ! R sneutrino mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL5(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(in) :: Mh(2)        ! soft Higgs masses squared
  Complex(dp), Intent(in) :: mue, B    ! mu, mu*B
  Real(dp), Intent(out) :: g1(285)            ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG3'

  g1(1:3) = gauge
  Do i1=1,3
   g1(74 + 2*i1) = Real(Mhlf(i1),dp)
   g1(75 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_nu(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_nu(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+74) = Real(Ae(i1,i2),dp)
    g1(SumI+75) = Aimag(Ae(i1,i2))
    g1(SumI+92) = Real(Anu(i1,i2),dp)
    g1(SumI+93) = Aimag(Anu(i1,i2))
    g1(SumI+110) = Real(Ad(i1,i2),dp)
    g1(SumI+111) = Aimag(Ad(i1,i2))
    g1(SumI+128) = Real(Au(i1,i2),dp)
    g1(SumI+129) = Aimag(Au(i1,i2))
    g1(SumI+146) = Real(Me(i1,i2),dp)
    g1(SumI+147) = Aimag(Me(i1,i2))
    g1(SumI+164) = Real(Ml(i1,i2),dp)
    g1(SumI+165) = Aimag(Ml(i1,i2))
    g1(SumI+182) = Real(Mr(i1,i2),dp)
    g1(SumI+183) = Aimag(Mr(i1,i2))
    g1(SumI+200) = Real(Md(i1,i2),dp)
    g1(SumI+201) = Aimag(Md(i1,i2))
    g1(SumI+218) = Real(Mq(i1,i2),dp)
    g1(SumI+219) = Aimag(Mq(i1,i2))
    g1(SumI+236) = Real(Mu(i1,i2),dp)
    g1(SumI+237) = Aimag(Mu(i1,i2))
    g1(SumI+260) = Real(MnuL5(i1,i2),dp)
    g1(SumI+261) = Aimag(MnuL5(i1,i2))
   End Do
  End Do
  g1(262) = mH(1)
  g1(263) = mH(2)
  g1(264) = Real(mue,dp)
  g1(265) = Aimag(mue)
  g1(266) = Real(B,dp)
  g1(267) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG3

 Subroutine ParametersToG4(gauge, yuk_l, yuk_T, yuk_d, yuk_u, lam1, lam2  &
         & ,  Mhlf, Ae, AT, Ad, Au, Alam1, Alam2, Me, Ml, Md, Mq, Mu, Mh  &
         & , MT, mue, B, MnuL5, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 08.01.07: including dim5 operator for neutrinos
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)  &  ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(in) :: Me(3,3) &  ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL5(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(in) :: Mh(2)         ! soft Higgs masses squared
  Real(dp), Intent(in) :: MT(2)         ! soft Higgs triplet masses squared
  Complex(dp), Intent(in) :: mue, B     ! mu, mu*B
  Real(dp), Intent(out) :: g1(277)           ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG4'

  g1(1:3) = gauge
  g1(76) = Real(lam1,dp)
  g1(77) = Aimag(lam1)
  g1(78) = Real(lam2,dp)
  g1(79) = Aimag(lam2)
  g1(158) = Real(Alam1,dp)
  g1(159) = Aimag(Alam1)
  g1(160) = Real(Alam2,dp)
  g1(161) = Aimag(Alam2)
  Do i1=1,3
   g1(78 + 2*i1) = Real(Mhlf(i1),dp)
   g1(79 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+78) = Real(Ae(i1,i2),dp)
    g1(SumI+79) = Aimag(Ae(i1,i2))
    g1(SumI+96) = Real(AT(i1,i2),dp)
    g1(SumI+97) = Aimag(AT(i1,i2))
    g1(SumI+114) = Real(Ad(i1,i2),dp)
    g1(SumI+115) = Aimag(Ad(i1,i2))
    g1(SumI+132) = Real(Au(i1,i2),dp)
    g1(SumI+133) = Aimag(Au(i1,i2))
    g1(SumI+154) = Real(Me(i1,i2),dp)
    g1(SumI+155) = Aimag(Me(i1,i2))
    g1(SumI+172) = Real(Ml(i1,i2),dp)
    g1(SumI+173) = Aimag(Ml(i1,i2))
    g1(SumI+190) = Real(Md(i1,i2),dp)
    g1(SumI+191) = Aimag(Md(i1,i2))
    g1(SumI+208) = Real(Mq(i1,i2),dp)
    g1(SumI+209) = Aimag(Mq(i1,i2))
    g1(SumI+226) = Real(Mu(i1,i2),dp)
    g1(SumI+227) = Aimag(Mu(i1,i2))
    g1(SumI+252) = Real(MnuL5(i1,i2),dp)
    g1(SumI+253) = Aimag(MnuL5(i1,i2))
   End Do
  End Do
  g1(252) = mH(1)
  g1(253) = mH(2)
  g1(254) = mT(1)
  g1(255) = mT(2)
  g1(256) = Real(mue,dp)
  g1(257) = Aimag(mue)
  g1(258) = Real(B,dp)
  g1(259) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG4

 Subroutine ParametersToG5(gauge, yuk_l, yuk_T, yuk_d, yuk_u, yuk_Z, yuk_S  &
         & , lam1, lam2, Mhlf, Ae, AT, Ad, Au, AZ, AS, Alam1, Alam2, Me, Ml &
         & , Md, Mq, Mu, Mh, MT, MZ, MS, MT15, MZ15, MS15, mue, B, MnuL5, g1)
 !-----------------------------------------------------------------------
 ! transform the couplings and susy mass parameters to a vector, 
 ! splitting complex parameters into real and imaginary part
 ! written by Werner Porod, 20.8.99 
 ! 07.03.2001: including right handed neutrinos
 ! 08.01.07: including dim5 operator for neutrinos
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: gauge(3)          ! gauge couplings U(1), SU(2), SU(3)
  Complex(dp), Intent(in) :: yuk_l(3,3)   & ! lepton Yukawa couplings
                         & , yuk_T(3,3)   & ! triplet Yukawa couplings
                         & , yuk_d(3,3)   & ! d-quark Yukawa couplings
                         & , yuk_u(3,3)   & ! u-quark Yukawa couplings
                         & , yuk_Z(3,3)   & ! triplet Z Yukawa couplings
                         & , yuk_S(3,3)   & ! triplet S Yukawa couplings
                         & , lam1         & ! coupling triplet H_d
                         & , lam2           ! coupling triplet H_u
  Complex(dp), Intent(in) :: Mhlf(3)       ! gaugino masses
  Complex(dp), Intent(in) :: Ae(3,3)  &  ! lepton A parameters
                         & , AT(3,3)   & ! triplet A parameters
                         & , Ad(3,3)   & ! d-quark A parameters
                         & , Au(3,3)   & ! u-quark A parameters
                         & , AZ(3,3)   & ! triplet Z A parameters
                         & , AS(3,3)   & ! triplet S A parameters
                         & , Alam1     & ! A-parameter triplet H_d
                         & , Alam2       ! A-parameter triplet H_u
  Complex(dp), Intent(in) :: Me(3,3) &  ! R-slepton mass parameters squared
                         &  , Ml(3,3) & ! L-slepton mass parameters squared
                         &  , Md(3,3) & ! R d-squark mass parameters squared
                         &  , Mq(3,3) & ! L squark mass parameters squared
                         &  , Mu(3,3) & ! R u-squark mass parameters squared
                         &  , MnuL5(3,3) ! dim-5 neutrino mass operator
  Real(dp), Intent(in) :: Mh(2)         ! soft Higgs masses squared
  Real(dp), Intent(in) :: MT(2)         ! soft Higgs triplet masses squared
  Real(dp), Intent(in) :: MZ(2)         ! soft Higgs Z triplet masses squared
  Real(dp), Intent(in) :: MS(2)         ! soft Higgs S triplet masses squared
  Real(dp), Intent(in) :: MT15          ! Higgs triplet masses squared
  Real(dp), Intent(in) :: MZ15          ! Higgs Z triplet masses squared
  Real(dp), Intent(in) :: MS15          ! Higgs S triplet masses squared
  Complex(dp), Intent(in) :: mue, B     ! mu, mu*B
  Real(dp), Intent(out) :: g1(356)           ! vector containing the parameters

  Integer i1,i2,SumI
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'ParametersToG5'

  g1(1:3) = gauge
  g1(112) = Real(lam1,dp)
  g1(113) = Aimag(lam1)
  g1(114) = Real(lam2,dp)
  g1(115) = Aimag(lam2)

  g1(230) = Real(Alam1,dp)
  g1(231) = Aimag(Alam1)
  g1(232) = Real(Alam2,dp)
  g1(233) = Aimag(Alam2)
  Do i1=1,3
   g1(114 + 2*i1) = Real(Mhlf(i1),dp)
   g1(115 + 2*i1) = Aimag(Mhlf(i1))
   Do i2=1,3
    SumI = 6*i1+2*i2
    g1(SumI-4) = Real(yuk_l(i1,i2),dp)
    g1(SumI-3) = Aimag(yuk_l(i1,i2))
    g1(SumI+14) = Real(yuk_T(i1,i2),dp)
    g1(SumI+15) = Aimag(yuk_T(i1,i2))
    g1(SumI+32) = Real(yuk_d(i1,i2),dp)
    g1(SumI+33) = Aimag(yuk_d(i1,i2))
    g1(SumI+50) = Real(yuk_u(i1,i2),dp)
    g1(SumI+51) = Aimag(yuk_u(i1,i2))
    g1(SumI+68) = Real(yuk_Z(i1,i2),dp)
    g1(SumI+69) = Aimag(yuk_Z(i1,i2))
    g1(SumI+86) = Real(yuk_S(i1,i2),dp)
    g1(SumI+87) = Aimag(yuk_S(i1,i2))

    g1(SumI+114) = Real(Ae(i1,i2),dp)
    g1(SumI+115) = Aimag(Ae(i1,i2))
    g1(SumI+132) = Real(AT(i1,i2),dp)
    g1(SumI+133) = Aimag(AT(i1,i2))
    g1(SumI+150) = Real(Ad(i1,i2),dp)
    g1(SumI+151) = Aimag(Ad(i1,i2))
    g1(SumI+168) = Real(Au(i1,i2),dp)
    g1(SumI+169) = Aimag(Au(i1,i2))
    g1(SumI+186) = Real(AZ(i1,i2),dp)
    g1(SumI+187) = Aimag(AZ(i1,i2))
    g1(SumI+204) = Real(AS(i1,i2),dp)
    g1(SumI+205) = Aimag(AS(i1,i2))

    g1(SumI+226) = Real(Me(i1,i2),dp)
    g1(SumI+227) = Aimag(Me(i1,i2))
    g1(SumI+244) = Real(Ml(i1,i2),dp)
    g1(SumI+245) = Aimag(Ml(i1,i2))
    g1(SumI+262) = Real(Md(i1,i2),dp)
    g1(SumI+263) = Aimag(Md(i1,i2))
    g1(SumI+280) = Real(Mq(i1,i2),dp)
    g1(SumI+281) = Aimag(Mq(i1,i2))
    g1(SumI+298) = Real(Mu(i1,i2),dp)
    g1(SumI+299) = Aimag(Mu(i1,i2))

    g1(SumI+331) = Real(MnuL5(i1,i2),dp)
    g1(SumI+332) = Aimag(MnuL5(i1,i2))
   End Do
  End Do
  g1(324:325) = mH
  g1(326:327) = mT
  g1(328:329) = mZ
  g1(330:331) = mS
  g1(332) = MT15
  g1(333) = MZ15
  g1(334) = MS15
  g1(335) = Real(mue,dp)
  g1(336) = Aimag(mue)
  g1(337) = Real(B,dp)
  g1(338) = Aimag(B)

  Iname = Iname - 1

 End Subroutine ParametersToG5

 Subroutine rge3(len, t,gy,f)
 !--------------------------------------------------------
 ! 2-loop rge's for g_3, e and m_b in the SM, valid in the
 ! range m_b<Q<m_Z
 ! written by Werner Porod, 3.1.01
 ! 25.09.01: Portation to f90
 !--------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: t, gy(len)
  Real(dp), Intent(out) :: f(len)

  Real(dp) :: g32, g34, g36, e2, e4, g32e2, beta1, beta2, beta3, q

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge3'

  q = t

  g32 = gy(1)**2
  g34 = gy(1)**4
  g36 = gy(1)**6
  e2 = gy(2)**2
  e4 = gy(2)**4
  g32e2 = g32 * e2 
 !--------
 ! g_3
 !--------
  beta1 = - 23._dp * g32 / 3._dp
  beta2 = 22._dp * g32e2 / 9._dp - 116._dp * g34 / 3._dp
  beta3 = - 9769._dp * g36 / 54._dp
  f(1) = oo16pi2 * gy(1) * ( beta1 + oo16pi2 * ( beta2 + oo16pi2 * beta3) )
 !--------
 ! e
 !--------
  beta1 = 80._dp * e2 / 9._dp
  beta2 = 464._dp * e4 / 27._dp + 176._dp * g32e2 / 9._dp
  f(2) = oo16pi2 * gy(2) * ( beta1 + oo16pi2 * beta2)
 !--------
 ! m_b
 !--------
  beta1 = - 8._dp * g32 - 2._dp * e2 / 3._dp
  beta2 = - 1012._dp * g34 / 9._dp - 8._dp * g32e2 / 9._dp   &
    &    + 397._dp * e4 / 81._dp
  beta3 = - 949.742_dp * g36
  f(3) = oo16pi2 * gy(3) * ( beta1 + oo16pi2 * ( beta2 + oo16pi2 * beta3) )

  Iname = Iname - 1

 End Subroutine rge3


 Subroutine rge6(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  6.10.2000: changing to f90
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T,GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(len),sumI,beta2(3), q

  q = t

  gy2 = gy**2

  Do i=1,3       ! gauge couplings two loop
   sumI = 0._dp
   Do j=1,3
    sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
   Enddo
   f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
  Enddo

  beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)   &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2)            &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

  beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp    &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

  beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) + 136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do

  End Subroutine rge6


 Subroutine rge7(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  10.01.00: including tan(beta)
 !  6.10.2000: changing to f90
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T, GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(6), sumI, beta2(3), gamma1, gamma2, q

  q = t

  gy2 = gy(1:6)**2

  If (TwoLoopRGE) Then
   Do i=1,3       ! gauge couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
    Enddo
    f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
   Enddo

   beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)  &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2 )           &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

   beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp   &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

   beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) +136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do
  !---------------
  ! Ln(tan(beta))
  !---------------
   gamma1 = 3._dp * (gy2(5) - gy2(6)) + gy2(4)

   gamma2 = 0.75_dp * ( 3._dp * (gy2(6)**2 - gy2(5)**2) - gy2(4)**2)  &
   &    - (1.9_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(6)  &
   &   + (0.4_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(5)  &
   &   + (1.8_dp * gy2(1) + 1.5_dp * gy2(2) ) * gy2(4)

   f(7) = oo16pi2 * (gamma1 + oo16pi2 * gamma2 )

  Else ! Everything at 1-loop
   f(1:3) = oo16pi2 * gy(1:3) * gy2(1:3) * b_1   ! gauge couplings
   Do i=1,3       ! yukawa couplings one loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * sumI 
   End Do
   f(7) = oo16pi2 * ( 3._dp * (gy2(5) - gy2(6)) + gy2(4) )
  End If

  End Subroutine rge7


 Subroutine rge8_NMSSM(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  6.10.2000: changing to f90
 !  4.12.2008: extension to include NMSSM couplings
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T,GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(len),sumI,beta2(3), q

  q = t

  gy2 = gy**2

  Do i=1,3       ! gauge couplings two loop
   sumI = 0._dp
   Do j=1,3
    sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
   Enddo
   f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
  Enddo

  beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)   &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2)            &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

  beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp    &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

  beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) + 136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    sumI = sumI + gy2(7) ! adding lambda**2
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do
   f(7) = oo16pi2 * gy(7) * (-gy2(1) - 3._dp * gy2(2) + gy2(4)         &
        &                   + 3._dp * (gy2(5)+gy2(6)) + 4._dp * gy2(7) &
        &                   + 2._dp* gy2(8) )
   f(8) = oo16pi2 * gy(8) * (  gy2(7) +  gy2(8) )

  End Subroutine rge8_NMSSM


 Subroutine rge9_NMSSM(len, T,GY,F)
 !-----------------------------------------------------------------
 !     Right hand side of truncated renormalization group equations
 !          dGY_i/dT = F_i(G)
 !     for the determination of M_GUT and the value of alpha_GUT
 !     and values of the Yukawas
 !  written by Werner Porod, 28.12.99
 !  10.01.00: including tan(beta)
 !  6.10.2000: changing to f90
 !  4.12.2008: extension to include NMSSM couplings
 !-----------------------------------------------------------------
 Implicit None
  Integer, Intent(in) :: len
  Real(Dp), Intent(in) :: T, GY(len)
  Real(Dp), Intent(out) :: F(len)

  Integer :: j,i
  Real(Dp) :: gy2(8), sumI, beta2(3), gamma1, gamma2, q

  q = t

  gy2 = gy(1:8)**2

  If (TwoLoopRGE) Then
   Do i=1,3       ! gauge couplings two loop
    sumI = 0._dp
    Do j=1,3
     sumI = sumI + b_2(i,j) * gy2(j) - a_2(i,j) * gy2(3+j)
    Enddo
    f(i) = oo16pi2 * gy(i) * gy2(i) * ( b_1(i) + oo16pi2 * sumI )
   Enddo

   beta2(1) = ( 13.5_dp* gy2(1) + 1.8_dp * gy2(2) + 1.2_dp * gy2(4)  &
   &         - 0.4_dp * gy2(5) ) * gy2(1)                            &
   &       + (7.5_dp * gy2(2) + 6._dp * gy2(4) ) * gy2(2 )           &
   &       + 16._dp * gy2(5) * gy2(3)                                &
   &       - (10._dp * gy2(4) + 9._dp * gy2(5) ) * gy2(4)            &
   &       - (9._dp * gy2(5) +  3._dp * gy2(6) ) * gy2(5)

   beta2(2) = ( 287._dp * gy2(1) / 90._dp + gy2(2) + 8._dp * gy2(3) / 9._dp   &
   &         + 1.2_dp * gy2(4) + 0.4_dp * gy2(5) + 0.8_dp * gy2(6) ) * gy2(1) &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(5) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(5) ) * gy2(3)                    &
   &       - 3._dp * (gy2(4) + gy2(5) ) * gy2(4)                              &
   &       - (22._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)                     &
   &       - 5._dp * gy2(6)**2

   beta2(3) = ( 2743._dp * gy2(1) / 450._dp + gy2(2) +136._dp * gy2(3)/45._dp &
   &         + 0.4_dp * gy2(5) + 1.2_dp * gy2(6) ) * gy2(1)                   &
   &       + (7.5_dp * gy2(2) + 8._dp * gy2(3) + 6._dp * gy2(6) ) * gy2(2)    &
   &       + 16._dp * (-gy2(3) / 9._dp + gy2(6) ) * gy2(3)                    &
   &       - (gy2(4) + 5._dp * gy2(5) + 5._dp * gy2(6) ) * gy2(5)             &
   &       - 22._dp * gy2(6)**2

   Do i=1,3       ! yukawa couplings two loop
    sumI = gy2(7)
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * (sumI + oo16pi2 * beta2(i) )
   End Do
  !---------------
  ! Ln(tan(beta))
  !---------------
   gamma1 = 3._dp * (gy2(5) - gy2(6)) + gy2(4)

   gamma2 = 0.75_dp * ( 3._dp * (gy2(6)**2 - gy2(5)**2) - gy2(4)**2)  &
   &    - (1.9_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(6)  &
   &   + (0.4_dp * gy2(1) + 4.5_dp * gy2(2) + 20._dp * gy2(3) ) * gy2(5)  &
   &   + (1.8_dp * gy2(1) + 1.5_dp * gy2(2) ) * gy2(4)

   f(7) = oo16pi2 * gy(7) * (-gy2(1) - 3._dp * gy2(2) + gy2(4)         &
        &                   + 3._dp * (gy2(5)+gy2(6)) + 4._dp * gy2(7) &
        &                   + 2._dp* gy2(8) )
   f(8) = oo16pi2 * gy(8) * (  gy2(7) +  gy2(8) )
   f(9) = oo16pi2 * (gamma1 + oo16pi2 * gamma2 )

  Else ! Everything at 1-loop

   f(1:3) = oo16pi2 * gy(1:3) * gy2(1:3) * b_1   ! gauge couplings
   Do i=1,3       ! yukawa couplings one loop
    sumI = gy2(7)
    Do j=1,3
     sumI = sumI + c1_1(i,j) * gy2(j)  + c2_1(j,i) * gy2(j+3) 
    End Do
    f(i+3) = oo16pi2 * gy(i+3) * sumI 
   End Do
   f(7) = oo16pi2 * gy(7) * (-gy2(1) - 3._dp * gy2(2) + gy2(4)         &
        &                   + 3._dp * (gy2(5)+gy2(6)) + 4._dp * gy2(7) &
        &                   + 2._dp* gy2(8) )
   f(8) = oo16pi2 * gy(8) * (  gy2(7) +  gy2(8) )
   f(9) = oo16pi2 * ( 3._dp * (gy2(5) - gy2(6)) + gy2(4) )
  End If

  End Subroutine rge9_NMSSM


 Subroutine RGE10_SM(len,t,gy,f)
 !--------------------------------------------------------
 ! RGEs within the SM assuming the MSbar scheme
 ! 2-loop RGEs for e
 ! 4-loop RGEs for g_3
 ! 2-loop RGEs for lepton masses
 ! 4-loop QCD and 2-loop QED RGES for quark masses
 ! Assumption: the only threhold to be checked is m_b
 ! input: t = Log(Q^2)
 !        gy(i) ... i=1  -> e(Q)
 !                  i=2  -> g_3
 !                  i=3  -> m_e
 !                  i=4  -> m_mu
 !                  i=5  -> m_tau
 !                  i=6  -> m_u
 !                  i=7  -> m_c
 !                  i=8  -> m_d
 !                  i=9  -> m_s
 !                  i=10 -> m_b, is optional
 ! output:
 !   f = d(gy)/d(t)
 ! written by Werner Porod, 03.12.03
 !--------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: t, gy(len)
  Real(dp), Intent(out) :: f(len)

  Integer :: i1
  Real(dp) :: g32, g34, g36, g38, e2, e4, g32e2, q
  Real(dp), Parameter :: b_e1(2) = (/ 76._dp / 9._dp , 80._dp / 9._dp /)    &
       & , b_e2(2) = (/ 460._dp / 27._dp , 464._dp / 27._dp /)              & 
       & , b_e3(2) = (/ 160._dp / 9._dp , 176._dp / 9._dp  /)               & 
       & , b_g1(2) = (/ -25._dp / 3._dp, -23._dp/3._dp /)                   &
       & , b_g2(2) = (/ -154._dp / 3._dp, -116._dp/3._dp /)                 &
       & , b_g3(2) = (/ 20._dp / 3._dp, 22._dp/3._dp /)                     &
       & , b_g4(2) = (/ -21943._dp/54._dp, 9769._dp/54._dp /)               &
       & , b_g5(2) = (/ -4918247._dp/1458._dp-414140._dp*zeta3/81._dp       &
       &             , 598391._dp/1458._dp - 352864._dp*zeta3/81._dp /)     &
       & , g_el1(2) = (/ -6._dp, -6._dp /)                                  &
       & , g_el2(2) = (/ 353._dp / 9._dp,  373._dp / 9._dp /)               & 
       & , g_eu1(2) = (/ -8._dp/3._dp, -8._dp/3._dp /)                      &
       & , g_eu2(2) = (/ 1472._dp / 81._dp, 1552._dp / 81._dp/)             & 
       & , g_eu3(2) = (/ -32._dp / 9._dp,  -32._dp / 9._dp/)                & 
       & , g_ed1(2) = (/ -2._dp/3._dp, -2._dp/3._dp /)                      &
       & , g_ed2(2) = (/ 377._dp / 81._dp,  397._dp / 81._dp /)             & 
       & , g_ed3(2) = (/ -8._dp / 9._dp,  -8._dp / 9._dp /)                 & 
       & , g_q1(2) = (/ - 8._dp , -8._dp /)                                 &
       & , g_q2(2) = (/ -1052._dp / 9._dp ,  -1012._dp / 9._dp /)           &
       & , g_q3(2) = (/ -144674._dp/81._dp + 1280._dp * zeta3 / 3._dp       &
       &              , -128858._dp/81._dp + 1600._dp * zeta3 / 3._dp /)    &
       & , g_q4(2) = (/ -7330357._dp/243._dp + 51584._dp* zeta3/3._dp       &
       &                - 16000._dp*zeta4 / 3._dp + 11200._dp* zeta5 /9._dp &
       &             , -1911065._dp/81._dp + 618400._dp* zeta3/27._dp       &
       &                - 18400._dp*zeta4 / 3._dp - 25600._dp* zeta5 /9._dp  /)

       
  Iname = Iname + 1
  NameOfUnit(Iname) = 'RGE10_SM'

  q = t

  If (len.Eq.9) Then ! check which beta function (anomalous dimension) to use
   i1 = 1
  Else If (len.Eq.10) Then
   i1 = 2
  Else
   Write(ErrCan,*) "Error in routine "//Trim(NameOfUnit(Iname))
   Write(ErrCan,*) "Length of the vector gy = ",len
   Call TerminateProgram
  End If

  g32 = gy(1)**2
  g34 = gy(1)**4
  g36 = gy(1)**6
  g38 = gy(1)**8
  e2 = gy(2)**2
  e4 = gy(2)**4
  g32e2 = g32 * e2 
 !--------
 ! g_3
 !--------
  f(1) = oo16pi2 * gy(1) * ( b_g1(i1)*g32                                     &
       &                   + oo16pi2 * ( b_g2(i1)*g34 + b_g3(i1)*g32e2        &
       &                               + oo16pi2 * ( b_g4(i1)*g36             &
       &                                           + oo16pi2 * b_g5(i1)*g38 )))
 !--------
 ! e
 !--------
  f(2) = oo16pi2 * gy(2) * ( b_e1(i1) * e2                                &
       &                   + oo16pi2 * (b_e2(i1) * e4 + b_e3(i1) * g32e2 ))
 !-----------------
 ! m_l, l=e,mu,tau
 !-----------------
  f(3:5) =  oo16pi2 * gy(3:5) * (g_el1(i1) * e2 + oo16pi2 *g_el2(i1) * e4)
 !---------
 ! m_u, m_c
 !---------
  f(6:7) = oo16pi2 * gy(6:7) * (g_eu1(i1) * e2 + g_q1(i1) * g32              &
         &                     + oo16pi2 * (g_eu2(i1)*e4 + g_eu3(i1) * g32e2 &
         &                                 + g_q2(i1) * g34                  &
         &                                 + oo16pi2 * (g_q3(i1) * g36       &
         &                                       + oo16pi2 * g_q4(i1) * g38 )))
 !---------------
 ! m_d, m_s, m_b
 !---------------
  f(8:len) = oo16pi2 * gy(8:len) * (g_ed1(i1) * e2 + g_q1(i1) * g32          &
         &                     + oo16pi2 * (g_ed2(i1)*e4 + g_ed3(i1) * g32e2 &
         &                                 + g_q2(i1) * g34                  &
         &                                 + oo16pi2 * (g_q3(i1) * g36       &
         &                                       + oo16pi2 * g_q4(i1) * g38 )))

  Iname = Iname - 1

 End Subroutine RGE10_SM


 Subroutine rge57(len,T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(3), Dgauge(3), q, TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge57'

  q = t

  OnlyDiagonal = .Not.GenerationMixing

  Call GToCouplings(gy,gauge,Ye,Yd,Yu)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYeYe(i1,i1) = Real(aYeYe(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
  End Do

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(2) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = (3._dp,0._dp) * TraceY(2) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3)              &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   !------------------------------------------------
   ! these are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdYdaYdYd(i1,i1) = Real(aYdYdaYdYd(i1,i1), dp)
    aYeYeaYeYe(i1,i1) = Real(aYeYeaYeYe(i1,i1), dp)
    aYuYuaYuYu(i1,i1) = Real(aYuYuaYuYu(i1,i1), dp)
   End Do

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(2)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(2) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(3)
   hd(2) = 9._dp * TraceY(2) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(3)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(3) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(2) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                  &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
  End If

  !--------------------------------------------
  ! This helps avoiding numerical instabilities
  !--------------------------------------------
!  Call Chop(DYe)
!  Call Chop(DYd)
!  Call Chop(DYu)

  Call CouplingsToG(Dgauge,DYe,DYd,DYu,f)

  Iname = Iname - 1

 End Subroutine rge57


 Subroutine rge59_SM(len,T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  !-------------------------------------------------------------------
  ! coefficients for gauge contribution in RGEs for gauge couplings
  !-------------------------------------------------------------------
  Real(dp), Parameter :: b_1(3) = (/ 4.1_dp, -19._dp/6._dp, -7._dp /)  &
    & , b_2(3,3) = Reshape( Source = (/ 199._dp / 50._dp, 0.9_dp, 1.1_dp    &
    &                                , 2.7_dp,   35._dp/5._dp ,  4.5_dp     &
    &                                , 8.8_dp,          12._dp, -26._dp /), &
    &                       Shape = (/3, 3/) )
  !-------------------------------------------------------------------
  ! coefficients for Yukawa contribution in RGEs for gauge couplings
  ! the matrix is revised compared to usual notation: tau,bottom,top
  !-------------------------------------------------------------------
  Real(dp), Parameter :: a_2(3,3) =     &
    &       Reshape( Source = (/  1.5_dp,  0.5_dp,  0._dp       &
    &                          ,  0.5_dp,  1.5_dp,  2._dp       &
    &                          ,  1.7_dp,  1.5_dp,  2._dp  /),  &
    &                Shape = (/3, 3/) )

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(3), Dgauge(3), q, TraceY2(4)      &
    & , lam, ln_v, lam2, Dlam, Dln_v, Y2, HS, Y4, chi4, gauge4(3), betaLam1 &
    & , betaLam2, beta_v1, beta_v2
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(3,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge59_SM'

  q = t

  OnlyDiagonal = .Not.GenerationMixing

  Call GToCouplings(gy(1:57),gauge,Ye,Yd,Yu)
  lam = gy(58)
  ln_v = gy(59)

  gauge2 = gauge**2
  gauge4 = gauge2**2
  lam2 = lam**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYeYe(i1,i1) = Real(aYeYe(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
  End Do

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  Y2 = TraceY(1) + 3._dp * (TraceY(2)+TraceY(3))
  diagonal(1:3,1) = Y2 - 2.25_dp * gauge2(2)

  diagonal(1,1) =   diagonal(1,1) - 2.25_dp * gauge2(1)
  sume1 = 1.5_dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = diagonal(2,1) - 0.25_dp * gauge2(1) - 8._dp * gauge2(3)
  sumd1  = 1.5_dp * (aYdYd - aYuYu)
  sumu1  =  - sumd1
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(3,1) = diagonal(3,1) - 0.85_dp * gauge2(1) - 8._dp * gauge2(3)
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  HS = TraceY(1)**2 + 3._dp * (TraceY(2)**2+TraceY(3)**2)

  betaLam1 = 12._dp * lam2 - 4._dp * HS                           &
         & - (1.8_dp * gauge2(1)+9._dp*gauge2(2))*lam             &
         & + 0.27_dp * gauge4(1) + 0.9_dp * gauge2(1) * gauge2(2) &
         & + 2.25_dp * gauge4(2)

  beta_v1 = 0.45_dp * gauge2(1) + 2.25_dp * gauge2(2) - Y2

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   !------------------------------------------------
   ! these are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdYdaYdYd(i1,i1) = Real(aYdYdaYdYd(i1,i1), dp)
    aYeYeaYeYe(i1,i1) = Real(aYeYeaYeYe(i1,i1), dp)
    aYuYuaYuYu(i1,i1) = Real(aYuYuaYuYu(i1,i1), dp)
   End Do

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   chi4 = 6.75_dp * (TraceY2(2)+TraceY2(3)) + 2.25_dp * TraceY2(1) &
        & - 1.5_dp * TraceY2(4)

   Y4 = (0.85_dp * gauge2(1) + 2.25_dp * gauge2(2) + 8._dp * gauge2(3)) &
    &    * TraceY(3)                                                    &
    & + (0.25_dp * gauge2(1) + 2.25_dp * gauge2(2) + 8._dp * gauge2(3)) &
    &    * TraceY(2)                                                    &
    & + 0.75_dp * (gauge2(1) + gauge2(2)) * TraceY(1)

   diagonal(2,1) = - chi4 + 1.5_dp * lam2 + 2.5_dp * Y4                       &
             &   + ( -5.75_dp * gauge2(2) + 1.35_dp * gauge2(1) ) * gauge2(2) &
             &   + 6.855_dp * gauge4(1)
   hd(1) = 2.25_dp *Y2 + 6._dp * lam &
       & - 4.8375_dp * gauge2(1) - 9._dp * gauge2(2)
   sume2 = 1.5_dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(2,1)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = - chi4 + 1.5_dp * lam2 + 2.5_dp * Y4      &
             &   - 5.75_dp * gauge4(2) - 108._dp * gauge4(3) &
             & + 9._dp * gauge2(2) * gauge2(3)
   diagonal(3,2) = diagonal(2,2)
   diagonal(2,2) = diagonal(2,2) - 1.27_dp * gauge4(1) / 6._dp &
     & + (31._dp * gauge2(3) / 15._dp - 1.35_dp * gauge2(2)) * gauge2(1)
   
   hd(1) = -2.25_dp * Y2 - 6._dp * lam &
       & + 8.4375_dp * gauge2(2) + 16._dp * gauge2(3)
   hd(2) =  1.25_dp * Y2 - 2._dp * lam   &
       & + 0.5625_dp * gauge2(2) - 16._dp * gauge2(3)
   sumd2 = 1.5_dp * aYdYdaYdYd - aYdYdaYuYu - 0.25_dp * aYuYuaYdYd         &
       & + 2.75_dp * aYuYuaYuYu + (hd(1) + 2.3375_dp * gauge2(1) ) * aYdYd &
       & + (hd(2) - 0.9875_dp * gauge2(1)) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = diagonal(3,2) + 12.41_dp * gauge4(1) / 18._dp &
     & + (19._dp * gauge2(3) / 15._dp - 0.45_dp * gauge2(2)) * gauge2(1)
   
   sumd2 = 1.5_dp * aYuYuaYuYu - aYuYuaYdYd - 0.25_dp * aYdYdaYuYu         &
       & + 2.75_dp * aYdYdaYdYd + (hd(2) - 0.5375_dp * gauge2(1) ) * aYdYd &
       & + (hd(1) + 2.7875_dp * gauge2(1)) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
   betaLam2 = (10.8_dp * gauge2(1) + 54._dp * gauge2(2) -78._dp * lam) * lam2 &
     & - (9.125_dp * gauge4(2) - 5.85_dp  * gauge2(2) * gauge2(1)             &
     &   + 2.661_dp *  gauge4(1) ) * lam                                      &
     & + (38.125_dp * gauge2(2) - 7.225_dp * gauge2(1)) * gauge4(2)           &
     & - (8.385_dp * gauge2(2) + 3.411_dp * gauge2(1)) * gauge4(1)            &
     & - 64._dp * gauge2(3) * (TraceY2(2) + TraceY2(3))                       &
     & - 1.8_dp * gauge2(1)                                                   &
     &          * ( 3._dp * TraceY2(1)-TraceY2(2) + 2._dp * TraceY2(3) )      &
     & - 1.5_dp * gauge4(2) * Y4 - 24._dp * lam2 * Y2 - lam * HS              &
     & + ( 7.5_dp * (gauge2(1)+gauge2(2)) * TraceY(1)                         &
     &   + (2.5_dp * gauge2(1) + 22.5_dp * gauge2(2) + 80._dp * gauge2(3) )   &
     &     * TraceY(2)                                                        &
     &   + (8.5_dp * gauge2(1) + 22.5_dp * gauge2(2) + 80._dp * gauge2(3) )   &
     &     * TraceY(3) ) * lam                                                &
     & + ( (-7.5_dp * gauge2(1) + 11._dp * gauge2(2)) * TraceY(1)             &
     &   + ( 1.5_dp * gauge2(1) + 9._dp * gauge2(2)) * TraceY(2)              &
     &   + (-5.7_dp * gauge2(1) + 21._dp * gauge2(2)) * TraceY(3)             &
     &   ) * 0.6_dp * gauge2(1)                                               &
     & - 12._dp * ( TraceY2(2) + TraceY2(3) + 2._dp * TraceY2(4) )            &
     & + 6._dp * lam  * TraceY2(4)                                            &
     & + 20._dp * ( 3._dp * cTrace(MatMul2(aYuYuaYuYu,aYuYu,OnlyDiagonal) )   &
     &            + 3._dp * cTrace(MatMul2(aYdYdaYdYd,aYdYd,OnlyDiagonal) )   &
     &            + cTrace(MatMul2(aYeYeaYeYe,aYeYe,OnlyDiagonal) ) )

   beta_v2 = chi4 - 1.5_dp * lam2 - 2.5_dp * Y4 - 1.61625_dp * gauge4(1) &
     & - 0.3375_dp * gauge2(2) * gauge2(1) + 8.46875_dp * gauge4(2)
     
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                  &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
  !--------------
  ! Higgs sector  
  !--------------
   Dlam = oo16pi2 * ( betaLam1 + oo16pi2 * betaLam2)
   Dln_v = oo16pi2 * ( beta_v1 + oo16pi2 * beta_v2)

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
  !--------------
  ! Higgs sector 
  !--------------
   Dlam = oo16pi2 * betaLam1
   Dln_v = oo16pi2 * beta_v1
  End If

  !--------------------------------------------
  ! This helps avoiding numerical instabilities
  !--------------------------------------------
!  Call Chop(DYe)
!  Call Chop(DYd)
!  Call Chop(DYu)

  Call CouplingsToG(Dgauge,DYe,DYd,DYu,f(1:57))
  f(58) = Dlam
  f(59) = Dln_v

  Iname = Iname - 1

 End Subroutine rge59_SM


 Subroutine rge75(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 04.03.2001: including neutrino Yukawas as given in the MSSM
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), q
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), TraceY2(4), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)             &
    & ,Ynu(3,3), aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3), DYnu(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge75'

  q = t

  OnlyDiagonal = .Not.GenerationMixing

  Call GToCouplings2(gy,gauge,Ye,Ynu,Yd,Yu)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYnuYnu = Matmul2(aYnu,Ynu,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = (3._dp,0._dp) * TraceY(4) + TraceY(2)     &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul2(Ynu,sumnu1,OnlyDiagonal)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(4,1) = (3._dp,0._dp) * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = cTrace(aYeYeaYeYe)
   TraceY2(2) = cTrace(aYdYdaYdYd)
   TraceY2(3) = cTrace(aYuYuaYuYu)
   TraceY2(4) = cTrace(aYdYdaYuYu)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
  End If

  Call CouplingsToG2(Dgauge,DYe,DYnu,DYd,DYu,f)

  Iname = Iname - 1

 End Subroutine rge75


 Subroutine rge79(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, DYT
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21
  Real(dp) :: lam12, lam22, b_1a(3), Q, b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge79'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToCouplings4(gy, gauge, Ye, YT, Yd, Yu, lam1, lam2)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYTYT = Matmul2(aYT,YT,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul2(YT,sumT1,OnlyDiagonal)  &
        & + Matmul2(Transpose(aYeYe),YT,OnlyDiagonal)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
  End If

  Call CouplingsToG4(Dgauge, DYe, DYT, DYd, DYu, Dlam1, Dlam2, f)

  Iname = Iname - 1

 End Subroutine rge79


 Subroutine rge93(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 
 ! 30.08.2002: taking rge75 as basis and adding the RGEs for the
 !             left-neutrino mass matrix
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(4), Dgauge(3), q
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), TraceY2(4), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)             &
    & ,Ynu(3,3), aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3), DYnu(3,3)
  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3), sumM1(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge93'

  OnlyDiagonal = .Not.GenerationMixing

  Call GToCouplings3(gy,gauge,Ye,Ynu,Yd,Yu,Mnu)
  
  gauge2 = gauge**2
  q = t
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYnuYnu = Matmul2(aYnu,Ynu,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = (3._dp,0._dp) * TraceY(4) + TraceY(2)     &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul2(Ynu,sumnu1,OnlyDiagonal)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(4,1) = (3._dp,0._dp) * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  sumM1 = aYeYe + aYnuYnu
  diagonal(5,1) = 2._dp * TraceY(2) + 6._dp * TraceY(4)   &
              & - 2._dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(sumM1), Mnu) + Matmul(Mnu, sumM1)  &
          & + diagonal(5,1) * Mnu
  

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = cTrace(aYeYeaYeYe)
   TraceY2(2) = cTrace(aYdYdaYdYd)
   TraceY2(3) = cTrace(aYuYuaYuYu)
   TraceY2(4) = cTrace(aYdYdaYuYu)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2                                &
        & * ( b_1 + oo16pi2 * ( Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If

  Call CouplingsToG3(Dgauge,DYe,DYnu,DYd,DYu,DMnu,f)

  Iname = Iname - 1

 End Subroutine rge93


 Subroutine rge118(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1
  Real(dp) :: gauge(3), gauge2(3), TraceY(6), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, DYT 
  Complex(dp), Dimension(3,3) :: YZ, aYZ, aYZYZ, YZaYZ, betaYZ1, sumZ1, DYZ  &
    & , YdaYd
  Complex(dp), Dimension(3,3) :: YS, aYS, aYSYS, YSaYS, betaYS1, sumS1, DYS 
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21
  Real(dp) :: lam12, lam22, b_1a(3), Q, M15(3), betaM15(3), DM15(3), b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge118'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToCouplings5(gy, gauge, Ye, YT, Yd, Yu, YZ, YS, lam1, lam2, M15)

  gauge2 = gauge**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)
  Call Adjungate(YZ,aYZ)
  Call Adjungate(YS,aYS)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYTYT = Matmul2(aYT,YT,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)
  aYZYZ = Matmul2(aYZ,YZ,OnlyDiagonal)
  aYSYS = Matmul2(aYS,YS,OnlyDiagonal)

  YdaYd = Matmul2(Yd,aYd,OnlyDiagonal)
  YZaYZ = Matmul2(YZ,aYZ,OnlyDiagonal)
  YSaYS = Matmul2(YS,aYS,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )
  TraceY(5) = Real( cTrace(aYZYZ),dp )
  TraceY(6) = Real( cTrace(aYSYS),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT + aYZYZ)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT + 3._dp * aYZYZ
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul2(YT,sumT1,OnlyDiagonal)                 &
        & + Matmul2(Transpose(aYeYe+3._dp * aYZYZ),YT,OnlyDiagonal)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)  &
        & + 2._dp * Matmul2(YZaYZ + 2._dp * YSaYS, Yd,OnlyDiagonal)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  diagonal(5,1) = TraceY(5) + c1_1(2,1) * gauge2(1) &
            &   + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumZ1 = aYeYe + 3._dp * aYTYT + 5._dp * aYZYZ
  Do i1=1,3
   sumZ1(i1,i1) = sumZ1(i1,i1) + diagonal(5,1)
  End Do

  betaYZ1 = Matmul2(YZ,sumZ1,OnlyDiagonal)                 &
        & + 2._dp * Matmul2(YdaYd+2._dp * aYSYS,YZ,OnlyDiagonal)

  diagonal(6,1) = TraceY(6) - 0.8_dp * gauge2(1) - 12._dp * gauge2(3)
  sumS1 = 2._dp * Transpose(YdaYd + YZaYZ) + 8._dp * aYSYS
  Do i1=1,3
   sumS1(i1,i1) = sumS1(i1,i1) + diagonal(6,1)
  End Do

  betaYS1 = Matmul2(YS,sumS1,OnlyDiagonal)                 &
        & + 2._dp * Matmul2(YdaYd + YZaYZ,YS,OnlyDiagonal)


  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  !---------------------------------------------
  ! beta function for the 15-plet, MT, M_Z, M_S
  !---------------------------------------------
   betaM15(1) = M15(1) * ( TraceY(2) + lam12 + lam22 &
            &        - 2.4_dp * gauge2(1) - 8._dp * gauge2(2)  )
   betaM15(2) = M15(2) * ( TraceY(5) - gauge2(1)/15._dp - 3._dp * gauge2(2) &
                     & - 16._dp * gauge2(3) / 3._dp  )
   betaM15(3) = M15(3) * ( TraceY(6)  &
                     & - (3.2_dp * gauge2(1) + 40._dp * gauge2(3))/3._dp  )

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2a,TraceY(1:4))))
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
  !----------
  ! 15-plet
  !----------
   DM15 = oo16pi2 * betaM15
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
  !----------
  ! 15-plet
  !----------
   DM15 = oo16pi2 * betaM15
  End If

  Call CouplingsToG5(Dgauge, DYe, DYT, DYd, DYu, DYZ, DYS, Dlam1, Dlam2 &
       & , DM15, f)

  Iname = Iname - 1

 End Subroutine rge118


 Subroutine rge213(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), sumI, TraceY(3), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)  &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)        &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)            &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)          &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                 &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(5,2)       &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(3), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3), betaAu2(3,3)
  Real(dp) :: TraceA(3)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2, g2Mi(3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge213'

  OnlyDiagonal = .Not.GenerationMixing

  q = t

  Call GToParameters(gy, gauge, Ye, Yd, Yu                            &
                  & , Mhlf, Ae, Ad, Au, Me, Ml, Md, Mq, Mu, Mh, mue, B)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYdYd = Matmul(aYd,Yd)
  aYeYe = Matmul(aYe,Ye)
  aYuYu = Matmul(aYu,Yu)

  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYeYe(i1,i1) = Real(aYeYe(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
  End Do
  
  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * TraceY(2) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul(Ye,sume1)

  diagonal(2,1) = 3._dp * TraceY(2) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul(Yd,sumd1)

  diagonal(3,1) = 3._dp * TraceY(3)              &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do

  betaYu1 = Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul(aYdYd,aYdYd)
   aYeYeaYeYe = Matmul(aYeYe,aYeYe)
   aYuYuaYuYu = Matmul(aYuYu,aYuYu)
   aYuYuaYdYd = Matmul(aYuYu,aYdYd)
!   aYdYdaYuYu = Matmul(aYdYd,aYuYu)
   !------------------------------------------------
   ! this are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdYdaYdYd(i1,i1) = Real(aYdYdaYdYd(i1,i1), dp)
    aYeYeaYeYe(i1,i1) = Real(aYeYeaYeYe(i1,i1), dp)
    aYuYuaYuYu(i1,i1) = Real(aYuYuaYuYu(i1,i1), dp)
   End Do

   Call Adjungate(aYuYuaYdYd, aYdYdaYuYu)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(2)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(2) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul(Ye,sume2)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(3)
   hd(2) = 9._dp * TraceY(2) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul(Yd,sumd2)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(3)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(3) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(2) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul(Yu,sumu2)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ad,aAd)
  Call Adjungate(Ae,aAe)
  Call Adjungate(Au,aAu)

  aAdAd = Matmul(aAd,Ad)
  aAeAe = Matmul(aAe,Ae)
  aAuAu = Matmul(aAu,Au)

  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aAdAd(i1,i1) = Real(aAdAd(i1,i1), dp)
   aAeAe(i1,i1) = Real(aAeAe(i1,i1), dp)
   aAuAu(i1,i1) = Real(aAuAu(i1,i1), dp)
  End Do

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aAdAd),dp )
  TraceA(3) = Real( cTrace(aAuAu),dp )

  aYdAd = Matmul(aYd,Ad)
  aYeAe = Matmul(aYe,Ae)
  aYuAu = Matmul(aYu,Au)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYdAd) 
  TraceaYA(3) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = Matmul(Ae,sume1)

  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + Matmul(Ye,sume1)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = Matmul(Ad,sumd1)
  
  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1) &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do
  betaAd1 = betaAd1 + Matmul(Yd,sumd1)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = Matmul(Au,sumu1)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3)              &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
  End Do
  betaAu1 = betaAu1 + Matmul(Yu,sumu1)

  If (TwoLoopRGE) Then
   aYdYdaYdAd = Matmul(aYdYd,aYdAd)
   aYdAdaYdYd = Matmul(aYdAd,aYdYd)
   aYeYeaYeAe = Matmul(aYeYe,aYeAe)
   aYeAeaYeYe = Matmul(aYeAe,aYeYe)
   aYuYuaYuAu = Matmul(aYuYu,aYuAu)
   aYuAuaYuYu = Matmul(aYuAu,aYuYu)
   aYuAuaYdYd = Matmul(aYuAu,aYdYd)
   aYuYuaYdAd = Matmul(aYuYu,aYdAd)
   aYdAdaYuYu = Matmul(aYdAd,aYuYu)
   aYdYdaYuAu = Matmul(aYdYd,aYuAu)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = Matmul(Ae,sume2)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(2)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(2)                &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1) - 2.4_dp * g2Mi(1) * TraceY(1)    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(2) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(2) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + Matmul(Ye,sume2)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = Matmul(Ad,sumd2)
    
   diagonal(2,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(2) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(3) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(2) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(3) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
   betaAd2 = betaAd2 + Matmul(Yd,sumd2)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(3) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = Matmul(Au,sumu2)
    
   diagonal(3,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(3)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(3)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(2) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(3) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(3) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
   betaAu2 = betaAu2 + Matmul(Yu,sumu2)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)
   Do i1=1,3
    S1 = S1 + Real(Me(i1,i1),dp) - Real(Ml(i1,i1),dp) &
       &    + Real(Md(i1,i1),dp) + Real(Mq(i1,i1),dp) &
       &    - 2._dp * Real(Mu(i1,i1),dp)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = Matmul(Yd,aYd)
   YeaYe = Matmul(Ye,aYe)
   YuaYu = Matmul(Yu,aYu)

   MdYdaYd = Matmul(Md,YdaYd)
   MeYeaYe = Matmul(Me,YeaYe)
   MlaYeYe = Matmul(Ml,aYeYe)
   MqaYdYd = Matmul(Mq,aYdYd)
   MqaYuYu = Matmul(Mq,aYuYu)
   MuYuaYu = Matmul(Mu,YuaYu)

   Call Adjungate(MdYdaYd, YdaYdMd) ! YdaYdMd = Matmul(YdaYd,Md,OnlyDiagonal)
   Call Adjungate(MeYeaYe, YeaYeMe) ! YeaYeMe = Matmul(YeaYe,Me,OnlyDiagonal)
   Call Adjungate(MlaYeYe, aYeYeMl) ! aYeYeMl = Matmul(aYeYe,Ml,OnlyDiagonal)
   Call Adjungate(MqaYdYd, aYdYdMq) ! aYdYdMq = Matmul(aYdYd,Mq,OnlyDiagonal)
   Call Adjungate(MqaYuYu, aYuYuMq) ! aYuYuMq = Matmul(aYuYu,Mq,OnlyDiagonal)
   Call Adjungate(MuYuaYu, YuaYuMu) ! YuaYuMu = Matmul(YuaYu,Mu,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AdaAd = Matmul(Ad,aAd)
   AeaAe = Matmul(Ae,aAe)
   AuaAu = Matmul(Au,aAu)

   !------------------------------------------------
   ! these are hermitian matrices, clean up to
   ! avoid numerical problems
   !------------------------------------------------
   Do i1=1,3
    aYdMdYd(i1,i1) = Real(aYdMdYd(i1,i1), dp)
    aYeMeYe(i1,i1) = Real(aYeMeYe(i1,i1), dp)
    aYuMuYu(i1,i1) = Real(aYuMuYu(i1,i1), dp)
    YdMqaYd(i1,i1) = Real(YdMqaYd(i1,i1), dp)
    YeMlaYe(i1,i1) = Real(YeMlaYe(i1,i1), dp)
    YuMqaYu(i1,i1) = Real(YuMqaYu(i1,i1), dp)
    AdaAd(i1,i1) = Real(AdaAd(i1,i1), dp)
    AeaAe(i1,i1) = Real(AeaAe(i1,i1), dp)
    AuaAu(i1,i1) = Real(AuaAu(i1,i1), dp)
   End Do

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(2,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(2,1)
   End Do

   diagonal(3,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(3,1)
   End Do

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(5,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    !------------------------------------------------
    ! these are hermitian matrices, clean up to
    ! avoid numerical problems
    !------------------------------------------------
    Do i1=1,3
     YdaYdYdaYd(i1,i1) = Real( YdaYdYdaYd(i1,i1), dp)
     YeaYeYeaYe(i1,i1) = Real( YeaYeYeaYe(i1,i1), dp)
     YuaYuYuaYu(i1,i1) = Real( YuaYuYuaYu(i1,i1), dp)
    End Do

    AdaYd = Matmul(Ad,aYd)
    AeaYe = Matmul(Ae,aYe)
    AuaYu = Matmul(Au,aYu)

    aAdYd = Matmul(aAd,Yd)
    aAeYe = Matmul(aAe,Ye)
    aAuYu = Matmul(aAu,Yu)

    Call Adjungate(AdaYd,YdaAd) ! YdaAd = Matmul(Yd,aAd)
    Call Adjungate(AeaYe,YeaAe) ! YeaAe = Matmul(Ye,aAe)
    Call Adjungate(AuaYu,YuaAu) ! YuaAu = Matmul(Yu,aAu)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = Matmul(Md,YdaYuYuaYd)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = Matmul(Mu,YuaYdYdaYu)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = Matmul(MeYeaYe,YeaYe)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = Matmul(aYeMeYe,aYeYe)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = Matmul(YeaYeMe,YeaYe)

    MlaYeYeaYeYe = Matmul(MlaYeYe,aYeYe)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = Matmul(YeMlaYe,YeaYe)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = Matmul(aYeYeMl,aYeYe)

    MdYdaYdYdaYd = Matmul(MdYdaYd,YdaYd)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = Matmul(aYdMdYd,aYdYd)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = Matmul(YdaYdMd,YdaYd)

    MqaYdYdaYdYd = Matmul(MqaYdYd,aYdYd)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = Matmul(YdMqaYd,YdaYd)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = Matmul(aYdYdMq,aYdYd)

    MqaYuYuaYuYu = Matmul(MqaYuYu,aYuYu)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = Matmul(YuMqaYu,YuaYu)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = Matmul(aYuYuMq,aYuYu)

    MuYuaYuYuaYu = Matmul(MuYuaYu,YuaYu)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = Matmul(aYuMuYu,aYuYu)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = Matmul(YuaYuMu,YuaYu)

    AdaAdYdaYd = Matmul(AdaAd,YdaYd)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = Matmul(AdaYd,YdaAd)
    YdaAdAdaYd = Matmul(YdaAd,AdaYd)

    aAdAdaYdYd = Matmul(aAdAd,aYdYd)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = Matmul(aAdYd,aYdAd)
    aYdAdaAdYd = Matmul(aYdAd,aAdYd)

    AeaAeYeaYe = Matmul(AeaAe,YeaYe)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = Matmul(AeaYe,YeaAe)
    YeaAeAeaYe = Matmul(YeaAe,AeaYe)

    aAeAeaYeYe = Matmul(aAeAe,aYeYe)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = Matmul(aAeYe,aYeAe)
    aYeAeaAeYe = Matmul(aYeAe,aAeYe)

    AuaAuYuaYu = Matmul(AuaAu,YuaYu)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = Matmul(AuaYu,YuaAu)
    YuaAuAuaYu = Matmul(YuaAu,AuaYu)

    aAuAuaYuYu = Matmul(aAuAu,aYuYu)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = Matmul(aAuYu,aYuAu)
    aYuAuaAuYu = Matmul(aYuAu,aAuYu)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(3) - MH(1) * TraceY(2) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - Real(YuMqaYu(i1,i1),dp) + 4._dp * Real(aYuMuYu(i1,i1),dp)   &
        &    - Real(YdMqaYd(i1,i1),dp) - 2._dp * Real(aYdMdYd(i1,i1),dp)   &
        &    + Real(YeMlaYe(i1,i1),dp) - 2._dp * Real(aYeMeYe(i1,i1),dp)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(2) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
              & + 3._dp * ( Real(cTrace(MqaYdYd),dp)+Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(2) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(2) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                     &         +Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(3)
    Tr3aYuAu = 3._dp * TraceaYA(3)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(2,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(2,2)
    End Do

    diagonal(3,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(3) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(3) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(5,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * TraceY(1) + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(2) + Real( cTrace(YdMqaYd),dp ) &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(2)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(1) = mH(2) * TraceY(3) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(3)

   betamH21 = 6._dp * TraceMH2(1) - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1)  &
          & + 0.6_dp * S1

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
           & + Real(cTrace(YdMqaYuYuaYd),dp) + Real(cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(2)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(2),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(2)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(2),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
      & + 6._dp*(Real(cTrace(MqaYuYuaYuYu),dp)+Real(cTrace(aYuMuYuaYuYu),dp)  &
      &        + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp))   &
      & + Real(cTrace(YuMqaYdYdaYu),dp)+Real(cTrace(YuaYdMdYdaYu),dp )        &
      & + Real(cTrace(YuaYdYdMqaYu),dp)+Real(cTrace(YdaYuMuYuaYd),dp )        &
      & + Real(cTrace(YuaAdAdaYu),dp) + Real( cTrace(AuaYdYdaAu),dp )         &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                                  &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(3)                                  &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(2)+TraceY(3)) + TraceY(1)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(2)+TraceaYA(3)) + 2._dp * TraceaYA(1)    &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(3)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(2)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * (Matmul(AuaYu,YuaYu)     &
              &                 + Matmul(AdaYd,YdaYd) )  &
              &       + Matmul(AeaYe,YeaYe)              &
              &       + Matmul(aYuAu,aYdYd)              &
              &       + Matmul(aYdAd,aYuYu) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(2)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(3)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(2)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1)                          &
      &   - 30._dp * gauge2(2)**2 * Mhlf(2)                                   &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1 + oo16pi2 * (Matmul(b_2,gauge2) - Matmul(a_2,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) )     &
          &      + a_2(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1 * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  End If

  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do

  !--------------------------------------------
  ! This helps avoiding numerical instabilities
  !--------------------------------------------
!   Call Chop(DYe)
!   Call Chop(DYd)
!   Call Chop(DYu)
!   Call Chop(DMhlf)
!   Call Chop(DAe)
!   Call Chop(DAd)
!   Call Chop(DAu)
!   Call Chop(DMe)
!   Call Chop(DMl)
!   Call Chop(DMd)
!   Call Chop(DMu)
!   Call Chop(DMq)
!   Call Chop(Dmue)
!   Call Chop(DB)


  Call ParametersToG(Dgauge, DYe, DYd, DYu, DMhlf, DAe, DAd, DAu &
                   &, DMe, DMl, DMd, DMq, DMu, DMh, Dmue, DB, f)

  Iname = Iname - 1

 End Subroutine rge213


 Subroutine rge267(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.03.2001: implementing right-handed neutrinos  at 1-loop
 !             up to now: Y_nu, A_nu,
 !             the parameters m_H1, M_H2, B and mu still need to be changed
 ! 07.10.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), sumI, TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4),Ynu(3,3) &
    & , aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3) ,DYnu(3,3)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(4), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3), Anu(3,3), aAnu(3,3), DAnu(3,3) ,aAnuAnu(3,3)       &
     &  , betaAnu1(3,3), aYnuAnu(3,3)!, betaAnu2(3,3)
  Real(dp) :: TraceA(4)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)
  Complex(dp) :: Mr(3,3), YnuaYnu(3,3), MlaYnuYnu(3,3), MrYnuaYnu(3,3)        &
     & , aYnuYnuMl(3,3), YnuaYnuMr(3,3), YnuMlaYnu(3,3), aYnuMrYnu(3,3)       &
     & , AnuaAnu(3,3), DMr(3,3), betaMr1(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2, g2Mi(3)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge267'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters2(gy, gauge, Ye, Ynu, Yd, Yu                            &
                & , Mhlf, Ae, Anu, Ad, Au, Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYnuYnu = Matmul2(aYnu,Ynu,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = 3._dp * TraceY(4) + TraceY(2)           &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul2(Ynu,sumnu1,OnlyDiagonal)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(4,1) = (3._dp,0._dp) * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(Anu,aAnu)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)

  aAdAd = MatMul2(aAd,Ad,OnlyDiagonal)
  aAnuAnu = MatMul2(aAnu,Anu,OnlyDiagonal)
  aAeAe = MatMul2(aAe,Ae,OnlyDiagonal)
  aAuAu = MatMul2(aAu,Au,OnlyDiagonal)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aAnuAnu),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )

  aYdAd = MatMul2(aYd,Ad,OnlyDiagonal)
  aYnuAnu = MatMul2(aYnu,Anu,OnlyDiagonal)
  aYeAe = MatMul2(aYe,Ae,OnlyDiagonal)
  aYuAu = MatMul2(aYu,Au,OnlyDiagonal)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYnuAnu) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = MatMul2(Ae,sume1,OnlyDiagonal)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe + 2._dp * aYnuAnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + MatMul2(Ye,sume1,OnlyDiagonal)

  !--------------
  ! A_nu
  !--------------
  sumnu1 = sumnu1 + 2._dp * aYnuYnu
  betaAnu1 = MatMul2(Anu,sumnu1,OnlyDiagonal)

  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2)     &
      &         + 0.6_dp * g2Mi(1) + 3._dp * g2Mi(2) )
  sumnu1 = 4._dp * aYnuAnu + 2._dp * aYeAe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do
  betaAnu1 = betaAnu1 + MatMul2(Ynu,sumnu1,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = MatMul2(Ad,sumd1,OnlyDiagonal)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)           &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + MatMul2(Yd,sumd1,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = MatMul2(Au,sumu1,OnlyDiagonal)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2)           &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + MatMul2(Yu,sumu1,OnlyDiagonal)

  If (TwoLoopRGE) Then
   aYdYdaYdAd = MatMul2(aYdYd,aYdAd,OnlyDiagonal)
   aYdAdaYdYd = MatMul2(aYdAd,aYdYd,OnlyDiagonal)
   aYeYeaYeAe = MatMul2(aYeYe,aYeAe,OnlyDiagonal)
   aYeAeaYeYe = MatMul2(aYeAe,aYeYe,OnlyDiagonal)
   aYuYuaYuAu = MatMul2(aYuYu,aYuAu,OnlyDiagonal)
   aYuAuaYuYu = MatMul2(aYuAu,aYuYu,OnlyDiagonal)
   aYuAuaYdYd = MatMul2(aYuAu,aYdYd,OnlyDiagonal)
   aYuYuaYdAd = MatMul2(aYuYu,aYdAd,OnlyDiagonal)
   aYdAdaYuYu = MatMul2(aYdAd,aYuYu,OnlyDiagonal)
   aYdYdaYuAu = MatMul2(aYdYd,aYuAu,OnlyDiagonal)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = MatMul2(Ae,sume2,OnlyDiagonal)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)         &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1)                            &
     &  - 2.4_dp * g2Mi(1) * TraceY(1)                    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + MatMul2(Ye,sume2,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = MatMul2(Ad,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + MatMul2(Yd,sumd2,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = MatMul2(Au,sumu2,OnlyDiagonal)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + MatMul2(Yu,sumu2,OnlyDiagonal)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = MatMul2(Yd,aYd,OnlyDiagonal)
   YnuaYnu = MatMul2(Ynu,aYnu,OnlyDiagonal)
   YeaYe = MatMul2(Ye,aYe,OnlyDiagonal)
   YuaYu = MatMul2(Yu,aYu,OnlyDiagonal)

   MeYeaYe = MatMul2(Me,YeaYe,OnlyDiagonal)
   MlaYeYe = MatMul2(Ml,aYeYe,OnlyDiagonal)
   MlaYnuYnu = MatMul2(Ml,aYnuYnu,OnlyDiagonal)
   MrYnuaYnu = MatMul2(Mr,YnuaYnu,OnlyDiagonal)

   MdYdaYd = MatMul2(Md,YdaYd,OnlyDiagonal)
   MqaYdYd = MatMul2(Mq,aYdYd,OnlyDiagonal)
   MqaYuYu = MatMul2(Mq,aYuYu,OnlyDiagonal)
   MuYuaYu = MatMul2(Mu,YuaYu,OnlyDiagonal)

   YeaYeMe = MatMul2(YeaYe,Me,OnlyDiagonal)
   aYeYeMl = MatMul2(aYeYe,Ml,OnlyDiagonal)
   aYnuYnuMl = MatMul2(aYnuYnu,Ml,OnlyDiagonal)
   YnuaYnuMr = MatMul2(YnuaYnu,Mr,OnlyDiagonal)

   YdaYdMd = MatMul2(YdaYd,Md,OnlyDiagonal)
   aYdYdMq = MatMul2(aYdYd,Mq,OnlyDiagonal)
   aYuYuMq = MatMul2(aYuYu,Mq,OnlyDiagonal)
   YuaYuMu = MatMul2(YuaYu,Mu,OnlyDiagonal)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YnuMlaYnu = MatMul3(Ynu,Ml,aYnu,OnlyDiagonal)
   aYnuMrYnu = MatMul3(aYnu,Mr,Ynu,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AeaAe = MatMul2(Ae,aAe,OnlyDiagonal)
   AnuaAnu = MatMul2(Anu,aAnu,OnlyDiagonal)
   AdaAd = MatMul2(Ad,aAd,OnlyDiagonal)
   AuaAu = MatMul2(Au,aAu,OnlyDiagonal)

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + MlaYnuYnu + aYnuYnuMl                                            &
         & + 2._dp * ( mH(2) * aYnuYnu + aYnuMrYnu + aAnuAnu )
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   betaMr1 = 2._dp * (MrYnuaYnu + YnuaYnuMr)                        &
         & + 4._dp * ( mH(2) * YnuaYnu + YnuMlaYnu + AnuaAnu )

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AdaYd = MatMul2(Ad,aYd,OnlyDiagonal)
    AeaYe = MatMul2(Ae,aYe,OnlyDiagonal)
    AuaYu = MatMul2(Au,aYu,OnlyDiagonal)

    aAdYd = MatMul2(aAd,Yd,OnlyDiagonal)
    aAeYe = MatMul2(aAe,Ye,OnlyDiagonal)
    aAuYu = MatMul2(aAu,Yu,OnlyDiagonal)

    YdaAd = MatMul2(Yd,aAd,OnlyDiagonal)
    YeaAe = MatMul2(Ye,aAe,OnlyDiagonal)
    YuaAu = MatMul2(Yu,aAu,OnlyDiagonal)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = MatMul2(Md,YdaYuYuaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = MatMul2(Mu,YuaYdYdaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = MatMul2(MeYeaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = MatMul2(aYeMeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = MatMul2(YeaYeMe,YeaYe,OnlyDiagonal)

    MlaYeYeaYeYe = MatMul2(MlaYeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = MatMul2(YeMlaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = MatMul2(aYeYeMl,aYeYe,OnlyDiagonal)

    MdYdaYdYdaYd = MatMul2(MdYdaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = MatMul2(aYdMdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = MatMul2(YdaYdMd,YdaYd,OnlyDiagonal)

    MqaYdYdaYdYd = MatMul2(MqaYdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = MatMul2(YdMqaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = MatMul2(aYdYdMq,aYdYd,OnlyDiagonal)

    MqaYuYuaYuYu = MatMul2(MqaYuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = MatMul2(YuMqaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = MatMul2(aYuYuMq,aYuYu,OnlyDiagonal)

    MuYuaYuYuaYu = MatMul2(MuYuaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = MatMul2(aYuMuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = MatMul2(YuaYuMu,YuaYu,OnlyDiagonal)

    AdaAdYdaYd = MatMul2(AdaAd,YdaYd,OnlyDiagonal)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = MatMul2(AdaYd,YdaAd,OnlyDiagonal)
    YdaAdAdaYd = MatMul2(YdaAd,AdaYd,OnlyDiagonal)

    aAdAdaYdYd = MatMul2(aAdAd,aYdYd,OnlyDiagonal)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = MatMul2(aAdYd,aYdAd,OnlyDiagonal)
    aYdAdaAdYd = MatMul2(aYdAd,aAdYd,OnlyDiagonal)

    AeaAeYeaYe = MatMul2(AeaAe,YeaYe,OnlyDiagonal)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = MatMul2(AeaYe,YeaAe,OnlyDiagonal)
    YeaAeAeaYe = MatMul2(YeaAe,AeaYe,OnlyDiagonal)

    aAeAeaYeYe = MatMul2(aAeAe,aYeYe,OnlyDiagonal)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = MatMul2(aAeYe,aYeAe,OnlyDiagonal)
    aYeAeaAeYe = MatMul2(aYeAe,aAeYe,OnlyDiagonal)

    AuaAuYuaYu = MatMul2(AuaAu,YuaYu,OnlyDiagonal)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = MatMul2(AuaYu,YuaAu,OnlyDiagonal)
    YuaAuAuaYu = MatMul2(YuaAu,AuaYu,OnlyDiagonal)

    aAuAuaYuYu = MatMul2(aAuAu,aYuYu,OnlyDiagonal)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = MatMul2(aAuYu,aYuAu,OnlyDiagonal)
    aYuAuaAuYu = MatMul2(aYuAu,aAuYu,OnlyDiagonal)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * TraceY(1) + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp ) &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(2) = mH(2) * TraceY(2) + Real( cTrace(YnuMlaYnu),dp )  &
             & + Real( cTrace(aYnuMrYnu),dp ) + TraceA(2)
   traceMH2(1) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(2) + 6._dp * TraceMH2(1)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)) + TraceY(1) + TraceY(2)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4))                          &
           & + 2._dp * (TraceaYA(1)+ TraceaYA(2))                         &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul2(AuaYu,YuaYu,OnlyDiagonal)     &
              &                 + Matmul2(AdaYd,YdaYd,OnlyDiagonal) )   &
              &       + MatMul2(AeaYe,YeaYe,OnlyDiagonal)               &
              &       + Matmul2(aYuAu,aYdYd,OnlyDiagonal)               &
              &       + MatMul2(aYdAd,aYuYu,OnlyDiagonal) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1) - 30._dp * gauge2(2)**2 * Mhlf(2)    &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1 + oo16pi2 * (Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,4
     sumI = sumI + a_2a(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMr = oo16pi2 * betaMr1 ! + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1 * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMr = oo16pi2 * betaMr1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMr(i1,i1) = Real(DMr(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do

  Call ParametersToG2(Dgauge, DYe, DYnu, DYd, DYu, DMhlf, DAe, DAnu, DAd, DAu &
                   &, DMe, DMl, DMr, DMd, DMq, DMu, DMh, Dmue, DB, f)

  Iname = Iname - 1

 End Subroutine rge267


 Subroutine rge277(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), sumI, TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(4), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3)
  Real(dp) :: TraceA(4)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2

  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, AT, aAT     &
      & , aATAT, betaAT1, aYTAT, MlaYTYT, aYTMlYT, aYTYTMl , DYT, DAT
  Complex(dp) :: g2Mi(3)
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21, lam1Alam1 &
      & , lam2Alam2, Alam1, Alam2, betaAlam11, betaAlam21, DAlam1, DAlam2
  Real(dp) :: MT(2), lam12, lam22, betaMT2(2), Alam12, Alam22, b_1a(3)     &
      & , DMT2(2), b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge277'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters4(gy, gauge, Ye, YT, Yd, Yu, lam1, lam2, Mhlf, Ae, AT, Ad &
                & , Au, Alam1, Alam2, Me, Ml, Md, Mq, Mu, Mh, mT, mue, B, Mnu)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYTYT = Matmul2(aYT,YT,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul2(YT,sumT1,OnlyDiagonal)  &
        & + Matmul2(Transpose(aYeYe),YT,OnlyDiagonal)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(AT,aAT)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)

  aAdAd = MatMul2(aAd,Ad,OnlyDiagonal)
  aATAT = MatMul2(aAT,AT,OnlyDiagonal)
  aAeAe = MatMul2(aAe,Ae,OnlyDiagonal)
  aAuAu = MatMul2(aAu,Au,OnlyDiagonal)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aATAT),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )

  aYdAd = MatMul2(aYd,Ad,OnlyDiagonal)
  aYTAT = MatMul2(aYT,AT,OnlyDiagonal)
  aYeAe = MatMul2(aYe,Ae,OnlyDiagonal)
  aYuAu = MatMul2(aYu,Au,OnlyDiagonal)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYTAT) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 
  lam1Alam1 = Conjg(lam1) * Alam1
  lam2Alam2 = Conjg(lam2) * Alam2
  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = MatMul2(Ae,sume1,OnlyDiagonal)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         + 3._dp * lam1Alam1                  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    ) 
  sume1 = 4._dp * aYeAe + 6._dp * aYTAT
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + MatMul2(Ye,sume1,OnlyDiagonal)

  !--------------
  ! A_T
  !--------------
   If (decoupling_heavy_states) then
    betaAT1 = 0._dp
   Else  
    diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
    sumT1 = aYeYe + 9._dp * aYTYT
    betaAT1 = MatMul2(AT,sumT1,OnlyDiagonal)
    betaAT1 = Transpose(betaAT1)
    Do i1=1,3
     sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
    End Do

    betaAT1 = betaAT1 + MatMul2(sumT1,AT,OnlyDiagonal)
  
    diagonal(2,1) = 2._dp * ( TraceaYA(2) + lam1Alam1    &
      &         + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) )
    betaAT1 = betaAT1 + diagonal(2,1) * YT
   End If
  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = MatMul2(Ad,sumd1,OnlyDiagonal)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1) &
                &         + 3._dp * lam1Alam1                 &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + MatMul2(Yd,sumd1,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = MatMul2(Au,sumu1,OnlyDiagonal)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4)                       &
                &         + 3._dp * lam2Alam2 - c1_1(3,1) * g2Mi(1)   &
                &         - c1_1(3,2) * g2Mi(2) - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + MatMul2(Yu,sumu1,OnlyDiagonal)

  !--------------
  ! A_1
  !--------------
  betaAlam11 = Alam1 * (21._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
          &           + 6._dp * TraceY(3)                              &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )       &
          & + 2._dp * lam1 * ( TraceaYA(2) + 2._dp * TraceaYA(1)       &
          &                  + 6._dp * TraceaYA(3)                     &
          &                  + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) ) 

  !--------------
  ! A_2
  !--------------
  betaAlam21 = Alam2 * (21._dp * lam22 + 6._dp * TraceY(4)          &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )    &
          & + 2._dp * lam2 * ( 6._dp * TraceaYA(3)                  &
          &                  + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdAd = MatMul2(aYdYd,aYdAd,OnlyDiagonal)
   aYdAdaYdYd = MatMul2(aYdAd,aYdYd,OnlyDiagonal)
   aYeYeaYeAe = MatMul2(aYeYe,aYeAe,OnlyDiagonal)
   aYeAeaYeYe = MatMul2(aYeAe,aYeYe,OnlyDiagonal)
   aYuYuaYuAu = MatMul2(aYuYu,aYuAu,OnlyDiagonal)
   aYuAuaYuYu = MatMul2(aYuAu,aYuYu,OnlyDiagonal)
   aYuAuaYdYd = MatMul2(aYuAu,aYdYd,OnlyDiagonal)
   aYuYuaYdAd = MatMul2(aYuYu,aYdAd,OnlyDiagonal)
   aYdAdaYuYu = MatMul2(aYdAd,aYuYu,OnlyDiagonal)
   aYdYdaYuAu = MatMul2(aYdYd,aYuAu,OnlyDiagonal)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = MatMul2(Ae,sume2,OnlyDiagonal)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)             &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1) - 2.4_dp * g2Mi(1) * TraceY(1) &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + MatMul2(Ye,sume2,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = MatMul2(Ad,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + MatMul2(Yd,sumd2,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = MatMul2(Au,sumu2,OnlyDiagonal)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + MatMul2(Yu,sumu2,OnlyDiagonal)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)  + 3._dp * (mT(1) - mT(2) )
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = MatMul2(Yd,aYd,OnlyDiagonal)
   YeaYe = MatMul2(Ye,aYe,OnlyDiagonal)
   YuaYu = MatMul2(Yu,aYu,OnlyDiagonal)

   MeYeaYe = MatMul2(Me,YeaYe,OnlyDiagonal)
   MlaYeYe = MatMul2(Ml,aYeYe,OnlyDiagonal)
   MlaYTYT = MatMul2(Ml,aYTYT,OnlyDiagonal)

   MdYdaYd = MatMul2(Md,YdaYd,OnlyDiagonal)
   MqaYdYd = MatMul2(Mq,aYdYd,OnlyDiagonal)
   MqaYuYu = MatMul2(Mq,aYuYu,OnlyDiagonal)
   MuYuaYu = MatMul2(Mu,YuaYu,OnlyDiagonal)

   YeaYeMe = MatMul2(YeaYe,Me,OnlyDiagonal)
   aYeYeMl = MatMul2(aYeYe,Ml,OnlyDiagonal)
   aYTYTMl = MatMul2(aYTYT,Ml,OnlyDiagonal)

   YdaYdMd = MatMul2(YdaYd,Md,OnlyDiagonal)
   aYdYdMq = MatMul2(aYdYd,Mq,OnlyDiagonal)
   aYuYuMq = MatMul2(aYuYu,Mq,OnlyDiagonal)
   YuaYuMu = MatMul2(YuaYu,Mu,OnlyDiagonal)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   aYTMlYT = MatMul3(aYT,Transpose(Ml),YT,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AeaAe = MatMul2(Ae,aAe,OnlyDiagonal)
   AdaAd = MatMul2(Ad,aAd,OnlyDiagonal)
   AuaAu = MatMul2(Au,aAu,OnlyDiagonal)
   Alam12 = Abs(Alam1)**2
   Alam22 = Abs(Alam2)**2

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + 3._dp * ( MlaYTYT + aYTYTMl )                                    &
         & + 6._dp * ( aYTMlYT + aATAT + MT(1) * aYTYT)
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AdaYd = MatMul2(Ad,aYd,OnlyDiagonal)
    AeaYe = MatMul2(Ae,aYe,OnlyDiagonal)
    AuaYu = MatMul2(Au,aYu,OnlyDiagonal)

    aAdYd = MatMul2(aAd,Yd,OnlyDiagonal)
    aAeYe = MatMul2(aAe,Ye,OnlyDiagonal)
    aAuYu = MatMul2(aAu,Yu,OnlyDiagonal)

    YdaAd = MatMul2(Yd,aAd,OnlyDiagonal)
    YeaAe = MatMul2(Ye,aAe,OnlyDiagonal)
    YuaAu = MatMul2(Yu,aAu,OnlyDiagonal)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = MatMul2(Md,YdaYuYuaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = MatMul2(Mu,YuaYdYdaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = MatMul2(MeYeaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = MatMul2(aYeMeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = MatMul2(YeaYeMe,YeaYe,OnlyDiagonal)

    MlaYeYeaYeYe = MatMul2(MlaYeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = MatMul2(YeMlaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = MatMul2(aYeYeMl,aYeYe,OnlyDiagonal)

    MdYdaYdYdaYd = MatMul2(MdYdaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = MatMul2(aYdMdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = MatMul2(YdaYdMd,YdaYd,OnlyDiagonal)

    MqaYdYdaYdYd = MatMul2(MqaYdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = MatMul2(YdMqaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = MatMul2(aYdYdMq,aYdYd,OnlyDiagonal)

    MqaYuYuaYuYu = MatMul2(MqaYuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = MatMul2(YuMqaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = MatMul2(aYuYuMq,aYuYu,OnlyDiagonal)

    MuYuaYuYuaYu = MatMul2(MuYuaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = MatMul2(aYuMuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = MatMul2(YuaYuMu,YuaYu,OnlyDiagonal)

    AdaAdYdaYd = MatMul2(AdaAd,YdaYd,OnlyDiagonal)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = MatMul2(AdaYd,YdaAd,OnlyDiagonal)
    YdaAdAdaYd = MatMul2(YdaAd,AdaYd,OnlyDiagonal)

    aAdAdaYdYd = MatMul2(aAdAd,aYdYd,OnlyDiagonal)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = MatMul2(aAdYd,aYdAd,OnlyDiagonal)
    aYdAdaAdYd = MatMul2(aYdAd,aAdYd,OnlyDiagonal)

    AeaAeYeaYe = MatMul2(AeaAe,YeaYe,OnlyDiagonal)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = MatMul2(AeaYe,YeaAe,OnlyDiagonal)
    YeaAeAeaYe = MatMul2(YeaAe,AeaYe,OnlyDiagonal)

    aAeAeaYeYe = MatMul2(aAeAe,aYeYe,OnlyDiagonal)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = MatMul2(aAeYe,aYeAe,OnlyDiagonal)
    aYeAeaAeYe = MatMul2(aYeAe,aAeYe,OnlyDiagonal)

    AuaAuYuaYu = MatMul2(AuaAu,YuaYu,OnlyDiagonal)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = MatMul2(AuaYu,YuaAu,OnlyDiagonal)
    YuaAuAuaYu = MatMul2(YuaAu,AuaYu,OnlyDiagonal)

    aAuAuaYuYu = MatMul2(aAuAu,aYuYu,OnlyDiagonal)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = MatMul2(aAuYu,aYuAu,OnlyDiagonal)
    aYuAuaAuYu = MatMul2(aYuAu,aAuYu,OnlyDiagonal)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * (TraceY(1) + 6._dp * lam12)            &
             & + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp )   &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3) + Alam12  &
             & + lam12 * MT(1)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(1) = 3._dp * (2._dp*mH(2)+ MT(2)) * lam22 + 3._dp * Alam22
   traceMH2(2) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(1) + 6._dp * TraceMH2(2)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   If (decoupling_heavy_states) then
    betaMT2 = 0._dp 

   else
    betaMT2(1) = MT(1) * (lam12 + TraceY(2)) + 2._dp * mH(1) * lam12   &
            & + 2._dp * Real( cTrace(aYTMlYT),dp ) + TraceA(2) + Alam12 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

    betaMT2(2) = MT(2) * lam22 + 2._dp * mH(2) * lam22 + Alam22 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

    betaMT2 = 2._dp * betaMT2
   End If

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(1) = traceMH2(2)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)+lam12+lam22) + TraceY(1)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4) + lam1Alam1 + lam2Alam2) &
           & + 2._dp * TraceaYA(1) + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul2(AuaYu,YuaYu,OnlyDiagonal)     &
              &                 + Matmul2(AdaYd,YdaYd,OnlyDiagonal) )   &
              &       + MatMul2(AeaYe,YeaYe,OnlyDiagonal)               &
              &       + Matmul2(aYuAu,aYdYd,OnlyDiagonal)               &
              &       + MatMul2(aYdAd,aYuYu,OnlyDiagonal) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1) - 30._dp * gauge2(2)**2 * Mhlf(2)    &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  diagonal(5,1) = 6._dp * TraceY(4) - 2._dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(aYeYe), Mnu) + Matmul(Mnu, aYeYe)  &
          & + diagonal(5,1) * Mnu
  

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2a(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,4
     sumI = sumI + a_2a(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1a(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
   DAT = oo16pi2 * betaAT1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
   DmT2 = oo16pi2 * betaMT2
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1a * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
   DAT = oo16pi2 * betaAT1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
   DmT2 = oo16pi2 * betaMT2
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  !----------------------
  ! left neutrino mass
  !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do

  Call ParametersToG4(Dgauge, DYe, DYT, DYd, DYu, Dlam1, Dlam2, DMhlf, DAe, DAT &
          & , DAd, DAu, DAlam1, DAlam2, DMe, DMl, DMd, DMq, DMu, DMh, DMT2      &
          & , Dmue, DB, DMnu, f)

  Iname = Iname - 1

 End Subroutine rge277


 Subroutine rge285(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.03.2001: implementing right-handed neutrinos  at 1-loop
 !             up to now: Y_nu, A_nu,
 !             the parameters m_H1, M_H2, B and mu still need to be changed
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), sumI, TraceY(4), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4),Ynu(3,3) &
    & , aYnu(3,3), aYnuYnu(3,3), sumnu1(3,3), betaYnu1(3,3) ,DYnu(3,3)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(4), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3), Anu(3,3), aAnu(3,3), DAnu(3,3) ,aAnuAnu(3,3)       &
     &  , betaAnu1(3,3), aYnuAnu(3,3)!, betaAnu2(3,3)
  Real(dp) :: TraceA(4)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)
  Complex(dp) :: Mr(3,3), YnuaYnu(3,3), MlaYnuYnu(3,3), MrYnuaYnu(3,3)        &
     & , aYnuYnuMl(3,3), YnuaYnuMr(3,3), YnuMlaYnu(3,3), aYnuMrYnu(3,3)       &
     & , AnuaAnu(3,3), DMr(3,3), betaMr1(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2, g2Mi(3)

  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3), sumM1(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge285'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters3(gy, gauge, Ye, Ynu, Yd, Yu, Mhlf, Ae, Anu, Ad, Au      &
                & , Me, Ml, Mr, Md, Mq, Mu, Mh, mue, B, Mnu)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  Call Adjungate(Yd,aYd)
  Call Adjungate(Ynu,aYnu)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYnuYnu = Matmul2(aYnu,Ynu,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYnuYnu),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )

  diagonal(1,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe + aYnuYnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = 3._dp * TraceY(4) + TraceY(2)           &
            &   - 0.6_dp * gauge2(1) - 3._dp * gauge2(2)
  sumnu1 = 3._dp * aYnuYnu + aYeYe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do

  betaYnu1 = Matmul2(Ynu,sumnu1,OnlyDiagonal)

  diagonal(3,1) = (3._dp,0._dp) * TraceY(3) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(4,1) = (3._dp,0._dp) * TraceY(4)  + TraceY(2)             &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(Anu,aAnu)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)

  aAdAd = MatMul2(aAd,Ad,OnlyDiagonal)
  aAnuAnu = MatMul2(aAnu,Anu,OnlyDiagonal)
  aAeAe = MatMul2(aAe,Ae,OnlyDiagonal)
  aAuAu = MatMul2(aAu,Au,OnlyDiagonal)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aAnuAnu),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )

  aYdAd = MatMul2(aYd,Ad,OnlyDiagonal)
  aYnuAnu = MatMul2(aYnu,Anu,OnlyDiagonal)
  aYeAe = MatMul2(aYe,Ae,OnlyDiagonal)
  aYuAu = MatMul2(aYu,Au,OnlyDiagonal)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYnuAnu) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = MatMul2(Ae,sume1,OnlyDiagonal)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe + 2._dp * aYnuAnu
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + MatMul2(Ye,sume1,OnlyDiagonal)

  !--------------
  ! A_nu
  !--------------
  sumnu1 = sumnu1 + 2._dp * aYnuYnu
  betaAnu1 = MatMul2(Anu,sumnu1,OnlyDiagonal)

  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2)     &
      &         + 0.6_dp * g2Mi(1) + 3._dp * g2Mi(2) )
  sumnu1 = 4._dp * aYnuAnu + 2._dp * aYeAe
  Do i1=1,3
   sumnu1(i1,i1) = sumnu1(i1,i1) + diagonal(2,1)
  End Do
  betaAnu1 = betaAnu1 + MatMul2(Ynu,sumnu1,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = MatMul2(Ad,sumd1,OnlyDiagonal)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1) &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + MatMul2(Yd,sumd1,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = MatMul2(Au,sumu1,OnlyDiagonal)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4) + TraceaYA(2) &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + MatMul2(Yu,sumu1,OnlyDiagonal)

  If (TwoLoopRGE) Then
   aYdYdaYdAd = MatMul2(aYdYd,aYdAd,OnlyDiagonal)
   aYdAdaYdYd = MatMul2(aYdAd,aYdYd,OnlyDiagonal)
   aYeYeaYeAe = MatMul2(aYeYe,aYeAe,OnlyDiagonal)
   aYeAeaYeYe = MatMul2(aYeAe,aYeYe,OnlyDiagonal)
   aYuYuaYuAu = MatMul2(aYuYu,aYuAu,OnlyDiagonal)
   aYuAuaYuYu = MatMul2(aYuAu,aYuYu,OnlyDiagonal)
   aYuAuaYdYd = MatMul2(aYuAu,aYdYd,OnlyDiagonal)
   aYuYuaYdAd = MatMul2(aYuYu,aYdAd,OnlyDiagonal)
   aYdAdaYuYu = MatMul2(aYdAd,aYuYu,OnlyDiagonal)
   aYdYdaYuAu = MatMul2(aYdYd,aYuAu,OnlyDiagonal)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = MatMul2(Ae,sume2,OnlyDiagonal)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)                &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1) - 2.4_dp * g2Mi(1) * TraceY(1)    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + MatMul2(Ye,sume2,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = MatMul2(Ad,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + MatMul2(Yd,sumd2,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = MatMul2(Au,sumu2,OnlyDiagonal)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)           &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + MatMul2(Yu,sumu2,OnlyDiagonal)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1)
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do
   S1 = S1 * gauge2(1)

   YdaYd = MatMul2(Yd,aYd,OnlyDiagonal)
   YnuaYnu = MatMul2(Ynu,aYnu,OnlyDiagonal)
   YeaYe = MatMul2(Ye,aYe,OnlyDiagonal)
   YuaYu = MatMul2(Yu,aYu,OnlyDiagonal)

   MeYeaYe = MatMul2(Me,YeaYe,OnlyDiagonal)
   MlaYeYe = MatMul2(Ml,aYeYe,OnlyDiagonal)
   MlaYnuYnu = MatMul2(Ml,aYnuYnu,OnlyDiagonal)
   MrYnuaYnu = MatMul2(Mr,YnuaYnu,OnlyDiagonal)

   MdYdaYd = MatMul2(Md,YdaYd,OnlyDiagonal)
   MqaYdYd = MatMul2(Mq,aYdYd,OnlyDiagonal)
   MqaYuYu = MatMul2(Mq,aYuYu,OnlyDiagonal)
   MuYuaYu = MatMul2(Mu,YuaYu,OnlyDiagonal)

   YeaYeMe = MatMul2(YeaYe,Me,OnlyDiagonal)
   aYeYeMl = MatMul2(aYeYe,Ml,OnlyDiagonal)
   aYnuYnuMl = MatMul2(aYnuYnu,Ml,OnlyDiagonal)
   YnuaYnuMr = MatMul2(YnuaYnu,Mr,OnlyDiagonal)

   YdaYdMd = MatMul2(YdaYd,Md,OnlyDiagonal)
   aYdYdMq = MatMul2(aYdYd,Mq,OnlyDiagonal)
   aYuYuMq = MatMul2(aYuYu,Mq,OnlyDiagonal)
   YuaYuMu = MatMul2(YuaYu,Mu,OnlyDiagonal)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YnuMlaYnu = MatMul3(Ynu,Ml,aYnu,OnlyDiagonal)
   aYnuMrYnu = MatMul3(aYnu,Mr,Ynu,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AeaAe = MatMul2(Ae,aAe,OnlyDiagonal)
   AnuaAnu = MatMul2(Anu,aAnu,OnlyDiagonal)
   AdaAd = MatMul2(Ad,aAd,OnlyDiagonal)
   AuaAu = MatMul2(Au,aAu,OnlyDiagonal)

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + MlaYnuYnu + aYnuYnuMl                                            &
         & + 2._dp * ( mH(2) * aYnuYnu + aYnuMrYnu + aAnuAnu )
   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   betaMr1 = 2._dp * (MrYnuaYnu + YnuaYnuMr)                        &
         & + 4._dp * ( mH(2) * YnuaYnu + YnuMlaYnu + AnuaAnu )

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd)             &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AdaYd = MatMul2(Ad,aYd,OnlyDiagonal)
    AeaYe = MatMul2(Ae,aYe,OnlyDiagonal)
    AuaYu = MatMul2(Au,aYu,OnlyDiagonal)

    aAdYd = MatMul2(aAd,Yd,OnlyDiagonal)
    aAeYe = MatMul2(aAe,Ye,OnlyDiagonal)
    aAuYu = MatMul2(aAu,Yu,OnlyDiagonal)

    YdaAd = MatMul2(Yd,aAd,OnlyDiagonal)
    YeaAe = MatMul2(Ye,aAe,OnlyDiagonal)
    YuaAu = MatMul2(Yu,aAu,OnlyDiagonal)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = MatMul2(Md,YdaYuYuaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = MatMul2(Mu,YuaYdYdaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = MatMul2(MeYeaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = MatMul2(aYeMeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = MatMul2(YeaYeMe,YeaYe,OnlyDiagonal)

    MlaYeYeaYeYe = MatMul2(MlaYeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = MatMul2(YeMlaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = MatMul2(aYeYeMl,aYeYe,OnlyDiagonal)

    MdYdaYdYdaYd = MatMul2(MdYdaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = MatMul2(aYdMdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = MatMul2(YdaYdMd,YdaYd,OnlyDiagonal)

    MqaYdYdaYdYd = MatMul2(MqaYdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = MatMul2(YdMqaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = MatMul2(aYdYdMq,aYdYd,OnlyDiagonal)

    MqaYuYuaYuYu = MatMul2(MqaYuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = MatMul2(YuMqaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = MatMul2(aYuYuMq,aYuYu,OnlyDiagonal)

    MuYuaYuYuaYu = MatMul2(MuYuaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = MatMul2(aYuMuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = MatMul2(YuaYuMu,YuaYu,OnlyDiagonal)

    AdaAdYdaYd = MatMul2(AdaAd,YdaYd,OnlyDiagonal)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = MatMul2(AdaYd,YdaAd,OnlyDiagonal)
    YdaAdAdaYd = MatMul2(YdaAd,AdaYd,OnlyDiagonal)

    aAdAdaYdYd = MatMul2(aAdAd,aYdYd,OnlyDiagonal)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = MatMul2(aAdYd,aYdAd,OnlyDiagonal)
    aYdAdaAdYd = MatMul2(aYdAd,aAdYd,OnlyDiagonal)

    AeaAeYeaYe = MatMul2(AeaAe,YeaYe,OnlyDiagonal)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = MatMul2(AeaYe,YeaAe,OnlyDiagonal)
    YeaAeAeaYe = MatMul2(YeaAe,AeaYe,OnlyDiagonal)

    aAeAeaYeYe = MatMul2(aAeAe,aYeYe,OnlyDiagonal)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = MatMul2(aAeYe,aYeAe,OnlyDiagonal)
    aYeAeaAeYe = MatMul2(aYeAe,aAeYe,OnlyDiagonal)

    AuaAuYuaYu = MatMul2(AuaAu,YuaYu,OnlyDiagonal)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = MatMul2(AuaYu,YuaAu,OnlyDiagonal)
    YuaAuAuaYu = MatMul2(YuaAu,AuaYu,OnlyDiagonal)

    aAuAuaYuYu = MatMul2(aAuAu,aYuYu,OnlyDiagonal)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = MatMul2(aAuYu,aYuAu,OnlyDiagonal)
    aYuAuaAuYu = MatMul2(aYuAu,aAuYu,OnlyDiagonal)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * TraceY(1) + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp ) &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(2) = mH(2) * TraceY(2) + Real( cTrace(YnuMlaYnu),dp )  &
             & + Real( cTrace(aYnuMrYnu),dp ) + TraceA(2)
   traceMH2(1) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(2) + 6._dp * TraceMH2(1)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)) + TraceY(1) + TraceY(2)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4))          &
           & + 2._dp * (TraceaYA(1)+ TraceaYA(2))         &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul2(AuaYu,YuaYu,OnlyDiagonal)     &
              &                 + Matmul2(AdaYd,YdaYd,OnlyDiagonal) )   &
              &       + MatMul2(AeaYe,YeaYe,OnlyDiagonal)               &
              &       + Matmul2(aYuAu,aYdYd,OnlyDiagonal)               &
              &       + MatMul2(aYdAd,aYuYu,OnlyDiagonal) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3) + 1.6_dp * g2Mi(1) ) * TraceY(4)               &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)               &
      &   - 2.4_dp * g2Mi(1) * TraceY(1) - 30._dp * gauge2(2)**2 * Mhlf(2)    &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  sumM1 = aYeYe + aYnuYnu
  diagonal(5,1) = 2._dp * TraceY(2) + 6._dp * TraceY(4)   &
              & - 2._dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(sumM1), Mnu) + Matmul(Mnu, sumM1)  &
          & + diagonal(5,1) * Mnu
  

 !---------------
 ! 2-loop RGEs
 !---------------
  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1 + oo16pi2 * (Matmul(b_2,gauge2) - Matmul(a_2a,TraceY) ) )
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYnu = oo16pi2 * betaYnu1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,4
     sumI = sumI + a_2a(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMr = oo16pi2 * betaMr1 ! + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYnu = oo16pi2 * betaYnu1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1 * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAnu = oo16pi2 * betaAnu1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMr = oo16pi2 * betaMr1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  !----------------------
  ! left neutrino mass
  !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMr(i1,i1) = Real(DMr(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do

  Call ParametersToG3(Dgauge, DYe, DYnu, DYd, DYu, DMhlf, DAe, DAnu, DAd, DAu &
                   &, DMe, DMl, DMr, DMd, DMq, DMu, DMh, Dmue, DB, DMnu, f)

  Iname = Iname - 1

 End Subroutine rge285


 Subroutine rge356(len, T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings.
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 08.01.06: including neutrino dim 5 operator
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2
  Real(dp) :: gauge(3), gauge2(3), sumI, TraceY(6), Dgauge(3), TraceY2(4)
  Complex(dp) :: Ye(3,3), Yd(3,3), Yu(3,3), aYe(3,3), aYd(3,3), aYu(3,3)      &
    & , aYdYd(3,3), aYeYe(3,3), aYuYu(3,3), sumd1(3,3), sume1(3,3)            &
    & , betaYd1(3,3), betaYd2(3,3), betaYe1(3,3), betaYe2(3,3)                &
    & , betaYu1(3,3), betaYu2(3,3), DYe(3,3), DYd(3,3), DYu(3,3)              &
    & , aYdYdaYdYd(3,3), aYeYeaYeYe(3,3), aYuYuaYuYu(3,3)                     &
    & , aYdYdaYuYu(3,3), aYuYuaYdYd(3,3), diagonal(6,2)           &
    & , hd(2), sumu1(3,3), sumd2(3,3), sume2(3,3), sumu2(3,3), hc(4)

  Complex(dp) :: Mhlf(3),DMhlf(3)

  Complex(dp) :: Ae(3,3), Ad(3,3), Au(3,3), aAe(3,3), aAd(3,3), aAu(3,3)   &
     &  , DAe(3,3), DAd(3,3), DAu(3,3), aAdAd(3,3), aAeAe(3,3), aAuAu(3,3) &
     &  , aYdAd(3,3), aYeAe(3,3), aYuAu(3,3), TraceaYA(6), betaAd1(3,3)    &
     &  , betaAd2(3,3), betaAe1(3,3), betaAe2(3,3), betaAu1(3,3)           &
     &  , betaAu2(3,3)
  Real(dp) :: TraceA(6)
  Complex(dp) :: aYdYdaYdAd(3,3), aYdAdaYdYd(3,3), TraceAY2(5)               &
     &  , aYeYeaYeAe(3,3), aYeAeaYeYe(3,3), aYuYuaYuAu(3,3), aYuAuaYuYu(3,3) &
     &  , aYuYuaYdAd(3,3), aYuAuaYdYd(3,3), aYdYdaYuAu(3,3), aYdAdaYuYu(3,3)
   
  Complex(dp) :: Me(3,3), Ml(3,3), Md(3,3), Mq(3,3), Mu(3,3), DMe(3,3)        &
     & , DMl(3,3), DMd(3,3), DMq(3,3), DMu(3,3), YdaYd(3,3), YeaYe(3,3)       &
     & , YuaYu(3,3), MdYdaYd(3,3), MeYeaYe(3,3), MuYuaYu(3,3), YdaYdMd(3,3)   &
     & , YeaYeMe(3,3), YuaYuMu(3,3), YdMqaYd(3,3), YeMlaYe(3,3)               &
     & ,  YuMqaYu(3,3), AdaAd(3,3), AeaAe(3,3), AuaAu(3,3), betaMd1(3,3)      &
     & , betaMd2(3,3), betaMe1(3,3), betaMe2(3,3), betaMl1(3,3), betaMl2(3,3) &
     & , betaMq1(3,3), betaMq2(3,3), betaMu1(3,3), betaMu2(3,3), MqaYdYd(3,3) &
     & , MqaYuYu(3,3), aYdYdMq(3,3), aYuYuMq(3,3), aYeYeMl(3,3), MlaYeYe(3,3) &
     & , aYeMeYe(3,3), aYdMdYd(3,3), aYuMuYu(3,3)                             &
     & , YdaYdYdaYd(3,3), YeaYeYeaYe(3,3), YuaYuYuaYu(3,3), MeYeaYeYeaYe(3,3) &
     & , YeaYeYeaYeMe(3,3), YeaYeMeYeaYe(3,3), YeMlaYeYeaYe(3,3), AeaYe(3,3)  &
     & , YeaYeYeMlaYe(3,3), AeaAeYeaYe(3,3), YeaYeAeaAe(3,3), YeaAe(3,3)      &
     & , AeaYeYeaAe(3,3), YeaAeAeaYe(3,3), Tr3aAdYdaAeYe

  Complex(dp) :: Tr3aYdAdaYeAe, AdaYd(3,3), YdaAd(3,3), MlaYeYeaYeYe(3,3)     &
     & , aYeYeaYeYeMl(3,3), aYeYeMlaYeYe(3,3), aYeYeaYeMeYe(3,3)              &
     & , aYeMeYeaYeYe(3,3), aAdYd(3,3),aAeYe(3,3), aAeAeaYeYe(3,3)            &
     & , aYeYeaAeAe(3,3), aAeYeaYeAe(3,3), aYeAeaAeYe(3,3), MdYdaYdYdaYd(3,3) &
     & , YdaYdYdaYdMd(3,3), YdMqaYdYdaYd(3,3),YdaYdMdYdaYd(3,3)               &
     & , YdaYdYdMqaYd(3,3), AdaAdYdaYd(3,3), YdaYDAdaAd(3,3), AdaYdYdaAd(3,3) &
     & , YdaAdAdaYd(3,3)
  Complex(dp) :: MdYdaYuYuaYd(3,3), YdaYuYuaYdMd(3,3), YdMqaYuYuaYd(3,3)      &
     & , YdaYuYuMqaYd(3,3), YdaYuMuYuaYd(3,3), AdaAuYuaYd(3,3)                &
     & , YdaYuAuaAd(3,3), AdaYuYuaAd(3,3), YdaAuAuaYd(3,3), YdaYuYuaYd(3,3)   &
     & , Tr3aYuAu, Tr3aAuYu, YuaAu(3,3)                                       &
     & , MqaYdYdaYdYd(3,3), aYdYdaYdYDMq(3,3), aYdMdYdaYdYd(3,3)              &
     & , aYdYdMqaYdYd(3,3), aYdYdaYdMdYd(3,3), aAdAdaYdYd(3,3)                &
     & , aYdYDaAdAd(3,3), aAdYdaYdAd(3,3), aYdAdaAdYd(3,3), MqaYuYuaYuYu(3,3) &
     & , aYuYuaYuYUMq(3,3), aYuMuYuaYuYu(3,3), aYuYuMqaYuYu(3,3)              &
     & , aYuYuaYuMuYu(3,3), aAuAuaYuYu(3,3), aYuYUaAuAu(3,3), aAuYuaYuAu(3,3) &
     & , aYuAuaAuYu(3,3), aAuYu(3,3), AuaYu(3,3), YuaYdYdaYu(3,3)             &
     & , AuaYdYdaAu(3,3), YuaAdAdaYu(3,3), AuaAdYdaYu(3,3), YuaYdAdaAu(3,3)   &
     & , YuMqaYuYuaYu(3,3), YuaYuYuMqaYu(3,3), MuYuaYuYuaYu(3,3)              &
     & , YuaYuYuaYuMu(3,3), YuaYuMuYuaYu(3,3), AuaAuYuaYu(3,3)                &
     & , YuaYuAuaAu(3,3), AuaYuYuaAu(3,3), YuaAuAuaYu(3,3), MuYuaYdYdaYu(3,3) &
     & , YuaYdYdaYuMu(3,3), YuMqaYdYdaYu(3,3), YuaYdYdMqaYu(3,3)              &
     & , YuaYdMdYdaYu(3,3)

  Real(dp) :: S1, S2, sig(3), Tr3aYdYdaYeYe, Tr3aAdAdaAeAe, AbsGM2(3)         &
     & , Tr3MqaYdYd3aYDMdYd , Tr3MqaYuYu3aYuMu, Tr3aAuAu

  Real(dp) :: Mh(2), DMh(2), TraceMH1(3), TraceMH2(2), betaMH11        &
     &  , betaMH12, betaMH21, betaMH22, q

  Complex(dp) :: mue, B, Dmue, DB, TraceMue(2), TraceB(2), betaMue1, betaMue2 &
     & , betaB1, betaB2

  Complex(dp) :: Mnu(3,3), DMnu(3,3), betamnu1(3,3)

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, betaYT1, sumT1, AT, aAT     &
      & , aATAT, betaAT1, aYTAT, MlaYTYT, aYTMlYT, aYTYTMl, DYT, DAT
  Complex(dp), Dimension(3,3) :: YZ, aYZ, aYZYZ, YZaYZ, betaYZ1, sumZ1, DYZ  &
    & , AS, aAS, ASaAS, betaAS1, aYSAS, DAS                                  &
    & , AZ, aAZ, aAZAZ, AZaAZ, betaAZ1, aYZAZ, DAZ
  Complex(dp), Dimension(3,3) :: YS, aYS, aYSYS, YSaYS, betaYS1, sumS1, DYS &
    & , AZaYZ, ASaYS, aYSMdYS, MlaYZYZ, aYZYZMl, aYZMdYZ, MdYZaYZ, YZaYZMd  &
    & , MdYSaYS, YSaYSMd, YZMlaYZ, YSMdaYS
  Complex(dp) :: g2Mi(3)
  Complex(dp) :: lam1, lam2, Dlam1, Dlam2, betalam11, betalam21, lam1Alam1 &
      & , lam2Alam2, Alam1, Alam2, betaAlam11, betaAlam21, DAlam1, DAlam2
  Real(dp) :: MT(2), lam12, lam22, betaMT2(2), Alam12, Alam22, b_1a(3), DMT2(2) &
      & , MS(2), betaMS2(2), DMS2(2), MZ(2), betaMZ2(2), DMZ2(2), MS15, MZ15    &
      & , MT15, DMS15, DMZ15, DMT15, betaMS15, betaMZ15, betaMT15, b_2a(3,3)

  Iname = Iname + 1
  NameOfUnit(Iname) = 'rge356'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters5(gy, gauge, Ye, YT, Yd, Yu, YZ, YS, lam1, lam2, Mhlf    &
                & , Ae, AT, Ad, Au, AZ, AS, Alam1, Alam2, Me, Ml, Md, Mq, Mu &
                & , Mh, mT, mZ, mS, MT15, MZ15, MS15, mue, B, Mnu)

  gauge2 = gauge**2
  AbsGM2 = gauge2 * Abs( Mhlf )**2
!-----------------
! beta functions
!-----------------
  lam12 = Abs(lam1)**2
  lam22 = Abs(lam2)**2

  Call Adjungate(Yd,aYd)
  Call Adjungate(YT,aYT)
  Call Adjungate(Ye,aYe)
  Call Adjungate(Yu,aYu)
  Call Adjungate(YZ,aYZ)
  Call Adjungate(YS,aYS)

  aYeYe = Matmul2(aYe,Ye,OnlyDiagonal)
  aYTYT = Matmul2(aYT,YT,OnlyDiagonal)
  aYdYd = Matmul2(aYd,Yd,OnlyDiagonal)
  aYuYu = Matmul2(aYu,Yu,OnlyDiagonal)
  aYZYZ = Matmul2(aYZ,YZ,OnlyDiagonal)
  aYSYS = Matmul2(aYS,YS,OnlyDiagonal)

  YdaYd = Matmul2(Yd,aYd,OnlyDiagonal)
  YZaYZ = Matmul2(YZ,aYZ,OnlyDiagonal)
  YSaYS = Matmul2(YS,aYS,OnlyDiagonal)

  TraceY(1) = Real( cTrace(aYeYe),dp )
  TraceY(2) = Real( cTrace(aYTYT),dp )
  TraceY(3) = Real( cTrace(aYdYd),dp )
  TraceY(4) = Real( cTrace(aYuYu),dp )
  TraceY(5) = Real( cTrace(aYZYZ),dp )
  TraceY(6) = Real( cTrace(aYSYS),dp )

  diagonal(1,1) = 3._dp * (TraceY(3) + lam12) + TraceY(1)    &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * (aYeYe + aYTYT + aYZYZ)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = aYeYe + 6._dp * aYTYT + 3._dp * aYZYZ
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul2(YT,sumT1,OnlyDiagonal)                 &
        & + Matmul2(Transpose(aYeYe+3._dp * aYZYZ),YT,OnlyDiagonal)

  diagonal(3,1) = 3._dp * (TraceY(3)  + lam12) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)  &
        & + 2._dp * Matmul2(YZaYZ + 2._dp * YSaYS, Yd,OnlyDiagonal)

  diagonal(4,1) = 3._dp * (TraceY(4) + lam22)                &
   &  + c1_1(3,1) * gauge2(1) + c1_1(3,2) * gauge2(2) + c1_1(3,3) * gauge2(3)
  sumu1  = 3._dp * aYuYu + aYdYd
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do

  betaYu1 = Matmul2(Yu,sumu1,OnlyDiagonal)

  diagonal(5,1) = TraceY(5) + c1_1(2,1) * gauge2(1) &
            &   + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumZ1 = aYeYe + 3._dp * aYTYT + 5._dp * aYZYZ
  Do i1=1,3
   sumZ1(i1,i1) = sumZ1(i1,i1) + diagonal(5,1)
  End Do

  betaYZ1 = Matmul2(YZ,sumZ1,OnlyDiagonal)                 &
        & + 2._dp * Matmul2(YdaYd+2._dp * aYSYS,YZ,OnlyDiagonal)

  diagonal(6,1) = TraceY(6) - 0.8_dp * gauge2(1) - 12._dp * gauge2(3)
  sumS1 = 2._dp * Transpose(YdaYd + YZaYZ) + 8._dp * aYSYS
  Do i1=1,3
   sumS1(i1,i1) = sumS1(i1,i1) + diagonal(6,1)
  End Do

  betaYS1 = Matmul2(YS,sumS1,OnlyDiagonal)                 &
        & + 2._dp * Matmul2(YdaYd + YZaYZ,YS,OnlyDiagonal)


  betalam11 = lam1 * (7._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
            &        + 6._dp * TraceY(3)                            &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  betalam21 = lam2 * (7._dp * lam22  + 6._dp * TraceY(4)     &
            &        - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdYd = Matmul2(aYdYd,aYdYd,OnlyDiagonal)
   aYeYeaYeYe = Matmul2(aYeYe,aYeYe,OnlyDiagonal)
   aYuYuaYuYu = Matmul2(aYuYu,aYuYu,OnlyDiagonal)
   aYuYuaYdYd = Matmul2(aYuYu,aYdYd,OnlyDiagonal)
   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)

   TraceY2(1) = Real( cTrace(aYeYeaYeYe), dp)
   TraceY2(2) = Real( cTrace(aYdYdaYdYd), dp)
   TraceY2(3) = Real( cTrace(aYuYuaYuYu), dp)
   TraceY2(4) = Real( cTrace(aYdYdaYuYu), dp)

   diagonal(1,2) = - 3._dp * (3._dp * TraceY2(2) + TraceY2(4) + TraceY2(1) ) &
             &   + ( 16._dp * gauge2(3) - 0.4_dp * gauge2(1) ) * TraceY(3)   &
             &   + 1.2_dp * gauge2(1) * TraceY(1)                            &
             &   + ( 7.5_dp * gauge2(2) + 1.8_dp * gauge2(1) ) * gauge2(2)   &
             &   + 13.5_dp * gauge2(1)**2
   hd(1) = 9._dp * TraceY(3) + 3._dp * TraceY(1) - 6._dp * gauge2(2)
   sume2 = - 4._dp * aYeYeaYeYe - hd(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
 
   betaYe2 = Matmul2(Ye,sume2,OnlyDiagonal)
    
   diagonal(2,2) = diagonal(1,2)                                        &
      &     + 8._dp * ( ( gauge2(1) - 2._dp * gauge2(3) ) / 9._dp       &
      &              + gauge2(2)  ) * gauge2(3)                         &
      &     - 0.8_dp * gauge2(1) * gauge2(2)                            &
      &     - 928._dp * gauge2(1)**2 / 90._dp
   hd(1) = 0.8_dp * gauge2(1) - 3._dp * TraceY(4)
   hd(2) = 9._dp * TraceY(3) + 3._dp * TraceY(1)     &
     &   - 6._dp * gauge2(2) - 0.8_dp * gauge2(1)
   sumd2 = - 4._dp * aYdYdaYdYd - 2._dp * aYuYuaYuYu - 2._dp * aYuYuaYdYd &
       & + hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(2,2)
   End Do
 
   betaYd2 = Matmul2(Yd,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = - 3._dp * (3._dp * TraceY2(3) + TraceY2(4) )            &
     &      + ( 16._dp * gauge2(3) + 0.8_dp * gauge2(1) ) * TraceY(4)      &
     &      + 8._dp * ( ( 3.4_dp * gauge2(1) - 2._dp* gauge2(3) ) / 9._dp  &
     &               + gauge2(2)  ) * gauge2(3)                            &
     &      + ( 7.5_dp * gauge2(2) + gauge2(1) ) * gauge2(2)               &
     &      + 2743._dp * gauge2(1)**2 / 450._dp
   hd(1) = 9._dp * TraceY(4) - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   hd(2) = 3._dp * TraceY(3) + TraceY(1) - 0.4_dp * gauge2(1)
   sumu2 = - 4._dp * aYuYuaYuYu - 2._dp * aYdYdaYdYd - 2._dp * aYdYdaYuYu  &
       & - hd(1) * aYuYu - hd(2) * aYdYd
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(3,2)
   End Do
 
   betaYu2 = Matmul2(Yu,sumu2,OnlyDiagonal)
    
  End If 

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ae,aAe)
  Call Adjungate(AT,aAT)
  Call Adjungate(Ad,aAd)
  Call Adjungate(Au,aAu)
  Call Adjungate(AZ,aAZ)
  Call Adjungate(AS,aAS)

  aAdAd = MatMul2(aAd,Ad,OnlyDiagonal)
  aATAT = MatMul2(aAT,AT,OnlyDiagonal)
  aAeAe = MatMul2(aAe,Ae,OnlyDiagonal)
  aAuAu = MatMul2(aAu,Au,OnlyDiagonal)
  aAZAZ = MatMul2(aAZ,AZ,OnlyDiagonal)

  aYdAd = MatMul2(aYd,Ad,OnlyDiagonal)
  aYTAT = MatMul2(aYT,AT,OnlyDiagonal)
  aYeAe = MatMul2(aYe,Ae,OnlyDiagonal)
  aYuAu = MatMul2(aYu,Au,OnlyDiagonal)
  aYZAZ = MatMul2(aYZ,AZ,OnlyDiagonal)
  aYSAS = MatMul2(aYS,AS,OnlyDiagonal)

  AeaYe = MatMul2(Ad,aYd,OnlyDiagonal)
  AdaYd = MatMul2(Ad,aYd,OnlyDiagonal)
  AZaYZ = MatMul2(AZ,aYZ,OnlyDiagonal)
  ASaYS = MatMul2(AS,aYS,OnlyDiagonal)

  ASaAS = MatMul2(AS,aAS,OnlyDiagonal)

  TraceA(1) = Real( cTrace(aAeAe),dp )
  TraceA(2) = Real( cTrace(aATAT),dp )
  TraceA(3) = Real( cTrace(aAdAd),dp )
  TraceA(4) = Real( cTrace(aAuAu),dp )
  TraceA(5) = Real( cTrace(aAZAZ),dp )
  TraceA(6) = Real( cTrace(ASaAS),dp )

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYTAT) 
  TraceaYA(3) = cTrace(aYdAd) 
  TraceaYA(4) = cTrace(aYuAu) 
  TraceaYA(5) = cTrace(aYZAZ) 
  TraceaYA(6) = cTrace(aYSAS) 
  lam1Alam1 = Conjg(lam1) * Alam1
  lam2Alam2 = Conjg(lam2) * Alam2
  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = MatMul2(Ae,sume1,OnlyDiagonal)
  
  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1)  &
                &         + 3._dp * lam1Alam1                  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    ) 
  sume1 = 4._dp * aYeAe + 6._dp * (aYTAT + aYZAZ)
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + MatMul2(Ye,sume1,OnlyDiagonal)

  !--------------
  ! A_T
  !--------------
  diagonal(2,1) = TraceY(2)  + lam12          &
            &   - 1.8_dp * gauge2(1) - 7._dp * gauge2(2)
  sumT1 = 2._dp * aYeYe + 9._dp * aYTYT + 3._dp * aYZYZ
  betaAT1 = MatMul2(AT,sumT1,OnlyDiagonal)
  betaAT1 = Transpose(betaAT1)
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaAT1 = betaAT1 + MatMul2(AT,sumT1,OnlyDiagonal)

  diagonal(2,1) = 2._dp * ( TraceaYA(2) + lam1Alam1    &
      &         + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) )
  betaAT1 = betaAT1 + diagonal(2,1) * YT                         &
        & + 6._dp * ( MatMul2(YT, aYZAZ,OnlyDiagonal)            &
        &           + MatMul2(Transpose(aYZAZ), YT,OnlyDiagonal) )

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = MatMul2(Ad,sumd1,OnlyDiagonal)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3) + TraceaYA(1) &
                &         + 3._dp * lam1Alam1                 &
                &         - c1_1(2,1) * g2Mi(1)   &
                &         - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(3,1)
  End Do
  betaAd1 = betaAd1 + MatMul2(Yd,sumd1,OnlyDiagonal)            &
        & + 2._dp * MatMul2(2._dp*YSaYS+YZaYZ, Ad,OnlyDiagonal) &
        & + 4._dp * MatMul2(2._dp*ASaYS+AZaYZ, Yd,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = MatMul2(Au,sumu1,OnlyDiagonal)
  
  diagonal(4,1) = 2._dp * ( 3._dp * TraceaYA(4) + 3._dp * lam2Alam2   &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2) &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(4,1)
  End Do
  betaAu1 = betaAu1 + MatMul2(Yu,sumu1,OnlyDiagonal)

  !--------------
  ! A_Z
  !--------------
  sumZ1 = sumZ1 + 3._dp * aYZYZ
  betaAZ1 = MatMul2(AZ,sumZ1,OnlyDiagonal)                                   &
    & + 2._dp * MatMul2(2._dp*YSaYs + YdaYd + 2._dp*YZaYZ, AZ,OnlyDiagonal)  &
    & + 6._dp * MatMul2(YZ, aYTAT, OnlyDiagonal)

  diagonal(5,1) = - 2._dp * ( c1_1(2,1) * g2Mi(1) + c1_1(2,2) * g2Mi(2) &
                &           + c1_1(2,3) * g2Mi(3)  ) + 2._dp * TraceaYA(5)

  betaAZ1 = betaAZ1 + diagonal(5,1) * YZ                                                &
        & +  2._dp * MatMul2(YZ, AeaYe,OnlyDiagonal)                              &
        & + 4._dp * MatMul2(AdaYd + 2._dp * ASaYS, YZ,OnlyDiagonal)

  !--------------
  ! A_S
  !--------------
  diagonal(6,1) = TraceY(6) - 0.8_dp * gauge2(1) - 12._dp * gauge2(3)
  sumS1 = 2._dp * Transpose(YdaYd + YZaYZ) + 12._dp * aYSYS
  betaAS1 = MatMul2(AS,sumS1,OnlyDiagonal)
  betaAS1 = Transpose(betaAS1)
  Do i1=1,3
   sumS1(i1,i1) = sumS1(i1,i1) + diagonal(6,1)
  End Do

  betaAS1 = betaAS1 + MatMul2(AS,sumS1,OnlyDiagonal)

  diagonal(6,1) = 2._dp * ( TraceaYA(6) + 0.8_dp * g2Mi(1) + 12._dp * g2Mi(3) )
  sumS1 = 4._dp * (AdaYd + AZaYZ) 
  betaAS1 = betaAS1 + diagonal(6,1) * YS               &
        & + MatMul2(sumS1, YS,OnlyDiagonal)            &
        & + MatMul2(YS,Transpose(sumS1), OnlyDiagonal) 

  !--------------
  ! A_1
  !--------------
  betaAlam11 = Alam1 * (21._dp * lam12 + TraceY(2) + 2._dp * TraceY(1) &
          &           + 6._dp * TraceY(3)                              &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )       &
          & + 2._dp * lam1 * ( TraceaYA(2) + 2._dp * TraceaYA(1)       &
          &                  + 6._dp * TraceaYA(3)                     &
          &                  + 1.8_dp * g2Mi(1)  + 7._dp * g2Mi(2) ) 

  !--------------
  ! A_2
  !--------------
  betaAlam21 = Alam2 * (21._dp * lam22 + 6._dp * TraceY(4)          &
          &           - 1.8_dp * gauge2(1) - 7._dp * gauge2(2) )    &
          & + 2._dp * lam2 * ( 6._dp * TraceaYA(4)                  &
          &                  + 1.8_dp * g2Mi(1) + 7._dp * g2Mi(2) ) 

  If (TwoLoopRGE) Then
   aYdYdaYdAd = MatMul2(aYdYd,aYdAd,OnlyDiagonal)
   aYdAdaYdYd = MatMul2(aYdAd,aYdYd,OnlyDiagonal)
   aYeYeaYeAe = MatMul2(aYeYe,aYeAe,OnlyDiagonal)
   aYeAeaYeYe = MatMul2(aYeAe,aYeYe,OnlyDiagonal)
   aYuYuaYuAu = MatMul2(aYuYu,aYuAu,OnlyDiagonal)
   aYuAuaYuYu = MatMul2(aYuAu,aYuYu,OnlyDiagonal)
   aYuAuaYdYd = MatMul2(aYuAu,aYdYd,OnlyDiagonal)
   aYuYuaYdAd = MatMul2(aYuYu,aYdAd,OnlyDiagonal)
   aYdAdaYuYu = MatMul2(aYdAd,aYuYu,OnlyDiagonal)
   aYdYdaYuAu = MatMul2(aYdYd,aYuAu,OnlyDiagonal)
   TraceAY2(1) = cTrace(aYeYeaYeAe)
   TraceAY2(2) = cTrace(aYdYdaYdAd)
   TraceAY2(3) = cTrace(aYuYuaYuAu)
   TraceAY2(4) = cTrace(aYuYuaYdAd)
   TraceAY2(5) = cTrace(aYdYdaYuAu)

  !--------------
  ! A_e
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = MatMul2(Ae,sume2,OnlyDiagonal)
    
   diagonal(1,2) = -6._dp * ( 6._dp * TraceAY2(2) + TraceAY2(4)       &
     &                      + TraceAY2(5) + 2._dp * TraceAY2(1)  )    &
     &  + ( 32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)   &
     &  - ( 32._dp * g2Mi(3)                              &
     &    - 0.8_dp * g2Mi(1) ) * TraceY(3)                &
     &  + 2.4_dp * gauge2(1) * TraceaYA(1)                            &
     &  - 2.4_dp * g2Mi(1) * TraceY(1)                    &
     &  - ( 30._dp * g2Mi(2)                              &
     &    + 3.6_dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)      &
     &  - 54._dp * gauge2(1)**2 * Mhlf(1)
   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)     &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1) &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2)
   sume2 = - 6._dp * aYeYeaYeAe - 8._dp * aYeAeaYeYe  &
         & - hd(1) * aYeAe - hc(1) * aYeYe
   Do i1=1,3
    sume2(i1,i1) = sume2(i1,i1) + diagonal(1,2)
   End Do
   betaAe2 = betaAe2 + MatMul2(Ye,sume2,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = MatMul2(Ad,sumd2,OnlyDiagonal)
    
   diagonal(3,2) = diagonal(1,2)                                   &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                 &
     &              - gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp      &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)    &
     &  + 1.6_dp * gauge2(1) * gauge2(2) * (Mhlf(1)+Mhlf(2))       & 
     &  + 1.856e3_dp * gauge2(1)**2 * Mhlf(1) / 4.5e1_dp

   hd(1) = 12._dp * TraceY(3) + 4._dp * TraceY(1)  &
       & - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hd(2) = 6._dp * TraceY(4) - 1.6_dp * gauge2(1) 
   hc(1) = 18._dp * TraceaYA(3) + 6._dp * TraceaYA(1)                 &
       & + 1.2e1_dp * gauge2(2) *  Mhlf(2) + 1.6_dp * gauge2(1) *  Mhlf(1)
   hc(2) = 6._dp * TraceaYA(4) + 1.6_dp * gauge2(1) *  Mhlf(1)
   sumd2 = - 6._dp * aYdYdaYdAd - 8._dp * aYdAdaYdYd                  &
       &   - 4._dp * ( aYuAuaYuYu + aYuYuaYuAu + aYuAuaYdYd )         &
       &   - 2._dp * aYuYuaYdAd - hd(1) * aYdAd - hc(1) * aYdYd       &
       &  - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumd2(i1,i1) = sumd2(i1,i1) + diagonal(3,2)
   End Do
   betaAd2 = betaAd2 + MatMul2(Yd,sumd2,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(4) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = MatMul2(Au,sumu2,OnlyDiagonal)
    
   diagonal(4,2) =  -6._dp * ( 6._dp * TraceAY2(3) + TraceAY2(4)        &
     &                        + TraceAY2(5)  )                          &
     &  + ( 32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)     &
     &  - ( 32._dp * g2Mi(3)                                &
     &    + 1.6_dp * g2Mi(1) ) * TraceY(4)                  &
     &  + 16._dp * ( ( 4._dp * g2Mi(3)                      &
     &              - 3.4_dp * gauge2(1) * (Mhlf(3)+Mhlf(1)) ) / 9._dp  &
     &            - gauge2(2) * (Mhlf(3)+Mhlf(2)) ) * gauge2(3)         &
     &  - ( 30._dp * g2Mi(2)                                &
     &    + 2._dp * gauge2(1) * (Mhlf(1)+Mhlf(2)) ) * gauge2(2)         &
     &  - 5486._dp * gauge2(1)**2 * Mhlf(1) / 225._dp
   hd(1) = 6._dp * TraceY(3) + 2._dp * TraceY(1) - 0.8_dp * gauge2(1)
   hc(1) = 6._dp * TraceaYA(3) + 2._dp * TraceaYA(1)   &
       &  + 0.8_dp * gauge2(1) *  Mhlf(1)
   hd(2) = 12._dp * TraceY(4) - 6._dp * gauge2(2) - 1.2_dp * gauge2(1)
   hc(2) = 18._dp * TraceaYA(4) + 1.2e1_dp * g2Mi(2) & 
       & + 0.8_dp * g2Mi(1)
   sumu2 = - 6._dp * aYuYuaYuAu - 8._dp * aYuAuaYuYu                     &
       &   - 4._dp * ( aYdAdaYdYd + aYdYdaYdAd + aYdAdaYuYu )            &
       &   - 2._dp * aYdYdaYuAu - hd(1) * aYdAd - hc(1) * aYdYd          &
       &   - hd(2) * aYuAu - hc(2) * aYuYu
   Do i1=1,3
    sumu2(i1,i1) = sumu2(i1,i1) + diagonal(4,2)
   End Do
   betaAu2 = betaAu2 + MatMul2(Yu,sumu2,OnlyDiagonal)

  End If 
!----------------------------------------------
! beta functions for Sfermion mass parameters
!----------------------------------------------
   S1 = mH(2) - mH(1) + 3._dp * (mT(1) - mT(2) ) &
    & + mZ(1) - mZ(2) + 4._dp * (mS(2) - mS(1))
   Do i1=1,3
    S1 = S1 + Me(i1,i1) - Ml(i1,i1) &
       &    + Md(i1,i1) + Mq(i1,i1) - 2._dp * Mu(i1,i1)
   End Do

   S1 = S1 * gauge2(1)

   YdaYd = MatMul2(Yd,aYd,OnlyDiagonal)
   YeaYe = MatMul2(Ye,aYe,OnlyDiagonal)
   YuaYu = MatMul2(Yu,aYu,OnlyDiagonal)

   MeYeaYe = MatMul2(Me,YeaYe,OnlyDiagonal)
   MlaYeYe = MatMul2(Ml,aYeYe,OnlyDiagonal)
   MlaYTYT = MatMul2(Ml,aYTYT,OnlyDiagonal)
   MlaYZYZ = MatMul2(Ml,aYZYZ,OnlyDiagonal)

   MdYdaYd = MatMul2(Md,YdaYd,OnlyDiagonal)
   MqaYdYd = MatMul2(Mq,aYdYd,OnlyDiagonal)
   MqaYuYu = MatMul2(Mq,aYuYu,OnlyDiagonal)
   MuYuaYu = MatMul2(Mu,YuaYu,OnlyDiagonal)
   MdYZaYZ = MatMul2(Md,YZaYZ,OnlyDiagonal)
   MdYSaYS = MatMul2(Md,YSaYS,OnlyDiagonal)

   YeaYeMe = MatMul2(YeaYe,Me,OnlyDiagonal)
   aYeYeMl = MatMul2(aYeYe,Ml,OnlyDiagonal)
   aYTYTMl = MatMul2(aYTYT,Ml,OnlyDiagonal)
   aYZYZMl = MatMul2(aYZYZ,Ml,OnlyDiagonal)

   YdaYdMd = MatMul2(YdaYd,Md,OnlyDiagonal)
   aYdYdMq = MatMul2(aYdYd,Mq,OnlyDiagonal)
   aYuYuMq = MatMul2(aYuYu,Mq,OnlyDiagonal)
   YuaYuMu = MatMul2(YuaYu,Mu,OnlyDiagonal)
   YZaYZMd = MatMul2(YZaYZ,Md,OnlyDiagonal)
   YSaYSMd = MatMul2(YSaYS,Md,OnlyDiagonal)

   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   aYTMlYT = MatMul3(aYT,Transpose(Ml),YT,OnlyDiagonal)
   YZMlaYZ = MatMul3(YZ,Ml,YZ,OnlyDiagonal)
   YSMdaYS = MatMul3(YS,Transpose(Md),aYS,OnlyDiagonal)
   aYSMdYS = MatMul3(aYS,Md,YS,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)
   aYZMdYZ = MatMul3(aYZ,Md,YZ,OnlyDiagonal)

   AeaAe = MatMul2(Ae,aAe,OnlyDiagonal)
   AdaAd = MatMul2(Ad,aAd,OnlyDiagonal)
   AuaAu = MatMul2(Au,aAu,OnlyDiagonal)
   AZaAZ = MatMul2(AZ,aAZ,OnlyDiagonal)
   Alam12 = Abs(Alam1)**2
   Alam22 = Abs(Alam2)**2

   diagonal(1,1) = - 4.8_dp * AbsGM2(1) + 1.2_dp * S1
   betaMe1 = 2._dp * (MeYeaYe + YeaYeMe)             &
         & + 4._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe )
   Do i1=1,3
    betaMe1(i1,i1) = betaMe1(i1,i1) + diagonal(1,1)
   End Do

   diagonal(3,1) = - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1
   betaMl1 = MlaYeYe + aYeYeMl + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe )  &
         & + 3._dp * ( MlaYTYT + MlaYZYZ + aYTYTMl + aYZYZMl )                &
         & + 6._dp * ( aYTMlYT + aATAT + MT(1) * aYTYT)                       &
         & + 6._dp * ( aYZMdYZ + aAZAZ + MZ(1) * aYZYZ)

   Do i1=1,3
    betaMl1(i1,i1) = betaMl1(i1,i1) + diagonal(3,1)
   End Do

   diagonal(4,1) = - ( 32._dp * AbsGM2(3) + 1.6_dp * AbsGM2(1) ) / 3._dp &
               & + 0.4_dp * S1
   betaMd1 = 2._dp * (MdYdaYd + YdaYdMd + MdYZaYZ + YZaYZMd )   &
         & + 4._dp * (MdYSaYS + YSaYSMd)                        &
         & + 4._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd )        &
         & + 4._dp * ( mZ(1) * YZaYZ + YZMlaYZ + AZaAZ )        &
         & + 8._dp * ( mS(1) * YSaYS + YSMdaYS + ASaAS )
   Do i1=1,3
    betaMd1(i1,i1) = betaMd1(i1,i1) + diagonal(4,1)
   End Do

   diagonal(5,1) = - ( 32._dp * AbsGM2(3) + 0.4_dp * AbsGM2(1) ) / 3._dp &
               & - 6._dp * AbsGM2(2) + 0.2_dp * S1
   betaMq1 = MqaYuYu + aYuYuMq + MqaYdYd + aYdYdMq             &
         & + 2._dp * ( mH(2) * aYuYu + mH(1) * aYdYd + aYuMuYu     &
         &           + aYdMdYd + aAuAu + aAdAd )
   Do i1=1,3
    betaMq1(i1,i1) = betaMq1(i1,i1) + diagonal(5,1)
   End Do

   diagonal(6,1) = - ( 32._dp * AbsGM2(3) + 6.4_dp * AbsGM2(1) ) / 3._dp &
               &   - 0.8_dp * S1
   betaMu1 = 2._dp * (MuYuaYu + YuaYuMu)             &
         & + 4._dp * ( mH(2) * YuaYu + YuMqaYu + AuaAu )
   Do i1=1,3
    betaMu1(i1,i1) = betaMu1(i1,i1) + diagonal(6,1)
   End Do

   If (TwoLoopRGE) Then
    YdaYdYdaYd = MatSquare(YdaYd,OnlyDiagonal)
    YeaYeYeaYe = MatSquare(YeaYe,OnlyDiagonal)
    YuaYuYuaYu = MatSquare(YuaYu,OnlyDiagonal)

    AuaYu = MatMul2(Au,aYu,OnlyDiagonal)

    aAdYd = MatMul2(aAd,Yd,OnlyDiagonal)
    aAeYe = MatMul2(aAe,Ye,OnlyDiagonal)
    aAuYu = MatMul2(aAu,Yu,OnlyDiagonal)

    YdaAd = MatMul2(Yd,aAd,OnlyDiagonal)
    YeaAe = MatMul2(Ye,aAe,OnlyDiagonal)
    YuaAu = MatMul2(Yu,aAu,OnlyDiagonal)

    YdaYuYuaYd = MatMul3(Yd,aYuYu,aYd,OnlyDiagonal)
    AdaYuYuaAd = MatMul3(Ad,aYuYu,aAd,OnlyDiagonal)
    YdaAuAuaYd = MatMul3(Yd,aAuAu,aYd,OnlyDiagonal)
    AdaAuYuaYd = MatMul4(Ad,aAu,Yu,aYd,OnlyDiagonal)
    YdaYuAuaAd = MatMul3(Yd,aYuAu,aAd,OnlyDiagonal)

    YuaYdYdaYu = MatMul3(Yu,aYdYd,aYu,OnlyDiagonal)
    AuaYdYdaAu = MatMul3(Au,aYdYd,aAu,OnlyDiagonal)
    YuaAdAdaYu = MatMul3(Yu,aAdAd,aYu,OnlyDiagonal)
    AuaAdYdaYu = MatMul4(Au,aAd,Yd,aYu,OnlyDiagonal)
    YuaYdAdaAu = MatMul3(Yu,aYdAd,aAu,OnlyDiagonal)

    MdYdaYuYuaYd = MatMul2(Md,YdaYuYuaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYuYuaYd, YdaYuYuaYdMd)
    YdMqaYuYuaYd = MatMul3(Yd,MqaYuYu,aYd,OnlyDiagonal)
    Call Adjungate(YdMqaYuYuaYd, YdaYuYuMqaYd)
    YdaYuMuYuaYd = MatMul3(Yd,aYuMuYu,aYd,OnlyDiagonal)

    MuYuaYdYdaYu = MatMul2(Mu,YuaYdYdaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYdYdaYu, YuaYdYdaYuMu)
    YuMqaYdYdaYu = MatMul3(Yu,MqaYdYd,aYu,OnlyDiagonal)
    Call Adjungate(YuMqaYdYdaYu, YuaYdYdMqaYu)
    YuaYdMdYdaYu = MatMul3(Yu,aYdMdYd,aYu,OnlyDiagonal)

    MeYeaYeYeaYe = MatMul2(MeYeaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(MeYeaYeYeaYe,YeaYeYeaYeMe)
    aYeMeYeaYeYe = MatMul2(aYeMeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(aYeMeYeaYeYe,aYeYeaYeMeYe)
    YeaYeMeYeaYe = MatMul2(YeaYeMe,YeaYe,OnlyDiagonal)

    MlaYeYeaYeYe = MatMul2(MlaYeYe,aYeYe,OnlyDiagonal)
    Call Adjungate(MlaYeYeaYeYe, aYeYeaYeYeMl)
    YeMlaYeYeaYe = MatMul2(YeMlaYe,YeaYe,OnlyDiagonal)
    Call Adjungate(YeMlaYeYeaYe, YeaYeYeMlaYe)
    aYeYeMlaYeYe = MatMul2(aYeYeMl,aYeYe,OnlyDiagonal)

    MdYdaYdYdaYd = MatMul2(MdYdaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(MdYdaYdYdaYd, YdaYdYdaYdMd)
    aYdMdYdaYdYd = MatMul2(aYdMdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(aYdMdYdaYdYd, aYdYdaYdMdYd)
    YdaYdMdYdaYd = MatMul2(YdaYdMd,YdaYd,OnlyDiagonal)

    MqaYdYdaYdYd = MatMul2(MqaYdYd,aYdYd,OnlyDiagonal)
    Call Adjungate(MqaYdYdaYdYd, aYdYdaYdYdMq)
    YdMqaYdYdaYd = MatMul2(YdMqaYd,YdaYd,OnlyDiagonal)
    Call Adjungate(YdMqaYdYdaYd, YdaYdYdMqaYd)
    aYdYdMqaYdYd = MatMul2(aYdYdMq,aYdYd,OnlyDiagonal)

    MqaYuYuaYuYu = MatMul2(MqaYuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(MqaYuYuaYuYu, aYuYuaYuYuMq)
    YuMqaYuYuaYu = MatMul2(YuMqaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(YuMqaYuYuaYu, YuaYuYuMqaYu)
    aYuYuMqaYuYu = MatMul2(aYuYuMq,aYuYu,OnlyDiagonal)

    MuYuaYuYuaYu = MatMul2(MuYuaYu,YuaYu,OnlyDiagonal)
    Call Adjungate(MuYuaYuYuaYu, YuaYuYuaYuMu)
    aYuMuYuaYuYu = MatMul2(aYuMuYu,aYuYu,OnlyDiagonal)
    Call Adjungate(aYuMuYuaYuYu, aYuYuaYuMuYu)
    YuaYuMuYuaYu = MatMul2(YuaYuMu,YuaYu,OnlyDiagonal)

    AdaAdYdaYd = MatMul2(AdaAd,YdaYd,OnlyDiagonal)
    Call Adjungate(AdaAdYdaYd, YdaYdAdaAd)
    AdaYdYdaAd = MatMul2(AdaYd,YdaAd,OnlyDiagonal)
    YdaAdAdaYd = MatMul2(YdaAd,AdaYd,OnlyDiagonal)

    aAdAdaYdYd = MatMul2(aAdAd,aYdYd,OnlyDiagonal)
    Call Adjungate(aAdAdaYdYd, aYdYdaAdAd)
    aAdYdaYdAd = MatMul2(aAdYd,aYdAd,OnlyDiagonal)
    aYdAdaAdYd = MatMul2(aYdAd,aAdYd,OnlyDiagonal)

    AeaAeYeaYe = MatMul2(AeaAe,YeaYe,OnlyDiagonal)
    Call Adjungate(AeaAeYeaYe, YeaYeAeaAe)
    AeaYeYeaAe = MatMul2(AeaYe,YeaAe,OnlyDiagonal)
    YeaAeAeaYe = MatMul2(YeaAe,AeaYe,OnlyDiagonal)

    aAeAeaYeYe = MatMul2(aAeAe,aYeYe,OnlyDiagonal)
    Call Adjungate(aAeAeaYeYe, aYeYeaAeAe)
    aAeYeaYeAe = MatMul2(aAeYe,aYeAe,OnlyDiagonal)
    aYeAeaAeYe = MatMul2(aYeAe,aAeYe,OnlyDiagonal)

    AuaAuYuaYu = MatMul2(AuaAu,YuaYu,OnlyDiagonal)
    Call Adjungate(AuaAuYuaYu, YuaYuAuaAu)
    AuaYuYuaAu = MatMul2(AuaYu,YuaAu,OnlyDiagonal)
    YuaAuAuaYu = MatMul2(YuaAu,AuaYu,OnlyDiagonal)

    aAuAuaYuYu = MatMul2(aAuAu,aYuYu,OnlyDiagonal)
    Call Adjungate(aAuAuaYuYu, aYuYuaAuAu)
    aAuYuaYuAu = MatMul2(aAuYu,aYuAu,OnlyDiagonal)
    aYuAuaAuYu = MatMul2(aYuAu,aAuYu,OnlyDiagonal)

    S2 = (1.5_dp * gauge2(2) + 0.3_dp * gauge2(1) )            &
     &      * (MH(2) - MH(1) - Real(cTrace(ML),dp) )              &
     & + ( (8._dp * gauge2(3) + 0.1_dp*gauge2(1)) / 3._dp      &
     &   + 1.5_dp * gauge2(2) ) * Real( cTrace(Mq),dp   )         &
     & - (16._dp * gauge2(3) + 3.2_dp*gauge2(1) )              &
     &    * Real(cTrace(Mu),dp) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md),dp) / 3._dp                           &
     & + 1.2_dp*gauge2(1) * Real(cTrace(Me),dp)                   &
     & - 3._dp * (MH(2)*TraceY(4) - MH(1) * TraceY(3) )        &
     & + mH(1) * TraceY(1)
    
    Do i1=1,3
     S2 = S2 - YuMqaYu(i1,i1) + 4._dp * aYuMuYu(i1,i1)   &
        &    - YdMqaYd(i1,i1) - 2._dp * aYdMdYd(i1,i1)   &
        &    + YeMlaYe(i1,i1) - 2._dp * aYeMeYe(i1,i1)
    End Do

    sig(1) = 3._dp * (MH(1) + MH(2) + Real(cTrace(Ml),dp) )          &
         & + Real(cTrace(Mq),dp) + 8._dp * Real(cTrace(Mu),dp)          &
         & + 2._dp * Real(cTrace(Md),dp) + 6._dp * Real(cTrace(Me),dp) 
    sig(1) = 0.2_dp * gauge2(1) * sig(1)
    sig(2) = gauge2(2) * ( MH(1) + MH(2) + Real( cTrace(Ml),dp )    &
           &             + 3._dp * Real( cTrace(Mq),dp ) )
    sig(3) = gauge2(3) * ( 2._dp * Real( cTrace(Mq),dp )            &
           &             + Real(cTrace(Mu),dp) + Real(cTrace(Md),dp) )

    Tr3aYdYdaYeYe = 3._dp * TraceY(3) + TraceY(1)
    Tr3MqaYdYd3aYDMdYd = Real(cTrace(MlaYeYe),dp) + Real(cTrace(aYeMeYe),dp)  &
            & + 3._dp * ( Real(cTrace(MqaYdYd),dp) + Real(cTrace(aYdMdYd),dp) )
    Tr3aAdAdaAeAe = 3._dp * TraceA(3) + TraceA(1)
    Tr3aYdAdaYeAe = 3._dp * TraceaYA(3) + TraceaYA(1)
    Tr3aAdYdaAeYe = Conjg( Tr3aYdAdaYeAe )
    Tr3MqaYuYu3aYuMu = 3._dp * (Real(cTrace(MqaYuYu),dp) &
                               + Real(cTrace(aYuMuYu),dp) )
    Tr3aAuAu = 3._dp * TraceA(4)
    Tr3aYuAu = 3._dp * TraceaYA(4)
    Tr3aAuYu = Conjg( Tr3aYuAu )

    diagonal(1,2) = 2.4_dp * gauge2(1) * (S2 + sig(1) )   &
                & + 112.32_dp * gauge2(1) * AbsGM2(1)
    hd(1) = 6._dp*gauge2(2)-1.2_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 4.8_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 2.4_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMe2 = -2._dp * ( MeYeaYeYeaYe + YeaYeYeaYeMe )                       &
      & - 4._dp * ( YeMlaYeYeaYe + YeaYeMeYeaYe + YeaYeYeMlaYe )             &
      & - 8._dp * mH(1) * YeaYeYeaYe                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MeYeaYe + YeaYeMe + 4._dp * MH(1) * YeaYe &
      &                           + 2._dp * (YeMlaYe + AeaAE)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YeaYe + AeaAeYeaYe + YeaYEAeaAe     &
      &           + AeaYeYeaAe + YeaAeAeaYe + Tr3aAdAdaAeAe * YeaYe          &
      &           + Tr3aAdYdaAeYe * AeaYe + Tr3aYdAdaYeAe * YeaAe )          &
      & + hd(1) * ( MeYeaYe + YeaYeMe                                        &
      &           + 2._dp * ( mH(1) * YeaYe + YeMlaYe + AeaAe ) )            &
      & + hd(2) * YeaYe + hc(1) * AeaYe + hc(2) * YeaAe
    Do i1=1,3
     betaMe2(i1,i1) = betaMe2(i1,i1) + diagonal(1,2)
    End Do

    diagonal(3,2) = gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)        &
     &     + 3._dp * gauge2(2) * sig(2)                                 &
     &     + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)      &
     &                   + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2       &
     &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )     &
     &     + 24.84_dp * gauge2(1) * AbsGM2(1)
    betaMl2 = -2._dp * ( MlaYeYeaYeYe + aYeYeaYeYeMl )                   &
      & - 4._dp * ( aYeMeYeaYeYe + aYeYeMlaYeYe + aYeYeaYeMeYe )         &
      & - 8._dp * mH(1) * aYeYeaYeYe                                     &
      & - Tr3aYdYdaYeYe * ( MlaYeYe + aYeYeML + 4._dp * MH(1) * aYeYe    &
      &                   + 2._dp * (aYeMeYe + aAeAE)  )                 &
      & - 4._dp * (aAeAeaYeYe + aYeYeaAeAe + aAeYeaYeAe + aYeAeaAeYe )   &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYeYe +  Tr3aAdAdaAeAe * aYeYe  &
      &           + Tr3aAdYdaAeYe * aYeAe + Tr3aYdAdaYeAe * aAeYe     )  &
      & + 1.2_dp*gauge2(1) * ( MlaYeYe + aYeYeMl                         &
      &                      + 2._dp * ( mH(1) * aYeYe + aYeMeYe + aAeAe &
      &                                - Mhlf(1) * aAeYe                 &
      &                                - Conjg(Mhlf(1)) * aYeAe ) )      &
      & + 4.8_dp * AbsGM2(1) * aYeYe
    Do i1=1,3
     betaMl2(i1,i1) = betaMl2(i1,i1) + diagonal(3,2)
    End Do

    diagonal(4,2) = 0.8_dp * gauge2(1) * (S2 + sig(1)/3._dp )             &
      &   + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp   &
      &   + 8.08e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                    &
      &   + 1.28e2_dp * ( gauge2(1) * AbsGM2(3)                           &
      &                 + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)     &
      &                   * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)+0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) + 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) - 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMd2 = -2._dp * ( MdYdaYdYdaYd + YdaYdYdaYdMd )                       &
      & - 4._dp * ( YdMqaYdYdaYd + YdaYdMdYdaYd + YdaYdYdMqaYd )             &
      & - 8._dp * mH(1) * YdaYdYdaYd                                         &
      & - 2._dp * Tr3aYdYdaYeYe * ( MdYdaYd + YdaYdMd + 4._dp * MH(1) * YdaYd &
      &                           + 2._dp * (YdMqaYd + AdaAD)  )             &
      & - 4._dp * ( Tr3MqaYdYd3aYDMdYd * YdaYd + AdaAdYdaYd + YdaYDAdaAd     &
      &           + AdaYdYdaAd + YdaAdAdaYd + Tr3aAdAdaAeAe * YdaYd          &
      &           + Tr3aAdYdaAeYe * AdaYd + Tr3aYdAdaYeAe * YdaAd )          &
      & + hd(1) * ( MdYdaYd + YdaYdMd                                        &
      &           + 2._dp * ( mH(1) * YdaYd + YdMqaYd + AdaAd ) )            &
      & + hd(2) * YdaYd + hc(1) * AdaYd + hc(2) * YdaAd                      &
      & - 2._dp * (MdYdaYuYuaYd + YdaYuYuaYdMd )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YdaYuYuaYd + YdMqaYuYuaYd+ YdaYuYuMqaYd  &
      &           + YdaYuMuYuaYd + AdaAuYuaYd + YdaYuAuaAd                   &
      &           + AdaYuYuaAd + YdaAuAuaYd )
    Do i1=1,3
     betaMd2(i1,i1) = betaMd2(i1,i1) + diagonal(4,2)
    End Do

    diagonal(5,2) = 0.2_dp * gauge2(1) * (2._dp * S2 + sig(1)/3._dp )     &
      &  + 3._dp * gauge2(2) * sig(2)                                     &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp    &
      &  + 1.99e2_dp * gauge2(1) * AbsGM2(1) / 75._dp                     &
      &  + 32._dp * ( gauge2(1) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)         &
      &              * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp      &
      &  + 33._dp * gauge2(2) * AbsGM2(2)                                 &
      &  + 32._dp * ( gauge2(2) * AbsGM2(3)                               &
      &             + gauge2(3) * AbsGM2(2)  + gauge2(3)*gauge2(2)        &
      &              * Real( Mhlf(2) * Conjg(Mhlf(3)),dp ) )                 &
      &  + 0.4_dp * ( gauge2(2) * AbsGM2(1)                               &
      &             + gauge2(1) * AbsGM2(2) + gauge2(1)*gauge2(2)         &
      &              * Real( Mhlf(2) * Conjg(Mhlf(1)),dp ) )
    hd(1) = 1.6_dp * AbsGM2(1)
    hc(1) = - 0.8_dp * g2Mi(1)
    hc(2) = Conjg( hc(1) )
    hd(2) = 3.2_dp * AbsGM2(1)
    hc(3) = - 1.6_dp * g2Mi(1)
    hc(4) = Conjg( hc(3) )
    betaMq2 = -2._dp * ( MqaYdYdaYdYd + aYdYdaYdYDMq )                       &
      & - 4._dp * ( aYdMdYdaYdYd + aYdYdMqaYdYd + aYdYdaYdMdYd )             &
      & - 8._dp * mH(1) * aYdYdaYdYd                                         &
      & - Tr3aYdYdaYeYe * ( MqaYdYd + aYdYdMq + 4._dp * MH(1) * aYdYd        &
      &                   + 2._dp * (aYdMdYd + aAdAd)  )                     &
      & - 2._dp * ( Tr3MqaYdYd3aYDMdYd * aYdYd                               &
      &           + 2._dp * ( aAdAdaYdYd + aYdYDaAdAd + aAdYdaYdAd           &
      &                     + aYdAdaAdYd )                                   &
      &           + Tr3aAdAdaAeAe * aYdYd + Tr3aAdYdaAeYe * aYdAd            &
      &           + Tr3aYdAdaYeAe * aAdYd )                                  &
      & + 0.4_dp*gauge2(1) * ( MqaYdYd + aYdYdMq                             &
      &                      + 2._dp * ( mH(1) * aYdYd + aYdMdYd + aAdAd ) ) &
      & + hd(1) * aYdYd + hc(1) * aAdYd + hc(2) * aYdAd                      &
      & - 2._dp * ( MqaYuYuaYuYu + aYuYuaYuYuMq )                            &
      & - 4._dp * ( aYuMuYuaYuYu + aYuYuMqaYuYu + aYuYuaYuMuYu )             &
      & - 8._dp * MH(2) * aYuYuaYuYu                                         &
      & - 3._dp * TraceY(4) * ( MqaYuYu + aYuYuMq + 4._dp * MH(2) * aYuYu    &
      &                       + 2._dp * (aYuMuYu + aAuAu) )                  &
      & - 2._dp * ( Tr3MqaYuYu3aYuMu * aYuYu                                 &
      &           + 2._dp * ( aAuAuaYuYu + aYuYuaAuAu                        &
      &                     + aAuYuaYuAu + aYuAuaAuYu )                      &
      &           + Tr3aAuAu * aYuYu + Tr3aAuYu * aYuAu + Tr3aYuAu * aAuYu ) &
      & + 0.8_dp*gauge2(1) * ( MqaYuYu + aYuYuMq                             &
      &                      + 2._dp * ( MH(2) * aYuYu + aYuMuYu + aAuAu ) ) &
      & + hd(2) * aYuYu + hc(3) * aAuYu + hc(4) * aYuAu
    Do i1=1,3
     betaMq2(i1,i1) = betaMq2(i1,i1) + diagonal(5,2)
    End Do

    diagonal(6,2) = 1.6_dp * gauge2(1) * (2._dp*sig(1)/3._dp - S2)           &
      &  + 16._dp * gauge2(3) * (-8._dp * AbsGM2(3) + sig(3) ) / 3._dp       &
      &  + 3424._dp * gauge2(1) * AbsGM2(1) / 75._dp                         &
      &  + 512._dp * ( gauge2(1) * AbsGM2(3)                                 &
      &              + gauge2(3) * AbsGM2(1) + gauge2(1)*gauge2(3)           &
      &                * Real( Mhlf(1) * Conjg(Mhlf(3)),dp ) ) / 4.5e1_dp

    hd(1) = 6._dp*gauge2(2)-0.4_dp*gauge2(1)
    hd(2) = 24._dp * AbsGM2(2) - 1.6_dp * AbsGM2(1)
    hc(2) = - 12._dp * g2Mi(2) + 0.8_dp * g2Mi(1)
    hc(1) = Conjg( hc(2) )
    betaMu2 = -2._dp * ( MuYuaYuYuaYu + YuaYuYuaYuMu )                       &
      & - 4._dp * ( YuMqaYuYuaYu + YuaYuMuYuaYu + YuaYuYuMqaYu )             &
      & - 8._dp * MH(2) * YuaYuYuaYu                                         &
      & - 6._dp * TraceY(4) * ( MuYuaYu + YuaYuMu + 4._dp * MH(2) * YuaYu    &
      &                       + 2._dp * (YuMqaYu + AuaAu)  )                 &
      & - 4._dp * ( Tr3MqaYuYu3aYuMu * YuaYu + AuaAuYuaYu + YuaYuAuaAu       &
      &           + AuaYuYuaAu + YuaAuAuaYu + Tr3aAuAu * YuaYu               &
      &           + Tr3aAuYu * AuaYu + Tr3aYuAu * YuaAu )                    &
      & + hd(1) * ( MuYuaYu + YuaYuMu                                        &
      &           + 2._dp * ( MH(2) * YuaYu + YuMqaYu + AuaAu ) )            &
      & + hd(2) * YuaYu + hc(1) * AuaYu + hc(2) * YuaAu                      &
      & - 2._dp * (MuYuaYdYdaYu + YuaYdYdaYuMu )                             &
      & - 4._dp * ( (MH(1)+MH(2)) * YuaYdYdaYu + YuMqaYdYdaYu + YuaYdYdMqaYu &
      &           + YuaYdMdYdaYu + AuaAdYdaYu + YuaYdAdaAu                   &
      &           + AuaYdYdaAu + YuaAdAdaYu )
    Do i1=1,3
     betaMu2(i1,i1) = betaMu2(i1,i1) + diagonal(6,2)
    End Do

   End If 

  !------------------------------------------
  ! beta functions for Higgs mass parameters
  !------------------------------------------
   traceMH1(1) = mH(1) * (TraceY(1) + 6._dp * lam12)            &
             & + Real( cTrace(YeMlaYe),dp ) &
             & + Real( cTrace(aYeMeYe),dp ) + TraceA(1)
   traceMH1(2) = mH(1) * TraceY(3) + Real( cTrace(YdMqaYd),dp )   &
             & + Real( cTrace(aYdMdYd),dp ) + TraceA(3) + Alam12  &
             & + lam12 * MT(1)
   betamH11 = 6._dp * TraceMH1(2) + 2._dp * TraceMH1(1)      &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) - 0.6_dp * S1

   traceMH2(1) = 3._dp * (2._dp* mH(2)+ MT(2)) * lam22 + 3._dp * Alam22
   traceMH2(2) = mH(2) * TraceY(4) + Real( cTrace(YuMqaYu),dp )  &
             & + Real( cTrace(aYuMuYu),dp ) + TraceA(4)
   betamH21 = 2._dp * TraceMH2(1) + 6._dp * TraceMH2(2)       &
          & - 6._dp * AbsGM2(2) - 1.2_dp * AbsGM2(1) + 0.6_dp * S1

   betaMT2(1) = MT(1) * (lam12 + TraceY(2)) + 2._dp * mH(1) * lam12   &
            & + 2._dp * Real( cTrace(aYTMlYT),dp ) + TraceA(2) + Alam12 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

   betaMT2(2) = MT(2) * lam22 + 2._dp * mH(2) * lam22 + Alam22 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

   betaMT2 = 2._dp * betaMT2
 
   betaMZ2(2) = - (0.4_dp * AbsGM2(1) + 32._dp * AbsGM2(3) ) / 3._dp  &
              & - 6._dp  * AbsGM2(2)

   betaMZ2(1) = 2._dp * MZ(1) * TraceY(5) + 2._dp * TraceA(5) + betaMZ2(2) &
            & + 2._dp * Real( cTrace(aYZMdYZ) + cTrace(YZMlaYZ) ,dp ) 

   betaMS2(2) = - (3.2_dp * AbsGM2(1) + 40._dp * AbsGM2(3) ) / 3._dp

   betaMS2(1) = MS(1) * TraceY(6) + 2._dp * Real( cTrace(aYSMdYS),dp ) &
            & + TraceA(6) + betaMS2(2)

   betaMS2 = 2._dp * betaMS2


   If (TwoLoopRGE) Then
    traceMH1(3) = MH(1) * (6._dp*TraceY2(2) + 2._dp*TraceY2(1) + TraceY2(4) ) &
              & + MH(2) * TraceY2(4)                                          &
              & + 6._dp * ( Real( cTrace(MqaYdYdaYdYd),dp )                   &
              &           + Real( cTrace(aYdMdYdaYdYd),dp )                   &
              &           + Real( cTrace(aAdAdaYdYd),dp )                     &
              &           + Real( cTrace(aAdYdaYdAd),dp )  )                  &
              & + 2._dp * ( Real( cTrace(MlaYeYeaYeYe),dp )                   &
              &           + Real( cTrace(aYeMeYeaYeYe),dp )                   &
              &           + Real( cTrace(aAeAeaYeYe),dp )                     &
              &           + Real( cTrace(aAeYeaYeAe),dp )  )                  &
       & + Real( cTrace(YdMqaYuYuaYd),dp ) + Real( cTrace(YdaYuMuYuaYd),dp ) &
       & + Real( cTrace(YdaYuYuMqaYd),dp ) + Real( cTrace(YuaYdMdYdaYu),dp ) &
       & + Real( cTrace(YdaAuAuaYd),dp ) + Real( cTrace(AdaYuYuaAd),dp )     &
              & + Real( cTrace(AdaAuYuaYd),dp ) + Real( cTrace(YdaYuAuaAd),dp )
    betaMH12 = - 6._dp * traceMH1(3)                                       &
      &   + (32._dp*gauge2(3) - 0.8_dp*gauge2(1) ) * traceMH1(2)           &
      &   + 64._dp * ( AbsGM2(3) * TraceY(3)                               &
      &             - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(3),dp ) )     &
      &   - 1.6_dp * ( AbsGM2(1) * TraceY(3)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(3),dp ) )     &
      &   + 2.4_dp*gauge2(1) * traceMH1(1)                                 &
      &   + 4.8_dp * ( AbsGM2(1) * TraceY(1)                               &
      &             - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(1),dp ) )     &
      &   + gauge2(1) * ( 0.6_dp * sig(1) - 1.2_dp * S2)                   &
      &   + 3._dp * gauge2(2) * sig(2)                                     &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)          &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2           &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

    traceMH2(2) = MH(2) * ( 6._dp * TraceY2(3) + TraceY2(4) )                 &
      &  + MH(1) * TraceY2(4)                                                 &
 &  + 6._dp * ( Real(cTrace(MqaYuYuaYuYu),dp) + Real(cTrace(aYuMuYuaYuYu),dp) &
 &            + Real(cTrace(aAuAuaYuYu),dp) + Real(cTrace(aAuYuaYuAu),dp)  )  &
 &  + Real( cTrace(YuMqaYdYdaYu),dp ) + Real( cTrace(YuaYdMdYdaYu),dp )       &
 &  + Real( cTrace(YuaYdYdMqaYu),dp ) + Real( cTrace(YdaYuMuYuaYd),dp )       &
 &  + Real( cTrace(YuaAdAdaYu),dp ) + Real( cTrace(AuaYdYdaAu),dp )           &
      &  + Real( cTrace(AuaAdYdaYu),dp ) + Real( cTrace(YuaYdAdaAu),dp )
    betaMH22 = - 6._dp * traceMH2(2)                                          &
      &   + (32._dp*gauge2(3) + 1.6_dp*gauge2(1) ) * traceMH2(1)              &
      &   + 64._dp * ( AbsGM2(3) * TraceY(4)                                  &
      &          - gauge2(3) * Real( Conjg(Mhlf(3))*TraceaYA(4),dp ) )        &
      &   + 3.2_dp * ( AbsGM2(1) * TraceY(4)                                  &
      &          - gauge2(1) * Real( Conjg(Mhlf(1))*TraceaYA(4),dp ) )        &
      &   + gauge2(1) * ( 0.6_dp * sig(1) + 1.2_dp * S2)                      &
      &   + 3._dp * gauge2(2) * sig(2)                                        &
      &   + gauge2(2) * ( 33._dp * AbsGM2(2) + 3.6_dp * AbsGM2(1)             &
      &                 + 3.6_dp * gauge2(1) * ( Abs(Mhlf(2))**2              &
      &                          +Real( Mhlf(1)*Conjg(Mhlf(2)),dp ) ) )       &
      &   + 24.84_dp * gauge2(1) * AbsGM2(1)

   End If
!------------------------------
! beta function for the 15-plet
!------------------------------
   betaMT15 = MT15 * ( TraceY(2) + lam12 + lam22 &
            &        - 2.4_dp * gauge2(1) - 8._dp * gauge2(2)  )
   betaMZ15 = MZ15 * ( TraceY(5) - gauge2(1)/15._dp - 3._dp * gauge2(2) &
                     & - 16._dp * gauge2(3) / 3._dp  )
   betaMS15 = MS15 * ( TraceY(6)  &
                     & - (3.2_dp * gauge2(1) + 40._dp * gauge2(3))/3._dp  )
!-----------------------------
! beta functions for mu and B
!-----------------------------
   TraceMue(1) = 3._dp * (TraceY(3)+TraceY(4)+lam12+lam22) + TraceY(1)  &
             & - 3._dp * gauge2(2) - 0.6_dp * gauge2(1)
   betaMue1 = mue * TraceMue(1)

   TraceB(1) = 6._dp * (TraceaYA(3)+TraceaYA(4) + lam1Alam1 + lam2Alam2) &
           & + 2._dp * TraceaYA(1)                         &
           & + 6._dp * g2Mi(2) + 1.2_dp * g2Mi(1)
   betaB1 = mue * TraceB(1) + B * TraceMue(1)

   If (TwoLoopRGE) Then
    TraceMue(2) = - 3._dp * ( 3._dp * (TraceY2(2) + TraceY2(3) )          &
      &                     + 2._dp * TraceY2(4) + TraceY2(1) )           &
      &         + (16._dp * gauge2(3) + 0.8_dp * gauge2(1)) * TraceY(4)   &
      &         + (16._dp * gauge2(3) - 0.4_dp * gauge2(1)) * TraceY(3)   &
      &         + 1.2_dp * gauge2(1) * TraceY(1)                          &
      &         + 7.5_dp * gauge2(2)**2                                   &
      &         + 1.8_dp * gauge2(2) * gauge2(1)                          &
      &         + 4.14_dp * gauge2(1)**2
    betaMue2 = mue * TraceMue(2)

    TraceB(1) = cTrace( 3._dp * ( Matmul2(AuaYu,YuaYu,OnlyDiagonal)     &
              &                 + Matmul2(AdaYd,YdaYd,OnlyDiagonal) )   &
              &       + MatMul2(AeaYe,YeaYe,OnlyDiagonal)               &
              &       + Matmul2(aYuAu,aYdYd,OnlyDiagonal)               &
              &       + MatMul2(aYdAd,aYuYu,OnlyDiagonal) ) 
    TraceB(2) = -12._dp * TraceB(1)                                           &
      &   + (32._dp * gauge2(3) + 1.6_dp * gauge2(1) ) * TraceaYA(4)          &
      &   + (32._dp * gauge2(3) - 0.8_dp * gauge2(1) ) * TraceaYA(3)          &
      &   + 2.4_dp * gauge2(1) * TraceaYA(1)                                  &
      &   - ( 32._dp * g2Mi(3)                                    &
      &     + 1.6_dp * g2Mi(1) ) * TraceY(4)                      &
      &   - ( 32._dp * g2Mi(3) - 0.8_dp * g2Mi(1) ) * TraceY(3)   &
      &   - 2.4_dp * g2Mi(1) * TraceY(1)                          &
      &   - 30._dp * gauge2(2)**2 * Mhlf(2)                                   &
      &   - 3.6_dp * gauge2(2) * gauge2(1) * (Mhlf(1) + Mhlf(2) )             &
      &   - 16.56_dp * gauge2(1)**2 * Mhlf(1) 
    betaB2 = mue * TraceB(2) + B * TraceMue(2)

   End If

  !--------------------------------
  ! neutrino dim. 5 operator
  !--------------------------------
  diagonal(5,1) = 6._dp * TraceY(4) - 2._dp * gauge2(1) - 6._dp * gauge2(2)
  betaMnu1 = Matmul( Transpose(aYeYe), Mnu) + Matmul(Mnu, aYeYe)  &
          & + diagonal(5,1) * Mnu
  

 !---------------
 ! 2-loop RGEs
 !---------------
  b_1a = b_1 + Delta_b_1
  b_2a = b_2 + Delta_b_2

  If (TwoLoopRGE) Then 
 !----------------------
 ! gauge couplings
 !----------------------
   If (MaxVal(Delta_b_2).eq.0._dp) then
    Dgauge = oo16pi2 * gauge * gauge2  &
         & * ( b_1a + oo16pi2 * ( Matmul(b_2a,gauge2) - a_2(:,1) * TraceY(1) &
         &                      - Matmul(a_2(:,2:3),TraceY(3:4))  ))
   else
    Dgauge = oo16pi2 * gauge * gauge2  &
         & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2b(:,1:6),TraceY) )  &
         &                      - a_2b(:,7)*lam12 - a_2b(:,8)*lam22 ) 
   End if
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * ( betaYe1 + oo16pi2 * betaYe2 )
   DYT = oo16pi2 * betaYT1 ! + oo16pi2 * betaYnu2 )
   DYd = oo16pi2 * ( betaYd1 + oo16pi2 * betaYd2 )
   DYu = oo16pi2 * ( betaYu1 + oo16pi2 * betaYu2 )
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   Do i1 = 1,3    
    sumI = 0._dp
    Do i2=1,3
     sumI = sumI + b_2a(i1,i2) * gauge2(i2) * (Mhlf(i1) + Mhlf(i2) ) 
    End Do
    Do i2=1,6
     sumI = sumI + a_2b(i1,i2) * ( TraceaYA(i2) - Mhlf(i1)*TraceY(i2) )
    End Do
    sumI = sumI + a_2b(i1,7) * lam1Alam1 + a_2b(i1,8) * lam2Alam2
    DMhlf(i1) = oo8pi2 * gauge2(i1) * ( b_1a(i1) * Mhlf(i1) + oo16pi2 * sumI)
   End Do
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * ( betaAe1 + oo16pi2 * betaAe2 )
   DAd = oo16pi2 * ( betaAd1 + oo16pi2 * betaAd2 )
   DAu = oo16pi2 * ( betaAu1 + oo16pi2 * betaAu2 )
   DAT = oo16pi2 * betaAT1
   DAZ = oo16pi2 * betaAZ1
   DAS = oo16pi2 * betaAS1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * ( betaMe1 + oo16pi2 * betaMe2 )
   DMl = oo16pi2 * ( betaMl1 + oo16pi2 * betaMl2 )
   DMd = oo16pi2 * ( betaMd1 + oo16pi2 * betaMd2 )
   DMq = oo16pi2 * ( betaMq1 + oo16pi2 * betaMq2 )
   DMu = oo16pi2 * ( betaMu1 + oo16pi2 * betaMu2 )
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * ( betaMH11 + oo16pi2 * betaMH12 )
   DmH(2) = oo16pi2 * ( betaMH21 + oo16pi2 * betaMH22 )
   DmT2 = oo16pi2 * betaMT2
   DmZ2 = oo16pi2 * betaMZ2
   DmS2 = oo16pi2 * betaMS2
  !----------
  ! 15-plet
  !----------
  DMT15 = oo16pi2 * betaMT15
  DMZ15 = oo16pi2 * betaMZ15
  DMS15 = oo16pi2 * betaMS15
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * ( betaMue1 + oo16pi2 * betaMue2 )
   DB = oo16pi2 * ( betaB1 + oo16pi2 * betaB2 )

 !----------------------
 ! left neutrino mass
 !----------------------
   DMnu = oo16pi2 * betaMnu1
 !---------------
 ! 1-loop RGEs
 !---------------
  Else 
 !----------------------
 ! gauge couplings
 !----------------------
   Dgauge = oo16pi2 * gauge * gauge2 * b_1a 
 !--------------------
 ! Yukawa couplings
 !--------------------
   DYe = oo16pi2 * betaYe1
   DYT = oo16pi2 * betaYT1 
   DYd = oo16pi2 * betaYd1
   DYu = oo16pi2 * betaYu1
   DYZ = oo16pi2 * betaYZ1 
   DYS = oo16pi2 * betaYS1 
   Dlam1 = oo16pi2 * betalam11
   Dlam2 = oo16pi2 * betalam21
 !--------------------------
 ! gaugino mass parameters
 !--------------------------
   DMhlf = oo8pi2 * gauge2 * b_1a * Mhlf
  !--------------------------
  ! trilinear parameters
  !--------------------------
   DAe = oo16pi2 * betaAe1
   DAd = oo16pi2 * betaAd1
   DAu = oo16pi2 * betaAu1
   DAT = oo16pi2 * betaAT1
   DAZ = oo16pi2 * betaAZ1
   DAS = oo16pi2 * betaAS1
   DAlam1 = oo16pi2 * betaAlam11
   DAlam2 = oo16pi2 * betaAlam21
  !---------------------------
  ! Sfermion mass parameters
  !---------------------------
   DMe = oo16pi2 * betaMe1
   DMl = oo16pi2 * betaMl1
   DMd = oo16pi2 * betaMd1
   DMq = oo16pi2 * betaMq1
   DMu = oo16pi2 * betaMu1
  !-----------------------
  ! Higgs mass parameters
  !-----------------------
   DmH(1) = oo16pi2 * betaMH11
   DmH(2) = oo16pi2 * betaMH21
   DmT2 = oo16pi2 * betaMT2
   DmZ2 = oo16pi2 * betaMZ2
   DmS2 = oo16pi2 * betaMS2
  !----------
  ! 15-plet
  !----------
  DMT15 = oo16pi2 * betaMT15
  DMZ15 = oo16pi2 * betaMZ15
  DMS15 = oo16pi2 * betaMS15
  !----------
  ! mu and B
  !----------
   DMue = oo16pi2 * betaMue1
   DB = oo16pi2 * betaB1
  !----------------------
  ! left neutrino mass
  !----------------------
   DMnu = oo16pi2 * betaMnu1
  End If


  !---------------------------------------
  ! to avoid numerical problems in odeint
  !---------------------------------------
  Do i1=1,3
   DMe(i1,i1) = Real(DMe(i1,i1),dp)
   DMl(i1,i1) = Real(DMl(i1,i1),dp)
   DMd(i1,i1) = Real(DMd(i1,i1),dp)
   DMu(i1,i1) = Real(DMu(i1,i1),dp)
   DMq(i1,i1) = Real(DMq(i1,i1),dp)
  End Do

  Call ParametersToG5(Dgauge, DYe, DYT, DYd, DYu, DYZ, DYS, Dlam1, Dlam2    &
          & , DMhlf, DAe, DAT, DAd, DAu, DAZ, DAS, DAlam1, DAlam2, DMe, DMl &
          & , DMd, DMq, DMu, DMh, DMT2, DmZ2, DmS2, DMT15, DMZ15, DMS15     &
          & , Dmue, DB, DMnu, f)

  Iname = Iname - 1

 End Subroutine rge356


 Subroutine RGE_SU5(len,T,GY,F)
 !-----------------------------------------------------------------------
 ! Right hand side of renormalization group equations dGY_i/dT = F_i(G) 
 ! of the gauge and Yukawa couplings, soft SUSY breaking parameters
 ! For the determination of M_GUT and the value of alpha_GUT
 ! and values of the Yukawas, all complex 3 times 3 matrices
 ! written by Werner Porod, 17.8.1999
 ! 25.09.01: portation to f90
 !-----------------------------------------------------------------------
 Implicit None

  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2, SumI
  Real(dp) :: g5, g52, b5, lam2, lamp2, GammaH, GammaHbar, GammaS, TraceY(3) &
    & , Q
  Complex(dp) :: lam, lamp
  Complex(dp), Dimension(3,3) :: Y_u, Y_d, Y_N, YdaYd, aYdYd, YnaYn, aYnYn &
    & , YuaYu, aYuYu, GammaT, GammaF, GammaN, aY_u, aY_d, aY_N 
  Real(dp) :: dG5
  Complex(dp) :: Dlam, Dlamp
  Complex(dp), Dimension(3,3) :: DYu1, DYn1, DYd1

  Complex(dp) :: MuH, DMuH, MuSig, DMuSig, M5, DM5, gM5, GammaAH, GammaAHbar &
    & , GammaAS
  Complex(dp), Dimension(3,3) :: A_u, A_d, A_N, AdaYd, aYdAd, aYnAn &
    & , AuaYu, aYuAu, GammaAT, GammaAF, GammaAN
  Complex(dp) :: Alam, Alamp, Alamlam, Alamplamp, DAlam, DAlamp, TraceAY(3)
  Complex(dp), Dimension(3,3) :: DAu1, DAn1, DAd1

  Complex(dp), Dimension(3,3) :: MT2, MF2, MN2, DMT2, DMF2, DMN2    &
    & , aA_u, aA_d, aA_N, AdaAd, aAdAd, AnaAn, aAnAn, AuaAu, aAuAu  &
    & , MT2YuaYu, YuaYuMT2, YuMT2aYu, MT2YdaYd, YdaYdMT2, aYdMT2Yd  &
    & , YdMF2aYd, MF2aYdYd, aYdYdMF2, MF2aYnYn, aYnYnMF2, YnMF2aYn  &
    & , aYnYnMN2, MN2aYnYn, aYnMN2Yn, aYnMF2Yn

  Complex(dp) :: TraceA(3)
  Real(dp) :: M52, Alam2, Alamp2, MH2, MHbar2, MS2, DMH2, DMHbar2, DMS2
 
  Iname = Iname + 1
  NameOfUnit(Iname) = 'RGE_SU5'

  q = t

  OnlyDiagonal = .Not.GenerationMixing

  !----------------------------------
  ! gauge couplings beta functions
  !----------------------------------
  g5 = gy(1)
  g52 = g5**2
  b5 = -3._dp

  !----------------------------
  ! Yukawa couplings
  !----------------------------
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    Y_u(i1,i2) = Cmplx( gy(SumI-6), gy(SumI-5),dp )
    Y_d(i1,i2) = Cmplx( gy(SumI+12), gy(SumI+13),dp )
    Y_n(i1,i2) = Cmplx( gy(SumI+30), gy(SumI+31),dp )
   End Do
  End Do

  lam = Cmplx( gy(56), gy(57),dp )
  lamp = Cmplx( gy(58), gy(59),dp )
  lam2 = Abs(lam)**2
  lamp2 = Abs(lamp)**2
  !-----------------
  ! beta functions
  !-----------------
  Call Adjungate(Y_d,aY_d)
  Call Adjungate(Y_n,aY_n)
  Call Adjungate(Y_u,aY_u)

  aYdYd = Matmul2(aY_d, Y_d,OnlyDiagonal)
  aYnYn = Matmul2(aY_n, Y_n,OnlyDiagonal)
  aYuYu = Matmul2(aY_u, Y_u,OnlyDiagonal)
  YdaYd = Matmul2(Y_d, aY_d,OnlyDiagonal)
  YnaYn = Matmul2(Y_n, aY_n,OnlyDiagonal)
  YuaYu = Matmul2(Y_u, aY_u,OnlyDiagonal)
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aYdYd(i1,i1) = Real(aYdYd(i1,i1), dp)
   aYnYn(i1,i1) = Real(aYnYn(i1,i1), dp)
   aYuYu(i1,i1) = Real(aYuYu(i1,i1), dp)
   YdaYd(i1,i1) = Real(YdaYd(i1,i1), dp)
   YnaYn(i1,i1) = Real(YnaYn(i1,i1), dp)
   YuaYu(i1,i1) = Real(YuaYu(i1,i1), dp)
  End Do

  TraceY(1) = Real( cTrace(aYnYn),dp )
  TraceY(2) = Real( cTrace(aYdYd),dp )
  TraceY(3) = Real( cTrace(aYuYu),dp )

  GammaT = 2._dp * YdaYd + 3._dp * YuaYu
  GammaF = 4._dp * aYdYd + aYnYn
  GammaN = 5._dp * aYnYn

  Do i1=1,3
   GammaT(i1,i1) = GammaT(i1,i1) - 7.2_dp * g52
   GammaF(i1,i1) = GammaF(i1,i1) - 4.8_dp * g52
  End Do

  GammaHbar = 4._dp * TraceY(2) + 4.8_dp * (lam2 - g52)
  GammaH = 3._dp * TraceY(3) + TraceY(1) + 4.8_dp * (lam2 - g52)
  GammaS = 1.05_dp * lamp2 + lam2 - 10._dp * g52

  Dg5 = oo16pi2 * b5 * g5 * g52
 
  DYu1 = Y_u * GammaH + Matmul2(GammaT,Y_u,OnlyDiagonal) &
       & + MatMul2(Y_u,Transpose(GammaT),OnlyDiagonal)
  DYd1 = Y_d * GammaHbar + Matmul2(GammaT,Y_d,OnlyDiagonal) &
       & + MatMul2(Y_d,GammaF,OnlyDiagonal)
  DYn1 = Y_n * GammaH + Matmul2(GammaT,Y_n,OnlyDiagonal) &
       & + MatMul2(Y_n,GammaN,OnlyDiagonal)

  DYn1 = oo16pi2 * DYn1 
  DYd1 = oo16pi2 * DYd1 
  DYu1 = oo16pi2 * DYu1 

  Dlam = oo16pi2 * lam * (GammaH + GammaHbar + GammaS)
  Dlamp = oo16pi2 * lamp * 3._dp * GammaS

  f(1) = Dg5
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    f(SumI-6) = Real(DYu1(i1,i2),dp)
    f(SumI-5) = Aimag(DYu1(i1,i2))
    f(SumI+12) = Real(DYd1(i1,i2),dp)
    f(SumI+13) = Aimag(DYd1(i1,i2))
    f(SumI+30) = Real(DYn1(i1,i2),dp)
    f(SumI+31) = Aimag(DYn1(i1,i2))
   End Do
  End Do
  f(56) = Real(Dlam, dp)
  f(57) = Aimag(Dlam)
  f(58) = Real(Dlamp, dp)
  f(59) = Aimag(Dlamp)

  If (len.Eq.59) Then
   Iname = Iname - 1
   Return
  End If

  !----------------------------------
  ! bilinear Higgs parameters
  !----------------------------------
  MuH = Cmplx( gy(60), gy(61),dp )
  MuSig = Cmplx( gy(62), gy(63),dp )

  DMuH = oo16pi2 * MuH * (GammaH + GammaHbar)
  DMuSig = oo16pi2 * MuSig * 2._dp * GammaS

  f(60) = Real(DMuH, dp)
  f(61) = Aimag(DMuH)
  f(62) = Real(DMuSig, dp)
  f(63) = Aimag(DMuSig)

  If (len.Eq.63) Then
   Iname = Iname - 1
   Return
  End If

  !----------------------------------
  ! gaugino beta function
  !----------------------------------
  M5 = Cmplx( gy(64), gy(65) , dp)
  gM5 = g52 * M5

  !----------------------------
  ! A-parameters
  !----------------------------
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    A_u(i1,i2) = Cmplx( gy(SumI+58), gy(SumI+59),dp )
    A_d(i1,i2) = Cmplx( gy(SumI+76), gy(SumI+77),dp )
    A_n(i1,i2) = Cmplx( gy(SumI+94), gy(SumI+95),dp )
   End Do
  End Do

  Alam = Cmplx( gy(120), gy(121),dp )
  Alamp = Cmplx( gy(122), gy(123),dp )
  Alamlam = Alam * Conjg(lam)
  Alamplamp = Alamp * Conjg(lamp)
  !-----------------
  ! beta functions
  !-----------------

  aYdAd = Matmul2(aY_d, A_d,OnlyDiagonal)
  aYnAn = Matmul2(aY_n, A_n,OnlyDiagonal)
  aYuAu = Matmul2(aY_u, A_u,OnlyDiagonal)
  AdaYd = Matmul2(A_d, aY_d,OnlyDiagonal)
  AuaYu = Matmul2(A_u, aY_u,OnlyDiagonal)

  TraceAY(1) = cTrace(aYnAn)
  TraceAY(2) = cTrace(aYdAd)
  TraceAY(3) = cTrace(aYuAu)

  GammaAT = -2._dp * AdaYd - 3._dp * AuaYu
  GammaAF = -4._dp * aYdAd - aYnAn
  GammaAN = -5._dp * aYnAn

  Do i1=1,3
   GammaAT(i1,i1) = GammaAT(i1,i1) - 7.2_dp * gM5
   GammaAF(i1,i1) = GammaAF(i1,i1) - 4.8_dp * gM5
  End Do

  GammaAHbar = - 4._dp * TraceAY(2) - 4.8_dp * (Alamlam + gM5)
  GammaAH = -3._dp * TraceAY(3) - TraceAY(1) - 4.8_dp * (Alamlam + gM5)
  GammaAS = - 1.05_dp * Alamplamp - Alamlam - 10._dp * gM5

  DM5 = oo8pi2 * b5 * gM5
 
  DAu1 = A_u * GammaH + Matmul2(GammaT,A_u,OnlyDiagonal)               &
       & + MatMul2(A_u,Transpose(GammaT),OnlyDiagonal)                 &
       & - 2._dp * ( Y_u * GammaAH + Matmul2(GammaAT,Y_u,OnlyDiagonal) &
       &           + MatMul2(Y_u,Transpose(GammaAT),OnlyDiagonal) ) 
  DAd1 = A_d * GammaHbar + Matmul2(GammaT,A_d,OnlyDiagonal)               &
       & + MatMul2(A_d,GammaF,OnlyDiagonal)                               &
       & - 2._dp * ( Y_d * GammaAHbar + Matmul2(GammaAT,Y_d,OnlyDiagonal) &
       &           + MatMul2(Y_d,GammaAF,OnlyDiagonal) )
  DAn1 = A_n * GammaH + Matmul2(GammaT,A_n,OnlyDiagonal)              &
       & + MatMul2(A_n,GammaN,OnlyDiagonal)                           &
       & - 2._dp * (Y_n * GammaAH + Matmul2(GammaAT,Y_n,OnlyDiagonal) &
       &           + MatMul2(Y_n,GammaAN,OnlyDiagonal) )

  DAn1 = oo16pi2 * DAn1 
  DAd1 = oo16pi2 * DAd1 
  DAu1 = oo16pi2 * DAu1 

  DAlam = oo16pi2 * ( Alam * (GammaH + GammaHbar + GammaS)       &
        &           - 2._dp * lam * (GammaAH + GammaAHbar + GammaAS) )
  DAlamp = oo16pi2 * ( Alamp * 3._dp * GammaS - lamp * 6._dp * GammaAS )

  f(64) = Real(DM5, dp)
  f(65) = Aimag( DM5 )

  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    f(SumI+58) = Real(DAu1(i1,i2),dp)
    f(SumI+59) = Aimag(DAu1(i1,i2))
    f(SumI+76) = Real(DAd1(i1,i2),dp)
    f(SumI+77) = Aimag(DAd1(i1,i2))
    f(SumI+94) = Real(DAn1(i1,i2),dp)
    f(SumI+95) = Aimag(DAn1(i1,i2))
   End Do
  End Do
  f(120) = Real(DAlam, dp)
  f(121) = Aimag(DAlam)
  f(122) = Real(DAlamp, dp)
  f(123) = Aimag(DAlamp)

  !----------------------------
  ! scalar masses squared
  !----------------------------
  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    MT2(i1,i2) = Cmplx( gy(SumI+116), gy(SumI+117),dp )
    MF2(i1,i2) = Cmplx( gy(SumI+134), gy(SumI+135),dp )
    MN2(i1,i2) = Cmplx( gy(SumI+152), gy(SumI+153),dp )
   End Do
  End Do

  MHbar2 = gy(178)
  MH2 = gy(179)
  mS2 = gy(180) 

  Call Adjungate(A_d,aA_d)
  Call Adjungate(A_n,aA_n)
  Call Adjungate(A_u,aA_u)

  aAdAd = Matmul2(aA_d, A_d,OnlyDiagonal)
  aAnAn = Matmul2(aA_n, A_n,OnlyDiagonal)
  aAuAu = Matmul2(aA_u, A_u,OnlyDiagonal)
  AdaAd = Matmul2(A_d, aA_d,OnlyDiagonal)
  AnaAn = Matmul2(A_n, aA_n,OnlyDiagonal)
  AuaAu = Matmul2(A_u, aA_u,OnlyDiagonal)
  
  !------------------------------------------------
  ! these are hermitian matrices, clean up to
  ! avoid numerical problems
  !------------------------------------------------
  Do i1=1,3
   aAdAd(i1,i1) = Real(aAdAd(i1,i1), dp)
   aAnAn(i1,i1) = Real(aAnAn(i1,i1), dp)
   aAuAu(i1,i1) = Real(aAuAu(i1,i1), dp)
   AdaAd(i1,i1) = Real(AdaAd(i1,i1), dp)
   AnaAn(i1,i1) = Real(AnaAn(i1,i1), dp)
   AuaAu(i1,i1) = Real(AuaAu(i1,i1), dp)
  End Do

  Alam2 = Abs(Alam)**2
  Alamp2 = Abs(Alamp)**2

  TraceA(1) = Real( cTrace(aAnAn),dp )
  TraceA(2) = Real( cTrace(aAdAd),dp )
  TraceA(3) = Real( cTrace(aAuAu),dp )

  MT2YuaYu = MatMul2( MT2, YuaYu, OnlyDiagonal )
  YuaYuMT2 = MatMul2( YuaYu, MT2, OnlyDiagonal )
  YuMT2aYu = MatMul3( Y_u, MT2, aY_u, OnlyDiagonal )
  MT2YdaYd = MatMul2( MT2, YdaYd, OnlyDiagonal )
  YdaYdMT2 = MatMul2( YdaYd, MT2, OnlyDiagonal )
  aYdMT2Yd = MatMul3( aY_d, MT2, Y_d, OnlyDiagonal )

  YdMF2aYd = MatMul3( Y_d, MF2, aY_d, OnlyDiagonal )
  MF2aYdYd = MatMul2( MF2, aYdYd, OnlyDiagonal )
  aYdYdMF2 = MatMul2( aYdYd, MF2, OnlyDiagonal )
  MF2aYnYn = MatMul2( MF2, aYnYn, OnlyDiagonal )
  aYnYnMF2 = MatMul2( aYnYn, MF2, OnlyDiagonal )
  YnMF2aYn = MatMul3( Y_n, MF2, aY_n, OnlyDiagonal )
  aYnMF2Yn = MatMul3( aY_n, MF2, Y_n, OnlyDiagonal )

  aYnYnMN2 = MatMul2( aYnYn, MN2, OnlyDiagonal )
  MN2aYnYn = MatMul2( MN2, aYnYn, OnlyDiagonal )
  aYnMN2Yn = MatMul3( aY_n, MN2, Y_n, OnlyDiagonal )

  DMT2 = 4._dp * AdaAd + 6._dp * AuaAu                      &
     & + 3._dp * ( MT2YuaYu + 2._dp * YuMT2aYu + YuaYuMT2)  &
     & + 2._dp * ( MT2YdaYd + 2._dp * YdMF2aYd + YdaYdMT2)  &
     & + 6._dp * MH2 * YuaYu + 4._dp * MHbar2 * YdaYd
  DMF2 = 8._dp * aAdAd + 2._dp * aAnAn                      &
     & + 4._dp * (MF2aYdYd + 2._dp * aYdMT2Yd + aYdYdMF2 )  &
     & + MF2aYnYn + aYnYnMF2 + 2._dp * aYnMN2Yn             &
     & + 2._dp * MHbar2 * aYdYd+ 2._dp * MH2 * aYnYn
  DMN2 = 10._dp * aAnAn                                    &
     & + 5._dp * ( MN2aYnYn + 2._dp * aYnMF2Yn + aYnYnMN2) &
     & + 10._dp * MH2 * aYnYn

  M52 = g52 * Abs(M5)**2

  Do i1=1,3
   DMT2(i1,i1) = DMT2(i1,i1) - 28.8_dp * M52
   DMF2(i1,i1) = DMF2(i1,i1) - 19.2_dp * M52
  End Do

  DMH2 = 6._dp * TraceA(3) + 2._dp * TraceA(1)  &
     & + 9.6_dp * (Alam2 - 2._dp * M52            &
     &            + lam2 * (MH2 + MHbar2 + MS2) )                &
     & + 2._dp * (3._dp * TraceY(3) + TraceY(1) ) * MH2       &
     & + 2._dp * Real( 6._dp * cTrace(YuMT2aYu)               &
     &               + cTrace(YnMF2aYn) + cTrace(aYnMN2Yn) ,dp)
  DMHbar2 = 9.6_dp * ( Alam2 - 2._dp * M52                           &
     &               + lam2 * (MH2 + MHbar2 + MS2) )                 &
     & + 8._dp * (TraceY(2) * MHbar2 + TraceA(2)                     &
     &           + Real( cTrace( aYdMT2Yd ) + cTrace(YdMF2aYd) , dp) ) 
  DMS2 = 2.1_dp * Alamp2 + 2._dp * Alamp2 - 40._dp * M52 &
     & + 6.3_dp * lamp2 * MS2 + 2._dp * lam2 * (MH2 + MHbar2 + MS2) 

  DMT2 = oo16pi2 * DMT2
  DMF2 = oo16pi2 * DMF2
  DMN2 = oo16pi2 * DMN2
  DMH2 = oo16pi2 * DMH2
  DMHbar2 = oo16pi2 * DMHbar2
  DMS2 = oo16pi2 * DMS2

  Do i1=1,3
   Do i2=1,3
    SumI = 6*i1+2*i2
    f(SumI+116) = Real(DMT2(i1,i2),dp)
    f(SumI+117) = Aimag(DMT2(i1,i2))
    f(SumI+134) = Real(DMF2(i1,i2),dp)
    f(SumI+135) = Aimag(DMF2(i1,i2))
    f(SumI+152) = Real(DMN2(i1,i2),dp)
    f(SumI+153) = Aimag(DMN2(i1,i2))
   End Do
  End Do

  f(178) = DMHbar2
  f(179) = DMH2
  f(180) = DMS2
  
  Iname = Iname - 1

 End Subroutine RGE_SU5


#ifdef SEESAWIII

 Subroutine GToParameters111(g,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3)

 Implicit None 
  Real(dp), Intent(in) :: g(111) 
  Real(dp),Intent(out) :: g1,g2,g3

  Complex(dp), Intent(out), Dimension(3,3) :: Yu, Yd, Ye, Yb3, Yw3, Yx3

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters111' 
 
  g1= g(1) 
  g2= g(2) 
  g3= g(3) 
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    Yb3(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    Yw3(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    Yx3(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
   End Do 
  End Do 
  
  Iname = Iname - 1 
 
 End Subroutine GToParameters111


 Subroutine ParametersToG111(g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,g)

 Implicit None 
  Real(dp), Intent(out) :: g(111) 
  Real(dp), Intent(in) :: g1,g2,g3

  Complex(dp), Intent(in), Dimension(3,3) :: Yu, Yd, Ye, Yb3, Yw3, Yx3

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG111' 
 
  g(1) = g1  
  g(2) = g2  
  g(3) = g3  

  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(Yb3(i1,i2), dp) 
    g(SumI+59) = Aimag(Yb3(i1,i2)) 
    g(SumI+76) = Real(Yw3(i1,i2), dp) 
    g(SumI+77) = Aimag(Yw3(i1,i2)) 
    g(SumI+94) = Real(Yx3(i1,i2), dp) 
    g(SumI+95) = Aimag(Yx3(i1,i2)) 
   End Do 
  End Do 

  Iname = Iname - 1 
 
 End Subroutine ParametersToG111


 Subroutine GToParameters353(g,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM,mue,      &
  & MZM,MSM,AYu,AYd,AYe,AYt,AYs,AYz,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2,ml2,mHd2, &
  & mHu2,md2,mu2,me2,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG)

 Implicit None 
  Real(dp), Intent(in) :: g(353) 
  Real(dp),Intent(out) :: g1,g2,g3,mHd2,mHu2

  Complex(dp),Intent(out) :: Yu(3,3),Yd(3,3),Ye(3,3),Yt(3,3),Ys(3,3),Yz(3,3)       &
    & ,L1,L2,MTM,mue,MZM,MSM,AYu(3,3),AYd(3,3),AYe(3,3),AYt(3,3),AYs(3,3),AYz(3,3) &
    & ,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3)    &
    & ,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG

  Integer i1, i2, SumI
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters353' 
 
  g1= g(1) 
  g2= g(2) 
  g3= g(3) 
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3)
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    Yt(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    Ys(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    Yz(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
    AYu(i1,i2) = Cmplx( g(SumI+124), g(SumI+125), dp) 
    AYd(i1,i2) = Cmplx( g(SumI+142), g(SumI+143), dp) 
    AYe(i1,i2) = Cmplx( g(SumI+160), g(SumI+161), dp) 
    AYt(i1,i2) = Cmplx( g(SumI+178), g(SumI+179), dp) 
    AYs(i1,i2) = Cmplx( g(SumI+196), g(SumI+197), dp) 
    AYz(i1,i2) = Cmplx( g(SumI+214), g(SumI+215), dp) 
    mq2(i1,i2) = Cmplx( g(SumI+244), g(SumI+245), dp) 
    ml2(i1,i2) = Cmplx( g(SumI+262), g(SumI+263), dp) 
    md2(i1,i2) = Cmplx( g(SumI+282), g(SumI+283), dp) 
    mu2(i1,i2) = Cmplx( g(SumI+300), g(SumI+301), dp) 
    me2(i1,i2) = Cmplx( g(SumI+318), g(SumI+319), dp) 
   End Do
  End Do
 
  L1= Cmplx(g(112),g(113),dp) 
  L2= Cmplx(g(114),g(115),dp) 
  MTM= Cmplx(g(116),g(117),dp) 
  mue= Cmplx(g(118),g(119),dp) 
  MZM= Cmplx(g(120),g(121),dp) 
  MSM= Cmplx(g(122),g(123),dp) 
  AL1= Cmplx(g(232),g(233),dp) 
  AL2= Cmplx(g(234),g(235),dp) 
  Amue= Cmplx(g(236),g(237),dp) 
  AMTM= Cmplx(g(238),g(239),dp) 
  AMZM= Cmplx(g(240),g(241),dp) 
  AMSM= Cmplx(g(242),g(243),dp) 
  mHd2= g(280) 
  mHu2= g(281)  
  mt2= Cmplx(g(336),g(337),dp) 
  mtb2= Cmplx(g(338),g(339),dp) 
  ms2= Cmplx(g(340),g(341),dp) 
  msb2= Cmplx(g(342),g(343),dp) 
  mz2= Cmplx(g(344),g(345),dp) 
  mzb2= Cmplx(g(346),g(347),dp) 
  MassB= Cmplx(g(348),g(349),dp) 
  MassWB= Cmplx(g(350),g(351),dp) 
  MassG= Cmplx(g(352),g(353),dp) 

  Iname = Iname - 1 
 
 End Subroutine GToParameters353


 Subroutine ParametersToG353(g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM,mue,MZM,     & 
  & MSM,AYu,AYd,AYe,AYt,AYs,AYz,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2,ml2,mHd2,mHu2, &
  & md2,mu2,me2,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG,g)

 Implicit None 
  Real(dp), Intent(out) :: g(353) 
  Real(dp), Intent(in) :: g1,g2,g3,mHd2,mHu2

  Complex(dp), Intent(in) :: Yu(3,3),Yd(3,3),Ye(3,3),Yt(3,3),Ys(3,3),Yz(3,3)       & 
    & ,L1,L2,MTM,mue,MZM,MSM,AYu(3,3),AYd(3,3),AYe(3,3),AYt(3,3),AYs(3,3),AYz(3,3) &
    & ,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3)    &
    & ,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG353' 
 
  g(1) = g1  
  g(2) = g2  
  g(3) = g3  
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(Yt(i1,i2), dp) 
    g(SumI+59) = Aimag(Yt(i1,i2)) 
    g(SumI+76) = Real(Ys(i1,i2), dp) 
    g(SumI+77) = Aimag(Ys(i1,i2)) 
    g(SumI+94) = Real(Yz(i1,i2), dp) 
    g(SumI+95) = Aimag(Yz(i1,i2)) 
    g(SumI+124) = Real(AYu(i1,i2), dp) 
    g(SumI+125) = Aimag(AYu(i1,i2)) 
    g(SumI+142) = Real(AYd(i1,i2), dp) 
    g(SumI+143) = Aimag(AYd(i1,i2)) 
    g(SumI+160) = Real(AYe(i1,i2), dp) 
    g(SumI+161) = Aimag(AYe(i1,i2)) 
    g(SumI+178) = Real(AYt(i1,i2), dp) 
    g(SumI+179) = Aimag(AYt(i1,i2)) 
    g(SumI+196) = Real(AYs(i1,i2), dp) 
    g(SumI+197) = Aimag(AYs(i1,i2)) 
    g(SumI+214) = Real(AYz(i1,i2), dp) 
    g(SumI+215) = Aimag(AYz(i1,i2)) 
    g(SumI+244) = Real(mq2(i1,i2), dp) 
    g(SumI+245) = Aimag(mq2(i1,i2)) 
    g(SumI+262) = Real(ml2(i1,i2), dp) 
    g(SumI+263) = Aimag(ml2(i1,i2)) 
    g(SumI+282) = Real(md2(i1,i2), dp) 
    g(SumI+283) = Aimag(md2(i1,i2)) 
    g(SumI+300) = Real(mu2(i1,i2), dp) 
    g(SumI+301) = Aimag(mu2(i1,i2)) 
    g(SumI+318) = Real(me2(i1,i2), dp) 
    g(SumI+319) = Aimag(me2(i1,i2)) 
   End Do
  End Do

  g(112) = Real(L1,dp)  
  g(113) = Aimag(L1)  
  g(114) = Real(L2,dp)  
  g(115) = Aimag(L2)  
  g(116) = Real(MTM,dp)  
  g(117) = Aimag(MTM)  
  g(118) = Real(mue,dp)  
  g(119) = Aimag(mue)  
  g(120) = Real(MZM,dp)  
  g(121) = Aimag(MZM)  
  g(122) = Real(MSM,dp)  
  g(123) = Aimag(MSM)
  g(232) = Real(AL1,dp)  
  g(233) = Aimag(AL1)  
  g(234) = Real(AL2,dp)  
  g(235) = Aimag(AL2)  
  g(236) = Real(Amue,dp)  
  g(237) = Aimag(Amue)  
  g(238) = Real(AMTM,dp)  
  g(239) = Aimag(AMTM)  
  g(240) = Real(AMZM,dp)  
  g(241) = Aimag(AMZM)  
  g(242) = Real(AMSM,dp)  
  g(243) = Aimag(AMSM)  
  g(280) = mHd2
  g(281) = mHu2
  g(336) = Real(mt2,dp)  
  g(337) = Aimag(mt2)  
  g(338) = Real(mtb2,dp)  
  g(339) = Aimag(mtb2)  
  g(340) = Real(ms2,dp)  
  g(341) = Aimag(ms2)  
  g(342) = Real(msb2,dp)  
  g(343) = Aimag(msb2)  
  g(344) = Real(mz2,dp)  
  g(345) = Aimag(mz2)  
  g(346) = Real(mzb2,dp)  
  g(347) = Aimag(mzb2)  
  g(348) = Real(MassB,dp)  
  g(349) = Aimag(MassB)  
  g(350) = Real(MassWB,dp)  
  g(351) = Aimag(MassWB)  
  g(352) = Real(MassG,dp)  
  g(353) = Aimag(MassG)  
  Iname = Iname - 1 
 
 End Subroutine ParametersToG353


 Subroutine rge117(len, T, GY, F)
 Implicit None
  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Real(dp) :: q
  Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,Dg3
  Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          &
  & ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       &
  & ,DYe(3,3),adjYe(3,3),Yt(3,3),betaYt1(3,3),betaYt2(3,3),DYt(3,3),adjYt(3,3)             &
  & ,Ys(3,3),betaYs1(3,3),betaYs2(3,3),DYs(3,3),adjYs(3,3),Yz(3,3),betaYz1(3,3)            &
  & ,betaYz2(3,3),DYz(3,3),adjYz(3,3),L1,betaL11,betaL12,DL1,L2,betaL21,betaL22,           &
  & DL2,MTM,betaMTM1,betaMTM2,DMTM
  Complex(dp) :: YdadjYd(3,3),YeadjYe(3,3),YsCYs(3,3),YtCYt(3,3),YuadjYu(3,3)           &
  & ,YzadjYz(3,3),adjYdYd(3,3),adjYdYs(3,3),adjYdYz(3,3),adjYeYe(3,3),adjYuYu(3,3)         &
  & ,adjYzYd(3,3),adjYzYs(3,3),adjYzYz(3,3),CYdTYd(3,3),CYeYt(3,3),CYsYd(3,3)              &
  & ,CYsYs(3,3),CYsYz(3,3),CYtYt(3,3),CYzYt(3,3),CYzTYz(3,3),YdadjYdYd(3,3),               &
  & YdadjYdYs(3,3),YdadjYdYz(3,3),YdadjYuYu(3,3),YeadjYeYe(3,3),YeadjYzYz(3,3)             &
  & ,YeCYtYt(3,3),YsCYdTYd(3,3),YsCYsYd(3,3),YsCYsYs(3,3),YsCYsYz(3,3),YsCYzTYz(3,3)       &
  & ,YtadjYeYe(3,3),YtadjYzYz(3,3),YtCYtYt(3,3),YuadjYdYd(3,3),YuadjYuYu(3,3)              &
  & ,YzadjYeYe(3,3),YzadjYzYd(3,3),YzadjYzYs(3,3),YzadjYzYz(3,3),YzCYtYt(3,3)              &
  & ,TYeCYeYt(3,3),TYzCYzYt(3,3)

  Complex(dp), Dimension(3,3) :: YdadjYu, YeadjYz, YeCYt, YtadjYe, YtadjYz     &
    & , YuadjYd, YzadjYe, YzCYt, CYeTYz, CYtTYz, CYuTYd, YeadjYzYd, YeadjYzYs  &
    & , YsCYzYt, YtadjYzYd, YtadjYzYs, YtCYtTYz, YuadjYdYs, YuadjYdYz          &
    & , adjYdYdadjYd, adjYdYdadjYu, adjYdYsCYs, adjYdYzadjYz, adjYeYeadjYe     &
    & , adjYeYeadjYz, adjYeYeCYt, adjYuYuadjYd, adjYuYuadjYu, adjYzYdadjYd     &
    & , adjYzYsCYs, adjYzYzadjYe, adjYzYzadjYz, adjYzYzCYt, CYsYdadjYd         &
    & , CYsYsCYs, CYsYzadjYz, CYtYtadjYe, CYtYtadjYz, CYtYtCYt, TYdCYdTYd      &
    & , TYdCYsYd, TYdCYsYs, TYdCYsYz, TYdCYzYt, TYeCYeTYz, TYuCYuTYd,TYzCYsYd  &
    & , TYzCYsYs, TYzCYsYz, TYzCYzTYz, YdadjYdYdadjYd, YdadjYdYsCYs            &
    & , YdadjYdYzadjYz, YdadjYuYuadjYd, YeadjYeYeadjYe, YeadjYzYzadjYe         &
    & , YeCYtYtadjYe, YsCYsYdadjYd, YsCYsYsCYs               &
  & , YsCYsYzadjYz, YtadjYzYzCYt, YtCYtYtCYt, YuadjYdYdadjYu&
  & , YuadjYuYuadjYu, YzadjYeYeadjYz, YzadjYzYdadjYd         &
  & , YzadjYzYzadjYz, adjYdYdadjYdYd, adjYdYdadjYdYs         &
  & , adjYdYdadjYdYz, adjYdYdadjYuYu, adjYdYsCYsYd, adjYdYzadjYzYd         &
  & , adjYeYeadjYeYe, adjYeYeadjYzYd, adjYeYeadjYzYs, adjYeYeadjYzYz       &
  & , adjYeYeCYtYt, adjYuYuadjYdYd, adjYuYuadjYdYs, adjYuYuadjYdYz         &
  & , adjYuYuadjYuYu, adjYzYdadjYdYz, adjYzYsCYsYz, adjYzYzadjYeYe         &
  & , adjYzYzadjYzYd, adjYzYzadjYzYs, adjYzYzadjYzYz, adjYzYzCYtYt         &
  & , CYdTYdCYdTYd, CYdTYdCYsYd, CYdTYdCYsYs, CYdTYdCYsYz, CYdTYdCYzYt &
  & , CYdTYuCYuTYd, CYeTYeCYeYt, CYsYdadjYdYs, CYsYsCYsYd, CYsYsCYsYs  &
  & , CYsYsCYsYz, CYsYsCYzYt, CYsYzadjYzYs, CYtYtadjYeYe, CYtYtadjYzYd &
  & , CYtYtadjYzYs, CYtYtadjYzYz, CYtYtCYtYt, CYtTYeCYeYt, CYtTYzCYzYt &
  & , CYzYtCYtTYz, CYzTYeCYeTYz, CYzTYzCYsYd, CYzTYzCYsYs, CYzTYzCYsYz &
  & , CYzTYzCYzYt, CYzTYzCYzTYz, YdadjYdYdadjYdYd, YdadjYdYdadjYdYs        &
  & , YdadjYdYdadjYdYz, YdadjYdYsCYsYd, YdadjYdYzadjYzYd, YdadjYuYuadjYdYd &
  & , YdadjYuYuadjYdYs, YdadjYuYuadjYdYz, YdadjYuYuadjYuYu, YeadjYeYeadjYeYe&
  & , YeadjYzYdadjYdYz, YeadjYzYsCYsYz, YeadjYzYzadjYeYe, YeadjYzYzadjYzYz &
  & , YeCYtYtadjYeYe, YeCYtYtCYtYt, YeCYtTYeCYeYt, YeCYtTYzCYzYt           &
  & , YsCYdTYdCYdTYd, YsCYdTYdCYsYd, YsCYdTYdCYsYs, YsCYdTYdCYsYz          &
  & , YsCYdTYuCYuTYd, YsCYsYdadjYdYs, YsCYsYsCYsYd, YsCYsYsCYsYs           &
  & , YsCYsYsCYsYz, YsCYsYzadjYzYs, YsCYzYtCYtTYz, YsCYzTYeCYeTYz          &
  & , YsCYzTYzCYsYd, YsCYzTYzCYsYs, YsCYzTYzCYsYz, YsCYzTYzCYzTYz          &
  & , YtadjYeYeadjYeYe, YtadjYeYeCYtYt, YtadjYzYdadjYdYz, YtadjYzYsCYsYz   &
  & , YtadjYzYzadjYzYz, YtadjYzYzCYtYt, YtCYtYtCYtYt, YtCYtTYeCYeYt        &
  & , YtCYtTYzCYzYt, YuadjYdYdadjYdYd, YuadjYdYdadjYuYu, YuadjYdYsCYsYd    &
  & , YuadjYdYzadjYzYd, YuadjYuYuadjYuYu, YzadjYeYeadjYeYe, YzadjYeYeadjYzYd&
  & , YzadjYeYeadjYzYs, YzadjYeYeadjYzYz, YzadjYzYdadjYdYz, YzadjYzYsCYsYz &
  & , YzadjYzYzadjYzYd, YzadjYzYzadjYzYs, YzadjYzYzadjYzYz, YzCYtYtadjYzYd &
  & , YzCYtYtadjYzYs, YzCYtYtadjYzYz, YzCYtYtCYtYt, YzCYtTYeCYeYt          &
  & , YzCYtTYzCYzYt, TYeCYeTYeCYeYt, TYzCYdTYdCYzYt, TYzCYsYsCYzYt         &
  & , TYzCYzTYzCYzYt

  Real(dp) :: TrCYtYt, TrCYsYs, TrYdadjYd, TrYeadjYe, TrYsCYs, TrYtCYt     &
   & , TrYuadjYu ,TrYzadjYz

  Complex(dp) :: TrYdadjYdYdadjYd,TrYdadjYdYsCYs,TrYdadjYdYzadjYz,TrYdadjYuYuadjYd,     &
  & TrYeadjYeYeadjYe,TrYeadjYzYzadjYe,TrYeCYtYtadjYe,TrYsCYsYdadjYd,          &
  & TrYsCYsYzadjYz,TrYtadjYzYzCYt,TrYtCYtYtCYt,TrYuadjYdYdadjYu,            &
  & TrYuadjYuYuadjYu,TrYzadjYeYeadjYz,TrYzadjYzYdadjYd,TrYzadjYzYzadjYz,    &
  & TrYzCYtYtadjYz,TrCYsYdadjYdYs,TrCYsYsCYsYs,TrCYsYzadjYzYs,TrCYtYtadjYeYe,              &
  & TrCYtYtadjYzYz,TrCYtYtCYtYt

  Real(dp) :: g12, g13, g14, g22, g23, g24, g32, g33, g34, AbsL1sq, AbsL2sq

  Iname = Iname +1
  NameOfUnit(Iname) = 'rge117'

  OnlyDiagonal = .Not.GenerationMixing

  q = t

  Call GToParameters117(gy,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM)

  g12 = g1**2
  g13 = g1*g12
  g14 = g12*g12
  g22 = g2**2
  g23 = g2*g22
  g24 = g22*g22
  g32 = g3**2
  g33 = g3*g32
  g34 = g32*g32
  AbsL1sq = Abs(L1)**2
  AbsL2sq = Abs(L2)**2

  Call Adjungate(Yu,adjYu)
  Call Adjungate(Yd,adjYd)
  Call Adjungate(Ye,adjYe)
  adjYt = Conjg(Yt)
  adjYs = Conjg(Ys)
  Call Adjungate(Yz,adjYz)
  YdadjYd = Matmul2(Yd,adjYd,OnlyDiagonal)
  YeadjYe = Matmul2(Ye,adjYe,OnlyDiagonal)
  YsCYs = Matmul2(Ys,adjYs,OnlyDiagonal)
  YtCYt = Matmul2(Yt,adjYt,OnlyDiagonal)
  YuadjYu = Matmul2(Yu,adjYu,OnlyDiagonal)
  YzadjYz = Matmul2(Yz,adjYz,OnlyDiagonal)
  adjYdYd = Matmul2(adjYd,Yd,OnlyDiagonal)
  adjYdYs = Matmul2(adjYd,Ys,OnlyDiagonal)
  adjYdYz = Matmul2(adjYd,Yz,OnlyDiagonal)
  adjYeYe = Matmul2(adjYe,Ye,OnlyDiagonal)
  adjYuYu = Matmul2(adjYu,Yu,OnlyDiagonal)
  adjYzYd = Matmul2(adjYz,Yd,OnlyDiagonal)
  adjYzYs = Matmul2(adjYz,Ys,OnlyDiagonal)
  adjYzYz = Matmul2(adjYz,Yz,OnlyDiagonal)

  CYdTYd = Matmul2(Conjg(Yd),Transpose(Yd),OnlyDiagonal)
  CYeYt = Matmul2(Conjg(Ye),Yt,OnlyDiagonal)
  CYsYd = Matmul2(adjYs,Yd,OnlyDiagonal)
  CYsYs = Matmul2(adjYs,Ys,OnlyDiagonal)
  CYsYz = Matmul2(adjYs,Yz,OnlyDiagonal)
  CYtYt = Matmul2(adjYt,Yt,OnlyDiagonal)
  CYzYt = Matmul2(Conjg(Yz),Yt,OnlyDiagonal)
  CYzTYz = Matmul2(Conjg(Yz),Transpose(Yz),OnlyDiagonal)

  YdadjYdYd = Matmul2(Yd,adjYdYd,OnlyDiagonal)
  YdadjYdYs = Matmul2(Yd,adjYdYs,OnlyDiagonal)
  YdadjYdYz = Matmul2(Yd,adjYdYz,OnlyDiagonal)
  YdadjYuYu = Matmul2(Yd,adjYuYu,OnlyDiagonal)
  YeadjYeYe = Matmul2(Ye,adjYeYe,OnlyDiagonal)
  YeadjYzYz = Matmul2(Ye,adjYzYz,OnlyDiagonal)
  YeCYtYt = Matmul2(Ye,CYtYt,OnlyDiagonal)
  YsCYdTYd = Matmul2(Ys,CYdTYd,OnlyDiagonal)
  YsCYsYd = Matmul2(Ys,CYsYd,OnlyDiagonal)
  YsCYsYs = Matmul2(Ys,CYsYs,OnlyDiagonal)
  YsCYsYz = Matmul2(Ys,CYsYz,OnlyDiagonal)
  YsCYzTYz = Matmul2(Ys,CYzTYz,OnlyDiagonal)
  YtadjYeYe = Matmul2(Yt,adjYeYe,OnlyDiagonal)
  YtadjYzYz = Matmul2(Yt,adjYzYz,OnlyDiagonal)
  YtCYtYt = Matmul2(Yt,CYtYt,OnlyDiagonal)
  YuadjYdYd = Matmul2(Yu,adjYdYd,OnlyDiagonal)
  YuadjYuYu = Matmul2(Yu,adjYuYu,OnlyDiagonal)
  YzadjYeYe = Matmul2(Yz,adjYeYe,OnlyDiagonal)
  YzadjYzYd = Matmul2(Yz,adjYzYd,OnlyDiagonal)
  YzadjYzYs = Matmul2(Yz,adjYzYs,OnlyDiagonal)
  YzadjYzYz = Matmul2(Yz,adjYzYz,OnlyDiagonal)
  YzCYtYt = Matmul2(Yz,CYtYt,OnlyDiagonal)
  TYeCYeYt = Matmul2(Transpose(Ye),CYeYt,OnlyDiagonal)
  TYzCYzYt = Matmul2(Transpose(Yz),CYzYt,OnlyDiagonal)
  TrYdadjYd = Real(cTrace(YdadjYd),dp)
  TrYeadjYe = Real(cTrace(YeadjYe),dp)
  TrYsCYs = Real(cTrace(YsCYs),dp)
  TrYtCYt = Real(cTrace(YtCYt),dp)
  TrYuadjYu = Real(cTrace(YuadjYu),dp)
  TrYzadjYz = Real(cTrace(YzadjYz),dp)
  TrCYsYs = TrYsCYs ! Real(cTrace(CYsYs),dp)
  TrCYtYt = TrYtCYt ! Real(cTrace(CYtYt),dp)


  If (TwoLoopRGE) Then
  YdadjYu = Matmul2(Yd,adjYu,OnlyDiagonal)
  YeadjYz = Matmul2(Ye,adjYz,OnlyDiagonal)
  YeCYt = Matmul2(Ye,adjYt,OnlyDiagonal)
  YtadjYe = Matmul2(Yt,adjYe,OnlyDiagonal)
  YtadjYz = Matmul2(Yt,adjYz,OnlyDiagonal)
  YuadjYd = Matmul2(Yu,adjYd,OnlyDiagonal)
  YzadjYe = Matmul2(Yz,adjYe,OnlyDiagonal)
  YzCYt = Matmul2(Yz,adjYt,OnlyDiagonal)
  CYeTYz = Matmul2(Conjg(Ye),Transpose(Yz),OnlyDiagonal)
  CYtTYz = Matmul2(adjYt,Transpose(Yz),OnlyDiagonal)
  CYuTYd = Matmul2(Conjg(Yu),Transpose(Yd),OnlyDiagonal)
  YeadjYzYd = Matmul2(Ye,adjYzYd,OnlyDiagonal)
  YeadjYzYs = Matmul2(Ye,adjYzYs,OnlyDiagonal)
  YsCYzYt = Matmul2(Ys,CYzYt,OnlyDiagonal)
  YtadjYzYd = Matmul2(Yt,adjYzYd,OnlyDiagonal)
  YtadjYzYs = Matmul2(Yt,adjYzYs,OnlyDiagonal)
  YtCYtTYz = Matmul2(Yt,CYtTYz,OnlyDiagonal)
  YuadjYdYs = Matmul2(Yu,adjYdYs,OnlyDiagonal)
  YuadjYdYz = Matmul2(Yu,adjYdYz,OnlyDiagonal)
  adjYdYdadjYd = Matmul2(adjYd,YdadjYd,OnlyDiagonal)
  adjYdYdadjYu = Matmul2(adjYd,YdadjYu,OnlyDiagonal)
  adjYdYsCYs = Matmul2(adjYd,YsCYs,OnlyDiagonal)
  adjYdYzadjYz = Matmul2(adjYd,YzadjYz,OnlyDiagonal)
  adjYeYeadjYe = Matmul2(adjYe,YeadjYe,OnlyDiagonal)
  adjYeYeadjYz = Matmul2(adjYe,YeadjYz,OnlyDiagonal)
  adjYeYeCYt = Matmul2(adjYe,YeCYt,OnlyDiagonal)
  adjYuYuadjYd = Matmul2(adjYu,YuadjYd,OnlyDiagonal)
  adjYuYuadjYu = Matmul2(adjYu,YuadjYu,OnlyDiagonal)
  adjYzYdadjYd = Matmul2(adjYz,YdadjYd,OnlyDiagonal)
  adjYzYsCYs = Matmul2(adjYz,YsCYs,OnlyDiagonal)
  adjYzYzadjYe = Matmul2(adjYz,YzadjYe,OnlyDiagonal)
  adjYzYzadjYz = Matmul2(adjYz,YzadjYz,OnlyDiagonal)
  adjYzYzCYt = Matmul2(adjYz,YzCYt,OnlyDiagonal)
  CYsYdadjYd = Matmul2(adjYs,YdadjYd,OnlyDiagonal)
  CYsYsCYs = Matmul2(adjYs,YsCYs,OnlyDiagonal)
  CYsYzadjYz = Matmul2(adjYs,YzadjYz,OnlyDiagonal)
  CYtYtadjYe = Matmul2(adjYt,YtadjYe,OnlyDiagonal)
  CYtYtadjYz = Matmul2(adjYt,YtadjYz,OnlyDiagonal)
  CYtYtCYt = Matmul2(adjYt,YtCYt,OnlyDiagonal)
  TYdCYdTYd = Matmul2(Transpose(Yd),CYdTYd,OnlyDiagonal)
  TYdCYsYd = Matmul2(Transpose(Yd),CYsYd,OnlyDiagonal)
  TYdCYsYs = Matmul2(Transpose(Yd),CYsYs,OnlyDiagonal)
  TYdCYsYz = Matmul2(Transpose(Yd),CYsYz,OnlyDiagonal)
  TYdCYzYt = Matmul2(Transpose(Yd),CYzYt,OnlyDiagonal)
  TYeCYeTYz = Matmul2(Transpose(Ye),CYeTYz,OnlyDiagonal)
  TYuCYuTYd = Matmul2(Transpose(Yu),CYuTYd,OnlyDiagonal)
  TYzCYsYd = Matmul2(Transpose(Yz),CYsYd,OnlyDiagonal)
  TYzCYsYs = Matmul2(Transpose(Yz),CYsYs,OnlyDiagonal)
  TYzCYsYz = Matmul2(Transpose(Yz),CYsYz,OnlyDiagonal)
  TYzCYzTYz = Matmul2(Transpose(Yz),CYzTYz,OnlyDiagonal)
  YdadjYdYdadjYd = Matmul2(Yd,adjYdYdadjYd,OnlyDiagonal)
  YdadjYdYsCYs = Matmul2(Yd,adjYdYsCYs,OnlyDiagonal)
  YdadjYdYzadjYz = Matmul2(Yd,adjYdYzadjYz,OnlyDiagonal)
  YdadjYuYuadjYd = Matmul2(Yd,adjYuYuadjYd,OnlyDiagonal)
  YeadjYeYeadjYe = Matmul2(Ye,adjYeYeadjYe,OnlyDiagonal)
  YeadjYzYzadjYe = Matmul2(Ye,adjYzYzadjYe,OnlyDiagonal)
  YeCYtYtadjYe = Matmul2(Ye,CYtYtadjYe,OnlyDiagonal)
  YsCYsYdadjYd = Matmul2(Ys,CYsYdadjYd,OnlyDiagonal)
  YsCYsYsCYs = Matmul2(Ys,CYsYsCYs,OnlyDiagonal)
  YsCYsYzadjYz = Matmul2(Ys,CYsYzadjYz,OnlyDiagonal)
  YtadjYzYzCYt = Matmul2(Yt,adjYzYzCYt,OnlyDiagonal)
  YtCYtYtCYt = Matmul2(Yt,CYtYtCYt,OnlyDiagonal)
  YuadjYdYdadjYu = Matmul2(Yu,adjYdYdadjYu,OnlyDiagonal)
  YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal)
  YzadjYeYeadjYz = Matmul2(Yz,adjYeYeadjYz,OnlyDiagonal)
  YzadjYzYdadjYd = Matmul2(Yz,adjYzYdadjYd,OnlyDiagonal)
  YzadjYzYzadjYz = Matmul2(Yz,adjYzYzadjYz,OnlyDiagonal)
  adjYdYdadjYdYd = Matmul2(adjYd,YdadjYdYd,OnlyDiagonal)
  adjYdYdadjYdYs = Matmul2(adjYd,YdadjYdYs,OnlyDiagonal)
  adjYdYdadjYdYz = Matmul2(adjYd,YdadjYdYz,OnlyDiagonal)
  adjYdYdadjYuYu = Matmul2(adjYd,YdadjYuYu,OnlyDiagonal)
  adjYdYsCYsYd = Matmul2(adjYd,YsCYsYd,OnlyDiagonal)
  adjYdYzadjYzYd = Matmul2(adjYd,YzadjYzYd,OnlyDiagonal)
  adjYeYeadjYeYe = Matmul2(adjYe,YeadjYeYe,OnlyDiagonal)
  adjYeYeadjYzYd = Matmul2(adjYe,YeadjYzYd,OnlyDiagonal)
  adjYeYeadjYzYs = Matmul2(adjYe,YeadjYzYs,OnlyDiagonal)
  adjYeYeadjYzYz = Matmul2(adjYe,YeadjYzYz,OnlyDiagonal)
  adjYeYeCYtYt = Matmul2(adjYe,YeCYtYt,OnlyDiagonal)
  adjYuYuadjYdYd = Matmul2(adjYu,YuadjYdYd,OnlyDiagonal)
  adjYuYuadjYdYs = Matmul2(adjYu,YuadjYdYs,OnlyDiagonal)
  adjYuYuadjYdYz = Matmul2(adjYu,YuadjYdYz,OnlyDiagonal)
  adjYuYuadjYuYu = Matmul2(adjYu,YuadjYuYu,OnlyDiagonal)
  adjYzYdadjYdYz = Matmul2(adjYz,YdadjYdYz,OnlyDiagonal)
  adjYzYsCYsYz = Matmul2(adjYz,YsCYsYz,OnlyDiagonal)
  adjYzYzadjYeYe = Matmul2(adjYz,YzadjYeYe,OnlyDiagonal)
  adjYzYzadjYzYd = Matmul2(adjYz,YzadjYzYd,OnlyDiagonal)
  adjYzYzadjYzYs = Matmul2(adjYz,YzadjYzYs,OnlyDiagonal)
  adjYzYzadjYzYz = Matmul2(adjYz,YzadjYzYz,OnlyDiagonal)
  adjYzYzCYtYt = Matmul2(adjYz,YzCYtYt,OnlyDiagonal)
  CYdTYdCYdTYd = Matmul2(Conjg(Yd),TYdCYdTYd,OnlyDiagonal)
  CYdTYdCYsYd = Matmul2(Conjg(Yd),TYdCYsYd,OnlyDiagonal)
  CYdTYdCYsYs = Matmul2(Conjg(Yd),TYdCYsYs,OnlyDiagonal)
  CYdTYdCYsYz = Matmul2(Conjg(Yd),TYdCYsYz,OnlyDiagonal)
  CYdTYdCYzYt = Matmul2(Conjg(Yd),TYdCYzYt,OnlyDiagonal)
  CYdTYuCYuTYd = Matmul2(Conjg(Yd),TYuCYuTYd,OnlyDiagonal)
  CYeTYeCYeYt = Matmul2(Conjg(Ye),TYeCYeYt,OnlyDiagonal)
  CYsYdadjYdYs = Matmul2(adjYs,YdadjYdYs,OnlyDiagonal)
  CYsYsCYsYd = Matmul2(adjYs,YsCYsYd,OnlyDiagonal)
  CYsYsCYsYs = Matmul2(adjYs,YsCYsYs,OnlyDiagonal)
  CYsYsCYsYz = Matmul2(adjYs,YsCYsYz,OnlyDiagonal)
  CYsYsCYzYt = Matmul2(adjYs,YsCYzYt,OnlyDiagonal)
  CYsYzadjYzYs = Matmul2(adjYs,YzadjYzYs,OnlyDiagonal)
  CYtYtadjYeYe = Matmul2(adjYt,YtadjYeYe,OnlyDiagonal)
  CYtYtadjYzYd = Matmul2(adjYt,YtadjYzYd,OnlyDiagonal)
  CYtYtadjYzYs = Matmul2(adjYt,YtadjYzYs,OnlyDiagonal)
  CYtYtadjYzYz = Matmul2(adjYt,YtadjYzYz,OnlyDiagonal)
  CYtYtCYtYt = Matmul2(adjYt,YtCYtYt,OnlyDiagonal)
  CYtTYeCYeYt = Matmul2(adjYt,TYeCYeYt,OnlyDiagonal)
  CYtTYzCYzYt = Matmul2(adjYt,TYzCYzYt,OnlyDiagonal)
  CYzYtCYtTYz = Matmul2(Conjg(Yz),YtCYtTYz,OnlyDiagonal)
  CYzTYeCYeTYz = Matmul2(Conjg(Yz),TYeCYeTYz,OnlyDiagonal)
  CYzTYzCYsYd = Matmul2(Conjg(Yz),TYzCYsYd,OnlyDiagonal)
  CYzTYzCYsYs = Matmul2(Conjg(Yz),TYzCYsYs,OnlyDiagonal)
  CYzTYzCYsYz = Matmul2(Conjg(Yz),TYzCYsYz,OnlyDiagonal)
  CYzTYzCYzYt = Matmul2(Conjg(Yz),TYzCYzYt,OnlyDiagonal)
  CYzTYzCYzTYz = Matmul2(Conjg(Yz),TYzCYzTYz,OnlyDiagonal)
  YdadjYdYdadjYdYd = Matmul2(Yd,adjYdYdadjYdYd,OnlyDiagonal)
  YdadjYdYdadjYdYs = Matmul2(Yd,adjYdYdadjYdYs,OnlyDiagonal)
  YdadjYdYdadjYdYz = Matmul2(Yd,adjYdYdadjYdYz,OnlyDiagonal)
  YdadjYdYsCYsYd = Matmul2(Yd,adjYdYsCYsYd,OnlyDiagonal)
  YdadjYdYzadjYzYd = Matmul2(Yd,adjYdYzadjYzYd,OnlyDiagonal)
  YdadjYuYuadjYdYd = Matmul2(Yd,adjYuYuadjYdYd,OnlyDiagonal)
  YdadjYuYuadjYdYs = Matmul2(Yd,adjYuYuadjYdYs,OnlyDiagonal)
  YdadjYuYuadjYdYz = Matmul2(Yd,adjYuYuadjYdYz,OnlyDiagonal)
  YdadjYuYuadjYuYu = Matmul2(Yd,adjYuYuadjYuYu,OnlyDiagonal)
  YeadjYeYeadjYeYe = Matmul2(Ye,adjYeYeadjYeYe,OnlyDiagonal)
  YeadjYzYdadjYdYz = Matmul2(Ye,adjYzYdadjYdYz,OnlyDiagonal)
  YeadjYzYsCYsYz = Matmul2(Ye,adjYzYsCYsYz,OnlyDiagonal)
  YeadjYzYzadjYeYe = Matmul2(Ye,adjYzYzadjYeYe,OnlyDiagonal)
  YeadjYzYzadjYzYz = Matmul2(Ye,adjYzYzadjYzYz,OnlyDiagonal)
  YeCYtYtadjYeYe = Matmul2(Ye,CYtYtadjYeYe,OnlyDiagonal)
  YeCYtYtCYtYt = Matmul2(Ye,CYtYtCYtYt,OnlyDiagonal)
  YeCYtTYeCYeYt = Matmul2(Ye,CYtTYeCYeYt,OnlyDiagonal)
  YeCYtTYzCYzYt = Matmul2(Ye,CYtTYzCYzYt,OnlyDiagonal)
  YsCYdTYdCYdTYd = Matmul2(Ys,CYdTYdCYdTYd,OnlyDiagonal)
  YsCYdTYdCYsYd = Matmul2(Ys,CYdTYdCYsYd,OnlyDiagonal)
  YsCYdTYdCYsYs = Matmul2(Ys,CYdTYdCYsYs,OnlyDiagonal)
  YsCYdTYdCYsYz = Matmul2(Ys,CYdTYdCYsYz,OnlyDiagonal)
  YsCYdTYuCYuTYd = Matmul2(Ys,CYdTYuCYuTYd,OnlyDiagonal)
  YsCYsYdadjYdYs = Matmul2(Ys,CYsYdadjYdYs,OnlyDiagonal)
  YsCYsYsCYsYd = Matmul2(Ys,CYsYsCYsYd,OnlyDiagonal)
  YsCYsYsCYsYs = Matmul2(Ys,CYsYsCYsYs,OnlyDiagonal)
  YsCYsYsCYsYz = Matmul2(Ys,CYsYsCYsYz,OnlyDiagonal)
  YsCYsYzadjYzYs = Matmul2(Ys,CYsYzadjYzYs,OnlyDiagonal)
  YsCYzYtCYtTYz = Matmul2(Ys,CYzYtCYtTYz,OnlyDiagonal)
  YsCYzTYeCYeTYz = Matmul2(Ys,CYzTYeCYeTYz,OnlyDiagonal)
  YsCYzTYzCYsYd = Matmul2(Ys,CYzTYzCYsYd,OnlyDiagonal)
  YsCYzTYzCYsYs = Matmul2(Ys,CYzTYzCYsYs,OnlyDiagonal)
  YsCYzTYzCYsYz = Matmul2(Ys,CYzTYzCYsYz,OnlyDiagonal)
  YsCYzTYzCYzTYz = Matmul2(Ys,CYzTYzCYzTYz,OnlyDiagonal)
  YtadjYeYeadjYeYe = Matmul2(Yt,adjYeYeadjYeYe,OnlyDiagonal)
  YtadjYeYeCYtYt = Matmul2(Yt,adjYeYeCYtYt,OnlyDiagonal)
  YtadjYzYdadjYdYz = Matmul2(Yt,adjYzYdadjYdYz,OnlyDiagonal)
  YtadjYzYsCYsYz = Matmul2(Yt,adjYzYsCYsYz,OnlyDiagonal)
  YtadjYzYzadjYzYz = Matmul2(Yt,adjYzYzadjYzYz,OnlyDiagonal)
  YtadjYzYzCYtYt = Matmul2(Yt,adjYzYzCYtYt,OnlyDiagonal)
  YtCYtYtCYtYt = Matmul2(Yt,CYtYtCYtYt,OnlyDiagonal)
  YtCYtTYeCYeYt = Matmul2(Yt,CYtTYeCYeYt,OnlyDiagonal)
  YtCYtTYzCYzYt = Matmul2(Yt,CYtTYzCYzYt,OnlyDiagonal)
  YuadjYdYdadjYdYd = Matmul2(Yu,adjYdYdadjYdYd,OnlyDiagonal)
  YuadjYdYdadjYuYu = Matmul2(Yu,adjYdYdadjYuYu,OnlyDiagonal)
  YuadjYdYsCYsYd = Matmul2(Yu,adjYdYsCYsYd,OnlyDiagonal)
  YuadjYdYzadjYzYd = Matmul2(Yu,adjYdYzadjYzYd,OnlyDiagonal)
  YuadjYuYuadjYuYu = Matmul2(Yu,adjYuYuadjYuYu,OnlyDiagonal)
  YzadjYeYeadjYeYe = Matmul2(Yz,adjYeYeadjYeYe,OnlyDiagonal)
  YzadjYeYeadjYzYd = Matmul2(Yz,adjYeYeadjYzYd,OnlyDiagonal)
  YzadjYeYeadjYzYs = Matmul2(Yz,adjYeYeadjYzYs,OnlyDiagonal)
  YzadjYeYeadjYzYz = Matmul2(Yz,adjYeYeadjYzYz,OnlyDiagonal)
  YzadjYzYdadjYdYz = Matmul2(Yz,adjYzYdadjYdYz,OnlyDiagonal)
  YzadjYzYsCYsYz = Matmul2(Yz,adjYzYsCYsYz,OnlyDiagonal)
  YzadjYzYzadjYzYd = Matmul2(Yz,adjYzYzadjYzYd,OnlyDiagonal)
  YzadjYzYzadjYzYs = Matmul2(Yz,adjYzYzadjYzYs,OnlyDiagonal)
  YzadjYzYzadjYzYz = Matmul2(Yz,adjYzYzadjYzYz,OnlyDiagonal)
  YzCYtYtadjYzYd = Matmul2(Yz,CYtYtadjYzYd,OnlyDiagonal)
  YzCYtYtadjYzYs = Matmul2(Yz,CYtYtadjYzYs,OnlyDiagonal)
  YzCYtYtadjYzYz = Matmul2(Yz,CYtYtadjYzYz,OnlyDiagonal)
  YzCYtYtCYtYt = Matmul2(Yz,CYtYtCYtYt,OnlyDiagonal)
  YzCYtTYeCYeYt = Matmul2(Yz,CYtTYeCYeYt,OnlyDiagonal)
  YzCYtTYzCYzYt = Matmul2(Yz,CYtTYzCYzYt,OnlyDiagonal)
  TYeCYeTYeCYeYt = Matmul2(Transpose(Ye),CYeTYeCYeYt,OnlyDiagonal)
  TYzCYdTYdCYzYt = Matmul2(Transpose(Yz),CYdTYdCYzYt,OnlyDiagonal)
  TYzCYsYsCYzYt = Matmul2(Transpose(Yz),CYsYsCYzYt,OnlyDiagonal)
  TYzCYzTYzCYzYt = Matmul2(Transpose(Yz),CYzTYzCYzYt,OnlyDiagonal)

  TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd)

  TrYdadjYdYsCYs = cTrace(YdadjYdYsCYs)
  TrCYsYdadjYdYs = TrYdadjYdYsCYs

  TrYdadjYdYzadjYz = cTrace(YdadjYdYzadjYz)
  TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd)
  TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe)
  TrYeadjYzYzadjYe = cTrace(YeadjYzYzadjYe)

  TrYeCYtYtadjYe = cTrace(YeCYtYtadjYe)
  TrCYtYtadjYeYe = TrYeCYtYtadjYe

  TrYsCYsYdadjYd = cTrace(YsCYsYdadjYd)

  TrYsCYsYzadjYz = cTrace(YsCYsYzadjYz)
  TrCYsYzadjYzYs = TrYsCYsYzadjYz

  TrYtadjYzYzCYt = cTrace(YtadjYzYzCYt)
  TrCYtYtadjYzYz = TrYtadjYzYzCYt
  TrYzCYtYtadjYz = TrYtadjYzYzCYt

  TrYtCYtYtCYt = cTrace(YtCYtYtCYt)
  TrCYtYtCYtYt = TrYtCYtYtCYt

  TrYuadjYdYdadjYu = cTrace(YuadjYdYdadjYu)
  TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu)
  TrYzadjYeYeadjYz = cTrace(YzadjYeYeadjYz)
  TrYzadjYzYdadjYd = cTrace(YzadjYzYdadjYd)

  TrYzadjYzYzadjYz = cTrace(YzadjYzYzadjYz)

  TrCYsYsCYsYs = cTrace(CYsYsCYsYs)

  End If


  !--------------------
  ! g1
  !--------------------

  betag11 = 13.6_dp

  If (TwoLoopRGE) Then
   betag12 = ( 1502._dp*g12/75._dp + 34.8_dp * g22 + 184._dp * g32/3._dp  &
           & - 4.8_dp * TrCYsYs - 5.4_dp * (TrCYtYt + AbsL1sq+AbsL2sq)    &
           & - 2.8_dp * (TrYdadjYd+TrYzadjYz)  - 3.6_dp * TrYeadjYe       &
           & - 5.2_dp * TrYuadjYu )

   Dg1 = oo16pi2*g13*( betag11 + oo16pi2 * betag12 )

  Else
   Dg1 = oo16pi2* betag11*g13
  End If

  !--------------------
  ! g2
  !--------------------

  betag21 = 8._dp

  If (TwoLoopRGE) Then
   betag22 = ( 11.6_dp * g12 + 94._dp * g22 + 40._dp * g32                     &
           & - 6._dp * (TrYdadjYd + TrYuadjYu + TrYzadjYz) - 2._dp * TrYeadjYe &
           & - 7._dp*(TrCYtYt + AbsL1sq + AbsL2sq) )

   Dg2 = oo16pi2*g23*( betag21 + oo16pi2 * betag22 )

  Else
   Dg2 = oo16pi2*g23* betag21
  End If

  !--------------------
  ! g3
  !--------------------

  betag31 = 4._dp

  If (TwoLoopRGE) Then
   betag32 = ( 23._dp*g12/3._dp + 15._dp*g22 + 400._dp*g32/3._dp       &
           & - 9._dp*TrCYsYs - 4._dp*(TrYdadjYd+TrYuadjYu+TrYzadjYz) )

   Dg3 = oo16pi2*g33*( betag31 + oo16pi2 * betag32 )

  Else
   Dg3 = oo16pi2*g33* betag31
  End If

  !--------------------
  ! Yu
  !--------------------

  betaYu1 = ( -13._dp * g12 / 15._dp - 3._dp * g22 - 16._dp * g32 /3._dp &
          & + 3._dp * (TrYuadjYu + AbsL2sq ) ) * Yu                      &
          & + YuadjYdYd + 3._dp * YuadjYuYu

  If (TwoLoopRGE) Then
   betaYu2 = ( 5473._dp*g14/450._dp + (g12+8._dp*g32)*g22+ 28.5_dp*g24        &
      & + 136._dp*g12*g32/45._dp + 320._dp*g34/9._dp - 3._dp*TrYuadjYdYdadjYu &
      & + 0.8_dp*g12*TrYuadjYu + 16._dp*g32*TrYuadjYu- 9._dp*TrYuadjYuYuadjYu &
      & + (3.6_dp*g12+ 12._dp*g22 - 9._dp*TrYuadjYu - 12._dp*AbsL2sq)*AbsL2sq &
      & ) * Yu                                                                &
      & + (0.4_dp*g12 - 3._dp*(TrYdadjYd+AbsL1sq) - TrYeadjYe ) * YuadjYdYd   &
      & + (0.4_dp*g12 + 6._dp*g22  - 9._dp*(TrYuadjYu+AbsL2sq) ) * YuadjYuYu  &
      & - 2._dp * (YuadjYdYdadjYdYd + YuadjYdYdadjYuYu + YuadjYdYzadjYzYd )   &
      & - 4._dp * (YuadjYdYsCYsYd + YuadjYuYuadjYuYu)

   DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 )

  Else
   DYu = oo16pi2* betaYu1
  End If


  !--------------------
  ! Yd
  !--------------------

  betaYd1 = (-7._dp*g12/15._dp - 3._dp*g22 - 16._dp*g32/3._dp   &
          & + 3._dp*(TrYdadjYd+AbsL1sq) + TrYeadjYe )*Yd        &
          & + 3._dp*YdadjYdYd + YdadjYuYu + 4._dp*YsCYsYd + 2._dp*YzadjYzYd

  If (TwoLoopRGE) Then

   betaYd2 = ( 581._dp*g14/90._dp + (8._dp*g32+g12)*g22 + 28.5_dp*g24          &
          & + (8._dp*g12*g32+320._dp*g34)/9._dp - 9._dp*TrYdadjYdYdadjYd       &
          & + (16._dp*g32 - 0.4_dp*g12 - 9._dp*AbsL1sq ) *TrYdadjYd            &
          & - 3._dp *(TrYdadjYuYuadjYd+TrYeadjYeYeadjYe) + 1.2_dp*g12*TrYeadjYe&
          & - 3._dp*(TrYeadjYzYzadjYe+TrYeCYtYtadjYe)  - 12._dp*TrYsCYsYdadjYd &
          & + (3.6_dp*g12+ 12._dp*g22 - 3._dp*(TrCYtYt+TrYeadjYe)              &
          &    - 12._dp * AbsL1sq) * AbsL1sq - 6._dp*TrYzadjYzYdadjYd          &
          & ) * Yd                                                             &
          & + ( 0.8_dp*g12 + 6._dp*g22 - 3._dp*TrYeadjYe                       &
          &   - 9._dp*(TrYdadjYd+AbsL1sq)  ) * YdadjYdYd                       &
          & + ( 0.8_dp*g12  - 3._dp * (TrYuadjYu + AbsL2sq) ) * YdadjYuYu      &
          & + ( 32._dp*g12/15._dp + 80._dp*g32/3._dp -4._dp*TrCYsYs) * YsCYsYd &
          & + ( 0.4_dp*g12  + 6._dp*g22 - 2._dp*TrYzadjYz ) * YzadjYzYd        &
          & - 4._dp * (YdadjYdYdadjYdYd + YdadjYdYsCYsYd)                      &
          & - 2._dp * (YdadjYdYzadjYzYd + YdadjYuYuadjYdYd + YdadjYuYuadjYuYu  &
          &           +YzadjYeYeadjYzYd )   - 16._dp*YsCYsYsCYsYd              &
          & - 8._dp * (YsCYdTYdCYsYd + YsCYzTYzCYsYd)                          &
          & - 6._dp * (YzadjYzYzadjYzYd + YzCYtYtadjYzYd)

   DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 )

  Else
   DYd = oo16pi2* betaYd1
  End If


  !--------------------
  ! Ye
  !--------------------

  betaYe1 = ( 3._dp*(TrYdadjYd+AbsL1sq) + TrYeadjYe    &
          & - 1.8_dp*g12 - 3._dp*g22 ) * Ye            &
          & + 3._dp * (YeadjYeYe + YeadjYzYz + YeCYtYt)

  If (TwoLoopRGE) Then
   betaYe2 = ( 26.1_dp*g14 + 1.8_dp*g12*g22  + 28.5_dp*g24                   &
           & + (16._dp*g32 - 0.4_dp*g12) * TrYdadjYd  + 1.2_dp*g12*TrYeadjYe &
           & - 9._dp*TrYdadjYdYdadjYd - 12._dp*TrYsCYsYdadjYd                &
           & - 3._dp * (TrYdadjYuYuadjYd+TrYeadjYeYeadjYe+TrYeadjYzYzadjYe   &
           &           + TrYeCYtYtadjYe ) - 6._dp*TrYzadjYzYdadjYd           &
           & + (3.6_dp*g12 + 12._dp*g22  - 3._dp*(TrCYtYt+TrYeadjYe)         &
           &   - 9._dp*TrYdadjYd - 12._dp*AbsL1sq) * AbsL1sq                 &
           & ) * Ye                                                          &
           & + ( 6._dp*g22 - 3._dp*TrYeadjYe - 9._dp*(TrYdadjYd+AbsL1sq)     &
           &   ) * YeadjYeYe                                                 &
           & + ( 16._dp*g32 - 0.4_dp*g12  - 3._dp*TrYzadjYz ) * YeadjYzYz    &
           & + ( 3.6_dp*g12 + 12._dp*g22  - 3._dp*(TrYtCYt+AbsL1sq)          &
           &   ) * YeCYtYt                                                   &
           & - 4._dp*YeadjYeYeadjYeYe - 12._dp*YeadjYzYsCYsYz                &
           & - 6._dp * ( YeadjYzYdadjYdYz + YeadjYzYzadjYeYe                 &
           &           + YeadjYzYzadjYzYz + YeCYtYtadjYeYe )                 &
           & - 3._dp*YeCYtTYeCYeYt - 9._dp * (YeCYtTYzCYzYt+YeCYtYtCYtYt)

   DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 )

  Else
   DYe = oo16pi2* betaYe1
  End If

  !--------------------
  ! Yt
  !--------------------

  betaYt1 = (TrCYtYt + AbsL1sq - 1.8_dp*g12  - 7._dp*g22 ) * Yt         &
        & + TYeCYeYt + 3._dp * (TYzCYzYt+YtadjYzYz)+ YtadjYeYe + 6._dp*YtCYtYt

  If (TwoLoopRGE) Then
   betaYt2 = ( 26.1_dp*g14 + 11.4_dp*g12*g22 + 76.5_dp *g24  - TrCYtYtadjYeYe &
         &   - (0.6_dp*g12 + g22)* TrCYtYt - 3._dp * TrCYtYtadjYzYz           &
         &   - 6._dp*TrCYtYtCYtYt - TrYeCYtYtadjYe - 3._dp*TrYtadjYzYzCYt     &
         &   - ( 0.6_dp *g12 + g22 + 6._dp * (TrYdadjYd + AbsL1sq)            &
         &     - 2._dp * TrYeadjYe ) * AbsL1sq                                &
         &   ) * Yt                                                           &
         &   + ( 1.2_dp *g12  - TrYeadjYe - 3._dp*(TrYdadjYd+AbsL1sq)         &
         &     ) * (TYeCYeYt + YtadjYeYe)                                     &
         &   + (16._dp*g32 -0.4_dp*g12 -3._dp*TrYzadjYz )*(TYzCYzYt+YtadjYzYz)&
         &   + ( 7.2_dp*g12 + 24._dp*g22 - 6._dp*(TrCYtYt+AbsL1sq) ) *YtCYtYt &
         &   - 2._dp * (TYeCYeTYeCYeYt + YtadjYeYeadjYeYe)                    &
         &   - 6._dp * (TYzCYdTYdCYzYt + TYzCYzTYzCYzYt                       &
         &             +YtadjYzYdadjYdYz + YtadjYzYzadjYzYz)                  &
         &   - 12._dp *(TYzCYsYsCYzYt + YtadjYzYsCYsYz)                       &
         &   - 3._dp * (YtadjYeYeCYtYt + YtCYtTYeCYeYt)                       &
         &   - 9._dp * (YtadjYzYzCYtYt + YtCYtTYzCYzYt) - 18._dp *YtCYtYtCYtYt

   DYt = oo16pi2*( betaYt1 + oo16pi2 * betaYt2 )

  Else
   DYt = oo16pi2* betaYt1
  End If


  !--------------------
  ! Ys
  !--------------------

  betaYs1 = ( TrCYsYs - 0.8_dp *g12 - 12._dp*g32 ) * Ys              &
          & + 2._dp * ( YdadjYdYs + YsCYdTYd +YsCYzTYz + YzadjYzYs ) &
          & + 8._dp * YsCYsYs

  If (TwoLoopRGE) Then
   betaYs2 = ( 11.2_dp*g14 + 128._dp*g12*g32/15._dp + 320._dp*g34/3._dp       &
         &   - 4._dp * (g12/15._dp + g32/3._dp )* TrCYsYs                     &
         &   - 4._dp * (TrCYsYdadjYdYs+TrCYsYzadjYzYs) - 8._dp * TrCYsYsCYsYs &
         &   ) * Ys                                                           &
         &   + (0.4_dp*g12 + 6._dp*(g22-TrYdadjYd-AbsL1sq) - 2._dp*TrYeadjYe  &
         &     ) * (YdadjYdYs+YsCYdTYd)                                       &
         &   + (64._dp*g12/15._dp + 160._dp*g32/3._dp - 8._dp*TrCYsYs)*YsCYsYs&
         &   + (0.4_dp*g12 +6._dp*g22 -2._dp*TrYzadjYz) *(YsCYzTYz+YzadjYzYs) &
         &   - 2._dp * (YdadjYdYdadjYdYs + YdadjYuYuadjYdYs + YsCYdTYdCYdTYd  &
         &             +YsCYdTYuCYuTYd + YsCYzTYeCYeTYz + YzadjYeYeadjYzYs)   &
         &   - 8._dp * ( YsCYdTYdCYsYs + YsCYsYdadjYdYs                       &
         &             + YsCYsYzadjYzYs + YsCYzTYzCYsYs )                     &
         &   - 6._dp * ( YsCYzTYzCYzTYz + YsCYzYtCYtTYz                       &
         &             + YzadjYzYzadjYzYs + YzCYtYtadjYzYs)                   &
         &   - 32._dp*YsCYsYsCYsYs

   DYs = oo16pi2*( betaYs1 + oo16pi2 * betaYs2 )

  Else
   DYs = oo16pi2* betaYs1
  End If


  !--------------------
  ! Yz
  !--------------------
  betaYz1 = (TrYzadjYz -7._dp*g12/15._dp -3._dp*g22 -16._dp*g32/3._dp ) * Yz &
        & + 2._dp*YdadjYdYz + 4._dp*YsCYsYz + YzadjYeYe + 5._dp*YzadjYzYz    &
        & + 3._dp*YzCYtYt


  If (TwoLoopRGE) Then
   betaYz2 = ( 581._dp*g14/90._dp + g12*g22 + 28.5_dp*g24  + 320._dp*g34/9._dp &
         &   + 8._dp*g32*(g12/9._dp + g22) + 0.4_dp*g12*TrYzadjYz              &
         &   - TrYzadjYeYeadjYz- 2._dp*TrYdadjYdYzadjYz- 4._dp*TrYsCYsYzadjYz  &
         &   - 5._dp*TrYzadjYzYzadjYz - 3._dp*TrYzCYtYtadjYz                   &
         &   ) * Yz                                                            &
         &   + (0.4_dp*g12 + 6._dp*(g22 - TrYdadjYd -AbsL1sq) -2._dp*TrYeadjYe &
         &      ) * YdadjYdYz                                                  &
         &   + (32._dp*g12/15._dp + 80._dp*g32/3._dp - 4._dp*TrCYsYs) *YsCYsYz &
         &   + (1.2_dp*g12 - 3._dp*(TrYdadjYd+AbsL1sq) - TrYeadjYe) *YzadjYeYe &
         &   + (3.6_dp*g12 + 12._dp*g22 - 3._dp*(TrCYtYt+AbsL1sq) ) * YzCYtYt  &
         &   + (6._dp*g22 + 16._dp*g32  - 5._dp*TrYzadjYz) * YzadjYzYz         &
         &   - 2._dp * ( YdadjYdYdadjYdYz + YdadjYuYuadjYdYz                   &
         &             + YzadjYeYeadjYeYe + YzadjYeYeadjYzYz )                 &
         &   - 9._dp * ( YzCYtTYzCYzYt + YzCYtYtCYtYt) - 3._dp*YzCYtTYeCYeYt   &
         &   - 8._dp * ( YsCYdTYdCYsYz +YsCYzTYzCYsYz) - 16._dp*YsCYsYsCYsYz   &
         &   - 6._dp * ( YzadjYzYdadjYdYz + YzCYtYtadjYzYz)                    &
         &   - 12._dp * (YzadjYzYsCYsYz + YzadjYzYzadjYzYz)

   DYz = oo16pi2*( betaYz1 + oo16pi2 * betaYz2 )

  Else
   DYz = oo16pi2* betaYz1
  End If


  !--------------------
  ! L1
  !--------------------

  betaL11 = -1.8_dp*g12 - 7._dp * g22 + TrCYtYt + 6._dp * TrYdadjYd &
          & + 2._dp*TrYeadjYe + 7._dp*AbsL1sq

  If (TwoLoopRGE) Then
   betaL12 = 26.1_dp*g14 + 11.4_dp*g12*g22 + 76.5_dp*g24                  &
         &   - (0.6_dp*g12+g22 + 6._dp* AbsL1sq) *TrCYtYt                 &
         &   + (32._dp*g32 - 0.8_dp*g12 - 24._dp* AbsL1sq) * TrYdadjYd    &
         &   + (2.4_dp*g12 - 8._dp* AbsL1sq) * TrYeadjYe                  &
         &   + (6.6_dp*g12 + 23._dp*g22 - 30._dp* AbsL1sq) * AbsL1sq      &
         &   - 8._dp*TrYeCYtYtadjYe  - 12._dp * TrYzadjYzYdadjYd          &
         &   - 6._dp * (TrCYtYtCYtYt + TrYdadjYuYuadjYd + TrCYtYtadjYzYz  &
         &             +TrYeadjYeYeadjYe + TrYeadjYzYzadjYe)              &
         &   - 18._dp*TrYdadjYdYdadjYd - 24._dp * TrYsCYsYdadjYd

   DL1 = oo16pi2 * L1 * ( betaL11 + oo16pi2 * betaL12 )

  Else
   DL1 = oo16pi2 * L1 * betaL11
  End If


  !--------------------
  ! L2
  !--------------------

  betaL21 = -1.8_dp*g12 - 7._dp * g22 + 6._dp * TrYuadjYu + 7._dp*AbsL2sq

  If (TwoLoopRGE) Then
   betaL22 = 26.1_dp*g14 + 11.4_dp*g12*g22 + 76.5_dp*g24             &
         & + (1.6_dp*g12 + 32._dp*g32 - 24._dp*AbsL2sq ) * TrYuadjYu &
         & + (6.6_dp*g12 + 23._dp*g22 - 30._dp*AbsL2sq ) * AbsL2sq   &
         & - 6._dp * TrYuadjYdYdadjYu - 18._dp*TrYuadjYuYuadjYu

   DL2 = oo16pi2 * L2 * ( betaL21 + oo16pi2 * betaL22 )

  Else
   DL2 = oo16pi2 * L2 * betaL21
  End If


  !--------------------
  ! MTM
  !--------------------


  betaMTM1 = -2.4_dp*g12 - 8._dp*g22 + TrCYtYt + AbsL1sq + AbsL2sq

  If (TwoLoopRGE) Then
   betaMTM2 = 35.52_dp*g14 + 19.2_dp*g12*g22 + 96._dp*g24           &
          & - (0.6_dp*g12 + g22) * TrCYtYt - 2._dp * TrCYtYtadjYeYe &
          & - 6._dp*TrCYtYtadjYzYz - 6._dp*TrCYtYtCYtYt             &
          & - (0.6_dp*g12 + g22 + 6._dp*(TrYdadjYd+AbsL1sq)         &
          &   + 2._dp*TrYeadjYe ) * AbsL1sq                         &
          & - (0.6_dp*g12 + g22 + 6._dp*(TrYuadjYu+AbsL2sq) ) * AbsL2sq

   DMTM = oo16pi2 * MTM * ( betaMTM1 + oo16pi2 * betaMTM2 )

  Else
   DMTM = oo16pi2 * MTM * betaMTM1
  End If

  Call ParametersToG117(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYt,DYs,DYz,DL1,DL2,DMTM,f)

  Iname = Iname - 1

 End Subroutine rge117

 Subroutine rge353(len, T, GY, F)
 Implicit None
  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)
  Integer :: i1,i2
  Real(dp) :: q
  Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,         &
  & Dg3,mHd2,betamHd21,betamHd22,DmHd2,mHu2,betamHu21,betamHu22,DmHu2
  Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          &
  & ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       &
  & ,DYe(3,3),adjYe(3,3),Yt(3,3),betaYt1(3,3),betaYt2(3,3),DYt(3,3),adjYt(3,3)             &
  & ,Ys(3,3),betaYs1(3,3),betaYs2(3,3),DYs(3,3),adjYs(3,3),Yz(3,3),betaYz1(3,3)            &
  & ,betaYz2(3,3),DYz(3,3),adjYz(3,3),L1,betaL11,betaL12,DL1,L2,betaL21,betaL22,           &
  & DL2,MTM,betaMTM1,betaMTM2,DMTM,mue,betamue1,betamue2,Dmue,MZM,betaMZM1,betaMZM2,       &
  & DMZM,MSM,betaMSM1,betaMSM2,DMSM,AYu(3,3),betaAYu1(3,3),betaAYu2(3,3),DAYu(3,3)         &
  & ,adjAYu(3,3),AYd(3,3),betaAYd1(3,3),betaAYd2(3,3),DAYd(3,3),adjAYd(3,3),               &
  & AYe(3,3),betaAYe1(3,3),betaAYe2(3,3),DAYe(3,3),adjAYe(3,3),AYt(3,3),betaAYt1(3,3)      &
  & ,betaAYt2(3,3),DAYt(3,3),adjAYt(3,3),AYs(3,3),betaAYs1(3,3),betaAYs2(3,3)              &
  & ,DAYs(3,3),adjAYs(3,3),AYz(3,3),betaAYz1(3,3),betaAYz2(3,3),DAYz(3,3),adjAYz(3,3)      &
  & ,AL1,betaAL11,betaAL12,DAL1,AL2,betaAL21,betaAL22,DAL2,Amue,betaAmue1,betaAmue2,       &
  & DAmue,AMTM,betaAMTM1,betaAMTM2,DAMTM,AMZM,betaAMZM1,betaAMZM2,DAMZM,AMSM,              &
  & betaAMSM1,betaAMSM2,DAMSM,mq2(3,3),betamq21(3,3),betamq22(3,3),Dmq2(3,3)               &
  & ,adjmq2(3,3),ml2(3,3),betaml21(3,3),betaml22(3,3),Dml2(3,3),adjml2(3,3),               &
  & md2(3,3),betamd21(3,3),betamd22(3,3),Dmd2(3,3),adjmd2(3,3),mu2(3,3),betamu21(3,3)      &
  & ,betamu22(3,3),Dmu2(3,3),adjmu2(3,3),me2(3,3),betame21(3,3),betame22(3,3)              &
  & ,Dme2(3,3),adjme2(3,3),mt2,betamt21,betamt22,Dmt2,mtb2,betamtb21,betamtb22,            &
  & Dmtb2,ms2,betams21,betams22,Dms2,msb2,betamsb21,betamsb22,Dmsb2,mz2,betamz21,          &
  & betamz22,Dmz2,mzb2,betamzb21,betamzb22,Dmzb2,MassB,betaMassB1,betaMassB2,              &
  & DMassB,MassWB,betaMassWB1,betaMassWB2,DMassWB,MassG,betaMassG1,betaMassG2,DMassG
  Complex(dp) :: Tr1(3),Tr2(3),Tr3(3)
  Complex(dp) :: md2Ys(3,3),md2CYd(3,3),md2CYs(3,3),md2CYz(3,3),me2CYe(3,3)             &
  & ,ml2Yt(3,3),ml2adjYe(3,3),ml2adjYz(3,3),ml2CYt(3,3),mq2adjYd(3,3),mq2adjYu(3,3)        &
  & ,mu2CYu(3,3),YdadjYd(3,3),YeadjYe(3,3),Ysmd2(3,3),YsCYs(3,3),Ytml2(3,3),               &
  & YtCYt(3,3),YuadjYu(3,3),YzadjYz(3,3),AYdadjYd(3,3),AYdadjAYd(3,3),AYeadjYe(3,3)        &
  & ,AYeadjAYe(3,3),AYsCAYs(3,3),AYtCAYt(3,3),AYuadjYu(3,3),AYuadjAYu(3,3),AYzadjYz(3,3)   &
  & ,AYzadjAYz(3,3),adjYdmd2(3,3),adjYdYd(3,3),adjYdYs(3,3),adjYdYz(3,3),adjYdAYd(3,3)     &
  & ,adjYdAYs(3,3),adjYdAYz(3,3),adjYeme2(3,3),adjYeYe(3,3),adjYeAYe(3,3),adjYumu2(3,3)    &
  & ,adjYuYu(3,3),adjYuAYu(3,3),adjYzmd2(3,3),adjYzYd(3,3),adjYzYs(3,3),adjYzYz(3,3)       &
  & ,adjYzAYd(3,3),adjYzAYs(3,3),adjYzAYz(3,3),CYdmq2(3,3),CYdTYd(3,3),CYdTAYd(3,3)        &
  & ,CYeml2(3,3),CYeYt(3,3),CYeAYt(3,3),CYsmd2(3,3),CYsYd(3,3),CYsYs(3,3),CYsYz(3,3)       &
  & ,CYsAYd(3,3),CYsAYs(3,3),CYsAYz(3,3),CYtml2(3,3),CYtYt(3,3),CYtAYt(3,3),               &
  & CYumq2(3,3),CYzml2(3,3),CYzYt(3,3),CYzAYt(3,3),CYzTYz(3,3),CYzTAYz(3,3),               &
  & CAYsAYs(3,3),CAYtAYt(3,3),TYdCYd(3,3),TYeCYe(3,3),TYuCYu(3,3),TYzCYz(3,3)              &
  & ,TAYdCAYd(3,3),TAYeCAYe(3,3),TAYuCAYu(3,3),TAYzCAYz(3,3),md2YdadjYd(3,3)               &
  & ,md2YsCYs(3,3),md2YzadjYz(3,3),me2YeadjYe(3,3),ml2YtCYt(3,3),ml2TYeCYe(3,3)            &
  & ,ml2TYzCYz(3,3),mq2TYdCYd(3,3),mq2TYuCYu(3,3),mu2YuadjYu(3,3),Ydmq2adjYd(3,3)          &
  & ,YdadjYdmd2(3,3),YdadjYdYd(3,3),YdadjYdYs(3,3),YdadjYdYz(3,3),YdadjYdAYd(3,3)          &
  & ,YdadjYdAYs(3,3),YdadjYdAYz(3,3),YdadjYuYu(3,3),YdadjYuAYu(3,3),Yeml2adjYe(3,3)        &
  & ,YeadjYeme2(3,3),YeadjYeYe(3,3),YeadjYeAYe(3,3),YeadjYzYz(3,3),YeadjYzAYz(3,3)         &
  & ,YeCYtYt(3,3),YeCYtAYt(3,3),Ysmd2CYs(3,3),YsCYdTYd(3,3),YsCYdTAYd(3,3),YsCYsmd2(3,3)   &
  & ,YsCYsYd(3,3),YsCYsYs(3,3),YsCYsYz(3,3),YsCYsAYd(3,3),YsCYsAYs(3,3),YsCYsAYz(3,3)      &
  & ,YsCYzTYz(3,3),YsCYzTAYz(3,3),Ytml2CYt(3,3),YtadjYeYe(3,3),YtadjYeAYe(3,3)             &
  & ,YtadjYzYz(3,3),YtadjYzAYz(3,3),YtCYtml2(3,3),YtCYtYt(3,3),YtCYtAYt(3,3)               &
  & ,Yumq2adjYu(3,3),YuadjYdYd(3,3),YuadjYdAYd(3,3),YuadjYumu2(3,3),YuadjYuYu(3,3)         &
  & ,YuadjYuAYu(3,3),Yzml2adjYz(3,3),YzadjYeYe(3,3),YzadjYeAYe(3,3),YzadjYzmd2(3,3)        &
  & ,YzadjYzYd(3,3),YzadjYzYs(3,3),YzadjYzYz(3,3),YzadjYzAYd(3,3),YzadjYzAYs(3,3)          &
  & ,YzadjYzAYz(3,3),YzCYtYt(3,3),YzCYtAYt(3,3),AYdadjYdYd(3,3),AYdadjYdYs(3,3)            &
  & ,AYdadjYdYz(3,3),AYdadjYuYu(3,3),AYeadjYeYe(3,3),AYeadjYzYz(3,3),AYeCYtYt(3,3)         &
  & ,AYsCYdTYd(3,3),AYsCYsYd(3,3),AYsCYsYs(3,3),AYsCYsYz(3,3),AYsCYzTYz(3,3)               &
  & ,AYtadjYeYe(3,3),AYtadjYzYz(3,3),AYtCYtYt(3,3),AYuadjYdYd(3,3),AYuadjYuYu(3,3)         &
  & ,AYzadjYeYe(3,3),AYzadjYzYd(3,3),AYzadjYzYs(3,3),AYzadjYzYz(3,3),AYzCYtYt(3,3)         &
  & ,CYsmd2Ys(3,3),CYsYsmd2(3,3),CYtml2Yt(3,3),CYtYtml2(3,3),TYdmd2CYd(3,3),               &
  & TYdCYdmq2(3,3),TYeme2CYe(3,3),TYeCYeml2(3,3),TYeCYeYt(3,3),TYeCYeAYt(3,3)              &
  & ,TYumu2CYu(3,3),TYuCYumq2(3,3),TYzmd2CYz(3,3),TYzCYzml2(3,3),TYzCYzYt(3,3)             &
  & ,TYzCYzAYt(3,3),TAYeCYeYt(3,3),TAYzCYzYt(3,3)

  Complex(dp) :: md2Yz(3,3),me2Ye(3,3),YdadjYu(3,3),YdadjAYd(3,3),YdadjAYu(3,3)         &
  & ,Yeml2(3,3),YeadjYz(3,3),YeadjAYe(3,3),YeadjAYz(3,3),YeCYt(3,3),YeCAYt(3,3)            &
  & ,YsCYd(3,3),YsCYz(3,3),YsCAYd(3,3),YsCAYs(3,3),YsCAYz(3,3),YtadjYe(3,3),               &
  & YtadjYz(3,3),YtadjAYe(3,3),YtadjAYz(3,3),YtCAYt(3,3),YuadjYd(3,3),YuadjAYd(3,3)        &
  & ,YuadjAYu(3,3),Yzml2(3,3),YzadjYe(3,3),YzadjAYe(3,3),YzadjAYz(3,3),YzCYt(3,3)          &
  & ,YzCAYt(3,3),AYdadjYu(3,3),AYdadjAYu(3,3),AYeadjYz(3,3),AYeadjAYz(3,3),AYeCYt(3,3)     &
  & ,AYeCAYt(3,3),AYsCYd(3,3),AYsCYs(3,3),AYsCYz(3,3),AYsCAYd(3,3),AYsCAYz(3,3)            &
  & ,AYtadjYe(3,3),AYtadjYz(3,3),AYtadjAYe(3,3),AYtadjAYz(3,3),AYtCYt(3,3),AYuadjYd(3,3)   &
  & ,AYuadjAYd(3,3),AYzadjYe(3,3),AYzadjAYe(3,3),AYzCYt(3,3),AYzCAYt(3,3),adjAYdYs(3,3)    &
  & ,adjAYdAYs(3,3),adjAYeYe(3,3),adjAYeAYe(3,3),adjAYzYs(3,3),adjAYzYz(3,3)               &
  & ,adjAYzAYs(3,3),adjAYzAYz(3,3),CYeTYz(3,3),CYeTAYz(3,3),CYtTYz(3,3),CYtTAYz(3,3)       &
  & ,CYuTYd(3,3),CYuTAYd(3,3),CAYsYs(3,3),CAYtYt(3,3),TYdCYs(3,3),TYdCYz(3,3)              &
  & ,TYdCAYd(3,3),TYdCAYs(3,3),TYdCAYz(3,3),TYeCAYe(3,3),TYuCAYu(3,3),TYzCYd(3,3)          &
  & ,TYzCYs(3,3),TYzCAYd(3,3),TYzCAYs(3,3),TYzCAYz(3,3),TAYdCYd(3,3),TAYdCYs(3,3)          &
  & ,TAYdCYz(3,3),TAYdCAYs(3,3),TAYdCAYz(3,3),TAYeCYe(3,3),TAYuCYu(3,3),TAYzCYd(3,3)       &
  & ,TAYzCYs(3,3),TAYzCYz(3,3),TAYzCAYd(3,3),TAYzCAYs(3,3),md2YdadjYu(3,3),md2YsCYd(3,3)   &
  & ,md2YsCYz(3,3),md2YzadjYe(3,3),md2YzCYt(3,3),md2CYsYs(3,3),me2YeadjYz(3,3)             &
  & ,me2YeCYt(3,3),ml2YtadjYe(3,3),ml2YtadjYz(3,3),ml2adjYeYe(3,3),ml2adjYzYs(3,3)         &
  & ,ml2adjYzYz(3,3),ml2CYtYt(3,3),ml2TYzCYd(3,3),ml2TYzCYs(3,3),mq2adjYdYs(3,3)           &
  & ,mq2TYdCYs(3,3),mq2TYdCYz(3,3),mu2YuadjYd(3,3),Ydmq2adjYu(3,3),YdadjYumu2(3,3)         &
  & ,YdadjAYdAYs(3,3),Yeml2adjYz(3,3),Yeml2CYt(3,3),YeadjYzmd2(3,3),YeadjYzYd(3,3)         &
  & ,YeadjYzYs(3,3),YeadjYzAYd(3,3),YeadjYzAYs(3,3),YeCYtml2(3,3),Ysmd2CYd(3,3)            &
  & ,Ysmd2CYz(3,3),YsCYdmq2(3,3),YsCYzml2(3,3),YsCYzYt(3,3),YsCYzAYt(3,3),YsCAYsAYs(3,3)   &
  & ,Ytml2adjYe(3,3),Ytml2adjYz(3,3),YtadjYeme2(3,3),YtadjYzmd2(3,3),YtadjYzYd(3,3)        &
  & ,YtadjYzYs(3,3),YtadjYzAYd(3,3),YtadjYzAYs(3,3),YtadjAYeAYe(3,3),YtadjAYzAYz(3,3)      &
  & ,YtCYtTYz(3,3),YtCYtTAYz(3,3),YtCAYtAYt(3,3),Yumq2adjYd(3,3),YuadjYdmd2(3,3)           &
  & ,YuadjYdYs(3,3),YuadjYdYz(3,3),YuadjYdAYs(3,3),YuadjYdAYz(3,3),Yzml2adjYe(3,3)         &
  & ,Yzml2CYt(3,3),YzadjYeme2(3,3),YzadjAYzAYs(3,3),YzCYtml2(3,3),AYdadjAYdYs(3,3)         &
  & ,AYeadjYzYd(3,3),AYeadjYzYs(3,3),AYsCYzYt(3,3),AYsCAYsYs(3,3),AYtadjYzYd(3,3)          &
  & ,AYtadjYzYs(3,3),AYtadjAYeYe(3,3),AYtadjAYzYz(3,3),AYtCYtTYz(3,3),AYtCAYtYt(3,3)       &
  & ,AYuadjYdYs(3,3),AYuadjYdYz(3,3),AYzadjAYzYs(3,3),adjYdmd2Ys(3,3),adjYdYdadjYd(3,3)    &
  & ,adjYdYdadjYu(3,3),adjYdYdadjAYd(3,3),adjYdYdadjAYu(3,3),adjYdYsmd2(3,3)               &
  & ,adjYdYsCYs(3,3),adjYdYzadjYz(3,3),adjYdAYdadjYd(3,3),adjYdAYdadjYu(3,3)               &
  & ,adjYdAYdadjAYd(3,3),adjYdAYdadjAYu(3,3),adjYdAYsCYs(3,3),adjYdAYsCAYs(3,3)            &
  & ,adjYdAYzadjYz(3,3),adjYdAYzadjAYz(3,3),adjYeme2Ye(3,3),adjYeYeml2(3,3),               &
  & adjYeYeadjYe(3,3),adjYeYeadjYz(3,3),adjYeYeadjAYe(3,3),adjYeYeadjAYz(3,3)              &
  & ,adjYeYeCYt(3,3),adjYeYeCAYt(3,3),adjYeAYeadjYe(3,3),adjYeAYeadjYz(3,3),               &
  & adjYeAYeadjAYe(3,3),adjYeAYeadjAYz(3,3),adjYeAYeCYt(3,3),adjYeAYeCAYt(3,3)             &
  & ,adjYuYuadjYd(3,3),adjYuYuadjYu(3,3),adjYuYuadjAYd(3,3),adjYuYuadjAYu(3,3)             &
  & ,adjYuAYuadjYd(3,3),adjYuAYuadjYu(3,3),adjYuAYuadjAYd(3,3),adjYuAYuadjAYu(3,3)         &
  & ,adjYzmd2Ys(3,3),adjYzmd2Yz(3,3),adjYzYdadjYd(3,3),adjYzYsmd2(3,3),adjYzYsCYs(3,3)     &
  & ,adjYzYzml2(3,3),adjYzYzadjYe(3,3),adjYzYzadjYz(3,3),adjYzYzadjAYe(3,3),               &
  & adjYzYzadjAYz(3,3),adjYzYzCYt(3,3),adjYzYzCAYt(3,3),adjYzAYdadjYd(3,3),adjYzAYdadjAYd(3,3)&
  & ,adjYzAYsCYs(3,3),adjYzAYsCAYs(3,3),adjYzAYzadjYe(3,3),adjYzAYzadjYz(3,3)              &
  & ,adjYzAYzadjAYe(3,3),adjYzAYzadjAYz(3,3),adjYzAYzCYt(3,3),adjYzAYzCAYt(3,3)            &
  & ,adjAYdYdadjYd(3,3),adjAYdYdadjYu(3,3),adjAYdAYdadjYd(3,3),adjAYdAYdadjYu(3,3)         &
  & ,adjAYdAYsCYs(3,3),adjAYdAYzadjYz(3,3),adjAYeYeadjYe(3,3),adjAYeYeadjYz(3,3)           &
  & ,adjAYeYeCYt(3,3),adjAYeAYeadjYe(3,3),adjAYeAYeadjYz(3,3),adjAYeAYeCYt(3,3)            &
  & ,adjAYuYuadjYd(3,3),adjAYuYuadjYu(3,3),adjAYuAYuadjYd(3,3),adjAYuAYuadjYu(3,3)         &
  & ,adjAYzYzadjYe(3,3),adjAYzYzadjYz(3,3),adjAYzYzCYt(3,3),adjAYzAYdadjYd(3,3)            &
  & ,adjAYzAYsCYs(3,3),adjAYzAYzadjYe(3,3),adjAYzAYzadjYz(3,3),adjAYzAYzCYt(3,3)           &
  & ,CYdTYdCYd(3,3),CYdTYdCYs(3,3),CYdTYdCYz(3,3),CYdTYdCAYd(3,3),CYdTYdCAYs(3,3)          &
  & ,CYdTYdCAYz(3,3),CYdTAYdCAYd(3,3),CYdTAYdCAYs(3,3),CYdTAYdCAYz(3,3),CYeTYeCYe(3,3)     &
  & ,CYeTYeCAYe(3,3),CYeTAYeCAYe(3,3),CYsYdadjYd(3,3),CYsYsCYd(3,3),CYsYsCYs(3,3)          &
  & ,CYsYsCYz(3,3),CYsYsCAYd(3,3),CYsYsCAYs(3,3),CYsYsCAYz(3,3),CYsYzadjYz(3,3)            &
  & ,CYsAYdadjYd(3,3),CYsAYdadjAYd(3,3),CYsAYsCYs(3,3),CYsAYsCAYd(3,3),CYsAYsCAYs(3,3)     &
  & ,CYsAYsCAYz(3,3),CYsAYzadjYz(3,3),CYsAYzadjAYz(3,3),CYtYtadjYe(3,3),CYtYtadjYz(3,3)    &
  & ,CYtYtadjAYe(3,3),CYtYtadjAYz(3,3),CYtYtCYt(3,3),CYtYtCAYt(3,3),CYtAYtadjYe(3,3)       &
  & ,CYtAYtadjYz(3,3),CYtAYtadjAYe(3,3),CYtAYtadjAYz(3,3),CYtAYtCYt(3,3),CYtAYtCAYt(3,3)   &
  & ,CYuTYuCYu(3,3),CYuTYuCAYu(3,3),CYuTAYuCAYu(3,3),CYzTYzCYd(3,3),CYzTYzCYs(3,3)         &
  & ,CYzTYzCYz(3,3),CYzTYzCAYd(3,3),CYzTYzCAYs(3,3),CYzTYzCAYz(3,3),CYzTAYzCAYd(3,3)       &
  & ,CYzTAYzCAYs(3,3),CYzTAYzCAYz(3,3),CAYdTYdCYd(3,3),CAYdTYdCYs(3,3),CAYdTYdCYz(3,3)     &
  & ,CAYdTAYdCYd(3,3),CAYdTAYdCYs(3,3),CAYdTAYdCYz(3,3),CAYeTYeCYe(3,3),CAYeTAYeCYe(3,3)   &
  & ,CAYsYsCYd(3,3),CAYsYsCYs(3,3),CAYsYsCYz(3,3),CAYsAYdadjYd(3,3),CAYsAYsCYd(3,3)        &
  & ,CAYsAYsCYs(3,3),CAYsAYsCYz(3,3),CAYsAYzadjYz(3,3),CAYtYtadjYe(3,3),CAYtYtadjYz(3,3)   &
  & ,CAYtYtCYt(3,3),CAYtAYtadjYe(3,3),CAYtAYtadjYz(3,3),CAYtAYtCYt(3,3),CAYuTYuCYu(3,3)    &
  & ,CAYuTAYuCYu(3,3),CAYzTYzCYd(3,3),CAYzTYzCYs(3,3),CAYzTYzCYz(3,3),CAYzTAYzCYd(3,3)     &
  & ,CAYzTAYzCYs(3,3),CAYzTAYzCYz(3,3),TYdmd2CYs(3,3),TYdmd2CYz(3,3),TYdCYdTYd(3,3)        &
  & ,TYdCYdTAYd(3,3),TYdCYsmd2(3,3),TYdCYsYd(3,3),TYdCYsYs(3,3),TYdCYsYz(3,3)              &
  & ,TYdCYsAYd(3,3),TYdCYsAYs(3,3),TYdCYsAYz(3,3),TYdCYzml2(3,3),TYdCYzYt(3,3)             &
  & ,TYdCYzAYt(3,3),TYeCYeTYz(3,3),TYeCYeTAYz(3,3),TYuCYuTYd(3,3),TYuCYuTAYd(3,3)          &
  & ,TYzmd2CYd(3,3),TYzmd2CYs(3,3),TYzCYdmq2(3,3),TYzCYsmd2(3,3),TYzCYsYd(3,3)             &
  & ,TYzCYsYs(3,3),TYzCYsYz(3,3),TYzCYsAYd(3,3),TYzCYsAYs(3,3),TYzCYsAYz(3,3)              &
  & ,TYzCYzTYz(3,3),TYzCYzTAYz(3,3),TAYdCYdTYd(3,3),TAYdCYsYd(3,3),TAYdCYsYs(3,3)          &
  & ,TAYdCYsYz(3,3),TAYdCYzYt(3,3),TAYeCYeTYz(3,3),TAYuCYuTYd(3,3),TAYzCYsYd(3,3)          &
  & ,TAYzCYsYs(3,3),TAYzCYsYz(3,3),TAYzCYzTYz(3,3),md2YsCYsYs(3,3),md2CYdTYdCYd(3,3)       &
  & ,md2CYdTYdCYs(3,3),md2CYdTYdCYz(3,3),md2CYsYdadjYd(3,3),md2CYsYsCYd(3,3)               &
  & ,md2CYsYsCYs(3,3),md2CYsYsCYz(3,3),md2CYsYzadjYz(3,3),md2CYzTYzCYd(3,3),               &
  & md2CYzTYzCYs(3,3),md2CYzTYzCYz(3,3),me2CYeTYeCYe(3,3),ml2YtCYtYt(3,3),ml2adjYeYeadjYe(3,3)&
  & ,ml2adjYeYeadjYz(3,3),ml2adjYeYeCYt(3,3),ml2adjYzYdadjYd(3,3),ml2adjYzYsCYs(3,3)       &
  & ,ml2adjYzYzadjYe(3,3),ml2adjYzYzadjYz(3,3),ml2adjYzYzCYt(3,3),ml2CYtYtadjYe(3,3)       &
  & ,ml2CYtYtadjYz(3,3),ml2CYtYtCYt(3,3),mq2adjYdYdadjYd(3,3),mq2adjYdYdadjYu(3,3)         &
  & ,mq2adjYdYsCYs(3,3),mq2adjYdYzadjYz(3,3),mq2adjYuYuadjYd(3,3),mq2adjYuYuadjYu(3,3)     &
  & ,mu2CYuTYuCYu(3,3),Ydmq2adjYdYs(3,3),YdadjYdmd2Ys(3,3),YdadjYdYdadjYd(3,3)             &
  & ,YdadjYdYsmd2(3,3),YdadjYdYsCYs(3,3),YdadjYdYzadjYz(3,3),YdadjYdAYdadjYd(3,3)          &
  & ,YdadjYdAYdadjAYd(3,3),YdadjYdAYsCYs(3,3),YdadjYdAYsCAYs(3,3),YdadjYdAYzadjYz(3,3)     &
  & ,YdadjYdAYzadjAYz(3,3),YdadjYuYuadjYd(3,3),YdadjYuAYuadjYd(3,3),YdadjYuAYuadjAYd(3,3)  &
  & ,YdadjAYdAYdadjYd(3,3),YdadjAYdAYsCYs(3,3),YdadjAYdAYzadjYz(3,3),YdadjAYuAYuadjYd(3,3) &
  & ,YeadjYeYeadjYe(3,3),YeadjYeAYeadjYe(3,3),YeadjYeAYeadjAYe(3,3),YeadjYzYzadjYe(3,3)    &
  & ,YeadjYzAYzadjYe(3,3),YeadjYzAYzadjAYe(3,3),YeadjAYeAYeadjYe(3,3),YeadjAYzAYzadjYe(3,3)&
  & ,YeCYtYtadjYe(3,3),YeCYtAYtadjYe(3,3),YeCYtAYtadjAYe(3,3),YeCAYtAYtadjYe(3,3)          &
  & ,Ysmd2CYsYs(3,3),YsCYdTYdCYs(3,3),YsCYdTAYdCAYs(3,3),YsCYsmd2Ys(3,3),YsCYsYdadjYd(3,3) &
  & ,YsCYsYsmd2(3,3),YsCYsYsCYs(3,3),YsCYsYzadjYz(3,3),YsCYsAYdadjYd(3,3),YsCYsAYdadjAYd(3,3)&
  & ,YsCYsAYsCYs(3,3),YsCYsAYsCAYs(3,3),YsCYsAYzadjYz(3,3),YsCYsAYzadjAYz(3,3)             &
  & ,YsCYzTYzCYs(3,3),YsCYzTAYzCAYs(3,3),YsCAYdTAYdCYs(3,3),YsCAYsAYdadjYd(3,3)            &
  & ,YsCAYsAYsCYs(3,3),YsCAYsAYzadjYz(3,3),YsCAYzTAYzCYs(3,3),Ytml2adjYeYe(3,3)            &
  & ,Ytml2adjYzYz(3,3),Ytml2CYtYt(3,3),YtadjYeme2Ye(3,3),YtadjYeYeml2(3,3),YtadjYeYeCYt(3,3)&
  & ,YtadjYeAYeCYt(3,3),YtadjYeAYeCAYt(3,3),YtadjYzmd2Yz(3,3),YtadjYzYzml2(3,3)            &
  & ,YtadjYzYzCYt(3,3),YtadjYzAYzCYt(3,3),YtadjYzAYzCAYt(3,3),YtadjAYeAYeCYt(3,3)          &
  & ,YtadjAYzAYzCYt(3,3),YtCYtml2Yt(3,3),YtCYtYtml2(3,3),YtCYtYtCYt(3,3),YtCYtAYtCYt(3,3)  &
  & ,YtCYtAYtCAYt(3,3),YtCAYtAYtCYt(3,3),YuadjYdYdadjYu(3,3),YuadjYdAYdadjYu(3,3)          &
  & ,YuadjYdAYdadjAYu(3,3),YuadjYuYuadjYu(3,3),YuadjYuAYuadjYu(3,3),YuadjYuAYuadjAYu(3,3)  &
  & ,YuadjAYdAYdadjYu(3,3),YuadjAYuAYuadjYu(3,3),Yzml2adjYzYs(3,3),YzadjYeYeadjYz(3,3)     &
  & ,YzadjYeAYeadjYz(3,3),YzadjYeAYeadjAYz(3,3),YzadjYzmd2Ys(3,3),YzadjYzYdadjYd(3,3)      &
  & ,YzadjYzYsmd2(3,3),YzadjYzYsCYs(3,3),YzadjYzYzadjYz(3,3),YzadjYzAYdadjYd(3,3)          &
  & ,YzadjYzAYdadjAYd(3,3),YzadjYzAYsCYs(3,3),YzadjYzAYsCAYs(3,3),YzadjYzAYzadjYz(3,3)     &
  & ,YzadjYzAYzadjAYz(3,3),YzadjAYeAYeadjYz(3,3),YzadjAYzAYdadjYd(3,3),YzadjAYzAYsCYs(3,3) &
  & ,YzadjAYzAYzadjYz(3,3),YzCYtYtadjYz(3,3),YzCYtAYtadjYz(3,3),YzCYtAYtadjAYz(3,3)        &
  & ,YzCAYtAYtadjYz(3,3),AYdadjYdYdadjAYd(3,3),AYdadjYuYuadjAYd(3,3),AYdadjAYdYdadjYd(3,3) &
  & ,AYdadjAYuYuadjYd(3,3),AYeadjYeYeadjAYe(3,3),AYeadjYzYzadjAYe(3,3),AYeadjAYeYeadjYe(3,3)&
  & ,AYeadjAYzYzadjYe(3,3),AYeCYtYtadjAYe(3,3),AYeCAYtYtadjYe(3,3),AYsCYdTYdCAYs(3,3)      &
  & ,AYsCYsYsCAYs(3,3),AYsCYzTYzCAYs(3,3),AYsCAYdTYdCYs(3,3),AYsCAYsYsCYs(3,3)             &
  & ,AYsCAYzTYzCYs(3,3),AYtadjYeYeCAYt(3,3),AYtadjYzYzCAYt(3,3),AYtadjAYeYeCYt(3,3)        &
  & ,AYtadjAYzYzCYt(3,3),AYtCYtYtCAYt(3,3),AYtCAYtYtCYt(3,3),AYuadjYdYdadjAYu(3,3)         &
  & ,AYuadjYuYuadjAYu(3,3),AYuadjAYdYdadjYu(3,3),AYuadjAYuYuadjYu(3,3),AYzadjYeYeadjAYz(3,3)&
  & ,AYzadjYzYzadjAYz(3,3),AYzadjAYeYeadjYz(3,3),AYzadjAYzYzadjYz(3,3),AYzCYtYtadjAYz(3,3) &
  & ,AYzCAYtYtadjYz(3,3),adjYdmd2YdadjYd(3,3),adjYdmd2YdadjYu(3,3),adjYdmd2YsCYs(3,3)      &
  & ,adjYdmd2YzadjYz(3,3),adjYdYdmq2adjYd(3,3),adjYdYdmq2adjYu(3,3),adjYdYdadjYdmd2(3,3)   &
  & ,adjYdYdadjYdYd(3,3),adjYdYdadjYdYs(3,3),adjYdYdadjYdYz(3,3),adjYdYdadjYdAYd(3,3)      &
  & ,adjYdYdadjYdAYs(3,3),adjYdYdadjYdAYz(3,3),adjYdYdadjYumu2(3,3),adjYdYdadjYuYu(3,3)    &
  & ,adjYdYdadjYuAYu(3,3),adjYdYsCYsYd(3,3),adjYdYsCYsAYd(3,3),adjYdYzml2adjYz(3,3)        &
  & ,adjYdYzadjYzYd(3,3),adjYdYzadjYzAYd(3,3),adjYdAYdadjYdYd(3,3),adjYdAYdadjYdYs(3,3)    &
  & ,adjYdAYdadjYdYz(3,3),adjYdAYdadjYuYu(3,3),adjYdAYsCYsYd(3,3),adjYdAYzadjYzYd(3,3)     &
  & ,adjYeme2YeadjYe(3,3),adjYeme2YeadjYz(3,3),adjYeme2YeCYt(3,3),adjYeYeml2adjYe(3,3)     &
  & ,adjYeYeml2adjYz(3,3),adjYeYeml2CYt(3,3),adjYeYeadjYeme2(3,3),adjYeYeadjYeYe(3,3)      &
  & ,adjYeYeadjYeAYe(3,3),adjYeYeadjYzmd2(3,3),adjYeYeadjYzYd(3,3),adjYeYeadjYzYs(3,3)     &
  & ,adjYeYeadjYzYz(3,3),adjYeYeadjYzAYd(3,3),adjYeYeadjYzAYs(3,3),adjYeYeadjYzAYz(3,3)    &
  & ,adjYeYeCYtml2(3,3),adjYeYeCYtYt(3,3),adjYeYeCYtAYt(3,3),adjYeAYeadjYeYe(3,3)          &
  & ,adjYeAYeadjYzYd(3,3),adjYeAYeadjYzYs(3,3),adjYeAYeadjYzYz(3,3),adjYeAYeCYtYt(3,3)     &
  & ,adjYumu2YuadjYd(3,3),adjYumu2YuadjYu(3,3),adjYuYumq2adjYd(3,3),adjYuYumq2adjYu(3,3)   &
  & ,adjYuYuadjYdmd2(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdYs(3,3),adjYuYuadjYdYz(3,3)      &
  & ,adjYuYuadjYdAYd(3,3),adjYuYuadjYdAYs(3,3),adjYuYuadjYdAYz(3,3),adjYuYuadjYumu2(3,3)   &
  & ,adjYuYuadjYuYu(3,3),adjYuYuadjYuAYu(3,3),adjYuAYuadjYdYd(3,3),adjYuAYuadjYdYs(3,3)    &
  & ,adjYuAYuadjYdYz(3,3),adjYuAYuadjYuYu(3,3),adjYzmd2YdadjYd(3,3),adjYzmd2YsCYs(3,3)     &
  & ,adjYzmd2YzadjYe(3,3),adjYzmd2YzadjYz(3,3),adjYzmd2YzCYt(3,3),adjYzYdmq2adjYd(3,3)     &
  & ,adjYzYdadjYdYz(3,3),adjYzYdadjYdAYz(3,3),adjYzYsCYsYz(3,3),adjYzYsCYsAYz(3,3)         &
  & ,adjYzYzml2adjYe(3,3),adjYzYzml2adjYz(3,3),adjYzYzml2CYt(3,3),adjYzYzadjYeme2(3,3)     &
  & ,adjYzYzadjYeYe(3,3),adjYzYzadjYeAYe(3,3),adjYzYzadjYzmd2(3,3),adjYzYzadjYzYd(3,3)     &
  & ,adjYzYzadjYzYs(3,3),adjYzYzadjYzYz(3,3),adjYzYzadjYzAYd(3,3),adjYzYzadjYzAYs(3,3)     &
  & ,adjYzYzadjYzAYz(3,3),adjYzYzCYtml2(3,3),adjYzYzCYtYt(3,3),adjYzYzCYtAYt(3,3)          &
  & ,adjYzAYdadjYdYz(3,3),adjYzAYsCYsYz(3,3),adjYzAYzadjYeYe(3,3),adjYzAYzadjYzYd(3,3)     &
  & ,adjYzAYzadjYzYs(3,3),adjYzAYzadjYzYz(3,3),adjYzAYzCYtYt(3,3),CYdmq2TYdCYd(3,3)        &
  & ,CYdmq2TYdCYs(3,3),CYdmq2TYdCYz(3,3),CYdTYdmd2CYd(3,3),CYdTYdmd2CYs(3,3)               &
  & ,CYdTYdmd2CYz(3,3),CYdTYdCYdmq2(3,3),CYdTYdCYdTYd(3,3),CYdTYdCYdTAYd(3,3)              &
  & ,CYdTYdCYsmd2(3,3),CYdTYdCYsYd(3,3),CYdTYdCYsYs(3,3),CYdTYdCYsYz(3,3),CYdTYdCYsAYd(3,3)&
  & ,CYdTYdCYsAYs(3,3),CYdTYdCYsAYz(3,3),CYdTYdCYzml2(3,3),CYdTYdCYzYt(3,3),               &
  & CYdTYdCYzAYt(3,3),CYdTYuCYuTYd(3,3),CYdTYuCYuTAYd(3,3),CYdTAYdCYdTYd(3,3)              &
  & ,CYdTAYdCYsYd(3,3),CYdTAYdCYsYs(3,3),CYdTAYdCYsYz(3,3),CYdTAYdCYzYt(3,3)               &
  & ,CYdTAYuCYuTYd(3,3),CYeml2TYeCYe(3,3),CYeTYeme2CYe(3,3),CYeTYeCYeml2(3,3)              &
  & ,CYeTYeCYeYt(3,3),CYeTYeCYeAYt(3,3),CYeTAYeCYeYt(3,3),CYsmd2YdadjYd(3,3)               &
  & ,CYsmd2YsCYd(3,3),CYsmd2YsCYs(3,3),CYsmd2YsCYz(3,3),CYsmd2YzadjYz(3,3),CYsYdmq2adjYd(3,3)&
  & ,CYsYdadjYdmd2(3,3),CYsYdadjYdYs(3,3),CYsYdadjYdAYs(3,3),CYsYdadjAYdAYs(3,3)           &
  & ,CYsYsmd2CYd(3,3),CYsYsmd2CYs(3,3),CYsYsmd2CYz(3,3),CYsYsCYdmq2(3,3),CYsYsCYsmd2(3,3)  &
  & ,CYsYsCYsYd(3,3),CYsYsCYsYs(3,3),CYsYsCYsYz(3,3),CYsYsCYsAYd(3,3),CYsYsCYsAYs(3,3)     &
  & ,CYsYsCYsAYz(3,3),CYsYsCYzml2(3,3),CYsYsCYzYt(3,3),CYsYsCYzAYt(3,3),CYsYsCAYsAYs(3,3)  &
  & ,CYsYzml2adjYz(3,3),CYsYzadjYzmd2(3,3),CYsYzadjYzYs(3,3),CYsYzadjYzAYs(3,3)            &
  & ,CYsYzadjAYzAYs(3,3),CYsAYdadjYdYs(3,3),CYsAYdadjAYdYs(3,3),CYsAYsCYsYd(3,3)           &
  & ,CYsAYsCYsYs(3,3),CYsAYsCYsYz(3,3),CYsAYsCYzYt(3,3),CYsAYsCAYsYs(3,3),CYsAYzadjYzYs(3,3)&
  & ,CYsAYzadjAYzYs(3,3),CYtml2YtadjYe(3,3),CYtml2YtadjYz(3,3),CYtml2YtCYt(3,3)            &
  & ,CYtYtml2adjYe(3,3),CYtYtml2adjYz(3,3),CYtYtml2CYt(3,3),CYtYtadjYeme2(3,3)             &
  & ,CYtYtadjYeYe(3,3),CYtYtadjYeAYe(3,3),CYtYtadjYzmd2(3,3),CYtYtadjYzYd(3,3)             &
  & ,CYtYtadjYzYs(3,3),CYtYtadjYzYz(3,3),CYtYtadjYzAYd(3,3),CYtYtadjYzAYs(3,3)             &
  & ,CYtYtadjYzAYz(3,3),CYtYtadjAYeAYe(3,3),CYtYtadjAYzAYz(3,3),CYtYtCYtml2(3,3)           &
  & ,CYtYtCYtYt(3,3),CYtYtCYtAYt(3,3),CYtYtCAYtAYt(3,3),CYtAYtadjYeYe(3,3),CYtAYtadjYzYd(3,3)&
  & ,CYtAYtadjYzYs(3,3),CYtAYtadjYzYz(3,3),CYtAYtadjAYeYe(3,3),CYtAYtadjAYzYz(3,3)         &
  & ,CYtAYtCYtYt(3,3),CYtAYtCAYtYt(3,3),CYtTYeCYeYt(3,3),CYtTYeCYeAYt(3,3),CYtTYzCYzYt(3,3)&
  & ,CYtTYzCYzAYt(3,3),CYtTAYeCYeYt(3,3),CYtTAYzCYzYt(3,3),CYumq2TYuCYu(3,3)               &
  & ,CYuTYumu2CYu(3,3),CYuTYuCYumq2(3,3),CYzml2TYzCYd(3,3),CYzml2TYzCYs(3,3)               &
  & ,CYzml2TYzCYz(3,3),CYzYtCYtTYz(3,3),CYzYtCYtTAYz(3,3),CYzAYtCYtTYz(3,3),               &
  & CYzTYeCYeTYz(3,3),CYzTYeCYeTAYz(3,3),CYzTYzmd2CYd(3,3),CYzTYzmd2CYs(3,3)               &
  & ,CYzTYzmd2CYz(3,3),CYzTYzCYdmq2(3,3),CYzTYzCYsmd2(3,3),CYzTYzCYsYd(3,3),               &
  & CYzTYzCYsYs(3,3),CYzTYzCYsYz(3,3),CYzTYzCYsAYd(3,3),CYzTYzCYsAYs(3,3),CYzTYzCYsAYz(3,3)&
  & ,CYzTYzCYzml2(3,3),CYzTYzCYzYt(3,3),CYzTYzCYzAYt(3,3),CYzTYzCYzTYz(3,3),               &
  & CYzTYzCYzTAYz(3,3),CYzTAYeCYeTYz(3,3),CYzTAYzCYsYd(3,3),CYzTAYzCYsYs(3,3)              &
  & ,CYzTAYzCYsYz(3,3),CYzTAYzCYzYt(3,3),CYzTAYzCYzTYz(3,3),TYdCYdTYdCYd(3,3)              &
  & ,TYdCYdTAYdCAYd(3,3),TYdCYsYsCYd(3,3),TYdCYsAYsCAYd(3,3),TYdCYzTYzCYd(3,3)             &
  & ,TYdCYzTAYzCAYd(3,3),TYdCAYdTAYdCYd(3,3),TYdCAYsAYsCYd(3,3),TYdCAYzTAYzCYd(3,3)        &
  & ,TYeCYeTYeCYe(3,3),TYeCYeTAYeCAYe(3,3),TYeCAYeTAYeCYe(3,3),TYuCYuTYuCYu(3,3)           &
  & ,TYuCYuTAYuCAYu(3,3),TYuCAYuTAYuCYu(3,3),TYzCYdTYdCYz(3,3),TYzCYdTAYdCAYz(3,3)         &
  & ,TYzCYsYsCYz(3,3),TYzCYsAYsCAYz(3,3),TYzCYzTYzCYz(3,3),TYzCYzTAYzCAYz(3,3)
  Complex(dp) :: TYzCAYdTAYdCYz(3,3),TYzCAYsAYsCYz(3,3),TYzCAYzTAYzCYz(3,3),TAYdCYdTYdCAYd(3,3)        &
  & ,TAYdCYsYsCAYd(3,3),TAYdCYzTYzCAYd(3,3),TAYdCAYdTYdCYd(3,3),TAYdCAYsYsCYd(3,3)         &
  & ,TAYdCAYzTYzCYd(3,3),TAYeCYeTYeCAYe(3,3),TAYeCAYeTYeCYe(3,3),TAYuCYuTYuCAYu(3,3)       &
  & ,TAYuCAYuTYuCYu(3,3),TAYzCYdTYdCAYz(3,3),TAYzCYsYsCAYz(3,3),TAYzCYzTYzCAYz(3,3)        &
  & ,TAYzCAYdTYdCYz(3,3),TAYzCAYsYsCYz(3,3),TAYzCAYzTYzCYz(3,3),md2YdadjYdYdadjYd(3,3)     &
  & ,md2YdadjYdYsCYs(3,3),md2YdadjYuYuadjYd(3,3),md2YsCYdTYdCYs(3,3),md2YsCYsYdadjYd(3,3)  &
  & ,md2YsCYsYsCYs(3,3),md2YsCYsYzadjYz(3,3),md2YsCYzTYzCYs(3,3),md2YzadjYeYeadjYz(3,3)    &
  & ,md2YzadjYzYsCYs(3,3),md2YzadjYzYzadjYz(3,3),md2YzCYtYtadjYz(3,3),md2CYsYsCYsYs(3,3)   &
  & ,me2YeadjYeYeadjYe(3,3),me2YeadjYzYzadjYe(3,3),me2YeCYtYtadjYe(3,3),ml2YtadjYeYeCYt(3,3)&
  & ,ml2YtadjYzYzCYt(3,3),ml2YtCYtYtCYt(3,3),ml2adjYzYsCYsYz(3,3),ml2CYtYtCYtYt(3,3)       &
  & ,ml2TYeCYeTYeCYe(3,3),ml2TYzCYdTYdCYz(3,3),ml2TYzCYsYsCYz(3,3),ml2TYzCYzTYzCYz(3,3)    &
  & ,mq2adjYdYsCYsYd(3,3),mq2TYdCYdTYdCYd(3,3),mq2TYdCYsYsCYd(3,3),mq2TYdCYzTYzCYd(3,3)    &
  & ,mq2TYuCYuTYuCYu(3,3),mu2YuadjYdYdadjYu(3,3),mu2YuadjYuYuadjYu(3,3),Ydmq2adjYdYdadjYd(3,3)&
  & ,Ydmq2adjYdYsCYs(3,3),Ydmq2adjYdYzadjYz(3,3),Ydmq2adjYuYuadjYd(3,3),YdadjYdmd2YdadjYd(3,3)&
  & ,YdadjYdmd2YsCYs(3,3),YdadjYdmd2YzadjYz(3,3),YdadjYdYdmq2adjYd(3,3),YdadjYdYdadjYdmd2(3,3)&
  & ,YdadjYdYdadjYdYd(3,3),YdadjYdYdadjYdYs(3,3),YdadjYdYdadjYdYz(3,3),YdadjYdYdadjYdAYd(3,3)&
  & ,YdadjYdYdadjYdAYs(3,3),YdadjYdYdadjYdAYz(3,3),YdadjYdYsCYsYd(3,3),YdadjYdYsCYsAYd(3,3)&
  & ,YdadjYdYzml2adjYz(3,3),YdadjYdYzadjYzYd(3,3),YdadjYdYzadjYzAYd(3,3),YdadjYdAYdadjYdYd(3,3)&
  & ,YdadjYdAYdadjYdYs(3,3),YdadjYdAYdadjYdYz(3,3),YdadjYdAYsCYsYd(3,3),YdadjYdAYzadjYzYd(3,3)&
  & ,YdadjYumu2YuadjYd(3,3),YdadjYuYumq2adjYd(3,3),YdadjYuYuadjYdmd2(3,3),YdadjYuYuadjYdYd(3,3)&
  & ,YdadjYuYuadjYdYs(3,3),YdadjYuYuadjYdYz(3,3),YdadjYuYuadjYdAYd(3,3),YdadjYuYuadjYdAYs(3,3)&
  & ,YdadjYuYuadjYdAYz(3,3),YdadjYuYuadjYuYu(3,3),YdadjYuYuadjYuAYu(3,3),YdadjYuAYuadjYdYd(3,3)&
  & ,YdadjYuAYuadjYdYs(3,3),YdadjYuAYuadjYdYz(3,3),YdadjYuAYuadjYuYu(3,3),Yeml2adjYeYeadjYe(3,3)&
  & ,Yeml2adjYzYzadjYe(3,3),Yeml2CYtYtadjYe(3,3),YeadjYeme2YeadjYe(3,3),YeadjYeYeml2adjYe(3,3)&
  & ,YeadjYeYeadjYeme2(3,3),YeadjYeYeadjYeYe(3,3),YeadjYeYeadjYeAYe(3,3),YeadjYeAYeadjYeYe(3,3)&
  & ,YeadjYzmd2YzadjYe(3,3),YeadjYzYdadjYdYz(3,3),YeadjYzYdadjYdAYz(3,3),YeadjYzYsCYsYz(3,3)&
  & ,YeadjYzYsCYsAYz(3,3),YeadjYzYzml2adjYe(3,3),YeadjYzYzadjYeme2(3,3),YeadjYzYzadjYeYe(3,3)&
  & ,YeadjYzYzadjYeAYe(3,3),YeadjYzYzadjYzYz(3,3),YeadjYzYzadjYzAYz(3,3),YeadjYzAYdadjYdYz(3,3)&
  & ,YeadjYzAYsCYsYz(3,3),YeadjYzAYzadjYeYe(3,3),YeadjYzAYzadjYzYz(3,3),YeCYtml2YtadjYe(3,3)&
  & ,YeCYtYtml2adjYe(3,3),YeCYtYtadjYeme2(3,3),YeCYtYtadjYeYe(3,3),YeCYtYtadjYeAYe(3,3)    &
  & ,YeCYtYtCYtYt(3,3),YeCYtYtCYtAYt(3,3),YeCYtAYtadjYeYe(3,3),YeCYtAYtCYtYt(3,3)          &
  & ,YeCYtTYeCYeYt(3,3),YeCYtTYeCYeAYt(3,3),YeCYtTYzCYzYt(3,3),YeCYtTYzCYzAYt(3,3)         &
  & ,YeCYtTAYeCYeYt(3,3),YeCYtTAYzCYzYt(3,3),Ysmd2CYdTYdCYs(3,3),Ysmd2CYsYdadjYd(3,3)      &
  & ,Ysmd2CYsYsCYs(3,3),Ysmd2CYsYzadjYz(3,3),Ysmd2CYzTYzCYs(3,3),YsCYdmq2TYdCYs(3,3)       &
  & ,YsCYdTYdmd2CYs(3,3),YsCYdTYdCYdTYd(3,3),YsCYdTYdCYdTAYd(3,3),YsCYdTYdCYsmd2(3,3)      &
  & ,YsCYdTYdCYsYd(3,3),YsCYdTYdCYsYs(3,3),YsCYdTYdCYsYz(3,3),YsCYdTYdCYsAYd(3,3)          &
  & ,YsCYdTYdCYsAYs(3,3),YsCYdTYdCYsAYz(3,3),YsCYdTYuCYuTYd(3,3),YsCYdTYuCYuTAYd(3,3)      &
  & ,YsCYdTAYdCYdTYd(3,3),YsCYdTAYdCYsYd(3,3),YsCYdTAYdCYsYs(3,3),YsCYdTAYdCYsYz(3,3)      &
  & ,YsCYdTAYuCYuTYd(3,3),YsCYsmd2YdadjYd(3,3),YsCYsmd2YsCYs(3,3),YsCYsmd2YzadjYz(3,3)     &
  & ,YsCYsYdmq2adjYd(3,3),YsCYsYdadjYdmd2(3,3),YsCYsYdadjYdYs(3,3),YsCYsYdadjYdAYs(3,3)    &
  & ,YsCYsYsmd2CYs(3,3),YsCYsYsCYsmd2(3,3),YsCYsYsCYsYd(3,3),YsCYsYsCYsYs(3,3)             &
  & ,YsCYsYsCYsYz(3,3),YsCYsYsCYsAYd(3,3),YsCYsYsCYsAYs(3,3),YsCYsYsCYsAYz(3,3)            &
  & ,YsCYsYzml2adjYz(3,3),YsCYsYzadjYzmd2(3,3),YsCYsYzadjYzYs(3,3),YsCYsYzadjYzAYs(3,3)    &
  & ,YsCYsAYdadjYdYs(3,3),YsCYsAYsCYsYd(3,3),YsCYsAYsCYsYs(3,3),YsCYsAYsCYsYz(3,3)         &
  & ,YsCYsAYzadjYzYs(3,3),YsCYzml2TYzCYs(3,3),YsCYzYtCYtTYz(3,3),YsCYzYtCYtTAYz(3,3)       &
  & ,YsCYzAYtCYtTYz(3,3),YsCYzTYeCYeTYz(3,3),YsCYzTYeCYeTAYz(3,3),YsCYzTYzmd2CYs(3,3)      &
  & ,YsCYzTYzCYsmd2(3,3),YsCYzTYzCYsYd(3,3),YsCYzTYzCYsYs(3,3),YsCYzTYzCYsYz(3,3)          &
  & ,YsCYzTYzCYsAYd(3,3),YsCYzTYzCYsAYs(3,3),YsCYzTYzCYsAYz(3,3),YsCYzTYzCYzTYz(3,3)       &
  & ,YsCYzTYzCYzTAYz(3,3),YsCYzTAYeCYeTYz(3,3),YsCYzTAYzCYsYd(3,3),YsCYzTAYzCYsYs(3,3)     &
  & ,YsCYzTAYzCYsYz(3,3),YsCYzTAYzCYzTYz(3,3),Ytml2adjYeYeCYt(3,3),Ytml2adjYzYzCYt(3,3)    &
  & ,Ytml2CYtYtCYt(3,3),YtadjYeme2YeCYt(3,3),YtadjYeYeml2CYt(3,3),YtadjYeYeadjYeYe(3,3)    &
  & ,YtadjYeYeadjYeAYe(3,3),YtadjYeYeCYtml2(3,3),YtadjYeYeCYtYt(3,3),YtadjYeYeCYtAYt(3,3)  &
  & ,YtadjYeAYeadjYeYe(3,3),YtadjYeAYeCYtYt(3,3),YtadjYzmd2YzCYt(3,3),YtadjYzYdadjYdYz(3,3)&
  & ,YtadjYzYdadjYdAYz(3,3),YtadjYzYsCYsYz(3,3),YtadjYzYsCYsAYz(3,3),YtadjYzYzml2CYt(3,3)  &
  & ,YtadjYzYzadjYzYz(3,3),YtadjYzYzadjYzAYz(3,3),YtadjYzYzCYtml2(3,3),YtadjYzYzCYtYt(3,3) &
  & ,YtadjYzYzCYtAYt(3,3),YtadjYzAYdadjYdYz(3,3),YtadjYzAYsCYsYz(3,3),YtadjYzAYzadjYzYz(3,3)&
  & ,YtadjYzAYzCYtYt(3,3),YtCYtml2YtCYt(3,3),YtCYtYtml2CYt(3,3),YtCYtYtCYtml2(3,3)         &
  & ,YtCYtYtCYtYt(3,3),YtCYtYtCYtAYt(3,3),YtCYtAYtCYtYt(3,3),YtCYtTYeCYeYt(3,3)            &
  & ,YtCYtTYeCYeAYt(3,3),YtCYtTYzCYzYt(3,3),YtCYtTYzCYzAYt(3,3),YtCYtTAYeCYeYt(3,3)        &
  & ,YtCYtTAYzCYzYt(3,3),Yumq2adjYdYdadjYu(3,3),Yumq2adjYuYuadjYu(3,3),YuadjYdmd2YdadjYu(3,3)&
  & ,YuadjYdYdmq2adjYu(3,3),YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYdAYd(3,3),YuadjYdYdadjYumu2(3,3)&
  & ,YuadjYdYdadjYuYu(3,3),YuadjYdYdadjYuAYu(3,3),YuadjYdYsCYsYd(3,3),YuadjYdYsCYsAYd(3,3) &
  & ,YuadjYdYzadjYzYd(3,3),YuadjYdYzadjYzAYd(3,3),YuadjYdAYdadjYdYd(3,3),YuadjYdAYdadjYuYu(3,3)&
  & ,YuadjYdAYsCYsYd(3,3),YuadjYdAYzadjYzYd(3,3),YuadjYumu2YuadjYu(3,3),YuadjYuYumq2adjYu(3,3)&
  & ,YuadjYuYuadjYumu2(3,3),YuadjYuYuadjYuYu(3,3),YuadjYuYuadjYuAYu(3,3),YuadjYuAYuadjYuYu(3,3)&
  & ,Yzml2adjYeYeadjYz(3,3),Yzml2adjYzYdadjYd(3,3),Yzml2adjYzYsCYs(3,3),Yzml2adjYzYzadjYz(3,3)&
  & ,Yzml2CYtYtadjYz(3,3),YzadjYeme2YeadjYz(3,3),YzadjYeYeml2adjYz(3,3),YzadjYeYeadjYeYe(3,3)&
  & ,YzadjYeYeadjYeAYe(3,3),YzadjYeYeadjYzmd2(3,3),YzadjYeYeadjYzYd(3,3),YzadjYeYeadjYzYs(3,3)&
  & ,YzadjYeYeadjYzYz(3,3),YzadjYeYeadjYzAYd(3,3),YzadjYeYeadjYzAYs(3,3),YzadjYeYeadjYzAYz(3,3)&
  & ,YzadjYeAYeadjYeYe(3,3),YzadjYeAYeadjYzYd(3,3),YzadjYeAYeadjYzYs(3,3),YzadjYeAYeadjYzYz(3,3)&
  & ,YzadjYzmd2YdadjYd(3,3),YzadjYzmd2YsCYs(3,3),YzadjYzmd2YzadjYz(3,3),YzadjYzYdmq2adjYd(3,3)&
  & ,YzadjYzYdadjYdYz(3,3),YzadjYzYdadjYdAYz(3,3),YzadjYzYsCYsYz(3,3),YzadjYzYsCYsAYz(3,3) &
  & ,YzadjYzYzml2adjYz(3,3),YzadjYzYzadjYzmd2(3,3),YzadjYzYzadjYzYd(3,3),YzadjYzYzadjYzYs(3,3)&
  & ,YzadjYzYzadjYzYz(3,3),YzadjYzYzadjYzAYd(3,3),YzadjYzYzadjYzAYs(3,3),YzadjYzYzadjYzAYz(3,3)&
  & ,YzadjYzAYdadjYdYz(3,3),YzadjYzAYsCYsYz(3,3),YzadjYzAYzadjYzYd(3,3),YzadjYzAYzadjYzYs(3,3)&
  & ,YzadjYzAYzadjYzYz(3,3),YzCYtml2YtadjYz(3,3),YzCYtYtml2adjYz(3,3),YzCYtYtadjYzmd2(3,3) &
  & ,YzCYtYtadjYzYd(3,3),YzCYtYtadjYzYs(3,3),YzCYtYtadjYzYz(3,3),YzCYtYtadjYzAYd(3,3)      &
  & ,YzCYtYtadjYzAYs(3,3),YzCYtYtadjYzAYz(3,3),YzCYtYtCYtYt(3,3),YzCYtYtCYtAYt(3,3)        &
  & ,YzCYtAYtadjYzYd(3,3),YzCYtAYtadjYzYs(3,3),YzCYtAYtadjYzYz(3,3),YzCYtAYtCYtYt(3,3)     &
  & ,YzCYtTYeCYeYt(3,3),YzCYtTYeCYeAYt(3,3),YzCYtTYzCYzYt(3,3),YzCYtTYzCYzAYt(3,3)         &
  & ,YzCYtTAYeCYeYt(3,3),YzCYtTAYzCYzYt(3,3),AYdadjYdYdadjYdYd(3,3),AYdadjYdYdadjYdYs(3,3) &
  & ,AYdadjYdYdadjYdYz(3,3),AYdadjYdYsCYsYd(3,3),AYdadjYdYzadjYzYd(3,3),AYdadjYuYuadjYdYd(3,3)&
  & ,AYdadjYuYuadjYdYs(3,3),AYdadjYuYuadjYdYz(3,3),AYdadjYuYuadjYuYu(3,3),AYeadjYeYeadjYeYe(3,3)&
  & ,AYeadjYzYdadjYdYz(3,3),AYeadjYzYsCYsYz(3,3),AYeadjYzYzadjYeYe(3,3),AYeadjYzYzadjYzYz(3,3)&
  & ,AYeCYtYtadjYeYe(3,3),AYeCYtYtCYtYt(3,3),AYeCYtTYeCYeYt(3,3),AYeCYtTYzCYzYt(3,3)       &
  & ,AYsCYdTYdCYdTYd(3,3),AYsCYdTYdCYsYd(3,3),AYsCYdTYdCYsYs(3,3),AYsCYdTYdCYsYz(3,3)      &
  & ,AYsCYdTYuCYuTYd(3,3),AYsCYsYdadjYdYs(3,3),AYsCYsYsCYsYd(3,3),AYsCYsYsCYsYs(3,3)       &
  & ,AYsCYsYsCYsYz(3,3),AYsCYsYzadjYzYs(3,3),AYsCYzYtCYtTYz(3,3),AYsCYzTYeCYeTYz(3,3)      &
  & ,AYsCYzTYzCYsYd(3,3),AYsCYzTYzCYsYs(3,3),AYsCYzTYzCYsYz(3,3),AYsCYzTYzCYzTYz(3,3)      &
  & ,AYtadjYeYeadjYeYe(3,3),AYtadjYeYeCYtYt(3,3),AYtadjYzYdadjYdYz(3,3),AYtadjYzYsCYsYz(3,3)&
  & ,AYtadjYzYzadjYzYz(3,3),AYtadjYzYzCYtYt(3,3),AYtCYtYtCYtYt(3,3),AYtCYtTYeCYeYt(3,3)    &
  & ,AYtCYtTYzCYzYt(3,3),AYuadjYdYdadjYdYd(3,3),AYuadjYdYdadjYuYu(3,3),AYuadjYdYsCYsYd(3,3)&
  & ,AYuadjYdYzadjYzYd(3,3),AYuadjYuYuadjYuYu(3,3),AYzadjYeYeadjYeYe(3,3),AYzadjYeYeadjYzYd(3,3)&
  & ,AYzadjYeYeadjYzYs(3,3),AYzadjYeYeadjYzYz(3,3),AYzadjYzYdadjYdYz(3,3),AYzadjYzYsCYsYz(3,3)&
  & ,AYzadjYzYzadjYzYd(3,3),AYzadjYzYzadjYzYs(3,3),AYzadjYzYzadjYzYz(3,3),AYzCYtYtadjYzYd(3,3)&
  & ,AYzCYtYtadjYzYs(3,3),AYzCYtYtadjYzYz(3,3),AYzCYtYtCYtYt(3,3),AYzCYtTYeCYeYt(3,3)      &
  & ,AYzCYtTYzCYzYt(3,3),CYsmd2YsCYsYs(3,3),CYsYdmq2adjYdYs(3,3),CYsYdadjYdmd2Ys(3,3)      &
  & ,CYsYdadjYdYsmd2(3,3),CYsYsmd2CYsYs(3,3),CYsYsCYsmd2Ys(3,3),CYsYsCYsYsmd2(3,3)         &
  & ,CYsYzml2adjYzYs(3,3),CYsYzadjYzmd2Ys(3,3),CYsYzadjYzYsmd2(3,3),CYtml2YtCYtYt(3,3)     &
  & ,CYtYtml2adjYeYe(3,3),CYtYtml2adjYzYz(3,3),CYtYtml2CYtYt(3,3),CYtYtadjYeme2Ye(3,3)     &
  & ,CYtYtadjYeYeml2(3,3),CYtYtadjYzmd2Yz(3,3),CYtYtadjYzYzml2(3,3),CYtYtCYtml2Yt(3,3)     &
  & ,CYtYtCYtYtml2(3,3),TYdmd2CYdTYdCYd(3,3),TYdmd2CYsYsCYd(3,3),TYdmd2CYzTYzCYd(3,3)      &
  & ,TYdCYdmq2TYdCYd(3,3),TYdCYdTYdmd2CYd(3,3),TYdCYdTYdCYdmq2(3,3),TYdCYsmd2YsCYd(3,3)    &
  & ,TYdCYsYsmd2CYd(3,3),TYdCYsYsCYdmq2(3,3),TYdCYzml2TYzCYd(3,3),TYdCYzTYzmd2CYd(3,3)     &
  & ,TYdCYzTYzCYdmq2(3,3),TYeme2CYeTYeCYe(3,3),TYeCYeml2TYeCYe(3,3),TYeCYeTYeme2CYe(3,3)   &
  & ,TYeCYeTYeCYeml2(3,3),TYeCYeTYeCYeYt(3,3),TYeCYeTYeCYeAYt(3,3),TYeCYeTAYeCYeYt(3,3)    &
  & ,TYumu2CYuTYuCYu(3,3),TYuCYumq2TYuCYu(3,3),TYuCYuTYumu2CYu(3,3),TYuCYuTYuCYumq2(3,3)   &
  & ,TYzmd2CYdTYdCYz(3,3),TYzmd2CYsYsCYz(3,3),TYzmd2CYzTYzCYz(3,3),TYzCYdmq2TYdCYz(3,3)    &
  & ,TYzCYdTYdmd2CYz(3,3),TYzCYdTYdCYzml2(3,3),TYzCYdTYdCYzYt(3,3),TYzCYdTYdCYzAYt(3,3)    &
  & ,TYzCYdTAYdCYzYt(3,3),TYzCYsmd2YsCYz(3,3),TYzCYsYsmd2CYz(3,3),TYzCYsYsCYzml2(3,3)      &
  & ,TYzCYsYsCYzYt(3,3),TYzCYsYsCYzAYt(3,3),TYzCYsAYsCYzYt(3,3),TYzCYzml2TYzCYz(3,3)       &
  & ,TYzCYzTYzmd2CYz(3,3),TYzCYzTYzCYzml2(3,3),TYzCYzTYzCYzYt(3,3),TYzCYzTYzCYzAYt(3,3)    &
  & ,TYzCYzTAYzCYzYt(3,3),TAYeCYeTYeCYeYt(3,3),TAYzCYdTYdCYzYt(3,3),TAYzCYsYsCYzYt(3,3)    &
  & ,TAYzCYzTYzCYzYt(3,3)

  Complex(dp) :: Trmd2,Trme2,Trml2,Trmq2,Trmu2,TrYdadjYd,TrYeadjYe,TrYsCYs,             &
  & TrYtCYt,TrYuadjYu,TrYzadjYz,TrAYdadjYd,TrAYdadjAYd,TrAYeadjYe,TrAYeadjAYe,             &
  & TrAYsCAYs,TrAYtCAYt,TrAYuadjYu,TrAYuadjAYu,TrAYzadjYz,TrAYzadjAYz,TrCYsYs,             &
  & TrCYsAYs,TrCYtYt,TrCYtAYt,TrCAYsAYs,TrCAYtAYt,Trmd2YdadjYd,Trmd2YzadjYz,               &
  & Trme2YeadjYe,Trmu2YuadjYu,TrYdmq2adjYd,TrYeml2adjYe,TrYsCYsmd2,TrYtCYtml2,             &
  & TrYumq2adjYu,TrYzml2adjYz,TrCYsmd2Ys,TrCYsYsmd2,TrCYtml2Yt,TrCYtYtml2

  Complex(dp) :: TrYdadjAYd,TrYeadjAYe,TrYsCAYs,TrYtCAYt,TrYuadjAYu,TrYzadjAYz,         &
  & Trmd2YsCYs,Trml2YtCYt,TrYdadjYdYdadjYd,TrYdadjYdYsCYs,TrYdadjYdYzadjYz,TrYdadjYdAYdadjYd,&
  & TrYdadjYdAYdadjAYd,TrYdadjYdAYsCYs,TrYdadjYdAYsCAYs,TrYdadjYdAYzadjYz,TrYdadjYdAYzadjAYz,&
  & TrYdadjYuYuadjYd,TrYdadjYuAYuadjYd,TrYdadjYuAYuadjAYd,TrYdadjAYdAYdadjYd,              &
  & TrYdadjAYdAYsCYs,TrYdadjAYdAYzadjYz,TrYdadjAYuAYuadjYd,TrYeadjYeYeadjYe,               &
  & TrYeadjYeAYeadjYe,TrYeadjYeAYeadjAYe,TrYeadjYzYzadjYe,TrYeadjYzAYzadjYe,               &
  & TrYeadjYzAYzadjAYe,TrYeadjAYeAYeadjYe,TrYeadjAYzAYzadjYe,TrYeCYtYtadjYe,               &
  & TrYeCYtAYtadjYe,TrYeCYtAYtadjAYe,TrYeCAYtAYtadjYe,TrYsCYsYdadjYd,TrYsCYsYsCYs,         &
  & TrYsCYsYzadjYz,TrYsCYsAYdadjYd,TrYsCYsAYdadjAYd,TrYsCYsAYsCYs,TrYsCYsAYsCAYs,          &
  & TrYsCYsAYzadjYz,TrYsCYsAYzadjAYz,TrYsCAYsAYdadjYd,TrYsCAYsAYsCYs,TrYsCAYsAYzadjYz,     &
  & TrYtadjYeYeCYt,TrYtadjYeAYeCYt,TrYtadjYeAYeCAYt,TrYtadjYzYzCYt,TrYtadjYzAYzCYt,        &
  & TrYtadjYzAYzCAYt,TrYtadjAYeAYeCYt,TrYtadjAYzAYzCYt,TrYtCYtYtCYt,TrYtCYtAYtCYt,         &
  & TrYtCYtAYtCAYt,TrYtCAYtAYtCYt,TrYuadjYdYdadjYu,TrYuadjYdAYdadjYu,TrYuadjYdAYdadjAYu,   &
  & TrYuadjYuYuadjYu,TrYuadjYuAYuadjYu,TrYuadjYuAYuadjAYu,TrYuadjAYdAYdadjYu,              &
  & TrYuadjAYuAYuadjYu,TrYzadjYeYeadjYz,TrYzadjYeAYeadjYz,TrYzadjYeAYeadjAYz,              &
  & TrYzadjYzYdadjYd,TrYzadjYzYsCYs,TrYzadjYzYzadjYz,TrYzadjYzAYdadjYd,TrYzadjYzAYdadjAYd, &
  & TrYzadjYzAYsCYs,TrYzadjYzAYsCAYs,TrYzadjYzAYzadjYz,TrYzadjYzAYzadjAYz,TrYzadjAYeAYeadjYz,&
  & TrYzadjAYzAYdadjYd,TrYzadjAYzAYsCYs,TrYzadjAYzAYzadjYz,TrYzCYtYtadjYz,TrYzCYtAYtadjYz, &
  & TrYzCYtAYtadjAYz,TrYzCAYtAYtadjYz,TrCYsYdadjYdYs,TrCYsYdadjYdAYs,TrCYsYdadjAYdAYs,     &
  & TrCYsYsCYsYs,TrCYsYsCYsAYs,TrCYsYsCAYsAYs,TrCYsYzadjYzYs,TrCYsYzadjYzAYs,              &
  & TrCYsYzadjAYzAYs,TrCYsAYdadjYdYs,TrCYsAYdadjAYdYs,TrCYsAYsCYsYs,TrCYsAYsCAYsYs,        &
  & TrCYsAYzadjYzYs,TrCYsAYzadjAYzYs,TrCYtYtadjYeYe,TrCYtYtadjYeAYe,TrCYtYtadjYzYz,        &
  & TrCYtYtadjYzAYz,TrCYtYtadjAYeAYe,TrCYtYtadjAYzAYz,TrCYtYtCYtYt,TrCYtYtCYtAYt,          &
  & TrCYtYtCAYtAYt,TrCYtAYtadjYeYe,TrCYtAYtadjYzYz,TrCYtAYtadjAYeYe,TrCYtAYtadjAYzYz,      &
  & TrCYtAYtCYtYt,TrCYtAYtCAYtYt,Trmd2YdadjYdYsCYs,Trmd2YsCYsYdadjYd,Trmd2YsCYsYsCYs,      &
  & Trmd2YsCYsYzadjYz,Trmd2YzadjYzYsCYs,Trmd2YzCYtYtadjYz,Trmd2CYsYsCYsYs,Trme2YeCYtYtadjYe,&
  & Trml2YtCYtYtCYt,Trml2adjYzYsCYsYz,Trml2CYtYtCYtYt,Trmq2adjYdYsCYsYd,TrYdmq2adjYdYdadjYd,&
  & TrYdmq2adjYdYsCYs,TrYdmq2adjYdYzadjYz,TrYdmq2adjYuYuadjYd,TrYdadjYdmd2YdadjYd,         &
  & TrYdadjYdmd2YsCYs,TrYdadjYdmd2YzadjYz,TrYdadjYdYdmq2adjYd,TrYdadjYdYzml2adjYz,         &
  & TrYdadjYumu2YuadjYd,TrYdadjYuYumq2adjYd,TrYeml2adjYeYeadjYe,TrYeml2adjYzYzadjYe,       &
  & TrYeml2CYtYtadjYe,TrYeadjYeme2YeadjYe,TrYeadjYeYeml2adjYe,TrYeadjYzmd2YzadjYe,         &
  & TrYeadjYzYzml2adjYe,TrYeCYtml2YtadjYe,TrYeCYtYtml2adjYe,TrYsmd2CYsYdadjYd,             &
  & TrYsmd2CYsYsCYs,TrYsmd2CYsYzadjYz,TrYsCYsmd2YdadjYd,TrYsCYsmd2YsCYs,TrYsCYsmd2YzadjYz, &
  & TrYsCYsYdmq2adjYd,TrYsCYsYdadjYdmd2,TrYsCYsYsmd2CYs,TrYsCYsYsCYsmd2,TrYsCYsYzml2adjYz, &
  & TrYsCYsYzadjYzmd2,TrYtml2adjYeYeCYt,TrYtml2adjYzYzCYt,TrYtml2CYtYtCYt,TrYtadjYeme2YeCYt,&
  & TrYtadjYeYeml2CYt,TrYtadjYeYeCYtml2,TrYtadjYzmd2YzCYt,TrYtadjYzYzml2CYt,               &
  & TrYtadjYzYzCYtml2,TrYtCYtml2YtCYt,TrYtCYtYtml2CYt,TrYtCYtYtCYtml2,TrYumq2adjYdYdadjYu, &
  & TrYumq2adjYuYuadjYu,TrYuadjYdmd2YdadjYu,TrYuadjYdYdmq2adjYu,TrYuadjYumu2YuadjYu,       &
  & TrYuadjYuYumq2adjYu,TrYzml2adjYeYeadjYz,TrYzml2adjYzYdadjYd,TrYzml2adjYzYsCYs,         &
  & TrYzml2adjYzYzadjYz,TrYzml2CYtYtadjYz,TrYzadjYeme2YeadjYz,TrYzadjYeYeml2adjYz,         &
  & TrYzadjYzmd2YdadjYd,TrYzadjYzmd2YsCYs,TrYzadjYzmd2YzadjYz,TrYzadjYzYdmq2adjYd,         &
  & TrYzadjYzYzml2adjYz,TrYzCYtml2YtadjYz,TrYzCYtYtml2adjYz,TrCYsmd2YsCYsYs,               &
  & TrCYsYdmq2adjYdYs,TrCYsYdadjYdmd2Ys,TrCYsYdadjYdYsmd2,TrCYsYsmd2CYsYs,TrCYsYsCYsmd2Ys, &
  & TrCYsYsCYsYsmd2,TrCYsYzml2adjYzYs,TrCYsYzadjYzmd2Ys,TrCYsYzadjYzYsmd2,TrCYtml2YtCYtYt, &
  & TrCYtYtml2adjYeYe,TrCYtYtml2adjYzYz,TrCYtYtml2CYtYt,TrCYtYtadjYeme2Ye,TrCYtYtadjYeYeml2,&
  & TrCYtYtadjYzmd2Yz,TrCYtYtadjYzYzml2,TrCYtYtCYtml2Yt,TrCYtYtCYtYtml2

  Real(dp) :: g12, g13, g14, g22, g23, g24, g32, g33, g34, AbsL1sq, AbsL2sq

  Iname = Iname +1
  NameOfUnit(Iname) = 'rge353'

  OnlyDiagonal = .Not.GenerationMixing
  q = t

  Call GToParameters353(gy,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM,mue,MZM,MSM,      &
     & AYu,AYd,AYe,AYt,AYs,AYz,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2,ml2,mHd2,mHu2,md2,   &
     & mu2,me2,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG)

  g12 = g1**2
  g13 = g1*g12
  g14 = g12*g12
  g22 = g2**2
  g23 = g2*g22
  g24 = g22*g22
  g32 = g3**2
  g33 = g3*g32
  g34 = g32*g32
  AbsL1sq = Abs(L1)**2
  AbsL2sq = Abs(L2)**2

    Call Adjungate(Yu,adjYu)
    Call Adjungate(Yd,adjYd)
    Call Adjungate(Ye,adjYe)
    Call Adjungate(Yt,adjYt)
    Call Adjungate(Ys,adjYs)
    Call Adjungate(Yz,adjYz)
    Call Adjungate(AYu,adjAYu)
    Call Adjungate(AYd,adjAYd)
    Call Adjungate(AYe,adjAYe)
    Call Adjungate(AYt,adjAYt)
    Call Adjungate(AYs,adjAYs)
    Call Adjungate(AYz,adjAYz)
    Call Adjungate(mq2,adjmq2)
    Call Adjungate(ml2,adjml2)
    Call Adjungate(md2,adjmd2)
    Call Adjungate(mu2,adjmu2)
    Call Adjungate(me2,adjme2)
   md2Ys = Matmul(md2,Ys)
   md2CYd = Matmul(md2,Conjg(Yd))
   md2CYs = Matmul(md2,adjYs)
   md2CYz = Matmul(md2,Conjg(Yz))
   me2CYe = Matmul(me2,Conjg(Ye))
   ml2Yt = Matmul(ml2,Yt)
   ml2adjYe = Matmul(ml2,adjYe)
   ml2adjYz = Matmul(ml2,adjYz)
   ml2CYt = Matmul(ml2,adjYt)
   mq2adjYd = Matmul(mq2,adjYd)
   mq2adjYu = Matmul(mq2,adjYu)
   mu2CYu = Matmul(mu2,Conjg(Yu))
   YdadjYd = Matmul(Yd,adjYd)
   YeadjYe = Matmul(Ye,adjYe)
   Ysmd2 = Matmul(Ys,md2)
   YsCYs = Matmul(Ys,adjYs)
   Ytml2 = Matmul(Yt,ml2)
   YtCYt = Matmul(Yt,adjYt)
   YuadjYu = Matmul(Yu,adjYu)
   YzadjYz = Matmul(Yz,adjYz)
   AYdadjYd = Matmul(AYd,adjYd)
   AYdadjAYd = Matmul(AYd,adjAYd)
   AYeadjYe = Matmul(AYe,adjYe)
   AYeadjAYe = Matmul(AYe,adjAYe)
   AYsCAYs = Matmul(AYs,adjAYs)
   AYtCAYt = Matmul(AYt,adjAYt)
   AYuadjYu = Matmul(AYu,adjYu)
   AYuadjAYu = Matmul(AYu,adjAYu)
   AYzadjYz = Matmul(AYz,adjYz)
   AYzadjAYz = Matmul(AYz,adjAYz)
   adjYdmd2 = Matmul(adjYd,md2)
   adjYdYd = Matmul(adjYd,Yd)
   adjYdYs = Matmul(adjYd,Ys)
   adjYdYz = Matmul(adjYd,Yz)
   adjYdAYd = Matmul(adjYd,AYd)
   adjYdAYs = Matmul(adjYd,AYs)
   adjYdAYz = Matmul(adjYd,AYz)
   adjYeme2 = Matmul(adjYe,me2)
   adjYeYe = Matmul(adjYe,Ye)
   adjYeAYe = Matmul(adjYe,AYe)
   adjYumu2 = Matmul(adjYu,mu2)
   adjYuYu = Matmul(adjYu,Yu)
   adjYuAYu = Matmul(adjYu,AYu)
   adjYzmd2 = Matmul(adjYz,md2)
   adjYzYd = Matmul(adjYz,Yd)
   adjYzYs = Matmul(adjYz,Ys)
   adjYzYz = Matmul(adjYz,Yz)
   adjYzAYd = Matmul(adjYz,AYd)
   adjYzAYs = Matmul(adjYz,AYs)
   adjYzAYz = Matmul(adjYz,AYz)
   CYdmq2 = Matmul(Conjg(Yd),mq2)
   CYdTYd = Matmul(Conjg(Yd),Transpose(Yd))
   CYdTAYd = Matmul(Conjg(Yd),Transpose(AYd))
   CYeml2 = Matmul(Conjg(Ye),ml2)
   CYeYt = Matmul(Conjg(Ye),Yt)
   CYeAYt = Matmul(Conjg(Ye),AYt)
   CYsmd2 = Matmul(adjYs,md2)
   CYsYd = Matmul(adjYs,Yd)
   CYsYs = Matmul(adjYs,Ys)
   CYsYz = Matmul(adjYs,Yz)
   CYsAYd = Matmul(adjYs,AYd)
   CYsAYs = Matmul(adjYs,AYs)
   CYsAYz = Matmul(adjYs,AYz)
   CYtml2 = Matmul(adjYt,ml2)
   CYtYt = Matmul(adjYt,Yt)
   CYtAYt = Matmul(adjYt,AYt)
   CYumq2 = Matmul(Conjg(Yu),mq2)
   CYzml2 = Matmul(Conjg(Yz),ml2)
   CYzYt = Matmul(Conjg(Yz),Yt)
   CYzAYt = Matmul(Conjg(Yz),AYt)
   CYzTYz = Matmul(Conjg(Yz),Transpose(Yz))
   CYzTAYz = Matmul(Conjg(Yz),Transpose(AYz))
   CAYsAYs = Matmul(adjAYs,AYs)
   CAYtAYt = Matmul(adjAYt,AYt)
   TYdCYd = Matmul(Transpose(Yd),Conjg(Yd))
   TYeCYe = Matmul(Transpose(Ye),Conjg(Ye))
   TYuCYu = Matmul(Transpose(Yu),Conjg(Yu))
   TYzCYz = Matmul(Transpose(Yz),Conjg(Yz))
   TAYdCAYd = Matmul(Transpose(AYd),Conjg(AYd))
   TAYeCAYe = Matmul(Transpose(AYe),Conjg(AYe))
   TAYuCAYu = Matmul(Transpose(AYu),Conjg(AYu))
   TAYzCAYz = Matmul(Transpose(AYz),Conjg(AYz))
   md2YdadjYd = Matmul(md2,YdadjYd)
   md2YsCYs = Matmul(md2,YsCYs)
   md2YzadjYz = Matmul(md2,YzadjYz)
   me2YeadjYe = Matmul(me2,YeadjYe)
   ml2YtCYt = Matmul(ml2,YtCYt)
   ml2TYeCYe = Matmul(ml2,TYeCYe)
   ml2TYzCYz = Matmul(ml2,TYzCYz)
   mq2TYdCYd = Matmul(mq2,TYdCYd)
   mq2TYuCYu = Matmul(mq2,TYuCYu)
   mu2YuadjYu = Matmul(mu2,YuadjYu)
   Ydmq2adjYd = Matmul(Yd,mq2adjYd)
   YdadjYdmd2 = Matmul(Yd,adjYdmd2)
   YdadjYdYd = Matmul(Yd,adjYdYd)
   YdadjYdYs = Matmul(Yd,adjYdYs)
   YdadjYdYz = Matmul(Yd,adjYdYz)
   YdadjYdAYd = Matmul(Yd,adjYdAYd)
   YdadjYdAYs = Matmul(Yd,adjYdAYs)
   YdadjYdAYz = Matmul(Yd,adjYdAYz)
   YdadjYuYu = Matmul(Yd,adjYuYu)
   YdadjYuAYu = Matmul(Yd,adjYuAYu)
   Yeml2adjYe = Matmul(Ye,ml2adjYe)
   YeadjYeme2 = Matmul(Ye,adjYeme2)
   YeadjYeYe = Matmul(Ye,adjYeYe)
   YeadjYeAYe = Matmul(Ye,adjYeAYe)
   YeadjYzYz = Matmul(Ye,adjYzYz)
   YeadjYzAYz = Matmul(Ye,adjYzAYz)
   YeCYtYt = Matmul(Ye,CYtYt)
   YeCYtAYt = Matmul(Ye,CYtAYt)
   Ysmd2CYs = Matmul(Ys,md2CYs)
   YsCYdTYd = Matmul(Ys,CYdTYd)
   YsCYdTAYd = Matmul(Ys,CYdTAYd)
   YsCYsmd2 = Matmul(Ys,CYsmd2)
   YsCYsYd = Matmul(Ys,CYsYd)
   YsCYsYs = Matmul(Ys,CYsYs)
   YsCYsYz = Matmul(Ys,CYsYz)
   YsCYsAYd = Matmul(Ys,CYsAYd)
   YsCYsAYs = Matmul(Ys,CYsAYs)
   YsCYsAYz = Matmul(Ys,CYsAYz)
   YsCYzTYz = Matmul(Ys,CYzTYz)
   YsCYzTAYz = Matmul(Ys,CYzTAYz)
   Ytml2CYt = Matmul(Yt,ml2CYt)
   YtadjYeYe = Matmul(Yt,adjYeYe)
   YtadjYeAYe = Matmul(Yt,adjYeAYe)
   YtadjYzYz = Matmul(Yt,adjYzYz)
   YtadjYzAYz = Matmul(Yt,adjYzAYz)
   YtCYtml2 = Matmul(Yt,CYtml2)
   YtCYtYt = Matmul(Yt,CYtYt)
   YtCYtAYt = Matmul(Yt,CYtAYt)
   Yumq2adjYu = Matmul(Yu,mq2adjYu)
   YuadjYdYd = Matmul(Yu,adjYdYd)
   YuadjYdAYd = Matmul(Yu,adjYdAYd)
   YuadjYumu2 = Matmul(Yu,adjYumu2)
   YuadjYuYu = Matmul(Yu,adjYuYu)
   YuadjYuAYu = Matmul(Yu,adjYuAYu)
   Yzml2adjYz = Matmul(Yz,ml2adjYz)
   YzadjYeYe = Matmul(Yz,adjYeYe)
   YzadjYeAYe = Matmul(Yz,adjYeAYe)
   YzadjYzmd2 = Matmul(Yz,adjYzmd2)
   YzadjYzYd = Matmul(Yz,adjYzYd)
   YzadjYzYs = Matmul(Yz,adjYzYs)
   YzadjYzYz = Matmul(Yz,adjYzYz)
   YzadjYzAYd = Matmul(Yz,adjYzAYd)
   YzadjYzAYs = Matmul(Yz,adjYzAYs)
   YzadjYzAYz = Matmul(Yz,adjYzAYz)
   YzCYtYt = Matmul(Yz,CYtYt)
   YzCYtAYt = Matmul(Yz,CYtAYt)
   AYdadjYdYd = Matmul(AYd,adjYdYd)
   AYdadjYdYs = Matmul(AYd,adjYdYs)
   AYdadjYdYz = Matmul(AYd,adjYdYz)
   AYdadjYuYu = Matmul(AYd,adjYuYu)
   AYeadjYeYe = Matmul(AYe,adjYeYe)
   AYeadjYzYz = Matmul(AYe,adjYzYz)
   AYeCYtYt = Matmul(AYe,CYtYt)
   AYsCYdTYd = Matmul(AYs,CYdTYd)
   AYsCYsYd = Matmul(AYs,CYsYd)
   AYsCYsYs = Matmul(AYs,CYsYs)
   AYsCYsYz = Matmul(AYs,CYsYz)
   AYsCYzTYz = Matmul(AYs,CYzTYz)
   AYtadjYeYe = Matmul(AYt,adjYeYe)
   AYtadjYzYz = Matmul(AYt,adjYzYz)
   AYtCYtYt = Matmul(AYt,CYtYt)
   AYuadjYdYd = Matmul(AYu,adjYdYd)
   AYuadjYuYu = Matmul(AYu,adjYuYu)
   AYzadjYeYe = Matmul(AYz,adjYeYe)
   AYzadjYzYd = Matmul(AYz,adjYzYd)
   AYzadjYzYs = Matmul(AYz,adjYzYs)
   AYzadjYzYz = Matmul(AYz,adjYzYz)
   AYzCYtYt = Matmul(AYz,CYtYt)
   CYsmd2Ys = Matmul(adjYs,md2Ys)
   CYsYsmd2 = Matmul(adjYs,Ysmd2)
   CYtml2Yt = Matmul(adjYt,ml2Yt)
   CYtYtml2 = Matmul(adjYt,Ytml2)
   TYdmd2CYd = Matmul(Transpose(Yd),md2CYd)
   TYdCYdmq2 = Matmul(Transpose(Yd),CYdmq2)
   TYeme2CYe = Matmul(Transpose(Ye),me2CYe)
   TYeCYeml2 = Matmul(Transpose(Ye),CYeml2)
   TYeCYeYt = Matmul(Transpose(Ye),CYeYt)
   TYeCYeAYt = Matmul(Transpose(Ye),CYeAYt)
   TYumu2CYu = Matmul(Transpose(Yu),mu2CYu)
   TYuCYumq2 = Matmul(Transpose(Yu),CYumq2)
   TYzmd2CYz = Matmul(Transpose(Yz),md2CYz)
   TYzCYzml2 = Matmul(Transpose(Yz),CYzml2)
   TYzCYzYt = Matmul(Transpose(Yz),CYzYt)
   TYzCYzAYt = Matmul(Transpose(Yz),CYzAYt)
   TAYeCYeYt = Matmul(Transpose(AYe),CYeYt)
   TAYzCYzYt = Matmul(Transpose(AYz),CYzYt)
   Trmd2 = Real(cTrace(md2),dp)
   Trme2 = Real(cTrace(me2),dp)
   Trml2 = Real(cTrace(ml2),dp)
   Trmq2 = Real(cTrace(mq2),dp)
   Trmu2 = Real(cTrace(mu2),dp)
   TrYdadjYd = Real(cTrace(YdadjYd),dp)
   TrYeadjYe = Real(cTrace(YeadjYe),dp)
   TrYsCYs = Real(cTrace(YsCYs),dp)
   TrYtCYt = Real(cTrace(YtCYt),dp)
   TrYuadjYu = Real(cTrace(YuadjYu),dp)
   TrYzadjYz = Real(cTrace(YzadjYz),dp)
   TrAYdadjYd = Real(cTrace(AYdadjYd),dp)
   TrAYdadjAYd = Real(cTrace(AYdadjAYd),dp)
   TrAYeadjYe = Real(cTrace(AYeadjYe),dp)
   TrAYeadjAYe = Real(cTrace(AYeadjAYe),dp)
   TrAYsCAYs = Real(cTrace(AYsCAYs),dp)
   TrAYtCAYt = Real(cTrace(AYtCAYt),dp)
   TrAYuadjYu = Real(cTrace(AYuadjYu),dp)
   TrAYuadjAYu = Real(cTrace(AYuadjAYu),dp)
   TrAYzadjYz = Real(cTrace(AYzadjYz),dp)
   TrAYzadjAYz = Real(cTrace(AYzadjAYz),dp)
   TrCYsYs = Real(cTrace(CYsYs),dp)
   TrCYsAYs = Real(cTrace(CYsAYs),dp)
   TrCYtYt = Real(cTrace(CYtYt),dp)
   TrCYtAYt = Real(cTrace(CYtAYt),dp)
   TrCAYsAYs = Real(cTrace(CAYsAYs),dp)
   TrCAYtAYt = Real(cTrace(CAYtAYt),dp)
   Trmd2YdadjYd = Real(cTrace(md2YdadjYd),dp)
   Trmd2YzadjYz = Real(cTrace(md2YzadjYz),dp)
   Trme2YeadjYe = Real(cTrace(me2YeadjYe),dp)
   Trmu2YuadjYu = Real(cTrace(mu2YuadjYu),dp)
   TrYdmq2adjYd = Real(cTrace(Ydmq2adjYd),dp)
   TrYeml2adjYe = Real(cTrace(Yeml2adjYe),dp)
   TrYsCYsmd2 = Real(cTrace(YsCYsmd2),dp)
   TrYtCYtml2 = Real(cTrace(YtCYtml2),dp)
   TrYumq2adjYu = Real(cTrace(Yumq2adjYu),dp)
   TrYzml2adjYz = Real(cTrace(Yzml2adjYz),dp)
   TrCYsmd2Ys = Real(cTrace(CYsmd2Ys),dp)
   TrCYsYsmd2 = Real(cTrace(CYsYsmd2),dp)
   TrCYtml2Yt = Real(cTrace(CYtml2Yt),dp)
   TrCYtYtml2 = Real(cTrace(CYtYtml2),dp)


  If (TwoLoopRGE) Then
   md2Yz = Matmul(md2,Yz)
   me2Ye = Matmul(me2,Ye)
   YdadjYu = Matmul(Yd,adjYu)
   YdadjAYd = Matmul(Yd,adjAYd)
   YdadjAYu = Matmul(Yd,adjAYu)
   Yeml2 = Matmul(Ye,ml2)
   YeadjYz = Matmul(Ye,adjYz)
   YeadjAYe = Matmul(Ye,adjAYe)
   YeadjAYz = Matmul(Ye,adjAYz)
   YeCYt = Matmul(Ye,adjYt)
   YeCAYt = Matmul(Ye,adjAYt)
   YsCYd = Matmul(Ys,Conjg(Yd))
   YsCYz = Matmul(Ys,Conjg(Yz))
   YsCAYd = Matmul(Ys,Conjg(AYd))
   YsCAYs = Matmul(Ys,adjAYs)
   YsCAYz = Matmul(Ys,Conjg(AYz))
   YtadjYe = Matmul(Yt,adjYe)
   YtadjYz = Matmul(Yt,adjYz)
   YtadjAYe = Matmul(Yt,adjAYe)
   YtadjAYz = Matmul(Yt,adjAYz)
   YtCAYt = Matmul(Yt,adjAYt)
   YuadjYd = Matmul(Yu,adjYd)
   YuadjAYd = Matmul(Yu,adjAYd)
   YuadjAYu = Matmul(Yu,adjAYu)
   Yzml2 = Matmul(Yz,ml2)
   YzadjYe = Matmul(Yz,adjYe)
   YzadjAYe = Matmul(Yz,adjAYe)
   YzadjAYz = Matmul(Yz,adjAYz)
   YzCYt = Matmul(Yz,adjYt)
   YzCAYt = Matmul(Yz,adjAYt)
   AYdadjYu = Matmul(AYd,adjYu)
   AYdadjAYu = Matmul(AYd,adjAYu)
   AYeadjYz = Matmul(AYe,adjYz)
   AYeadjAYz = Matmul(AYe,adjAYz)
   AYeCYt = Matmul(AYe,adjYt)
   AYeCAYt = Matmul(AYe,adjAYt)
   AYsCYd = Matmul(AYs,Conjg(Yd))
   AYsCYs = Matmul(AYs,adjYs)
   AYsCYz = Matmul(AYs,Conjg(Yz))
   AYsCAYd = Matmul(AYs,Conjg(AYd))
   AYsCAYz = Matmul(AYs,Conjg(AYz))
   AYtadjYe = Matmul(AYt,adjYe)
   AYtadjYz = Matmul(AYt,adjYz)
   AYtadjAYe = Matmul(AYt,adjAYe)
   AYtadjAYz = Matmul(AYt,adjAYz)
   AYtCYt = Matmul(AYt,adjYt)
   AYuadjYd = Matmul(AYu,adjYd)
   AYuadjAYd = Matmul(AYu,adjAYd)
   AYzadjYe = Matmul(AYz,adjYe)
   AYzadjAYe = Matmul(AYz,adjAYe)
   AYzCYt = Matmul(AYz,adjYt)
   AYzCAYt = Matmul(AYz,adjAYt)
   adjAYdYs = Matmul(adjAYd,Ys)
   adjAYdAYs = Matmul(adjAYd,AYs)
   adjAYeYe = Matmul(adjAYe,Ye)
   adjAYeAYe = Matmul(adjAYe,AYe)
   adjAYzYs = Matmul(adjAYz,Ys)
   adjAYzYz = Matmul(adjAYz,Yz)
   adjAYzAYs = Matmul(adjAYz,AYs)
   adjAYzAYz = Matmul(adjAYz,AYz)
   CYeTYz = Matmul(Conjg(Ye),Transpose(Yz))
   CYeTAYz = Matmul(Conjg(Ye),Transpose(AYz))
   CYtTYz = Matmul(adjYt,Transpose(Yz))
   CYtTAYz = Matmul(adjYt,Transpose(AYz))
   CYuTYd = Matmul(Conjg(Yu),Transpose(Yd))
   CYuTAYd = Matmul(Conjg(Yu),Transpose(AYd))
   CAYsYs = Matmul(adjAYs,Ys)
   CAYtYt = Matmul(adjAYt,Yt)
   TYdCYs = Matmul(Transpose(Yd),adjYs)
   TYdCYz = Matmul(Transpose(Yd),Conjg(Yz))
   TYdCAYd = Matmul(Transpose(Yd),Conjg(AYd))
   TYdCAYs = Matmul(Transpose(Yd),adjAYs)
   TYdCAYz = Matmul(Transpose(Yd),Conjg(AYz))
   TYeCAYe = Matmul(Transpose(Ye),Conjg(AYe))
   TYuCAYu = Matmul(Transpose(Yu),Conjg(AYu))
   TYzCYd = Matmul(Transpose(Yz),Conjg(Yd))
   TYzCYs = Matmul(Transpose(Yz),adjYs)
   TYzCAYd = Matmul(Transpose(Yz),Conjg(AYd))
   TYzCAYs = Matmul(Transpose(Yz),adjAYs)
   TYzCAYz = Matmul(Transpose(Yz),Conjg(AYz))
   TAYdCYd = Matmul(Transpose(AYd),Conjg(Yd))
   TAYdCYs = Matmul(Transpose(AYd),adjYs)
   TAYdCYz = Matmul(Transpose(AYd),Conjg(Yz))
   TAYdCAYs = Matmul(Transpose(AYd),adjAYs)
   TAYdCAYz = Matmul(Transpose(AYd),Conjg(AYz))
   TAYeCYe = Matmul(Transpose(AYe),Conjg(Ye))
   TAYuCYu = Matmul(Transpose(AYu),Conjg(Yu))
   TAYzCYd = Matmul(Transpose(AYz),Conjg(Yd))
   TAYzCYs = Matmul(Transpose(AYz),adjYs)
   TAYzCYz = Matmul(Transpose(AYz),Conjg(Yz))
   TAYzCAYd = Matmul(Transpose(AYz),Conjg(AYd))
   TAYzCAYs = Matmul(Transpose(AYz),adjAYs)
   md2YdadjYu = Matmul(md2,YdadjYu)
   md2YsCYd = Matmul(md2,YsCYd)
   md2YsCYz = Matmul(md2,YsCYz)
   md2YzadjYe = Matmul(md2,YzadjYe)
   md2YzCYt = Matmul(md2,YzCYt)
   md2CYsYs = Matmul(md2,CYsYs)
   me2YeadjYz = Matmul(me2,YeadjYz)
   me2YeCYt = Matmul(me2,YeCYt)
   ml2YtadjYe = Matmul(ml2,YtadjYe)
   ml2YtadjYz = Matmul(ml2,YtadjYz)
   ml2adjYeYe = Matmul(ml2,adjYeYe)
   ml2adjYzYs = Matmul(ml2,adjYzYs)
   ml2adjYzYz = Matmul(ml2,adjYzYz)
   ml2CYtYt = Matmul(ml2,CYtYt)
   ml2TYzCYd = Matmul(ml2,TYzCYd)
   ml2TYzCYs = Matmul(ml2,TYzCYs)
   mq2adjYdYs = Matmul(mq2,adjYdYs)
   mq2TYdCYs = Matmul(mq2,TYdCYs)
   mq2TYdCYz = Matmul(mq2,TYdCYz)
   mu2YuadjYd = Matmul(mu2,YuadjYd)
   Ydmq2adjYu = Matmul(Yd,mq2adjYu)
   YdadjYumu2 = Matmul(Yd,adjYumu2)
   YdadjAYdAYs = Matmul(Yd,adjAYdAYs)
   Yeml2adjYz = Matmul(Ye,ml2adjYz)
   Yeml2CYt = Matmul(Ye,ml2CYt)
   YeadjYzmd2 = Matmul(Ye,adjYzmd2)
   YeadjYzYd = Matmul(Ye,adjYzYd)
   YeadjYzYs = Matmul(Ye,adjYzYs)
   YeadjYzAYd = Matmul(Ye,adjYzAYd)
   YeadjYzAYs = Matmul(Ye,adjYzAYs)
   YeCYtml2 = Matmul(Ye,CYtml2)
   Ysmd2CYd = Matmul(Ys,md2CYd)
   Ysmd2CYz = Matmul(Ys,md2CYz)
   YsCYdmq2 = Matmul(Ys,CYdmq2)
   YsCYzml2 = Matmul(Ys,CYzml2)
   YsCYzYt = Matmul(Ys,CYzYt)
   YsCYzAYt = Matmul(Ys,CYzAYt)
   YsCAYsAYs = Matmul(Ys,CAYsAYs)
   Ytml2adjYe = Matmul(Yt,ml2adjYe)
   Ytml2adjYz = Matmul(Yt,ml2adjYz)
   YtadjYeme2 = Matmul(Yt,adjYeme2)
   YtadjYzmd2 = Matmul(Yt,adjYzmd2)
   YtadjYzYd = Matmul(Yt,adjYzYd)
   YtadjYzYs = Matmul(Yt,adjYzYs)
   YtadjYzAYd = Matmul(Yt,adjYzAYd)
   YtadjYzAYs = Matmul(Yt,adjYzAYs)
   YtadjAYeAYe = Matmul(Yt,adjAYeAYe)
   YtadjAYzAYz = Matmul(Yt,adjAYzAYz)
   YtCYtTYz = Matmul(Yt,CYtTYz)
   YtCYtTAYz = Matmul(Yt,CYtTAYz)
   YtCAYtAYt = Matmul(Yt,CAYtAYt)
   Yumq2adjYd = Matmul(Yu,mq2adjYd)
   YuadjYdmd2 = Matmul(Yu,adjYdmd2)
   YuadjYdYs = Matmul(Yu,adjYdYs)
   YuadjYdYz = Matmul(Yu,adjYdYz)
   YuadjYdAYs = Matmul(Yu,adjYdAYs)
   YuadjYdAYz = Matmul(Yu,adjYdAYz)
   Yzml2adjYe = Matmul(Yz,ml2adjYe)
   Yzml2CYt = Matmul(Yz,ml2CYt)
   YzadjYeme2 = Matmul(Yz,adjYeme2)
   YzadjAYzAYs = Matmul(Yz,adjAYzAYs)
   YzCYtml2 = Matmul(Yz,CYtml2)
   AYdadjAYdYs = Matmul(AYd,adjAYdYs)
   AYeadjYzYd = Matmul(AYe,adjYzYd)
   AYeadjYzYs = Matmul(AYe,adjYzYs)
   AYsCYzYt = Matmul(AYs,CYzYt)
   AYsCAYsYs = Matmul(AYs,CAYsYs)
   AYtadjYzYd = Matmul(AYt,adjYzYd)
   AYtadjYzYs = Matmul(AYt,adjYzYs)
   AYtadjAYeYe = Matmul(AYt,adjAYeYe)
   AYtadjAYzYz = Matmul(AYt,adjAYzYz)
   AYtCYtTYz = Matmul(AYt,CYtTYz)
   AYtCAYtYt = Matmul(AYt,CAYtYt)
   AYuadjYdYs = Matmul(AYu,adjYdYs)
   AYuadjYdYz = Matmul(AYu,adjYdYz)
   AYzadjAYzYs = Matmul(AYz,adjAYzYs)
   adjYdmd2Ys = Matmul(adjYd,md2Ys)
   adjYdYdadjYd = Matmul(adjYd,YdadjYd)
   adjYdYdadjYu = Matmul(adjYd,YdadjYu)
   adjYdYdadjAYd = Matmul(adjYd,YdadjAYd)
   adjYdYdadjAYu = Matmul(adjYd,YdadjAYu)
   adjYdYsmd2 = Matmul(adjYd,Ysmd2)
   adjYdYsCYs = Matmul(adjYd,YsCYs)
   adjYdYzadjYz = Matmul(adjYd,YzadjYz)
   adjYdAYdadjYd = Matmul(adjYd,AYdadjYd)
   adjYdAYdadjYu = Matmul(adjYd,AYdadjYu)
   adjYdAYdadjAYd = Matmul(adjYd,AYdadjAYd)
   adjYdAYdadjAYu = Matmul(adjYd,AYdadjAYu)
   adjYdAYsCYs = Matmul(adjYd,AYsCYs)
   adjYdAYsCAYs = Matmul(adjYd,AYsCAYs)
   adjYdAYzadjYz = Matmul(adjYd,AYzadjYz)
   adjYdAYzadjAYz = Matmul(adjYd,AYzadjAYz)
   adjYeme2Ye = Matmul(adjYe,me2Ye)
   adjYeYeml2 = Matmul(adjYe,Yeml2)
   adjYeYeadjYe = Matmul(adjYe,YeadjYe)
   adjYeYeadjYz = Matmul(adjYe,YeadjYz)
   adjYeYeadjAYe = Matmul(adjYe,YeadjAYe)
   adjYeYeadjAYz = Matmul(adjYe,YeadjAYz)
   adjYeYeCYt = Matmul(adjYe,YeCYt)
   adjYeYeCAYt = Matmul(adjYe,YeCAYt)
   adjYeAYeadjYe = Matmul(adjYe,AYeadjYe)
   adjYeAYeadjYz = Matmul(adjYe,AYeadjYz)
   adjYeAYeadjAYe = Matmul(adjYe,AYeadjAYe)
   adjYeAYeadjAYz = Matmul(adjYe,AYeadjAYz)
   adjYeAYeCYt = Matmul(adjYe,AYeCYt)
   adjYeAYeCAYt = Matmul(adjYe,AYeCAYt)
   adjYuYuadjYd = Matmul(adjYu,YuadjYd)
   adjYuYuadjYu = Matmul(adjYu,YuadjYu)
   adjYuYuadjAYd = Matmul(adjYu,YuadjAYd)
   adjYuYuadjAYu = Matmul(adjYu,YuadjAYu)
   adjYuAYuadjYd = Matmul(adjYu,AYuadjYd)
   adjYuAYuadjYu = Matmul(adjYu,AYuadjYu)
   adjYuAYuadjAYd = Matmul(adjYu,AYuadjAYd)
   adjYuAYuadjAYu = Matmul(adjYu,AYuadjAYu)
   adjYzmd2Ys = Matmul(adjYz,md2Ys)
   adjYzmd2Yz = Matmul(adjYz,md2Yz)
   adjYzYdadjYd = Matmul(adjYz,YdadjYd)
   adjYzYsmd2 = Matmul(adjYz,Ysmd2)
   adjYzYsCYs = Matmul(adjYz,YsCYs)
   adjYzYzml2 = Matmul(adjYz,Yzml2)
   adjYzYzadjYe = Matmul(adjYz,YzadjYe)
   adjYzYzadjYz = Matmul(adjYz,YzadjYz)
   adjYzYzadjAYe = Matmul(adjYz,YzadjAYe)
   adjYzYzadjAYz = Matmul(adjYz,YzadjAYz)
   adjYzYzCYt = Matmul(adjYz,YzCYt)
   adjYzYzCAYt = Matmul(adjYz,YzCAYt)
   adjYzAYdadjYd = Matmul(adjYz,AYdadjYd)
   adjYzAYdadjAYd = Matmul(adjYz,AYdadjAYd)
   adjYzAYsCYs = Matmul(adjYz,AYsCYs)
   adjYzAYsCAYs = Matmul(adjYz,AYsCAYs)
   adjYzAYzadjYe = Matmul(adjYz,AYzadjYe)
   adjYzAYzadjYz = Matmul(adjYz,AYzadjYz)
   adjYzAYzadjAYe = Matmul(adjYz,AYzadjAYe)
   adjYzAYzadjAYz = Matmul(adjYz,AYzadjAYz)
   adjYzAYzCYt = Matmul(adjYz,AYzCYt)
   adjYzAYzCAYt = Matmul(adjYz,AYzCAYt)
   adjAYdYdadjYd = Matmul(adjAYd,YdadjYd)
   adjAYdYdadjYu = Matmul(adjAYd,YdadjYu)
   adjAYdAYdadjYd = Matmul(adjAYd,AYdadjYd)
   adjAYdAYdadjYu = Matmul(adjAYd,AYdadjYu)
   adjAYdAYsCYs = Matmul(adjAYd,AYsCYs)
   adjAYdAYzadjYz = Matmul(adjAYd,AYzadjYz)
   adjAYeYeadjYe = Matmul(adjAYe,YeadjYe)
   adjAYeYeadjYz = Matmul(adjAYe,YeadjYz)
   adjAYeYeCYt = Matmul(adjAYe,YeCYt)
   adjAYeAYeadjYe = Matmul(adjAYe,AYeadjYe)
   adjAYeAYeadjYz = Matmul(adjAYe,AYeadjYz)
   adjAYeAYeCYt = Matmul(adjAYe,AYeCYt)
   adjAYuYuadjYd = Matmul(adjAYu,YuadjYd)
   adjAYuYuadjYu = Matmul(adjAYu,YuadjYu)
   adjAYuAYuadjYd = Matmul(adjAYu,AYuadjYd)
   adjAYuAYuadjYu = Matmul(adjAYu,AYuadjYu)
   adjAYzYzadjYe = Matmul(adjAYz,YzadjYe)
   adjAYzYzadjYz = Matmul(adjAYz,YzadjYz)
   adjAYzYzCYt = Matmul(adjAYz,YzCYt)
   adjAYzAYdadjYd = Matmul(adjAYz,AYdadjYd)
   adjAYzAYsCYs = Matmul(adjAYz,AYsCYs)
   adjAYzAYzadjYe = Matmul(adjAYz,AYzadjYe)
   adjAYzAYzadjYz = Matmul(adjAYz,AYzadjYz)
   adjAYzAYzCYt = Matmul(adjAYz,AYzCYt)
   CYdTYdCYd = Matmul(Conjg(Yd),TYdCYd)
   CYdTYdCYs = Matmul(Conjg(Yd),TYdCYs)
   CYdTYdCYz = Matmul(Conjg(Yd),TYdCYz)
   CYdTYdCAYd = Matmul(Conjg(Yd),TYdCAYd)
   CYdTYdCAYs = Matmul(Conjg(Yd),TYdCAYs)
   CYdTYdCAYz = Matmul(Conjg(Yd),TYdCAYz)
   CYdTAYdCAYd = Matmul(Conjg(Yd),TAYdCAYd)
   CYdTAYdCAYs = Matmul(Conjg(Yd),TAYdCAYs)
   CYdTAYdCAYz = Matmul(Conjg(Yd),TAYdCAYz)
   CYeTYeCYe = Matmul(Conjg(Ye),TYeCYe)
   CYeTYeCAYe = Matmul(Conjg(Ye),TYeCAYe)
   CYeTAYeCAYe = Matmul(Conjg(Ye),TAYeCAYe)
   CYsYdadjYd = Matmul(adjYs,YdadjYd)
   CYsYsCYd = Matmul(adjYs,YsCYd)
   CYsYsCYs = Matmul(adjYs,YsCYs)
   CYsYsCYz = Matmul(adjYs,YsCYz)
   CYsYsCAYd = Matmul(adjYs,YsCAYd)
   CYsYsCAYs = Matmul(adjYs,YsCAYs)
   CYsYsCAYz = Matmul(adjYs,YsCAYz)
   CYsYzadjYz = Matmul(adjYs,YzadjYz)
   CYsAYdadjYd = Matmul(adjYs,AYdadjYd)
   CYsAYdadjAYd = Matmul(adjYs,AYdadjAYd)
   CYsAYsCYs = Matmul(adjYs,AYsCYs)
   CYsAYsCAYd = Matmul(adjYs,AYsCAYd)
   CYsAYsCAYs = Matmul(adjYs,AYsCAYs)
   CYsAYsCAYz = Matmul(adjYs,AYsCAYz)
   CYsAYzadjYz = Matmul(adjYs,AYzadjYz)
   CYsAYzadjAYz = Matmul(adjYs,AYzadjAYz)
   CYtYtadjYe = Matmul(adjYt,YtadjYe)
   CYtYtadjYz = Matmul(adjYt,YtadjYz)
   CYtYtadjAYe = Matmul(adjYt,YtadjAYe)
   CYtYtadjAYz = Matmul(adjYt,YtadjAYz)
   CYtYtCYt = Matmul(adjYt,YtCYt)
   CYtYtCAYt = Matmul(adjYt,YtCAYt)
   CYtAYtadjYe = Matmul(adjYt,AYtadjYe)
   CYtAYtadjYz = Matmul(adjYt,AYtadjYz)
   CYtAYtadjAYe = Matmul(adjYt,AYtadjAYe)
   CYtAYtadjAYz = Matmul(adjYt,AYtadjAYz)
   CYtAYtCYt = Matmul(adjYt,AYtCYt)
   CYtAYtCAYt = Matmul(adjYt,AYtCAYt)
   CYuTYuCYu = Matmul(Conjg(Yu),TYuCYu)
   CYuTYuCAYu = Matmul(Conjg(Yu),TYuCAYu)
   CYuTAYuCAYu = Matmul(Conjg(Yu),TAYuCAYu)
   CYzTYzCYd = Matmul(Conjg(Yz),TYzCYd)
   CYzTYzCYs = Matmul(Conjg(Yz),TYzCYs)
   CYzTYzCYz = Matmul(Conjg(Yz),TYzCYz)
   CYzTYzCAYd = Matmul(Conjg(Yz),TYzCAYd)
   CYzTYzCAYs = Matmul(Conjg(Yz),TYzCAYs)
   CYzTYzCAYz = Matmul(Conjg(Yz),TYzCAYz)
   CYzTAYzCAYd = Matmul(Conjg(Yz),TAYzCAYd)
   CYzTAYzCAYs = Matmul(Conjg(Yz),TAYzCAYs)
   CYzTAYzCAYz = Matmul(Conjg(Yz),TAYzCAYz)
   CAYdTYdCYd = Matmul(Conjg(AYd),TYdCYd)
   CAYdTYdCYs = Matmul(Conjg(AYd),TYdCYs)
   CAYdTYdCYz = Matmul(Conjg(AYd),TYdCYz)
   CAYdTAYdCYd = Matmul(Conjg(AYd),TAYdCYd)
   CAYdTAYdCYs = Matmul(Conjg(AYd),TAYdCYs)
   CAYdTAYdCYz = Matmul(Conjg(AYd),TAYdCYz)
   CAYeTYeCYe = Matmul(Conjg(AYe),TYeCYe)
   CAYeTAYeCYe = Matmul(Conjg(AYe),TAYeCYe)
   CAYsYsCYd = Matmul(adjAYs,YsCYd)
   CAYsYsCYs = Matmul(adjAYs,YsCYs)
   CAYsYsCYz = Matmul(adjAYs,YsCYz)
   CAYsAYdadjYd = Matmul(adjAYs,AYdadjYd)
   CAYsAYsCYd = Matmul(adjAYs,AYsCYd)
   CAYsAYsCYs = Matmul(adjAYs,AYsCYs)
   CAYsAYsCYz = Matmul(adjAYs,AYsCYz)
   CAYsAYzadjYz = Matmul(adjAYs,AYzadjYz)
   CAYtYtadjYe = Matmul(adjAYt,YtadjYe)
   CAYtYtadjYz = Matmul(adjAYt,YtadjYz)
   CAYtYtCYt = Matmul(adjAYt,YtCYt)
   CAYtAYtadjYe = Matmul(adjAYt,AYtadjYe)
   CAYtAYtadjYz = Matmul(adjAYt,AYtadjYz)
   CAYtAYtCYt = Matmul(adjAYt,AYtCYt)
   CAYuTYuCYu = Matmul(Conjg(AYu),TYuCYu)
   CAYuTAYuCYu = Matmul(Conjg(AYu),TAYuCYu)
   CAYzTYzCYd = Matmul(Conjg(AYz),TYzCYd)
   CAYzTYzCYs = Matmul(Conjg(AYz),TYzCYs)
   CAYzTYzCYz = Matmul(Conjg(AYz),TYzCYz)
   CAYzTAYzCYd = Matmul(Conjg(AYz),TAYzCYd)
   CAYzTAYzCYs = Matmul(Conjg(AYz),TAYzCYs)
   CAYzTAYzCYz = Matmul(Conjg(AYz),TAYzCYz)
   TYdmd2CYs = Matmul(Transpose(Yd),md2CYs)
   TYdmd2CYz = Matmul(Transpose(Yd),md2CYz)
   TYdCYdTYd = Matmul(Transpose(Yd),CYdTYd)
   TYdCYdTAYd = Matmul(Transpose(Yd),CYdTAYd)
   TYdCYsmd2 = Matmul(Transpose(Yd),CYsmd2)
   TYdCYsYd = Matmul(Transpose(Yd),CYsYd)
   TYdCYsYs = Matmul(Transpose(Yd),CYsYs)
   TYdCYsYz = Matmul(Transpose(Yd),CYsYz)
   TYdCYsAYd = Matmul(Transpose(Yd),CYsAYd)
   TYdCYsAYs = Matmul(Transpose(Yd),CYsAYs)
   TYdCYsAYz = Matmul(Transpose(Yd),CYsAYz)
   TYdCYzml2 = Matmul(Transpose(Yd),CYzml2)
   TYdCYzYt = Matmul(Transpose(Yd),CYzYt)
   TYdCYzAYt = Matmul(Transpose(Yd),CYzAYt)
   TYeCYeTYz = Matmul(Transpose(Ye),CYeTYz)
   TYeCYeTAYz = Matmul(Transpose(Ye),CYeTAYz)
   TYuCYuTYd = Matmul(Transpose(Yu),CYuTYd)
   TYuCYuTAYd = Matmul(Transpose(Yu),CYuTAYd)
   TYzmd2CYd = Matmul(Transpose(Yz),md2CYd)
   TYzmd2CYs = Matmul(Transpose(Yz),md2CYs)
   TYzCYdmq2 = Matmul(Transpose(Yz),CYdmq2)
   TYzCYsmd2 = Matmul(Transpose(Yz),CYsmd2)
   TYzCYsYd = Matmul(Transpose(Yz),CYsYd)
   TYzCYsYs = Matmul(Transpose(Yz),CYsYs)
   TYzCYsYz = Matmul(Transpose(Yz),CYsYz)
   TYzCYsAYd = Matmul(Transpose(Yz),CYsAYd)
   TYzCYsAYs = Matmul(Transpose(Yz),CYsAYs)
   TYzCYsAYz = Matmul(Transpose(Yz),CYsAYz)
   TYzCYzTYz = Matmul(Transpose(Yz),CYzTYz)
   TYzCYzTAYz = Matmul(Transpose(Yz),CYzTAYz)
   TAYdCYdTYd = Matmul(Transpose(AYd),CYdTYd)
   TAYdCYsYd = Matmul(Transpose(AYd),CYsYd)
   TAYdCYsYs = Matmul(Transpose(AYd),CYsYs)
   TAYdCYsYz = Matmul(Transpose(AYd),CYsYz)
   TAYdCYzYt = Matmul(Transpose(AYd),CYzYt)
   TAYeCYeTYz = Matmul(Transpose(AYe),CYeTYz)
   TAYuCYuTYd = Matmul(Transpose(AYu),CYuTYd)
   TAYzCYsYd = Matmul(Transpose(AYz),CYsYd)
   TAYzCYsYs = Matmul(Transpose(AYz),CYsYs)
   TAYzCYsYz = Matmul(Transpose(AYz),CYsYz)
   TAYzCYzTYz = Matmul(Transpose(AYz),CYzTYz)
   md2YsCYsYs = Matmul(md2,YsCYsYs)
   md2CYdTYdCYd = Matmul(md2,CYdTYdCYd)
   md2CYdTYdCYs = Matmul(md2,CYdTYdCYs)
   md2CYdTYdCYz = Matmul(md2,CYdTYdCYz)
   md2CYsYdadjYd = Matmul(md2,CYsYdadjYd)
   md2CYsYsCYd = Matmul(md2,CYsYsCYd)
   md2CYsYsCYs = Matmul(md2,CYsYsCYs)
   md2CYsYsCYz = Matmul(md2,CYsYsCYz)
   md2CYsYzadjYz = Matmul(md2,CYsYzadjYz)
   md2CYzTYzCYd = Matmul(md2,CYzTYzCYd)
   md2CYzTYzCYs = Matmul(md2,CYzTYzCYs)
   md2CYzTYzCYz = Matmul(md2,CYzTYzCYz)
   me2CYeTYeCYe = Matmul(me2,CYeTYeCYe)
   ml2YtCYtYt = Matmul(ml2,YtCYtYt)
   ml2adjYeYeadjYe = Matmul(ml2,adjYeYeadjYe)
   ml2adjYeYeadjYz = Matmul(ml2,adjYeYeadjYz)
   ml2adjYeYeCYt = Matmul(ml2,adjYeYeCYt)
   ml2adjYzYdadjYd = Matmul(ml2,adjYzYdadjYd)
   ml2adjYzYsCYs = Matmul(ml2,adjYzYsCYs)
   ml2adjYzYzadjYe = Matmul(ml2,adjYzYzadjYe)
   ml2adjYzYzadjYz = Matmul(ml2,adjYzYzadjYz)
   ml2adjYzYzCYt = Matmul(ml2,adjYzYzCYt)
   ml2CYtYtadjYe = Matmul(ml2,CYtYtadjYe)
   ml2CYtYtadjYz = Matmul(ml2,CYtYtadjYz)
   ml2CYtYtCYt = Matmul(ml2,CYtYtCYt)
   mq2adjYdYdadjYd = Matmul(mq2,adjYdYdadjYd)
   mq2adjYdYdadjYu = Matmul(mq2,adjYdYdadjYu)
   mq2adjYdYsCYs = Matmul(mq2,adjYdYsCYs)
   mq2adjYdYzadjYz = Matmul(mq2,adjYdYzadjYz)
   mq2adjYuYuadjYd = Matmul(mq2,adjYuYuadjYd)
   mq2adjYuYuadjYu = Matmul(mq2,adjYuYuadjYu)
   mu2CYuTYuCYu = Matmul(mu2,CYuTYuCYu)
   Ydmq2adjYdYs = Matmul(Yd,mq2adjYdYs)
   YdadjYdmd2Ys = Matmul(Yd,adjYdmd2Ys)
   YdadjYdYdadjYd = Matmul(Yd,adjYdYdadjYd)
   YdadjYdYsmd2 = Matmul(Yd,adjYdYsmd2)
   YdadjYdYsCYs = Matmul(Yd,adjYdYsCYs)
   YdadjYdYzadjYz = Matmul(Yd,adjYdYzadjYz)
   YdadjYdAYdadjYd = Matmul(Yd,adjYdAYdadjYd)
   YdadjYdAYdadjAYd = Matmul(Yd,adjYdAYdadjAYd)
   YdadjYdAYsCYs = Matmul(Yd,adjYdAYsCYs)
   YdadjYdAYsCAYs = Matmul(Yd,adjYdAYsCAYs)
   YdadjYdAYzadjYz = Matmul(Yd,adjYdAYzadjYz)
   YdadjYdAYzadjAYz = Matmul(Yd,adjYdAYzadjAYz)
   YdadjYuYuadjYd = Matmul(Yd,adjYuYuadjYd)
   YdadjYuAYuadjYd = Matmul(Yd,adjYuAYuadjYd)
   YdadjYuAYuadjAYd = Matmul(Yd,adjYuAYuadjAYd)
   YdadjAYdAYdadjYd = Matmul(Yd,adjAYdAYdadjYd)
   YdadjAYdAYsCYs = Matmul(Yd,adjAYdAYsCYs)
   YdadjAYdAYzadjYz = Matmul(Yd,adjAYdAYzadjYz)
   YdadjAYuAYuadjYd = Matmul(Yd,adjAYuAYuadjYd)
   YeadjYeYeadjYe = Matmul(Ye,adjYeYeadjYe)
   YeadjYeAYeadjYe = Matmul(Ye,adjYeAYeadjYe)
   YeadjYeAYeadjAYe = Matmul(Ye,adjYeAYeadjAYe)
   YeadjYzYzadjYe = Matmul(Ye,adjYzYzadjYe)
   YeadjYzAYzadjYe = Matmul(Ye,adjYzAYzadjYe)
   YeadjYzAYzadjAYe = Matmul(Ye,adjYzAYzadjAYe)
   YeadjAYeAYeadjYe = Matmul(Ye,adjAYeAYeadjYe)
   YeadjAYzAYzadjYe = Matmul(Ye,adjAYzAYzadjYe)
   YeCYtYtadjYe = Matmul(Ye,CYtYtadjYe)
   YeCYtAYtadjYe = Matmul(Ye,CYtAYtadjYe)
   YeCYtAYtadjAYe = Matmul(Ye,CYtAYtadjAYe)
   YeCAYtAYtadjYe = Matmul(Ye,CAYtAYtadjYe)
   Ysmd2CYsYs = Matmul(Ys,md2CYsYs)
   YsCYdTYdCYs = Matmul(Ys,CYdTYdCYs)
   YsCYdTAYdCAYs = Matmul(Ys,CYdTAYdCAYs)
   YsCYsmd2Ys = Matmul(Ys,CYsmd2Ys)
   YsCYsYdadjYd = Matmul(Ys,CYsYdadjYd)
   YsCYsYsmd2 = Matmul(Ys,CYsYsmd2)
   YsCYsYsCYs = Matmul(Ys,CYsYsCYs)
   YsCYsYzadjYz = Matmul(Ys,CYsYzadjYz)
   YsCYsAYdadjYd = Matmul(Ys,CYsAYdadjYd)
   YsCYsAYdadjAYd = Matmul(Ys,CYsAYdadjAYd)
   YsCYsAYsCYs = Matmul(Ys,CYsAYsCYs)
   YsCYsAYsCAYs = Matmul(Ys,CYsAYsCAYs)
   YsCYsAYzadjYz = Matmul(Ys,CYsAYzadjYz)
   YsCYsAYzadjAYz = Matmul(Ys,CYsAYzadjAYz)
   YsCYzTYzCYs = Matmul(Ys,CYzTYzCYs)
   YsCYzTAYzCAYs = Matmul(Ys,CYzTAYzCAYs)
   YsCAYdTAYdCYs = Matmul(Ys,CAYdTAYdCYs)
   YsCAYsAYdadjYd = Matmul(Ys,CAYsAYdadjYd)
   YsCAYsAYsCYs = Matmul(Ys,CAYsAYsCYs)
   YsCAYsAYzadjYz = Matmul(Ys,CAYsAYzadjYz)
   YsCAYzTAYzCYs = Matmul(Ys,CAYzTAYzCYs)
   Ytml2adjYeYe = Matmul(Yt,ml2adjYeYe)
   Ytml2adjYzYz = Matmul(Yt,ml2adjYzYz)
   Ytml2CYtYt = Matmul(Yt,ml2CYtYt)
   YtadjYeme2Ye = Matmul(Yt,adjYeme2Ye)
   YtadjYeYeml2 = Matmul(Yt,adjYeYeml2)
   YtadjYeYeCYt = Matmul(Yt,adjYeYeCYt)
   YtadjYeAYeCYt = Matmul(Yt,adjYeAYeCYt)
   YtadjYeAYeCAYt = Matmul(Yt,adjYeAYeCAYt)
   YtadjYzmd2Yz = Matmul(Yt,adjYzmd2Yz)
   YtadjYzYzml2 = Matmul(Yt,adjYzYzml2)
   YtadjYzYzCYt = Matmul(Yt,adjYzYzCYt)
   YtadjYzAYzCYt = Matmul(Yt,adjYzAYzCYt)
   YtadjYzAYzCAYt = Matmul(Yt,adjYzAYzCAYt)
   YtadjAYeAYeCYt = Matmul(Yt,adjAYeAYeCYt)
   YtadjAYzAYzCYt = Matmul(Yt,adjAYzAYzCYt)
   YtCYtml2Yt = Matmul(Yt,CYtml2Yt)
   YtCYtYtml2 = Matmul(Yt,CYtYtml2)
   YtCYtYtCYt = Matmul(Yt,CYtYtCYt)
   YtCYtAYtCYt = Matmul(Yt,CYtAYtCYt)
   YtCYtAYtCAYt = Matmul(Yt,CYtAYtCAYt)
   YtCAYtAYtCYt = Matmul(Yt,CAYtAYtCYt)
   YuadjYdYdadjYu = Matmul(Yu,adjYdYdadjYu)
   YuadjYdAYdadjYu = Matmul(Yu,adjYdAYdadjYu)
   YuadjYdAYdadjAYu = Matmul(Yu,adjYdAYdadjAYu)
   YuadjYuYuadjYu = Matmul(Yu,adjYuYuadjYu)
   YuadjYuAYuadjYu = Matmul(Yu,adjYuAYuadjYu)
   YuadjYuAYuadjAYu = Matmul(Yu,adjYuAYuadjAYu)
   YuadjAYdAYdadjYu = Matmul(Yu,adjAYdAYdadjYu)
   YuadjAYuAYuadjYu = Matmul(Yu,adjAYuAYuadjYu)
   Yzml2adjYzYs = Matmul(Yz,ml2adjYzYs)
   YzadjYeYeadjYz = Matmul(Yz,adjYeYeadjYz)
   YzadjYeAYeadjYz = Matmul(Yz,adjYeAYeadjYz)
   YzadjYeAYeadjAYz = Matmul(Yz,adjYeAYeadjAYz)
   YzadjYzmd2Ys = Matmul(Yz,adjYzmd2Ys)
   YzadjYzYdadjYd = Matmul(Yz,adjYzYdadjYd)
   YzadjYzYsmd2 = Matmul(Yz,adjYzYsmd2)
   YzadjYzYsCYs = Matmul(Yz,adjYzYsCYs)
   YzadjYzYzadjYz = Matmul(Yz,adjYzYzadjYz)
   YzadjYzAYdadjYd = Matmul(Yz,adjYzAYdadjYd)
   YzadjYzAYdadjAYd = Matmul(Yz,adjYzAYdadjAYd)
   YzadjYzAYsCYs = Matmul(Yz,adjYzAYsCYs)
   YzadjYzAYsCAYs = Matmul(Yz,adjYzAYsCAYs)
   YzadjYzAYzadjYz = Matmul(Yz,adjYzAYzadjYz)
   YzadjYzAYzadjAYz = Matmul(Yz,adjYzAYzadjAYz)
   YzadjAYeAYeadjYz = Matmul(Yz,adjAYeAYeadjYz)
   YzadjAYzAYdadjYd = Matmul(Yz,adjAYzAYdadjYd)
   YzadjAYzAYsCYs = Matmul(Yz,adjAYzAYsCYs)
   YzadjAYzAYzadjYz = Matmul(Yz,adjAYzAYzadjYz)
   YzCYtYtadjYz = Matmul(Yz,CYtYtadjYz)
   YzCYtAYtadjYz = Matmul(Yz,CYtAYtadjYz)
   YzCYtAYtadjAYz = Matmul(Yz,CYtAYtadjAYz)
   YzCAYtAYtadjYz = Matmul(Yz,CAYtAYtadjYz)
   AYdadjYdYdadjAYd = Matmul(AYd,adjYdYdadjAYd)
   AYdadjYuYuadjAYd = Matmul(AYd,adjYuYuadjAYd)
   AYdadjAYdYdadjYd = Matmul(AYd,adjAYdYdadjYd)
   AYdadjAYuYuadjYd = Matmul(AYd,adjAYuYuadjYd)
   AYeadjYeYeadjAYe = Matmul(AYe,adjYeYeadjAYe)
   AYeadjYzYzadjAYe = Matmul(AYe,adjYzYzadjAYe)
   AYeadjAYeYeadjYe = Matmul(AYe,adjAYeYeadjYe)
   AYeadjAYzYzadjYe = Matmul(AYe,adjAYzYzadjYe)
   AYeCYtYtadjAYe = Matmul(AYe,CYtYtadjAYe)
   AYeCAYtYtadjYe = Matmul(AYe,CAYtYtadjYe)
   AYsCYdTYdCAYs = Matmul(AYs,CYdTYdCAYs)
   AYsCYsYsCAYs = Matmul(AYs,CYsYsCAYs)
   AYsCYzTYzCAYs = Matmul(AYs,CYzTYzCAYs)
   AYsCAYdTYdCYs = Matmul(AYs,CAYdTYdCYs)
   AYsCAYsYsCYs = Matmul(AYs,CAYsYsCYs)
   AYsCAYzTYzCYs = Matmul(AYs,CAYzTYzCYs)
   AYtadjYeYeCAYt = Matmul(AYt,adjYeYeCAYt)
   AYtadjYzYzCAYt = Matmul(AYt,adjYzYzCAYt)
   AYtadjAYeYeCYt = Matmul(AYt,adjAYeYeCYt)
   AYtadjAYzYzCYt = Matmul(AYt,adjAYzYzCYt)
   AYtCYtYtCAYt = Matmul(AYt,CYtYtCAYt)
   AYtCAYtYtCYt = Matmul(AYt,CAYtYtCYt)
   AYuadjYdYdadjAYu = Matmul(AYu,adjYdYdadjAYu)
   AYuadjYuYuadjAYu = Matmul(AYu,adjYuYuadjAYu)
   AYuadjAYdYdadjYu = Matmul(AYu,adjAYdYdadjYu)
   AYuadjAYuYuadjYu = Matmul(AYu,adjAYuYuadjYu)
   AYzadjYeYeadjAYz = Matmul(AYz,adjYeYeadjAYz)
   AYzadjYzYzadjAYz = Matmul(AYz,adjYzYzadjAYz)
   AYzadjAYeYeadjYz = Matmul(AYz,adjAYeYeadjYz)
   AYzadjAYzYzadjYz = Matmul(AYz,adjAYzYzadjYz)
   AYzCYtYtadjAYz = Matmul(AYz,CYtYtadjAYz)
   AYzCAYtYtadjYz = Matmul(AYz,CAYtYtadjYz)
   adjYdmd2YdadjYd = Matmul(adjYd,md2YdadjYd)
   adjYdmd2YdadjYu = Matmul(adjYd,md2YdadjYu)
   adjYdmd2YsCYs = Matmul(adjYd,md2YsCYs)
   adjYdmd2YzadjYz = Matmul(adjYd,md2YzadjYz)
   adjYdYdmq2adjYd = Matmul(adjYd,Ydmq2adjYd)
   adjYdYdmq2adjYu = Matmul(adjYd,Ydmq2adjYu)
   adjYdYdadjYdmd2 = Matmul(adjYd,YdadjYdmd2)
   adjYdYdadjYdYd = Matmul(adjYd,YdadjYdYd)
   adjYdYdadjYdYs = Matmul(adjYd,YdadjYdYs)
   adjYdYdadjYdYz = Matmul(adjYd,YdadjYdYz)
   adjYdYdadjYdAYd = Matmul(adjYd,YdadjYdAYd)
   adjYdYdadjYdAYs = Matmul(adjYd,YdadjYdAYs)
   adjYdYdadjYdAYz = Matmul(adjYd,YdadjYdAYz)
   adjYdYdadjYumu2 = Matmul(adjYd,YdadjYumu2)
   adjYdYdadjYuYu = Matmul(adjYd,YdadjYuYu)
   adjYdYdadjYuAYu = Matmul(adjYd,YdadjYuAYu)
   adjYdYsCYsYd = Matmul(adjYd,YsCYsYd)
   adjYdYsCYsAYd = Matmul(adjYd,YsCYsAYd)
   adjYdYzml2adjYz = Matmul(adjYd,Yzml2adjYz)
   adjYdYzadjYzYd = Matmul(adjYd,YzadjYzYd)
   adjYdYzadjYzAYd = Matmul(adjYd,YzadjYzAYd)
   adjYdAYdadjYdYd = Matmul(adjYd,AYdadjYdYd)
   adjYdAYdadjYdYs = Matmul(adjYd,AYdadjYdYs)
   adjYdAYdadjYdYz = Matmul(adjYd,AYdadjYdYz)
   adjYdAYdadjYuYu = Matmul(adjYd,AYdadjYuYu)
   adjYdAYsCYsYd = Matmul(adjYd,AYsCYsYd)
   adjYdAYzadjYzYd = Matmul(adjYd,AYzadjYzYd)
   adjYeme2YeadjYe = Matmul(adjYe,me2YeadjYe)
   adjYeme2YeadjYz = Matmul(adjYe,me2YeadjYz)
   adjYeme2YeCYt = Matmul(adjYe,me2YeCYt)
   adjYeYeml2adjYe = Matmul(adjYe,Yeml2adjYe)
   adjYeYeml2adjYz = Matmul(adjYe,Yeml2adjYz)
   adjYeYeml2CYt = Matmul(adjYe,Yeml2CYt)
   adjYeYeadjYeme2 = Matmul(adjYe,YeadjYeme2)
   adjYeYeadjYeYe = Matmul(adjYe,YeadjYeYe)
   adjYeYeadjYeAYe = Matmul(adjYe,YeadjYeAYe)
   adjYeYeadjYzmd2 = Matmul(adjYe,YeadjYzmd2)
   adjYeYeadjYzYd = Matmul(adjYe,YeadjYzYd)
   adjYeYeadjYzYs = Matmul(adjYe,YeadjYzYs)
   adjYeYeadjYzYz = Matmul(adjYe,YeadjYzYz)
   adjYeYeadjYzAYd = Matmul(adjYe,YeadjYzAYd)
   adjYeYeadjYzAYs = Matmul(adjYe,YeadjYzAYs)
   adjYeYeadjYzAYz = Matmul(adjYe,YeadjYzAYz)
   adjYeYeCYtml2 = Matmul(adjYe,YeCYtml2)
   adjYeYeCYtYt = Matmul(adjYe,YeCYtYt)
   adjYeYeCYtAYt = Matmul(adjYe,YeCYtAYt)
   adjYeAYeadjYeYe = Matmul(adjYe,AYeadjYeYe)
   adjYeAYeadjYzYd = Matmul(adjYe,AYeadjYzYd)
   adjYeAYeadjYzYs = Matmul(adjYe,AYeadjYzYs)
   adjYeAYeadjYzYz = Matmul(adjYe,AYeadjYzYz)
   adjYeAYeCYtYt = Matmul(adjYe,AYeCYtYt)
   adjYumu2YuadjYd = Matmul(adjYu,mu2YuadjYd)
   adjYumu2YuadjYu = Matmul(adjYu,mu2YuadjYu)
   adjYuYumq2adjYd = Matmul(adjYu,Yumq2adjYd)
   adjYuYumq2adjYu = Matmul(adjYu,Yumq2adjYu)
   adjYuYuadjYdmd2 = Matmul(adjYu,YuadjYdmd2)
   adjYuYuadjYdYd = Matmul(adjYu,YuadjYdYd)
   adjYuYuadjYdYs = Matmul(adjYu,YuadjYdYs)
   adjYuYuadjYdYz = Matmul(adjYu,YuadjYdYz)
   adjYuYuadjYdAYd = Matmul(adjYu,YuadjYdAYd)
   adjYuYuadjYdAYs = Matmul(adjYu,YuadjYdAYs)
   adjYuYuadjYdAYz = Matmul(adjYu,YuadjYdAYz)
   adjYuYuadjYumu2 = Matmul(adjYu,YuadjYumu2)
   adjYuYuadjYuYu = Matmul(adjYu,YuadjYuYu)
   adjYuYuadjYuAYu = Matmul(adjYu,YuadjYuAYu)
   adjYuAYuadjYdYd = Matmul(adjYu,AYuadjYdYd)
   adjYuAYuadjYdYs = Matmul(adjYu,AYuadjYdYs)
   adjYuAYuadjYdYz = Matmul(adjYu,AYuadjYdYz)
   adjYuAYuadjYuYu = Matmul(adjYu,AYuadjYuYu)
   adjYzmd2YdadjYd = Matmul(adjYz,md2YdadjYd)
   adjYzmd2YsCYs = Matmul(adjYz,md2YsCYs)
   adjYzmd2YzadjYe = Matmul(adjYz,md2YzadjYe)
   adjYzmd2YzadjYz = Matmul(adjYz,md2YzadjYz)
   adjYzmd2YzCYt = Matmul(adjYz,md2YzCYt)
   adjYzYdmq2adjYd = Matmul(adjYz,Ydmq2adjYd)
   adjYzYdadjYdYz = Matmul(adjYz,YdadjYdYz)
   adjYzYdadjYdAYz = Matmul(adjYz,YdadjYdAYz)
   adjYzYsCYsYz = Matmul(adjYz,YsCYsYz)
   adjYzYsCYsAYz = Matmul(adjYz,YsCYsAYz)
   adjYzYzml2adjYe = Matmul(adjYz,Yzml2adjYe)
   adjYzYzml2adjYz = Matmul(adjYz,Yzml2adjYz)
   adjYzYzml2CYt = Matmul(adjYz,Yzml2CYt)
   adjYzYzadjYeme2 = Matmul(adjYz,YzadjYeme2)
   adjYzYzadjYeYe = Matmul(adjYz,YzadjYeYe)
   adjYzYzadjYeAYe = Matmul(adjYz,YzadjYeAYe)
   adjYzYzadjYzmd2 = Matmul(adjYz,YzadjYzmd2)
   adjYzYzadjYzYd = Matmul(adjYz,YzadjYzYd)
   adjYzYzadjYzYs = Matmul(adjYz,YzadjYzYs)
   adjYzYzadjYzYz = Matmul(adjYz,YzadjYzYz)
   adjYzYzadjYzAYd = Matmul(adjYz,YzadjYzAYd)
   adjYzYzadjYzAYs = Matmul(adjYz,YzadjYzAYs)
   adjYzYzadjYzAYz = Matmul(adjYz,YzadjYzAYz)
   adjYzYzCYtml2 = Matmul(adjYz,YzCYtml2)
   adjYzYzCYtYt = Matmul(adjYz,YzCYtYt)
   adjYzYzCYtAYt = Matmul(adjYz,YzCYtAYt)
   adjYzAYdadjYdYz = Matmul(adjYz,AYdadjYdYz)
   adjYzAYsCYsYz = Matmul(adjYz,AYsCYsYz)
   adjYzAYzadjYeYe = Matmul(adjYz,AYzadjYeYe)
   adjYzAYzadjYzYd = Matmul(adjYz,AYzadjYzYd)
   adjYzAYzadjYzYs = Matmul(adjYz,AYzadjYzYs)
   adjYzAYzadjYzYz = Matmul(adjYz,AYzadjYzYz)
   adjYzAYzCYtYt = Matmul(adjYz,AYzCYtYt)
   CYdmq2TYdCYd = Matmul(Conjg(Yd),mq2TYdCYd)
   CYdmq2TYdCYs = Matmul(Conjg(Yd),mq2TYdCYs)
   CYdmq2TYdCYz = Matmul(Conjg(Yd),mq2TYdCYz)
   CYdTYdmd2CYd = Matmul(Conjg(Yd),TYdmd2CYd)
   CYdTYdmd2CYs = Matmul(Conjg(Yd),TYdmd2CYs)
   CYdTYdmd2CYz = Matmul(Conjg(Yd),TYdmd2CYz)
   CYdTYdCYdmq2 = Matmul(Conjg(Yd),TYdCYdmq2)
   CYdTYdCYdTYd = Matmul(Conjg(Yd),TYdCYdTYd)
   CYdTYdCYdTAYd = Matmul(Conjg(Yd),TYdCYdTAYd)
   CYdTYdCYsmd2 = Matmul(Conjg(Yd),TYdCYsmd2)
   CYdTYdCYsYd = Matmul(Conjg(Yd),TYdCYsYd)
   CYdTYdCYsYs = Matmul(Conjg(Yd),TYdCYsYs)
   CYdTYdCYsYz = Matmul(Conjg(Yd),TYdCYsYz)
   CYdTYdCYsAYd = Matmul(Conjg(Yd),TYdCYsAYd)
   CYdTYdCYsAYs = Matmul(Conjg(Yd),TYdCYsAYs)
   CYdTYdCYsAYz = Matmul(Conjg(Yd),TYdCYsAYz)
   CYdTYdCYzml2 = Matmul(Conjg(Yd),TYdCYzml2)
   CYdTYdCYzYt = Matmul(Conjg(Yd),TYdCYzYt)
   CYdTYdCYzAYt = Matmul(Conjg(Yd),TYdCYzAYt)
   CYdTYuCYuTYd = Matmul(Conjg(Yd),TYuCYuTYd)
   CYdTYuCYuTAYd = Matmul(Conjg(Yd),TYuCYuTAYd)
   CYdTAYdCYdTYd = Matmul(Conjg(Yd),TAYdCYdTYd)
   CYdTAYdCYsYd = Matmul(Conjg(Yd),TAYdCYsYd)
   CYdTAYdCYsYs = Matmul(Conjg(Yd),TAYdCYsYs)
   CYdTAYdCYsYz = Matmul(Conjg(Yd),TAYdCYsYz)
   CYdTAYdCYzYt = Matmul(Conjg(Yd),TAYdCYzYt)
   CYdTAYuCYuTYd = Matmul(Conjg(Yd),TAYuCYuTYd)
   CYeml2TYeCYe = Matmul(Conjg(Ye),ml2TYeCYe)
   CYeTYeme2CYe = Matmul(Conjg(Ye),TYeme2CYe)
   CYeTYeCYeml2 = Matmul(Conjg(Ye),TYeCYeml2)
   CYeTYeCYeYt = Matmul(Conjg(Ye),TYeCYeYt)
   CYeTYeCYeAYt = Matmul(Conjg(Ye),TYeCYeAYt)
   CYeTAYeCYeYt = Matmul(Conjg(Ye),TAYeCYeYt)
   CYsmd2YdadjYd = Matmul(adjYs,md2YdadjYd)
   CYsmd2YsCYd = Matmul(adjYs,md2YsCYd)
   CYsmd2YsCYs = Matmul(adjYs,md2YsCYs)
   CYsmd2YsCYz = Matmul(adjYs,md2YsCYz)
   CYsmd2YzadjYz = Matmul(adjYs,md2YzadjYz)
   CYsYdmq2adjYd = Matmul(adjYs,Ydmq2adjYd)
   CYsYdadjYdmd2 = Matmul(adjYs,YdadjYdmd2)
   CYsYdadjYdYs = Matmul(adjYs,YdadjYdYs)
   CYsYdadjYdAYs = Matmul(adjYs,YdadjYdAYs)
   CYsYdadjAYdAYs = Matmul(adjYs,YdadjAYdAYs)
   CYsYsmd2CYd = Matmul(adjYs,Ysmd2CYd)
   CYsYsmd2CYs = Matmul(adjYs,Ysmd2CYs)
   CYsYsmd2CYz = Matmul(adjYs,Ysmd2CYz)
   CYsYsCYdmq2 = Matmul(adjYs,YsCYdmq2)
   CYsYsCYsmd2 = Matmul(adjYs,YsCYsmd2)
   CYsYsCYsYd = Matmul(adjYs,YsCYsYd)
   CYsYsCYsYs = Matmul(adjYs,YsCYsYs)
   CYsYsCYsYz = Matmul(adjYs,YsCYsYz)
   CYsYsCYsAYd = Matmul(adjYs,YsCYsAYd)
   CYsYsCYsAYs = Matmul(adjYs,YsCYsAYs)
   CYsYsCYsAYz = Matmul(adjYs,YsCYsAYz)
   CYsYsCYzml2 = Matmul(adjYs,YsCYzml2)
   CYsYsCYzYt = Matmul(adjYs,YsCYzYt)
   CYsYsCYzAYt = Matmul(adjYs,YsCYzAYt)
   CYsYsCAYsAYs = Matmul(adjYs,YsCAYsAYs)
   CYsYzml2adjYz = Matmul(adjYs,Yzml2adjYz)
   CYsYzadjYzmd2 = Matmul(adjYs,YzadjYzmd2)
   CYsYzadjYzYs = Matmul(adjYs,YzadjYzYs)
   CYsYzadjYzAYs = Matmul(adjYs,YzadjYzAYs)
   CYsYzadjAYzAYs = Matmul(adjYs,YzadjAYzAYs)
   CYsAYdadjYdYs = Matmul(adjYs,AYdadjYdYs)
   CYsAYdadjAYdYs = Matmul(adjYs,AYdadjAYdYs)
   CYsAYsCYsYd = Matmul(adjYs,AYsCYsYd)
   CYsAYsCYsYs = Matmul(adjYs,AYsCYsYs)
   CYsAYsCYsYz = Matmul(adjYs,AYsCYsYz)
   CYsAYsCYzYt = Matmul(adjYs,AYsCYzYt)
   CYsAYsCAYsYs = Matmul(adjYs,AYsCAYsYs)
   CYsAYzadjYzYs = Matmul(adjYs,AYzadjYzYs)
   CYsAYzadjAYzYs = Matmul(adjYs,AYzadjAYzYs)
   CYtml2YtadjYe = Matmul(adjYt,ml2YtadjYe)
   CYtml2YtadjYz = Matmul(adjYt,ml2YtadjYz)
   CYtml2YtCYt = Matmul(adjYt,ml2YtCYt)
   CYtYtml2adjYe = Matmul(adjYt,Ytml2adjYe)
   CYtYtml2adjYz = Matmul(adjYt,Ytml2adjYz)
   CYtYtml2CYt = Matmul(adjYt,Ytml2CYt)
   CYtYtadjYeme2 = Matmul(adjYt,YtadjYeme2)
   CYtYtadjYeYe = Matmul(adjYt,YtadjYeYe)
   CYtYtadjYeAYe = Matmul(adjYt,YtadjYeAYe)
   CYtYtadjYzmd2 = Matmul(adjYt,YtadjYzmd2)
   CYtYtadjYzYd = Matmul(adjYt,YtadjYzYd)
   CYtYtadjYzYs = Matmul(adjYt,YtadjYzYs)
   CYtYtadjYzYz = Matmul(adjYt,YtadjYzYz)
   CYtYtadjYzAYd = Matmul(adjYt,YtadjYzAYd)
   CYtYtadjYzAYs = Matmul(adjYt,YtadjYzAYs)
   CYtYtadjYzAYz = Matmul(adjYt,YtadjYzAYz)
   CYtYtadjAYeAYe = Matmul(adjYt,YtadjAYeAYe)
   CYtYtadjAYzAYz = Matmul(adjYt,YtadjAYzAYz)
   CYtYtCYtml2 = Matmul(adjYt,YtCYtml2)
   CYtYtCYtYt = Matmul(adjYt,YtCYtYt)
   CYtYtCYtAYt = Matmul(adjYt,YtCYtAYt)
   CYtYtCAYtAYt = Matmul(adjYt,YtCAYtAYt)
   CYtAYtadjYeYe = Matmul(adjYt,AYtadjYeYe)
   CYtAYtadjYzYd = Matmul(adjYt,AYtadjYzYd)
   CYtAYtadjYzYs = Matmul(adjYt,AYtadjYzYs)
   CYtAYtadjYzYz = Matmul(adjYt,AYtadjYzYz)
   CYtAYtadjAYeYe = Matmul(adjYt,AYtadjAYeYe)
   CYtAYtadjAYzYz = Matmul(adjYt,AYtadjAYzYz)
   CYtAYtCYtYt = Matmul(adjYt,AYtCYtYt)
   CYtAYtCAYtYt = Matmul(adjYt,AYtCAYtYt)
   CYtTYeCYeYt = Matmul(adjYt,TYeCYeYt)
   CYtTYeCYeAYt = Matmul(adjYt,TYeCYeAYt)
   CYtTYzCYzYt = Matmul(adjYt,TYzCYzYt)
   CYtTYzCYzAYt = Matmul(adjYt,TYzCYzAYt)
   CYtTAYeCYeYt = Matmul(adjYt,TAYeCYeYt)
   CYtTAYzCYzYt = Matmul(adjYt,TAYzCYzYt)
   CYumq2TYuCYu = Matmul(Conjg(Yu),mq2TYuCYu)
   CYuTYumu2CYu = Matmul(Conjg(Yu),TYumu2CYu)
   CYuTYuCYumq2 = Matmul(Conjg(Yu),TYuCYumq2)
   CYzml2TYzCYd = Matmul(Conjg(Yz),ml2TYzCYd)
   CYzml2TYzCYs = Matmul(Conjg(Yz),ml2TYzCYs)
   CYzml2TYzCYz = Matmul(Conjg(Yz),ml2TYzCYz)
   CYzYtCYtTYz = Matmul(Conjg(Yz),YtCYtTYz)
   CYzYtCYtTAYz = Matmul(Conjg(Yz),YtCYtTAYz)
   CYzAYtCYtTYz = Matmul(Conjg(Yz),AYtCYtTYz)
   CYzTYeCYeTYz = Matmul(Conjg(Yz),TYeCYeTYz)
   CYzTYeCYeTAYz = Matmul(Conjg(Yz),TYeCYeTAYz)
   CYzTYzmd2CYd = Matmul(Conjg(Yz),TYzmd2CYd)
   CYzTYzmd2CYs = Matmul(Conjg(Yz),TYzmd2CYs)
   CYzTYzmd2CYz = Matmul(Conjg(Yz),TYzmd2CYz)
   CYzTYzCYdmq2 = Matmul(Conjg(Yz),TYzCYdmq2)
   CYzTYzCYsmd2 = Matmul(Conjg(Yz),TYzCYsmd2)
   CYzTYzCYsYd = Matmul(Conjg(Yz),TYzCYsYd)
   CYzTYzCYsYs = Matmul(Conjg(Yz),TYzCYsYs)
   CYzTYzCYsYz = Matmul(Conjg(Yz),TYzCYsYz)
   CYzTYzCYsAYd = Matmul(Conjg(Yz),TYzCYsAYd)
   CYzTYzCYsAYs = Matmul(Conjg(Yz),TYzCYsAYs)
   CYzTYzCYsAYz = Matmul(Conjg(Yz),TYzCYsAYz)
   CYzTYzCYzml2 = Matmul(Conjg(Yz),TYzCYzml2)
   CYzTYzCYzYt = Matmul(Conjg(Yz),TYzCYzYt)
   CYzTYzCYzAYt = Matmul(Conjg(Yz),TYzCYzAYt)
   CYzTYzCYzTYz = Matmul(Conjg(Yz),TYzCYzTYz)
   CYzTYzCYzTAYz = Matmul(Conjg(Yz),TYzCYzTAYz)
   CYzTAYeCYeTYz = Matmul(Conjg(Yz),TAYeCYeTYz)
   CYzTAYzCYsYd = Matmul(Conjg(Yz),TAYzCYsYd)
   CYzTAYzCYsYs = Matmul(Conjg(Yz),TAYzCYsYs)
   CYzTAYzCYsYz = Matmul(Conjg(Yz),TAYzCYsYz)
   CYzTAYzCYzYt = Matmul(Conjg(Yz),TAYzCYzYt)
   CYzTAYzCYzTYz = Matmul(Conjg(Yz),TAYzCYzTYz)
   TYdCYdTYdCYd = Matmul(Transpose(Yd),CYdTYdCYd)
   TYdCYdTAYdCAYd = Matmul(Transpose(Yd),CYdTAYdCAYd)
   TYdCYsYsCYd = Matmul(Transpose(Yd),CYsYsCYd)
   TYdCYsAYsCAYd = Matmul(Transpose(Yd),CYsAYsCAYd)
   TYdCYzTYzCYd = Matmul(Transpose(Yd),CYzTYzCYd)
   TYdCYzTAYzCAYd = Matmul(Transpose(Yd),CYzTAYzCAYd)
   TYdCAYdTAYdCYd = Matmul(Transpose(Yd),CAYdTAYdCYd)
   TYdCAYsAYsCYd = Matmul(Transpose(Yd),CAYsAYsCYd)
   TYdCAYzTAYzCYd = Matmul(Transpose(Yd),CAYzTAYzCYd)
   TYeCYeTYeCYe = Matmul(Transpose(Ye),CYeTYeCYe)
   TYeCYeTAYeCAYe = Matmul(Transpose(Ye),CYeTAYeCAYe)
   TYeCAYeTAYeCYe = Matmul(Transpose(Ye),CAYeTAYeCYe)
   TYuCYuTYuCYu = Matmul(Transpose(Yu),CYuTYuCYu)
   TYuCYuTAYuCAYu = Matmul(Transpose(Yu),CYuTAYuCAYu)
   TYuCAYuTAYuCYu = Matmul(Transpose(Yu),CAYuTAYuCYu)
   TYzCYdTYdCYz = Matmul(Transpose(Yz),CYdTYdCYz)
   TYzCYdTAYdCAYz = Matmul(Transpose(Yz),CYdTAYdCAYz)
   TYzCYsYsCYz = Matmul(Transpose(Yz),CYsYsCYz)
   TYzCYsAYsCAYz = Matmul(Transpose(Yz),CYsAYsCAYz)
   TYzCYzTYzCYz = Matmul(Transpose(Yz),CYzTYzCYz)
   TYzCYzTAYzCAYz = Matmul(Transpose(Yz),CYzTAYzCAYz)
   TYzCAYdTAYdCYz = Matmul(Transpose(Yz),CAYdTAYdCYz)
   TYzCAYsAYsCYz = Matmul(Transpose(Yz),CAYsAYsCYz)
   TYzCAYzTAYzCYz = Matmul(Transpose(Yz),CAYzTAYzCYz)
   TAYdCYdTYdCAYd = Matmul(Transpose(AYd),CYdTYdCAYd)
   TAYdCYsYsCAYd = Matmul(Transpose(AYd),CYsYsCAYd)
   TAYdCYzTYzCAYd = Matmul(Transpose(AYd),CYzTYzCAYd)
   TAYdCAYdTYdCYd = Matmul(Transpose(AYd),CAYdTYdCYd)
   TAYdCAYsYsCYd = Matmul(Transpose(AYd),CAYsYsCYd)
   TAYdCAYzTYzCYd = Matmul(Transpose(AYd),CAYzTYzCYd)
   TAYeCYeTYeCAYe = Matmul(Transpose(AYe),CYeTYeCAYe)
   TAYeCAYeTYeCYe = Matmul(Transpose(AYe),CAYeTYeCYe)
   TAYuCYuTYuCAYu = Matmul(Transpose(AYu),CYuTYuCAYu)
   TAYuCAYuTYuCYu = Matmul(Transpose(AYu),CAYuTYuCYu)
   TAYzCYdTYdCAYz = Matmul(Transpose(AYz),CYdTYdCAYz)
   TAYzCYsYsCAYz = Matmul(Transpose(AYz),CYsYsCAYz)
   TAYzCYzTYzCAYz = Matmul(Transpose(AYz),CYzTYzCAYz)
   TAYzCAYdTYdCYz = Matmul(Transpose(AYz),CAYdTYdCYz)
   TAYzCAYsYsCYz = Matmul(Transpose(AYz),CAYsYsCYz)
   TAYzCAYzTYzCYz = Matmul(Transpose(AYz),CAYzTYzCYz)
   md2YdadjYdYdadjYd = Matmul(md2,YdadjYdYdadjYd)
   md2YdadjYdYsCYs = Matmul(md2,YdadjYdYsCYs)
   md2YdadjYuYuadjYd = Matmul(md2,YdadjYuYuadjYd)
   md2YsCYdTYdCYs = Matmul(md2,YsCYdTYdCYs)
   md2YsCYsYdadjYd = Matmul(md2,YsCYsYdadjYd)
   md2YsCYsYsCYs = Matmul(md2,YsCYsYsCYs)
   md2YsCYsYzadjYz = Matmul(md2,YsCYsYzadjYz)
   md2YsCYzTYzCYs = Matmul(md2,YsCYzTYzCYs)
   md2YzadjYeYeadjYz = Matmul(md2,YzadjYeYeadjYz)
   md2YzadjYzYsCYs = Matmul(md2,YzadjYzYsCYs)
   md2YzadjYzYzadjYz = Matmul(md2,YzadjYzYzadjYz)
   md2YzCYtYtadjYz = Matmul(md2,YzCYtYtadjYz)
   md2CYsYsCYsYs = Matmul(md2,CYsYsCYsYs)
   me2YeadjYeYeadjYe = Matmul(me2,YeadjYeYeadjYe)
   me2YeadjYzYzadjYe = Matmul(me2,YeadjYzYzadjYe)
   me2YeCYtYtadjYe = Matmul(me2,YeCYtYtadjYe)
   ml2YtadjYeYeCYt = Matmul(ml2,YtadjYeYeCYt)
   ml2YtadjYzYzCYt = Matmul(ml2,YtadjYzYzCYt)
   ml2YtCYtYtCYt = Matmul(ml2,YtCYtYtCYt)
   ml2adjYzYsCYsYz = Matmul(ml2,adjYzYsCYsYz)
   ml2CYtYtCYtYt = Matmul(ml2,CYtYtCYtYt)
   ml2TYeCYeTYeCYe = Matmul(ml2,TYeCYeTYeCYe)
   ml2TYzCYdTYdCYz = Matmul(ml2,TYzCYdTYdCYz)
   ml2TYzCYsYsCYz = Matmul(ml2,TYzCYsYsCYz)
   ml2TYzCYzTYzCYz = Matmul(ml2,TYzCYzTYzCYz)
   mq2adjYdYsCYsYd = Matmul(mq2,adjYdYsCYsYd)
   mq2TYdCYdTYdCYd = Matmul(mq2,TYdCYdTYdCYd)
   mq2TYdCYsYsCYd = Matmul(mq2,TYdCYsYsCYd)
   mq2TYdCYzTYzCYd = Matmul(mq2,TYdCYzTYzCYd)
   mq2TYuCYuTYuCYu = Matmul(mq2,TYuCYuTYuCYu)
   mu2YuadjYdYdadjYu = Matmul(mu2,YuadjYdYdadjYu)
   mu2YuadjYuYuadjYu = Matmul(mu2,YuadjYuYuadjYu)
   Ydmq2adjYdYdadjYd = Matmul(Yd,mq2adjYdYdadjYd)
   Ydmq2adjYdYsCYs = Matmul(Yd,mq2adjYdYsCYs)
   Ydmq2adjYdYzadjYz = Matmul(Yd,mq2adjYdYzadjYz)
   Ydmq2adjYuYuadjYd = Matmul(Yd,mq2adjYuYuadjYd)
   YdadjYdmd2YdadjYd = Matmul(Yd,adjYdmd2YdadjYd)
   YdadjYdmd2YsCYs = Matmul(Yd,adjYdmd2YsCYs)
   YdadjYdmd2YzadjYz = Matmul(Yd,adjYdmd2YzadjYz)
   YdadjYdYdmq2adjYd = Matmul(Yd,adjYdYdmq2adjYd)
   YdadjYdYdadjYdmd2 = Matmul(Yd,adjYdYdadjYdmd2)
   YdadjYdYdadjYdYd = Matmul(Yd,adjYdYdadjYdYd)
   YdadjYdYdadjYdYs = Matmul(Yd,adjYdYdadjYdYs)
   YdadjYdYdadjYdYz = Matmul(Yd,adjYdYdadjYdYz)
   YdadjYdYdadjYdAYd = Matmul(Yd,adjYdYdadjYdAYd)
   YdadjYdYdadjYdAYs = Matmul(Yd,adjYdYdadjYdAYs)
   YdadjYdYdadjYdAYz = Matmul(Yd,adjYdYdadjYdAYz)
   YdadjYdYsCYsYd = Matmul(Yd,adjYdYsCYsYd)
   YdadjYdYsCYsAYd = Matmul(Yd,adjYdYsCYsAYd)
   YdadjYdYzml2adjYz = Matmul(Yd,adjYdYzml2adjYz)
   YdadjYdYzadjYzYd = Matmul(Yd,adjYdYzadjYzYd)
   YdadjYdYzadjYzAYd = Matmul(Yd,adjYdYzadjYzAYd)
   YdadjYdAYdadjYdYd = Matmul(Yd,adjYdAYdadjYdYd)
   YdadjYdAYdadjYdYs = Matmul(Yd,adjYdAYdadjYdYs)
   YdadjYdAYdadjYdYz = Matmul(Yd,adjYdAYdadjYdYz)
   YdadjYdAYsCYsYd = Matmul(Yd,adjYdAYsCYsYd)
   YdadjYdAYzadjYzYd = Matmul(Yd,adjYdAYzadjYzYd)
   YdadjYumu2YuadjYd = Matmul(Yd,adjYumu2YuadjYd)
   YdadjYuYumq2adjYd = Matmul(Yd,adjYuYumq2adjYd)
   YdadjYuYuadjYdmd2 = Matmul(Yd,adjYuYuadjYdmd2)
   YdadjYuYuadjYdYd = Matmul(Yd,adjYuYuadjYdYd)
   YdadjYuYuadjYdYs = Matmul(Yd,adjYuYuadjYdYs)
   YdadjYuYuadjYdYz = Matmul(Yd,adjYuYuadjYdYz)
   YdadjYuYuadjYdAYd = Matmul(Yd,adjYuYuadjYdAYd)
   YdadjYuYuadjYdAYs = Matmul(Yd,adjYuYuadjYdAYs)
   YdadjYuYuadjYdAYz = Matmul(Yd,adjYuYuadjYdAYz)
   YdadjYuYuadjYuYu = Matmul(Yd,adjYuYuadjYuYu)
   YdadjYuYuadjYuAYu = Matmul(Yd,adjYuYuadjYuAYu)
   YdadjYuAYuadjYdYd = Matmul(Yd,adjYuAYuadjYdYd)
   YdadjYuAYuadjYdYs = Matmul(Yd,adjYuAYuadjYdYs)
   YdadjYuAYuadjYdYz = Matmul(Yd,adjYuAYuadjYdYz)
   YdadjYuAYuadjYuYu = Matmul(Yd,adjYuAYuadjYuYu)
   Yeml2adjYeYeadjYe = Matmul(Ye,ml2adjYeYeadjYe)
   Yeml2adjYzYzadjYe = Matmul(Ye,ml2adjYzYzadjYe)
   Yeml2CYtYtadjYe = Matmul(Ye,ml2CYtYtadjYe)
   YeadjYeme2YeadjYe = Matmul(Ye,adjYeme2YeadjYe)
   YeadjYeYeml2adjYe = Matmul(Ye,adjYeYeml2adjYe)
   YeadjYeYeadjYeme2 = Matmul(Ye,adjYeYeadjYeme2)
   YeadjYeYeadjYeYe = Matmul(Ye,adjYeYeadjYeYe)
   YeadjYeYeadjYeAYe = Matmul(Ye,adjYeYeadjYeAYe)
   YeadjYeAYeadjYeYe = Matmul(Ye,adjYeAYeadjYeYe)
   YeadjYzmd2YzadjYe = Matmul(Ye,adjYzmd2YzadjYe)
   YeadjYzYdadjYdYz = Matmul(Ye,adjYzYdadjYdYz)
   YeadjYzYdadjYdAYz = Matmul(Ye,adjYzYdadjYdAYz)
   YeadjYzYsCYsYz = Matmul(Ye,adjYzYsCYsYz)
   YeadjYzYsCYsAYz = Matmul(Ye,adjYzYsCYsAYz)
   YeadjYzYzml2adjYe = Matmul(Ye,adjYzYzml2adjYe)
   YeadjYzYzadjYeme2 = Matmul(Ye,adjYzYzadjYeme2)
   YeadjYzYzadjYeYe = Matmul(Ye,adjYzYzadjYeYe)
   YeadjYzYzadjYeAYe = Matmul(Ye,adjYzYzadjYeAYe)
   YeadjYzYzadjYzYz = Matmul(Ye,adjYzYzadjYzYz)
   YeadjYzYzadjYzAYz = Matmul(Ye,adjYzYzadjYzAYz)
   YeadjYzAYdadjYdYz = Matmul(Ye,adjYzAYdadjYdYz)
   YeadjYzAYsCYsYz = Matmul(Ye,adjYzAYsCYsYz)
   YeadjYzAYzadjYeYe = Matmul(Ye,adjYzAYzadjYeYe)
   YeadjYzAYzadjYzYz = Matmul(Ye,adjYzAYzadjYzYz)
   YeCYtml2YtadjYe = Matmul(Ye,CYtml2YtadjYe)
   YeCYtYtml2adjYe = Matmul(Ye,CYtYtml2adjYe)
   YeCYtYtadjYeme2 = Matmul(Ye,CYtYtadjYeme2)
   YeCYtYtadjYeYe = Matmul(Ye,CYtYtadjYeYe)
   YeCYtYtadjYeAYe = Matmul(Ye,CYtYtadjYeAYe)
   YeCYtYtCYtYt = Matmul(Ye,CYtYtCYtYt)
   YeCYtYtCYtAYt = Matmul(Ye,CYtYtCYtAYt)
   YeCYtAYtadjYeYe = Matmul(Ye,CYtAYtadjYeYe)
   YeCYtAYtCYtYt = Matmul(Ye,CYtAYtCYtYt)
   YeCYtTYeCYeYt = Matmul(Ye,CYtTYeCYeYt)
   YeCYtTYeCYeAYt = Matmul(Ye,CYtTYeCYeAYt)
   YeCYtTYzCYzYt = Matmul(Ye,CYtTYzCYzYt)
   YeCYtTYzCYzAYt = Matmul(Ye,CYtTYzCYzAYt)
   YeCYtTAYeCYeYt = Matmul(Ye,CYtTAYeCYeYt)
   YeCYtTAYzCYzYt = Matmul(Ye,CYtTAYzCYzYt)
   Ysmd2CYdTYdCYs = Matmul(Ys,md2CYdTYdCYs)
   Ysmd2CYsYdadjYd = Matmul(Ys,md2CYsYdadjYd)
   Ysmd2CYsYsCYs = Matmul(Ys,md2CYsYsCYs)
   Ysmd2CYsYzadjYz = Matmul(Ys,md2CYsYzadjYz)
   Ysmd2CYzTYzCYs = Matmul(Ys,md2CYzTYzCYs)
   YsCYdmq2TYdCYs = Matmul(Ys,CYdmq2TYdCYs)
   YsCYdTYdmd2CYs = Matmul(Ys,CYdTYdmd2CYs)
   YsCYdTYdCYdTYd = Matmul(Ys,CYdTYdCYdTYd)
   YsCYdTYdCYdTAYd = Matmul(Ys,CYdTYdCYdTAYd)
   YsCYdTYdCYsmd2 = Matmul(Ys,CYdTYdCYsmd2)
   YsCYdTYdCYsYd = Matmul(Ys,CYdTYdCYsYd)
   YsCYdTYdCYsYs = Matmul(Ys,CYdTYdCYsYs)
   YsCYdTYdCYsYz = Matmul(Ys,CYdTYdCYsYz)
   YsCYdTYdCYsAYd = Matmul(Ys,CYdTYdCYsAYd)
   YsCYdTYdCYsAYs = Matmul(Ys,CYdTYdCYsAYs)
   YsCYdTYdCYsAYz = Matmul(Ys,CYdTYdCYsAYz)
   YsCYdTYuCYuTYd = Matmul(Ys,CYdTYuCYuTYd)
   YsCYdTYuCYuTAYd = Matmul(Ys,CYdTYuCYuTAYd)
   YsCYdTAYdCYdTYd = Matmul(Ys,CYdTAYdCYdTYd)
   YsCYdTAYdCYsYd = Matmul(Ys,CYdTAYdCYsYd)
   YsCYdTAYdCYsYs = Matmul(Ys,CYdTAYdCYsYs)
   YsCYdTAYdCYsYz = Matmul(Ys,CYdTAYdCYsYz)
   YsCYdTAYuCYuTYd = Matmul(Ys,CYdTAYuCYuTYd)
   YsCYsmd2YdadjYd = Matmul(Ys,CYsmd2YdadjYd)
   YsCYsmd2YsCYs = Matmul(Ys,CYsmd2YsCYs)
   YsCYsmd2YzadjYz = Matmul(Ys,CYsmd2YzadjYz)
   YsCYsYdmq2adjYd = Matmul(Ys,CYsYdmq2adjYd)
   YsCYsYdadjYdmd2 = Matmul(Ys,CYsYdadjYdmd2)
   YsCYsYdadjYdYs = Matmul(Ys,CYsYdadjYdYs)
   YsCYsYdadjYdAYs = Matmul(Ys,CYsYdadjYdAYs)
   YsCYsYsmd2CYs = Matmul(Ys,CYsYsmd2CYs)
   YsCYsYsCYsmd2 = Matmul(Ys,CYsYsCYsmd2)
   YsCYsYsCYsYd = Matmul(Ys,CYsYsCYsYd)
   YsCYsYsCYsYs = Matmul(Ys,CYsYsCYsYs)
   YsCYsYsCYsYz = Matmul(Ys,CYsYsCYsYz)
   YsCYsYsCYsAYd = Matmul(Ys,CYsYsCYsAYd)
   YsCYsYsCYsAYs = Matmul(Ys,CYsYsCYsAYs)
   YsCYsYsCYsAYz = Matmul(Ys,CYsYsCYsAYz)
   YsCYsYzml2adjYz = Matmul(Ys,CYsYzml2adjYz)
   YsCYsYzadjYzmd2 = Matmul(Ys,CYsYzadjYzmd2)
   YsCYsYzadjYzYs = Matmul(Ys,CYsYzadjYzYs)
   YsCYsYzadjYzAYs = Matmul(Ys,CYsYzadjYzAYs)
   YsCYsAYdadjYdYs = Matmul(Ys,CYsAYdadjYdYs)
   YsCYsAYsCYsYd = Matmul(Ys,CYsAYsCYsYd)
   YsCYsAYsCYsYs = Matmul(Ys,CYsAYsCYsYs)
   YsCYsAYsCYsYz = Matmul(Ys,CYsAYsCYsYz)
   YsCYsAYzadjYzYs = Matmul(Ys,CYsAYzadjYzYs)
   YsCYzml2TYzCYs = Matmul(Ys,CYzml2TYzCYs)
   YsCYzYtCYtTYz = Matmul(Ys,CYzYtCYtTYz)
   YsCYzYtCYtTAYz = Matmul(Ys,CYzYtCYtTAYz)
   YsCYzAYtCYtTYz = Matmul(Ys,CYzAYtCYtTYz)
   YsCYzTYeCYeTYz = Matmul(Ys,CYzTYeCYeTYz)
   YsCYzTYeCYeTAYz = Matmul(Ys,CYzTYeCYeTAYz)
   YsCYzTYzmd2CYs = Matmul(Ys,CYzTYzmd2CYs)
   YsCYzTYzCYsmd2 = Matmul(Ys,CYzTYzCYsmd2)
   YsCYzTYzCYsYd = Matmul(Ys,CYzTYzCYsYd)
   YsCYzTYzCYsYs = Matmul(Ys,CYzTYzCYsYs)
   YsCYzTYzCYsYz = Matmul(Ys,CYzTYzCYsYz)
   YsCYzTYzCYsAYd = Matmul(Ys,CYzTYzCYsAYd)
   YsCYzTYzCYsAYs = Matmul(Ys,CYzTYzCYsAYs)
   YsCYzTYzCYsAYz = Matmul(Ys,CYzTYzCYsAYz)
   YsCYzTYzCYzTYz = Matmul(Ys,CYzTYzCYzTYz)
   YsCYzTYzCYzTAYz = Matmul(Ys,CYzTYzCYzTAYz)
   YsCYzTAYeCYeTYz = Matmul(Ys,CYzTAYeCYeTYz)
   YsCYzTAYzCYsYd = Matmul(Ys,CYzTAYzCYsYd)
   YsCYzTAYzCYsYs = Matmul(Ys,CYzTAYzCYsYs)
   YsCYzTAYzCYsYz = Matmul(Ys,CYzTAYzCYsYz)
   YsCYzTAYzCYzTYz = Matmul(Ys,CYzTAYzCYzTYz)
   Ytml2adjYeYeCYt = Matmul(Yt,ml2adjYeYeCYt)
   Ytml2adjYzYzCYt = Matmul(Yt,ml2adjYzYzCYt)
   Ytml2CYtYtCYt = Matmul(Yt,ml2CYtYtCYt)
   YtadjYeme2YeCYt = Matmul(Yt,adjYeme2YeCYt)
   YtadjYeYeml2CYt = Matmul(Yt,adjYeYeml2CYt)
   YtadjYeYeadjYeYe = Matmul(Yt,adjYeYeadjYeYe)
   YtadjYeYeadjYeAYe = Matmul(Yt,adjYeYeadjYeAYe)
   YtadjYeYeCYtml2 = Matmul(Yt,adjYeYeCYtml2)
   YtadjYeYeCYtYt = Matmul(Yt,adjYeYeCYtYt)
   YtadjYeYeCYtAYt = Matmul(Yt,adjYeYeCYtAYt)
   YtadjYeAYeadjYeYe = Matmul(Yt,adjYeAYeadjYeYe)
   YtadjYeAYeCYtYt = Matmul(Yt,adjYeAYeCYtYt)
   YtadjYzmd2YzCYt = Matmul(Yt,adjYzmd2YzCYt)
   YtadjYzYdadjYdYz = Matmul(Yt,adjYzYdadjYdYz)
   YtadjYzYdadjYdAYz = Matmul(Yt,adjYzYdadjYdAYz)
   YtadjYzYsCYsYz = Matmul(Yt,adjYzYsCYsYz)
   YtadjYzYsCYsAYz = Matmul(Yt,adjYzYsCYsAYz)
   YtadjYzYzml2CYt = Matmul(Yt,adjYzYzml2CYt)
   YtadjYzYzadjYzYz = Matmul(Yt,adjYzYzadjYzYz)
   YtadjYzYzadjYzAYz = Matmul(Yt,adjYzYzadjYzAYz)
   YtadjYzYzCYtml2 = Matmul(Yt,adjYzYzCYtml2)
   YtadjYzYzCYtYt = Matmul(Yt,adjYzYzCYtYt)
   YtadjYzYzCYtAYt = Matmul(Yt,adjYzYzCYtAYt)
   YtadjYzAYdadjYdYz = Matmul(Yt,adjYzAYdadjYdYz)
   YtadjYzAYsCYsYz = Matmul(Yt,adjYzAYsCYsYz)
   YtadjYzAYzadjYzYz = Matmul(Yt,adjYzAYzadjYzYz)
   YtadjYzAYzCYtYt = Matmul(Yt,adjYzAYzCYtYt)
   YtCYtml2YtCYt = Matmul(Yt,CYtml2YtCYt)
   YtCYtYtml2CYt = Matmul(Yt,CYtYtml2CYt)
   YtCYtYtCYtml2 = Matmul(Yt,CYtYtCYtml2)
   YtCYtYtCYtYt = Matmul(Yt,CYtYtCYtYt)
   YtCYtYtCYtAYt = Matmul(Yt,CYtYtCYtAYt)
   YtCYtAYtCYtYt = Matmul(Yt,CYtAYtCYtYt)
   YtCYtTYeCYeYt = Matmul(Yt,CYtTYeCYeYt)
   YtCYtTYeCYeAYt = Matmul(Yt,CYtTYeCYeAYt)
   YtCYtTYzCYzYt = Matmul(Yt,CYtTYzCYzYt)
   YtCYtTYzCYzAYt = Matmul(Yt,CYtTYzCYzAYt)
   YtCYtTAYeCYeYt = Matmul(Yt,CYtTAYeCYeYt)
   YtCYtTAYzCYzYt = Matmul(Yt,CYtTAYzCYzYt)
   Yumq2adjYdYdadjYu = Matmul(Yu,mq2adjYdYdadjYu)
   Yumq2adjYuYuadjYu = Matmul(Yu,mq2adjYuYuadjYu)
   YuadjYdmd2YdadjYu = Matmul(Yu,adjYdmd2YdadjYu)
   YuadjYdYdmq2adjYu = Matmul(Yu,adjYdYdmq2adjYu)
   YuadjYdYdadjYdYd = Matmul(Yu,adjYdYdadjYdYd)
   YuadjYdYdadjYdAYd = Matmul(Yu,adjYdYdadjYdAYd)
   YuadjYdYdadjYumu2 = Matmul(Yu,adjYdYdadjYumu2)
   YuadjYdYdadjYuYu = Matmul(Yu,adjYdYdadjYuYu)
   YuadjYdYdadjYuAYu = Matmul(Yu,adjYdYdadjYuAYu)
   YuadjYdYsCYsYd = Matmul(Yu,adjYdYsCYsYd)
   YuadjYdYsCYsAYd = Matmul(Yu,adjYdYsCYsAYd)
   YuadjYdYzadjYzYd = Matmul(Yu,adjYdYzadjYzYd)
   YuadjYdYzadjYzAYd = Matmul(Yu,adjYdYzadjYzAYd)
   YuadjYdAYdadjYdYd = Matmul(Yu,adjYdAYdadjYdYd)
   YuadjYdAYdadjYuYu = Matmul(Yu,adjYdAYdadjYuYu)
   YuadjYdAYsCYsYd = Matmul(Yu,adjYdAYsCYsYd)
   YuadjYdAYzadjYzYd = Matmul(Yu,adjYdAYzadjYzYd)
   YuadjYumu2YuadjYu = Matmul(Yu,adjYumu2YuadjYu)
   YuadjYuYumq2adjYu = Matmul(Yu,adjYuYumq2adjYu)
   YuadjYuYuadjYumu2 = Matmul(Yu,adjYuYuadjYumu2)
   YuadjYuYuadjYuYu = Matmul(Yu,adjYuYuadjYuYu)
   YuadjYuYuadjYuAYu = Matmul(Yu,adjYuYuadjYuAYu)
   YuadjYuAYuadjYuYu = Matmul(Yu,adjYuAYuadjYuYu)
   Yzml2adjYeYeadjYz = Matmul(Yz,ml2adjYeYeadjYz)
   Yzml2adjYzYdadjYd = Matmul(Yz,ml2adjYzYdadjYd)
   Yzml2adjYzYsCYs = Matmul(Yz,ml2adjYzYsCYs)
   Yzml2adjYzYzadjYz = Matmul(Yz,ml2adjYzYzadjYz)
   Yzml2CYtYtadjYz = Matmul(Yz,ml2CYtYtadjYz)
   YzadjYeme2YeadjYz = Matmul(Yz,adjYeme2YeadjYz)
   YzadjYeYeml2adjYz = Matmul(Yz,adjYeYeml2adjYz)
   YzadjYeYeadjYeYe = Matmul(Yz,adjYeYeadjYeYe)
   YzadjYeYeadjYeAYe = Matmul(Yz,adjYeYeadjYeAYe)
   YzadjYeYeadjYzmd2 = Matmul(Yz,adjYeYeadjYzmd2)
   YzadjYeYeadjYzYd = Matmul(Yz,adjYeYeadjYzYd)
   YzadjYeYeadjYzYs = Matmul(Yz,adjYeYeadjYzYs)
   YzadjYeYeadjYzYz = Matmul(Yz,adjYeYeadjYzYz)
   YzadjYeYeadjYzAYd = Matmul(Yz,adjYeYeadjYzAYd)
   YzadjYeYeadjYzAYs = Matmul(Yz,adjYeYeadjYzAYs)
   YzadjYeYeadjYzAYz = Matmul(Yz,adjYeYeadjYzAYz)
   YzadjYeAYeadjYeYe = Matmul(Yz,adjYeAYeadjYeYe)
   YzadjYeAYeadjYzYd = Matmul(Yz,adjYeAYeadjYzYd)
   YzadjYeAYeadjYzYs = Matmul(Yz,adjYeAYeadjYzYs)
   YzadjYeAYeadjYzYz = Matmul(Yz,adjYeAYeadjYzYz)
   YzadjYzmd2YdadjYd = Matmul(Yz,adjYzmd2YdadjYd)
   YzadjYzmd2YsCYs = Matmul(Yz,adjYzmd2YsCYs)
   YzadjYzmd2YzadjYz = Matmul(Yz,adjYzmd2YzadjYz)
   YzadjYzYdmq2adjYd = Matmul(Yz,adjYzYdmq2adjYd)
   YzadjYzYdadjYdYz = Matmul(Yz,adjYzYdadjYdYz)
   YzadjYzYdadjYdAYz = Matmul(Yz,adjYzYdadjYdAYz)
   YzadjYzYsCYsYz = Matmul(Yz,adjYzYsCYsYz)
   YzadjYzYsCYsAYz = Matmul(Yz,adjYzYsCYsAYz)
   YzadjYzYzml2adjYz = Matmul(Yz,adjYzYzml2adjYz)
   YzadjYzYzadjYzmd2 = Matmul(Yz,adjYzYzadjYzmd2)
   YzadjYzYzadjYzYd = Matmul(Yz,adjYzYzadjYzYd)
   YzadjYzYzadjYzYs = Matmul(Yz,adjYzYzadjYzYs)
   YzadjYzYzadjYzYz = Matmul(Yz,adjYzYzadjYzYz)
   YzadjYzYzadjYzAYd = Matmul(Yz,adjYzYzadjYzAYd)
   YzadjYzYzadjYzAYs = Matmul(Yz,adjYzYzadjYzAYs)
   YzadjYzYzadjYzAYz = Matmul(Yz,adjYzYzadjYzAYz)
   YzadjYzAYdadjYdYz = Matmul(Yz,adjYzAYdadjYdYz)
   YzadjYzAYsCYsYz = Matmul(Yz,adjYzAYsCYsYz)
   YzadjYzAYzadjYzYd = Matmul(Yz,adjYzAYzadjYzYd)
   YzadjYzAYzadjYzYs = Matmul(Yz,adjYzAYzadjYzYs)
   YzadjYzAYzadjYzYz = Matmul(Yz,adjYzAYzadjYzYz)
   YzCYtml2YtadjYz = Matmul(Yz,CYtml2YtadjYz)
   YzCYtYtml2adjYz = Matmul(Yz,CYtYtml2adjYz)
   YzCYtYtadjYzmd2 = Matmul(Yz,CYtYtadjYzmd2)
   YzCYtYtadjYzYd = Matmul(Yz,CYtYtadjYzYd)
   YzCYtYtadjYzYs = Matmul(Yz,CYtYtadjYzYs)
   YzCYtYtadjYzYz = Matmul(Yz,CYtYtadjYzYz)
   YzCYtYtadjYzAYd = Matmul(Yz,CYtYtadjYzAYd)
   YzCYtYtadjYzAYs = Matmul(Yz,CYtYtadjYzAYs)
   YzCYtYtadjYzAYz = Matmul(Yz,CYtYtadjYzAYz)
   YzCYtYtCYtYt = Matmul(Yz,CYtYtCYtYt)
   YzCYtYtCYtAYt = Matmul(Yz,CYtYtCYtAYt)
   YzCYtAYtadjYzYd = Matmul(Yz,CYtAYtadjYzYd)
   YzCYtAYtadjYzYs = Matmul(Yz,CYtAYtadjYzYs)
   YzCYtAYtadjYzYz = Matmul(Yz,CYtAYtadjYzYz)
   YzCYtAYtCYtYt = Matmul(Yz,CYtAYtCYtYt)
   YzCYtTYeCYeYt = Matmul(Yz,CYtTYeCYeYt)
   YzCYtTYeCYeAYt = Matmul(Yz,CYtTYeCYeAYt)
   YzCYtTYzCYzYt = Matmul(Yz,CYtTYzCYzYt)
   YzCYtTYzCYzAYt = Matmul(Yz,CYtTYzCYzAYt)
   YzCYtTAYeCYeYt = Matmul(Yz,CYtTAYeCYeYt)
   YzCYtTAYzCYzYt = Matmul(Yz,CYtTAYzCYzYt)
   AYdadjYdYdadjYdYd = Matmul(AYd,adjYdYdadjYdYd)
   AYdadjYdYdadjYdYs = Matmul(AYd,adjYdYdadjYdYs)
   AYdadjYdYdadjYdYz = Matmul(AYd,adjYdYdadjYdYz)
   AYdadjYdYsCYsYd = Matmul(AYd,adjYdYsCYsYd)
   AYdadjYdYzadjYzYd = Matmul(AYd,adjYdYzadjYzYd)
   AYdadjYuYuadjYdYd = Matmul(AYd,adjYuYuadjYdYd)
   AYdadjYuYuadjYdYs = Matmul(AYd,adjYuYuadjYdYs)
   AYdadjYuYuadjYdYz = Matmul(AYd,adjYuYuadjYdYz)
   AYdadjYuYuadjYuYu = Matmul(AYd,adjYuYuadjYuYu)
   AYeadjYeYeadjYeYe = Matmul(AYe,adjYeYeadjYeYe)
   AYeadjYzYdadjYdYz = Matmul(AYe,adjYzYdadjYdYz)
   AYeadjYzYsCYsYz = Matmul(AYe,adjYzYsCYsYz)
   AYeadjYzYzadjYeYe = Matmul(AYe,adjYzYzadjYeYe)
   AYeadjYzYzadjYzYz = Matmul(AYe,adjYzYzadjYzYz)
   AYeCYtYtadjYeYe = Matmul(AYe,CYtYtadjYeYe)
   AYeCYtYtCYtYt = Matmul(AYe,CYtYtCYtYt)
   AYeCYtTYeCYeYt = Matmul(AYe,CYtTYeCYeYt)
   AYeCYtTYzCYzYt = Matmul(AYe,CYtTYzCYzYt)
   AYsCYdTYdCYdTYd = Matmul(AYs,CYdTYdCYdTYd)
   AYsCYdTYdCYsYd = Matmul(AYs,CYdTYdCYsYd)
   AYsCYdTYdCYsYs = Matmul(AYs,CYdTYdCYsYs)
   AYsCYdTYdCYsYz = Matmul(AYs,CYdTYdCYsYz)
   AYsCYdTYuCYuTYd = Matmul(AYs,CYdTYuCYuTYd)
   AYsCYsYdadjYdYs = Matmul(AYs,CYsYdadjYdYs)
   AYsCYsYsCYsYd = Matmul(AYs,CYsYsCYsYd)
   AYsCYsYsCYsYs = Matmul(AYs,CYsYsCYsYs)
   AYsCYsYsCYsYz = Matmul(AYs,CYsYsCYsYz)
   AYsCYsYzadjYzYs = Matmul(AYs,CYsYzadjYzYs)
   AYsCYzYtCYtTYz = Matmul(AYs,CYzYtCYtTYz)
   AYsCYzTYeCYeTYz = Matmul(AYs,CYzTYeCYeTYz)
   AYsCYzTYzCYsYd = Matmul(AYs,CYzTYzCYsYd)
   AYsCYzTYzCYsYs = Matmul(AYs,CYzTYzCYsYs)
   AYsCYzTYzCYsYz = Matmul(AYs,CYzTYzCYsYz)
   AYsCYzTYzCYzTYz = Matmul(AYs,CYzTYzCYzTYz)
   AYtadjYeYeadjYeYe = Matmul(AYt,adjYeYeadjYeYe)
   AYtadjYeYeCYtYt = Matmul(AYt,adjYeYeCYtYt)
   AYtadjYzYdadjYdYz = Matmul(AYt,adjYzYdadjYdYz)
   AYtadjYzYsCYsYz = Matmul(AYt,adjYzYsCYsYz)
   AYtadjYzYzadjYzYz = Matmul(AYt,adjYzYzadjYzYz)
   AYtadjYzYzCYtYt = Matmul(AYt,adjYzYzCYtYt)
   AYtCYtYtCYtYt = Matmul(AYt,CYtYtCYtYt)
   AYtCYtTYeCYeYt = Matmul(AYt,CYtTYeCYeYt)
   AYtCYtTYzCYzYt = Matmul(AYt,CYtTYzCYzYt)
   AYuadjYdYdadjYdYd = Matmul(AYu,adjYdYdadjYdYd)
   AYuadjYdYdadjYuYu = Matmul(AYu,adjYdYdadjYuYu)
   AYuadjYdYsCYsYd = Matmul(AYu,adjYdYsCYsYd)
   AYuadjYdYzadjYzYd = Matmul(AYu,adjYdYzadjYzYd)
   AYuadjYuYuadjYuYu = Matmul(AYu,adjYuYuadjYuYu)
   AYzadjYeYeadjYeYe = Matmul(AYz,adjYeYeadjYeYe)
   AYzadjYeYeadjYzYd = Matmul(AYz,adjYeYeadjYzYd)
   AYzadjYeYeadjYzYs = Matmul(AYz,adjYeYeadjYzYs)
   AYzadjYeYeadjYzYz = Matmul(AYz,adjYeYeadjYzYz)
   AYzadjYzYdadjYdYz = Matmul(AYz,adjYzYdadjYdYz)
   AYzadjYzYsCYsYz = Matmul(AYz,adjYzYsCYsYz)
   AYzadjYzYzadjYzYd = Matmul(AYz,adjYzYzadjYzYd)
   AYzadjYzYzadjYzYs = Matmul(AYz,adjYzYzadjYzYs)
   AYzadjYzYzadjYzYz = Matmul(AYz,adjYzYzadjYzYz)
   AYzCYtYtadjYzYd = Matmul(AYz,CYtYtadjYzYd)
   AYzCYtYtadjYzYs = Matmul(AYz,CYtYtadjYzYs)
   AYzCYtYtadjYzYz = Matmul(AYz,CYtYtadjYzYz)
   AYzCYtYtCYtYt = Matmul(AYz,CYtYtCYtYt)
   AYzCYtTYeCYeYt = Matmul(AYz,CYtTYeCYeYt)
   AYzCYtTYzCYzYt = Matmul(AYz,CYtTYzCYzYt)
   CYsmd2YsCYsYs = Matmul(adjYs,md2YsCYsYs)
   CYsYdmq2adjYdYs = Matmul(adjYs,Ydmq2adjYdYs)
   CYsYdadjYdmd2Ys = Matmul(adjYs,YdadjYdmd2Ys)
   CYsYdadjYdYsmd2 = Matmul(adjYs,YdadjYdYsmd2)
   CYsYsmd2CYsYs = Matmul(adjYs,Ysmd2CYsYs)
   CYsYsCYsmd2Ys = Matmul(adjYs,YsCYsmd2Ys)
   CYsYsCYsYsmd2 = Matmul(adjYs,YsCYsYsmd2)
   CYsYzml2adjYzYs = Matmul(adjYs,Yzml2adjYzYs)
   CYsYzadjYzmd2Ys = Matmul(adjYs,YzadjYzmd2Ys)
   CYsYzadjYzYsmd2 = Matmul(adjYs,YzadjYzYsmd2)
   CYtml2YtCYtYt = Matmul(adjYt,ml2YtCYtYt)
   CYtYtml2adjYeYe = Matmul(adjYt,Ytml2adjYeYe)
   CYtYtml2adjYzYz = Matmul(adjYt,Ytml2adjYzYz)
   CYtYtml2CYtYt = Matmul(adjYt,Ytml2CYtYt)
   CYtYtadjYeme2Ye = Matmul(adjYt,YtadjYeme2Ye)
   CYtYtadjYeYeml2 = Matmul(adjYt,YtadjYeYeml2)
   CYtYtadjYzmd2Yz = Matmul(adjYt,YtadjYzmd2Yz)
   CYtYtadjYzYzml2 = Matmul(adjYt,YtadjYzYzml2)
   CYtYtCYtml2Yt = Matmul(adjYt,YtCYtml2Yt)
   CYtYtCYtYtml2 = Matmul(adjYt,YtCYtYtml2)
   TYdmd2CYdTYdCYd = Matmul(Transpose(Yd),md2CYdTYdCYd)
   TYdmd2CYsYsCYd = Matmul(Transpose(Yd),md2CYsYsCYd)
   TYdmd2CYzTYzCYd = Matmul(Transpose(Yd),md2CYzTYzCYd)
   TYdCYdmq2TYdCYd = Matmul(Transpose(Yd),CYdmq2TYdCYd)
   TYdCYdTYdmd2CYd = Matmul(Transpose(Yd),CYdTYdmd2CYd)
   TYdCYdTYdCYdmq2 = Matmul(Transpose(Yd),CYdTYdCYdmq2)
   TYdCYsmd2YsCYd = Matmul(Transpose(Yd),CYsmd2YsCYd)
   TYdCYsYsmd2CYd = Matmul(Transpose(Yd),CYsYsmd2CYd)
   TYdCYsYsCYdmq2 = Matmul(Transpose(Yd),CYsYsCYdmq2)
   TYdCYzml2TYzCYd = Matmul(Transpose(Yd),CYzml2TYzCYd)
   TYdCYzTYzmd2CYd = Matmul(Transpose(Yd),CYzTYzmd2CYd)
   TYdCYzTYzCYdmq2 = Matmul(Transpose(Yd),CYzTYzCYdmq2)
   TYeme2CYeTYeCYe = Matmul(Transpose(Ye),me2CYeTYeCYe)
   TYeCYeml2TYeCYe = Matmul(Transpose(Ye),CYeml2TYeCYe)
   TYeCYeTYeme2CYe = Matmul(Transpose(Ye),CYeTYeme2CYe)
   TYeCYeTYeCYeml2 = Matmul(Transpose(Ye),CYeTYeCYeml2)
   TYeCYeTYeCYeYt = Matmul(Transpose(Ye),CYeTYeCYeYt)
   TYeCYeTYeCYeAYt = Matmul(Transpose(Ye),CYeTYeCYeAYt)
   TYeCYeTAYeCYeYt = Matmul(Transpose(Ye),CYeTAYeCYeYt)
   TYumu2CYuTYuCYu = Matmul(Transpose(Yu),mu2CYuTYuCYu)
   TYuCYumq2TYuCYu = Matmul(Transpose(Yu),CYumq2TYuCYu)
   TYuCYuTYumu2CYu = Matmul(Transpose(Yu),CYuTYumu2CYu)
   TYuCYuTYuCYumq2 = Matmul(Transpose(Yu),CYuTYuCYumq2)
   TYzmd2CYdTYdCYz = Matmul(Transpose(Yz),md2CYdTYdCYz)
   TYzmd2CYsYsCYz = Matmul(Transpose(Yz),md2CYsYsCYz)
   TYzmd2CYzTYzCYz = Matmul(Transpose(Yz),md2CYzTYzCYz)
   TYzCYdmq2TYdCYz = Matmul(Transpose(Yz),CYdmq2TYdCYz)
   TYzCYdTYdmd2CYz = Matmul(Transpose(Yz),CYdTYdmd2CYz)
   TYzCYdTYdCYzml2 = Matmul(Transpose(Yz),CYdTYdCYzml2)
   TYzCYdTYdCYzYt = Matmul(Transpose(Yz),CYdTYdCYzYt)
   TYzCYdTYdCYzAYt = Matmul(Transpose(Yz),CYdTYdCYzAYt)
   TYzCYdTAYdCYzYt = Matmul(Transpose(Yz),CYdTAYdCYzYt)
   TYzCYsmd2YsCYz = Matmul(Transpose(Yz),CYsmd2YsCYz)
   TYzCYsYsmd2CYz = Matmul(Transpose(Yz),CYsYsmd2CYz)
   TYzCYsYsCYzml2 = Matmul(Transpose(Yz),CYsYsCYzml2)
   TYzCYsYsCYzYt = Matmul(Transpose(Yz),CYsYsCYzYt)
   TYzCYsYsCYzAYt = Matmul(Transpose(Yz),CYsYsCYzAYt)
   TYzCYsAYsCYzYt = Matmul(Transpose(Yz),CYsAYsCYzYt)
   TYzCYzml2TYzCYz = Matmul(Transpose(Yz),CYzml2TYzCYz)
   TYzCYzTYzmd2CYz = Matmul(Transpose(Yz),CYzTYzmd2CYz)
   TYzCYzTYzCYzml2 = Matmul(Transpose(Yz),CYzTYzCYzml2)
   TYzCYzTYzCYzYt = Matmul(Transpose(Yz),CYzTYzCYzYt)
   TYzCYzTYzCYzAYt = Matmul(Transpose(Yz),CYzTYzCYzAYt)
   TYzCYzTAYzCYzYt = Matmul(Transpose(Yz),CYzTAYzCYzYt)
   TAYeCYeTYeCYeYt = Matmul(Transpose(AYe),CYeTYeCYeYt)
   TAYzCYdTYdCYzYt = Matmul(Transpose(AYz),CYdTYdCYzYt)
   TAYzCYsYsCYzYt = Matmul(Transpose(AYz),CYsYsCYzYt)
   TAYzCYzTYzCYzYt = Matmul(Transpose(AYz),CYzTYzCYzYt)
   TrYdadjAYd = cTrace(YdadjAYd)
   TrYeadjAYe = cTrace(YeadjAYe)
   TrYsCAYs = cTrace(YsCAYs)
   TrYtCAYt = cTrace(YtCAYt)
   TrYuadjAYu = cTrace(YuadjAYu)
   TrYzadjAYz = cTrace(YzadjAYz)
   Trmd2YsCYs = cTrace(md2YsCYs)
   Trml2YtCYt = cTrace(ml2YtCYt)
   TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd)
   TrYdadjYdYsCYs = cTrace(YdadjYdYsCYs)
   TrYdadjYdYzadjYz = cTrace(YdadjYdYzadjYz)
   TrYdadjYdAYdadjYd = cTrace(YdadjYdAYdadjYd)
   TrYdadjYdAYdadjAYd = cTrace(YdadjYdAYdadjAYd)
   TrYdadjYdAYsCYs = cTrace(YdadjYdAYsCYs)
   TrYdadjYdAYsCAYs = cTrace(YdadjYdAYsCAYs)
   TrYdadjYdAYzadjYz = cTrace(YdadjYdAYzadjYz)
   TrYdadjYdAYzadjAYz = cTrace(YdadjYdAYzadjAYz)
   TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd)
   TrYdadjYuAYuadjYd = cTrace(YdadjYuAYuadjYd)
   TrYdadjYuAYuadjAYd = cTrace(YdadjYuAYuadjAYd)
   TrYdadjAYdAYdadjYd = cTrace(YdadjAYdAYdadjYd)
   TrYdadjAYdAYsCYs = cTrace(YdadjAYdAYsCYs)
   TrYdadjAYdAYzadjYz = cTrace(YdadjAYdAYzadjYz)
   TrYdadjAYuAYuadjYd = cTrace(YdadjAYuAYuadjYd)
   TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe)
   TrYeadjYeAYeadjYe = cTrace(YeadjYeAYeadjYe)
   TrYeadjYeAYeadjAYe = cTrace(YeadjYeAYeadjAYe)
   TrYeadjYzYzadjYe = cTrace(YeadjYzYzadjYe)
   TrYeadjYzAYzadjYe = cTrace(YeadjYzAYzadjYe)
   TrYeadjYzAYzadjAYe = cTrace(YeadjYzAYzadjAYe)
   TrYeadjAYeAYeadjYe = cTrace(YeadjAYeAYeadjYe)
   TrYeadjAYzAYzadjYe = cTrace(YeadjAYzAYzadjYe)
   TrYeCYtYtadjYe = cTrace(YeCYtYtadjYe)
   TrYeCYtAYtadjYe = cTrace(YeCYtAYtadjYe)
   TrYeCYtAYtadjAYe = cTrace(YeCYtAYtadjAYe)
   TrYeCAYtAYtadjYe = cTrace(YeCAYtAYtadjYe)
   TrYsCYsYdadjYd = cTrace(YsCYsYdadjYd)
   TrYsCYsYsCYs = cTrace(YsCYsYsCYs)
   TrYsCYsYzadjYz = cTrace(YsCYsYzadjYz)
   TrYsCYsAYdadjYd = cTrace(YsCYsAYdadjYd)
   TrYsCYsAYdadjAYd = cTrace(YsCYsAYdadjAYd)
   TrYsCYsAYsCYs = cTrace(YsCYsAYsCYs)
   TrYsCYsAYsCAYs = cTrace(YsCYsAYsCAYs)
   TrYsCYsAYzadjYz = cTrace(YsCYsAYzadjYz)
   TrYsCYsAYzadjAYz = cTrace(YsCYsAYzadjAYz)
   TrYsCAYsAYdadjYd = cTrace(YsCAYsAYdadjYd)
   TrYsCAYsAYsCYs = cTrace(YsCAYsAYsCYs)
   TrYsCAYsAYzadjYz = cTrace(YsCAYsAYzadjYz)
   TrYtadjYeYeCYt = cTrace(YtadjYeYeCYt)
   TrYtadjYeAYeCYt = cTrace(YtadjYeAYeCYt)
   TrYtadjYeAYeCAYt = cTrace(YtadjYeAYeCAYt)
   TrYtadjYzYzCYt = cTrace(YtadjYzYzCYt)
   TrYtadjYzAYzCYt = cTrace(YtadjYzAYzCYt)
   TrYtadjYzAYzCAYt = cTrace(YtadjYzAYzCAYt)
   TrYtadjAYeAYeCYt = cTrace(YtadjAYeAYeCYt)
   TrYtadjAYzAYzCYt = cTrace(YtadjAYzAYzCYt)
   TrYtCYtYtCYt = cTrace(YtCYtYtCYt)
   TrYtCYtAYtCYt = cTrace(YtCYtAYtCYt)
   TrYtCYtAYtCAYt = cTrace(YtCYtAYtCAYt)
   TrYtCAYtAYtCYt = cTrace(YtCAYtAYtCYt)
   TrYuadjYdYdadjYu = cTrace(YuadjYdYdadjYu)
   TrYuadjYdAYdadjYu = cTrace(YuadjYdAYdadjYu)
   TrYuadjYdAYdadjAYu = cTrace(YuadjYdAYdadjAYu)
   TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu)
   TrYuadjYuAYuadjYu = cTrace(YuadjYuAYuadjYu)
   TrYuadjYuAYuadjAYu = cTrace(YuadjYuAYuadjAYu)
   TrYuadjAYdAYdadjYu = cTrace(YuadjAYdAYdadjYu)
   TrYuadjAYuAYuadjYu = cTrace(YuadjAYuAYuadjYu)
   TrYzadjYeYeadjYz = cTrace(YzadjYeYeadjYz)
   TrYzadjYeAYeadjYz = cTrace(YzadjYeAYeadjYz)
   TrYzadjYeAYeadjAYz = cTrace(YzadjYeAYeadjAYz)
   TrYzadjYzYdadjYd = cTrace(YzadjYzYdadjYd)
   TrYzadjYzYsCYs = cTrace(YzadjYzYsCYs)
   TrYzadjYzYzadjYz = cTrace(YzadjYzYzadjYz)
   TrYzadjYzAYdadjYd = cTrace(YzadjYzAYdadjYd)
   TrYzadjYzAYdadjAYd = cTrace(YzadjYzAYdadjAYd)
   TrYzadjYzAYsCYs = cTrace(YzadjYzAYsCYs)
   TrYzadjYzAYsCAYs = cTrace(YzadjYzAYsCAYs)
   TrYzadjYzAYzadjYz = cTrace(YzadjYzAYzadjYz)
   TrYzadjYzAYzadjAYz = cTrace(YzadjYzAYzadjAYz)
   TrYzadjAYeAYeadjYz = cTrace(YzadjAYeAYeadjYz)
   TrYzadjAYzAYdadjYd = cTrace(YzadjAYzAYdadjYd)
   TrYzadjAYzAYsCYs = cTrace(YzadjAYzAYsCYs)
   TrYzadjAYzAYzadjYz = cTrace(YzadjAYzAYzadjYz)
   TrYzCYtYtadjYz = cTrace(YzCYtYtadjYz)
   TrYzCYtAYtadjYz = cTrace(YzCYtAYtadjYz)
   TrYzCYtAYtadjAYz = cTrace(YzCYtAYtadjAYz)
   TrYzCAYtAYtadjYz = cTrace(YzCAYtAYtadjYz)
   TrCYsYdadjYdYs = cTrace(CYsYdadjYdYs)
   TrCYsYdadjYdAYs = cTrace(CYsYdadjYdAYs)
   TrCYsYdadjAYdAYs = cTrace(CYsYdadjAYdAYs)
   TrCYsYsCYsYs = cTrace(CYsYsCYsYs)
   TrCYsYsCYsAYs = cTrace(CYsYsCYsAYs)
   TrCYsYsCAYsAYs = cTrace(CYsYsCAYsAYs)
   TrCYsYzadjYzYs = cTrace(CYsYzadjYzYs)
   TrCYsYzadjYzAYs = cTrace(CYsYzadjYzAYs)
   TrCYsYzadjAYzAYs = cTrace(CYsYzadjAYzAYs)
   TrCYsAYdadjYdYs = cTrace(CYsAYdadjYdYs)
   TrCYsAYdadjAYdYs = cTrace(CYsAYdadjAYdYs)
   TrCYsAYsCYsYs = cTrace(CYsAYsCYsYs)
   TrCYsAYsCAYsYs = cTrace(CYsAYsCAYsYs)
   TrCYsAYzadjYzYs = cTrace(CYsAYzadjYzYs)
   TrCYsAYzadjAYzYs = cTrace(CYsAYzadjAYzYs)
   TrCYtYtadjYeYe = cTrace(CYtYtadjYeYe)
   TrCYtYtadjYeAYe = cTrace(CYtYtadjYeAYe)
   TrCYtYtadjYzYz = cTrace(CYtYtadjYzYz)
   TrCYtYtadjYzAYz = cTrace(CYtYtadjYzAYz)
   TrCYtYtadjAYeAYe = cTrace(CYtYtadjAYeAYe)
   TrCYtYtadjAYzAYz = cTrace(CYtYtadjAYzAYz)
   TrCYtYtCYtYt = cTrace(CYtYtCYtYt)
   TrCYtYtCYtAYt = cTrace(CYtYtCYtAYt)
   TrCYtYtCAYtAYt = cTrace(CYtYtCAYtAYt)
   TrCYtAYtadjYeYe = cTrace(CYtAYtadjYeYe)
   TrCYtAYtadjYzYz = cTrace(CYtAYtadjYzYz)
   TrCYtAYtadjAYeYe = cTrace(CYtAYtadjAYeYe)
   TrCYtAYtadjAYzYz = cTrace(CYtAYtadjAYzYz)
   TrCYtAYtCYtYt = cTrace(CYtAYtCYtYt)
   TrCYtAYtCAYtYt = cTrace(CYtAYtCAYtYt)
   Trmd2YdadjYdYsCYs = cTrace(md2YdadjYdYsCYs)
   Trmd2YsCYsYdadjYd = cTrace(md2YsCYsYdadjYd)
   Trmd2YsCYsYsCYs = cTrace(md2YsCYsYsCYs)
   Trmd2YsCYsYzadjYz = cTrace(md2YsCYsYzadjYz)
   Trmd2YzadjYzYsCYs = cTrace(md2YzadjYzYsCYs)
   Trmd2YzCYtYtadjYz = cTrace(md2YzCYtYtadjYz)
   Trmd2CYsYsCYsYs = cTrace(md2CYsYsCYsYs)
   Trme2YeCYtYtadjYe = cTrace(me2YeCYtYtadjYe)
   Trml2YtCYtYtCYt = cTrace(ml2YtCYtYtCYt)
   Trml2adjYzYsCYsYz = cTrace(ml2adjYzYsCYsYz)
   Trml2CYtYtCYtYt = cTrace(ml2CYtYtCYtYt)
   Trmq2adjYdYsCYsYd = cTrace(mq2adjYdYsCYsYd)
   TrYdmq2adjYdYdadjYd = cTrace(Ydmq2adjYdYdadjYd)
   TrYdmq2adjYdYsCYs = cTrace(Ydmq2adjYdYsCYs)
   TrYdmq2adjYdYzadjYz = cTrace(Ydmq2adjYdYzadjYz)
   TrYdmq2adjYuYuadjYd = cTrace(Ydmq2adjYuYuadjYd)
   TrYdadjYdmd2YdadjYd = cTrace(YdadjYdmd2YdadjYd)
   TrYdadjYdmd2YsCYs = cTrace(YdadjYdmd2YsCYs)
   TrYdadjYdmd2YzadjYz = cTrace(YdadjYdmd2YzadjYz)
   TrYdadjYdYdmq2adjYd = cTrace(YdadjYdYdmq2adjYd)
   TrYdadjYdYzml2adjYz = cTrace(YdadjYdYzml2adjYz)
   TrYdadjYumu2YuadjYd = cTrace(YdadjYumu2YuadjYd)
   TrYdadjYuYumq2adjYd = cTrace(YdadjYuYumq2adjYd)
   TrYeml2adjYeYeadjYe = cTrace(Yeml2adjYeYeadjYe)
   TrYeml2adjYzYzadjYe = cTrace(Yeml2adjYzYzadjYe)
   TrYeml2CYtYtadjYe = cTrace(Yeml2CYtYtadjYe)
   TrYeadjYeme2YeadjYe = cTrace(YeadjYeme2YeadjYe)
   TrYeadjYeYeml2adjYe = cTrace(YeadjYeYeml2adjYe)
   TrYeadjYzmd2YzadjYe = cTrace(YeadjYzmd2YzadjYe)
   TrYeadjYzYzml2adjYe = cTrace(YeadjYzYzml2adjYe)
   TrYeCYtml2YtadjYe = cTrace(YeCYtml2YtadjYe)
   TrYeCYtYtml2adjYe = cTrace(YeCYtYtml2adjYe)
   TrYsmd2CYsYdadjYd = cTrace(Ysmd2CYsYdadjYd)
   TrYsmd2CYsYsCYs = cTrace(Ysmd2CYsYsCYs)
   TrYsmd2CYsYzadjYz = cTrace(Ysmd2CYsYzadjYz)
   TrYsCYsmd2YdadjYd = cTrace(YsCYsmd2YdadjYd)
   TrYsCYsmd2YsCYs = cTrace(YsCYsmd2YsCYs)
   TrYsCYsmd2YzadjYz = cTrace(YsCYsmd2YzadjYz)
   TrYsCYsYdmq2adjYd = cTrace(YsCYsYdmq2adjYd)
   TrYsCYsYdadjYdmd2 = cTrace(YsCYsYdadjYdmd2)
   TrYsCYsYsmd2CYs = cTrace(YsCYsYsmd2CYs)
   TrYsCYsYsCYsmd2 = cTrace(YsCYsYsCYsmd2)
   TrYsCYsYzml2adjYz = cTrace(YsCYsYzml2adjYz)
   TrYsCYsYzadjYzmd2 = cTrace(YsCYsYzadjYzmd2)
   TrYtml2adjYeYeCYt = cTrace(Ytml2adjYeYeCYt)
   TrYtml2adjYzYzCYt = cTrace(Ytml2adjYzYzCYt)
   TrYtml2CYtYtCYt = cTrace(Ytml2CYtYtCYt)
   TrYtadjYeme2YeCYt = cTrace(YtadjYeme2YeCYt)
   TrYtadjYeYeml2CYt = cTrace(YtadjYeYeml2CYt)
   TrYtadjYeYeCYtml2 = cTrace(YtadjYeYeCYtml2)
   TrYtadjYzmd2YzCYt = cTrace(YtadjYzmd2YzCYt)
   TrYtadjYzYzml2CYt = cTrace(YtadjYzYzml2CYt)
   TrYtadjYzYzCYtml2 = cTrace(YtadjYzYzCYtml2)
   TrYtCYtml2YtCYt = cTrace(YtCYtml2YtCYt)
   TrYtCYtYtml2CYt = cTrace(YtCYtYtml2CYt)
   TrYtCYtYtCYtml2 = cTrace(YtCYtYtCYtml2)
   TrYumq2adjYdYdadjYu = cTrace(Yumq2adjYdYdadjYu)
   TrYumq2adjYuYuadjYu = cTrace(Yumq2adjYuYuadjYu)
   TrYuadjYdmd2YdadjYu = cTrace(YuadjYdmd2YdadjYu)
   TrYuadjYdYdmq2adjYu = cTrace(YuadjYdYdmq2adjYu)
   TrYuadjYumu2YuadjYu = cTrace(YuadjYumu2YuadjYu)
   TrYuadjYuYumq2adjYu = cTrace(YuadjYuYumq2adjYu)
   TrYzml2adjYeYeadjYz = cTrace(Yzml2adjYeYeadjYz)
   TrYzml2adjYzYdadjYd = cTrace(Yzml2adjYzYdadjYd)
   TrYzml2adjYzYsCYs = cTrace(Yzml2adjYzYsCYs)
   TrYzml2adjYzYzadjYz = cTrace(Yzml2adjYzYzadjYz)
   TrYzml2CYtYtadjYz = cTrace(Yzml2CYtYtadjYz)
   TrYzadjYeme2YeadjYz = cTrace(YzadjYeme2YeadjYz)
   TrYzadjYeYeml2adjYz = cTrace(YzadjYeYeml2adjYz)
   TrYzadjYzmd2YdadjYd = cTrace(YzadjYzmd2YdadjYd)
   TrYzadjYzmd2YsCYs = cTrace(YzadjYzmd2YsCYs)
   TrYzadjYzmd2YzadjYz = cTrace(YzadjYzmd2YzadjYz)
   TrYzadjYzYdmq2adjYd = cTrace(YzadjYzYdmq2adjYd)
   TrYzadjYzYzml2adjYz = cTrace(YzadjYzYzml2adjYz)
   TrYzCYtml2YtadjYz = cTrace(YzCYtml2YtadjYz)
   TrYzCYtYtml2adjYz = cTrace(YzCYtYtml2adjYz)
   TrCYsmd2YsCYsYs = cTrace(CYsmd2YsCYsYs)
   TrCYsYdmq2adjYdYs = cTrace(CYsYdmq2adjYdYs)
   TrCYsYdadjYdmd2Ys = cTrace(CYsYdadjYdmd2Ys)
   TrCYsYdadjYdYsmd2 = cTrace(CYsYdadjYdYsmd2)
   TrCYsYsmd2CYsYs = cTrace(CYsYsmd2CYsYs)
   TrCYsYsCYsmd2Ys = cTrace(CYsYsCYsmd2Ys)
   TrCYsYsCYsYsmd2 = cTrace(CYsYsCYsYsmd2)
   TrCYsYzml2adjYzYs = cTrace(CYsYzml2adjYzYs)
   TrCYsYzadjYzmd2Ys = cTrace(CYsYzadjYzmd2Ys)
   TrCYsYzadjYzYsmd2 = cTrace(CYsYzadjYzYsmd2)
   TrCYtml2YtCYtYt = cTrace(CYtml2YtCYtYt)
   TrCYtYtml2adjYeYe = cTrace(CYtYtml2adjYeYe)
   TrCYtYtml2adjYzYz = cTrace(CYtYtml2adjYzYz)
   TrCYtYtml2CYtYt = cTrace(CYtYtml2CYtYt)
   TrCYtYtadjYeme2Ye = cTrace(CYtYtadjYeme2Ye)
   TrCYtYtadjYeYeml2 = cTrace(CYtYtadjYeYeml2)
   TrCYtYtadjYzmd2Yz = cTrace(CYtYtadjYzmd2Yz)
   TrCYtYtadjYzYzml2 = cTrace(CYtYtadjYzYzml2)
   TrCYtYtCYtml2Yt = cTrace(CYtYtCYtml2Yt)
   TrCYtYtCYtYtml2 = cTrace(CYtYtCYtYtml2)
  End If



  Tr1(1) = -3._dp*mHd2/5._dp + 3._dp*mHu2/5._dp - 12._dp*ms2/5._dp + 12._dp*msb2/5._dp +&
  &  9._dp*mt2/5._dp - 9._dp*mtb2/5._dp + 3._dp*mz2/5._dp - 3._dp*mzb2/5._dp +           &
  &  3._dp*Trmd2/5._dp + 3._dp*Trme2/5._dp - 3._dp*Trml2/5._dp + 3._dp*Trmq2/5._dp -       &
  &  6._dp*Trmu2/5._dp


  If (TwoLoopRGE) Then
  Tr2(1) = 3._dp*mHd2/10._dp + 3._dp*mHu2/10._dp + 12._dp*ms2/5._dp + 12._dp*msb2/5._dp +&
  &  12._dp*mt2/5._dp + 12._dp*mtb2/5._dp + mz2/10._dp + mzb2/10._dp + Trmd2        &
  & /5._dp + 3._dp*Trme2/5._dp + 3._dp*Trml2/10._dp + Trmq2           &
  & /10._dp + 4._dp*Trmu2/5._dp

  Tr3(1) = -9._dp*g12*mHd2/100._dp - 9._dp*g22*mHd2/20._dp + 9._dp*g12*mHu2/100._dp +&
  &  9._dp*g22*mHu2/20._dp - 24._dp*g12*ms2/25._dp - 12._dp*g32*ms2 + 24._dp*g12*msb2/25._dp +&
  &  12._dp*g32*msb2 + 36._dp*g12*mt2/25._dp + 24._dp*g22*mt2/5._dp - 36._dp*g12*mtb2/25._dp -&
  &  24._dp*g22*mtb2/5._dp + (g12*mz2)/100._dp + 9._dp*g22*mz2/20._dp +              &
  &  4._dp*g32*mz2/5._dp - (g12*mzb2)/100._dp - 9._dp*g22*mzb2/20._dp -              &
  &  4._dp*g32*mzb2/5._dp + 9._dp*ms2*TrCYsYs/20._dp - 3._dp*TrCYsYsmd2/5._dp -          &
  &  3._dp*mt2*TrCYtYt/10._dp + 9._dp*TrCYtYtml2/20._dp - 3._dp*Trmd2YdadjYd/5._dp -       &
  &  3._dp*Trmd2YsCYs/5._dp - 3._dp*Trmd2YzadjYz/5._dp - 3._dp*Trme2YeadjYe/5._dp +        &
  &  9._dp*Trml2YtCYt/20._dp + 6._dp*Trmu2YuadjYu/5._dp + 9._dp*mHd2*TrYdadjYd/10._dp -    &
  &  3._dp*TrYdmq2adjYd/10._dp + 3._dp*mHd2*TrYeadjYe/10._dp + 3._dp*TrYeml2adjYe/10._dp + &
  &  27._dp*ms2*TrYsCYs/20._dp - 9._dp*mt2*TrYtCYt/10._dp - 9._dp*mHu2*TrYuadjYu/10._dp -  &
  &  3._dp*TrYumq2adjYu/10._dp - 3._dp*mz2*TrYzadjYz/10._dp + 9._dp*TrYzml2adjYz/10._dp +  &
  &  9._dp*L1*mHd2*Conjg(L1)/10._dp - 6._dp*L1*mt2*Conjg(L1)/5._dp - 9._dp*L2*mHu2*Conjg(L2)&
  & /10._dp + 6._dp*L2*mtb2*Conjg(L2)/5._dp + (g12*Trmd2)/25._dp + 4._dp*g32*Trmd2&
  & /5._dp + 9._dp*g12*Trme2/25._dp - 9._dp*g12*Trml2/100._dp -          &
  &  9._dp*g22*Trml2/20._dp + (g12*Trmq2)/100._dp + 9._dp*g22*Trmq2&
  & /20._dp + 4._dp*g32*Trmq2/5._dp - 8._dp*g12*Trmu2/25._dp -           &
  &  8._dp*g32*Trmu2/5._dp

  Tr2(2) = mHd2/2._dp + mHu2/2._dp + 4._dp*mt2 + 4._dp*mtb2 + 3._dp*mz2/2._dp +         &
  &  3._dp*mzb2/2._dp + Trml2/2._dp + 3._dp*Trmq2/2._dp

  Tr2(3) = 15._dp*ms2/2._dp + 15._dp*msb2/2._dp + mz2 + mzb2 + Trmd2             &
  & /2._dp + Trmq2 + Trmu2/2._dp

  Tr3(3) = -3._dp*ms2*TrCYsYs/4._dp - 3._dp*ms2*TrYsCYs/4._dp

  End If


  !--------------------
  ! g1
  !--------------------

  betag11 = 13.6_dp

  If (TwoLoopRGE) Then
   betag12 = ( 1502._dp*g12/75._dp + 34.8_dp * g22 + 184._dp * g32/3._dp  &
           & - 4.8_dp * TrCYsYs - 5.4_dp * (TrCYtYt + AbsL1sq+AbsL2sq)    &
           & - 2.8_dp * (TrYdadjYd+TrYzadjYz)  - 3.6_dp * TrYeadjYe       &
           & - 5.2_dp * TrYuadjYu )

   Dg1 = oo16pi2*g13*( betag11 + oo16pi2 * betag12 )

  Else
   Dg1 = oo16pi2* betag11*g13
  End If

  !--------------------
  ! g2
  !--------------------

  betag21 = 8._dp

  If (TwoLoopRGE) Then
   betag22 = ( 11.6_dp * g12 + 94._dp * g22 + 40._dp * g32                     &
           & - 6._dp * (TrYdadjYd + TrYuadjYu + TrYzadjYz) - 2._dp * TrYeadjYe &
           & - 7._dp*(TrCYtYt + AbsL1sq + AbsL2sq) )

   Dg2 = oo16pi2*g23*( betag21 + oo16pi2 * betag22 )

  Else
   Dg2 = oo16pi2*g23* betag21
  End If

  !--------------------
  ! g3
  !--------------------

  betag31 = 4._dp

  If (TwoLoopRGE) Then
   betag32 = ( 23._dp*g12/3._dp + 15._dp*g22 + 400._dp*g32/3._dp       &
           & - 9._dp*TrCYsYs - 4._dp*(TrYdadjYd+TrYuadjYu+TrYzadjYz) )

   Dg3 = oo16pi2*g33*( betag31 + oo16pi2 * betag32 )

  Else
   Dg3 = oo16pi2*g33* betag31
  End If

  !--------------------
  ! Yu
  !--------------------

  betaYu1 = ( -13._dp * g12 / 15._dp - 3._dp * g22 - 16._dp * g32 /3._dp &
          & + 3._dp * (TrYuadjYu + AbsL2sq ) ) * Yu                      &
          & + YuadjYdYd + 3._dp * YuadjYuYu

  If (TwoLoopRGE) Then
   betaYu2 = ( 5473._dp*g14/450._dp + (g12+8._dp*g32)*g22+ 28.5_dp*g24        &
      & + 136._dp*g12*g32/45._dp + 320._dp*g34/9._dp - 3._dp*TrYuadjYdYdadjYu &
      & + 0.8_dp*g12*TrYuadjYu + 16._dp*g32*TrYuadjYu- 9._dp*TrYuadjYuYuadjYu &
      & + (3.6_dp*g12+ 12._dp*g22 - 9._dp*TrYuadjYu - 12._dp*AbsL2sq)*AbsL2sq &
      & ) * Yu                                                                &
      & + (0.4_dp*g12 - 3._dp*(TrYdadjYd+AbsL1sq) - TrYeadjYe ) * YuadjYdYd   &
      & + (0.4_dp*g12 + 6._dp*g22  - 9._dp*(TrYuadjYu+AbsL2sq) ) * YuadjYuYu  &
      & - 2._dp * (YuadjYdYdadjYdYd + YuadjYdYdadjYuYu + YuadjYdYzadjYzYd )   &
      & - 4._dp * (YuadjYdYsCYsYd + YuadjYuYuadjYuYu)

   DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 )

  Else
   DYu = oo16pi2* betaYu1
  End If


  !--------------------
  ! Yd
  !--------------------

  betaYd1 = (-7._dp*g12/15._dp - 3._dp*g22 - 16._dp*g32/3._dp   &
          & + 3._dp*(TrYdadjYd+AbsL1sq) + TrYeadjYe )*Yd        &
          & + 3._dp*YdadjYdYd + YdadjYuYu + 4._dp*YsCYsYd + 2._dp*YzadjYzYd

  If (TwoLoopRGE) Then

   betaYd2 = ( 581._dp*g14/90._dp + (8._dp*g32+g12)*g22 + 28.5_dp*g24          &
          & + (8._dp*g12*g32+320._dp*g34)/9._dp - 9._dp*TrYdadjYdYdadjYd       &
          & + (16._dp*g32 - 0.4_dp*g12 - 9._dp*AbsL1sq ) *TrYdadjYd            &
          & - 3._dp *(TrYdadjYuYuadjYd+TrYeadjYeYeadjYe) + 1.2_dp*g12*TrYeadjYe&
          & - 3._dp*(TrYeadjYzYzadjYe+TrYeCYtYtadjYe)  - 12._dp*TrYsCYsYdadjYd &
          & + (3.6_dp*g12+ 12._dp*g22 - 3._dp*(TrCYtYt+TrYeadjYe)              &
          &    - 12._dp * AbsL1sq) * AbsL1sq - 6._dp*TrYzadjYzYdadjYd          &
          & ) * Yd                                                             &
          & + ( 0.8_dp*g12 + 6._dp*g22 - 3._dp*TrYeadjYe                       &
          &   - 9._dp*(TrYdadjYd+AbsL1sq)  ) * YdadjYdYd                       &
          & + ( 0.8_dp*g12  - 3._dp * (TrYuadjYu + AbsL2sq) ) * YdadjYuYu      &
          & + ( 32._dp*g12/15._dp + 80._dp*g32/3._dp -4._dp*TrCYsYs) * YsCYsYd &
          & + ( 0.4_dp*g12  + 6._dp*g22 - 2._dp*TrYzadjYz ) * YzadjYzYd        &
          & - 4._dp * (YdadjYdYdadjYdYd + YdadjYdYsCYsYd)                      &
          & - 2._dp * (YdadjYdYzadjYzYd + YdadjYuYuadjYdYd + YdadjYuYuadjYuYu  &
          &           +YzadjYeYeadjYzYd )   - 16._dp*YsCYsYsCYsYd              &
          & - 8._dp * (YsCYdTYdCYsYd + YsCYzTYzCYsYd)                          &
          & - 6._dp * (YzadjYzYzadjYzYd + YzCYtYtadjYzYd)

   DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 )

  Else
   DYd = oo16pi2* betaYd1
  End If


  !--------------------
  ! Ye
  !--------------------

  betaYe1 = ( 3._dp*(TrYdadjYd+AbsL1sq) + TrYeadjYe    &
          & - 1.8_dp*g12 - 3._dp*g22 ) * Ye            &
          & + 3._dp * (YeadjYeYe + YeadjYzYz + YeCYtYt)

  If (TwoLoopRGE) Then
   betaYe2 = ( 26.1_dp*g14 + 1.8_dp*g12*g22  + 28.5_dp*g24                   &
           & + (16._dp*g32 - 0.4_dp*g12) * TrYdadjYd  + 1.2_dp*g12*TrYeadjYe &
           & - 9._dp*TrYdadjYdYdadjYd - 12._dp*TrYsCYsYdadjYd                &
           & - 3._dp * (TrYdadjYuYuadjYd+TrYeadjYeYeadjYe+TrYeadjYzYzadjYe   &
           &           + TrYeCYtYtadjYe ) - 6._dp*TrYzadjYzYdadjYd           &
           & + (3.6_dp*g12 + 12._dp*g22  - 3._dp*(TrCYtYt+TrYeadjYe)         &
           &   - 9._dp*TrYdadjYd - 12._dp*AbsL1sq) * AbsL1sq                 &
           & ) * Ye                                                          &
           & + ( 6._dp*g22 - 3._dp*TrYeadjYe - 9._dp*(TrYdadjYd+AbsL1sq)     &
           &   ) * YeadjYeYe                                                 &
           & + ( 16._dp*g32 - 0.4_dp*g12  - 3._dp*TrYzadjYz ) * YeadjYzYz    &
           & + ( 3.6_dp*g12 + 12._dp*g22  - 3._dp*(TrYtCYt+AbsL1sq)          &
           &   ) * YeCYtYt                                                   &
           & - 4._dp*YeadjYeYeadjYeYe - 12._dp*YeadjYzYsCYsYz                &
           & - 6._dp * ( YeadjYzYdadjYdYz + YeadjYzYzadjYeYe                 &
           &           + YeadjYzYzadjYzYz + YeCYtYtadjYeYe )                 &
           & - 3._dp*YeCYtTYeCYeYt - 9._dp * (YeCYtTYzCYzYt+YeCYtYtCYtYt)

   DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 )

  Else
   DYe = oo16pi2* betaYe1
  End If

  !--------------------
  ! Yt
  !--------------------

  betaYt1 = (TrCYtYt + AbsL1sq - 1.8_dp*g12  - 7._dp*g22 ) * Yt         &
        & + TYeCYeYt + 3._dp * (TYzCYzYt+YtadjYzYz)+ YtadjYeYe + 6._dp*YtCYtYt

  If (TwoLoopRGE) Then
   betaYt2 = ( 26.1_dp*g14 + 11.4_dp*g12*g22 + 76.5_dp *g24  - TrCYtYtadjYeYe &
         &   - (0.6_dp*g12 + g22)* TrCYtYt - 3._dp * TrCYtYtadjYzYz           &
         &   - 6._dp*TrCYtYtCYtYt - TrYeCYtYtadjYe - 3._dp*TrYtadjYzYzCYt     &
         &   - ( 0.6_dp *g12 + g22 + 6._dp * (TrYdadjYd + AbsL1sq)            &
         &     - 2._dp * TrYeadjYe ) * AbsL1sq                                &
         &   ) * Yt                                                           &
         &   + ( 1.2_dp *g12  - TrYeadjYe - 3._dp*(TrYdadjYd+AbsL1sq)         &
         &     ) * (TYeCYeYt + YtadjYeYe)                                     &
         &   + (16._dp*g32 -0.4_dp*g12 -3._dp*TrYzadjYz )*(TYzCYzYt+YtadjYzYz)&
         &   + ( 7.2_dp*g12 + 24._dp*g22 - 6._dp*(TrCYtYt+AbsL1sq) ) *YtCYtYt &
         &   - 2._dp * (TYeCYeTYeCYeYt + YtadjYeYeadjYeYe)                    &
         &   - 6._dp * (TYzCYdTYdCYzYt + TYzCYzTYzCYzYt                       &
         &             +YtadjYzYdadjYdYz + YtadjYzYzadjYzYz)                  &
         &   - 12._dp *(TYzCYsYsCYzYt + YtadjYzYsCYsYz)                       &
         &   - 3._dp * (YtadjYeYeCYtYt + YtCYtTYeCYeYt)                       &
         &   - 9._dp * (YtadjYzYzCYtYt + YtCYtTYzCYzYt) - 18._dp *YtCYtYtCYtYt

   DYt = oo16pi2*( betaYt1 + oo16pi2 * betaYt2 )

  Else
   DYt = oo16pi2* betaYt1
  End If


  !--------------------
  ! Ys
  !--------------------

  betaYs1 = ( TrCYsYs - 0.8_dp *g12 - 12._dp*g32 ) * Ys              &
          & + 2._dp * ( YdadjYdYs + YsCYdTYd +YsCYzTYz + YzadjYzYs ) &
          & + 8._dp * YsCYsYs

  If (TwoLoopRGE) Then
   betaYs2 = ( 11.2_dp*g14 + 128._dp*g12*g32/15._dp + 320._dp*g34/3._dp       &
         &   - 4._dp * (g12/15._dp + g32/3._dp )* TrCYsYs                     &
         &   - 4._dp * (TrCYsYdadjYdYs+TrCYsYzadjYzYs) - 8._dp * TrCYsYsCYsYs &
         &   ) * Ys                                                           &
         &   + (0.4_dp*g12 + 6._dp*(g22-TrYdadjYd-AbsL1sq) - 2._dp*TrYeadjYe  &
         &     ) * (YdadjYdYs+YsCYdTYd)                                       &
         &   + (64._dp*g12/15._dp + 160._dp*g32/3._dp - 8._dp*TrCYsYs)*YsCYsYs&
         &   + (0.4_dp*g12 +6._dp*g22 -2._dp*TrYzadjYz) *(YsCYzTYz+YzadjYzYs) &
         &   - 2._dp * (YdadjYdYdadjYdYs + YdadjYuYuadjYdYs + YsCYdTYdCYdTYd  &
         &             +YsCYdTYuCYuTYd + YsCYzTYeCYeTYz + YzadjYeYeadjYzYs)   &
         &   - 8._dp * ( YsCYdTYdCYsYs + YsCYsYdadjYdYs                       &
         &             + YsCYsYzadjYzYs + YsCYzTYzCYsYs )                     &
         &   - 6._dp * ( YsCYzTYzCYzTYz + YsCYzYtCYtTYz                       &
         &             + YzadjYzYzadjYzYs + YzCYtYtadjYzYs)                   &
         &   - 32._dp*YsCYsYsCYsYs

   DYs = oo16pi2*( betaYs1 + oo16pi2 * betaYs2 )

  Else
   DYs = oo16pi2* betaYs1
  End If


  !--------------------
  ! Yz
  !--------------------
  betaYz1 = (TrYzadjYz -7._dp*g12/15._dp -3._dp*g22 -16._dp*g32/3._dp ) * Yz &
        & + 2._dp*YdadjYdYz + 4._dp*YsCYsYz + YzadjYeYe + 5._dp*YzadjYzYz    &
        & + 3._dp*YzCYtYt


  If (TwoLoopRGE) Then
   betaYz2 = ( 581._dp*g14/90._dp + g12*g22 + 28.5_dp*g24  + 320._dp*g34/9._dp &
         &   + 8._dp*g32*(g12/9._dp + g22) + 0.4_dp*g12*TrYzadjYz              &
         &   - TrYzadjYeYeadjYz- 2._dp*TrYdadjYdYzadjYz- 4._dp*TrYsCYsYzadjYz  &
         &   - 5._dp*TrYzadjYzYzadjYz - 3._dp*TrYzCYtYtadjYz                   &
         &   ) * Yz                                                            &
         &   + (0.4_dp*g12 + 6._dp*(g22 - TrYdadjYd -AbsL1sq) -2._dp*TrYeadjYe &
         &      ) * YdadjYdYz                                                  &
         &   + (32._dp*g12/15._dp + 80._dp*g32/3._dp - 4._dp*TrCYsYs) *YsCYsYz &
         &   + (1.2_dp*g12 - 3._dp*(TrYdadjYd+AbsL1sq) - TrYeadjYe) *YzadjYeYe &
         &   + (3.6_dp*g12 + 12._dp*g22 - 3._dp*(TrCYtYt+AbsL1sq) ) * YzCYtYt  &
         &   + (6._dp*g22 + 16._dp*g32  - 5._dp*TrYzadjYz) * YzadjYzYz         &
         &   - 2._dp * ( YdadjYdYdadjYdYz + YdadjYuYuadjYdYz                   &
         &             + YzadjYeYeadjYeYe + YzadjYeYeadjYzYz )                 &
         &   - 9._dp * ( YzCYtTYzCYzYt + YzCYtYtCYtYt) - 3._dp*YzCYtTYeCYeYt   &
         &   - 8._dp * ( YsCYdTYdCYsYz +YsCYzTYzCYsYz) - 16._dp*YsCYsYsCYsYz   &
         &   - 6._dp * ( YzadjYzYdadjYdYz + YzCYtYtadjYzYz)                    &
         &   - 12._dp * (YzadjYzYsCYsYz + YzadjYzYzadjYzYz)

   DYz = oo16pi2*( betaYz1 + oo16pi2 * betaYz2 )

  Else
   DYz = oo16pi2* betaYz1
  End If


  !--------------------
  ! L1
  !--------------------

  betaL11 = -1.8_dp*g12 - 7._dp * g22 + TrCYtYt + 6._dp * TrYdadjYd &
          & + 2._dp*TrYeadjYe + 7._dp*AbsL1sq

  If (TwoLoopRGE) Then
   betaL12 = 26.1_dp*g14 + 11.4_dp*g12*g22 + 76.5_dp*g24                  &
         &   - (0.6_dp*g12+g22 + 6._dp* AbsL1sq) *TrCYtYt                 &
         &   + (32._dp*g32 - 0.8_dp*g12 - 24._dp* AbsL1sq) * TrYdadjYd    &
         &   + (2.4_dp*g12 - 8._dp* AbsL1sq) * TrYeadjYe                  &
         &   + (6.6_dp*g12 + 23._dp*g22 - 30._dp* AbsL1sq) * AbsL1sq      &
         &   - 8._dp*TrYeCYtYtadjYe  - 12._dp * TrYzadjYzYdadjYd          &
         &   - 6._dp * (TrCYtYtCYtYt + TrYdadjYuYuadjYd + TrCYtYtadjYzYz  &
         &             +TrYeadjYeYeadjYe + TrYeadjYzYzadjYe)              &
         &   - 18._dp*TrYdadjYdYdadjYd - 24._dp * TrYsCYsYdadjYd

   DL1 = oo16pi2 * L1 * ( betaL11 + oo16pi2 * betaL12 )

  Else
   DL1 = oo16pi2 * L1 * betaL11
  End If


  !--------------------
  ! L2
  !--------------------

  betaL21 = -1.8_dp*g12 - 7._dp * g22 + 6._dp * TrYuadjYu + 7._dp*AbsL2sq

  If (TwoLoopRGE) Then
   betaL22 = 26.1_dp*g14 + 11.4_dp*g12*g22 + 76.5_dp*g24             &
         & + (1.6_dp*g12 + 32._dp*g32 - 24._dp*AbsL2sq ) * TrYuadjYu &
         & + (6.6_dp*g12 + 23._dp*g22 - 30._dp*AbsL2sq ) * AbsL2sq   &
         & - 6._dp * TrYuadjYdYdadjYu - 18._dp*TrYuadjYuYuadjYu

   DL2 = oo16pi2 * L2 * ( betaL21 + oo16pi2 * betaL22 )

  Else
   DL2 = oo16pi2 * L2 * betaL21
  End If


  !--------------------
  ! MTM
  !--------------------

  betaMTM1 = -2.4_dp*g12 - 8._dp*g22 + TrCYtYt + AbsL1sq + AbsL2sq

  If (TwoLoopRGE) Then
   betaMTM2 = 35.52_dp*g14 + 19.2_dp*g12*g22 + 96._dp*g24           &
          & - (0.6_dp*g12 + g22) * TrCYtYt - 2._dp * TrCYtYtadjYeYe &
          & - 6._dp*TrCYtYtadjYzYz - 6._dp*TrCYtYtCYtYt             &
          & - (0.6_dp*g12 + g22 + 6._dp*(TrYdadjYd+AbsL1sq)         &
          &   + 2._dp*TrYeadjYe ) * AbsL1sq                         &
          & - (0.6_dp*g12 + g22 + 6._dp*(TrYuadjYu+AbsL2sq) ) * AbsL2sq

   DMTM = oo16pi2 * MTM * ( betaMTM1 + oo16pi2 * betaMTM2 )

  Else
   DMTM = oo16pi2 * MTM * betaMTM1
  End If

  !--------------------
  ! mue
  !--------------------

  betamue1 = -0.6_dp*g12 + TrYeadjYe &
         & + 3._dp * ( TrYdadjYd + TrYuadjYu + AbsL1sq + AbsL2sq - g22)

  If (TwoLoopRGE) Then
   betamue2 = 8.34_dp*g14 + 1.8_dp*g12*g22 + 28.5_dp*g24 &
          & + ( 16._dp*g32 - 0.4_dp*g12 - 9._dp*AbsL1sq) * TrYdadjYd &
          & + ( 16._dp*g32 + 0.8_dp*g12 - 9._dp*AbsL2sq) * TrYuadjYu &
          & + ( 1.2_dp*g12 - 3._dp*AbsL1sq) * TrYeadjYe &
          & + ( 3.6_dp*g12 + 12._dp*g22 - 3._dp*TrCYtYt - 12._dp*AbsL1sq) * AbsL1sq &
          & + ( 3.6_dp*g12 + 12._dp*g22 - 12._dp*AbsL2sq) * AbsL2sq &
          & - 9._dp * ( TrYdadjYdYdadjYd + TrYuadjYuYuadjYu ) &
          & - 3._dp * ( TrYdadjYuYuadjYd + TrYeadjYeYeadjYe + TrYuadjYdYdadjYu  &
          &           + TrYeadjYzYzadjYe + TrYeCYtYtadjYe )   &
          & - 12._dp*TrYsCYsYdadjYd - 6._dp*TrYzadjYzYdadjYd

   Dmue = oo16pi2*mue*( betamue1 + oo16pi2 * betamue2 )

  Else
   Dmue = oo16pi2*mue* betamue1
  End If


  !--------------------
  ! MZM
  !--------------------

  betaMZM1 = -g12/15._dp - 3._dp*g22 - 16._dp*g32/3._dp +  TrYzadjYz

  If (TwoLoopRGE) Then
   betaMZM2 = 409._dp*g14/450._dp + (g12*g22)/5._dp + 28.5_dp*g24 &
        &  + 16._dp*g32*(g12/45._dp + g22) + 320._dp*g34/9._dp    &
        &  - 2._dp*TrYdadjYdYzadjYz - 4._dp*TrYsCYsYzadjYz - TrYzadjYeYeadjYz  &
        &  + 0.4_dp*g12*TrYzadjYz - 5._dp*TrYzadjYzYzadjYz - 3._dp*TrYzCYtYtadjYz

   DMZM = oo16pi2*MZM*( betaMZM1 + oo16pi2 * betaMZM2 )

  Else
   DMZM = oo16pi2*MZM* betaMZM1
  End If


  !--------------------
  ! MSM
  !--------------------

  betaMSM1 = -16._dp*g12/15._dp - 40._dp*g32/3._dp + 0.75_dp*TrCYsYs &
         & + 0.25_dp * TrYsCYs

  If (TwoLoopRGE) Then
   betaMSM2 = 3392._dp*g14/225._dp +128._dp*g12*g32/9._dp +1280._dp*g34/9._dp &
          & - 4._dp * (g12/15._dp + g32/3._dp) * TrCYsYs                      &
          & - 8._dp*TrCYsYsCYsYs - 4._dp*(TrCYsYzadjYzYs + TrCYsYdadjYdYs)

   DMSM = oo16pi2*MSM*( betaMSM1 + oo16pi2 * betaMSM2 )

  Else
   DMSM = oo16pi2*MSM* betaMSM1
  End If


  !--------------------
  ! AYu
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betaAYu1(i1,i2) = AYuadjYdYd(i1,i2) + 5._dp*AYuadjYuYu(i1,i2) + 26._dp*g12*MassB*Yu(i1,i2)&
  & /15._dp + 32._dp*g32*MassG*Yu(i1,i2)/3._dp + 6._dp*g22*MassWB*Yu(i1,i2)            &
  &  + 6._dp*TrAYuadjYu*Yu(i1,i2) + 6._dp*AL2*Conjg(L2)*Yu(i1,i2) + 2._dp*YuadjYdAYd(i1,i2)&
  &  + 4._dp*YuadjYuAYu(i1,i2) - 13._dp*g12*AYu(i1,i2)/15._dp - 3._dp*g22*AYu(i1,i2)   &
  &  - 16._dp*g32*AYu(i1,i2)/3._dp + 3._dp*TrYuadjYu*AYu(i1,i2) + 3._dp*AbsL2sq     &
  & *AYu(i1,i2)
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betaAYu2(i1,i2) = 2._dp*g12*AYuadjYdYd(i1,i2)/5._dp - 3._dp*TrYdadjYd*AYuadjYdYd(i1,i2)&
  &  - TrYeadjYe*AYuadjYdYd(i1,i2) - 2._dp*AYuadjYdYdadjYdYd(i1,i2) - 4._dp*AYuadjYdYdadjYuYu(i1,i2)&
  &  - 4._dp*AYuadjYdYsCYsYd(i1,i2) - 2._dp*AYuadjYdYzadjYzYd(i1,i2) + 12._dp*g22*AYuadjYuYu(i1,i2)&
  &  - 15._dp*TrYuadjYu*AYuadjYuYu(i1,i2) - 6._dp*AYuadjYuYuadjYuYu(i1,i2) -               &
  &  3._dp*L1*AYuadjYdYd(i1,i2)*Conjg(L1) - 15._dp*L2*AYuadjYuYu(i1,i2)*Conjg(L2)          &
  &  - 10946._dp*g14*MassB*Yu(i1,i2)/225._dp - 2._dp*g12*g22*MassB*Yu(i1,i2)         &
  &  - 272._dp*g12*g32*MassB*Yu(i1,i2)/45._dp - 272._dp*g12*g32*MassG*Yu(i1,i2)    &
  & /45._dp - 16._dp*g22*g32*MassG*Yu(i1,i2) - 1280._dp*g34*MassG*Yu(i1,i2)          &
  & /9._dp - 2._dp*g12*g22*MassWB*Yu(i1,i2) - 114._dp*g24*MassWB*Yu(i1,i2)           &
  &  - 16._dp*g22*g32*MassWB*Yu(i1,i2) + 8._dp*g12*TrAYuadjYu*Yu(i1,i2)              &
  & /5._dp + 32._dp*g32*TrAYuadjYu*Yu(i1,i2) - 6._dp*TrYdadjYuAYuadjYd*Yu(i1,i2)         &
  &  - 6._dp*TrYuadjYdAYdadjYu*Yu(i1,i2) - 8._dp*g12*MassB*TrYuadjYu*Yu(i1,i2)           &
  & /5._dp - 32._dp*g32*MassG*TrYuadjYu*Yu(i1,i2) - 36._dp*TrYuadjYuAYuadjYu*Yu(i1,i2)   &
  &  - 36._dp*g12*L2*MassB*Conjg(L2)*Yu(i1,i2)/5._dp - 24._dp*g22*L2*MassWB*Conjg(L2)  &
  & *Yu(i1,i2) - 18._dp*L2*TrAYuadjYu*Conjg(L2)*Yu(i1,i2) + 36._dp*g12*AL2*Conjg(L2)     &
  & *Yu(i1,i2)/5._dp + 24._dp*g22*AL2*Conjg(L2)*Yu(i1,i2) - 18._dp*TrYuadjYu*AL2*Conjg(L2)&
  & *Yu(i1,i2) - 48._dp*L2*AL2*Conjg(L2)**2*Yu(i1,i2) + 4._dp*g12*YuadjYdAYd(i1,i2)      &
  & /5._dp - 6._dp*TrYdadjYd*YuadjYdAYd(i1,i2) - 2._dp*TrYeadjYe*YuadjYdAYd(i1,i2)         &
  &  - 6._dp*AbsL1sq*YuadjYdAYd(i1,i2) - 4._dp*YuadjYdAYdadjYdYd(i1,i2)               &
  &  - 4._dp*YuadjYdAYdadjYuYu(i1,i2) - 8._dp*YuadjYdAYsCYsYd(i1,i2) - 4._dp*YuadjYdAYzadjYzYd(i1,i2)&
  &  - 4._dp*g12*MassB*YuadjYdYd(i1,i2)/5._dp - 6._dp*TrAYdadjYd*YuadjYdYd(i1,i2)        &
  &  - 2._dp*TrAYeadjYe*YuadjYdYd(i1,i2) - 6._dp*AL1*Conjg(L1)*YuadjYdYd(i1,i2)            &
  &  - 4._dp*YuadjYdYdadjYdAYd(i1,i2) - 2._dp*YuadjYdYdadjYuAYu(i1,i2) - 8._dp*YuadjYdYsCYsAYd(i1,i2)&
  &  - 4._dp*YuadjYdYzadjYzAYd(i1,i2) + 6._dp*g12*YuadjYuAYu(i1,i2)/5._dp +              &
  &  6._dp*g22*YuadjYuAYu(i1,i2) - 12._dp*TrYuadjYu*YuadjYuAYu(i1,i2) - 12._dp*AbsL2sq&
  & *YuadjYuAYu(i1,i2) - 8._dp*YuadjYuAYuadjYuYu(i1,i2) - 4._dp*g12*MassB*YuadjYuYu(i1,i2)&
  & /5._dp - 12._dp*g22*MassWB*YuadjYuYu(i1,i2) - 18._dp*TrAYuadjYu*YuadjYuYu(i1,i2)     &
  &  - 18._dp*AL2*Conjg(L2)*YuadjYuYu(i1,i2) - 6._dp*YuadjYuYuadjYuAYu(i1,i2)              &
  &  + 5473._dp*g14*AYu(i1,i2)/450._dp + g12*g22*AYu(i1,i2) + 57._dp*g24*AYu(i1,i2)&
  & /2._dp + 136._dp*g12*g32*AYu(i1,i2)/45._dp + 8._dp*g22*g32*AYu(i1,i2)          &
  &  + 320._dp*g34*AYu(i1,i2)/9._dp - 3._dp*TrYdadjYuYuadjYd*AYu(i1,i2) + 4._dp*g12*TrYuadjYu*AYu(i1,i2)&
  & /5._dp + 16._dp*g32*TrYuadjYu*AYu(i1,i2) - 9._dp*TrYuadjYuYuadjYu*AYu(i1,i2)         &
  &  + 18._dp*g12*AbsL2sq*AYu(i1,i2)/5._dp + 12._dp*g22*AbsL2sq              &
  & *AYu(i1,i2) - 9._dp*L2*TrYuadjYu*Conjg(L2)*AYu(i1,i2) - 12._dp*L2**2*Conjg(L2)         &
  & **2*AYu(i1,i2)
   End Do
  End Do


  DAYu = oo16pi2*( betaAYu1 + oo16pi2 * betaAYu2 )


  Else
  DAYu = oo16pi2* betaAYu1
  End If


  !--------------------
  ! AYd
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betaAYd1(i1,i2) = 5._dp*AYdadjYdYd(i1,i2) + AYdadjYuYu(i1,i2) + 8._dp*AYsCYsYd(i1,i2) &
  &  + 4._dp*AYzadjYzYd(i1,i2) + 14._dp*g12*MassB*Yd(i1,i2)/15._dp + 32._dp*g32*MassG*Yd(i1,i2)&
  & /3._dp + 6._dp*g22*MassWB*Yd(i1,i2) + 6._dp*TrAYdadjYd*Yd(i1,i2) + 2._dp*TrAYeadjYe*Yd(i1,i2)&
  &  + 6._dp*AL1*Conjg(L1)*Yd(i1,i2) + 4._dp*YdadjYdAYd(i1,i2) + 2._dp*YdadjYuAYu(i1,i2)   &
  &  + 4._dp*YsCYsAYd(i1,i2) + 2._dp*YzadjYzAYd(i1,i2) - 7._dp*g12*AYd(i1,i2)            &
  & /15._dp - 3._dp*g22*AYd(i1,i2) - 16._dp*g32*AYd(i1,i2)/3._dp + 3._dp*TrYdadjYd*AYd(i1,i2)&
  &  + TrYeadjYe*AYd(i1,i2) + 3._dp*AbsL1sq*AYd(i1,i2)
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betaAYd2(i1,i2) = 6._dp*g12*AYdadjYdYd(i1,i2)/5._dp + 12._dp*g22*AYdadjYdYd(i1,i2)&
  &  - 15._dp*TrYdadjYd*AYdadjYdYd(i1,i2) - 5._dp*TrYeadjYe*AYdadjYdYd(i1,i2)              &
  &  - 6._dp*AYdadjYdYdadjYdYd(i1,i2) - 4._dp*AYdadjYdYsCYsYd(i1,i2) - 2._dp*AYdadjYdYzadjYzYd(i1,i2)&
  &  + 4._dp*g12*AYdadjYuYu(i1,i2)/5._dp - 3._dp*TrYuadjYu*AYdadjYuYu(i1,i2)             &
  &  - 4._dp*AYdadjYuYuadjYdYd(i1,i2) - 2._dp*AYdadjYuYuadjYuYu(i1,i2) - 16._dp*AYsCYdTYdCYsYd(i1,i2)&
  &  + 64._dp*g12*AYsCYsYd(i1,i2)/15._dp + 160._dp*g32*AYsCYsYd(i1,i2)/3._dp -         &
  &  6._dp*TrCYsYs*AYsCYsYd(i1,i2) - 2._dp*TrYsCYs*AYsCYsYd(i1,i2) - 32._dp*AYsCYsYsCYsYd(i1,i2)&
  &  - 16._dp*AYsCYzTYzCYsYd(i1,i2) - 4._dp*AYzadjYeYeadjYzYd(i1,i2) + 4._dp*g12*AYzadjYzYd(i1,i2)&
  & /5._dp + 12._dp*g22*AYzadjYzYd(i1,i2) - 4._dp*TrYzadjYz*AYzadjYzYd(i1,i2)            &
  &  - 12._dp*AYzadjYzYzadjYzYd(i1,i2) - 12._dp*AYzCYtYtadjYzYd(i1,i2) - 15._dp*L1*AYdadjYdYd(i1,i2)&
  & *Conjg(L1) - 3._dp*L2*AYdadjYuYu(i1,i2)*Conjg(L2) - 1162._dp*g14*MassB*Yd(i1,i2)     &
  & /45._dp - 2._dp*g12*g22*MassB*Yd(i1,i2) - 16._dp*g12*g32*MassB*Yd(i1,i2)       &
  & /9._dp - 16._dp*g12*g32*MassG*Yd(i1,i2)/9._dp - 16._dp*g22*g32*MassG*Yd(i1,i2) &
  &  - 1280._dp*g34*MassG*Yd(i1,i2)/9._dp - 2._dp*g12*g22*MassWB*Yd(i1,i2)           &
  &  - 114._dp*g24*MassWB*Yd(i1,i2) - 16._dp*g22*g32*MassWB*Yd(i1,i2) -              &
  &  4._dp*g12*TrAYdadjYd*Yd(i1,i2)/5._dp + 32._dp*g32*TrAYdadjYd*Yd(i1,i2)            &
  &  + 12._dp*g12*TrAYeadjYe*Yd(i1,i2)/5._dp - 12._dp*TrCYsAYdadjYdYs*Yd(i1,i2)          &
  &  - 15._dp*TrCYsYdadjYdAYs*Yd(i1,i2) - 3._dp*TrCYtAYtadjYeYe*Yd(i1,i2) - 4._dp*TrCYtYtadjYeAYe*Yd(i1,i2)&
  &  + 4._dp*g12*MassB*TrYdadjYd*Yd(i1,i2)/5._dp - 32._dp*g32*MassG*TrYdadjYd*Yd(i1,i2)&
  &  - 36._dp*TrYdadjYdAYdadjYd*Yd(i1,i2) - 9._dp*TrYdadjYdAYsCYs*Yd(i1,i2) -              &
  &  12._dp*TrYdadjYdAYzadjYz*Yd(i1,i2) - 6._dp*TrYdadjYuAYuadjYd*Yd(i1,i2) -              &
  &  12._dp*g12*MassB*TrYeadjYe*Yd(i1,i2)/5._dp - 12._dp*TrYeadjYeAYeadjYe*Yd(i1,i2)     &
  &  - 6._dp*TrYeadjYzAYzadjYe*Yd(i1,i2) - 3._dp*TrYeCYtAYtadjYe*Yd(i1,i2) -               &
  &  12._dp*TrYsCYsAYdadjYd*Yd(i1,i2) - 2._dp*TrYtadjYeAYeCYt*Yd(i1,i2) - 6._dp*TrYuadjYdAYdadjYu*Yd(i1,i2)&
  &  - 6._dp*TrYzadjYeAYeadjYz*Yd(i1,i2) - 12._dp*TrYzadjYzAYdadjYd*Yd(i1,i2)              &
  &  - 36._dp*g12*L1*MassB*Conjg(L1)*Yd(i1,i2)/5._dp - 24._dp*g22*L1*MassWB*Conjg(L1)  &
  & *Yd(i1,i2) - 18._dp*L1*TrAYdadjYd*Conjg(L1)*Yd(i1,i2) - 6._dp*L1*TrAYeadjYe*Conjg(L1)  &
  & *Yd(i1,i2) - 6._dp*L1*TrCYtAYt*Conjg(L1)*Yd(i1,i2) + 36._dp*g12*AL1*Conjg(L1)        &
  & *Yd(i1,i2)/5._dp + 24._dp*g22*AL1*Conjg(L1)*Yd(i1,i2) - 9._dp*TrCYtYt*AL1*Conjg(L1)  &
  & *Yd(i1,i2)/2._dp - 18._dp*TrYdadjYd*AL1*Conjg(L1)*Yd(i1,i2) - 6._dp*TrYeadjYe*AL1*Conjg(L1)&
  & *Yd(i1,i2) - 3._dp*TrYtCYt*AL1*Conjg(L1)*Yd(i1,i2)/2._dp - 48._dp*L1*AL1*Conjg(L1)     &
  & **2*Yd(i1,i2) + 6._dp*g12*YdadjYdAYd(i1,i2)/5._dp + 6._dp*g22*YdadjYdAYd(i1,i2)    &
  &  - 12._dp*TrYdadjYd*YdadjYdAYd(i1,i2) - 4._dp*TrYeadjYe*YdadjYdAYd(i1,i2)              &
  &  - 12._dp*AbsL1sq*YdadjYdAYd(i1,i2) - 8._dp*YdadjYdAYdadjYdYd(i1,i2)              &
  &  - 8._dp*YdadjYdAYsCYsYd(i1,i2) - 4._dp*YdadjYdAYzadjYzYd(i1,i2) - 8._dp*g12*MassB*YdadjYdYd(i1,i2)&
  & /5._dp - 12._dp*g22*MassWB*YdadjYdYd(i1,i2) - 18._dp*TrAYdadjYd*YdadjYdYd(i1,i2)     &
  &  - 6._dp*TrAYeadjYe*YdadjYdYd(i1,i2) - 18._dp*AL1*Conjg(L1)*YdadjYdYd(i1,i2)           &
  &  - 6._dp*YdadjYdYdadjYdAYd(i1,i2) - 8._dp*YdadjYdYsCYsAYd(i1,i2) - 4._dp*YdadjYdYzadjYzAYd(i1,i2)&
  &  + 8._dp*g12*YdadjYuAYu(i1,i2)/5._dp - 6._dp*TrYuadjYu*YdadjYuAYu(i1,i2)             &
  &  - 6._dp*AbsL2sq*YdadjYuAYu(i1,i2) - 4._dp*YdadjYuAYuadjYdYd(i1,i2)               &
  &  - 4._dp*YdadjYuAYuadjYuYu(i1,i2) - 8._dp*g12*MassB*YdadjYuYu(i1,i2)/5._dp -         &
  &  6._dp*TrAYuadjYu*YdadjYuYu(i1,i2) - 6._dp*AL2*Conjg(L2)*YdadjYuYu(i1,i2)              &
  &  - 2._dp*YdadjYuYuadjYdAYd(i1,i2) - 4._dp*YdadjYuYuadjYuAYu(i1,i2) - 16._dp*YsCYdTAYdCYsYd(i1,i2)&
  &  - 8._dp*YsCYdTYdCYsAYd(i1,i2) + 32._dp*g12*YsCYsAYd(i1,i2)/15._dp + 80._dp*g32*YsCYsAYd(i1,i2)&
  & /3._dp - 3._dp*TrCYsYs*YsCYsAYd(i1,i2) - TrYsCYs*YsCYsAYd(i1,i2) - 32._dp*YsCYsAYsCYsYd(i1,i2)&
  &  - 64._dp*g12*MassB*YsCYsYd(i1,i2)/15._dp - 160._dp*g32*MassG*YsCYsYd(i1,i2)       &
  & /3._dp - 8._dp*TrCYsAYs*YsCYsYd(i1,i2) - 16._dp*YsCYsYsCYsAYd(i1,i2) - 16._dp*YsCYzTAYzCYsYd(i1,i2)&
  &  - 8._dp*YsCYzTYzCYsAYd(i1,i2) - 4._dp*YzadjYeAYeadjYzYd(i1,i2) - 2._dp*YzadjYeYeadjYzAYd(i1,i2)&
  &  + 2._dp*g12*YzadjYzAYd(i1,i2)/5._dp + 6._dp*g22*YzadjYzAYd(i1,i2) -               &
  &  2._dp*TrYzadjYz*YzadjYzAYd(i1,i2) - 12._dp*YzadjYzAYzadjYzYd(i1,i2) - 4._dp*g12*MassB*YzadjYzYd(i1,i2)&
  & /5._dp - 12._dp*g22*MassWB*YzadjYzYd(i1,i2) - 4._dp*TrAYzadjYz*YzadjYzYd(i1,i2)      &
  &  - 6._dp*YzadjYzYzadjYzAYd(i1,i2) - 12._dp*YzCYtAYtadjYzYd(i1,i2) - 6._dp*YzCYtYtadjYzAYd(i1,i2)&
  &  + 581._dp*g14*AYd(i1,i2)/90._dp + g12*g22*AYd(i1,i2) + 57._dp*g24*AYd(i1,i2)  &
  & /2._dp + 8._dp*g12*g32*AYd(i1,i2)/9._dp + 8._dp*g22*g32*AYd(i1,i2)             &
  &  + 320._dp*g34*AYd(i1,i2)/9._dp - 6._dp*TrCYsYdadjYdYs*AYd(i1,i2) - 3._dp*TrCYtYtadjYeYe*AYd(i1,i2)&
  & /2._dp - 2._dp*g12*TrYdadjYd*AYd(i1,i2)/5._dp + 16._dp*g32*TrYdadjYd*AYd(i1,i2)    &
  &  - 9._dp*TrYdadjYdYdadjYd*AYd(i1,i2) - 9._dp*TrYdadjYdYsCYs*AYd(i1,i2)/2._dp -         &
  &  6._dp*TrYdadjYdYzadjYz*AYd(i1,i2) + 6._dp*g12*TrYeadjYe*AYd(i1,i2)/5._dp -          &
  &  3._dp*TrYeadjYeYeadjYe*AYd(i1,i2) - (TrYeCYtYtadjYe*AYd(i1,i2))/2._dp -               &
  &  3._dp*TrYsCYsYdadjYd*AYd(i1,i2)/2._dp - TrYtadjYeYeCYt*AYd(i1,i2) - 3._dp*TrYuadjYdYdadjYu*AYd(i1,i2)&
  &  - 3._dp*TrYzadjYeYeadjYz*AYd(i1,i2) + 18._dp*g12*AbsL1sq*AYd(i1,i2)            &
  & /5._dp + 12._dp*g22*AbsL1sq*AYd(i1,i2) - 9._dp*L1*TrCYtYt*Conjg(L1)             &
  & *AYd(i1,i2)/4._dp - 9._dp*L1*TrYdadjYd*Conjg(L1)*AYd(i1,i2) - 3._dp*L1*TrYeadjYe*Conjg(L1)&
  & *AYd(i1,i2) - 3._dp*L1*TrYtCYt*Conjg(L1)*AYd(i1,i2)/4._dp - 12._dp*L1**2*Conjg(L1)     &
  & **2*AYd(i1,i2)
   End Do
  End Do


  DAYd = oo16pi2*( betaAYd1 + oo16pi2 * betaAYd2 )


  Else
  DAYd = oo16pi2* betaAYd1
  End If


  !--------------------
  ! AYe
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betaAYe1(i1,i2) = 5._dp*AYeadjYeYe(i1,i2) + 3._dp*AYeadjYzYz(i1,i2) + 3._dp*AYeCYtYt(i1,i2)&
  &  + 18._dp*g12*MassB*Ye(i1,i2)/5._dp + 6._dp*g22*MassWB*Ye(i1,i2) + 6._dp*TrAYdadjYd*Ye(i1,i2)&
  &  + 2._dp*TrAYeadjYe*Ye(i1,i2) + 6._dp*AL1*Conjg(L1)*Ye(i1,i2) + 4._dp*YeadjYeAYe(i1,i2)&
  &  + 6._dp*YeadjYzAYz(i1,i2) + 6._dp*YeCYtAYt(i1,i2) - 9._dp*g12*AYe(i1,i2)            &
  & /5._dp - 3._dp*g22*AYe(i1,i2) + 3._dp*TrYdadjYd*AYe(i1,i2) + TrYeadjYe*AYe(i1,i2)    &
  &  + 3._dp*AbsL1sq*AYe(i1,i2)
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betaAYe2(i1,i2) = -6._dp*g12*AYeadjYeYe(i1,i2)/5._dp + 12._dp*g22*AYeadjYeYe(i1,i2)&
  &  - 15._dp*TrYdadjYd*AYeadjYeYe(i1,i2) - 5._dp*TrYeadjYe*AYeadjYeYe(i1,i2)              &
  &  - 6._dp*AYeadjYeYeadjYeYe(i1,i2) - 6._dp*AYeadjYzYdadjYdYz(i1,i2) - 12._dp*AYeadjYzYsCYsYz(i1,i2)&
  &  - 2._dp*g12*AYeadjYzYz(i1,i2)/5._dp + 16._dp*g32*AYeadjYzYz(i1,i2) -              &
  &  3._dp*TrYzadjYz*AYeadjYzYz(i1,i2) - 12._dp*AYeadjYzYzadjYeYe(i1,i2) - 6._dp*AYeadjYzYzadjYzYz(i1,i2)&
  &  - 3._dp*AYeCYtTYeCYeYt(i1,i2) - 9._dp*AYeCYtTYzCYzYt(i1,i2) + 18._dp*g12*AYeCYtYt(i1,i2)&
  & /5._dp + 12._dp*g22*AYeCYtYt(i1,i2) - 9._dp*TrCYtYt*AYeCYtYt(i1,i2)/4._dp -          &
  &  3._dp*TrYtCYt*AYeCYtYt(i1,i2)/4._dp - 12._dp*AYeCYtYtadjYeYe(i1,i2) - 9._dp*AYeCYtYtCYtYt(i1,i2)&
  &  - 15._dp*L1*AYeadjYeYe(i1,i2)*Conjg(L1) - 3._dp*L1*AYeCYtYt(i1,i2)*Conjg(L1)          &
  &  - 522._dp*g14*MassB*Ye(i1,i2)/5._dp - 18._dp*g12*g22*MassB*Ye(i1,i2)            &
  & /5._dp - 18._dp*g12*g22*MassWB*Ye(i1,i2)/5._dp - 114._dp*g24*MassWB*Ye(i1,i2)    &
  &  - 4._dp*g12*TrAYdadjYd*Ye(i1,i2)/5._dp + 32._dp*g32*TrAYdadjYd*Ye(i1,i2)          &
  &  + 12._dp*g12*TrAYeadjYe*Ye(i1,i2)/5._dp - 12._dp*TrCYsAYdadjYdYs*Ye(i1,i2)          &
  &  - 15._dp*TrCYsYdadjYdAYs*Ye(i1,i2) - 3._dp*TrCYtAYtadjYeYe*Ye(i1,i2) - 4._dp*TrCYtYtadjYeAYe*Ye(i1,i2)&
  &  + 4._dp*g12*MassB*TrYdadjYd*Ye(i1,i2)/5._dp - 32._dp*g32*MassG*TrYdadjYd*Ye(i1,i2)&
  &  - 36._dp*TrYdadjYdAYdadjYd*Ye(i1,i2) - 9._dp*TrYdadjYdAYsCYs*Ye(i1,i2) -              &
  &  12._dp*TrYdadjYdAYzadjYz*Ye(i1,i2) - 6._dp*TrYdadjYuAYuadjYd*Ye(i1,i2) -              &
  &  12._dp*g12*MassB*TrYeadjYe*Ye(i1,i2)/5._dp - 12._dp*TrYeadjYeAYeadjYe*Ye(i1,i2)     &
  &  - 6._dp*TrYeadjYzAYzadjYe*Ye(i1,i2) - 3._dp*TrYeCYtAYtadjYe*Ye(i1,i2) -               &
  &  12._dp*TrYsCYsAYdadjYd*Ye(i1,i2) - 2._dp*TrYtadjYeAYeCYt*Ye(i1,i2) - 6._dp*TrYuadjYdAYdadjYu*Ye(i1,i2)&
  &  - 6._dp*TrYzadjYeAYeadjYz*Ye(i1,i2) - 12._dp*TrYzadjYzAYdadjYd*Ye(i1,i2)              &
  &  - 36._dp*g12*L1*MassB*Conjg(L1)*Ye(i1,i2)/5._dp - 24._dp*g22*L1*MassWB*Conjg(L1)  &
  & *Ye(i1,i2) - 18._dp*L1*TrAYdadjYd*Conjg(L1)*Ye(i1,i2) - 6._dp*L1*TrAYeadjYe*Conjg(L1)  &
  & *Ye(i1,i2) - 6._dp*L1*TrCYtAYt*Conjg(L1)*Ye(i1,i2) + 36._dp*g12*AL1*Conjg(L1)        &
  & *Ye(i1,i2)/5._dp + 24._dp*g22*AL1*Conjg(L1)*Ye(i1,i2) - 9._dp*TrCYtYt*AL1*Conjg(L1)  &
  & *Ye(i1,i2)/2._dp - 18._dp*TrYdadjYd*AL1*Conjg(L1)*Ye(i1,i2) - 6._dp*TrYeadjYe*AL1*Conjg(L1)&
  & *Ye(i1,i2) - 3._dp*TrYtCYt*AL1*Conjg(L1)*Ye(i1,i2)/2._dp - 48._dp*L1*AL1*Conjg(L1)     &
  & **2*Ye(i1,i2) + 6._dp*g12*YeadjYeAYe(i1,i2)/5._dp + 6._dp*g22*YeadjYeAYe(i1,i2)    &
  &  - 12._dp*TrYdadjYd*YeadjYeAYe(i1,i2) - 4._dp*TrYeadjYe*YeadjYeAYe(i1,i2)              &
  &  - 12._dp*AbsL1sq*YeadjYeAYe(i1,i2) - 8._dp*YeadjYeAYeadjYeYe(i1,i2)              &
  &  - 12._dp*g22*MassWB*YeadjYeYe(i1,i2) - 18._dp*TrAYdadjYd*YeadjYeYe(i1,i2)           &
  &  - 6._dp*TrAYeadjYe*YeadjYeYe(i1,i2) - 18._dp*AL1*Conjg(L1)*YeadjYeYe(i1,i2)           &
  &  - 6._dp*YeadjYeYeadjYeAYe(i1,i2) - 12._dp*YeadjYzAYdadjYdYz(i1,i2) - 24._dp*YeadjYzAYsCYsYz(i1,i2)&
  &  - 4._dp*g12*YeadjYzAYz(i1,i2)/5._dp + 32._dp*g32*YeadjYzAYz(i1,i2) -              &
  &  6._dp*TrYzadjYz*YeadjYzAYz(i1,i2) - 12._dp*YeadjYzAYzadjYeYe(i1,i2) - 12._dp*YeadjYzAYzadjYzYz(i1,i2)&
  &  - 12._dp*YeadjYzYdadjYdAYz(i1,i2) - 24._dp*YeadjYzYsCYsAYz(i1,i2) + 4._dp*g12*MassB*YeadjYzYz(i1,i2)&
  & /5._dp - 32._dp*g32*MassG*YeadjYzYz(i1,i2) - 6._dp*TrAYzadjYz*YeadjYzYz(i1,i2)       &
  &  - 6._dp*YeadjYzYzadjYeAYe(i1,i2) - 12._dp*YeadjYzYzadjYzAYz(i1,i2) + 36._dp*g12*YeCYtAYt(i1,i2)&
  & /5._dp + 24._dp*g22*YeCYtAYt(i1,i2) - 9._dp*TrCYtYt*YeCYtAYt(i1,i2)/2._dp -          &
  &  3._dp*TrYtCYt*YeCYtAYt(i1,i2)/2._dp - 6._dp*AbsL1sq*YeCYtAYt(i1,i2)              &
  &  - 12._dp*YeCYtAYtadjYeYe(i1,i2) - 18._dp*YeCYtAYtCYtYt(i1,i2) - 6._dp*YeCYtTAYeCYeYt(i1,i2)&
  &  - 18._dp*YeCYtTAYzCYzYt(i1,i2) - 6._dp*YeCYtTYeCYeAYt(i1,i2) - 18._dp*YeCYtTYzCYzAYt(i1,i2)&
  &  - 36._dp*g12*MassB*YeCYtYt(i1,i2)/5._dp - 24._dp*g22*MassWB*YeCYtYt(i1,i2)        &
  &  - 6._dp*TrCYtAYt*YeCYtYt(i1,i2) - 6._dp*AL1*Conjg(L1)*YeCYtYt(i1,i2) - 6._dp*YeCYtYtadjYeAYe(i1,i2)&
  &  - 18._dp*YeCYtYtCYtAYt(i1,i2) + 261._dp*g14*AYe(i1,i2)/10._dp + 9._dp*g12*g22*AYe(i1,i2)&
  & /5._dp + 57._dp*g24*AYe(i1,i2)/2._dp - 6._dp*TrCYsYdadjYdYs*AYe(i1,i2)               &
  &  - 3._dp*TrCYtYtadjYeYe*AYe(i1,i2)/2._dp - 2._dp*g12*TrYdadjYd*AYe(i1,i2)            &
  & /5._dp + 16._dp*g32*TrYdadjYd*AYe(i1,i2) - 9._dp*TrYdadjYdYdadjYd*AYe(i1,i2)         &
  &  - 9._dp*TrYdadjYdYsCYs*AYe(i1,i2)/2._dp - 6._dp*TrYdadjYdYzadjYz*AYe(i1,i2)           &
  &  + 6._dp*g12*TrYeadjYe*AYe(i1,i2)/5._dp - 3._dp*TrYeadjYeYeadjYe*AYe(i1,i2)          &
  &  - (TrYeCYtYtadjYe*AYe(i1,i2))/2._dp - 3._dp*TrYsCYsYdadjYd*AYe(i1,i2)/2._dp -         &
  &  TrYtadjYeYeCYt*AYe(i1,i2) - 3._dp*TrYuadjYdYdadjYu*AYe(i1,i2) - 3._dp*TrYzadjYeYeadjYz*AYe(i1,i2)&
  &  + 18._dp*g12*AbsL1sq*AYe(i1,i2)/5._dp + 12._dp*g22*AbsL1sq              &
  & *AYe(i1,i2) - 9._dp*L1*TrCYtYt*Conjg(L1)*AYe(i1,i2)/4._dp - 9._dp*L1*TrYdadjYd*Conjg(L1)&
  & *AYe(i1,i2) - 3._dp*L1*TrYeadjYe*Conjg(L1)*AYe(i1,i2) - 3._dp*L1*TrYtCYt*Conjg(L1)     &
  & *AYe(i1,i2)/4._dp - 12._dp*L1**2*Conjg(L1)**2*AYe(i1,i2)
   End Do
  End Do


  DAYe = oo16pi2*( betaAYe1 + oo16pi2 * betaAYe2 )


  Else
  DAYe = oo16pi2* betaAYe1
  End If


  !--------------------
  ! AYt
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betaAYt1(i1,i2) = (sqrt2*AYtadjYeYe(i1,i2) + 3._dp*sqrt2*AYtadjYzYz(i1,i2)&
  &  + 9._dp*sqrt2*AYtCYtYt(i1,i2) + 2._dp*sqrt2*TAYeCYeYt(i1,i2)              &
  &  + 6._dp*sqrt2*TAYzCYzYt(i1,i2) + sqrt2*TYeCYeAYt(i1,i2) + 3._dp*sqrt2&
  & *TYzCYzAYt(i1,i2) + 18._dp*sqrt2*g12*MassB*Yt(i1,i2)/5._dp + 14._dp*sqrt2&
  & *g22*MassWB*Yt(i1,i2) + 2._dp*sqrt2*TrCYtAYt*Yt(i1,i2) + 2._dp*sqrt2     &
  & *AL1*Conjg(L1)*Yt(i1,i2) + 2._dp*sqrt2*YtadjYeAYe(i1,i2) + 6._dp*sqrt2     &
  & *YtadjYzAYz(i1,i2) + 9._dp*sqrt2*YtCYtAYt(i1,i2) - 9._dp*sqrt2             &
  & *g12*AYt(i1,i2)/5._dp - 7._dp*sqrt2*g22*AYt(i1,i2) + 3._dp*TrCYtYt*AYt(i1,i2)&
  & /(2._dp*sqrt2) + (TrYtCYt*AYt(i1,i2))/(2._dp*sqrt2) + sqrt2          &
  & *AbsL1sq*AYt(i1,i2))/sqrt2
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betaAYt2(i1,i2) = (6._dp*sqrt2*g12*AYtadjYeYe(i1,i2)/5._dp - 3._dp*sqrt2&
  & *TrYdadjYd*AYtadjYeYe(i1,i2) - sqrt2*TrYeadjYe*AYtadjYeYe(i1,i2) -               &
  &  2._dp*sqrt2*AYtadjYeYeadjYeYe(i1,i2) - 6._dp*sqrt2*AYtadjYeYeCYtYt(i1,i2) &
  &  - 6._dp*sqrt2*AYtadjYzYdadjYdYz(i1,i2) - 12._dp*sqrt2*AYtadjYzYsCYsYz(i1,i2)&
  &  - 2._dp*sqrt2*g12*AYtadjYzYz(i1,i2)/5._dp + 16._dp*sqrt2*g32*AYtadjYzYz(i1,i2)&
  &  - 3._dp*sqrt2*TrYzadjYz*AYtadjYzYz(i1,i2) - 6._dp*sqrt2*AYtadjYzYzadjYzYz(i1,i2)&
  &  - 18._dp*sqrt2*AYtadjYzYzCYtYt(i1,i2) - 3._dp*sqrt2*AYtCYtTYeCYeYt(i1,i2) &
  &  - 9._dp*sqrt2*AYtCYtTYzCYzYt(i1,i2) + 54._dp*sqrt2*g12*AYtCYtYt(i1,i2)  &
  & /5._dp + 36._dp*sqrt2*g22*AYtCYtYt(i1,i2) - 7._dp*TrCYtYt*AYtCYtYt(i1,i2)      &
  & /(2._dp*sqrt2) - 5._dp*sqrt2*TrCYtYt*AYtCYtYt(i1,i2) - 9._dp*TrYtCYt*AYtCYtYt(i1,i2)&
  & /(2._dp*sqrt2) - 27._dp*sqrt2*AYtCYtYtCYtYt(i1,i2) - 3._dp*sqrt2     &
  & *L1*AYtadjYeYe(i1,i2)*Conjg(L1) - 9._dp*sqrt2*L1*AYtCYtYt(i1,i2)*Conjg(L1)       &
  &  - 4._dp*sqrt2*TAYeCYeTYeCYeYt(i1,i2) + 12._dp*sqrt2*g12*TAYeCYeYt(i1,i2)&
  & /5._dp - 6._dp*sqrt2*TrYdadjYd*TAYeCYeYt(i1,i2) - 2._dp*sqrt2              &
  & *TrYeadjYe*TAYeCYeYt(i1,i2) - 6._dp*sqrt2*AbsL1sq*TAYeCYeYt(i1,i2)          &
  &  - 12._dp*sqrt2*TAYzCYdTYdCYzYt(i1,i2) - 24._dp*sqrt2*TAYzCYsYsCYzYt(i1,i2)&
  &  - 12._dp*sqrt2*TAYzCYzTYzCYzYt(i1,i2) - 4._dp*sqrt2*g12*TAYzCYzYt(i1,i2)&
  & /5._dp + 32._dp*sqrt2*g32*TAYzCYzYt(i1,i2) - 6._dp*sqrt2*TrYzadjYz*TAYzCYzYt(i1,i2)&
  &  + 6._dp*sqrt2*g12*TYeCYeAYt(i1,i2)/5._dp - 3._dp*sqrt2*TrYdadjYd*TYeCYeAYt(i1,i2)&
  &  - sqrt2*TrYeadjYe*TYeCYeAYt(i1,i2) - 3._dp*sqrt2*AbsL1sq             &
  & *TYeCYeAYt(i1,i2) - 4._dp*sqrt2*TYeCYeTAYeCYeYt(i1,i2) - 2._dp*sqrt2       &
  & *TYeCYeTYeCYeAYt(i1,i2) - 12._dp*sqrt2*g12*MassB*TYeCYeYt(i1,i2)               &
  & /5._dp - 6._dp*sqrt2*TrAYdadjYd*TYeCYeYt(i1,i2) - 2._dp*sqrt2              &
  & *TrAYeadjYe*TYeCYeYt(i1,i2) - 6._dp*sqrt2*AL1*Conjg(L1)*TYeCYeYt(i1,i2)          &
  &  - 12._dp*sqrt2*TYzCYdTAYdCYzYt(i1,i2) - 6._dp*sqrt2*TYzCYdTYdCYzAYt(i1,i2)&
  &  - 24._dp*sqrt2*TYzCYsAYsCYzYt(i1,i2) - 12._dp*sqrt2*TYzCYsYsCYzAYt(i1,i2) &
  &  - 2._dp*sqrt2*g12*TYzCYzAYt(i1,i2)/5._dp + 16._dp*sqrt2*g32*TYzCYzAYt(i1,i2)&
  &  - 3._dp*sqrt2*TrYzadjYz*TYzCYzAYt(i1,i2) - 12._dp*sqrt2*TYzCYzTAYzCYzYt(i1,i2)&
  &  - 6._dp*sqrt2*TYzCYzTYzCYzAYt(i1,i2) + 4._dp*sqrt2*g12*MassB*TYzCYzYt(i1,i2)&
  & /5._dp - 32._dp*sqrt2*g32*MassG*TYzCYzYt(i1,i2) - 6._dp*sqrt2            &
  & *TrAYzadjYz*TYzCYzYt(i1,i2) - 522._dp*sqrt2*g14*MassB*Yt(i1,i2)/5._dp -        &
  &  114._dp*sqrt2*g12*g22*MassB*Yt(i1,i2)/5._dp - 114._dp*sqrt2           &
  & *g12*g22*MassWB*Yt(i1,i2)/5._dp - 306._dp*sqrt2*g24*MassWB*Yt(i1,i2)       &
  &  - 6._dp*sqrt2*g12*TrCYtAYt*Yt(i1,i2)/5._dp - 2._dp*sqrt2*g22*TrCYtAYt*Yt(i1,i2)&
  &  - 9._dp*sqrt2*TrCYtAYtCYtYt*Yt(i1,i2) + 9._dp*g12*MassB*TrCYtYt*Yt(i1,i2)     &
  & /(5._dp*sqrt2) + 3._dp*(g22*MassWB*TrCYtYt*Yt(i1,i2))/sqrt2              &
  &  - 4._dp*sqrt2*TrCYtYtadjYeAYe*Yt(i1,i2) - 12._dp*sqrt2*TrCYtYtadjYzAYz*Yt(i1,i2)&
  &  - 12._dp*sqrt2*TrCYtYtCYtAYt*Yt(i1,i2) - 4._dp*sqrt2*TrYeCYtAYtadjYe*Yt(i1,i2)&
  &  + 3._dp*g12*MassB*TrYtCYt*Yt(i1,i2)/(5._dp*sqrt2) + (g22*MassWB*TrYtCYt*Yt(i1,i2))&
  & /sqrt2 - 3._dp*sqrt2*TrYtCYtAYtCYt*Yt(i1,i2) - 12._dp*sqrt2          &
  & *TrYzCYtAYtadjYz*Yt(i1,i2) + 6._dp*sqrt2*g12*L1*MassB*Conjg(L1)*Yt(i1,i2)      &
  & /5._dp + 2._dp*sqrt2*g22*L1*MassWB*Conjg(L1)*Yt(i1,i2) - 12._dp*sqrt2    &
  & *L1*TrAYdadjYd*Conjg(L1)*Yt(i1,i2) - 4._dp*sqrt2*L1*TrAYeadjYe*Conjg(L1)         &
  & *Yt(i1,i2) - 6._dp*sqrt2*g12*AL1*Conjg(L1)*Yt(i1,i2)/5._dp - 2._dp*sqrt2 &
  & *g22*AL1*Conjg(L1)*Yt(i1,i2) - 12._dp*sqrt2*TrYdadjYd*AL1*Conjg(L1)            &
  & *Yt(i1,i2) - 4._dp*sqrt2*TrYeadjYe*AL1*Conjg(L1)*Yt(i1,i2) - 24._dp*sqrt2  &
  & *L1*AL1*Conjg(L1)**2*Yt(i1,i2) + 12._dp*sqrt2*g12*YtadjYeAYe(i1,i2)            &
  & /5._dp - 6._dp*sqrt2*TrYdadjYd*YtadjYeAYe(i1,i2) - 2._dp*sqrt2             &
  & *TrYeadjYe*YtadjYeAYe(i1,i2) - 6._dp*sqrt2*AbsL1sq*YtadjYeAYe(i1,i2)        &
  &  - 4._dp*sqrt2*YtadjYeAYeadjYeYe(i1,i2) - 6._dp*sqrt2*YtadjYeAYeCYtYt(i1,i2)&
  &  - 12._dp*sqrt2*g12*MassB*YtadjYeYe(i1,i2)/5._dp - 6._dp*sqrt2           &
  & *TrAYdadjYd*YtadjYeYe(i1,i2) - 2._dp*sqrt2*TrAYeadjYe*YtadjYeYe(i1,i2)           &
  &  - 6._dp*sqrt2*AL1*Conjg(L1)*YtadjYeYe(i1,i2) - 4._dp*sqrt2*YtadjYeYeadjYeAYe(i1,i2)&
  &  - 3._dp*sqrt2*YtadjYeYeCYtAYt(i1,i2) - 12._dp*sqrt2*YtadjYzAYdadjYdYz(i1,i2)&
  &  - 24._dp*sqrt2*YtadjYzAYsCYsYz(i1,i2) - 4._dp*sqrt2*g12*YtadjYzAYz(i1,i2)&
  & /5._dp + 32._dp*sqrt2*g32*YtadjYzAYz(i1,i2) - 6._dp*sqrt2*TrYzadjYz*YtadjYzAYz(i1,i2)&
  &  - 12._dp*sqrt2*YtadjYzAYzadjYzYz(i1,i2) - 18._dp*sqrt2*YtadjYzAYzCYtYt(i1,i2)&
  &  - 12._dp*sqrt2*YtadjYzYdadjYdAYz(i1,i2) - 24._dp*sqrt2*YtadjYzYsCYsAYz(i1,i2)&
  &  + 4._dp*sqrt2*g12*MassB*YtadjYzYz(i1,i2)/5._dp - 32._dp*sqrt2           &
  & *g32*MassG*YtadjYzYz(i1,i2) - 6._dp*sqrt2*TrAYzadjYz*YtadjYzYz(i1,i2)          &
  &  - 12._dp*sqrt2*YtadjYzYzadjYzAYz(i1,i2) - 9._dp*sqrt2*YtadjYzYzCYtAYt(i1,i2)&
  &  + 54._dp*sqrt2*g12*YtCYtAYt(i1,i2)/5._dp + 36._dp*sqrt2*g22*YtCYtAYt(i1,i2)&
  &  - 27._dp*TrCYtYt*YtCYtAYt(i1,i2)/(2._dp*sqrt2) - 9._dp*TrYtCYt*YtCYtAYt(i1,i2)  &
  & /(2._dp*sqrt2) - 9._dp*sqrt2*AbsL1sq*YtCYtAYt(i1,i2) - 36._dp*sqrt2&
  & *YtCYtAYtCYtYt(i1,i2) - 6._dp*sqrt2*YtCYtTAYeCYeYt(i1,i2) - 18._dp*sqrt2   &
  & *YtCYtTAYzCYzYt(i1,i2) - 6._dp*sqrt2*YtCYtTYeCYeAYt(i1,i2) - 18._dp*sqrt2  &
  & *YtCYtTYzCYzAYt(i1,i2) - 72._dp*sqrt2*g12*MassB*YtCYtYt(i1,i2)/5._dp -         &
  &  48._dp*sqrt2*g22*MassWB*YtCYtYt(i1,i2) - 12._dp*sqrt2*TrCYtAYt*YtCYtYt(i1,i2)&
  &  - 12._dp*sqrt2*AL1*Conjg(L1)*YtCYtYt(i1,i2) - 27._dp*sqrt2*YtCYtYtCYtAYt(i1,i2)&
  &  + 261._dp*g14*AYt(i1,i2)/(5._dp*sqrt2) + 57._dp*sqrt2*g12*g22*AYt(i1,i2)&
  & /5._dp + 153._dp*(g24*AYt(i1,i2))/sqrt2 - 9._dp*g12*TrCYtYt*AYt(i1,i2)       &
  & /(10._dp*sqrt2) - 3._dp*g22*TrCYtYt*AYt(i1,i2)/(2._dp*sqrt2)             &
  &  - 9._dp*(TrCYtYtCYtYt*AYt(i1,i2))/sqrt2 - 2._dp*sqrt2*TrYeCYtYtadjYe*AYt(i1,i2)&
  &  - 3._dp*g12*TrYtCYt*AYt(i1,i2)/(10._dp*sqrt2) - (g22*TrYtCYt*AYt(i1,i2))    &
  & /(2._dp*sqrt2) - 3._dp*(TrYtCYtYtCYt*AYt(i1,i2))/sqrt2 - 6._dp*sqrt2 &
  & *TrYzCYtYtadjYz*AYt(i1,i2) - 3._dp*sqrt2*g12*AbsL1sq*AYt(i1,i2)           &
  & /5._dp - sqrt2*g22*AbsL1sq*AYt(i1,i2) - 6._dp*sqrt2*L1*TrYdadjYd*Conjg(L1)&
  & *AYt(i1,i2) - 2._dp*sqrt2*L1*TrYeadjYe*Conjg(L1)*AYt(i1,i2) - 6._dp*sqrt2  &
  & *L1**2*Conjg(L1)**2*AYt(i1,i2))/sqrt2
   End Do
  End Do


  DAYt = oo16pi2*( betaAYt1 + oo16pi2 * betaAYt2 )


  Else
  DAYt = oo16pi2* betaAYt1
  End If


  !--------------------
  ! AYs
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betaAYs1(i1,i2) = (4._dp*sqrt2*AYdadjYdYs(i1,i2) + 2._dp*sqrt2            &
  & *AYsCYdTYd(i1,i2) + 12._dp*sqrt2*AYsCYsYs(i1,i2) + 2._dp*sqrt2             &
  & *AYsCYzTYz(i1,i2) + 4._dp*sqrt2*AYzadjYzYs(i1,i2) + 2._dp*sqrt2            &
  & *YdadjYdAYs(i1,i2) + 8._dp*sqrt2*g12*MassB*Ys(i1,i2)/5._dp + 24._dp*sqrt2&
  & *g32*MassG*Ys(i1,i2) + 2._dp*sqrt2*TrCYsAYs*Ys(i1,i2) + 4._dp*sqrt2      &
  & *YsCYdTAYd(i1,i2) + 12._dp*sqrt2*YsCYsAYs(i1,i2) + 4._dp*sqrt2             &
  & *YsCYzTAYz(i1,i2) + 2._dp*sqrt2*YzadjYzAYs(i1,i2) - 4._dp*sqrt2            &
  & *g12*AYs(i1,i2)/5._dp - 12._dp*sqrt2*g32*AYs(i1,i2) + 3._dp*TrCYsYs*AYs(i1,i2)&
  & /(2._dp*sqrt2) + (TrYsCYs*AYs(i1,i2))/(2._dp*sqrt2))/sqrt2
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betaAYs2(i1,i2) = (-4._dp*sqrt2*AYdadjYdYdadjYdYs(i1,i2) + 4._dp*sqrt2    &
  & *g12*AYdadjYdYs(i1,i2)/5._dp + 12._dp*sqrt2*g22*AYdadjYdYs(i1,i2)            &
  &  - 12._dp*sqrt2*TrYdadjYd*AYdadjYdYs(i1,i2) - 4._dp*sqrt2*TrYeadjYe*AYdadjYdYs(i1,i2)&
  &  - 4._dp*sqrt2*AYdadjYuYuadjYdYs(i1,i2) + 2._dp*sqrt2*g12*AYsCYdTYd(i1,i2)&
  & /5._dp + 6._dp*sqrt2*g22*AYsCYdTYd(i1,i2) - 6._dp*sqrt2*TrYdadjYd*AYsCYdTYd(i1,i2)&
  &  - 2._dp*sqrt2*TrYeadjYe*AYsCYdTYd(i1,i2) - 2._dp*sqrt2*AYsCYdTYdCYdTYd(i1,i2)&
  &  - 16._dp*sqrt2*AYsCYdTYdCYsYs(i1,i2) - 2._dp*sqrt2*AYsCYdTYuCYuTYd(i1,i2) &
  &  - 8._dp*sqrt2*AYsCYsYdadjYdYs(i1,i2) + 32._dp*sqrt2*g12*AYsCYsYs(i1,i2) &
  & /5._dp + 80._dp*sqrt2*g32*AYsCYsYs(i1,i2) - 9._dp*sqrt2*TrCYsYs*AYsCYsYs(i1,i2)&
  &  - 3._dp*sqrt2*TrYsCYs*AYsCYsYs(i1,i2) - 48._dp*sqrt2*AYsCYsYsCYsYs(i1,i2) &
  &  - 8._dp*sqrt2*AYsCYsYzadjYzYs(i1,i2) - 2._dp*sqrt2*AYsCYzTYeCYeTYz(i1,i2) &
  &  + 2._dp*sqrt2*g12*AYsCYzTYz(i1,i2)/5._dp + 6._dp*sqrt2*g22*AYsCYzTYz(i1,i2)&
  &  - 2._dp*sqrt2*TrYzadjYz*AYsCYzTYz(i1,i2) - 16._dp*sqrt2*AYsCYzTYzCYsYs(i1,i2)&
  &  - 6._dp*sqrt2*AYsCYzTYzCYzTYz(i1,i2) - 6._dp*sqrt2*AYsCYzYtCYtTYz(i1,i2)  &
  &  - 4._dp*sqrt2*AYzadjYeYeadjYzYs(i1,i2) + 4._dp*sqrt2*g12*AYzadjYzYs(i1,i2)&
  & /5._dp + 12._dp*sqrt2*g22*AYzadjYzYs(i1,i2) - 4._dp*sqrt2*TrYzadjYz*AYzadjYzYs(i1,i2)&
  &  - 12._dp*sqrt2*AYzadjYzYzadjYzYs(i1,i2) - 12._dp*sqrt2*AYzCYtYtadjYzYs(i1,i2)&
  &  - 12._dp*sqrt2*L1*AYdadjYdYs(i1,i2)*Conjg(L1) - 6._dp*sqrt2               &
  & *L1*AYsCYdTYd(i1,i2)*Conjg(L1) - 4._dp*sqrt2*YdadjYdAYdadjYdYs(i1,i2)            &
  &  + 2._dp*sqrt2*g12*YdadjYdAYs(i1,i2)/5._dp + 6._dp*sqrt2*g22*YdadjYdAYs(i1,i2)&
  &  - 6._dp*sqrt2*TrYdadjYd*YdadjYdAYs(i1,i2) - 2._dp*sqrt2*TrYeadjYe*YdadjYdAYs(i1,i2)&
  &  - 6._dp*sqrt2*AbsL1sq*YdadjYdAYs(i1,i2) - 2._dp*sqrt2*YdadjYdYdadjYdAYs(i1,i2)&
  &  - 4._dp*sqrt2*g12*MassB*YdadjYdYs(i1,i2)/5._dp - 12._dp*sqrt2           &
  & *g22*MassWB*YdadjYdYs(i1,i2) - 12._dp*sqrt2*TrAYdadjYd*YdadjYdYs(i1,i2)        &
  &  - 4._dp*sqrt2*TrAYeadjYe*YdadjYdYs(i1,i2) - 12._dp*sqrt2*AL1*Conjg(L1)    &
  & *YdadjYdYs(i1,i2) - 4._dp*sqrt2*YdadjYuAYuadjYdYs(i1,i2) - 2._dp*sqrt2     &
  & *YdadjYuYuadjYdAYs(i1,i2) - 224._dp*sqrt2*g14*MassB*Ys(i1,i2)/5._dp -          &
  &  256._dp*sqrt2*g12*g32*MassB*Ys(i1,i2)/15._dp - 256._dp*sqrt2          &
  & *g12*g32*MassG*Ys(i1,i2)/15._dp - 1280._dp*sqrt2*g34*MassG*Ys(i1,i2)       &
  & /3._dp - 8._dp*sqrt2*g12*TrCYsAYs*Ys(i1,i2)/15._dp - 8._dp*sqrt2         &
  & *g32*TrCYsAYs*Ys(i1,i2)/3._dp - 12._dp*sqrt2*TrCYsAYsCYsYs*Ys(i1,i2)           &
  &  - 6._dp*sqrt2*TrCYsYdadjYdAYs*Ys(i1,i2) + 2._dp*sqrt2*g12*MassB*TrCYsYs*Ys(i1,i2)&
  & /5._dp + 2._dp*sqrt2*g32*MassG*TrCYsYs*Ys(i1,i2) - 16._dp*sqrt2          &
  & *TrCYsYsCYsAYs*Ys(i1,i2) - 6._dp*sqrt2*TrCYsYzadjYzAYs*Ys(i1,i2) -               &
  &  2._dp*sqrt2*TrYdadjYdAYsCYs*Ys(i1,i2) + 2._dp*sqrt2*g12*MassB*TrYsCYs*Ys(i1,i2)&
  & /15._dp + 2._dp*sqrt2*g32*MassG*TrYsCYs*Ys(i1,i2)/3._dp - 8._dp*sqrt2    &
  & *TrYsCYsAYdadjYd*Ys(i1,i2) - 4._dp*sqrt2*TrYsCYsAYsCYs*Ys(i1,i2) -               &
  &  8._dp*sqrt2*TrYsCYsAYzadjYz*Ys(i1,i2) - 2._dp*sqrt2*TrYzadjYzAYsCYs*Ys(i1,i2)&
  &  + 4._dp*sqrt2*g12*YsCYdTAYd(i1,i2)/5._dp + 12._dp*sqrt2*g22*YsCYdTAYd(i1,i2)&
  &  - 12._dp*sqrt2*TrYdadjYd*YsCYdTAYd(i1,i2) - 4._dp*sqrt2*TrYeadjYe*YsCYdTAYd(i1,i2)&
  &  - 12._dp*sqrt2*AbsL1sq*YsCYdTAYd(i1,i2) - 4._dp*sqrt2*YsCYdTAYdCYdTYd(i1,i2)&
  &  - 16._dp*sqrt2*YsCYdTAYdCYsYs(i1,i2) - 4._dp*sqrt2*YsCYdTAYuCYuTYd(i1,i2) &
  &  - 4._dp*sqrt2*g12*MassB*YsCYdTYd(i1,i2)/5._dp - 12._dp*sqrt2            &
  & *g22*MassWB*YsCYdTYd(i1,i2) - 12._dp*sqrt2*TrAYdadjYd*YsCYdTYd(i1,i2)          &
  &  - 4._dp*sqrt2*TrAYeadjYe*YsCYdTYd(i1,i2) - 12._dp*sqrt2*AL1*Conjg(L1)     &
  & *YsCYdTYd(i1,i2) - 4._dp*sqrt2*YsCYdTYdCYdTAYd(i1,i2) - 8._dp*sqrt2        &
  & *YsCYdTYdCYsAYs(i1,i2) - 4._dp*sqrt2*YsCYdTYuCYuTAYd(i1,i2) - 16._dp*sqrt2 &
  & *YsCYsAYdadjYdYs(i1,i2) + 32._dp*sqrt2*g12*YsCYsAYs(i1,i2)/5._dp +             &
  &  80._dp*sqrt2*g32*YsCYsAYs(i1,i2) - 9._dp*sqrt2*TrCYsYs*YsCYsAYs(i1,i2)  &
  &  - 3._dp*sqrt2*TrYsCYs*YsCYsAYs(i1,i2) - 64._dp*sqrt2*YsCYsAYsCYsYs(i1,i2) &
  &  - 16._dp*sqrt2*YsCYsAYzadjYzYs(i1,i2) - 16._dp*sqrt2*YsCYsYdadjYdAYs(i1,i2)&
  &  - 128._dp*sqrt2*g12*MassB*YsCYsYs(i1,i2)/15._dp - 320._dp*sqrt2         &
  & *g32*MassG*YsCYsYs(i1,i2)/3._dp - 16._dp*sqrt2*TrCYsAYs*YsCYsYs(i1,i2)         &
  &  - 48._dp*sqrt2*YsCYsYsCYsAYs(i1,i2) - 16._dp*sqrt2*YsCYsYzadjYzAYs(i1,i2) &
  &  - 12._dp*sqrt2*YsCYzAYtCYtTYz(i1,i2) - 4._dp*sqrt2*YsCYzTAYeCYeTYz(i1,i2) &
  &  + 4._dp*sqrt2*g12*YsCYzTAYz(i1,i2)/5._dp + 12._dp*sqrt2*g22*YsCYzTAYz(i1,i2)&
  &  - 4._dp*sqrt2*TrYzadjYz*YsCYzTAYz(i1,i2) - 16._dp*sqrt2*YsCYzTAYzCYsYs(i1,i2)&
  &  - 12._dp*sqrt2*YsCYzTAYzCYzTYz(i1,i2) - 4._dp*sqrt2*YsCYzTYeCYeTAYz(i1,i2)&
  &  - 4._dp*sqrt2*g12*MassB*YsCYzTYz(i1,i2)/5._dp - 12._dp*sqrt2            &
  & *g22*MassWB*YsCYzTYz(i1,i2) - 4._dp*sqrt2*TrAYzadjYz*YsCYzTYz(i1,i2)           &
  &  - 8._dp*sqrt2*YsCYzTYzCYsAYs(i1,i2) - 12._dp*sqrt2*YsCYzTYzCYzTAYz(i1,i2) &
  &  - 12._dp*sqrt2*YsCYzYtCYtTAYz(i1,i2) - 4._dp*sqrt2*YzadjYeAYeadjYzYs(i1,i2)&
  &  - 2._dp*sqrt2*YzadjYeYeadjYzAYs(i1,i2) + 2._dp*sqrt2*g12*YzadjYzAYs(i1,i2)&
  & /5._dp + 6._dp*sqrt2*g22*YzadjYzAYs(i1,i2) - 2._dp*sqrt2*TrYzadjYz*YzadjYzAYs(i1,i2)&
  &  - 12._dp*sqrt2*YzadjYzAYzadjYzYs(i1,i2) - 4._dp*sqrt2*g12*MassB*YzadjYzYs(i1,i2)&
  & /5._dp - 12._dp*sqrt2*g22*MassWB*YzadjYzYs(i1,i2) - 4._dp*sqrt2          &
  & *TrAYzadjYz*YzadjYzYs(i1,i2) - 6._dp*sqrt2*YzadjYzYzadjYzAYs(i1,i2)              &
  &  - 12._dp*sqrt2*YzCYtAYtadjYzYs(i1,i2) - 6._dp*sqrt2*YzCYtYtadjYzAYs(i1,i2)&
  &  + 56._dp*sqrt2*g14*AYs(i1,i2)/5._dp + 128._dp*sqrt2*g12*g32*AYs(i1,i2)&
  & /15._dp + 320._dp*sqrt2*g34*AYs(i1,i2)/3._dp - (sqrt2*g12*TrCYsYs*AYs(i1,i2))&
  & /5._dp - sqrt2*g32*TrCYsYs*AYs(i1,i2) - 6._dp*sqrt2*TrCYsYsCYsYs*AYs(i1,i2)&
  &  - (sqrt2*g12*TrYsCYs*AYs(i1,i2))/15._dp - (sqrt2*g32*TrYsCYs*AYs(i1,i2))&
  & /3._dp - 4._dp*sqrt2*TrYsCYsYdadjYd*AYs(i1,i2) - 2._dp*sqrt2               &
  & *TrYsCYsYsCYs*AYs(i1,i2) - 4._dp*sqrt2*TrYsCYsYzadjYz*AYs(i1,i2))/sqrt2
   End Do
  End Do


  DAYs = oo16pi2*( betaAYs1 + oo16pi2 * betaAYs2 )


  Else
  DAYs = oo16pi2* betaAYs1
  End If


  !--------------------
  ! AYz
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betaAYz1(i1,i2) = 4._dp*AYdadjYdYz(i1,i2) + 8._dp*AYsCYsYz(i1,i2) + AYzadjYeYe(i1,i2) &
  &  + 7._dp*AYzadjYzYz(i1,i2) + 3._dp*AYzCYtYt(i1,i2) + 2._dp*YdadjYdAYz(i1,i2)           &
  &  + 4._dp*YsCYsAYz(i1,i2) + 14._dp*g12*MassB*Yz(i1,i2)/15._dp + 32._dp*g32*MassG*Yz(i1,i2)&
  & /3._dp + 6._dp*g22*MassWB*Yz(i1,i2) + 2._dp*TrAYzadjYz*Yz(i1,i2) + 2._dp*YzadjYeAYe(i1,i2)&
  &  + 8._dp*YzadjYzAYz(i1,i2) + 6._dp*YzCYtAYt(i1,i2) - 7._dp*g12*AYz(i1,i2)            &
  & /15._dp - 3._dp*g22*AYz(i1,i2) - 16._dp*g32*AYz(i1,i2)/3._dp + TrYzadjYz*AYz(i1,i2)
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betaAYz2(i1,i2) = -4._dp*AYdadjYdYdadjYdYz(i1,i2) + 4._dp*g12*AYdadjYdYz(i1,i2)     &
  & /5._dp + 12._dp*g22*AYdadjYdYz(i1,i2) - 12._dp*TrYdadjYd*AYdadjYdYz(i1,i2)           &
  &  - 4._dp*TrYeadjYe*AYdadjYdYz(i1,i2) - 4._dp*AYdadjYuYuadjYdYz(i1,i2) - 16._dp*AYsCYdTYdCYsYz(i1,i2)&
  &  - 32._dp*AYsCYsYsCYsYz(i1,i2) + 64._dp*g12*AYsCYsYz(i1,i2)/15._dp + 160._dp*g32*AYsCYsYz(i1,i2)&
  & /3._dp - 6._dp*TrCYsYs*AYsCYsYz(i1,i2) - 2._dp*TrYsCYs*AYsCYsYz(i1,i2) -               &
  &  16._dp*AYsCYzTYzCYsYz(i1,i2) + 6._dp*g12*AYzadjYeYe(i1,i2)/5._dp - 3._dp*TrYdadjYd*AYzadjYeYe(i1,i2)&
  &  - TrYeadjYe*AYzadjYeYe(i1,i2) - 2._dp*AYzadjYeYeadjYeYe(i1,i2) - 4._dp*AYzadjYeYeadjYzYz(i1,i2)&
  &  - 6._dp*AYzadjYzYdadjYdYz(i1,i2) - 12._dp*AYzadjYzYsCYsYz(i1,i2) + 2._dp*g12*AYzadjYzYz(i1,i2)&
  & /5._dp + 12._dp*g22*AYzadjYzYz(i1,i2) + 16._dp*g32*AYzadjYzYz(i1,i2)               &
  &  - 7._dp*TrYzadjYz*AYzadjYzYz(i1,i2) - 18._dp*AYzadjYzYzadjYzYz(i1,i2) -               &
  &  3._dp*AYzCYtTYeCYeYt(i1,i2) - 9._dp*AYzCYtTYzCYzYt(i1,i2) + 18._dp*g12*AYzCYtYt(i1,i2)&
  & /5._dp + 12._dp*g22*AYzCYtYt(i1,i2) - 9._dp*TrCYtYt*AYzCYtYt(i1,i2)/4._dp -          &
  &  3._dp*TrYtCYt*AYzCYtYt(i1,i2)/4._dp - 12._dp*AYzCYtYtadjYzYz(i1,i2) - 9._dp*AYzCYtYtCYtYt(i1,i2)&
  &  - 12._dp*L1*AYdadjYdYz(i1,i2)*Conjg(L1) - 3._dp*L1*AYzadjYeYe(i1,i2)*Conjg(L1)        &
  &  - 3._dp*L1*AYzCYtYt(i1,i2)*Conjg(L1) - 4._dp*YdadjYdAYdadjYdYz(i1,i2) +               &
  &  2._dp*g12*YdadjYdAYz(i1,i2)/5._dp + 6._dp*g22*YdadjYdAYz(i1,i2) - 6._dp*TrYdadjYd*YdadjYdAYz(i1,i2)&
  &  - 2._dp*TrYeadjYe*YdadjYdAYz(i1,i2) - 6._dp*AbsL1sq*YdadjYdAYz(i1,i2)            &
  &  - 2._dp*YdadjYdYdadjYdAYz(i1,i2) - 4._dp*g12*MassB*YdadjYdYz(i1,i2)/5._dp -         &
  &  12._dp*g22*MassWB*YdadjYdYz(i1,i2) - 12._dp*TrAYdadjYd*YdadjYdYz(i1,i2)             &
  &  - 4._dp*TrAYeadjYe*YdadjYdYz(i1,i2) - 12._dp*AL1*Conjg(L1)*YdadjYdYz(i1,i2)           &
  &  - 4._dp*YdadjYuAYuadjYdYz(i1,i2) - 2._dp*YdadjYuYuadjYdAYz(i1,i2) - 16._dp*YsCYdTAYdCYsYz(i1,i2)&
  &  - 8._dp*YsCYdTYdCYsAYz(i1,i2) - 32._dp*YsCYsAYsCYsYz(i1,i2) + 32._dp*g12*YsCYsAYz(i1,i2)&
  & /15._dp + 80._dp*g32*YsCYsAYz(i1,i2)/3._dp - 3._dp*TrCYsYs*YsCYsAYz(i1,i2)           &
  &  - TrYsCYs*YsCYsAYz(i1,i2) - 16._dp*YsCYsYsCYsAYz(i1,i2) - 64._dp*g12*MassB*YsCYsYz(i1,i2)&
  & /15._dp - 160._dp*g32*MassG*YsCYsYz(i1,i2)/3._dp - 8._dp*TrCYsAYs*YsCYsYz(i1,i2)     &
  &  - 16._dp*YsCYzTAYzCYsYz(i1,i2) - 8._dp*YsCYzTYzCYsAYz(i1,i2) - 1162._dp*g14*MassB*Yz(i1,i2)&
  & /45._dp - 2._dp*g12*g22*MassB*Yz(i1,i2) - 16._dp*g12*g32*MassB*Yz(i1,i2)       &
  & /9._dp - 16._dp*g12*g32*MassG*Yz(i1,i2)/9._dp - 16._dp*g22*g32*MassG*Yz(i1,i2) &
  &  - 1280._dp*g34*MassG*Yz(i1,i2)/9._dp - 2._dp*g12*g22*MassWB*Yz(i1,i2)           &
  &  - 114._dp*g24*MassWB*Yz(i1,i2) - 16._dp*g22*g32*MassWB*Yz(i1,i2) +              &
  &  4._dp*g12*TrAYzadjYz*Yz(i1,i2)/5._dp - 4._dp*TrCYsAYzadjYzYs*Yz(i1,i2)              &
  &  - 5._dp*TrCYsYzadjYzAYs*Yz(i1,i2) - 3._dp*TrCYtAYtadjYzYz*Yz(i1,i2) - 4._dp*TrCYtYtadjYzAYz*Yz(i1,i2)&
  &  - 4._dp*TrYdadjYdAYzadjYz*Yz(i1,i2) - 2._dp*TrYeadjYzAYzadjYe*Yz(i1,i2)               &
  &  - 4._dp*TrYsCYsAYzadjYz*Yz(i1,i2) - 2._dp*TrYtadjYzAYzCYt*Yz(i1,i2) - 2._dp*TrYzadjYeAYeadjYz*Yz(i1,i2)&
  &  - 4._dp*g12*MassB*TrYzadjYz*Yz(i1,i2)/5._dp - 4._dp*TrYzadjYzAYdadjYd*Yz(i1,i2)     &
  &  - 3._dp*TrYzadjYzAYsCYs*Yz(i1,i2) - 20._dp*TrYzadjYzAYzadjYz*Yz(i1,i2) -              &
  &  3._dp*TrYzCYtAYtadjYz*Yz(i1,i2) + 12._dp*g12*YzadjYeAYe(i1,i2)/5._dp -              &
  &  6._dp*TrYdadjYd*YzadjYeAYe(i1,i2) - 2._dp*TrYeadjYe*YzadjYeAYe(i1,i2) -               &
  &  6._dp*AbsL1sq*YzadjYeAYe(i1,i2) - 4._dp*YzadjYeAYeadjYeYe(i1,i2) -               &
  &  4._dp*YzadjYeAYeadjYzYz(i1,i2) - 12._dp*g12*MassB*YzadjYeYe(i1,i2)/5._dp -          &
  &  6._dp*TrAYdadjYd*YzadjYeYe(i1,i2) - 2._dp*TrAYeadjYe*YzadjYeYe(i1,i2) -               &
  &  6._dp*AL1*Conjg(L1)*YzadjYeYe(i1,i2) - 4._dp*YzadjYeYeadjYeAYe(i1,i2) -               &
  &  2._dp*YzadjYeYeadjYzAYz(i1,i2) - 12._dp*YzadjYzAYdadjYdYz(i1,i2) - 24._dp*YzadjYzAYsCYsYz(i1,i2)&
  &  - 2._dp*g12*YzadjYzAYz(i1,i2)/5._dp + 6._dp*g22*YzadjYzAYz(i1,i2) +               &
  &  32._dp*g32*YzadjYzAYz(i1,i2) - 8._dp*TrYzadjYz*YzadjYzAYz(i1,i2) - 24._dp*YzadjYzAYzadjYzYz(i1,i2)&
  &  - 12._dp*YzadjYzYdadjYdAYz(i1,i2) - 24._dp*YzadjYzYsCYsAYz(i1,i2) - 32._dp*g32*MassG*YzadjYzYz(i1,i2)&
  &  - 12._dp*g22*MassWB*YzadjYzYz(i1,i2) - 10._dp*TrAYzadjYz*YzadjYzYz(i1,i2)           &
  &  - 18._dp*YzadjYzYzadjYzAYz(i1,i2) + 36._dp*g12*YzCYtAYt(i1,i2)/5._dp +              &
  &  24._dp*g22*YzCYtAYt(i1,i2) - 9._dp*TrCYtYt*YzCYtAYt(i1,i2)/2._dp - 3._dp*TrYtCYt*YzCYtAYt(i1,i2)&
  & /2._dp - 6._dp*AbsL1sq*YzCYtAYt(i1,i2) - 12._dp*YzCYtAYtadjYzYz(i1,i2)            &
  &  - 18._dp*YzCYtAYtCYtYt(i1,i2) - 6._dp*YzCYtTAYeCYeYt(i1,i2) - 18._dp*YzCYtTAYzCYzYt(i1,i2)&
  &  - 6._dp*YzCYtTYeCYeAYt(i1,i2) - 18._dp*YzCYtTYzCYzAYt(i1,i2) - 36._dp*g12*MassB*YzCYtYt(i1,i2)&
  & /5._dp - 24._dp*g22*MassWB*YzCYtYt(i1,i2) - 6._dp*TrCYtAYt*YzCYtYt(i1,i2)            &
  &  - 6._dp*AL1*Conjg(L1)*YzCYtYt(i1,i2) - 6._dp*YzCYtYtadjYzAYz(i1,i2) - 18._dp*YzCYtYtCYtAYt(i1,i2)&
  &  + 581._dp*g14*AYz(i1,i2)/90._dp + g12*g22*AYz(i1,i2) + 57._dp*g24*AYz(i1,i2)  &
  & /2._dp + 8._dp*g12*g32*AYz(i1,i2)/9._dp + 8._dp*g22*g32*AYz(i1,i2)             &
  &  + 320._dp*g34*AYz(i1,i2)/9._dp - 2._dp*TrCYsYzadjYzYs*AYz(i1,i2) - 3._dp*TrCYtYtadjYzYz*AYz(i1,i2)&
  & /2._dp - TrYeadjYzYzadjYe*AYz(i1,i2) - (TrYsCYsYzadjYz*AYz(i1,i2))/2._dp -             &
  &  TrYtadjYzYzCYt*AYz(i1,i2) + 2._dp*g12*TrYzadjYz*AYz(i1,i2)/5._dp - 2._dp*TrYzadjYzYdadjYd*AYz(i1,i2)&
  &  - 3._dp*TrYzadjYzYsCYs*AYz(i1,i2)/2._dp - 5._dp*TrYzadjYzYzadjYz*AYz(i1,i2)           &
  &  - (TrYzCYtYtadjYz*AYz(i1,i2))/2._dp
   End Do
  End Do


  DAYz = oo16pi2*( betaAYz1 + oo16pi2 * betaAYz2 )


  Else
  DAYz = oo16pi2* betaAYz1
  End If


  !--------------------
  ! AL1
  !--------------------

  betaAL11 = (18._dp*sqrt2*g12*L1*MassB/5._dp + 14._dp*sqrt2              &
  & *g22*L1*MassWB + 12._dp*sqrt2*L1*TrAYdadjYd + 4._dp*sqrt2*L1*TrAYeadjYe +&
  &  2._dp*sqrt2*L1*TrCYtAYt - 9._dp*sqrt2*g12*AL1/5._dp - 7._dp*sqrt2 &
  & *g22*AL1 + 3._dp*TrCYtYt*AL1/(2._dp*sqrt2) + 6._dp*sqrt2*TrYdadjYd*AL1 + &
  &  2._dp*sqrt2*TrYeadjYe*AL1 + (TrYtCYt*AL1)/(2._dp*sqrt2) + 21._dp*sqrt2&
  & *L1*AL1*Conjg(L1))/sqrt2


  If (TwoLoopRGE) Then
  betaAL12 = (-522._dp*sqrt2*g14*L1*MassB/5._dp - 114._dp*sqrt2           &
  & *g12*g22*L1*MassB/5._dp - 114._dp*sqrt2*g12*g22*L1*MassWB/5._dp -        &
  &  306._dp*sqrt2*g24*L1*MassWB - 8._dp*sqrt2*g12*L1*TrAYdadjYd/5._dp +   &
  &  64._dp*sqrt2*g32*L1*TrAYdadjYd + 24._dp*sqrt2*g12*L1*TrAYeadjYe/5._dp -&
  &  24._dp*sqrt2*L1*TrCYsAYdadjYdYs - 30._dp*sqrt2*L1*TrCYsYdadjYdAYs -       &
  &  6._dp*sqrt2*g12*L1*TrCYtAYt/5._dp - 2._dp*sqrt2*g22*L1*TrCYtAYt -     &
  &  6._dp*sqrt2*L1*TrCYtAYtadjYeYe - 9._dp*sqrt2*L1*TrCYtAYtCYtYt +           &
  &  9._dp*g12*L1*MassB*TrCYtYt/(5._dp*sqrt2) + 3._dp*(g22*L1*MassWB*TrCYtYt)    &
  & /sqrt2 - 12._dp*sqrt2*L1*TrCYtYtadjYeAYe - 12._dp*sqrt2              &
  & *L1*TrCYtYtadjYzAYz - 12._dp*sqrt2*L1*TrCYtYtCYtAYt + 8._dp*sqrt2          &
  & *g12*L1*MassB*TrYdadjYd/5._dp - 64._dp*sqrt2*g32*L1*MassG*TrYdadjYd -        &
  &  72._dp*sqrt2*L1*TrYdadjYdAYdadjYd - 18._dp*sqrt2*L1*TrYdadjYdAYsCYs -     &
  &  24._dp*sqrt2*L1*TrYdadjYdAYzadjYz - 12._dp*sqrt2*L1*TrYdadjYuAYuadjYd -   &
  &  24._dp*sqrt2*g12*L1*MassB*TrYeadjYe/5._dp - 24._dp*sqrt2*L1*TrYeadjYeAYeadjYe -&
  &  12._dp*sqrt2*L1*TrYeadjYzAYzadjYe - 10._dp*sqrt2*L1*TrYeCYtAYtadjYe -     &
  &  24._dp*sqrt2*L1*TrYsCYsAYdadjYd - 4._dp*sqrt2*L1*TrYtadjYeAYeCYt +        &
  &  3._dp*g12*L1*MassB*TrYtCYt/(5._dp*sqrt2) + (g22*L1*MassWB*TrYtCYt)          &
  & /sqrt2 - 3._dp*sqrt2*L1*TrYtCYtAYtCYt - 12._dp*sqrt2*L1*TrYuadjYdAYdadjYu -&
  &  12._dp*sqrt2*L1*TrYzadjYeAYeadjYz - 24._dp*sqrt2*L1*TrYzadjYzAYdadjYd -   &
  &  12._dp*sqrt2*L1*TrYzCYtAYtadjYz + 261._dp*g14*AL1/(5._dp*sqrt2)         &
  &  + 57._dp*sqrt2*g12*g22*AL1/5._dp + 153._dp*(g24*AL1)/sqrt2          &
  &  - 12._dp*sqrt2*TrCYsYdadjYdYs*AL1 - 9._dp*g12*TrCYtYt*AL1/(10._dp*sqrt2)&
  &  - 3._dp*g22*TrCYtYt*AL1/(2._dp*sqrt2) - 3._dp*sqrt2*TrCYtYtadjYeYe*AL1 -&
  &  9._dp*(TrCYtYtCYtYt*AL1)/sqrt2 - 4._dp*sqrt2*g12*TrYdadjYd*AL1/5._dp +  &
  &  32._dp*sqrt2*g32*TrYdadjYd*AL1 - 18._dp*sqrt2*TrYdadjYdYdadjYd*AL1 -    &
  &  9._dp*sqrt2*TrYdadjYdYsCYs*AL1 - 12._dp*sqrt2*TrYdadjYdYzadjYz*AL1 +      &
  &  12._dp*sqrt2*g12*TrYeadjYe*AL1/5._dp - 6._dp*sqrt2*TrYeadjYeYeadjYe*AL1 -&
  &  3._dp*sqrt2*TrYeCYtYtadjYe*AL1 - 3._dp*sqrt2*TrYsCYsYdadjYd*AL1 -         &
  &  2._dp*sqrt2*TrYtadjYeYeCYt*AL1 - 3._dp*g12*TrYtCYt*AL1/(10._dp*sqrt2)   &
  &  - (g22*TrYtCYt*AL1)/(2._dp*sqrt2) - 3._dp*(TrYtCYtYtCYt*AL1)/sqrt2      &
  &  - 6._dp*sqrt2*TrYuadjYdYdadjYu*AL1 - 6._dp*sqrt2*TrYzadjYeYeadjYz*AL1 -   &
  &  6._dp*sqrt2*TrYzCYtYtadjYz*AL1 - 66._dp*sqrt2*g12*L1**2*MassB*Conjg(L1) &
  & /5._dp - 46._dp*sqrt2*g22*L1**2*MassWB*Conjg(L1) - 48._dp*sqrt2          &
  & *L1**2*TrAYdadjYd*Conjg(L1) - 16._dp*sqrt2*L1**2*TrAYeadjYe*Conjg(L1)            &
  &  - 12._dp*sqrt2*L1**2*TrCYtAYt*Conjg(L1) + 99._dp*sqrt2*g12*L1*AL1*Conjg(L1)&
  & /5._dp + 69._dp*sqrt2*g22*L1*AL1*Conjg(L1) - 27._dp*(L1*TrCYtYt*AL1*Conjg(L1)) &
  & /sqrt2 - 72._dp*sqrt2*L1*TrYdadjYd*AL1*Conjg(L1) - 24._dp*sqrt2      &
  & *L1*TrYeadjYe*AL1*Conjg(L1) - 9._dp*(L1*TrYtCYt*AL1*Conjg(L1))/sqrt2             &
  &  - 150._dp*sqrt2*L1**2*AL1*Conjg(L1)**2)/sqrt2


  DAL1 = oo16pi2*( betaAL11 + oo16pi2 * betaAL12 )


  Else
  DAL1 = oo16pi2* betaAL11
  End If


  !--------------------
  ! AL2
  !--------------------

  betaAL21 = (18._dp*g12*L2*MassB/5._dp + 14._dp              &
  & *g22*L2*MassWB + 12._dp*L2*TrAYuadjYu - 9._dp*g12*AL2/5._dp -&
  &  7._dp*g22*AL2 + 6._dp*TrYuadjYu*AL2 + 21._dp    &
  & *L2*AL2*Conjg(L2))


  If (TwoLoopRGE) Then
  betaAL22 = (-522._dp*sqrt2*g14*L2*MassB/5._dp - 114._dp*sqrt2           &
  & *g12*g22*L2*MassB/5._dp - 114._dp*sqrt2*g12*g22*L2*MassWB/5._dp -        &
  &  306._dp*sqrt2*g24*L2*MassWB + 16._dp*sqrt2*g12*L2*TrAYuadjYu/5._dp +  &
  &  64._dp*sqrt2*g32*L2*TrAYuadjYu - 12._dp*sqrt2*L2*TrYdadjYuAYuadjYd -    &
  &  12._dp*sqrt2*L2*TrYuadjYdAYdadjYu - 16._dp*sqrt2*g12*L2*MassB*TrYuadjYu/5._dp -&
  &  64._dp*sqrt2*g32*L2*MassG*TrYuadjYu - 72._dp*sqrt2*L2*TrYuadjYuAYuadjYu +&
  &  261._dp*g14*AL2/(5._dp*sqrt2) + 57._dp*sqrt2*g12*g22*AL2/5._dp +    &
  &  153._dp*(g24*AL2)/sqrt2 - 6._dp*sqrt2*TrYdadjYuYuadjYd*AL2 +            &
  &  8._dp*sqrt2*g12*TrYuadjYu*AL2/5._dp + 32._dp*sqrt2*g32*TrYuadjYu*AL2 -&
  &  18._dp*sqrt2*TrYuadjYuYuadjYu*AL2 - 66._dp*sqrt2*g12*L2**2*MassB*Conjg(L2)&
  & /5._dp - 46._dp*sqrt2*g22*L2**2*MassWB*Conjg(L2) - 48._dp*sqrt2          &
  & *L2**2*TrAYuadjYu*Conjg(L2) + 99._dp*sqrt2*g12*L2*AL2*Conjg(L2)/5._dp +        &
  &  69._dp*sqrt2*g22*L2*AL2*Conjg(L2) - 72._dp*sqrt2*L2*TrYuadjYu*AL2*Conjg(L2)&
  &  - 150._dp*sqrt2*L2**2*AL2*Conjg(L2)**2)/sqrt2


  DAL2 = oo16pi2*( betaAL21 + oo16pi2 * betaAL22 )


  Else
  DAL2 = oo16pi2* betaAL21
  End If


  !--------------------
  ! Amue
  !--------------------

  betaAmue1 = 6._dp*g12*MassB*mue/5._dp + 6._dp*g22*MassWB*mue + 6._dp*TrAYdadjYd*mue +&
  &  2._dp*TrAYeadjYe*mue + 6._dp*TrAYuadjYu*mue - 3._dp*g12*Amue/5._dp - 3._dp*g22*Amue +&
  &  3._dp*TrYdadjYd*Amue + TrYeadjYe*Amue + 3._dp*TrYuadjYu*Amue + 6._dp*mue*AL1*Conjg(L1)&
  &  + 3._dp*L1*Amue*Conjg(L1) + 6._dp*mue*AL2*Conjg(L2) + 3._dp*L2*Amue*Conjg(L2)


  If (TwoLoopRGE) Then
  betaAmue2 = -834._dp*g14*MassB*mue/25._dp - 18._dp*g12*g22*MassB*mue/5._dp -    &
  &  18._dp*g12*g22*MassWB*mue/5._dp - 114._dp*g24*MassWB*mue - 4._dp*g12*TrAYdadjYd*mue/5._dp +&
  &  32._dp*g32*TrAYdadjYd*mue + 12._dp*g12*TrAYeadjYe*mue/5._dp + 8._dp*g12*TrAYuadjYu*mue/5._dp +&
  &  32._dp*g32*TrAYuadjYu*mue - 12._dp*TrCYsAYdadjYdYs*mue - 15._dp*TrCYsYdadjYdAYs*mue -&
  &  3._dp*TrCYtAYtadjYeYe*mue - 4._dp*TrCYtYtadjYeAYe*mue + 4._dp*g12*MassB*TrYdadjYd*mue/5._dp -&
  &  32._dp*g32*MassG*TrYdadjYd*mue - 36._dp*TrYdadjYdAYdadjYd*mue - 9._dp*TrYdadjYdAYsCYs*mue -&
  &  12._dp*TrYdadjYdAYzadjYz*mue - 12._dp*TrYdadjYuAYuadjYd*mue - 12._dp*g12*MassB*TrYeadjYe*mue/5._dp -&
  &  12._dp*TrYeadjYeAYeadjYe*mue - 6._dp*TrYeadjYzAYzadjYe*mue - 3._dp*TrYeCYtAYtadjYe*mue -&
  &  12._dp*TrYsCYsAYdadjYd*mue - 2._dp*TrYtadjYeAYeCYt*mue - 12._dp*TrYuadjYdAYdadjYu*mue -&
  &  8._dp*g12*MassB*TrYuadjYu*mue/5._dp - 32._dp*g32*MassG*TrYuadjYu*mue -            &
  &  36._dp*TrYuadjYuAYuadjYu*mue - 6._dp*TrYzadjYeAYeadjYz*mue - 12._dp*TrYzadjYzAYdadjYd*mue +&
  &  417._dp*g14*Amue/50._dp + 9._dp*g12*g22*Amue/5._dp + 57._dp*g24*Amue/2._dp -  &
  &  6._dp*TrCYsYdadjYdYs*Amue - 3._dp*TrCYtYtadjYeYe*Amue/2._dp - 2._dp*g12*TrYdadjYd*Amue/5._dp +&
  &  16._dp*g32*TrYdadjYd*Amue - 9._dp*TrYdadjYdYdadjYd*Amue - 9._dp*TrYdadjYdYsCYs*Amue/2._dp -&
  &  6._dp*TrYdadjYdYzadjYz*Amue - 3._dp*TrYdadjYuYuadjYd*Amue + 6._dp*g12*TrYeadjYe*Amue/5._dp -&
  &  3._dp*TrYeadjYeYeadjYe*Amue - (TrYeCYtYtadjYe*Amue)/2._dp - 3._dp*TrYsCYsYdadjYd*Amue/2._dp -&
  &  TrYtadjYeYeCYt*Amue - 3._dp*TrYuadjYdYdadjYu*Amue + 4._dp*g12*TrYuadjYu*Amue/5._dp +&
  &  16._dp*g32*TrYuadjYu*Amue - 9._dp*TrYuadjYuYuadjYu*Amue - 3._dp*TrYzadjYeYeadjYz*Amue -&
  &  36._dp*g12*L1*MassB*mue*Conjg(L1)/5._dp - 24._dp*g22*L1*MassWB*mue*Conjg(L1)      &
  &  - 18._dp*L1*TrAYdadjYd*mue*Conjg(L1) - 6._dp*L1*TrAYeadjYe*mue*Conjg(L1)              &
  &  - 6._dp*L1*TrCYtAYt*mue*Conjg(L1) + 36._dp*g12*mue*AL1*Conjg(L1)/5._dp +            &
  &  24._dp*g22*mue*AL1*Conjg(L1) - 9._dp*TrCYtYt*mue*AL1*Conjg(L1)/2._dp -              &
  &  18._dp*TrYdadjYd*mue*AL1*Conjg(L1) - 6._dp*TrYeadjYe*mue*AL1*Conjg(L1) -              &
  &  3._dp*TrYtCYt*mue*AL1*Conjg(L1)/2._dp + 18._dp*g12*L1*Amue*Conjg(L1)/5._dp +        &
  &  12._dp*g22*L1*Amue*Conjg(L1) - 9._dp*L1*TrCYtYt*Amue*Conjg(L1)/4._dp -              &
  &  9._dp*L1*TrYdadjYd*Amue*Conjg(L1) - 3._dp*L1*TrYeadjYe*Amue*Conjg(L1) -               &
  &  3._dp*L1*TrYtCYt*Amue*Conjg(L1)/4._dp - 48._dp*L1*mue*AL1*Conjg(L1)**2 -              &
  &  12._dp*L1**2*Amue*Conjg(L1)**2 - 36._dp*g12*L2*MassB*mue*Conjg(L2)/5._dp -          &
  &  24._dp*g22*L2*MassWB*mue*Conjg(L2) - 18._dp*L2*TrAYuadjYu*mue*Conjg(L2)             &
  &  + 36._dp*g12*mue*AL2*Conjg(L2)/5._dp + 24._dp*g22*mue*AL2*Conjg(L2)               &
  &  - 18._dp*TrYuadjYu*mue*AL2*Conjg(L2) + 18._dp*g12*L2*Amue*Conjg(L2)/5._dp +         &
  &  12._dp*g22*L2*Amue*Conjg(L2) - 9._dp*L2*TrYuadjYu*Amue*Conjg(L2) - 48._dp*L2*mue*AL2*Conjg(L2)&
  & **2 - 12._dp*L2**2*Amue*Conjg(L2)**2


  DAmue = oo16pi2*( betaAmue1 + oo16pi2 * betaAmue2 )


  Else
  DAmue = oo16pi2* betaAmue1
  End If


  !--------------------
  ! AMTM
  !--------------------

  betaAMTM1 = 24._dp*g12*MassB*MTM/5._dp + 16._dp*g22*MassWB*MTM + 2._dp*MTM*TrCYtAYt -&
  &  12._dp*g12*AMTM/5._dp - 8._dp*g22*AMTM + 3._dp*TrCYtYt*AMTM/4._dp +               &
  &  (TrYtCYt*AMTM)/4._dp + 2._dp*MTM*AL1*Conjg(L1) + L1*AMTM*Conjg(L1) + 2._dp*MTM*AL2*Conjg(L2)&
  &  + L2*AMTM*Conjg(L2)


  If (TwoLoopRGE) Then
  betaAMTM2 = -3552._dp*g14*MassB*MTM/25._dp - 192._dp*g12*g22*MassB*MTM/5._dp -  &
  &  192._dp*g12*g22*MassWB*MTM/5._dp - 384._dp*g24*MassWB*MTM - 6._dp*g12*MTM*TrCYtAYt/5._dp -&
  &  2._dp*g22*MTM*TrCYtAYt - 9._dp*MTM*TrCYtAYtCYtYt + 9._dp*g12*MassB*MTM*TrCYtYt/10._dp +&
  &  3._dp*g22*MassWB*MTM*TrCYtYt/2._dp - 4._dp*MTM*TrCYtYtadjYeAYe - 12._dp*MTM*TrCYtYtadjYzAYz -&
  &  12._dp*MTM*TrCYtYtCYtAYt - 4._dp*MTM*TrYeCYtAYtadjYe + 3._dp*g12*MassB*MTM*TrYtCYt/10._dp +&
  &  (g22*MassWB*MTM*TrYtCYt)/2._dp - 3._dp*MTM*TrYtCYtAYtCYt - 12._dp*MTM*TrYzCYtAYtadjYz +&
  &  888._dp*g14*AMTM/25._dp + 96._dp*g12*g22*AMTM/5._dp + 96._dp*g24*AMTM -       &
  &  9._dp*g12*TrCYtYt*AMTM/20._dp - 3._dp*g22*TrCYtYt*AMTM/4._dp - 9._dp*TrCYtYtCYtYt*AMTM/2._dp -&
  &  2._dp*TrYeCYtYtadjYe*AMTM - 3._dp*g12*TrYtCYt*AMTM/20._dp - (g22*TrYtCYt*AMTM)    &
  & /4._dp - 3._dp*TrYtCYtYtCYt*AMTM/2._dp - 6._dp*TrYzCYtYtadjYz*AMTM + 6._dp*g12*L1*MassB*MTM*Conjg(L1)&
  & /5._dp + 2._dp*g22*L1*MassWB*MTM*Conjg(L1) - 12._dp*L1*MTM*TrAYdadjYd*Conjg(L1)      &
  &  - 4._dp*L1*MTM*TrAYeadjYe*Conjg(L1) - 6._dp*g12*MTM*AL1*Conjg(L1)/5._dp -           &
  &  2._dp*g22*MTM*AL1*Conjg(L1) - 12._dp*MTM*TrYdadjYd*AL1*Conjg(L1) - 4._dp*MTM*TrYeadjYe*AL1*Conjg(L1)&
  &  - 3._dp*g12*L1*AMTM*Conjg(L1)/5._dp - g22*L1*AMTM*Conjg(L1) - 6._dp*L1*TrYdadjYd*AMTM*Conjg(L1)&
  &  - 2._dp*L1*TrYeadjYe*AMTM*Conjg(L1) - 24._dp*L1*MTM*AL1*Conjg(L1)**2 - 6._dp*L1**2*AMTM*Conjg(L1)&
  & **2 + 6._dp*g12*L2*MassB*MTM*Conjg(L2)/5._dp + 2._dp*g22*L2*MassWB*MTM*Conjg(L2)   &
  &  - 12._dp*L2*MTM*TrAYuadjYu*Conjg(L2) - 6._dp*g12*MTM*AL2*Conjg(L2)/5._dp -          &
  &  2._dp*g22*MTM*AL2*Conjg(L2) - 12._dp*MTM*TrYuadjYu*AL2*Conjg(L2) - 3._dp*g12*L2*AMTM*Conjg(L2)&
  & /5._dp - g22*L2*AMTM*Conjg(L2) - 6._dp*L2*TrYuadjYu*AMTM*Conjg(L2) - 24._dp*L2*MTM*AL2*Conjg(L2)&
  & **2 - 6._dp*L2**2*AMTM*Conjg(L2)**2


  DAMTM = oo16pi2*( betaAMTM1 + oo16pi2 * betaAMTM2 )


  Else
  DAMTM = oo16pi2* betaAMTM1
  End If


  !--------------------
  ! AMZM
  !--------------------

  betaAMZM1 = 2._dp*g12*MassB*MZM/15._dp + 32._dp*g32*MassG*MZM/3._dp +             &
  &  6._dp*g22*MassWB*MZM + 2._dp*MZM*TrAYzadjYz - (g12*AMZM)/15._dp - 3._dp*g22*AMZM -&
  &  16._dp*g32*AMZM/3._dp + TrYzadjYz*AMZM


  If (TwoLoopRGE) Then
  betaAMZM2 = -818._dp*g14*MassB*MZM/225._dp - 2._dp*g12*g22*MassB*MZM/5._dp -    &
  &  32._dp*g12*g32*MassB*MZM/45._dp - 32._dp*g12*g32*MassG*MZM/45._dp -           &
  &  32._dp*g22*g32*MassG*MZM - 1280._dp*g34*MassG*MZM/9._dp - 2._dp*g12*g22*MassWB*MZM/5._dp -&
  &  114._dp*g24*MassWB*MZM - 32._dp*g22*g32*MassWB*MZM + 4._dp*g12*MZM*TrAYzadjYz/5._dp -&
  &  4._dp*MZM*TrCYsAYzadjYzYs - 5._dp*MZM*TrCYsYzadjYzAYs - 3._dp*MZM*TrCYtAYtadjYzYz -   &
  &  4._dp*MZM*TrCYtYtadjYzAYz - 4._dp*MZM*TrYdadjYdAYzadjYz - 2._dp*MZM*TrYeadjYzAYzadjYe -&
  &  4._dp*MZM*TrYsCYsAYzadjYz - 2._dp*MZM*TrYtadjYzAYzCYt - 2._dp*MZM*TrYzadjYeAYeadjYz - &
  &  4._dp*g12*MassB*MZM*TrYzadjYz/5._dp - 4._dp*MZM*TrYzadjYzAYdadjYd - 3._dp*MZM*TrYzadjYzAYsCYs -&
  &  20._dp*MZM*TrYzadjYzAYzadjYz - 3._dp*MZM*TrYzCYtAYtadjYz + 409._dp*g14*AMZM/450._dp +&
  &  (g12*g22*AMZM)/5._dp + 57._dp*g24*AMZM/2._dp + 16._dp*g12*g32*AMZM/45._dp + &
  &  16._dp*g22*g32*AMZM + 320._dp*g34*AMZM/9._dp - 2._dp*TrCYsYzadjYzYs*AMZM -      &
  &  3._dp*TrCYtYtadjYzYz*AMZM/2._dp - TrYeadjYzYzadjYe*AMZM - (TrYsCYsYzadjYz*AMZM)       &
  & /2._dp - TrYtadjYzYzCYt*AMZM + 2._dp*g12*TrYzadjYz*AMZM/5._dp - 2._dp*TrYzadjYzYdadjYd*AMZM -&
  &  3._dp*TrYzadjYzYsCYs*AMZM/2._dp - 5._dp*TrYzadjYzYzadjYz*AMZM - (TrYzCYtYtadjYz*AMZM)/2._dp


  DAMZM = oo16pi2*( betaAMZM1 + oo16pi2 * betaAMZM2 )


  Else
  DAMZM = oo16pi2* betaAMZM1
  End If


  !--------------------
  ! AMSM
  !--------------------

  betaAMSM1 = 32._dp*g12*MassB*MSM/15._dp + 80._dp*g32*MassG*MSM/3._dp +            &
  &  2._dp*MSM*TrCYsAYs - 16._dp*g12*AMSM/15._dp - 40._dp*g32*AMSM/3._dp +             &
  &  3._dp*TrCYsYs*AMSM/4._dp + (TrYsCYs*AMSM)/4._dp


  If (TwoLoopRGE) Then
  betaAMSM2 = -13568._dp*g14*MassB*MSM/225._dp - 256._dp*g12*g32*MassB*MSM/9._dp -&
  &  256._dp*g12*g32*MassG*MSM/9._dp - 5120._dp*g34*MassG*MSM/9._dp - 8._dp*g12*MSM*TrCYsAYs/15._dp -&
  &  8._dp*g32*MSM*TrCYsAYs/3._dp - 12._dp*MSM*TrCYsAYsCYsYs - 6._dp*MSM*TrCYsYdadjYdAYs +&
  &  2._dp*g12*MassB*MSM*TrCYsYs/5._dp + 2._dp*g32*MassG*MSM*TrCYsYs - 16._dp*MSM*TrCYsYsCYsAYs -&
  &  6._dp*MSM*TrCYsYzadjYzAYs - 2._dp*MSM*TrYdadjYdAYsCYs + 2._dp*g12*MassB*MSM*TrYsCYs/15._dp +&
  &  2._dp*g32*MassG*MSM*TrYsCYs/3._dp - 8._dp*MSM*TrYsCYsAYdadjYd - 4._dp*MSM*TrYsCYsAYsCYs -&
  &  8._dp*MSM*TrYsCYsAYzadjYz - 2._dp*MSM*TrYzadjYzAYsCYs + 3392._dp*g14*AMSM/225._dp + &
  &  128._dp*g12*g32*AMSM/9._dp + 1280._dp*g34*AMSM/9._dp - (g12*TrCYsYs*AMSM)     &
  & /5._dp - g32*TrCYsYs*AMSM - 6._dp*TrCYsYsCYsYs*AMSM - (g12*TrYsCYs*AMSM)           &
  & /15._dp - (g32*TrYsCYs*AMSM)/3._dp - 4._dp*TrYsCYsYdadjYd*AMSM - 2._dp*TrYsCYsYsCYs*AMSM -&
  &  4._dp*TrYsCYsYzadjYz*AMSM


  DAMSM = oo16pi2*( betaAMSM1 + oo16pi2 * betaAMSM2 )


  Else
  DAMSM = oo16pi2* betaAMSM1
  End If


  !--------------------
  ! mq2
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betamq21(i1,i2) = mq2TYdCYd(i1,i2) + mq2TYuCYu(i1,i2) + 2._dp*TAYdCAYd(i1,i2)         &
  &  + 2._dp*TAYuCAYu(i1,i2) + 2._dp*mHd2*TYdCYd(i1,i2) + TYdCYdmq2(i1,i2) +               &
  &  2._dp*TYdmd2CYd(i1,i2) + 2._dp*mHu2*TYuCYu(i1,i2) + TYuCYumq2(i1,i2) + 2._dp*TYumu2CYu(i1,i2)

  If (i1.eq.i2) Then
  betamq21(i1,i2) = betamq21(i1,i2)+(-2._dp*g12*MassB*Conjg(MassB)/15._dp - 32._dp*g32*MassG*Conjg(MassG)&
  & /3._dp - 6._dp*g22*MassWB*Conjg(MassWB) + (g12*Tr1(1))/3._dp)
  End If
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betamq22(i1,i2) = 2._dp*g12*mq2TYdCYd(i1,i2)/5._dp - 3._dp*TrYdadjYd*mq2TYdCYd(i1,i2)&
  &  - TrYeadjYe*mq2TYdCYd(i1,i2) - 3._dp*AbsL1sq*mq2TYdCYd(i1,i2) - 2._dp*mq2TYdCYdTYdCYd(i1,i2)&
  &  - 4._dp*mq2TYdCYsYsCYd(i1,i2) - 2._dp*mq2TYdCYzTYzCYd(i1,i2) + 4._dp*g12*mq2TYuCYu(i1,i2)&
  & /5._dp - 3._dp*TrYuadjYu*mq2TYuCYu(i1,i2) - 3._dp*AbsL2sq*mq2TYuCYu(i1,i2)        &
  &  - 2._dp*mq2TYuCYuTYuCYu(i1,i2) + 4._dp*g12*TAYdCAYd(i1,i2)/5._dp - 6._dp*TrYdadjYd*TAYdCAYd(i1,i2)&
  &  - 2._dp*TrYeadjYe*TAYdCAYd(i1,i2) - 6._dp*AbsL1sq*TAYdCAYd(i1,i2) -              &
  &  4._dp*TAYdCAYdTYdCYd(i1,i2) - 8._dp*TAYdCAYsYsCYd(i1,i2) - 4._dp*TAYdCAYzTYzCYd(i1,i2)&
  &  - 6._dp*TrYdadjAYd*TAYdCYd(i1,i2) - 2._dp*TrYeadjAYe*TAYdCYd(i1,i2) - 4._dp*g12*Conjg(MassB)&
  & *TAYdCYd(i1,i2)/5._dp - 6._dp*L1*Conjg(AL1)*TAYdCYd(i1,i2) - 4._dp*TAYdCYdTYdCAYd(i1,i2)&
  &  - 8._dp*TAYdCYsYsCAYd(i1,i2) - 4._dp*TAYdCYzTYzCAYd(i1,i2) + 8._dp*g12*TAYuCAYu(i1,i2)&
  & /5._dp - 6._dp*TrYuadjYu*TAYuCAYu(i1,i2) - 6._dp*AbsL2sq*TAYuCAYu(i1,i2)          &
  &  - 4._dp*TAYuCAYuTYuCYu(i1,i2) - 6._dp*TrYuadjAYu*TAYuCYu(i1,i2) - 8._dp*g12*Conjg(MassB)&
  & *TAYuCYu(i1,i2)/5._dp - 6._dp*L2*Conjg(AL2)*TAYuCYu(i1,i2) - 4._dp*TAYuCYuTYuCAYu(i1,i2)&
  &  - 4._dp*g12*MassB*TYdCAYd(i1,i2)/5._dp - 6._dp*TrAYdadjYd*TYdCAYd(i1,i2)            &
  &  - 2._dp*TrAYeadjYe*TYdCAYd(i1,i2) - 6._dp*AL1*Conjg(L1)*TYdCAYd(i1,i2) -              &
  &  4._dp*TYdCAYdTAYdCYd(i1,i2) - 8._dp*TYdCAYsAYsCYd(i1,i2) - 4._dp*TYdCAYzTAYzCYd(i1,i2)&
  &  + 4._dp*g12*mHd2*TYdCYd(i1,i2)/5._dp - 6._dp*TrAYdadjAYd*TYdCYd(i1,i2)              &
  &  - 2._dp*TrAYeadjAYe*TYdCYd(i1,i2) - 6._dp*Trmd2YdadjYd*TYdCYd(i1,i2) - 2._dp*Trme2YeadjYe*TYdCYd(i1,i2)&
  &  - 12._dp*mHd2*TrYdadjYd*TYdCYd(i1,i2) - 6._dp*TrYdmq2adjYd*TYdCYd(i1,i2)              &
  &  - 4._dp*mHd2*TrYeadjYe*TYdCYd(i1,i2) - 2._dp*TrYeml2adjYe*TYdCYd(i1,i2)               &
  &  - 18._dp*L1*mHd2*Conjg(L1)*TYdCYd(i1,i2) - 6._dp*L1*mt2*Conjg(L1)*TYdCYd(i1,i2)       &
  &  + 8._dp*g12*MassB*Conjg(MassB)*TYdCYd(i1,i2)/5._dp - 6._dp*AL1*Conjg(AL1)           &
  & *TYdCYd(i1,i2) + 2._dp*g12*TYdCYdmq2(i1,i2)/5._dp - 3._dp*TrYdadjYd*TYdCYdmq2(i1,i2) &
  &  - TrYeadjYe*TYdCYdmq2(i1,i2) - 3._dp*AbsL1sq*TYdCYdmq2(i1,i2) - 4._dp*TYdCYdmq2TYdCYd(i1,i2)&
  &  - 4._dp*TYdCYdTAYdCAYd(i1,i2) - 8._dp*mHd2*TYdCYdTYdCYd(i1,i2) - 2._dp*TYdCYdTYdCYdmq2(i1,i2)&
  &  - 4._dp*TYdCYdTYdmd2CYd(i1,i2) - 8._dp*TYdCYsAYsCAYd(i1,i2) - 8._dp*TYdCYsmd2YsCYd(i1,i2)&
  &  - 8._dp*mHd2*TYdCYsYsCYd(i1,i2) - 8._dp*ms2*TYdCYsYsCYd(i1,i2) - 4._dp*TYdCYsYsCYdmq2(i1,i2)&
  &  - 8._dp*TYdCYsYsmd2CYd(i1,i2) - 4._dp*TYdCYzml2TYzCYd(i1,i2) - 4._dp*TYdCYzTAYzCAYd(i1,i2)&
  &  - 4._dp*mHd2*TYdCYzTYzCYd(i1,i2) - 4._dp*mz2*TYdCYzTYzCYd(i1,i2) - 2._dp*TYdCYzTYzCYdmq2(i1,i2)&
  &  - 4._dp*TYdCYzTYzmd2CYd(i1,i2) + 4._dp*g12*TYdmd2CYd(i1,i2)/5._dp - 6._dp*TrYdadjYd*TYdmd2CYd(i1,i2)&
  &  - 2._dp*TrYeadjYe*TYdmd2CYd(i1,i2) - 6._dp*AbsL1sq*TYdmd2CYd(i1,i2)              &
  &  - 4._dp*TYdmd2CYdTYdCYd(i1,i2) - 8._dp*TYdmd2CYsYsCYd(i1,i2) - 4._dp*TYdmd2CYzTYzCYd(i1,i2)&
  &  - 8._dp*g12*MassB*TYuCAYu(i1,i2)/5._dp - 6._dp*TrAYuadjYu*TYuCAYu(i1,i2)            &
  &  - 6._dp*AL2*Conjg(L2)*TYuCAYu(i1,i2) - 4._dp*TYuCAYuTAYuCYu(i1,i2) + 8._dp*g12*mHu2*TYuCYu(i1,i2)&
  & /5._dp - 6._dp*TrAYuadjAYu*TYuCYu(i1,i2) - 6._dp*Trmu2YuadjYu*TYuCYu(i1,i2)            &
  &  - 12._dp*mHu2*TrYuadjYu*TYuCYu(i1,i2) - 6._dp*TrYumq2adjYu*TYuCYu(i1,i2)              &
  &  - 18._dp*L2*mHu2*Conjg(L2)*TYuCYu(i1,i2) - 6._dp*L2*mtb2*Conjg(L2)*TYuCYu(i1,i2)      &
  &  + 16._dp*g12*MassB*Conjg(MassB)*TYuCYu(i1,i2)/5._dp - 6._dp*AL2*Conjg(AL2)          &
  & *TYuCYu(i1,i2) + 4._dp*g12*TYuCYumq2(i1,i2)/5._dp - 3._dp*TrYuadjYu*TYuCYumq2(i1,i2) &
  &  - 3._dp*AbsL2sq*TYuCYumq2(i1,i2) - 4._dp*TYuCYumq2TYuCYu(i1,i2) - 4._dp*TYuCYuTAYuCAYu(i1,i2)&
  &  - 8._dp*mHu2*TYuCYuTYuCYu(i1,i2) - 2._dp*TYuCYuTYuCYumq2(i1,i2) - 4._dp*TYuCYuTYumu2CYu(i1,i2)&
  &  + 8._dp*g12*TYumu2CYu(i1,i2)/5._dp - 6._dp*TrYuadjYu*TYumu2CYu(i1,i2)               &
  &  - 6._dp*AbsL2sq*TYumu2CYu(i1,i2) - 4._dp*TYumu2CYuTYuCYu(i1,i2)

  If (i1.eq.i2) Then
  betamq22(i1,i2) = betamq22(i1,i2)+(409._dp*g14*MassB*Conjg(MassB)/75._dp + 2._dp*g12*g22*MassB*Conjg(MassB)&
  & /5._dp + 32._dp*g12*g32*MassB*Conjg(MassB)/45._dp + 16._dp*g12*g32*MassG*Conjg(MassB)&
  & /45._dp + (g12*g22*MassWB*Conjg(MassB))/5._dp + 16._dp*g12*g32*MassB*Conjg(MassG)&
  & /45._dp + 32._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 32._dp*g22*g32*MassG*Conjg(MassG)&
  &  + 544._dp*g34*MassG*Conjg(MassG)/3._dp + 16._dp*g22*g32*MassWB*Conjg(MassG)     &
  &  + (g12*g22*MassB*Conjg(MassWB))/5._dp + 16._dp*g22*g32*MassG*Conjg(MassWB)    &
  &  + 2._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g24*MassWB*Conjg(MassWB)   &
  &  + 32._dp*g22*g32*MassWB*Conjg(MassWB) + 2._dp*g14*Tr2(1)/15._dp +               &
  &  6._dp*g24*Tr2(2) + 32._dp*g34*Tr2(3)/3._dp + 4._dp*g12*Tr3(1)/3._dp +           &
  &  4._dp*g32*Tr3(3) + 4._dp*(g32*Tr3(3))/sqrt3)
  End If
   End Do
  End Do


  Dmq2 = oo16pi2*( betamq21 + oo16pi2 * betamq22 )


  Else
  Dmq2 = oo16pi2* betamq21
  End If


  !--------------------
  ! ml2
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betaml21(i1,i2) = 6._dp*AYtCAYt(i1,i2) + ml2TYeCYe(i1,i2) + 3._dp*ml2TYzCYz(i1,i2)    &
  &  + 3._dp*ml2YtCYt(i1,i2) + 2._dp*TAYeCAYe(i1,i2) + 6._dp*TAYzCAYz(i1,i2)               &
  &  + 2._dp*mHd2*TYeCYe(i1,i2) + TYeCYeml2(i1,i2) + 2._dp*TYeme2CYe(i1,i2) +              &
  &  6._dp*mz2*TYzCYz(i1,i2) + 3._dp*TYzCYzml2(i1,i2) + 6._dp*TYzmd2CYz(i1,i2)             &
  &  + 6._dp*mt2*YtCYt(i1,i2) + 3._dp*YtCYtml2(i1,i2) + 6._dp*Ytml2CYt(i1,i2)

  If (i1.eq.i2) Then
  betaml21(i1,i2) = betaml21(i1,i2)+(-6._dp*g12*MassB*Conjg(MassB)/5._dp - 6._dp*g22*MassWB*Conjg(MassWB)&
  &  - g12*Tr1(1))
  End If
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betaml22(i1,i2) = -6._dp*AYtadjAYeYeCYt(i1,i2) - 18._dp*AYtadjAYzYzCYt(i1,i2)         &
  &  - 6._dp*AYtadjYeYeCAYt(i1,i2) - 18._dp*AYtadjYzYzCAYt(i1,i2) + 36._dp*g12*AYtCAYt(i1,i2)&
  & /5._dp + 24._dp*g22*AYtCAYt(i1,i2) - 9._dp*TrCYtYt*AYtCAYt(i1,i2)/2._dp -            &
  &  3._dp*TrYtCYt*AYtCAYt(i1,i2)/2._dp - 18._dp*AYtCAYtYtCYt(i1,i2) - 6._dp*TrYtCAYt*AYtCYt(i1,i2)&
  &  - 18._dp*AYtCYtYtCAYt(i1,i2) - 6._dp*L1*AYtCAYt(i1,i2)*Conjg(L1) - 36._dp*g12*AYtCYt(i1,i2)&
  & *Conjg(MassB)/5._dp - 24._dp*g22*AYtCYt(i1,i2)*Conjg(MassWB) - 6._dp*L1*AYtCYt(i1,i2)&
  & *Conjg(AL1) + 6._dp*g12*ml2TYeCYe(i1,i2)/5._dp - 3._dp*TrYdadjYd*ml2TYeCYe(i1,i2)    &
  &  - TrYeadjYe*ml2TYeCYe(i1,i2) - 3._dp*AbsL1sq*ml2TYeCYe(i1,i2) - 2._dp*ml2TYeCYeTYeCYe(i1,i2)&
  &  - 6._dp*ml2TYzCYdTYdCYz(i1,i2) - 12._dp*ml2TYzCYsYsCYz(i1,i2) - 2._dp*g12*ml2TYzCYz(i1,i2)&
  & /5._dp + 16._dp*g32*ml2TYzCYz(i1,i2) - 3._dp*TrYzadjYz*ml2TYzCYz(i1,i2)              &
  &  - 6._dp*ml2TYzCYzTYzCYz(i1,i2) - 3._dp*ml2YtadjYeYeCYt(i1,i2) - 9._dp*ml2YtadjYzYzCYt(i1,i2)&
  &  + 18._dp*g12*ml2YtCYt(i1,i2)/5._dp + 12._dp*g22*ml2YtCYt(i1,i2) - 9._dp*TrCYtYt*ml2YtCYt(i1,i2)&
  & /4._dp - 3._dp*TrYtCYt*ml2YtCYt(i1,i2)/4._dp - 3._dp*AbsL1sq*ml2YtCYt(i1,i2)      &
  &  - 9._dp*ml2YtCYtYtCYt(i1,i2) + 12._dp*g12*TAYeCAYe(i1,i2)/5._dp - 6._dp*TrYdadjYd*TAYeCAYe(i1,i2)&
  &  - 2._dp*TrYeadjYe*TAYeCAYe(i1,i2) - 6._dp*AbsL1sq*TAYeCAYe(i1,i2) -              &
  &  4._dp*TAYeCAYeTYeCYe(i1,i2) - 6._dp*TrYdadjAYd*TAYeCYe(i1,i2) - 2._dp*TrYeadjAYe*TAYeCYe(i1,i2)&
  &  - 12._dp*g12*Conjg(MassB)*TAYeCYe(i1,i2)/5._dp - 6._dp*L1*Conjg(AL1)*TAYeCYe(i1,i2) &
  &  - 4._dp*TAYeCYeTYeCAYe(i1,i2) - 12._dp*TAYzCAYdTYdCYz(i1,i2) - 24._dp*TAYzCAYsYsCYz(i1,i2)&
  &  - 4._dp*g12*TAYzCAYz(i1,i2)/5._dp + 32._dp*g32*TAYzCAYz(i1,i2) - 6._dp*TrYzadjYz*TAYzCAYz(i1,i2)&
  &  - 12._dp*TAYzCAYzTYzCYz(i1,i2) - 12._dp*TAYzCYdTYdCAYz(i1,i2) - 24._dp*TAYzCYsYsCAYz(i1,i2)&
  &  - 6._dp*TrYzadjAYz*TAYzCYz(i1,i2) + 4._dp*g12*Conjg(MassB)*TAYzCYz(i1,i2)           &
  & /5._dp - 32._dp*g32*Conjg(MassG)*TAYzCYz(i1,i2) - 12._dp*TAYzCYzTYzCAYz(i1,i2)       &
  &  - 12._dp*g12*MassB*TYeCAYe(i1,i2)/5._dp - 6._dp*TrAYdadjYd*TYeCAYe(i1,i2)           &
  &  - 2._dp*TrAYeadjYe*TYeCAYe(i1,i2) - 6._dp*AL1*Conjg(L1)*TYeCAYe(i1,i2) -              &
  &  4._dp*TYeCAYeTAYeCYe(i1,i2) + 12._dp*g12*mHd2*TYeCYe(i1,i2)/5._dp - 6._dp*TrAYdadjAYd*TYeCYe(i1,i2)&
  &  - 2._dp*TrAYeadjAYe*TYeCYe(i1,i2) - 6._dp*Trmd2YdadjYd*TYeCYe(i1,i2) - 2._dp*Trme2YeadjYe*TYeCYe(i1,i2)&
  &  - 12._dp*mHd2*TrYdadjYd*TYeCYe(i1,i2) - 6._dp*TrYdmq2adjYd*TYeCYe(i1,i2)              &
  &  - 4._dp*mHd2*TrYeadjYe*TYeCYe(i1,i2) - 2._dp*TrYeml2adjYe*TYeCYe(i1,i2)               &
  &  - 18._dp*L1*mHd2*Conjg(L1)*TYeCYe(i1,i2) - 6._dp*L1*mt2*Conjg(L1)*TYeCYe(i1,i2)       &
  &  + 24._dp*g12*MassB*Conjg(MassB)*TYeCYe(i1,i2)/5._dp - 6._dp*AL1*Conjg(AL1)          &
  & *TYeCYe(i1,i2) + 6._dp*g12*TYeCYeml2(i1,i2)/5._dp - 3._dp*TrYdadjYd*TYeCYeml2(i1,i2) &
  &  - TrYeadjYe*TYeCYeml2(i1,i2) - 3._dp*AbsL1sq*TYeCYeml2(i1,i2) - 4._dp*TYeCYeml2TYeCYe(i1,i2)&
  &  - 4._dp*TYeCYeTAYeCAYe(i1,i2) - 8._dp*mHd2*TYeCYeTYeCYe(i1,i2) - 2._dp*TYeCYeTYeCYeml2(i1,i2)&
  &  - 4._dp*TYeCYeTYeme2CYe(i1,i2) + 12._dp*g12*TYeme2CYe(i1,i2)/5._dp - 6._dp*TrYdadjYd*TYeme2CYe(i1,i2)&
  &  - 2._dp*TrYeadjYe*TYeme2CYe(i1,i2) - 6._dp*AbsL1sq*TYeme2CYe(i1,i2)              &
  &  - 4._dp*TYeme2CYeTYeCYe(i1,i2) - 12._dp*TYzCAYdTAYdCYz(i1,i2) - 24._dp*TYzCAYsAYsCYz(i1,i2)&
  &  + 4._dp*g12*MassB*TYzCAYz(i1,i2)/5._dp - 32._dp*g32*MassG*TYzCAYz(i1,i2)          &
  &  - 6._dp*TrAYzadjYz*TYzCAYz(i1,i2) - 12._dp*TYzCAYzTAYzCYz(i1,i2) - 12._dp*TYzCYdmq2TYdCYz(i1,i2)&
  &  - 12._dp*TYzCYdTAYdCAYz(i1,i2) - 12._dp*mHd2*TYzCYdTYdCYz(i1,i2) - 12._dp*mz2*TYzCYdTYdCYz(i1,i2)&
  &  - 6._dp*TYzCYdTYdCYzml2(i1,i2) - 12._dp*TYzCYdTYdmd2CYz(i1,i2) - 24._dp*TYzCYsAYsCAYz(i1,i2)&
  &  - 24._dp*TYzCYsmd2YsCYz(i1,i2) - 24._dp*ms2*TYzCYsYsCYz(i1,i2) - 24._dp*mz2*TYzCYsYsCYz(i1,i2)&
  &  - 12._dp*TYzCYsYsCYzml2(i1,i2) - 24._dp*TYzCYsYsmd2CYz(i1,i2) - 4._dp*g12*mz2*TYzCYz(i1,i2)&
  & /5._dp + 32._dp*g32*mz2*TYzCYz(i1,i2) - 6._dp*TrAYzadjAYz*TYzCYz(i1,i2)              &
  &  - 6._dp*Trmd2YzadjYz*TYzCYz(i1,i2) - 12._dp*mz2*TrYzadjYz*TYzCYz(i1,i2)               &
  &  - 6._dp*TrYzml2adjYz*TYzCYz(i1,i2) - 8._dp*g12*MassB*Conjg(MassB)*TYzCYz(i1,i2)     &
  & /5._dp + 64._dp*g32*MassG*Conjg(MassG)*TYzCYz(i1,i2) - 2._dp*g12*TYzCYzml2(i1,i2)  &
  & /5._dp + 16._dp*g32*TYzCYzml2(i1,i2) - 3._dp*TrYzadjYz*TYzCYzml2(i1,i2)              &
  &  - 12._dp*TYzCYzml2TYzCYz(i1,i2) - 12._dp*TYzCYzTAYzCAYz(i1,i2) - 24._dp*mz2*TYzCYzTYzCYz(i1,i2)&
  &  - 6._dp*TYzCYzTYzCYzml2(i1,i2) - 12._dp*TYzCYzTYzmd2CYz(i1,i2) - 12._dp*TYzmd2CYdTYdCYz(i1,i2)&
  &  - 24._dp*TYzmd2CYsYsCYz(i1,i2) - 4._dp*g12*TYzmd2CYz(i1,i2)/5._dp + 32._dp*g32*TYzmd2CYz(i1,i2)&
  &  - 6._dp*TrYzadjYz*TYzmd2CYz(i1,i2) - 12._dp*TYzmd2CYzTYzCYz(i1,i2) - 6._dp*YtadjAYeAYeCYt(i1,i2)&
  &  - 18._dp*YtadjAYzAYzCYt(i1,i2) - 6._dp*YtadjYeAYeCAYt(i1,i2) - 6._dp*YtadjYeme2YeCYt(i1,i2)&
  &  - 6._dp*mHd2*YtadjYeYeCYt(i1,i2) - 6._dp*mt2*YtadjYeYeCYt(i1,i2) - 3._dp*YtadjYeYeCYtml2(i1,i2)&
  &  - 6._dp*YtadjYeYeml2CYt(i1,i2) - 18._dp*YtadjYzAYzCAYt(i1,i2) - 18._dp*YtadjYzmd2YzCYt(i1,i2)&
  &  - 18._dp*mt2*YtadjYzYzCYt(i1,i2) - 18._dp*mz2*YtadjYzYzCYt(i1,i2) - 9._dp*YtadjYzYzCYtml2(i1,i2)&
  &  - 18._dp*YtadjYzYzml2CYt(i1,i2) - 36._dp*g12*MassB*YtCAYt(i1,i2)/5._dp -            &
  &  24._dp*g22*MassWB*YtCAYt(i1,i2) - 6._dp*TrCYtAYt*YtCAYt(i1,i2) - 6._dp*AL1*Conjg(L1)&
  & *YtCAYt(i1,i2) - 18._dp*YtCAYtAYtCYt(i1,i2) + 36._dp*g12*mt2*YtCYt(i1,i2)            &
  & /5._dp + 24._dp*g22*mt2*YtCYt(i1,i2) - 3._dp*TrAYtCAYt*YtCYt(i1,i2)/2._dp -          &
  &  9._dp*TrCAYtAYt*YtCYt(i1,i2)/2._dp - 3._dp*TrCYtml2Yt*YtCYt(i1,i2) - 9._dp*mt2*TrCYtYt*YtCYt(i1,i2)&
  &  - 6._dp*TrCYtYtml2*YtCYt(i1,i2) - 3._dp*mt2*TrYtCYt*YtCYt(i1,i2) - 3._dp*TrYtCYtml2*YtCYt(i1,i2)&
  &  - 12._dp*L1*mHd2*Conjg(L1)*YtCYt(i1,i2) - 12._dp*L1*mt2*Conjg(L1)*YtCYt(i1,i2)        &
  &  + 72._dp*g12*MassB*Conjg(MassB)*YtCYt(i1,i2)/5._dp + 48._dp*g22*MassWB*Conjg(MassWB)&
  & *YtCYt(i1,i2) - 6._dp*AL1*Conjg(AL1)*YtCYt(i1,i2) - 18._dp*YtCYtAYtCAYt(i1,i2)         &
  &  + 18._dp*g12*YtCYtml2(i1,i2)/5._dp + 12._dp*g22*YtCYtml2(i1,i2) - 9._dp*TrCYtYt*YtCYtml2(i1,i2)&
  & /4._dp - 3._dp*TrYtCYt*YtCYtml2(i1,i2)/4._dp - 3._dp*AbsL1sq*YtCYtml2(i1,i2)      &
  &  - 18._dp*YtCYtml2YtCYt(i1,i2) - 36._dp*mt2*YtCYtYtCYt(i1,i2) - 9._dp*YtCYtYtCYtml2(i1,i2)&
  &  - 18._dp*YtCYtYtml2CYt(i1,i2) - 6._dp*Ytml2adjYeYeCYt(i1,i2) - 18._dp*Ytml2adjYzYzCYt(i1,i2)&
  &  + 36._dp*g12*Ytml2CYt(i1,i2)/5._dp + 24._dp*g22*Ytml2CYt(i1,i2) - 9._dp*TrCYtYt*Ytml2CYt(i1,i2)&
  & /2._dp - 3._dp*TrYtCYt*Ytml2CYt(i1,i2)/2._dp - 6._dp*AbsL1sq*Ytml2CYt(i1,i2)      &
  &  - 18._dp*Ytml2CYtYtCYt(i1,i2)

  If (i1.eq.i2) Then
  betaml22(i1,i2) = betaml22(i1,i2)+(1251._dp*g14*MassB*Conjg(MassB)/25._dp + 18._dp*g12*g22*MassB*Conjg(MassB)&
  & /5._dp + 9._dp*g12*g22*MassWB*Conjg(MassB)/5._dp + 9._dp*g12*g22*MassB*Conjg(MassWB)&
  & /5._dp + 18._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g24*MassWB*Conjg(MassWB)&
  &  + 6._dp*g14*Tr2(1)/5._dp + 6._dp*g24*Tr2(2) - 4._dp*g12*Tr3(1))
  End If
   End Do
  End Do


  Dml2 = oo16pi2*( betaml21 + oo16pi2 * betaml22 )


  Else
  Dml2 = oo16pi2* betaml21
  End If


  !--------------------
  ! mHd2
  !--------------------

  betamHd21 = 6._dp*TrAYdadjAYd + 2._dp*TrAYeadjAYe + 6._dp*Trmd2YdadjYd +              &
  &  2._dp*Trme2YeadjYe + 6._dp*mHd2*TrYdadjYd + 6._dp*TrYdmq2adjYd + 2._dp*mHd2*TrYeadjYe +&
  &  2._dp*TrYeml2adjYe + 12._dp*L1*mHd2*Conjg(L1) + 6._dp*L1*mt2*Conjg(L1) -              &
  &  6._dp*g12*MassB*Conjg(MassB)/5._dp - 6._dp*g22*MassWB*Conjg(MassWB)               &
  &  + 6._dp*AL1*Conjg(AL1) - g12*Tr1(1)


  If (TwoLoopRGE) Then
  betamHd22 = -4._dp*g12*TrAYdadjAYd/5._dp + 32._dp*g32*TrAYdadjAYd +               &
  &  12._dp*g12*TrAYeadjAYe/5._dp - 12._dp*TrCYsAYdadjAYdYs - 15._dp*TrCYsYdadjAYdAYs -  &
  &  12._dp*TrCYsYdadjYdmd2Ys - 21._dp*mHd2*TrCYsYdadjYdYs/2._dp - 12._dp*ms2*TrCYsYdadjYdYs -&
  &  12._dp*TrCYsYdadjYdYsmd2 - 12._dp*TrCYsYdmq2adjYdYs - 2._dp*TrCYtAYtadjAYeYe -        &
  &  3._dp*TrCYtYtadjAYeAYe - 3._dp*TrCYtYtadjYeme2Ye - 3._dp*mHd2*TrCYtYtadjYeYe/2._dp -  &
  &  2._dp*mt2*TrCYtYtadjYeYe - 3._dp*TrCYtYtadjYeYeml2 - 2._dp*TrCYtYtml2adjYeYe -        &
  &  4._dp*g12*Trmd2YdadjYd/5._dp + 32._dp*g32*Trmd2YdadjYd - 9._dp*Trmd2YdadjYdYsCYs +&
  &  12._dp*g12*Trme2YeadjYe/5._dp - Trme2YeCYtYtadjYe - 3._dp*Trmq2adjYdYsCYsYd +       &
  &  4._dp*g12*MassB*TrYdadjAYd/5._dp - 32._dp*g32*MassG*TrYdadjAYd - 36._dp*TrYdadjAYdAYdadjYd -&
  &  9._dp*TrYdadjAYdAYsCYs - 12._dp*TrYdadjAYdAYzadjYz - 6._dp*TrYdadjAYuAYuadjYd -       &
  &  4._dp*g12*mHd2*TrYdadjYd/5._dp + 32._dp*g32*mHd2*TrYdadjYd - 36._dp*TrYdadjYdAYdadjAYd -&
  &  24._dp*TrYdadjYdAYsCAYs - 12._dp*TrYdadjYdAYzadjAYz - 36._dp*TrYdadjYdmd2YdadjYd -    &
  &  9._dp*TrYdadjYdmd2YsCYs - 12._dp*TrYdadjYdmd2YzadjYz - 36._dp*mHd2*TrYdadjYdYdadjYd - &
  &  18._dp*TrYdadjYdYdmq2adjYd - 9._dp*mHd2*TrYdadjYdYsCYs - 9._dp*ms2*TrYdadjYdYsCYs -   &
  &  12._dp*mHd2*TrYdadjYdYzadjYz - 12._dp*mz2*TrYdadjYdYzadjYz - 12._dp*TrYdadjYdYzml2adjYz -&
  &  6._dp*TrYdadjYuAYuadjAYd - 6._dp*TrYdadjYumu2YuadjYd - 3._dp*mHd2*TrYdadjYuYuadjYd -  &
  &  4._dp*g12*TrYdmq2adjYd/5._dp + 32._dp*g32*TrYdmq2adjYd - 18._dp*TrYdmq2adjYdYdadjYd -&
  &  9._dp*TrYdmq2adjYdYsCYs - 12._dp*TrYdmq2adjYdYzadjYz - 12._dp*g12*MassB*TrYeadjAYe/5._dp -&
  &  12._dp*TrYeadjAYeAYeadjYe - 6._dp*TrYeadjAYzAYzadjYe + 12._dp*g12*mHd2*TrYeadjYe/5._dp -&
  &  12._dp*TrYeadjYeAYeadjAYe - 12._dp*TrYeadjYeme2YeadjYe - 12._dp*mHd2*TrYeadjYeYeadjYe -&
  &  6._dp*TrYeadjYeYeml2adjYe - 6._dp*TrYeadjYzAYzadjAYe - 6._dp*TrYeadjYzmd2YzadjYe -    &
  &  3._dp*mHd2*TrYeadjYzYzadjYe - 6._dp*TrYeCAYtAYtadjYe - 4._dp*TrYeCYtAYtadjAYe -       &
  &  3._dp*TrYeCYtml2YtadjYe - 7._dp*mHd2*TrYeCYtYtadjYe/2._dp - mt2*TrYeCYtYtadjYe -      &
  &  TrYeCYtYtml2adjYe + 12._dp*g12*TrYeml2adjYe/5._dp - 6._dp*TrYeml2adjYeYeadjYe -     &
  &  TrYeml2CYtYtadjYe - 24._dp*TrYsCAYsAYdadjYd - 12._dp*TrYsCYsAYdadjAYd -               &
  &  15._dp*TrYsCYsmd2YdadjYd - 9._dp*mHd2*TrYsCYsYdadjYd/2._dp - 3._dp*ms2*TrYsCYsYdadjYd -&
  &  3._dp*TrYsCYsYdadjYdmd2 - 12._dp*TrYsmd2CYsYdadjYd - 3._dp*TrYtadjAYeAYeCYt -         &
  &  6._dp*TrYtadjYeAYeCAYt - 2._dp*TrYtadjYeme2YeCYt - mHd2*TrYtadjYeYeCYt -              &
  &  3._dp*mt2*TrYtadjYeYeCYt - 3._dp*TrYtadjYeYeCYtml2 - 2._dp*TrYtadjYeYeml2CYt -        &
  &  3._dp*TrYtml2adjYeYeCYt - 6._dp*TrYuadjAYdAYdadjYu - 6._dp*TrYuadjYdAYdadjAYu -       &
  &  6._dp*TrYuadjYdmd2YdadjYu - 3._dp*mHd2*TrYuadjYdYdadjYu - 6._dp*mHu2*TrYuadjYdYdadjYu -&
  &  6._dp*TrYuadjYdYdmq2adjYu - 6._dp*TrYumq2adjYdYdadjYu - 6._dp*TrYzadjAYeAYeadjYz -    &
  &  12._dp*TrYzadjAYzAYdadjYd - 6._dp*TrYzadjYeAYeadjAYz - 6._dp*TrYzadjYeme2YeadjYz -    &
  &  3._dp*mHd2*TrYzadjYeYeadjYz - 6._dp*mz2*TrYzadjYeYeadjYz - 6._dp*TrYzadjYeYeml2adjYz -&
  &  12._dp*TrYzadjYzAYdadjAYd - 12._dp*TrYzadjYzmd2YdadjYd - 6._dp*TrYzml2adjYeYeadjYz +  &
  &  72._dp*g12*L1*mHd2*Conjg(L1)/5._dp + 48._dp*g22*L1*mHd2*Conjg(L1) +               &
  &  36._dp*g12*L1*mt2*Conjg(L1)/5._dp + 24._dp*g22*L1*mt2*Conjg(L1) - 18._dp*L1*TrAYdadjAYd*Conjg(L1)&
  &  - 6._dp*L1*TrAYeadjAYe*Conjg(L1) - 3._dp*L1*TrAYtCAYt*Conjg(L1)/2._dp -               &
  &  9._dp*L1*TrCAYtAYt*Conjg(L1)/2._dp - 3._dp*L1*TrCYtml2Yt*Conjg(L1) - 9._dp*L1*mHd2*TrCYtYt*Conjg(L1)&
  &  - 9._dp*L1*mt2*TrCYtYt*Conjg(L1) - 6._dp*L1*TrCYtYtml2*Conjg(L1) - 18._dp*L1*Trmd2YdadjYd*Conjg(L1)&
  &  - 6._dp*L1*Trme2YeadjYe*Conjg(L1) - 54._dp*L1*mHd2*TrYdadjYd*Conjg(L1) -              &
  &  18._dp*L1*mt2*TrYdadjYd*Conjg(L1) - 18._dp*L1*TrYdmq2adjYd*Conjg(L1) - 18._dp*L1*mHd2*TrYeadjYe*Conjg(L1)&
  &  - 6._dp*L1*mt2*TrYeadjYe*Conjg(L1) - 6._dp*L1*TrYeml2adjYe*Conjg(L1) - 3._dp*L1*mHd2*TrYtCYt*Conjg(L1)&
  &  - 3._dp*L1*mt2*TrYtCYt*Conjg(L1) - 3._dp*L1*TrYtCYtml2*Conjg(L1) - 18._dp*TrYdadjAYd*AL1*Conjg(L1)&
  &  - 6._dp*TrYeadjAYe*AL1*Conjg(L1) - 6._dp*TrYtCAYt*AL1*Conjg(L1) - 96._dp*L1**2*mHd2*Conjg(L1)&
  & **2 - 48._dp*L1**2*mt2*Conjg(L1)**2 + 1251._dp*g14*MassB*Conjg(MassB)/25._dp +       &
  &  18._dp*g12*g22*MassB*Conjg(MassB)/5._dp + 9._dp*g12*g22*MassWB*Conjg(MassB)   &
  & /5._dp + 4._dp*g12*TrAYdadjYd*Conjg(MassB)/5._dp - 12._dp*g12*TrAYeadjYe*Conjg(MassB)&
  & /5._dp - 8._dp*g12*MassB*TrYdadjYd*Conjg(MassB)/5._dp + 24._dp*g12*MassB*TrYeadjYe*Conjg(MassB)&
  & /5._dp + 72._dp*g12*L1*MassB*Conjg(L1)*Conjg(MassB)/5._dp - 36._dp*g12*AL1*Conjg(L1)&
  & *Conjg(MassB)/5._dp - 32._dp*g32*TrAYdadjYd*Conjg(MassG) + 64._dp*g32*MassG*TrYdadjYd*Conjg(MassG)&
  &  + 9._dp*g12*g22*MassB*Conjg(MassWB)/5._dp + 18._dp*g12*g22*MassWB*Conjg(MassWB)&
  & /5._dp + 159._dp*g24*MassWB*Conjg(MassWB) + 48._dp*g22*L1*MassWB*Conjg(L1)         &
  & *Conjg(MassWB) - 24._dp*g22*AL1*Conjg(L1)*Conjg(MassWB) - 36._dp*g12*L1*MassB*Conjg(AL1)&
  & /5._dp - 24._dp*g22*L1*MassWB*Conjg(AL1) - 18._dp*L1*TrAYdadjYd*Conjg(AL1)           &
  &  - 6._dp*L1*TrAYeadjYe*Conjg(AL1) - 6._dp*L1*TrCYtAYt*Conjg(AL1) + 36._dp*g12*AL1*Conjg(AL1)&
  & /5._dp + 24._dp*g22*AL1*Conjg(AL1) - 9._dp*TrCYtYt*AL1*Conjg(AL1)/2._dp -            &
  &  18._dp*TrYdadjYd*AL1*Conjg(AL1) - 6._dp*TrYeadjYe*AL1*Conjg(AL1) - 3._dp*TrYtCYt*AL1*Conjg(AL1)&
  & /2._dp - 96._dp*L1*AL1*Conjg(L1)*Conjg(AL1) + 6._dp*g14*Tr2(1)/5._dp +               &
  &  6._dp*g24*Tr2(2) - 4._dp*g12*Tr3(1)


  DmHd2 = oo16pi2*( betamHd21 + oo16pi2 * betamHd22 )


  Else
  DmHd2 = oo16pi2* betamHd21
  End If


  !--------------------
  ! mHu2
  !--------------------

  betamHu21 = 6._dp*TrAYuadjAYu + 6._dp*Trmu2YuadjYu + 6._dp*mHu2*TrYuadjYu +           &
  &  6._dp*TrYumq2adjYu + 12._dp*L2*mHu2*Conjg(L2) + 6._dp*L2*mtb2*Conjg(L2)               &
  &  - 6._dp*g12*MassB*Conjg(MassB)/5._dp - 6._dp*g22*MassWB*Conjg(MassWB)             &
  &  + 6._dp*AL2*Conjg(AL2) + g12*Tr1(1)


  If (TwoLoopRGE) Then
  betamHu22 = 8._dp*g12*TrAYuadjAYu/5._dp + 32._dp*g32*TrAYuadjAYu + 8._dp*g12*Trmu2YuadjYu/5._dp +&
  &  32._dp*g32*Trmu2YuadjYu - 6._dp*TrYdadjAYuAYuadjYd - 6._dp*TrYdadjYuAYuadjAYd -     &
  &  6._dp*TrYdadjYumu2YuadjYd - 6._dp*mHd2*TrYdadjYuYuadjYd - 3._dp*mHu2*TrYdadjYuYuadjYd -&
  &  6._dp*TrYdadjYuYumq2adjYd - 6._dp*TrYdmq2adjYuYuadjYd - 6._dp*TrYuadjAYdAYdadjYu -    &
  &  8._dp*g12*MassB*TrYuadjAYu/5._dp - 32._dp*g32*MassG*TrYuadjAYu - 36._dp*TrYuadjAYuAYuadjYu -&
  &  6._dp*TrYuadjYdAYdadjAYu - 6._dp*TrYuadjYdmd2YdadjYu - 3._dp*mHu2*TrYuadjYdYdadjYu +  &
  &  8._dp*g12*mHu2*TrYuadjYu/5._dp + 32._dp*g32*mHu2*TrYuadjYu - 36._dp*TrYuadjYuAYuadjAYu -&
  &  36._dp*TrYuadjYumu2YuadjYu - 36._dp*mHu2*TrYuadjYuYuadjYu - 18._dp*TrYuadjYuYumq2adjYu +&
  &  8._dp*g12*TrYumq2adjYu/5._dp + 32._dp*g32*TrYumq2adjYu - 18._dp*TrYumq2adjYuYuadjYu +&
  &  72._dp*g12*L2*mHu2*Conjg(L2)/5._dp + 48._dp*g22*L2*mHu2*Conjg(L2) +               &
  &  36._dp*g12*L2*mtb2*Conjg(L2)/5._dp + 24._dp*g22*L2*mtb2*Conjg(L2) -               &
  &  18._dp*L2*TrAYuadjAYu*Conjg(L2) - 18._dp*L2*Trmu2YuadjYu*Conjg(L2) - 54._dp*L2*mHu2*TrYuadjYu*Conjg(L2)&
  &  - 18._dp*L2*mtb2*TrYuadjYu*Conjg(L2) - 18._dp*L2*TrYumq2adjYu*Conjg(L2)               &
  &  - 18._dp*TrYuadjAYu*AL2*Conjg(L2) - 96._dp*L2**2*mHu2*Conjg(L2)**2 - 48._dp*L2**2*mtb2*Conjg(L2)&
  & **2 + 1251._dp*g14*MassB*Conjg(MassB)/25._dp + 18._dp*g12*g22*MassB*Conjg(MassB) &
  & /5._dp + 9._dp*g12*g22*MassWB*Conjg(MassB)/5._dp - 8._dp*g12*TrAYuadjYu*Conjg(MassB)&
  & /5._dp + 16._dp*g12*MassB*TrYuadjYu*Conjg(MassB)/5._dp + 72._dp*g12*L2*MassB*Conjg(L2)&
  & *Conjg(MassB)/5._dp - 36._dp*g12*AL2*Conjg(L2)*Conjg(MassB)/5._dp - 32._dp*g32*TrAYuadjYu*Conjg(MassG)&
  &  + 64._dp*g32*MassG*TrYuadjYu*Conjg(MassG) + 9._dp*g12*g22*MassB*Conjg(MassWB)   &
  & /5._dp + 18._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g24*MassWB*Conjg(MassWB)&
  &  + 48._dp*g22*L2*MassWB*Conjg(L2)*Conjg(MassWB) - 24._dp*g22*AL2*Conjg(L2)         &
  & *Conjg(MassWB) - 36._dp*g12*L2*MassB*Conjg(AL2)/5._dp - 24._dp*g22*L2*MassWB*Conjg(AL2)&
  &  - 18._dp*L2*TrAYuadjYu*Conjg(AL2) + 36._dp*g12*AL2*Conjg(AL2)/5._dp +               &
  &  24._dp*g22*AL2*Conjg(AL2) - 18._dp*TrYuadjYu*AL2*Conjg(AL2) - 96._dp*L2*AL2*Conjg(L2)&
  & *Conjg(AL2) + 6._dp*g14*Tr2(1)/5._dp + 6._dp*g24*Tr2(2) + 4._dp*g12*Tr3(1)


  DmHu2 = oo16pi2*( betamHu21 + oo16pi2 * betamHu22 )


  Else
  DmHu2 = oo16pi2* betamHu21
  End If


  !--------------------
  ! md2
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betamd21(i1,i2) = 4._dp*AYdadjAYd(i1,i2) + 8._dp*AYsCAYs(i1,i2) + 4._dp*AYzadjAYz(i1,i2)&
  &  + 2._dp*md2YdadjYd(i1,i2) + 4._dp*md2YsCYs(i1,i2) + 2._dp*md2YzadjYz(i1,i2)           &
  &  + 4._dp*mHd2*YdadjYd(i1,i2) + 2._dp*YdadjYdmd2(i1,i2) + 4._dp*Ydmq2adjYd(i1,i2)       &
  &  + 8._dp*ms2*YsCYs(i1,i2) + 4._dp*YsCYsmd2(i1,i2) + 8._dp*Ysmd2CYs(i1,i2)              &
  &  + 4._dp*mz2*YzadjYz(i1,i2) + 2._dp*YzadjYzmd2(i1,i2) + 4._dp*Yzml2adjYz(i1,i2)

  If (i1.eq.i2) Then
  betamd21(i1,i2) = betamd21(i1,i2)+(-8._dp*g12*MassB*Conjg(MassB)/15._dp - 32._dp*g32*MassG*Conjg(MassG)&
  & /3._dp + 2._dp*g12*Tr1(1)/3._dp)
  End If
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betamd22(i1,i2) = 4._dp*g12*AYdadjAYd(i1,i2)/5._dp + 12._dp*g22*AYdadjAYd(i1,i2)  &
  &  - 12._dp*TrYdadjYd*AYdadjAYd(i1,i2) - 4._dp*TrYeadjYe*AYdadjAYd(i1,i2) -              &
  &  4._dp*AYdadjAYdYdadjYd(i1,i2) - 4._dp*AYdadjAYuYuadjYd(i1,i2) - 12._dp*TrYdadjAYd*AYdadjYd(i1,i2)&
  &  - 4._dp*TrYeadjAYe*AYdadjYd(i1,i2) - 4._dp*AYdadjYdYdadjAYd(i1,i2) - 4._dp*AYdadjYuYuadjAYd(i1,i2)&
  &  - 16._dp*AYsCAYdTYdCYs(i1,i2) + 64._dp*g12*AYsCAYs(i1,i2)/15._dp + 160._dp*g32*AYsCAYs(i1,i2)&
  & /3._dp - 6._dp*TrCYsYs*AYsCAYs(i1,i2) - 2._dp*TrYsCYs*AYsCAYs(i1,i2) - 32._dp*AYsCAYsYsCYs(i1,i2)&
  &  - 16._dp*AYsCAYzTYzCYs(i1,i2) - 16._dp*AYsCYdTYdCAYs(i1,i2) - 8._dp*TrYsCAYs*AYsCYs(i1,i2)&
  &  - 32._dp*AYsCYsYsCAYs(i1,i2) - 16._dp*AYsCYzTYzCAYs(i1,i2) - 4._dp*AYzadjAYeYeadjYz(i1,i2)&
  &  + 4._dp*g12*AYzadjAYz(i1,i2)/5._dp + 12._dp*g22*AYzadjAYz(i1,i2) - 4._dp*TrYzadjYz*AYzadjAYz(i1,i2)&
  &  - 12._dp*AYzadjAYzYzadjYz(i1,i2) - 4._dp*AYzadjYeYeadjAYz(i1,i2) - 4._dp*TrYzadjAYz*AYzadjYz(i1,i2)&
  &  - 12._dp*AYzadjYzYzadjAYz(i1,i2) - 12._dp*AYzCAYtYtadjYz(i1,i2) - 12._dp*AYzCYtYtadjAYz(i1,i2)&
  &  - 12._dp*L1*AYdadjAYd(i1,i2)*Conjg(L1) - 4._dp*g12*AYdadjYd(i1,i2)*Conjg(MassB)     &
  & /5._dp - 64._dp*g12*AYsCYs(i1,i2)*Conjg(MassB)/15._dp - 4._dp*g12*AYzadjYz(i1,i2)  &
  & *Conjg(MassB)/5._dp - 160._dp*g32*AYsCYs(i1,i2)*Conjg(MassG)/3._dp - 12._dp*g22*AYdadjYd(i1,i2)&
  & *Conjg(MassWB) - 12._dp*g22*AYzadjYz(i1,i2)*Conjg(MassWB) - 12._dp*L1*AYdadjYd(i1,i2)&
  & *Conjg(AL1) + 2._dp*g12*md2YdadjYd(i1,i2)/5._dp + 6._dp*g22*md2YdadjYd(i1,i2)      &
  &  - 6._dp*TrYdadjYd*md2YdadjYd(i1,i2) - 2._dp*TrYeadjYe*md2YdadjYd(i1,i2)               &
  &  - 6._dp*AbsL1sq*md2YdadjYd(i1,i2) - 2._dp*md2YdadjYdYdadjYd(i1,i2)               &
  &  - 2._dp*md2YdadjYuYuadjYd(i1,i2) - 8._dp*md2YsCYdTYdCYs(i1,i2) + 32._dp*g12*md2YsCYs(i1,i2)&
  & /15._dp + 80._dp*g32*md2YsCYs(i1,i2)/3._dp - 3._dp*TrCYsYs*md2YsCYs(i1,i2)           &
  &  - TrYsCYs*md2YsCYs(i1,i2) - 16._dp*md2YsCYsYsCYs(i1,i2) - 8._dp*md2YsCYzTYzCYs(i1,i2) &
  &  - 2._dp*md2YzadjYeYeadjYz(i1,i2) + 2._dp*g12*md2YzadjYz(i1,i2)/5._dp +              &
  &  6._dp*g22*md2YzadjYz(i1,i2) - 2._dp*TrYzadjYz*md2YzadjYz(i1,i2) - 6._dp*md2YzadjYzYzadjYz(i1,i2)&
  &  - 6._dp*md2YzCYtYtadjYz(i1,i2) - 4._dp*g12*MassB*YdadjAYd(i1,i2)/5._dp -            &
  &  12._dp*g22*MassWB*YdadjAYd(i1,i2) - 12._dp*TrAYdadjYd*YdadjAYd(i1,i2)               &
  &  - 4._dp*TrAYeadjYe*YdadjAYd(i1,i2) - 12._dp*AL1*Conjg(L1)*YdadjAYd(i1,i2)             &
  &  - 4._dp*YdadjAYdAYdadjYd(i1,i2) - 4._dp*YdadjAYuAYuadjYd(i1,i2) + 4._dp*g12*mHd2*YdadjYd(i1,i2)&
  & /5._dp + 12._dp*g22*mHd2*YdadjYd(i1,i2) - 12._dp*TrAYdadjAYd*YdadjYd(i1,i2)          &
  &  - 4._dp*TrAYeadjAYe*YdadjYd(i1,i2) - 12._dp*Trmd2YdadjYd*YdadjYd(i1,i2)               &
  &  - 4._dp*Trme2YeadjYe*YdadjYd(i1,i2) - 24._dp*mHd2*TrYdadjYd*YdadjYd(i1,i2)            &
  &  - 12._dp*TrYdmq2adjYd*YdadjYd(i1,i2) - 8._dp*mHd2*TrYeadjYe*YdadjYd(i1,i2)            &
  &  - 4._dp*TrYeml2adjYe*YdadjYd(i1,i2) - 36._dp*L1*mHd2*Conjg(L1)*YdadjYd(i1,i2)         &
  &  - 12._dp*L1*mt2*Conjg(L1)*YdadjYd(i1,i2) + 8._dp*g12*MassB*Conjg(MassB)             &
  & *YdadjYd(i1,i2)/5._dp + 24._dp*g22*MassWB*Conjg(MassWB)*YdadjYd(i1,i2)               &
  &  - 12._dp*AL1*Conjg(AL1)*YdadjYd(i1,i2) - 4._dp*YdadjYdAYdadjAYd(i1,i2) +              &
  &  2._dp*g12*YdadjYdmd2(i1,i2)/5._dp + 6._dp*g22*YdadjYdmd2(i1,i2) - 6._dp*TrYdadjYd*YdadjYdmd2(i1,i2)&
  &  - 2._dp*TrYeadjYe*YdadjYdmd2(i1,i2) - 6._dp*AbsL1sq*YdadjYdmd2(i1,i2)            &
  &  - 4._dp*YdadjYdmd2YdadjYd(i1,i2) - 8._dp*mHd2*YdadjYdYdadjYd(i1,i2) - 2._dp*YdadjYdYdadjYdmd2(i1,i2)&
  &  - 4._dp*YdadjYdYdmq2adjYd(i1,i2) - 4._dp*YdadjYuAYuadjAYd(i1,i2) - 4._dp*YdadjYumu2YuadjYd(i1,i2)&
  &  - 4._dp*mHd2*YdadjYuYuadjYd(i1,i2) - 4._dp*mHu2*YdadjYuYuadjYd(i1,i2) -               &
  &  2._dp*YdadjYuYuadjYdmd2(i1,i2) - 4._dp*YdadjYuYumq2adjYd(i1,i2) + 4._dp*g12*Ydmq2adjYd(i1,i2)&
  & /5._dp + 12._dp*g22*Ydmq2adjYd(i1,i2) - 12._dp*TrYdadjYd*Ydmq2adjYd(i1,i2)           &
  &  - 4._dp*TrYeadjYe*Ydmq2adjYd(i1,i2) - 12._dp*AbsL1sq*Ydmq2adjYd(i1,i2)           &
  &  - 4._dp*Ydmq2adjYdYdadjYd(i1,i2) - 4._dp*Ydmq2adjYuYuadjYd(i1,i2) - 16._dp*YsCAYdTAYdCYs(i1,i2)&
  &  - 64._dp*g12*MassB*YsCAYs(i1,i2)/15._dp - 160._dp*g32*MassG*YsCAYs(i1,i2)         &
  & /3._dp - 8._dp*TrCYsAYs*YsCAYs(i1,i2) - 32._dp*YsCAYsAYsCYs(i1,i2) - 16._dp*YsCAYzTAYzCYs(i1,i2)&
  &  - 16._dp*YsCYdmq2TYdCYs(i1,i2) - 16._dp*YsCYdTAYdCAYs(i1,i2) - 16._dp*mHd2*YsCYdTYdCYs(i1,i2)&
  &  - 16._dp*ms2*YsCYdTYdCYs(i1,i2) - 8._dp*YsCYdTYdCYsmd2(i1,i2) - 16._dp*YsCYdTYdmd2CYs(i1,i2)&
  &  + 64._dp*g12*ms2*YsCYs(i1,i2)/15._dp + 160._dp*g32*ms2*YsCYs(i1,i2)               &
  & /3._dp - 2._dp*TrAYsCAYs*YsCYs(i1,i2) - 6._dp*TrCAYsAYs*YsCYs(i1,i2) - 4._dp*TrCYsmd2Ys*YsCYs(i1,i2)&
  &  - 12._dp*ms2*TrCYsYs*YsCYs(i1,i2) - 8._dp*TrCYsYsmd2*YsCYs(i1,i2) - 4._dp*ms2*TrYsCYs*YsCYs(i1,i2)&
  &  - 4._dp*TrYsCYsmd2*YsCYs(i1,i2) + 128._dp*g12*MassB*Conjg(MassB)*YsCYs(i1,i2)       &
  & /15._dp + 320._dp*g32*MassG*Conjg(MassG)*YsCYs(i1,i2)/3._dp - 32._dp*YsCYsAYsCAYs(i1,i2)&
  &  + 32._dp*g12*YsCYsmd2(i1,i2)/15._dp + 80._dp*g32*YsCYsmd2(i1,i2)/3._dp -          &
  &  3._dp*TrCYsYs*YsCYsmd2(i1,i2) - TrYsCYs*YsCYsmd2(i1,i2) - 32._dp*YsCYsmd2YsCYs(i1,i2) &
  &  - 64._dp*ms2*YsCYsYsCYs(i1,i2) - 16._dp*YsCYsYsCYsmd2(i1,i2) - 32._dp*YsCYsYsmd2CYs(i1,i2)&
  &  - 16._dp*YsCYzml2TYzCYs(i1,i2) - 16._dp*YsCYzTAYzCAYs(i1,i2) - 16._dp*ms2*YsCYzTYzCYs(i1,i2)&
  &  - 16._dp*mz2*YsCYzTYzCYs(i1,i2) - 8._dp*YsCYzTYzCYsmd2(i1,i2) - 16._dp*YsCYzTYzmd2CYs(i1,i2)&
  &  - 16._dp*Ysmd2CYdTYdCYs(i1,i2) + 64._dp*g12*Ysmd2CYs(i1,i2)/15._dp + 160._dp*g32*Ysmd2CYs(i1,i2)&
  & /3._dp - 6._dp*TrCYsYs*Ysmd2CYs(i1,i2) - 2._dp*TrYsCYs*Ysmd2CYs(i1,i2) -               &
  &  32._dp*Ysmd2CYsYsCYs(i1,i2) - 16._dp*Ysmd2CYzTYzCYs(i1,i2) - 4._dp*YzadjAYeAYeadjYz(i1,i2)&
  &  - 4._dp*g12*MassB*YzadjAYz(i1,i2)/5._dp - 12._dp*g22*MassWB*YzadjAYz(i1,i2)       &
  &  - 4._dp*TrAYzadjYz*YzadjAYz(i1,i2) - 12._dp*YzadjAYzAYzadjYz(i1,i2) - 4._dp*YzadjYeAYeadjAYz(i1,i2)&
  &  - 4._dp*YzadjYeme2YeadjYz(i1,i2) - 4._dp*mHd2*YzadjYeYeadjYz(i1,i2) - 4._dp*mz2*YzadjYeYeadjYz(i1,i2)&
  &  - 2._dp*YzadjYeYeadjYzmd2(i1,i2) - 4._dp*YzadjYeYeml2adjYz(i1,i2) + 4._dp*g12*mz2*YzadjYz(i1,i2)&
  & /5._dp + 12._dp*g22*mz2*YzadjYz(i1,i2) - 4._dp*TrAYzadjAYz*YzadjYz(i1,i2)            &
  &  - 4._dp*Trmd2YzadjYz*YzadjYz(i1,i2) - 8._dp*mz2*TrYzadjYz*YzadjYz(i1,i2)              &
  &  - 4._dp*TrYzml2adjYz*YzadjYz(i1,i2) + 8._dp*g12*MassB*Conjg(MassB)*YzadjYz(i1,i2)   &
  & /5._dp + 24._dp*g22*MassWB*Conjg(MassWB)*YzadjYz(i1,i2) - 12._dp*YzadjYzAYzadjAYz(i1,i2)&
  &  + 2._dp*g12*YzadjYzmd2(i1,i2)/5._dp + 6._dp*g22*YzadjYzmd2(i1,i2) -               &
  &  2._dp*TrYzadjYz*YzadjYzmd2(i1,i2) - 12._dp*YzadjYzmd2YzadjYz(i1,i2) - 24._dp*mz2*YzadjYzYzadjYz(i1,i2)&
  &  - 6._dp*YzadjYzYzadjYzmd2(i1,i2) - 12._dp*YzadjYzYzml2adjYz(i1,i2) - 12._dp*YzCAYtAYtadjYz(i1,i2)&
  &  - 12._dp*YzCYtAYtadjAYz(i1,i2) - 12._dp*YzCYtml2YtadjYz(i1,i2) - 12._dp*mt2*YzCYtYtadjYz(i1,i2)&
  &  - 12._dp*mz2*YzCYtYtadjYz(i1,i2) - 6._dp*YzCYtYtadjYzmd2(i1,i2) - 12._dp*YzCYtYtml2adjYz(i1,i2)&
  &  - 4._dp*Yzml2adjYeYeadjYz(i1,i2) + 4._dp*g12*Yzml2adjYz(i1,i2)/5._dp +              &
  &  12._dp*g22*Yzml2adjYz(i1,i2) - 4._dp*TrYzadjYz*Yzml2adjYz(i1,i2) - 12._dp*Yzml2adjYzYzadjYz(i1,i2)&
  &  - 12._dp*Yzml2CYtYtadjYz(i1,i2)

  If (i1.eq.i2) Then
  betamd22(i1,i2) = betamd22(i1,i2)+(1648._dp*g14*MassB*Conjg(MassB)/75._dp + 128._dp*g12*g32*MassB*Conjg(MassB)&
  & /45._dp + 64._dp*g12*g32*MassG*Conjg(MassB)/45._dp + 64._dp*g12*g32*MassB*Conjg(MassG)&
  & /45._dp + 128._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 544._dp*g34*MassG*Conjg(MassG)&
  & /3._dp + 8._dp*g14*Tr2(1)/15._dp + 32._dp*g34*Tr2(3)/3._dp + 8._dp*g12*Tr3(1)    &
  & /3._dp - 4._dp*g32*Tr3(3) - 4._dp*(g32*Tr3(3))/sqrt3)
  End If
   End Do
  End Do


  Dmd2 = oo16pi2*( betamd21 + oo16pi2 * betamd22 )


  Else
  Dmd2 = oo16pi2* betamd21
  End If


  !--------------------
  ! mu2
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betamu21(i1,i2) = 4._dp*AYuadjAYu(i1,i2) + 2._dp*mu2YuadjYu(i1,i2) + 4._dp*mHu2*YuadjYu(i1,i2)&
  &  + 2._dp*YuadjYumu2(i1,i2) + 4._dp*Yumq2adjYu(i1,i2)

  If (i1.eq.i2) Then
  betamu21(i1,i2) = betamu21(i1,i2)+(-32._dp*g12*MassB*Conjg(MassB)/15._dp - 32._dp*g32*MassG*Conjg(MassG)&
  & /3._dp - 4._dp*g12*Tr1(1)/3._dp)
  End If
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betamu22(i1,i2) = -4._dp*AYuadjAYdYdadjYu(i1,i2) - 4._dp*g12*AYuadjAYu(i1,i2)       &
  & /5._dp + 12._dp*g22*AYuadjAYu(i1,i2) - 12._dp*TrYuadjYu*AYuadjAYu(i1,i2)             &
  &  - 4._dp*AYuadjAYuYuadjYu(i1,i2) - 4._dp*AYuadjYdYdadjAYu(i1,i2) - 12._dp*TrYuadjAYu*AYuadjYu(i1,i2)&
  &  - 4._dp*AYuadjYuYuadjAYu(i1,i2) - 12._dp*L2*AYuadjAYu(i1,i2)*Conjg(L2) +              &
  &  4._dp*g12*AYuadjYu(i1,i2)*Conjg(MassB)/5._dp - 12._dp*g22*AYuadjYu(i1,i2)         &
  & *Conjg(MassWB) - 12._dp*L2*AYuadjYu(i1,i2)*Conjg(AL2) - 2._dp*mu2YuadjYdYdadjYu(i1,i2) &
  &  - 2._dp*g12*mu2YuadjYu(i1,i2)/5._dp + 6._dp*g22*mu2YuadjYu(i1,i2) -               &
  &  6._dp*TrYuadjYu*mu2YuadjYu(i1,i2) - 6._dp*AbsL2sq*mu2YuadjYu(i1,i2)              &
  &  - 2._dp*mu2YuadjYuYuadjYu(i1,i2) - 4._dp*YuadjAYdAYdadjYu(i1,i2) + 4._dp*g12*MassB*YuadjAYu(i1,i2)&
  & /5._dp - 12._dp*g22*MassWB*YuadjAYu(i1,i2) - 12._dp*TrAYuadjYu*YuadjAYu(i1,i2)       &
  &  - 12._dp*AL2*Conjg(L2)*YuadjAYu(i1,i2) - 4._dp*YuadjAYuAYuadjYu(i1,i2) -              &
  &  4._dp*YuadjYdAYdadjAYu(i1,i2) - 4._dp*YuadjYdmd2YdadjYu(i1,i2) - 4._dp*mHd2*YuadjYdYdadjYu(i1,i2)&
  &  - 4._dp*mHu2*YuadjYdYdadjYu(i1,i2) - 2._dp*YuadjYdYdadjYumu2(i1,i2) - 4._dp*YuadjYdYdmq2adjYu(i1,i2)&
  &  - 4._dp*g12*mHu2*YuadjYu(i1,i2)/5._dp + 12._dp*g22*mHu2*YuadjYu(i1,i2)            &
  &  - 12._dp*TrAYuadjAYu*YuadjYu(i1,i2) - 12._dp*Trmu2YuadjYu*YuadjYu(i1,i2)              &
  &  - 24._dp*mHu2*TrYuadjYu*YuadjYu(i1,i2) - 12._dp*TrYumq2adjYu*YuadjYu(i1,i2)           &
  &  - 36._dp*L2*mHu2*Conjg(L2)*YuadjYu(i1,i2) - 12._dp*L2*mtb2*Conjg(L2)*YuadjYu(i1,i2)   &
  &  - 8._dp*g12*MassB*Conjg(MassB)*YuadjYu(i1,i2)/5._dp + 24._dp*g22*MassWB*Conjg(MassWB)&
  & *YuadjYu(i1,i2) - 12._dp*AL2*Conjg(AL2)*YuadjYu(i1,i2) - 4._dp*YuadjYuAYuadjAYu(i1,i2) &
  &  - 2._dp*g12*YuadjYumu2(i1,i2)/5._dp + 6._dp*g22*YuadjYumu2(i1,i2) -               &
  &  6._dp*TrYuadjYu*YuadjYumu2(i1,i2) - 6._dp*AbsL2sq*YuadjYumu2(i1,i2)              &
  &  - 4._dp*YuadjYumu2YuadjYu(i1,i2) - 8._dp*mHu2*YuadjYuYuadjYu(i1,i2) - 2._dp*YuadjYuYuadjYumu2(i1,i2)&
  &  - 4._dp*YuadjYuYumq2adjYu(i1,i2) - 4._dp*Yumq2adjYdYdadjYu(i1,i2) - 4._dp*g12*Yumq2adjYu(i1,i2)&
  & /5._dp + 12._dp*g22*Yumq2adjYu(i1,i2) - 12._dp*TrYuadjYu*Yumq2adjYu(i1,i2)           &
  &  - 12._dp*AbsL2sq*Yumq2adjYu(i1,i2) - 4._dp*Yumq2adjYuYuadjYu(i1,i2)

  If (i1.eq.i2) Then
  betamu22(i1,i2) = betamu22(i1,i2)+(6784._dp*g14*MassB*Conjg(MassB)/75._dp + 512._dp*g12*g32*MassB*Conjg(MassB)&
  & /45._dp + 256._dp*g12*g32*MassG*Conjg(MassB)/45._dp + 256._dp*g12*g32*MassB*Conjg(MassG)&
  & /45._dp + 512._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 544._dp*g34*MassG*Conjg(MassG)&
  & /3._dp + 32._dp*g14*Tr2(1)/15._dp + 32._dp*g34*Tr2(3)/3._dp - 16._dp*g12*Tr3(1)  &
  & /3._dp - 4._dp*g32*Tr3(3) - 4._dp*(g32*Tr3(3))/sqrt3)
  End If
   End Do
  End Do


  Dmu2 = oo16pi2*( betamu21 + oo16pi2 * betamu22 )


  Else
  Dmu2 = oo16pi2* betamu21
  End If


  !--------------------
  ! me2
  !--------------------

  Do i1 = 1,3
  Do i2 = 1,3
  betame21(i1,i2) = 4._dp*AYeadjAYe(i1,i2) + 2._dp*me2YeadjYe(i1,i2) + 4._dp*mHd2*YeadjYe(i1,i2)&
  &  + 2._dp*YeadjYeme2(i1,i2) + 4._dp*Yeml2adjYe(i1,i2)

  If (i1.eq.i2) Then
  betame21(i1,i2) = betame21(i1,i2)+(-24._dp*g12*MassB*Conjg(MassB)/5._dp + 2._dp*g12*Tr1(1))
  End If
   End Do
  End Do


  If (TwoLoopRGE) Then
  Do i1 = 1,3
  Do i2 = 1,3
  betame22(i1,i2) = -12._dp*g12*AYeadjAYe(i1,i2)/5._dp + 12._dp*g22*AYeadjAYe(i1,i2)&
  &  - 12._dp*TrYdadjYd*AYeadjAYe(i1,i2) - 4._dp*TrYeadjYe*AYeadjAYe(i1,i2) -              &
  &  4._dp*AYeadjAYeYeadjYe(i1,i2) - 12._dp*AYeadjAYzYzadjYe(i1,i2) - 12._dp*TrYdadjAYd*AYeadjYe(i1,i2)&
  &  - 4._dp*TrYeadjAYe*AYeadjYe(i1,i2) - 4._dp*AYeadjYeYeadjAYe(i1,i2) - 12._dp*AYeadjYzYzadjAYe(i1,i2)&
  &  - 12._dp*AYeCAYtYtadjYe(i1,i2) - 12._dp*AYeCYtYtadjAYe(i1,i2) - 12._dp*L1*AYeadjAYe(i1,i2)&
  & *Conjg(L1) + 12._dp*g12*AYeadjYe(i1,i2)*Conjg(MassB)/5._dp - 12._dp*g22*AYeadjYe(i1,i2)&
  & *Conjg(MassWB) - 12._dp*L1*AYeadjYe(i1,i2)*Conjg(AL1) - 6._dp*g12*me2YeadjYe(i1,i2)  &
  & /5._dp + 6._dp*g22*me2YeadjYe(i1,i2) - 6._dp*TrYdadjYd*me2YeadjYe(i1,i2)             &
  &  - 2._dp*TrYeadjYe*me2YeadjYe(i1,i2) - 6._dp*AbsL1sq*me2YeadjYe(i1,i2)            &
  &  - 2._dp*me2YeadjYeYeadjYe(i1,i2) - 6._dp*me2YeadjYzYzadjYe(i1,i2) - 6._dp*me2YeCYtYtadjYe(i1,i2)&
  &  + 12._dp*g12*MassB*YeadjAYe(i1,i2)/5._dp - 12._dp*g22*MassWB*YeadjAYe(i1,i2)      &
  &  - 12._dp*TrAYdadjYd*YeadjAYe(i1,i2) - 4._dp*TrAYeadjYe*YeadjAYe(i1,i2) -              &
  &  12._dp*AL1*Conjg(L1)*YeadjAYe(i1,i2) - 4._dp*YeadjAYeAYeadjYe(i1,i2) - 12._dp*YeadjAYzAYzadjYe(i1,i2)&
  &  - 12._dp*g12*mHd2*YeadjYe(i1,i2)/5._dp + 12._dp*g22*mHd2*YeadjYe(i1,i2)           &
  &  - 12._dp*TrAYdadjAYd*YeadjYe(i1,i2) - 4._dp*TrAYeadjAYe*YeadjYe(i1,i2) -              &
  &  12._dp*Trmd2YdadjYd*YeadjYe(i1,i2) - 4._dp*Trme2YeadjYe*YeadjYe(i1,i2) -              &
  &  24._dp*mHd2*TrYdadjYd*YeadjYe(i1,i2) - 12._dp*TrYdmq2adjYd*YeadjYe(i1,i2)             &
  &  - 8._dp*mHd2*TrYeadjYe*YeadjYe(i1,i2) - 4._dp*TrYeml2adjYe*YeadjYe(i1,i2)             &
  &  - 36._dp*L1*mHd2*Conjg(L1)*YeadjYe(i1,i2) - 12._dp*L1*mt2*Conjg(L1)*YeadjYe(i1,i2)    &
  &  - 24._dp*g12*MassB*Conjg(MassB)*YeadjYe(i1,i2)/5._dp + 24._dp*g22*MassWB*Conjg(MassWB)&
  & *YeadjYe(i1,i2) - 12._dp*AL1*Conjg(AL1)*YeadjYe(i1,i2) - 4._dp*YeadjYeAYeadjAYe(i1,i2) &
  &  - 6._dp*g12*YeadjYeme2(i1,i2)/5._dp + 6._dp*g22*YeadjYeme2(i1,i2) -               &
  &  6._dp*TrYdadjYd*YeadjYeme2(i1,i2) - 2._dp*TrYeadjYe*YeadjYeme2(i1,i2) -               &
  &  6._dp*AbsL1sq*YeadjYeme2(i1,i2) - 4._dp*YeadjYeme2YeadjYe(i1,i2) -               &
  &  8._dp*mHd2*YeadjYeYeadjYe(i1,i2) - 2._dp*YeadjYeYeadjYeme2(i1,i2) - 4._dp*YeadjYeYeml2adjYe(i1,i2)&
  &  - 12._dp*YeadjYzAYzadjAYe(i1,i2) - 12._dp*YeadjYzmd2YzadjYe(i1,i2) - 12._dp*mHd2*YeadjYzYzadjYe(i1,i2)&
  &  - 12._dp*mz2*YeadjYzYzadjYe(i1,i2) - 6._dp*YeadjYzYzadjYeme2(i1,i2) - 12._dp*YeadjYzYzml2adjYe(i1,i2)&
  &  - 12._dp*YeCAYtAYtadjYe(i1,i2) - 12._dp*YeCYtAYtadjAYe(i1,i2) - 12._dp*YeCYtml2YtadjYe(i1,i2)&
  &  - 12._dp*mHd2*YeCYtYtadjYe(i1,i2) - 12._dp*mt2*YeCYtYtadjYe(i1,i2) - 6._dp*YeCYtYtadjYeme2(i1,i2)&
  &  - 12._dp*YeCYtYtml2adjYe(i1,i2) - 12._dp*g12*Yeml2adjYe(i1,i2)/5._dp +              &
  &  12._dp*g22*Yeml2adjYe(i1,i2) - 12._dp*TrYdadjYd*Yeml2adjYe(i1,i2) - 4._dp*TrYeadjYe*Yeml2adjYe(i1,i2)&
  &  - 12._dp*AbsL1sq*Yeml2adjYe(i1,i2) - 4._dp*Yeml2adjYeYeadjYe(i1,i2)              &
  &  - 12._dp*Yeml2adjYzYzadjYe(i1,i2) - 12._dp*Yeml2CYtYtadjYe(i1,i2)

  If (i1.eq.i2) Then
  betame22(i1,i2) = betame22(i1,i2)+(5328._dp*g14*MassB*Conjg(MassB)/25._dp + 24._dp*g14*Tr2(1)&
  & /5._dp + 8._dp*g12*Tr3(1))
  End If
   End Do
  End Do


  Dme2 = oo16pi2*( betame21 + oo16pi2 * betame22 )


  Else
  Dme2 = oo16pi2* betame21
  End If


  !--------------------
  ! mt2
  !--------------------

  betamt21 = TrAYtCAYt/2._dp + 3._dp*TrCAYtAYt/2._dp + TrCYtml2Yt + 3._dp*mt2*TrCYtYt/2._dp +&
  &  2._dp*TrCYtYtml2 + (mt2*TrYtCYt)/2._dp + TrYtCYtml2 + 4._dp*L1*mHd2*Conjg(L1)         &
  &  + 2._dp*L1*mt2*Conjg(L1) - 24._dp*g12*MassB*Conjg(MassB)/5._dp - 16._dp*g22*MassWB*Conjg(MassWB)&
  &  + 2._dp*AL1*Conjg(AL1) + 2._dp*g12*Tr1(1)


  If (TwoLoopRGE) Then
  betamt22 = -3._dp*g12*TrAYtCAYt/10._dp - (g22*TrAYtCAYt)/2._dp - 9._dp*g12*TrCAYtAYt/10._dp -&
  &  3._dp*g22*TrCAYtAYt/2._dp - 19._dp*TrCYtAYtCAYtYt/2._dp - 3._dp*g12*TrCYtml2Yt/5._dp -&
  &  g22*TrCYtml2Yt - TrCYtml2YtCYtYt/2._dp - 9._dp*g12*mt2*TrCYtYt/10._dp -           &
  &  3._dp*g22*mt2*TrCYtYt/2._dp - 4._dp*TrCYtYtadjAYeAYe - 12._dp*TrCYtYtadjAYzAYz -    &
  &  3._dp*TrCYtYtadjYeme2Ye - mt2*TrCYtYtadjYeYe - 9._dp*TrCYtYtadjYzmd2Yz -              &
  &  3._dp*mt2*TrCYtYtadjYzYz - 43._dp*TrCYtYtCAYtAYt/2._dp - 25._dp*TrCYtYtCYtml2Yt/2._dp -&
  &  19._dp*mt2*TrCYtYtCYtYt - 10._dp*TrCYtYtCYtYtml2 - 6._dp*g12*TrCYtYtml2/5._dp -     &
  &  2._dp*g22*TrCYtYtml2 - 8._dp*TrCYtYtml2CYtYt - 3._dp*Trmd2YzCYtYtadjYz -            &
  &  Trme2YeCYtYtadjYe - Trml2CYtYtCYtYt/2._dp - 3._dp*g12*Trml2YtCYt/5._dp -            &
  &  g22*Trml2YtCYt - 2._dp*Trml2YtCYtYtCYt - 4._dp*TrYeCAYtAYtadjYe - 4._dp*TrYeCYtAYtadjAYe -&
  &  4._dp*TrYeCYtml2YtadjYe - 4._dp*mHd2*TrYeCYtYtadjYe - 3._dp*mt2*TrYeCYtYtadjYe -      &
  &  4._dp*TrYeCYtYtml2adjYe - 4._dp*TrYeml2CYtYtadjYe - 4._dp*TrYtadjYeAYeCAYt -          &
  &  12._dp*TrYtadjYzAYzCAYt + 6._dp*g12*MassB*TrYtCAYt/5._dp + 2._dp*g22*MassWB*TrYtCAYt -&
  &  5._dp*TrYtCAYtAYtCYt/2._dp - 3._dp*g12*mt2*TrYtCYt/10._dp - (g22*mt2*TrYtCYt)     &
  & /2._dp - 29._dp*TrYtCYtAYtCAYt/2._dp - TrYtCYtml2YtCYt - 5._dp*mt2*TrYtCYtYtCYt -      &
  &  7._dp*TrYtCYtYtCYtml2/2._dp - 7._dp*TrYtCYtYtml2CYt/2._dp - 13._dp*TrYtml2CYtYtCYt/2._dp -&
  &  12._dp*TrYzCAYtAYtadjYz - 12._dp*TrYzCYtAYtadjAYz - 12._dp*TrYzCYtml2YtadjYz -        &
  &  9._dp*mt2*TrYzCYtYtadjYz - 12._dp*mz2*TrYzCYtYtadjYz - 12._dp*TrYzCYtYtml2adjYz -     &
  &  12._dp*TrYzml2CYtYtadjYz - 12._dp*g12*L1*mHd2*Conjg(L1)/5._dp - 4._dp*g22*L1*mHd2*Conjg(L1)&
  &  - 6._dp*g12*L1*mt2*Conjg(L1)/5._dp - 2._dp*g22*L1*mt2*Conjg(L1) - 12._dp*L1*TrAYdadjAYd*Conjg(L1)&
  &  - 4._dp*L1*TrAYeadjAYe*Conjg(L1) - 12._dp*L1*Trmd2YdadjYd*Conjg(L1) - 4._dp*L1*Trme2YeadjYe*Conjg(L1)&
  &  - 36._dp*L1*mHd2*TrYdadjYd*Conjg(L1) - 12._dp*L1*mt2*TrYdadjYd*Conjg(L1)              &
  &  - 12._dp*L1*TrYdmq2adjYd*Conjg(L1) - 12._dp*L1*mHd2*TrYeadjYe*Conjg(L1)               &
  &  - 4._dp*L1*mt2*TrYeadjYe*Conjg(L1) - 4._dp*L1*TrYeml2adjYe*Conjg(L1) - 12._dp*TrYdadjAYd*AL1*Conjg(L1)&
  &  - 4._dp*TrYeadjAYe*AL1*Conjg(L1) - 48._dp*L1**2*mHd2*Conjg(L1)**2 - 24._dp*L1**2*mt2*Conjg(L1)&
  & **2 + 5328._dp*g14*MassB*Conjg(MassB)/25._dp + 192._dp*g12*g22*MassB*Conjg(MassB)&
  & /5._dp + 96._dp*g12*g22*MassWB*Conjg(MassB)/5._dp + 6._dp*g12*TrCYtAYt*Conjg(MassB)&
  & /5._dp - 9._dp*g12*MassB*TrCYtYt*Conjg(MassB)/5._dp - 3._dp*g12*MassB*TrYtCYt*Conjg(MassB)&
  & /5._dp - 12._dp*g12*L1*MassB*Conjg(L1)*Conjg(MassB)/5._dp + 6._dp*g12*AL1*Conjg(L1)&
  & *Conjg(MassB)/5._dp + 96._dp*g12*g22*MassB*Conjg(MassWB)/5._dp + 192._dp*g12*g22*MassWB*Conjg(MassWB)&
  & /5._dp + 544._dp*g24*MassWB*Conjg(MassWB) + 2._dp*g22*TrCYtAYt*Conjg(MassWB)       &
  &  - 3._dp*g22*MassWB*TrCYtYt*Conjg(MassWB) - g22*MassWB*TrYtCYt*Conjg(MassWB)       &
  &  - 4._dp*g22*L1*MassWB*Conjg(L1)*Conjg(MassWB) + 2._dp*g22*AL1*Conjg(L1)           &
  & *Conjg(MassWB) + 6._dp*g12*L1*MassB*Conjg(AL1)/5._dp + 2._dp*g22*L1*MassWB*Conjg(AL1)&
  &  - 12._dp*L1*TrAYdadjYd*Conjg(AL1) - 4._dp*L1*TrAYeadjYe*Conjg(AL1) - 6._dp*g12*AL1*Conjg(AL1)&
  & /5._dp - 2._dp*g22*AL1*Conjg(AL1) - 12._dp*TrYdadjYd*AL1*Conjg(AL1) - 4._dp*TrYeadjYe*AL1*Conjg(AL1)&
  &  - 48._dp*L1*AL1*Conjg(L1)*Conjg(AL1) + 24._dp*g14*Tr2(1)/5._dp + 16._dp*g24*Tr2(2)&
  &  + 8._dp*g12*Tr3(1)


  Dmt2 = oo16pi2*( betamt21 + oo16pi2 * betamt22 )


  Else
  Dmt2 = oo16pi2* betamt21
  End If


  !--------------------
  ! mtb2
  !--------------------

  betamtb21 = 4._dp*L2*mHu2*Conjg(L2) + 2._dp*L2*mtb2*Conjg(L2) - 24._dp*g12*MassB*Conjg(MassB)&
  & /5._dp - 16._dp*g22*MassWB*Conjg(MassWB) + 2._dp*AL2*Conjg(AL2) - 2._dp*g12*Tr1(1)


  If (TwoLoopRGE) Then
  betamtb22 = -12._dp*g12*L2*mHu2*Conjg(L2)/5._dp - 4._dp*g22*L2*mHu2*Conjg(L2)     &
  &  - 6._dp*g12*L2*mtb2*Conjg(L2)/5._dp - 2._dp*g22*L2*mtb2*Conjg(L2) -               &
  &  12._dp*L2*TrAYuadjAYu*Conjg(L2) - 12._dp*L2*Trmu2YuadjYu*Conjg(L2) - 36._dp*L2*mHu2*TrYuadjYu*Conjg(L2)&
  &  - 12._dp*L2*mtb2*TrYuadjYu*Conjg(L2) - 12._dp*L2*TrYumq2adjYu*Conjg(L2)               &
  &  - 12._dp*TrYuadjAYu*AL2*Conjg(L2) - 48._dp*L2**2*mHu2*Conjg(L2)**2 - 24._dp*L2**2*mtb2*Conjg(L2)&
  & **2 + 5328._dp*g14*MassB*Conjg(MassB)/25._dp + 192._dp*g12*g22*MassB*Conjg(MassB)&
  & /5._dp + 96._dp*g12*g22*MassWB*Conjg(MassB)/5._dp - 12._dp*g12*L2*MassB*Conjg(L2)&
  & *Conjg(MassB)/5._dp + 6._dp*g12*AL2*Conjg(L2)*Conjg(MassB)/5._dp + 96._dp*g12*g22*MassB*Conjg(MassWB)&
  & /5._dp + 192._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 544._dp*g24*MassWB*Conjg(MassWB)&
  &  - 4._dp*g22*L2*MassWB*Conjg(L2)*Conjg(MassWB) + 2._dp*g22*AL2*Conjg(L2)           &
  & *Conjg(MassWB) + 6._dp*g12*L2*MassB*Conjg(AL2)/5._dp + 2._dp*g22*L2*MassWB*Conjg(AL2)&
  &  - 12._dp*L2*TrAYuadjYu*Conjg(AL2) - 6._dp*g12*AL2*Conjg(AL2)/5._dp - 2._dp*g22*AL2*Conjg(AL2)&
  &  - 12._dp*TrYuadjYu*AL2*Conjg(AL2) - 48._dp*L2*AL2*Conjg(L2)*Conjg(AL2) +              &
  &  24._dp*g14*Tr2(1)/5._dp + 16._dp*g24*Tr2(2) - 8._dp*g12*Tr3(1)


  Dmtb2 = oo16pi2*( betamtb21 + oo16pi2 * betamtb22 )


  Else
  Dmtb2 = oo16pi2* betamtb21
  End If


  !--------------------
  ! ms2
  !--------------------

  betams21 = TrAYsCAYs/2._dp + 3._dp*TrCAYsAYs/2._dp + TrCYsmd2Ys + 3._dp*ms2*TrCYsYs/2._dp +&
  &  2._dp*TrCYsYsmd2 + (ms2*TrYsCYs)/2._dp + TrYsCYsmd2 - 32._dp*g12*MassB*Conjg(MassB) &
  & /15._dp - 80._dp*g32*MassG*Conjg(MassG)/3._dp - 4._dp*g12*Tr1(1)/3._dp


  If (TwoLoopRGE) Then
  betams22 = -2._dp*g12*TrAYsCAYs/15._dp - 2._dp*g32*TrAYsCAYs/3._dp -              &
  &  2._dp*g12*TrCAYsAYs/5._dp - 2._dp*g32*TrCAYsAYs - 13._dp*TrCYsAYsCAYsYs -         &
  &  4._dp*g12*TrCYsmd2Ys/15._dp - 4._dp*g32*TrCYsmd2Ys/3._dp - TrCYsmd2YsCYsYs/2._dp -&
  &  6._dp*TrCYsYdadjAYdAYs - 2._dp*TrCYsYdadjYdmd2Ys - 2._dp*g12*ms2*TrCYsYs/5._dp -    &
  &  2._dp*g32*ms2*TrCYsYs - 29._dp*TrCYsYsCAYsAYs - 35._dp*TrCYsYsCYsmd2Ys/2._dp -      &
  &  26._dp*ms2*TrCYsYsCYsYs - 27._dp*TrCYsYsCYsYsmd2/2._dp - 8._dp*g12*TrCYsYsmd2/15._dp -&
  &  8._dp*g32*TrCYsYsmd2/3._dp - 11._dp*TrCYsYsmd2CYsYs - 6._dp*TrCYsYzadjAYzAYs -      &
  &  2._dp*TrCYsYzadjYzmd2Ys - Trmd2CYsYsCYsYs/2._dp - 4._dp*g12*Trmd2YsCYs/15._dp -     &
  &  4._dp*g32*Trmd2YsCYs/3._dp - 6._dp*Trmd2YsCYsYdadjYd - 5._dp*Trmd2YsCYsYsCYs/2._dp -&
  &  6._dp*Trmd2YsCYsYzadjYz - 2._dp*Trml2adjYzYsCYsYz - 2._dp*Trmq2adjYdYsCYsYd -         &
  &  2._dp*TrYdadjAYdAYsCYs - 8._dp*TrYdadjYdAYsCAYs - ms2*TrYdadjYdYsCYs + 8._dp*g12*MassB*TrYsCAYs/15._dp +&
  &  8._dp*g32*MassG*TrYsCAYs/3._dp - 8._dp*TrYsCAYsAYdadjYd - 3._dp*TrYsCAYsAYsCYs -    &
  &  8._dp*TrYsCAYsAYzadjYz - 2._dp*g12*ms2*TrYsCYs/15._dp - 2._dp*g32*ms2*TrYsCYs/3._dp -&
  &  8._dp*TrYsCYsAYdadjAYd - 19._dp*TrYsCYsAYsCAYs - 8._dp*TrYsCYsAYzadjAYz -             &
  &  8._dp*TrYsCYsmd2YdadjYd - 3._dp*TrYsCYsmd2YsCYs/2._dp - 8._dp*TrYsCYsmd2YzadjYz -     &
  &  8._dp*mHd2*TrYsCYsYdadjYd - 7._dp*ms2*TrYsCYsYdadjYd - 6._dp*TrYsCYsYdmq2adjYd -      &
  &  6._dp*ms2*TrYsCYsYsCYs - 4._dp*TrYsCYsYsCYsmd2 - 5._dp*TrYsCYsYsmd2CYs -              &
  &  7._dp*ms2*TrYsCYsYzadjYz - 8._dp*mz2*TrYsCYsYzadjYz - 6._dp*TrYsCYsYzml2adjYz -       &
  &  8._dp*TrYsmd2CYsYdadjYd - 8._dp*TrYsmd2CYsYsCYs - 8._dp*TrYsmd2CYsYzadjYz -           &
  &  2._dp*TrYzadjAYzAYsCYs - 8._dp*TrYzadjYzAYsCAYs - ms2*TrYzadjYzYsCYs + 6784._dp*g14*MassB*Conjg(MassB)&
  & /75._dp + 256._dp*g12*g32*MassB*Conjg(MassB)/9._dp + 128._dp*g12*g32*MassG*Conjg(MassB)&
  & /9._dp + 8._dp*g12*TrCYsAYs*Conjg(MassB)/15._dp - 4._dp*g12*MassB*TrCYsYs*Conjg(MassB)&
  & /5._dp - 4._dp*g12*MassB*TrYsCYs*Conjg(MassB)/15._dp + 128._dp*g12*g32*MassB*Conjg(MassG)&
  & /9._dp + 256._dp*g12*g32*MassG*Conjg(MassG)/9._dp + 2320._dp*g34*MassG*Conjg(MassG)&
  & /3._dp + 8._dp*g32*TrCYsAYs*Conjg(MassG)/3._dp - 4._dp*g32*MassG*TrCYsYs*Conjg(MassG)&
  &  - 4._dp*g32*MassG*TrYsCYs*Conjg(MassG)/3._dp + 32._dp*g14*Tr2(1)/15._dp +         &
  &  80._dp*g34*Tr2(3)/3._dp - 16._dp*g12*Tr3(1)/3._dp + 4._dp*g32*Tr3(3)


  Dms2 = oo16pi2*( betams21 + oo16pi2 * betams22 )


  Else
  Dms2 = oo16pi2* betams21
  End If


  !--------------------
  ! msb2
  !--------------------

  betamsb21 = -32._dp*g12*MassB*Conjg(MassB)/15._dp - 80._dp*g32*MassG*Conjg(MassG) &
  & /3._dp + 4._dp*g12*Tr1(1)/3._dp


  If (TwoLoopRGE) Then
  betamsb22 = 6784._dp*g14*MassB*Conjg(MassB)/75._dp + 256._dp*g12*g32*MassB*Conjg(MassB)&
  & /9._dp + 128._dp*g12*g32*MassG*Conjg(MassB)/9._dp + 128._dp*g12*g32*MassB*Conjg(MassG)&
  & /9._dp + 256._dp*g12*g32*MassG*Conjg(MassG)/9._dp + 2320._dp*g34*MassG*Conjg(MassG)&
  & /3._dp + 32._dp*g14*Tr2(1)/15._dp + 80._dp*g34*Tr2(3)/3._dp + 16._dp*g12*Tr3(1)  &
  & /3._dp - 4._dp*g32*Tr3(3)


  Dmsb2 = oo16pi2*( betamsb21 + oo16pi2 * betamsb22 )


  Else
  Dmsb2 = oo16pi2* betamsb21
  End If


  !--------------------
  ! mz2
  !--------------------

  betamz21 = 2._dp*TrAYzadjAYz + 2._dp*Trmd2YzadjYz + 2._dp*mz2*TrYzadjYz +             &
  &  2._dp*TrYzml2adjYz - 2._dp*g12*MassB*Conjg(MassB)/15._dp - 32._dp*g32*MassG*Conjg(MassG)&
  & /3._dp - 6._dp*g22*MassWB*Conjg(MassWB) + (g12*Tr1(1))/3._dp


  If (TwoLoopRGE) Then
  betamz22 = 4._dp*g12*TrAYzadjAYz/5._dp - 4._dp*TrCYsAYzadjAYzYs - 5._dp*TrCYsYzadjAYzAYs -&
  &  4._dp*TrCYsYzadjYzmd2Ys - 4._dp*ms2*TrCYsYzadjYzYs - 7._dp*mz2*TrCYsYzadjYzYs/2._dp - &
  &  4._dp*TrCYsYzadjYzYsmd2 - 4._dp*TrCYsYzml2adjYzYs - 2._dp*TrCYtAYtadjAYzYz -          &
  &  3._dp*TrCYtYtadjAYzAYz - 3._dp*TrCYtYtadjYzmd2Yz - 2._dp*mt2*TrCYtYtadjYzYz -         &
  &  3._dp*mz2*TrCYtYtadjYzYz/2._dp - 3._dp*TrCYtYtadjYzYzml2 - 2._dp*TrCYtYtml2adjYzYz +  &
  &  4._dp*g12*Trmd2YzadjYz/5._dp - 3._dp*Trmd2YzadjYzYsCYs - Trmd2YzCYtYtadjYz -        &
  &  Trml2adjYzYsCYsYz - 4._dp*TrYdadjAYdAYzadjYz - 4._dp*TrYdadjYdAYzadjAYz -             &
  &  4._dp*TrYdadjYdmd2YzadjYz - 2._dp*TrYeadjAYzAYzadjYe - 2._dp*TrYeadjYzAYzadjAYe -     &
  &  2._dp*TrYeadjYzmd2YzadjYe - 2._dp*mHd2*TrYeadjYzYzadjYe - mz2*TrYeadjYzYzadjYe -      &
  &  2._dp*TrYeadjYzYzml2adjYe - 2._dp*TrYeml2adjYzYzadjYe - 8._dp*TrYsCAYsAYzadjYz -      &
  &  4._dp*TrYsCYsAYzadjAYz - 5._dp*TrYsCYsmd2YzadjYz - ms2*TrYsCYsYzadjYz -               &
  &  3._dp*mz2*TrYsCYsYzadjYz/2._dp - TrYsCYsYzadjYzmd2 - 4._dp*TrYsmd2CYsYzadjYz -        &
  &  3._dp*TrYtadjAYzAYzCYt - 6._dp*TrYtadjYzAYzCAYt - 2._dp*TrYtadjYzmd2YzCYt -           &
  &  3._dp*mt2*TrYtadjYzYzCYt - mz2*TrYtadjYzYzCYt - 3._dp*TrYtadjYzYzCYtml2 -             &
  &  2._dp*TrYtadjYzYzml2CYt - 3._dp*TrYtml2adjYzYzCYt - 2._dp*TrYzadjAYeAYeadjYz -        &
  &  4._dp*g12*MassB*TrYzadjAYz/5._dp - 4._dp*TrYzadjAYzAYdadjYd - 3._dp*TrYzadjAYzAYsCYs -&
  &  20._dp*TrYzadjAYzAYzadjYz - 2._dp*TrYzadjYeAYeadjAYz - 2._dp*TrYzadjYeme2YeadjYz -    &
  &  mz2*TrYzadjYeYeadjYz + 4._dp*g12*mz2*TrYzadjYz/5._dp - 4._dp*TrYzadjYzAYdadjAYd -   &
  &  8._dp*TrYzadjYzAYsCAYs - 20._dp*TrYzadjYzAYzadjAYz - 4._dp*TrYzadjYzmd2YdadjYd -      &
  &  3._dp*TrYzadjYzmd2YsCYs - 20._dp*TrYzadjYzmd2YzadjYz - 4._dp*mHd2*TrYzadjYzYdadjYd -  &
  &  4._dp*mz2*TrYzadjYzYdadjYd - 4._dp*TrYzadjYzYdmq2adjYd - 3._dp*ms2*TrYzadjYzYsCYs -   &
  &  3._dp*mz2*TrYzadjYzYsCYs - 20._dp*mz2*TrYzadjYzYzadjYz - 10._dp*TrYzadjYzYzml2adjYz - &
  &  6._dp*TrYzCAYtAYtadjYz - 4._dp*TrYzCYtAYtadjAYz - 3._dp*TrYzCYtml2YtadjYz -           &
  &  mt2*TrYzCYtYtadjYz - 7._dp*mz2*TrYzCYtYtadjYz/2._dp - TrYzCYtYtml2adjYz +             &
  &  4._dp*g12*TrYzml2adjYz/5._dp - 4._dp*TrYzml2adjYzYdadjYd - 3._dp*TrYzml2adjYzYsCYs -&
  &  10._dp*TrYzml2adjYzYzadjYz - TrYzml2CYtYtadjYz + 409._dp*g14*MassB*Conjg(MassB)     &
  & /75._dp + 2._dp*g12*g22*MassB*Conjg(MassB)/5._dp + 32._dp*g12*g32*MassB*Conjg(MassB)&
  & /45._dp + 16._dp*g12*g32*MassG*Conjg(MassB)/45._dp + (g12*g22*MassWB*Conjg(MassB))&
  & /5._dp - 4._dp*g12*TrAYzadjYz*Conjg(MassB)/5._dp + 8._dp*g12*MassB*TrYzadjYz*Conjg(MassB)&
  & /5._dp + 16._dp*g12*g32*MassB*Conjg(MassG)/45._dp + 32._dp*g12*g32*MassG*Conjg(MassG)&
  & /45._dp + 32._dp*g22*g32*MassG*Conjg(MassG) + 544._dp*g34*MassG*Conjg(MassG)     &
  & /3._dp + 16._dp*g22*g32*MassWB*Conjg(MassG) + (g12*g22*MassB*Conjg(MassWB))    &
  & /5._dp + 16._dp*g22*g32*MassG*Conjg(MassWB) + 2._dp*g12*g22*MassWB*Conjg(MassWB)&
  & /5._dp + 159._dp*g24*MassWB*Conjg(MassWB) + 32._dp*g22*g32*MassWB*Conjg(MassWB)  &
  &  + 2._dp*g14*Tr2(1)/15._dp + 6._dp*g24*Tr2(2) + 32._dp*g34*Tr2(3)/3._dp +        &
  &  4._dp*g12*Tr3(1)/3._dp + 4._dp*g32*Tr3(3) + 4._dp*(g32*Tr3(3))/sqrt3


  Dmz2 = oo16pi2*( betamz21 + oo16pi2 * betamz22 )


  Else
  Dmz2 = oo16pi2* betamz21
  End If


  !--------------------
  ! mzb2
  !--------------------

  betamzb21 = -2._dp*g12*MassB*Conjg(MassB)/15._dp - 32._dp*g32*MassG*Conjg(MassG)  &
  & /3._dp - 6._dp*g22*MassWB*Conjg(MassWB) - (g12*Tr1(1))/3._dp


  If (TwoLoopRGE) Then
  betamzb22 = 409._dp*g14*MassB*Conjg(MassB)/75._dp + 2._dp*g12*g22*MassB*Conjg(MassB)&
  & /5._dp + 32._dp*g12*g32*MassB*Conjg(MassB)/45._dp + 16._dp*g12*g32*MassG*Conjg(MassB)&
  & /45._dp + (g12*g22*MassWB*Conjg(MassB))/5._dp + 16._dp*g12*g32*MassB*Conjg(MassG)&
  & /45._dp + 32._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 32._dp*g22*g32*MassG*Conjg(MassG)&
  &  + 544._dp*g34*MassG*Conjg(MassG)/3._dp + 16._dp*g22*g32*MassWB*Conjg(MassG)     &
  &  + (g12*g22*MassB*Conjg(MassWB))/5._dp + 16._dp*g22*g32*MassG*Conjg(MassWB)    &
  &  + 2._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g24*MassWB*Conjg(MassWB)   &
  &  + 32._dp*g22*g32*MassWB*Conjg(MassWB) + 2._dp*g14*Tr2(1)/15._dp +               &
  &  6._dp*g24*Tr2(2) + 32._dp*g34*Tr2(3)/3._dp - 4._dp*g12*Tr3(1)/3._dp -           &
  &  4._dp*g32*Tr3(3) - 4._dp*(g32*Tr3(3))/sqrt3


  Dmzb2 = oo16pi2*( betamzb21 + oo16pi2 * betamzb22 )


  Else
  Dmzb2 = oo16pi2* betamzb21
  End If


  !--------------------
  ! MassB
  !--------------------

  betaMassB1 = 136._dp*MassB/5._dp


  If (TwoLoopRGE) Then
  betaMassB2 = 6008._dp*g12*MassB/75._dp + 348._dp*g22*MassB/5._dp +          &
  &  368._dp*g32*MassB/3._dp + 368._dp*g32*MassG/3._dp + 348._dp*g22*MassWB/5._dp +&
  &  28._dp*TrAYdadjYd/5._dp + 36._dp*TrAYeadjYe/5._dp + 52._dp*TrAYuadjYu/5._dp +&
  &  28._dp*TrAYzadjYz/5._dp + 48._dp*TrCYsAYs/5._dp - 42._dp*MassB*TrCYsYs/5._dp +&
  &  54._dp*TrCYtAYt/5._dp - 9._dp*MassB*TrCYtYt - 28._dp*MassB*TrYdadjYd/5._dp -&
  &  36._dp*MassB*TrYeadjYe/5._dp - 6._dp*MassB*TrYsCYs/5._dp - 9._dp*MassB*TrYtCYt/5._dp -&
  &  52._dp*MassB*TrYuadjYu/5._dp - 28._dp*MassB*TrYzadjYz/5._dp -             &
  &  54._dp*L1*MassB*Conjg(L1)/5._dp + 54._dp*AL1*Conjg(L1)/5._dp -            &
  &  54._dp*L2*MassB*Conjg(L2)/5._dp + 54._dp*AL2*Conjg(L2)/5._dp


  DMassB = oo16pi2*g12*( betaMassB1 + oo16pi2 * betaMassB2 )


  Else
  DMassB = oo16pi2*g12* betaMassB1
  End If


  !--------------------
  ! MassWB
  !--------------------

  betaMassWB1 = 16._dp*MassWB

  If (TwoLoopRGE) Then
  betaMassWB2 = 116._dp*g12*MassB/5._dp + 80._dp*g32*MassG +            &
  &  116._dp*g12*MassWB/5._dp + 376._dp*g22*MassWB + 80._dp*g32*MassWB + &
  &  12._dp*TrAYdadjYd + 4._dp*TrAYeadjYe + 12._dp*TrAYuadjYu +          &
  &  12._dp*TrAYzadjYz + 14._dp*TrCYtAYt - 35._dp*MassWB*TrCYtYt/3._dp - &
  &  12._dp*MassWB*TrYdadjYd - 4._dp*MassWB*TrYeadjYe - 7._dp*MassWB*TrYtCYt/3._dp -&
  &  12._dp*MassWB*TrYuadjYu - 12._dp*MassWB*TrYzadjYz - 14._dp*L1*MassWB*Conjg(L1)&
  &  + 14._dp*AL1*Conjg(L1) - 14._dp*L2*MassWB*Conjg(L2) + 14._dp*AL2*Conjg(L2)


  DMassWB = oo16pi2*g22*( betaMassWB1 + oo16pi2 * betaMassWB2 )


  Else
  DMassWB = oo16pi2*g22* betaMassWB1
  End If


  !--------------------
  ! MassG
  !--------------------

  betaMassG1 = 8._dp*MassG


  If (TwoLoopRGE) Then
  betaMassG2 = 46._dp*g12*MassB/3._dp + 46._dp*g12*MassG/3._dp +        &
  &  30._dp*g22*MassG + 1600._dp*g32*MassG/3._dp + 30._dp*g22*MassWB +   &
  &  8._dp*TrAYdadjYd + 8._dp*TrAYuadjYu + 8._dp*TrAYzadjYz +            &
  &  18._dp*TrCYsAYs - 63._dp*MassG*TrCYsYs/4._dp - 8._dp*MassG*TrYdadjYd -&
  &  9._dp*MassG*TrYsCYs/4._dp - 8._dp*MassG*TrYuadjYu - 8._dp*MassG*TrYzadjYz


  DMassG = oo16pi2*g32*( betaMassG1 + oo16pi2 * betaMassG2 )


  Else
  DMassG = oo16pi2*g32* betaMassG1
  End If


  Call ParametersToG353(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYt,DYs,DYz,DL1,DL2,DMTM,               &
  & Dmue,DMZM,DMSM,DAYu,DAYd,DAYe,DAYt,DAYs,DAYz,DAL1,DAL2,DAmue,DAMTM,DAMZM,              &
  & DAMSM,Dmq2,Dml2,DmHd2,DmHu2,Dmd2,Dmu2,Dme2,Dmt2,Dmtb2,Dms2,Dmsb2,Dmz2,Dmzb2,           &
  & DMassB,DMassWB,DMassG,f)

  Iname = Iname - 1

 End Subroutine rge353


 Subroutine GToParameters117(g,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM)

 Implicit None 
  Real(dp), Intent(in) :: g(117) 
  Real(dp),Intent(out) :: g1,g2,g3

  Complex(dp),Intent(out) :: Yu(3,3),Yd(3,3),Ye(3,3),Yt(3,3),Ys(3,3),Yz(3,3)     & 
    & ,L1,L2,MTM

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters117' 
 
  g1= g(1) 
  g2= g(2) 
  g3= g(3) 
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3)
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    Yt(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    Ys(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    Yz(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
   End Do 
  End Do 
 
  L1= Cmplx(g(112),g(113),dp) 
  L2= Cmplx(g(114),g(115),dp) 
  MTM= Cmplx(g(116),g(117),dp) 
  Iname = Iname - 1 
 
 End Subroutine GToParameters117


 Subroutine ParametersToG117(g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM,g)

 Implicit None 
  Real(dp), Intent(out) :: g(117) 
  Real(dp), Intent(in) :: g1,g2,g3

  Complex(dp), Intent(in) :: Yu(3,3),Yd(3,3),Ye(3,3),Yt(3,3),Ys(3,3),Yz(3,3)            & 
   & ,L1,L2,MTM

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG117' 
 
  g(1) = g1  
  g(2) = g2  
  g(3) = g3  
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(Yt(i1,i2), dp) 
    g(SumI+59) = Aimag(Yt(i1,i2)) 
    g(SumI+76) = Real(Ys(i1,i2), dp) 
    g(SumI+77) = Aimag(Ys(i1,i2)) 
    g(SumI+94) = Real(Yz(i1,i2), dp) 
    g(SumI+95) = Aimag(Yz(i1,i2)) 
   End Do 
  End Do 

  g(112) = Real(L1,dp)  
  g(113) = Aimag(L1)  
  g(114) = Real(L2,dp)  
  g(115) = Aimag(L2)  
  g(116) = Real(MTM,dp)  
  g(117) = Aimag(MTM)  

  Iname = Iname - 1 
 
 End Subroutine ParametersToG117


 Subroutine GToParameters555(g,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,mue,MXM3,MWM3,     &
 & MGM3,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,Bmue,BMXM3,BMWM3,BMGM3,BMBM3,mq2,ml2,   &
 & mHd2,mHu2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MassB,MassWB,MassG,MnuL)

 Implicit None 
  Real(dp), Intent(in) :: g(573) 
  Real(dp),Intent(out) :: g1,g2,g3,mHd2,mHu2

  Complex(dp),Intent(out) :: mue, Bmue, MassB, MassWB, MassG
 
  Complex(dp),Intent(out), Dimension(3,3) :: Yu,Yd,Ye,Yb3,Yw3,Yx3              &
     & ,MXM3,MWM3,MGM3,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,BMXM3,BMWM3,BMGM3,BMBM3 &
     & ,mq2,ml2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MnuL
 
  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'GToParameters555' 
 
  g1= g(1) 
  g2= g(2) 
  g3= g(3) 
  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
    Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp)
    Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
    Yb3(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
    Yw3(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
    Yx3(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
    MXM3(i1,i2) = Cmplx( g(SumI+114), g(SumI+115), dp) 
    MWM3(i1,i2) = Cmplx( g(SumI+132), g(SumI+133), dp) 
    MGM3(i1,i2) = Cmplx( g(SumI+150), g(SumI+151), dp) 
    MBM3(i1,i2) = Cmplx( g(SumI+168), g(SumI+169), dp) 
    TYu(i1,i2) = Cmplx( g(SumI+186), g(SumI+187), dp) 
    TYd(i1,i2) = Cmplx( g(SumI+204), g(SumI+205), dp) 
    TYe(i1,i2) = Cmplx( g(SumI+222), g(SumI+223), dp) 
    TYb3(i1,i2) = Cmplx( g(SumI+240), g(SumI+241), dp) 
    TYw3(i1,i2) = Cmplx( g(SumI+258), g(SumI+259), dp) 
    TYx3(i1,i2) = Cmplx( g(SumI+276), g(SumI+277), dp) 
    BMXM3(i1,i2) = Cmplx( g(SumI+296), g(SumI+297), dp) 
    BMWM3(i1,i2) = Cmplx( g(SumI+314), g(SumI+315), dp) 
    BMGM3(i1,i2) = Cmplx( g(SumI+332), g(SumI+333), dp) 
    BMBM3(i1,i2) = Cmplx( g(SumI+350), g(SumI+351), dp) 
    mq2(i1,i2) = Cmplx( g(SumI+368), g(SumI+369), dp) 
    ml2(i1,i2) = Cmplx( g(SumI+386), g(SumI+387), dp) 
    md2(i1,i2) = Cmplx( g(SumI+406), g(SumI+407), dp) 
    mu2(i1,i2) = Cmplx( g(SumI+424), g(SumI+425), dp) 
    me2(i1,i2) = Cmplx( g(SumI+442), g(SumI+443), dp) 
    mHw32(i1,i2) = Cmplx( g(SumI+460), g(SumI+461), dp) 
    mHg3p2(i1,i2) = Cmplx( g(SumI+478), g(SumI+479), dp) 
    mHb32(i1,i2) = Cmplx( g(SumI+496), g(SumI+497), dp) 
    mHx32(i1,i2) = Cmplx( g(SumI+514), g(SumI+515), dp) 
    mHxb32(i1,i2) = Cmplx( g(SumI+532), g(SumI+533), dp) 
    MnuL(i1,i2) = Cmplx( g(SumI+556), g(SumI+557), dp) 
   End Do 
  End Do 

  mue= Cmplx(g(112),g(113),dp) 
  Bmue= Cmplx(g(294),g(295),dp) 
  mHd2= g(404) 
  mHu2= g(405) 
  MassB= Cmplx(g(550),g(551),dp) 
  MassWB= Cmplx(g(552),g(553),dp) 
  MassG= Cmplx(g(554),g(555),dp) 

  Iname = Iname - 1 
 
 End Subroutine GToParameters555


 Subroutine ParametersToG555(g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,mue,MXM3,MWM3,MGM3  &
  & ,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,Bmue,BMXM3,BMWM3,BMGM3,BMBM3,mq2,ml2,mHd2 &
  & ,mHu2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MassB,MassWB,MassG,MnuL,g)

 Implicit None 
  Real(dp), Intent(out) :: g(573) 
  Real(dp), Intent(in) :: g1,g2,g3,mHd2,mHu2

  Complex(dp),Intent(in) :: mue, Bmue, MassB, MassWB, MassG
 
  Complex(dp),Intent(in), Dimension(3,3) :: Yu,Yd,Ye,Yb3,Yw3,Yx3               &
     & ,MXM3,MWM3,MGM3,MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,BMXM3,BMWM3,BMGM3,BMBM3 &
     & ,mq2,ml2,md2,mu2,me2,mHw32,mHg3p2,mHb32,mHx32,mHxb32,MnuL

  Integer i1, i2, SumI 
 
  Iname = Iname +1 
  NameOfUnit(Iname) = 'ParametersToG555' 
 
  g(1) = g1  
  g(2) = g2  
  g(3) = g3  

  Do i1 = 1,3
   Do i2 = 1,3
    SumI = 2 * ( (i2-1) + (i1-1)*3 )
    g(SumI+4) = Real(Yu(i1,i2), dp) 
    g(SumI+5) = Aimag(Yu(i1,i2)) 
    g(SumI+22) = Real(Yd(i1,i2), dp) 
    g(SumI+23) = Aimag(Yd(i1,i2)) 
    g(SumI+40) = Real(Ye(i1,i2), dp) 
    g(SumI+41) = Aimag(Ye(i1,i2)) 
    g(SumI+58) = Real(Yb3(i1,i2), dp) 
    g(SumI+59) = Aimag(Yb3(i1,i2)) 
    g(SumI+76) = Real(Yw3(i1,i2), dp) 
    g(SumI+77) = Aimag(Yw3(i1,i2)) 
    g(SumI+94) = Real(Yx3(i1,i2), dp) 
    g(SumI+95) = Aimag(Yx3(i1,i2)) 
    g(SumI+114) = Real(MXM3(i1,i2), dp) 
    g(SumI+115) = Aimag(MXM3(i1,i2)) 
    g(SumI+132) = Real(MWM3(i1,i2), dp) 
    g(SumI+133) = Aimag(MWM3(i1,i2)) 
    g(SumI+150) = Real(MGM3(i1,i2), dp) 
    g(SumI+151) = Aimag(MGM3(i1,i2)) 
    g(SumI+168) = Real(MBM3(i1,i2), dp) 
    g(SumI+169) = Aimag(MBM3(i1,i2)) 
    g(SumI+186) = Real(TYu(i1,i2), dp) 
    g(SumI+187) = Aimag(TYu(i1,i2)) 
    g(SumI+204) = Real(TYd(i1,i2), dp) 
    g(SumI+205) = Aimag(TYd(i1,i2)) 
    g(SumI+222) = Real(TYe(i1,i2), dp) 
    g(SumI+223) = Aimag(TYe(i1,i2)) 
    g(SumI+240) = Real(TYb3(i1,i2), dp) 
    g(SumI+241) = Aimag(TYb3(i1,i2)) 
    g(SumI+258) = Real(TYw3(i1,i2), dp) 
    g(SumI+259) = Aimag(TYw3(i1,i2)) 
    g(SumI+276) = Real(TYx3(i1,i2), dp) 
    g(SumI+277) = Aimag(TYx3(i1,i2)) 
    g(SumI+296) = Real(BMXM3(i1,i2), dp) 
    g(SumI+297) = Aimag(BMXM3(i1,i2)) 
    g(SumI+314) = Real(BMWM3(i1,i2), dp) 
    g(SumI+315) = Aimag(BMWM3(i1,i2)) 
    g(SumI+332) = Real(BMGM3(i1,i2), dp) 
    g(SumI+333) = Aimag(BMGM3(i1,i2)) 
    g(SumI+350) = Real(BMBM3(i1,i2), dp) 
    g(SumI+351) = Aimag(BMBM3(i1,i2)) 
    g(SumI+368) = Real(mq2(i1,i2), dp) 
    g(SumI+369) = Aimag(mq2(i1,i2)) 
    g(SumI+386) = Real(ml2(i1,i2), dp) 
    g(SumI+387) = Aimag(ml2(i1,i2)) 
    g(SumI+406) = Real(md2(i1,i2), dp) 
    g(SumI+407) = Aimag(md2(i1,i2)) 
    g(SumI+424) = Real(mu2(i1,i2), dp) 
    g(SumI+425) = Aimag(mu2(i1,i2)) 
    g(SumI+442) = Real(me2(i1,i2), dp) 
    g(SumI+443) = Aimag(me2(i1,i2)) 
    g(SumI+460) = Real(mHw32(i1,i2), dp) 
    g(SumI+461) = Aimag(mHw32(i1,i2)) 
    g(SumI+478) = Real(mHg3p2(i1,i2), dp) 
    g(SumI+479) = Aimag(mHg3p2(i1,i2)) 
    g(SumI+496) = Real(mHb32(i1,i2), dp) 
    g(SumI+497) = Aimag(mHb32(i1,i2)) 
    g(SumI+514) = Real(mHx32(i1,i2), dp) 
    g(SumI+515) = Aimag(mHx32(i1,i2)) 
    g(SumI+532) = Real(mHxb32(i1,i2), dp) 
    g(SumI+533) = Aimag(mHxb32(i1,i2)) 
    g(SumI+556) = Real(MnuL(i1,i2), dp) 
    g(SumI+557) = Aimag(MnuL(i1,i2)) 
   End Do 
  End Do 

  g(112) = Real(mue,dp)  
  g(113) = Aimag(mue)
  g(294) = Real(Bmue,dp)  
  g(295) = Aimag(Bmue)  
  g(404) = mHd2  
  g(405) = mHu2
  g(550) = Real(MassB,dp)  
  g(551) = Aimag(MassB)  
  g(552) = Real(MassWB,dp)  
  g(553) = Aimag(MassWB)  
  g(554) = Real(MassG,dp)  
  g(555) = Aimag(MassG)  

  Iname = Iname - 1 
 
 End Subroutine ParametersToG555


Subroutine rge111(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i2
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,Dg3
Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          & 
& ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       & 
& ,DYe(3,3),adjYe(3,3),Yb3(3,3),betaYb31(3,3),betaYb32(3,3),DYb3(3,3),adjYb3(3,3)        & 
& ,Yw3(3,3),betaYw31(3,3),betaYw32(3,3),DYw3(3,3),adjYw3(3,3),Yx3(3,3),betaYx31(3,3)     & 
& ,betaYx32(3,3),DYx3(3,3),adjYx3(3,3)
Complex(dp) :: Yb3adjYb3(3,3),YdadjYd(3,3),YeadjYe(3,3),YuadjYu(3,3),Yw3adjYw3(3,3),Yx3adjYx3(3,3),  & 
& adjYb3Yb3(3,3),adjYdYd(3,3),adjYeYe(3,3),adjYuYu(3,3),adjYw3Yw3(3,3),adjYx3Yx3(3,3),   & 
& CYdTpYd(3,3),CYx3Yd(3,3),Yb3adjYb3Yb3(3,3),Yb3adjYeYe(3,3),Yb3adjYw3Yw3(3,3),          & 
& YdadjYdYd(3,3),YdadjYuYu(3,3),YeadjYb3Yb3(3,3),YeadjYeYe(3,3),YeadjYw3Yw3(3,3),        & 
& YuadjYdYd(3,3),YuadjYuYu(3,3),Yw3adjYb3Yb3(3,3),Yw3adjYeYe(3,3),Yw3adjYw3Yw3(3,3),     & 
& Yx3adjYx3Yx3(3,3),Yx3CYdTpYd(3,3),TpYx3CYx3Yd(3,3)

Complex(dp) :: YeadjYb3(3,3),YuadjYd(3,3),Yw3adjYb3(3,3),Yw3adjYe(3,3),CYuTpYd(3,3),TpYx3CYx3(3,3),  & 
& adjYb3Yb3adjYb3(3,3),adjYdYdadjYd(3,3),adjYdTpYx3CYx3(3,3),adjYeYeadjYb3(3,3),         & 
& adjYeYeadjYe(3,3),adjYuYuadjYd(3,3),adjYuYuadjYu(3,3),adjYw3Yw3adjYb3(3,3),            & 
& adjYw3Yw3adjYe(3,3),adjYw3Yw3adjYw3(3,3),adjYx3Yx3adjYx3(3,3),TpYdadjYx3Yx3(3,3),      & 
& TpYdCYdTpYd(3,3),TpYuCYuTpYd(3,3),Yb3adjYb3Yb3adjYb3(3,3),Yb3adjYeYeadjYb3(3,3),       & 
& Yb3adjYw3Yw3adjYb3(3,3),YdadjYdYdadjYd(3,3),YdadjYdTpYx3CYx3(3,3),YdadjYuYuadjYd(3,3), & 
& YeadjYeYeadjYe(3,3),YeadjYw3Yw3adjYe(3,3),YuadjYuYuadjYu(3,3),Yw3adjYw3Yw3adjYw3(3,3), & 
& Yx3adjYx3Yx3adjYx3(3,3),adjYb3Yb3adjYb3Yb3(3,3),adjYb3Yb3adjYeYe(3,3),adjYb3Yb3adjYw3Yw3(3,3),& 
& adjYdYdadjYdYd(3,3),adjYdYdadjYuYu(3,3),adjYdTpYx3CYx3Yd(3,3),adjYeYeadjYb3Yb3(3,3),   & 
& adjYeYeadjYeYe(3,3),adjYeYeadjYw3Yw3(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYuYu(3,3),     & 
& adjYw3Yw3adjYb3Yb3(3,3),adjYw3Yw3adjYeYe(3,3),adjYw3Yw3adjYw3Yw3(3,3),adjYx3Yx3adjYx3Yx3(3,3),& 
& CYdTpYdadjYx3Yx3(3,3),CYdTpYdCYdTpYd(3,3),CYdTpYuCYuTpYd(3,3),CYx3TpYx3CYx3Yd(3,3),    & 
& Yb3adjYb3Yb3adjYb3Yb3(3,3),Yb3adjYb3Yb3adjYw3Yw3(3,3),Yb3adjYeYeadjYb3Yb3(3,3),        & 
& Yb3adjYeYeadjYeYe(3,3),Yb3adjYw3Yw3adjYb3Yb3(3,3),Yb3adjYw3Yw3adjYw3Yw3(3,3),          & 
& YdadjYdYdadjYdYd(3,3),YdadjYdTpYx3CYx3Yd(3,3),YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYuYu(3,3),& 
& YeadjYb3Yb3adjYb3Yb3(3,3),YeadjYb3Yb3adjYeYe(3,3),YeadjYb3Yb3adjYw3Yw3(3,3),           & 
& YeadjYeYeadjYeYe(3,3),YeadjYw3Yw3adjYb3Yb3(3,3),YeadjYw3Yw3adjYeYe(3,3),               & 
& YeadjYw3Yw3adjYw3Yw3(3,3),YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYuYu(3,3),YuadjYdTpYx3CYx3Yd(3,3),& 
& YuadjYuYuadjYuYu(3,3),Yw3adjYb3Yb3adjYb3Yb3(3,3),Yw3adjYb3Yb3adjYw3Yw3(3,3),           & 
& Yw3adjYeYeadjYeYe(3,3),Yw3adjYeYeadjYw3Yw3(3,3),Yw3adjYw3Yw3adjYb3Yb3(3,3),            & 
& Yw3adjYw3Yw3adjYw3Yw3(3,3),Yx3adjYx3Yx3adjYx3Yx3(3,3),Yx3CYdTpYdadjYx3Yx3(3,3),        & 
& Yx3CYdTpYdCYdTpYd(3,3),Yx3CYdTpYuCYuTpYd(3,3),TpYx3CYx3TpYx3CYx3Yd(3,3)

Complex(dp) :: TrYb3adjYb3,TrYdadjYd,TrYeadjYe,TrYuadjYu,TrYw3adjYw3,TrYx3adjYx3

Complex(dp) :: TrYb3adjYb3Yb3adjYb3,TrYb3adjYeYeadjYb3,TrYb3adjYw3Yw3adjYb3,TrYdadjYdYdadjYd,        & 
& TrYdadjYdTpYx3CYx3,TrYdadjYuYuadjYd,TrYeadjYeYeadjYe,TrYeadjYw3Yw3adjYe,               & 
& TrYuadjYuYuadjYu,TrYw3adjYw3Yw3adjYw3,TrYx3adjYx3Yx3adjYx3

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Real(dp) :: g1p4,g2p4,g3p4

Iname = Iname +1 
NameOfUnit(Iname) = 'rge111' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters111(gy,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3)

If (ThresholdCrossed.lt.1) Then 
Yx3(1,:) = 0._dp 
Yb3(1,:) = 0._dp 
Yw3(1,:) = 0._dp 
End if 

If (ThresholdCrossed.lt.2) Then 
Yx3(2,:) = 0._dp 
Yb3(2,:) = 0._dp 
Yw3(2,:) = 0._dp 
End if 

If (ThresholdCrossed.lt.3) Then 
Yx3(3,:) = 0._dp 
Yb3(3,:) = 0._dp 
Yw3(3,:) = 0._dp 
End if 

Call Adjungate(Yu,adjYu)
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(Yb3,adjYb3)
Call Adjungate(Yw3,adjYw3)
Call Adjungate(Yx3,adjYx3)
 Yb3adjYb3 = Matmul2(Yb3,adjYb3,OnlyDiagonal) 
Forall(i2=1:3)  Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3(i2,i2),dp) 
 YdadjYd = Matmul2(Yd,adjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 YeadjYe = Matmul2(Ye,adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 YuadjYu = Matmul2(Yu,adjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 Yw3adjYw3 = Matmul2(Yw3,adjYw3,OnlyDiagonal) 
Forall(i2=1:3)  Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3(i2,i2),dp) 
 Yx3adjYx3 = Matmul2(Yx3,adjYx3,OnlyDiagonal) 
Forall(i2=1:3)  Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3(i2,i2),dp) 
 adjYb3Yb3 = Matmul2(adjYb3,Yb3,OnlyDiagonal) 
Forall(i2=1:3)  adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3(i2,i2),dp) 
 adjYdYd = Matmul2(adjYd,Yd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYeYe = Matmul2(adjYe,Ye,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYuYu = Matmul2(adjYu,Yu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYw3Yw3 = Matmul2(adjYw3,Yw3,OnlyDiagonal) 
Forall(i2=1:3)  adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3(i2,i2),dp) 
 adjYx3Yx3 = Matmul2(adjYx3,Yx3,OnlyDiagonal) 
Forall(i2=1:3)  adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3(i2,i2),dp) 
 CYdTpYd = Matmul2(Conjg(Yd),Transpose(Yd),OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYd(i2,i2) =  Real(CYdTpYd(i2,i2),dp) 
 CYx3Yd = Matmul2(Conjg(Yx3),Yd,OnlyDiagonal) 
 Yb3adjYb3Yb3 = Matmul2(Yb3,adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYeYe = Matmul2(Yb3,adjYeYe,OnlyDiagonal) 
 Yb3adjYw3Yw3 = Matmul2(Yb3,adjYw3Yw3,OnlyDiagonal) 
 YdadjYdYd = Matmul2(Yd,adjYdYd,OnlyDiagonal) 
 YdadjYuYu = Matmul2(Yd,adjYuYu,OnlyDiagonal) 
 YeadjYb3Yb3 = Matmul2(Ye,adjYb3Yb3,OnlyDiagonal) 
 YeadjYeYe = Matmul2(Ye,adjYeYe,OnlyDiagonal) 
 YeadjYw3Yw3 = Matmul2(Ye,adjYw3Yw3,OnlyDiagonal) 
 YuadjYdYd = Matmul2(Yu,adjYdYd,OnlyDiagonal) 
 YuadjYuYu = Matmul2(Yu,adjYuYu,OnlyDiagonal) 
 Yw3adjYb3Yb3 = Matmul2(Yw3,adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYeYe = Matmul2(Yw3,adjYeYe,OnlyDiagonal) 
 Yw3adjYw3Yw3 = Matmul2(Yw3,adjYw3Yw3,OnlyDiagonal) 
 Yx3adjYx3Yx3 = Matmul2(Yx3,adjYx3Yx3,OnlyDiagonal) 
 Yx3CYdTpYd = Matmul2(Yx3,CYdTpYd,OnlyDiagonal) 
 TpYx3CYx3Yd = Matmul2(Transpose(Yx3),CYx3Yd,OnlyDiagonal) 
 TrYb3adjYb3 = Real(cTrace(Yb3adjYb3),dp) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYw3adjYw3 = Real(cTrace(Yw3adjYw3),dp) 
 TrYx3adjYx3 = Real(cTrace(Yx3adjYx3),dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 


If (TwoLoopRGE) Then 
 YeadjYb3 = Matmul2(Ye,adjYb3,OnlyDiagonal) 
 YuadjYd = Matmul2(Yu,adjYd,OnlyDiagonal) 
 Yw3adjYb3 = Matmul2(Yw3,adjYb3,OnlyDiagonal) 
 Yw3adjYe = Matmul2(Yw3,adjYe,OnlyDiagonal) 
 CYuTpYd = Matmul2(Conjg(Yu),Transpose(Yd),OnlyDiagonal) 
 TpYx3CYx3 = Matmul2(Transpose(Yx3),Conjg(Yx3),OnlyDiagonal) 
Forall(i2=1:3)  TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3(i2,i2),dp) 
 adjYb3Yb3adjYb3 = Matmul2(adjYb3,Yb3adjYb3,OnlyDiagonal) 
 adjYdYdadjYd = Matmul2(adjYd,YdadjYd,OnlyDiagonal) 
 adjYdTpYx3CYx3 = Matmul2(adjYd,TpYx3CYx3,OnlyDiagonal) 
 adjYeYeadjYb3 = Matmul2(adjYe,YeadjYb3,OnlyDiagonal) 
 adjYeYeadjYe = Matmul2(adjYe,YeadjYe,OnlyDiagonal) 
 adjYuYuadjYd = Matmul2(adjYu,YuadjYd,OnlyDiagonal) 
 adjYuYuadjYu = Matmul2(adjYu,YuadjYu,OnlyDiagonal) 
 adjYw3Yw3adjYb3 = Matmul2(adjYw3,Yw3adjYb3,OnlyDiagonal) 
 adjYw3Yw3adjYe = Matmul2(adjYw3,Yw3adjYe,OnlyDiagonal) 
 adjYw3Yw3adjYw3 = Matmul2(adjYw3,Yw3adjYw3,OnlyDiagonal) 
 adjYx3Yx3adjYx3 = Matmul2(adjYx3,Yx3adjYx3,OnlyDiagonal) 
 TpYdadjYx3Yx3 = Matmul2(Transpose(Yd),adjYx3Yx3,OnlyDiagonal) 
 TpYdCYdTpYd = Matmul2(Transpose(Yd),CYdTpYd,OnlyDiagonal) 
 TpYuCYuTpYd = Matmul2(Transpose(Yu),CYuTpYd,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3 = Matmul2(Yb3,adjYb3Yb3adjYb3,OnlyDiagonal) 
Forall(i2=1:3)  Yb3adjYb3Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3Yb3adjYb3(i2,i2),dp) 
 Yb3adjYeYeadjYb3 = Matmul2(Yb3,adjYeYeadjYb3,OnlyDiagonal) 
Forall(i2=1:3)  Yb3adjYeYeadjYb3(i2,i2) =  Real(Yb3adjYeYeadjYb3(i2,i2),dp) 
 Yb3adjYw3Yw3adjYb3 = Matmul2(Yb3,adjYw3Yw3adjYb3,OnlyDiagonal) 
Forall(i2=1:3)  Yb3adjYw3Yw3adjYb3(i2,i2) =  Real(Yb3adjYw3Yw3adjYb3(i2,i2),dp) 
 YdadjYdYdadjYd = Matmul2(Yd,adjYdYdadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdTpYx3CYx3 = Matmul2(Yd,adjYdTpYx3CYx3,OnlyDiagonal) 
 YdadjYuYuadjYd = Matmul2(Yd,adjYuYuadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YeadjYeYeadjYe = Matmul2(Ye,adjYeYeadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYw3Yw3adjYe = Matmul2(Ye,adjYw3Yw3adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYw3Yw3adjYe(i2,i2) =  Real(YeadjYw3Yw3adjYe(i2,i2),dp) 
 YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 Yw3adjYw3Yw3adjYw3 = Matmul2(Yw3,adjYw3Yw3adjYw3,OnlyDiagonal) 
Forall(i2=1:3)  Yw3adjYw3Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3Yw3adjYw3(i2,i2),dp) 
 Yx3adjYx3Yx3adjYx3 = Matmul2(Yx3,adjYx3Yx3adjYx3,OnlyDiagonal) 
Forall(i2=1:3)  Yx3adjYx3Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3Yx3adjYx3(i2,i2),dp) 
 adjYb3Yb3adjYb3Yb3 = Matmul2(adjYb3,Yb3adjYb3Yb3,OnlyDiagonal) 
Forall(i2=1:3)  adjYb3Yb3adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3adjYb3Yb3(i2,i2),dp) 
 adjYb3Yb3adjYeYe = Matmul2(adjYb3,Yb3adjYeYe,OnlyDiagonal) 
 adjYb3Yb3adjYw3Yw3 = Matmul2(adjYb3,Yb3adjYw3Yw3,OnlyDiagonal) 
 adjYdYdadjYdYd = Matmul2(adjYd,YdadjYdYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYuYu = Matmul2(adjYd,YdadjYuYu,OnlyDiagonal) 
 adjYdTpYx3CYx3Yd = Matmul2(adjYd,TpYx3CYx3Yd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdTpYx3CYx3Yd(i2,i2) =  Real(adjYdTpYx3CYx3Yd(i2,i2),dp) 
 adjYeYeadjYb3Yb3 = Matmul2(adjYe,YeadjYb3Yb3,OnlyDiagonal) 
 adjYeYeadjYeYe = Matmul2(adjYe,YeadjYeYe,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYw3Yw3 = Matmul2(adjYe,YeadjYw3Yw3,OnlyDiagonal) 
 adjYuYuadjYdYd = Matmul2(adjYu,YuadjYdYd,OnlyDiagonal) 
 adjYuYuadjYuYu = Matmul2(adjYu,YuadjYuYu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYw3Yw3adjYb3Yb3 = Matmul2(adjYw3,Yw3adjYb3Yb3,OnlyDiagonal) 
 adjYw3Yw3adjYeYe = Matmul2(adjYw3,Yw3adjYeYe,OnlyDiagonal) 
 adjYw3Yw3adjYw3Yw3 = Matmul2(adjYw3,Yw3adjYw3Yw3,OnlyDiagonal) 
Forall(i2=1:3)  adjYw3Yw3adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3adjYw3Yw3(i2,i2),dp) 
 adjYx3Yx3adjYx3Yx3 = Matmul2(adjYx3,Yx3adjYx3Yx3,OnlyDiagonal) 
Forall(i2=1:3)  adjYx3Yx3adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3adjYx3Yx3(i2,i2),dp) 
 CYdTpYdadjYx3Yx3 = Matmul2(Conjg(Yd),TpYdadjYx3Yx3,OnlyDiagonal) 
 CYdTpYdCYdTpYd = Matmul2(Conjg(Yd),TpYdCYdTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYdCYdTpYd(i2,i2) =  Real(CYdTpYdCYdTpYd(i2,i2),dp) 
 CYdTpYuCYuTpYd = Matmul2(Conjg(Yd),TpYuCYuTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYuCYuTpYd(i2,i2) =  Real(CYdTpYuCYuTpYd(i2,i2),dp) 
 CYx3TpYx3CYx3Yd = Matmul2(Conjg(Yx3),TpYx3CYx3Yd,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3Yb3 = Matmul2(Yb3,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYw3Yw3 = Matmul2(Yb3,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 Yb3adjYeYeadjYb3Yb3 = Matmul2(Yb3,adjYeYeadjYb3Yb3,OnlyDiagonal) 
 Yb3adjYeYeadjYeYe = Matmul2(Yb3,adjYeYeadjYeYe,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYb3Yb3 = Matmul2(Yb3,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYw3Yw3 = Matmul2(Yb3,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 YdadjYdYdadjYdYd = Matmul2(Yd,adjYdYdadjYdYd,OnlyDiagonal) 
 YdadjYdTpYx3CYx3Yd = Matmul2(Yd,adjYdTpYx3CYx3Yd,OnlyDiagonal) 
 YdadjYuYuadjYdYd = Matmul2(Yd,adjYuYuadjYdYd,OnlyDiagonal) 
 YdadjYuYuadjYuYu = Matmul2(Yd,adjYuYuadjYuYu,OnlyDiagonal) 
 YeadjYb3Yb3adjYb3Yb3 = Matmul2(Ye,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 YeadjYb3Yb3adjYeYe = Matmul2(Ye,adjYb3Yb3adjYeYe,OnlyDiagonal) 
 YeadjYb3Yb3adjYw3Yw3 = Matmul2(Ye,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 YeadjYeYeadjYeYe = Matmul2(Ye,adjYeYeadjYeYe,OnlyDiagonal) 
 YeadjYw3Yw3adjYb3Yb3 = Matmul2(Ye,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 YeadjYw3Yw3adjYeYe = Matmul2(Ye,adjYw3Yw3adjYeYe,OnlyDiagonal) 
 YeadjYw3Yw3adjYw3Yw3 = Matmul2(Ye,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 YuadjYdYdadjYdYd = Matmul2(Yu,adjYdYdadjYdYd,OnlyDiagonal) 
 YuadjYdYdadjYuYu = Matmul2(Yu,adjYdYdadjYuYu,OnlyDiagonal) 
 YuadjYdTpYx3CYx3Yd = Matmul2(Yu,adjYdTpYx3CYx3Yd,OnlyDiagonal) 
 YuadjYuYuadjYuYu = Matmul2(Yu,adjYuYuadjYuYu,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYb3Yb3 = Matmul2(Yw3,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYw3Yw3 = Matmul2(Yw3,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 Yw3adjYeYeadjYeYe = Matmul2(Yw3,adjYeYeadjYeYe,OnlyDiagonal) 
 Yw3adjYeYeadjYw3Yw3 = Matmul2(Yw3,adjYeYeadjYw3Yw3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYb3Yb3 = Matmul2(Yw3,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYw3Yw3 = Matmul2(Yw3,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 Yx3adjYx3Yx3adjYx3Yx3 = Matmul2(Yx3,adjYx3Yx3adjYx3Yx3,OnlyDiagonal) 
 Yx3CYdTpYdadjYx3Yx3 = Matmul2(Yx3,CYdTpYdadjYx3Yx3,OnlyDiagonal) 
 Yx3CYdTpYdCYdTpYd = Matmul2(Yx3,CYdTpYdCYdTpYd,OnlyDiagonal) 
 Yx3CYdTpYuCYuTpYd = Matmul2(Yx3,CYdTpYuCYuTpYd,OnlyDiagonal) 
 TpYx3CYx3TpYx3CYx3Yd = Matmul2(Transpose(Yx3),CYx3TpYx3CYx3Yd,OnlyDiagonal) 
 TrYb3adjYb3Yb3adjYb3 = cTrace(Yb3adjYb3Yb3adjYb3) 
 TrYb3adjYeYeadjYb3 = cTrace(Yb3adjYeYeadjYb3) 
 TrYb3adjYw3Yw3adjYb3 = cTrace(Yb3adjYw3Yw3adjYb3) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdTpYx3CYx3 = cTrace(YdadjYdTpYx3CYx3) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYw3Yw3adjYe = cTrace(YeadjYw3Yw3adjYe) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYw3adjYw3Yw3adjYw3 = cTrace(Yw3adjYw3Yw3adjYw3) 
 TrYx3adjYx3Yx3adjYx3 = cTrace(Yx3adjYx3Yx3adjYx3) 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = (g1p3*(66 + 25*NGHx3 + 25*NGHxb3))/10._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (g1p3*(6*(199._dp*(g1p2) + 135._dp*(g2p2) + 440._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -    & 
&  70._dp*(TrYdadjYd) - 90._dp*(TrYeadjYe) - 130._dp*(TrYuadjYu) - 60._dp*(TrYw3adjYw3) -& 
&  190._dp*(TrYx3adjYx3)) + 125*(5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 +& 
&  NGHxb3)))/150._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = (g2p3*(2 + 4*NGHw3 + 3*NGHx3 +               & 
&  3*NGHxb3))/2._dp

 
 
If (TwoLoopRGE) Then 
betag22 = (g2p3*(2*(27._dp*(g1p2) + 375._dp*(g2p2) + 360._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -     & 
&  90._dp*(TrYdadjYd) - 30._dp*(TrYeadjYe) - 90._dp*(TrYuadjYu) - 140._dp*(TrYw3adjYw3) -& 
&  90._dp*(TrYx3adjYx3)) + 720*g2p2*NGHw3 + 15*(5._dp*(g1p2) +            & 
&  21._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 + NGHxb3)))/30._dp

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = g3p3*(-3 + 3*NGHg3 + NGHx3 + NGHxb3)

 
 
If (TwoLoopRGE) Then 
betag32 = (g3p3*(3*(11._dp*(g1p2) + 45._dp*(g2p2) + 70._dp*(g3p2) - 20._dp*(TrYdadjYd) -        & 
&  20._dp*(TrYuadjYu) - 20._dp*(TrYx3adjYx3)) + 810*g3p2*NGHg3 +          & 
&  5*(5._dp*(g1p2) + 9._dp*(g2p2) + 34._dp*(g3p2))*(NGHx3 +               & 
&  NGHxb3)))/15._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = ((-26._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)    & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yu)               & 
& /30._dp + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = (2._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YuadjYdYd + ((4._dp*(g1p2) +     & 
&  60._dp*(g2p2) - 9._dp*(TrYb3adjYb3) - 90._dp*(TrYuadjYu) - 45._dp*(TrYw3adjYw3) -     & 
&  90._dp*(TrYx3adjYx3))*YuadjYuYu)/10._dp - 2*(YuadjYdTpYx3CYx3Yd + YuadjYdYdadjYdYd +  & 
&  YuadjYdYdadjYuYu + 2._dp*(YuadjYuYuadjYuYu)) + (Yu*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - & 
&  540._dp*(TrYb3adjYeYeadjYb3) - 2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(5486._dp*(g1p4) + & 
&  20*g1p2*(45._dp*(g2p2) + 136._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 -       & 
&  32._dp*(g3p4)) - 5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) -        & 
&  1350._dp*(TrYeadjYw3Yw3adjYe) - 8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  8100._dp*(TrYx3adjYx3Yx3adjYx3)) + 60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) +& 
&  65*g1p2*(NGHx3 + NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu +   & 
&  TrYx3adjYx3) + g3p2*(3*NGHg3 + NGHx3 + NGHxb3)) +& 
&  2700*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/1800._dp

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = 2._dp*(TpYx3CYx3Yd) + (-7._dp*(g1p2)/15._dp - 3._dp*(g2p2) -               & 
&  16._dp*(g3p2)/3._dp + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd)           & 
&  + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = (4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) - 3._dp*(TrYeadjYe))*YdadjYdYd +& 
&  (2*TpYx3CYx3Yd*(-3._dp*(TrYb3adjYb3) + 5*(2._dp*(g1p2) + 6._dp*(g2p2) -               & 
&  6._dp*(TrYuadjYu) - 3._dp*(TrYw3adjYw3) - 6._dp*(TrYx3adjYx3))) + (8._dp*(g1p2) -     & 
&  3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*YdadjYuYu)/10._dp -& 
&  2*(TpYx3CYx3TpYx3CYx3Yd + YdadjYdTpYx3CYx3Yd + 2._dp*(YdadjYdYdadjYdYd) +             & 
&  YdadjYuYuadjYdYd + YdadjYuYuadjYuYu) + (Yd*(287._dp*(g1p4) + 90*g1p2*g2p2 +           & 
&  675._dp*(g2p4) + 80*(g1p2 + 9._dp*(g2p2))*g3p2 - 160._dp*(g3p4) - 27._dp*(TrYb3adjYeYeadjYb3) -& 
&  540._dp*(TrYdadjYdTpYx3CYx3) - 810._dp*(TrYdadjYdYdadjYd) - 270._dp*(TrYdadjYuYuadjYd) -& 
&  270._dp*(TrYeadjYeYeadjYe) - 135._dp*(TrYeadjYw3Yw3adjYe) + 135*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3)) + 3*g1p2*(-12*(TrYdadjYd -          & 
&  3._dp*(TrYeadjYe)) + 35*g1p2*(NGHx3 + NGHxb3)) +        & 
&  480*g3p2*(3._dp*(TrYdadjYd) + g3p2*(3*NGHg3 + NGHx3 +   & 
&  NGHxb3))))/90._dp

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (-9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)       & 
& *Ye + (3*(YeadjYb3Yb3 + 5*(2._dp*(YeadjYeYe) + YeadjYw3Yw3)))/10._dp

 
 
If (TwoLoopRGE) Then 
betaYe2 = (-36._dp*(YeadjYb3Yb3adjYb3Yb3) - 18*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) +            & 
&  TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(YeadjYb3Yb3 + 5._dp*(YeadjYw3Yw3)) -             & 
&  600*((-2._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)*YeadjYeYe - 2*g2p2*YeadjYw3Yw3) -& 
&  5*(24._dp*(YeadjYb3Yb3adjYeYe) + 12._dp*(YeadjYb3Yb3adjYw3Yw3) + 160._dp*(YeadjYeYeadjYeYe) +& 
&  9._dp*(YeadjYw3Yw3adjYb3Yb3) + 60*(2._dp*(YeadjYw3Yw3adjYeYe) + YeadjYw3Yw3adjYw3Yw3)) +& 
&  20*Ye*(-4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 12*g1p2*TrYeadjYe + 3*(45._dp*(g1p4) +   & 
&  6*g1p2*g2p2 + 25._dp*(g2p4) - TrYb3adjYeYeadjYb3 - 20._dp*(TrYdadjYdTpYx3CYx3) -      & 
&  30._dp*(TrYdadjYdYdadjYd) - 10._dp*(TrYdadjYuYuadjYd) - 10._dp*(TrYeadjYeYeadjYe) -   & 
&  5._dp*(TrYeadjYw3Yw3adjYe)) + 45*g1p4*NGHx3 + 45*g1p4*NGHxb3 +& 
&  15*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))/200._dp

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yb3 
!-------------------- 
 
betaYb31  = (3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) -              & 
&  2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*Yb3)/10._dp + 9._dp*(Yb3adjYb3Yb3)     & 
& /10._dp + Yb3adjYeYe + 3._dp*(Yb3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYb32 = (18*(4._dp*(g1p2) + 20._dp*(g2p2) - 3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) -        & 
&  15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*Yb3adjYb3Yb3 - 72._dp*(Yb3adjYb3Yb3adjYb3Yb3) +& 
&  40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yb3adjYeYe -               & 
&  30*(3._dp*(TrYb3adjYb3) + 5*(-8._dp*(g2p2) + 6._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3) +& 
&  6._dp*(TrYx3adjYx3)))*Yb3adjYw3Yw3 - 5*(12._dp*(Yb3adjYb3Yb3adjYw3Yw3) +              & 
&  24._dp*(Yb3adjYeYeadjYb3Yb3) + 80._dp*(Yb3adjYeYeadjYeYe) + 45._dp*(Yb3adjYw3Yw3adjYb3Yb3) +& 
&  60._dp*(Yb3adjYw3Yw3adjYw3Yw3)) + Yb3*(3*(276._dp*(g1p4) + 120*g1p2*g2p2 +            & 
&  500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) - 20._dp*(TrYb3adjYeYeadjYb3) -        & 
&  95._dp*(TrYb3adjYw3Yw3adjYb3) - 400._dp*(TrYdadjYdTpYx3CYx3) - 200._dp*(TrYdadjYuYuadjYd) -& 
&  100._dp*(TrYeadjYw3Yw3adjYe) - 600._dp*(TrYuadjYuYuadjYu) - 250._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  600._dp*(TrYx3adjYx3Yx3adjYx3)) + 20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu +            & 
&  20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 + 15*g1p4*NGHx3 + 15*g1p4*NGHxb3) +& 
&  300*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/200._dp

 
DYb3 = oo16pi2*( betaYb31 + oo16pi2 * betaYb32 ) 

 
Else 
DYb3 = oo16pi2* betaYb31 
End If 
 
 
!-------------------- 
! Yw3 
!-------------------- 
 
betaYw31  = (-3._dp*(g1p2)/5._dp - 7._dp*(g2p2) + 3._dp*(TrYb3adjYb3)/10._dp +        & 
&  3._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3)/2._dp + 3._dp*(TrYx3adjYx3))*Yw3 +            & 
&  3._dp*(Yw3adjYb3Yb3)/10._dp + Yw3adjYeYe + 5._dp*(Yw3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYw32 = (-36._dp*(Yw3adjYb3Yb3adjYb3Yb3) + 40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) -            & 
&  5._dp*(TrYeadjYe))*Yw3adjYeYe + 40*(3._dp*(g1p2) + 25._dp*(g2p2))*Yw3adjYw3Yw3 -      & 
&  6*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(3._dp*(Yw3adjYb3Yb3) +& 
&  25._dp*(Yw3adjYw3Yw3)) - 5*(24._dp*(Yw3adjYb3Yb3adjYw3Yw3) + 80._dp*(Yw3adjYeYeadjYeYe) +& 
&  40._dp*(Yw3adjYeYeadjYw3Yw3) + 9._dp*(Yw3adjYw3Yw3adjYb3Yb3) + 120._dp*(Yw3adjYw3Yw3adjYw3Yw3)) +& 
&  Yw3*(828._dp*(g1p4) + 360*g1p2*g2p2 + 5500._dp*(g2p4) - 54._dp*(TrYb3adjYb3Yb3adjYb3) -& 
&  60._dp*(TrYb3adjYeYeadjYb3) - 285._dp*(TrYb3adjYw3Yw3adjYb3) - 1200._dp*(TrYdadjYdTpYx3CYx3) -& 
&  600._dp*(TrYdadjYuYuadjYd) - 300._dp*(TrYeadjYw3Yw3adjYe) - 1800._dp*(TrYuadjYuYuadjYu) +& 
&  1200*g2p2*TrYw3adjYw3 - 750._dp*(TrYw3adjYw3Yw3adjYw3) - 1800._dp*(TrYx3adjYx3Yx3adjYx3) +& 
&  20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 +       & 
&  15*g1p4*NGHx3 + 15*g1p4*NGHxb3) + 700*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3))))/200._dp

 
DYw3 = oo16pi2*( betaYw31 + oo16pi2 * betaYw32 ) 

 
Else 
DYw3 = oo16pi2* betaYw31 
End If 
 
 
!-------------------- 
! Yx3 
!-------------------- 
 
betaYx31  = ((-38._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)   & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yx3)              & 
& /30._dp + 3._dp*(Yx3adjYx3Yx3) + 2._dp*(Yx3CYdTpYd)

 
 
If (TwoLoopRGE) Then 
betaYx32 = (8._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYb3adjYb3)/10._dp - 9._dp*(TrYuadjYu) - & 
&  9._dp*(TrYw3adjYw3)/2._dp - 9._dp*(TrYx3adjYx3))*Yx3adjYx3Yx3 + (2*(g1p2 +            & 
&  15._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yx3CYdTpYd)/5._dp -           & 
&  2*(2._dp*(Yx3adjYx3Yx3adjYx3Yx3) + Yx3CYdTpYdadjYx3Yx3 + Yx3CYdTpYdCYdTpYd +          & 
&  Yx3CYdTpYuCYuTpYd) + (Yx3*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(8246._dp*(g1p4) + 20*g1p2*(153._dp*(g2p2) +      & 
&  232._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 - 32._dp*(g3p4)) -               & 
&  5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) - 1350._dp*(TrYeadjYw3Yw3adjYe) -& 
&  8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) - 8100._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 95*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/1800._dp

 
DYx3 = oo16pi2*( betaYx31 + oo16pi2 * betaYx32 ) 

 
Else 
DYx3 = oo16pi2* betaYx31 
End If 
 
 
If (ThresholdCrossed.lt.1) Then 
DYx3(1,:) = 0._dp 
DYb3(1,:) = 0._dp 
DYw3(1,:) = 0._dp 
End if 

If (ThresholdCrossed.lt.2) Then 
DYx3(2,:) = 0._dp 
DYb3(2,:) = 0._dp 
DYw3(2,:) = 0._dp 
End if 

If (ThresholdCrossed.lt.3) Then 
DYx3(3,:) = 0._dp 
DYb3(3,:) = 0._dp 
DYw3(3,:) = 0._dp 
End if 

Call ParametersToG111(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYb3,DYw3,DYx3,f)

Iname = Iname - 1 
 
End Subroutine rge111


Subroutine rge555(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,         & 
& Dg3,mHd2,betamHd21,betamHd22,DmHd2,mHu2,betamHu21,betamHu22,DmHu2
Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)          & 
& ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       & 
& ,DYe(3,3),adjYe(3,3),Yb3(3,3),betaYb31(3,3),betaYb32(3,3),DYb3(3,3),adjYb3(3,3)        & 
& ,Yw3(3,3),betaYw31(3,3),betaYw32(3,3),DYw3(3,3),adjYw3(3,3),Yx3(3,3),betaYx31(3,3)     & 
& ,betaYx32(3,3),DYx3(3,3),adjYx3(3,3),mue,betamue1,betamue2,Dmue,MXM3(3,3)              & 
& ,betaMXM31(3,3),betaMXM32(3,3),DMXM3(3,3),adjMXM3(3,3),MWM3(3,3),betaMWM31(3,3)        & 
& ,betaMWM32(3,3),DMWM3(3,3),adjMWM3(3,3),MGM3(3,3),betaMGM31(3,3),betaMGM32(3,3)        & 
& ,DMGM3(3,3),adjMGM3(3,3),MBM3(3,3),betaMBM31(3,3),betaMBM32(3,3),DMBM3(3,3)            & 
& ,adjMBM3(3,3),TYu(3,3),betaTYu1(3,3),betaTYu2(3,3),DTYu(3,3),adjTYu(3,3)               & 
& ,TYd(3,3),betaTYd1(3,3),betaTYd2(3,3),DTYd(3,3),adjTYd(3,3),TYe(3,3),betaTYe1(3,3)     & 
& ,betaTYe2(3,3),DTYe(3,3),adjTYe(3,3),TYb3(3,3),betaTYb31(3,3),betaTYb32(3,3)           & 
& ,DTYb3(3,3),adjTYb3(3,3),TYw3(3,3),betaTYw31(3,3),betaTYw32(3,3),DTYw3(3,3)            & 
& ,adjTYw3(3,3),TYx3(3,3),betaTYx31(3,3),betaTYx32(3,3),DTYx3(3,3),adjTYx3(3,3)          & 
& ,Bmue,betaBmue1,betaBmue2,DBmue,BMXM3(3,3),betaBMXM31(3,3),betaBMXM32(3,3)             & 
& ,DBMXM3(3,3),adjBMXM3(3,3),BMWM3(3,3),betaBMWM31(3,3),betaBMWM32(3,3),DBMWM3(3,3)      & 
& ,adjBMWM3(3,3),BMGM3(3,3),betaBMGM31(3,3),betaBMGM32(3,3),DBMGM3(3,3),adjBMGM3(3,3)    & 
& ,BMBM3(3,3),betaBMBM31(3,3),betaBMBM32(3,3),DBMBM3(3,3),adjBMBM3(3,3),mq2(3,3)         & 
& ,betamq21(3,3),betamq22(3,3),Dmq2(3,3),adjmq2(3,3),ml2(3,3),betaml21(3,3)              & 
& ,betaml22(3,3),Dml2(3,3),adjml2(3,3),md2(3,3),betamd21(3,3),betamd22(3,3)              & 
& ,Dmd2(3,3),adjmd2(3,3),mu2(3,3),betamu21(3,3),betamu22(3,3),Dmu2(3,3),adjmu2(3,3)      & 
& ,me2(3,3),betame21(3,3),betame22(3,3),Dme2(3,3),adjme2(3,3),mHw32(3,3),betamHw321(3,3) & 
& ,betamHw322(3,3),DmHw32(3,3),adjmHw32(3,3),mHg32(3,3),betamHg321(3,3),betamHg322(3,3)  & 
& ,DmHg32(3,3),adjmHg32(3,3),mHb32(3,3),betamHb321(3,3),betamHb322(3,3),DmHb32(3,3)      & 
& ,adjmHb32(3,3),mHx32(3,3),betamHx321(3,3),betamHx322(3,3),DmHx32(3,3),adjmHx32(3,3)    & 
& ,mHxb32(3,3),betamHxb321(3,3),betamHxb322(3,3),DmHxb32(3,3),adjmHxb32(3,3)             & 
& ,MassB,betaMassB1,betaMassB2,DMassB,MassWB,betaMassWB1,betaMassWB2,DMassWB,MassG,betaMassG1,betaMassG2,DMassG
Complex(dp) :: Tr1(3),Tr2(3),Tr3(3) 
Real(dp) :: AbsMassB,AbsMassWB,AbsMassG
Complex(dp) :: md2adjYx3(3,3),md2CYd(3,3),me2CYe(3,3),mHb32CYb3(3,3),mHw32CYw3(3,3),mHxb32CYx3(3,3), & 
& ml2adjYb3(3,3),ml2adjYe(3,3),ml2adjYw3(3,3),mq2adjYd(3,3),mq2adjYu(3,3),               & 
& mu2CYu(3,3),Yb3adjYb3(3,3),YdadjYd(3,3),YeadjYe(3,3),YuadjYu(3,3),Yw3adjYw3(3,3),      & 
& Yx3adjYx3(3,3),adjYb3MBM3(3,3),adjYb3mHb32(3,3),adjYb3Yb3(3,3),adjYb3BMBM3(3,3),       & 
& adjYb3TYb3(3,3),adjYdmd2(3,3),adjYdYd(3,3),adjYdTYd(3,3),adjYeme2(3,3),adjYeYe(3,3),   & 
& adjYeTYe(3,3),adjYumu2(3,3),adjYuYu(3,3),adjYuTYu(3,3),adjYw3mHw32(3,3),               & 
& adjYw3MWM3(3,3),adjYw3Yw3(3,3),adjYw3BMWM3(3,3),adjYw3TYw3(3,3),adjYx3mHxb32(3,3),     & 
& adjYx3Yx3(3,3),adjYx3TYx3(3,3),CYb3ml2(3,3),CYb3TpYb3(3,3),CYb3TpTYb3(3,3),            & 
& CYdmq2(3,3),CYdTpYd(3,3),CYdTpTYd(3,3),CYeml2(3,3),CYumq2(3,3),CYw3ml2(3,3),           & 
& CYw3TpYw3(3,3),CYw3TpTYw3(3,3),CYx3md2(3,3),CYx3Yd(3,3),CYx3TYd(3,3),CYx3TpYx3(3,3),   & 
& CYx3TpTYx3(3,3),CTYb3TpTYb3(3,3),CTYdTpTYd(3,3),CTYeTpTYe(3,3),CTYuTpTYu(3,3),         & 
& CTYw3TpTYw3(3,3),CTYx3TpTYx3(3,3),TYb3adjTYb3(3,3),TYdadjTYd(3,3),TYeadjTYe(3,3),      & 
& TYuadjTYu(3,3),TYw3adjTYw3(3,3),TYx3adjTYx3(3,3),TpYb3CYb3(3,3),TpYdCYd(3,3),          & 
& TpYeCYe(3,3),TpYuCYu(3,3),TpYw3CYw3(3,3),TpYx3CYx3(3,3),TpTYb3CTYb3(3,3),              & 
& TpTYdCTYd(3,3),TpTYeCTYe(3,3),TpTYuCTYu(3,3),TpTYw3CTYw3(3,3),TpTYx3CTYx3(3,3),        & 
& MBM3CYb3TpYb3(3,3),MBM3CYb3TpTYb3(3,3),md2YdadjYd(3,3),md2adjYx3Yx3(3,3),              & 
& md2TpYx3CYx3(3,3),me2YeadjYe(3,3),mHb32Yb3adjYb3(3,3),mHw32Yw3adjYw3(3,3),             & 
& mHxb32Yx3adjYx3(3,3),ml2adjYb3Yb3(3,3),ml2adjYeYe(3,3),ml2adjYw3Yw3(3,3),              & 
& ml2TpYb3CYb3(3,3),ml2TpYeCYe(3,3),ml2TpYw3CYw3(3,3),mq2adjYdYd(3,3),mq2adjYuYu(3,3),   & 
& mq2TpYdCYd(3,3),mq2TpYuCYu(3,3),mu2YuadjYu(3,3),MWM3CYw3TpYw3(3,3),MWM3CYw3TpTYw3(3,3),& 
& MXM3CYx3TpYx3(3,3),MXM3CYx3TpTYx3(3,3),Yb3ml2adjYb3(3,3),Yb3adjYb3MBM3(3,3),           & 
& Yb3adjYb3mHb32(3,3),Yb3adjYb3Yb3(3,3),Yb3adjYb3BMBM3(3,3),Yb3adjYb3TYb3(3,3),          & 
& Yb3adjYeYe(3,3),Yb3adjYeTYe(3,3),Yb3adjYw3Yw3(3,3),Yb3adjYw3TYw3(3,3),Ydmq2adjYd(3,3), & 
& YdadjYdmd2(3,3),YdadjYdYd(3,3),YdadjYdTYd(3,3),YdadjYuYu(3,3),YdadjYuTYu(3,3),         & 
& Yeml2adjYe(3,3),YeadjYb3Yb3(3,3),YeadjYb3TYb3(3,3),YeadjYeme2(3,3),YeadjYeYe(3,3),     & 
& YeadjYeTYe(3,3),YeadjYw3Yw3(3,3),YeadjYw3TYw3(3,3),Yumq2adjYu(3,3),YuadjYdYd(3,3),     & 
& YuadjYdTYd(3,3),YuadjYumu2(3,3),YuadjYuYu(3,3),YuadjYuTYu(3,3),Yw3ml2adjYw3(3,3),      & 
& Yw3adjYb3Yb3(3,3),Yw3adjYb3TYb3(3,3),Yw3adjYeYe(3,3),Yw3adjYeTYe(3,3),Yw3adjYw3mHw32(3,3),& 
& Yw3adjYw3MWM3(3,3),Yw3adjYw3Yw3(3,3),Yw3adjYw3BMWM3(3,3),Yw3adjYw3TYw3(3,3),           & 
& Yx3md2adjYx3(3,3),Yx3adjYx3mHxb32(3,3),Yx3adjYx3Yx3(3,3),Yx3adjYx3TYx3(3,3),           & 
& Yx3CYdTpYd(3,3),Yx3CYdTpTYd(3,3),BMBM3CYb3TpYb3(3,3),BMWM3CYw3TpYw3(3,3),              & 
& BMXM3CYx3TpYx3(3,3),TYb3adjYb3MBM3(3,3),TYb3adjYb3Yb3(3,3),TYb3adjYeYe(3,3),           & 
& TYb3adjYw3Yw3(3,3),TYdadjYdYd(3,3),TYdadjYuYu(3,3),TYeadjYb3Yb3(3,3),TYeadjYeYe(3,3),  & 
& TYeadjYw3Yw3(3,3),TYuadjYdYd(3,3),TYuadjYuYu(3,3),TYw3adjYb3Yb3(3,3),TYw3adjYeYe(3,3), & 
& TYw3adjYw3MWM3(3,3),TYw3adjYw3Yw3(3,3),TYx3adjYx3Yx3(3,3),TYx3CYdTpYd(3,3),            & 
& TpYb3mHb32CYb3(3,3),TpYb3CYb3ml2(3,3),TpYdmd2CYd(3,3),TpYdCYdmq2(3,3),TpYeme2CYe(3,3), & 
& TpYeCYeml2(3,3),TpYumu2CYu(3,3),TpYuCYumq2(3,3),TpYw3mHw32CYw3(3,3),TpYw3CYw3ml2(3,3)

Complex(dp) :: TpYx3mHxb32CYx3(3,3),TpYx3CYx3md2(3,3),TpYx3CYx3Yd(3,3),TpYx3CYx3TYd(3,3),             & 
& TpTYx3CYx3Yd(3,3)

Complex(dp) :: Yb3adjYe(3,3),Yb3adjYw3(3,3),Yb3adjTYb3(3,3),Yb3adjTYe(3,3),Yb3adjTYw3(3,3),          & 
& YdadjYu(3,3),YdadjTYd(3,3),YdadjTYu(3,3),YeadjYb3(3,3),YeadjYw3(3,3),YeadjTYb3(3,3),   & 
& YeadjTYe(3,3),YeadjTYw3(3,3),YuadjYd(3,3),YuadjTYd(3,3),YuadjTYu(3,3),Yw3adjYb3(3,3),  & 
& Yw3adjYe(3,3),Yw3adjTYb3(3,3),Yw3adjTYe(3,3),Yw3adjTYw3(3,3),Yx3adjTYx3(3,3),          & 
& Yx3CYd(3,3),Yx3CTYd(3,3),adjYdadjTYx3(3,3),adjYdTpYx3(3,3),adjYdTpTYx3(3,3),           & 
& adjTYdadjYx3(3,3),CYb3TpYw3(3,3),CYb3TpTYw3(3,3),CYeTpYb3(3,3),CYeTpYw3(3,3),          & 
& CYeTpTYb3(3,3),CYeTpTYw3(3,3),CYuTpYd(3,3),CYuTpTYd(3,3),CYw3TpYb3(3,3),               & 
& CYw3TpTYb3(3,3),CTYb3adjYb3(3,3),CTYb3adjYe(3,3),CTYb3adjYw3(3,3),CTYb3TpYb3(3,3),     & 
& CTYdadjYd(3,3),CTYdadjYu(3,3),CTYdTpYd(3,3),CTYeadjYb3(3,3),CTYeadjYe(3,3),            & 
& CTYeadjYw3(3,3),CTYeTpYe(3,3),CTYuadjYd(3,3),CTYuadjYu(3,3),CTYuTpYu(3,3),             & 
& CTYw3adjYb3(3,3),CTYw3adjYe(3,3),CTYw3adjYw3(3,3),CTYw3TpYw3(3,3),CTYx3adjYx3(3,3),    & 
& CTYx3TYd(3,3),CTYx3TpYd(3,3),CTYx3TpYx3(3,3),CTYx3TpTYd(3,3),TYb3adjYb3(3,3),          & 
& TYb3adjYe(3,3),TYb3adjYw3(3,3),TYb3adjTYe(3,3),TYb3adjTYw3(3,3),TYdadjYd(3,3),         & 
& TYdadjYu(3,3),TYdadjTYu(3,3),TYeadjYb3(3,3),TYeadjYe(3,3),TYeadjYw3(3,3),              & 
& TYeadjTYb3(3,3),TYeadjTYw3(3,3),TYuadjYd(3,3),TYuadjYu(3,3),TYuadjTYd(3,3),            & 
& TYw3adjYb3(3,3),TYw3adjYe(3,3),TYw3adjYw3(3,3),TYw3adjTYb3(3,3),TYw3adjTYe(3,3),       & 
& TYx3adjYx3(3,3),TYx3CTYd(3,3),TpYb3CTYb3(3,3),TpYdadjYx3(3,3),TpYdadjTYx3(3,3),        & 
& TpYdCTYd(3,3),TpYeCTYe(3,3),TpYuCTYu(3,3),TpYw3CTYw3(3,3),TpYx3CTYx3(3,3),             & 
& TpTYb3CYb3(3,3),TpTYdadjYx3(3,3),TpTYdCYd(3,3),TpTYeCYe(3,3),TpTYuCYu(3,3),            & 
& TpTYw3CYw3(3,3),TpTYx3CYx3(3,3),md2YdadjYu(3,3),me2YeadjYb3(3,3),me2YeadjYw3(3,3),     & 
& mHb32Yb3adjYe(3,3),mHb32Yb3adjYw3(3,3),mHw32Yw3adjYb3(3,3),mHw32Yw3adjYe(3,3),         & 
& mHxb32Yx3CYd(3,3),mq2TpYdadjYx3(3,3),mu2YuadjYd(3,3),Yb3ml2adjYe(3,3),Yb3ml2adjYw3(3,3),& 
& Yb3adjYeme2(3,3),Yb3adjYw3mHw32(3,3),Yb3adjYw3MWM3(3,3),Yb3adjYw3BMWM3(3,3),           & 
& Ydmq2adjYu(3,3),YdadjYdTpYx3(3,3),YdadjYdTpTYx3(3,3),YdadjYumu2(3,3),YdadjTYdadjYx3(3,3),& 
& Yeml2adjYb3(3,3),Yeml2adjYw3(3,3),YeadjYb3MBM3(3,3),YeadjYb3mHb32(3,3),YeadjYb3BMBM3(3,3),& 
& YeadjYw3mHw32(3,3),YeadjYw3MWM3(3,3),YeadjYw3BMWM3(3,3),Yumq2adjYd(3,3),               & 
& YuadjYdmd2(3,3),Yw3ml2adjYb3(3,3),Yw3ml2adjYe(3,3),Yw3adjYb3MBM3(3,3),Yw3adjYb3mHb32(3,3),& 
& Yw3adjYb3BMBM3(3,3),Yw3adjYeme2(3,3),Yx3md2CYd(3,3),Yx3CYdmq2(3,3),adjYb3Yb3adjYb3(3,3),& 
& adjYb3Yb3adjYe(3,3),adjYb3Yb3adjYw3(3,3),adjYb3Yb3adjTYb3(3,3),adjYb3Yb3adjTYe(3,3),   & 
& adjYb3Yb3adjTYw3(3,3),adjYb3TYb3adjYb3(3,3),adjYb3TYb3adjYe(3,3),adjYb3TYb3adjYw3(3,3),& 
& adjYb3TYb3adjTYb3(3,3),adjYb3TYb3adjTYe(3,3),adjYb3TYb3adjTYw3(3,3),adjYdYdadjYd(3,3), & 
& adjYdYdadjYu(3,3),adjYdYdadjTYd(3,3),adjYdYdadjTYu(3,3),adjYdadjYx3Yx3(3,3),           & 
& adjYdTYdadjYd(3,3),adjYdTYdadjYu(3,3),adjYdTYdadjTYd(3,3),adjYdTYdadjTYu(3,3),         & 
& adjYdTpYx3CYx3(3,3),adjYdTpYx3CTYx3(3,3),adjYdTpTYx3CTYx3(3,3),adjYeYeadjYb3(3,3),     & 
& adjYeYeadjYe(3,3),adjYeYeadjYw3(3,3),adjYeYeadjTYb3(3,3),adjYeYeadjTYe(3,3),           & 
& adjYeYeadjTYw3(3,3),adjYeTYeadjYb3(3,3),adjYeTYeadjYe(3,3),adjYeTYeadjYw3(3,3),        & 
& adjYeTYeadjTYb3(3,3),adjYeTYeadjTYe(3,3),adjYeTYeadjTYw3(3,3),adjYuYuadjYd(3,3)

Complex(dp) :: adjYuYuadjYu(3,3),adjYuYuadjTYd(3,3),adjYuYuadjTYu(3,3),adjYuTYuadjYd(3,3),            & 
& adjYuTYuadjYu(3,3),adjYuTYuadjTYd(3,3),adjYuTYuadjTYu(3,3),adjYw3Yw3adjYb3(3,3),       & 
& adjYw3Yw3adjYe(3,3),adjYw3Yw3adjYw3(3,3),adjYw3Yw3adjTYb3(3,3),adjYw3Yw3adjTYe(3,3),   & 
& adjYw3Yw3adjTYw3(3,3),adjYw3TYw3adjYb3(3,3),adjYw3TYw3adjYe(3,3),adjYw3TYw3adjYw3(3,3),& 
& adjYw3TYw3adjTYb3(3,3),adjYw3TYw3adjTYe(3,3),adjYw3TYw3adjTYw3(3,3),adjYx3Yx3adjYx3(3,3),& 
& adjYx3Yx3adjTYx3(3,3),adjYx3Yx3CYd(3,3),adjYx3Yx3CTYd(3,3),adjYx3TYx3adjYx3(3,3),      & 
& adjYx3TYx3adjTYx3(3,3),adjYx3TYx3CTYd(3,3),adjTYb3TYb3adjYb3(3,3),adjTYb3TYb3adjYe(3,3),& 
& adjTYb3TYb3adjYw3(3,3),adjTYdTYdadjYd(3,3),adjTYdTYdadjYu(3,3),adjTYeTYeadjYb3(3,3),   & 
& adjTYeTYeadjYe(3,3),adjTYeTYeadjYw3(3,3),adjTYuTYuadjYd(3,3),adjTYuTYuadjYu(3,3),      & 
& adjTYw3TYw3adjYb3(3,3),adjTYw3TYw3adjYe(3,3),adjTYw3TYw3adjYw3(3,3),adjTYx3TYx3adjYx3(3,3),& 
& CYb3TpYb3CYb3(3,3),CYb3TpYb3CTYb3(3,3),CYb3TpYw3CYw3(3,3),CYb3TpYw3CTYw3(3,3),         & 
& CYb3TpTYb3CTYb3(3,3),CYb3TpTYw3CTYw3(3,3),CYdTpYdadjYx3(3,3),CYdTpYdadjTYx3(3,3),      & 
& CYdTpYdCYd(3,3),CYdTpYdCTYd(3,3),CYdTpTYdCTYd(3,3),CYeTpYeCYe(3,3),CYeTpYeCTYe(3,3),   & 
& CYeTpTYeCTYe(3,3),CYuTpYuCYu(3,3),CYuTpYuCTYu(3,3),CYuTpTYuCTYu(3,3),CYw3TpYb3CYb3(3,3),& 
& CYw3TpYb3CTYb3(3,3),CYw3TpYw3CYw3(3,3),CYw3TpYw3CTYw3(3,3),CYw3TpTYb3CTYb3(3,3),       & 
& CYw3TpTYw3CTYw3(3,3),CYx3YdadjYd(3,3),CYx3TpYx3CYx3(3,3),CYx3TpYx3CTYx3(3,3),          & 
& CYx3TpTYx3CTYx3(3,3),CTYb3TpTYb3CYb3(3,3),CTYb3TpTYw3CYw3(3,3),CTYdTpTYdadjYx3(3,3),   & 
& CTYdTpTYdCYd(3,3),CTYeTpTYeCYe(3,3),CTYuTpTYuCYu(3,3),CTYw3TpTYb3CYb3(3,3),            & 
& CTYw3TpTYw3CYw3(3,3),CTYx3TpTYx3CYx3(3,3),TYb3adjYw3MWM3(3,3),TYb3TpYb3CYb3(3,3),      & 
& TYb3TpYw3CYw3(3,3),TYdadjYdadjTYx3(3,3),TYdadjYdTpYx3(3,3),TYdTpYdCYd(3,3),            & 
& TYeadjYb3MBM3(3,3),TYeadjYw3MWM3(3,3),TYeTpYeCYe(3,3),TYuTpYuCYu(3,3),TYw3adjYb3MBM3(3,3),& 
& TYw3TpYb3CYb3(3,3),TYw3TpYw3CYw3(3,3),TYx3CTYdTpYd(3,3),TYx3TpYx3CYx3(3,3),            & 
& TpYb3CYb3TpYb3(3,3),TpYb3CYb3TpYw3(3,3),TpYb3CYb3TpTYb3(3,3),TpYb3CYb3TpTYw3(3,3),     & 
& TpYb3CTYb3adjYb3(3,3),TpYb3CTYb3adjYe(3,3),TpYb3CTYb3adjYw3(3,3),TpYdmd2adjYx3(3,3),   & 
& TpYdadjYx3mHxb32(3,3),TpYdadjYx3Yx3(3,3),TpYdadjYx3TYx3(3,3),TpYdCYdTpYd(3,3),         & 
& TpYdCYdTpTYd(3,3),TpYdCTYdadjYd(3,3),TpYdCTYdadjYu(3,3),TpYeCYeTpYb3(3,3),             & 
& TpYeCYeTpYw3(3,3),TpYeCYeTpTYb3(3,3),TpYeCYeTpTYw3(3,3),TpYeCTYeadjYb3(3,3),           & 
& TpYeCTYeadjYe(3,3),TpYeCTYeadjYw3(3,3),TpYuCYuTpYd(3,3),TpYuCYuTpTYd(3,3),             & 
& TpYuCTYuadjYd(3,3),TpYuCTYuadjYu(3,3),TpYw3CYw3TpYb3(3,3),TpYw3CYw3TpYw3(3,3),         & 
& TpYw3CYw3TpTYb3(3,3),TpYw3CYw3TpTYw3(3,3),TpYw3CTYw3adjYb3(3,3),TpYw3CTYw3adjYe(3,3),  & 
& TpYw3CTYw3adjYw3(3,3),TpYx3CYx3TpYx3(3,3),TpYx3CYx3TpTYx3(3,3),TpYx3CTYx3adjYx3(3,3),  & 
& TpYx3CTYx3TYd(3,3),TpYx3CTYx3TpTYd(3,3),TpTYb3CYb3TpYb3(3,3),TpTYb3CYb3TpYw3(3,3),     & 
& TpTYb3CTYb3adjYb3(3,3),TpTYb3CTYb3adjYe(3,3),TpTYb3CTYb3adjYw3(3,3),TpTYdCYdTpYd(3,3), & 
& TpTYdCTYdadjYd(3,3),TpTYdCTYdadjYu(3,3),TpTYeCYeTpYb3(3,3),TpTYeCYeTpYw3(3,3),         & 
& TpTYeCTYeadjYb3(3,3),TpTYeCTYeadjYe(3,3),TpTYeCTYeadjYw3(3,3),TpTYuCYuTpYd(3,3),       & 
& TpTYuCTYuadjYd(3,3),TpTYuCTYuadjYu(3,3),TpTYw3CYw3TpYb3(3,3),TpTYw3CYw3TpYw3(3,3),     & 
& TpTYw3CTYw3adjYb3(3,3),TpTYw3CTYw3adjYe(3,3),TpTYw3CTYw3adjYw3(3,3),TpTYx3CYx3TpYx3(3,3)

Complex(dp) :: TpTYx3CTYx3adjYx3(3,3),TpTYx3CTYx3TpYd(3,3),md2adjYx3Yx3adjYx3(3,3),md2adjYx3Yx3CYd(3,3),& 
& md2CYdTpYdadjYx3(3,3),md2CYdTpYdCYd(3,3),me2CYeTpYeCYe(3,3),mHb32CYb3TpYb3CYb3(3,3),   & 
& mHb32CYb3TpYw3CYw3(3,3),mHw32CYw3TpYb3CYb3(3,3),mHw32CYw3TpYw3CYw3(3,3),               & 
& mHxb32CYx3TpYx3CYx3(3,3),ml2adjYb3Yb3adjYb3(3,3),ml2adjYb3Yb3adjYe(3,3),               & 
& ml2adjYb3Yb3adjYw3(3,3),ml2adjYeYeadjYb3(3,3),ml2adjYeYeadjYe(3,3),ml2adjYeYeadjYw3(3,3),& 
& ml2adjYw3Yw3adjYb3(3,3),ml2adjYw3Yw3adjYe(3,3),ml2adjYw3Yw3adjYw3(3,3),mq2adjYdYdadjYd(3,3),& 
& mq2adjYdYdadjYu(3,3),mq2adjYuYuadjYd(3,3),mq2adjYuYuadjYu(3,3),mq2adjYx3Yx3CYd(3,3),   & 
& mu2CYuTpYuCYu(3,3),Yb3adjYb3Yb3adjYb3(3,3),Yb3adjYb3TYb3adjYb3(3,3),Yb3adjYb3TYb3adjTYb3(3,3),& 
& Yb3adjYeYeadjYb3(3,3),Yb3adjYeTYeadjYb3(3,3),Yb3adjYeTYeadjTYb3(3,3),Yb3adjYw3Yw3adjYb3(3,3),& 
& Yb3adjYw3TYw3adjYb3(3,3),Yb3adjYw3TYw3adjTYb3(3,3),Yb3adjTYb3TYb3adjYb3(3,3),          & 
& Yb3adjTYeTYeadjYb3(3,3),Yb3adjTYw3TYw3adjYb3(3,3),Yb3TpTYb3CTYb3adjYb3(3,3),           & 
& Yb3TpTYeCTYeadjYb3(3,3),Yb3TpTYw3CTYw3adjYb3(3,3),YdadjYdYdadjYd(3,3),YdadjYdTYdadjYd(3,3),& 
& YdadjYdTYdadjTYd(3,3),YdadjYdTpYx3CYx3(3,3),YdadjYdTpTYx3CTYx3(3,3),YdadjYuYuadjYd(3,3),& 
& YdadjYuTYuadjYd(3,3),YdadjYuTYuadjTYd(3,3),YdadjTYdTYdadjYd(3,3),YdadjTYuTYuadjYd(3,3),& 
& YdTpTYdCTYdadjYd(3,3),YdTpTYuCTYuadjYd(3,3),YeadjYb3Yb3adjYe(3,3),YeadjYb3TYb3adjYe(3,3),& 
& YeadjYb3TYb3adjTYe(3,3),YeadjYeYeadjYe(3,3),YeadjYeTYeadjYe(3,3),YeadjYeTYeadjTYe(3,3),& 
& YeadjYw3Yw3adjYe(3,3),YeadjYw3TYw3adjYe(3,3),YeadjYw3TYw3adjTYe(3,3),YeadjTYb3TYb3adjYe(3,3),& 
& YeadjTYeTYeadjYe(3,3),YeadjTYw3TYw3adjYe(3,3),YeTpTYb3CTYb3adjYe(3,3),YeTpTYeCTYeadjYe(3,3),& 
& YeTpTYw3CTYw3adjYe(3,3),YuadjYdYdadjYu(3,3),YuadjYdTYdadjYu(3,3),YuadjYdTYdadjTYu(3,3),& 
& YuadjYuYuadjYu(3,3),YuadjYuTYuadjYu(3,3),YuadjYuTYuadjTYu(3,3),YuadjTYdTYdadjYu(3,3),  & 
& YuadjTYuTYuadjYu(3,3),YuTpTYdCTYdadjYu(3,3),YuTpTYuCTYuadjYu(3,3),Yw3adjYb3Yb3adjYw3(3,3),& 
& Yw3adjYb3TYb3adjYw3(3,3),Yw3adjYb3TYb3adjTYw3(3,3),Yw3adjYeYeadjYw3(3,3),              & 
& Yw3adjYeTYeadjYw3(3,3),Yw3adjYeTYeadjTYw3(3,3),Yw3adjYw3Yw3adjYw3(3,3),Yw3adjYw3TYw3adjYw3(3,3),& 
& Yw3adjYw3TYw3adjTYw3(3,3),Yw3adjTYb3TYb3adjYw3(3,3),Yw3adjTYeTYeadjYw3(3,3),           & 
& Yw3adjTYw3TYw3adjYw3(3,3),Yw3TpTYb3CTYb3adjYw3(3,3),Yw3TpTYeCTYeadjYw3(3,3),           & 
& Yw3TpTYw3CTYw3adjYw3(3,3),Yx3adjYx3Yx3adjYx3(3,3),Yx3adjYx3TYx3adjYx3(3,3),            & 
& Yx3adjYx3TYx3adjTYx3(3,3),Yx3adjTYx3TYx3adjYx3(3,3),Yx3CYdTpYdadjYx3(3,3),             & 
& Yx3CTYdTpTYdadjYx3(3,3),Yx3TYdadjYdadjTYx3(3,3),Yx3TpTYx3CTYx3adjYx3(3,3),             & 
& adjYb3mHb32Yb3adjYb3(3,3),adjYb3mHb32Yb3adjYe(3,3),adjYb3mHb32Yb3adjYw3(3,3),          & 
& adjYb3Yb3ml2adjYb3(3,3),adjYb3Yb3ml2adjYe(3,3),adjYb3Yb3ml2adjYw3(3,3),adjYb3Yb3adjYb3MBM3(3,3),& 
& adjYb3Yb3adjYb3mHb32(3,3),adjYb3Yb3adjYb3Yb3(3,3),adjYb3Yb3adjYb3BMBM3(3,3),           & 
& adjYb3Yb3adjYb3TYb3(3,3),adjYb3Yb3adjYeme2(3,3),adjYb3Yb3adjYeYe(3,3),adjYb3Yb3adjYeTYe(3,3),& 
& adjYb3Yb3adjYw3mHw32(3,3),adjYb3Yb3adjYw3MWM3(3,3),adjYb3Yb3adjYw3Yw3(3,3),            & 
& adjYb3Yb3adjYw3BMWM3(3,3),adjYb3Yb3adjYw3TYw3(3,3),adjYb3TYb3adjYb3MBM3(3,3),          & 
& adjYb3TYb3adjYb3Yb3(3,3),adjYb3TYb3adjYeYe(3,3),adjYb3TYb3adjYw3MWM3(3,3),             & 
& adjYb3TYb3adjYw3Yw3(3,3),adjYdmd2YdadjYd(3,3),adjYdmd2YdadjYu(3,3),adjYdYdmq2adjYd(3,3),& 
& adjYdYdmq2adjYu(3,3),adjYdYdadjYdmd2(3,3),adjYdYdadjYdYd(3,3),adjYdYdadjYdTYd(3,3)

Complex(dp) :: adjYdYdadjYumu2(3,3),adjYdYdadjYuYu(3,3),adjYdYdadjYuTYu(3,3),adjYdTYdadjYdYd(3,3),    & 
& adjYdTYdadjYuYu(3,3),adjYdTpYx3CYx3Yd(3,3),adjYdTpYx3CYx3TYd(3,3),adjYdTpYx3CTYx3TYd(3,3),& 
& adjYdTpYx3CTYx3TpTYd(3,3),adjYdTpTYx3CYx3Yd(3,3),adjYdTpTYx3CTYx3TpYd(3,3),            & 
& adjYeme2YeadjYb3(3,3),adjYeme2YeadjYe(3,3),adjYeme2YeadjYw3(3,3),adjYeYeml2adjYb3(3,3),& 
& adjYeYeml2adjYe(3,3),adjYeYeml2adjYw3(3,3),adjYeYeadjYb3MBM3(3,3),adjYeYeadjYb3mHb32(3,3),& 
& adjYeYeadjYb3Yb3(3,3),adjYeYeadjYb3BMBM3(3,3),adjYeYeadjYb3TYb3(3,3),adjYeYeadjYeme2(3,3),& 
& adjYeYeadjYeYe(3,3),adjYeYeadjYeTYe(3,3),adjYeYeadjYw3mHw32(3,3),adjYeYeadjYw3MWM3(3,3),& 
& adjYeYeadjYw3Yw3(3,3),adjYeYeadjYw3BMWM3(3,3),adjYeYeadjYw3TYw3(3,3),adjYeTYeadjYb3MBM3(3,3),& 
& adjYeTYeadjYb3Yb3(3,3),adjYeTYeadjYeYe(3,3),adjYeTYeadjYw3MWM3(3,3),adjYeTYeadjYw3Yw3(3,3),& 
& adjYumu2YuadjYd(3,3),adjYumu2YuadjYu(3,3),adjYuYumq2adjYd(3,3),adjYuYumq2adjYu(3,3),   & 
& adjYuYuadjYdmd2(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdTYd(3,3),adjYuYuadjYumu2(3,3),    & 
& adjYuYuadjYuYu(3,3),adjYuYuadjYuTYu(3,3),adjYuTYuadjYdYd(3,3),adjYuTYuadjYuYu(3,3),    & 
& adjYw3mHw32Yw3adjYb3(3,3),adjYw3mHw32Yw3adjYe(3,3),adjYw3mHw32Yw3adjYw3(3,3),          & 
& adjYw3Yw3ml2adjYb3(3,3),adjYw3Yw3ml2adjYe(3,3),adjYw3Yw3ml2adjYw3(3,3),adjYw3Yw3adjYb3MBM3(3,3),& 
& adjYw3Yw3adjYb3mHb32(3,3),adjYw3Yw3adjYb3Yb3(3,3),adjYw3Yw3adjYb3BMBM3(3,3),           & 
& adjYw3Yw3adjYb3TYb3(3,3),adjYw3Yw3adjYeme2(3,3),adjYw3Yw3adjYeYe(3,3),adjYw3Yw3adjYeTYe(3,3),& 
& adjYw3Yw3adjYw3mHw32(3,3),adjYw3Yw3adjYw3MWM3(3,3),adjYw3Yw3adjYw3Yw3(3,3),            & 
& adjYw3Yw3adjYw3BMWM3(3,3),adjYw3Yw3adjYw3TYw3(3,3),adjYw3TYw3adjYb3MBM3(3,3),          & 
& adjYw3TYw3adjYb3Yb3(3,3),adjYw3TYw3adjYeYe(3,3),adjYw3TYw3adjYw3MWM3(3,3),             & 
& adjYw3TYw3adjYw3Yw3(3,3),adjYx3mHxb32Yx3adjYx3(3,3),adjYx3mHxb32Yx3CYd(3,3),           & 
& adjYx3Yx3md2adjYx3(3,3),adjYx3Yx3md2CYd(3,3),adjYx3Yx3adjYx3mHxb32(3,3),               & 
& adjYx3Yx3adjYx3Yx3(3,3),adjYx3Yx3adjYx3TYx3(3,3),adjYx3Yx3CYdmq2(3,3),adjYx3TYx3adjYx3Yx3(3,3),& 
& adjYx3TYx3CYdTpYd(3,3),adjYx3TYx3CTYdTpYd(3,3),adjTYb3TYb3TpYb3CYb3(3,3),              & 
& adjTYb3TYb3TpYw3CYw3(3,3),adjTYdTYdTpYdCYd(3,3),adjTYeTYeTpYeCYe(3,3),adjTYuTYuTpYuCYu(3,3),& 
& adjTYw3TYw3TpYb3CYb3(3,3),adjTYw3TYw3TpYw3CYw3(3,3),adjTYx3TYx3TpYx3CYx3(3,3),         & 
& CYb3ml2TpYb3CYb3(3,3),CYb3ml2TpYw3CYw3(3,3),CYb3TpYb3mHb32CYb3(3,3),CYb3TpYb3CYb3ml2(3,3),& 
& CYb3TpYb3CYb3TpYb3(3,3),CYb3TpYb3CYb3TpTYb3(3,3),CYb3TpYeCYeTpYb3(3,3),CYb3TpYeCYeTpTYb3(3,3),& 
& CYb3TpYw3mHw32CYw3(3,3),CYb3TpYw3CYw3ml2(3,3),CYb3TpYw3CYw3TpYb3(3,3),CYb3TpYw3CYw3TpTYb3(3,3),& 
& CYb3TpTYb3CYb3TpYb3(3,3),CYb3TpTYeCYeTpYb3(3,3),CYb3TpTYw3CYw3TpYb3(3,3),              & 
& CYdmq2TpYdadjYx3(3,3),CYdmq2TpYdCYd(3,3),CYdTpYdmd2adjYx3(3,3),CYdTpYdmd2CYd(3,3),     & 
& CYdTpYdadjYx3mHxb32(3,3),CYdTpYdadjYx3Yx3(3,3),CYdTpYdadjYx3TYx3(3,3),CYdTpYdCYdmq2(3,3),& 
& CYdTpYdCYdTpYd(3,3),CYdTpYdCYdTpTYd(3,3),CYdTpYuCYuTpYd(3,3),CYdTpYuCYuTpTYd(3,3),     & 
& CYdTpTYdCYdTpYd(3,3),CYdTpTYuCYuTpYd(3,3),CYeml2TpYeCYe(3,3),CYeTpYeme2CYe(3,3),       & 
& CYeTpYeCYeml2(3,3),CYumq2TpYuCYu(3,3),CYuTpYumu2CYu(3,3),CYuTpYuCYumq2(3,3),           & 
& CYw3ml2TpYb3CYb3(3,3),CYw3ml2TpYw3CYw3(3,3),CYw3TpYb3mHb32CYb3(3,3),CYw3TpYb3CYb3ml2(3,3),& 
& CYw3TpYb3CYb3TpYw3(3,3),CYw3TpYb3CYb3TpTYw3(3,3),CYw3TpYeCYeTpYw3(3,3),CYw3TpYeCYeTpTYw3(3,3),& 
& CYw3TpYw3mHw32CYw3(3,3),CYw3TpYw3CYw3ml2(3,3),CYw3TpYw3CYw3TpYw3(3,3),CYw3TpYw3CYw3TpTYw3(3,3)

Complex(dp) :: CYw3TpTYb3CYb3TpYw3(3,3),CYw3TpTYeCYeTpYw3(3,3),CYw3TpTYw3CYw3TpYw3(3,3),              & 
& CYx3md2TpYx3CYx3(3,3),CYx3YdadjYdTpYx3(3,3),CYx3YdadjYdTpTYx3(3,3),CYx3TYdadjYdTpYx3(3,3),& 
& CYx3TpYx3mHxb32CYx3(3,3),CYx3TpYx3CYx3md2(3,3),CYx3TpYx3CYx3Yd(3,3),CYx3TpYx3CYx3TYd(3,3),& 
& CYx3TpYx3CYx3TpYx3(3,3),CYx3TpYx3CYx3TpTYx3(3,3),CYx3TpTYx3CYx3Yd(3,3),CYx3TpTYx3CYx3TpYx3(3,3),& 
& TYb3adjYb3Yb3adjTYb3(3,3),TYb3adjYeYeadjTYb3(3,3),TYb3adjYw3Yw3adjTYb3(3,3),           & 
& TYb3TpYb3CTYb3adjYb3(3,3),TYb3TpYeCTYeadjYb3(3,3),TYb3TpYw3CTYw3adjYb3(3,3),           & 
& TYdadjYdYdadjTYd(3,3),TYdadjYdadjYx3Yx3(3,3),TYdadjYuYuadjTYd(3,3),TYdTpYdCTYdadjYd(3,3),& 
& TYdTpYuCTYuadjYd(3,3),TYeadjYb3Yb3adjTYe(3,3),TYeadjYeYeadjTYe(3,3),TYeadjYw3Yw3adjTYe(3,3),& 
& TYeTpYb3CTYb3adjYe(3,3),TYeTpYeCTYeadjYe(3,3),TYeTpYw3CTYw3adjYe(3,3),TYuadjYdYdadjTYu(3,3),& 
& TYuadjYuYuadjTYu(3,3),TYuTpYdCTYdadjYu(3,3),TYuTpYuCTYuadjYu(3,3),TYw3adjYb3Yb3adjTYw3(3,3),& 
& TYw3adjYeYeadjTYw3(3,3),TYw3adjYw3Yw3adjTYw3(3,3),TYw3TpYb3CTYb3adjYw3(3,3),           & 
& TYw3TpYeCTYeadjYw3(3,3),TYw3TpYw3CTYw3adjYw3(3,3),TYx3YdadjTYdadjYx3(3,3),             & 
& TYx3adjYx3Yx3adjTYx3(3,3),TYx3CYdTpYdadjTYx3(3,3),TYx3TpYx3CTYx3adjYx3(3,3),           & 
& TpYb3CYb3TpYb3CYb3(3,3),TpYb3CYb3TpYw3CYw3(3,3),TpYb3CYb3TpTYb3CTYb3(3,3),             & 
& TpYb3CYb3TpTYw3CTYw3(3,3),TpYb3CTYb3TpTYb3CYb3(3,3),TpYb3CTYb3TpTYw3CYw3(3,3),         & 
& TpYdadjYdTpTYx3CTYx3(3,3),TpYdadjYx3Yx3CYd(3,3),TpYdadjYx3TYx3CTYd(3,3),               & 
& TpYdCYdTpYdCYd(3,3),TpYdCYdTpTYdCTYd(3,3),TpYdCTYdTpTYdCYd(3,3),TpYeCYeTpYeCYe(3,3),   & 
& TpYeCYeTpTYeCTYe(3,3),TpYeCTYeTpTYeCYe(3,3),TpYuCYuTpYuCYu(3,3),TpYuCYuTpTYuCTYu(3,3), & 
& TpYuCTYuTpTYuCYu(3,3),TpYw3CYw3TpYb3CYb3(3,3),TpYw3CYw3TpYw3CYw3(3,3),TpYw3CYw3TpTYb3CTYb3(3,3),& 
& TpYw3CYw3TpTYw3CTYw3(3,3),TpYw3CTYw3TpTYb3CYb3(3,3),TpYw3CTYw3TpTYw3CYw3(3,3),         & 
& TpYx3CYx3YdadjYd(3,3),TpYx3CYx3TpYx3CYx3(3,3),TpYx3CYx3TpTYx3CTYx3(3,3),               & 
& TpYx3CTYx3TpTYx3CYx3(3,3),TpTYb3CYb3TpYb3CTYb3(3,3),TpTYb3CYb3TpYw3CTYw3(3,3),         & 
& TpTYdadjYdTpYx3CTYx3(3,3),TpTYdadjYx3Yx3CTYd(3,3),TpTYdCYdTpYdCTYd(3,3),               & 
& TpTYeCYeTpYeCTYe(3,3),TpTYuCYuTpYuCTYu(3,3),TpTYw3CYw3TpYb3CTYb3(3,3),TpTYw3CYw3TpYw3CTYw3(3,3),& 
& TpTYx3CYx3TpYx3CTYx3(3,3),MBM3CYb3TpYb3CYb3TpYb3(3,3),MBM3CYb3TpYb3CYb3TpTYb3(3,3),    & 
& MBM3CYb3TpYeCYeTpYb3(3,3),MBM3CYb3TpYeCYeTpTYb3(3,3),MBM3CYb3TpYw3CYw3TpYb3(3,3),      & 
& MBM3CYb3TpYw3CYw3TpTYb3(3,3),MBM3CYb3TpTYb3CYb3TpYb3(3,3),MBM3CYb3TpTYeCYeTpYb3(3,3),  & 
& MBM3CYb3TpTYw3CYw3TpYb3(3,3),md2YdadjYdYdadjYd(3,3),md2YdadjYdTpYx3CYx3(3,3),          & 
& md2YdadjYuYuadjYd(3,3),md2adjYx3Yx3adjYx3Yx3(3,3),md2TpYx3CYx3YdadjYd(3,3),            & 
& md2TpYx3CYx3TpYx3CYx3(3,3),me2YeadjYb3Yb3adjYe(3,3),me2YeadjYeYeadjYe(3,3),            & 
& me2YeadjYw3Yw3adjYe(3,3),mHb32Yb3adjYb3Yb3adjYb3(3,3),mHb32Yb3adjYeYeadjYb3(3,3),      & 
& mHb32Yb3adjYw3Yw3adjYb3(3,3),mHw32Yw3adjYb3Yb3adjYw3(3,3),mHw32Yw3adjYeYeadjYw3(3,3),  & 
& mHw32Yw3adjYw3Yw3adjYw3(3,3),mHxb32Yx3adjYx3Yx3adjYx3(3,3),mHxb32Yx3CYdTpYdadjYx3(3,3),& 
& mHxb32CYx3YdadjYdTpYx3(3,3),ml2adjYb3Yb3adjYb3Yb3(3,3),ml2adjYb3Yb3adjYeYe(3,3),       & 
& ml2adjYb3Yb3adjYw3Yw3(3,3),ml2adjYeYeadjYb3Yb3(3,3),ml2adjYeYeadjYeYe(3,3),            & 
& ml2adjYeYeadjYw3Yw3(3,3),ml2adjYw3Yw3adjYb3Yb3(3,3),ml2adjYw3Yw3adjYeYe(3,3),          & 
& ml2adjYw3Yw3adjYw3Yw3(3,3),ml2TpYb3CYb3TpYb3CYb3(3,3),ml2TpYb3CYb3TpYw3CYw3(3,3)

Complex(dp) :: ml2TpYeCYeTpYeCYe(3,3),ml2TpYw3CYw3TpYb3CYb3(3,3),ml2TpYw3CYw3TpYw3CYw3(3,3),          & 
& mq2adjYdYdadjYdYd(3,3),mq2adjYdYdadjYuYu(3,3),mq2adjYdTpYx3CYx3Yd(3,3),mq2adjYuYuadjYdYd(3,3),& 
& mq2adjYuYuadjYuYu(3,3),mq2TpYdCYdTpYdCYd(3,3),mq2TpYuCYuTpYuCYu(3,3),mu2YuadjYdYdadjYu(3,3),& 
& mu2YuadjYuYuadjYu(3,3),MWM3CYw3TpYb3CYb3TpYw3(3,3),MWM3CYw3TpYb3CYb3TpTYw3(3,3),       & 
& MWM3CYw3TpYeCYeTpYw3(3,3),MWM3CYw3TpYeCYeTpTYw3(3,3),MWM3CYw3TpYw3CYw3TpYw3(3,3),      & 
& MWM3CYw3TpYw3CYw3TpTYw3(3,3),MWM3CYw3TpTYb3CYb3TpYw3(3,3),MWM3CYw3TpTYeCYeTpYw3(3,3),  & 
& MWM3CYw3TpTYw3CYw3TpYw3(3,3),MXM3CYx3YdadjYdTpYx3(3,3),MXM3CYx3YdadjYdTpTYx3(3,3),     & 
& MXM3CYx3TYdadjYdTpYx3(3,3),MXM3CYx3TpYx3CYx3TpYx3(3,3),MXM3CYx3TpYx3CYx3TpTYx3(3,3),   & 
& MXM3CYx3TpTYx3CYx3TpYx3(3,3),Yb3ml2adjYb3Yb3adjYb3(3,3),Yb3ml2adjYeYeadjYb3(3,3),      & 
& Yb3ml2adjYw3Yw3adjYb3(3,3),Yb3adjYb3mHb32Yb3adjYb3(3,3),Yb3adjYb3Yb3ml2adjYb3(3,3),    & 
& Yb3adjYb3Yb3adjYb3MBM3(3,3),Yb3adjYb3Yb3adjYb3mHb32(3,3),Yb3adjYb3Yb3adjYb3Yb3(3,3),   & 
& Yb3adjYb3Yb3adjYb3BMBM3(3,3),Yb3adjYb3Yb3adjYb3TYb3(3,3),Yb3adjYb3Yb3adjYw3Yw3(3,3),   & 
& Yb3adjYb3Yb3adjYw3TYw3(3,3),Yb3adjYb3TYb3adjYb3MBM3(3,3),Yb3adjYb3TYb3adjYb3Yb3(3,3),  & 
& Yb3adjYb3TYb3adjYw3Yw3(3,3),Yb3adjYeme2YeadjYb3(3,3),Yb3adjYeYeml2adjYb3(3,3),         & 
& Yb3adjYeYeadjYb3MBM3(3,3),Yb3adjYeYeadjYb3mHb32(3,3),Yb3adjYeYeadjYb3Yb3(3,3),         & 
& Yb3adjYeYeadjYb3BMBM3(3,3),Yb3adjYeYeadjYb3TYb3(3,3),Yb3adjYeYeadjYeYe(3,3),           & 
& Yb3adjYeYeadjYeTYe(3,3),Yb3adjYeYeadjYw3TYw3(3,3),Yb3adjYeTYeadjYb3MBM3(3,3),          & 
& Yb3adjYeTYeadjYb3Yb3(3,3),Yb3adjYeTYeadjYeYe(3,3),Yb3adjYeTYeadjYw3Yw3(3,3),           & 
& Yb3adjYw3mHw32Yw3adjYb3(3,3),Yb3adjYw3Yw3ml2adjYb3(3,3),Yb3adjYw3Yw3adjYb3MBM3(3,3),   & 
& Yb3adjYw3Yw3adjYb3mHb32(3,3),Yb3adjYw3Yw3adjYb3Yb3(3,3),Yb3adjYw3Yw3adjYb3BMBM3(3,3),  & 
& Yb3adjYw3Yw3adjYb3TYb3(3,3),Yb3adjYw3Yw3adjYw3Yw3(3,3),Yb3adjYw3Yw3adjYw3TYw3(3,3),    & 
& Yb3adjYw3TYw3adjYb3MBM3(3,3),Yb3adjYw3TYw3adjYb3Yb3(3,3),Yb3adjYw3TYw3adjYw3Yw3(3,3),  & 
& Ydmq2adjYdYdadjYd(3,3),Ydmq2adjYuYuadjYd(3,3),Ydmq2adjYx3Yx3CYd(3,3),YdadjYdmd2YdadjYd(3,3),& 
& YdadjYdYdmq2adjYd(3,3),YdadjYdYdadjYdmd2(3,3),YdadjYdYdadjYdYd(3,3),YdadjYdYdadjYdTYd(3,3),& 
& YdadjYdTYdadjYdYd(3,3),YdadjYdTpYx3CYx3Yd(3,3),YdadjYdTpYx3CYx3TYd(3,3),               & 
& YdadjYdTpTYx3CYx3Yd(3,3),YdadjYumu2YuadjYd(3,3),YdadjYuYumq2adjYd(3,3),YdadjYuYuadjYdmd2(3,3),& 
& YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYdTYd(3,3),YdadjYuYuadjYuYu(3,3),YdadjYuYuadjYuTYu(3,3),& 
& YdadjYuTYuadjYdYd(3,3),YdadjYuTYuadjYuYu(3,3),Yeml2adjYb3Yb3adjYe(3,3),Yeml2adjYeYeadjYe(3,3),& 
& Yeml2adjYw3Yw3adjYe(3,3),YeadjYb3mHb32Yb3adjYe(3,3),YeadjYb3Yb3ml2adjYe(3,3),          & 
& YeadjYb3Yb3adjYb3Yb3(3,3),YeadjYb3Yb3adjYb3TYb3(3,3),YeadjYb3Yb3adjYeme2(3,3),         & 
& YeadjYb3Yb3adjYeYe(3,3),YeadjYb3Yb3adjYeTYe(3,3),YeadjYb3Yb3adjYw3Yw3(3,3),            & 
& YeadjYb3Yb3adjYw3TYw3(3,3),YeadjYb3TYb3adjYb3Yb3(3,3),YeadjYb3TYb3adjYeYe(3,3),        & 
& YeadjYb3TYb3adjYw3Yw3(3,3),YeadjYeme2YeadjYe(3,3),YeadjYeYeml2adjYe(3,3),              & 
& YeadjYeYeadjYeme2(3,3),YeadjYeYeadjYeYe(3,3),YeadjYeYeadjYeTYe(3,3),YeadjYeTYeadjYeYe(3,3),& 
& YeadjYw3mHw32Yw3adjYe(3,3),YeadjYw3Yw3ml2adjYe(3,3),YeadjYw3Yw3adjYb3Yb3(3,3),         & 
& YeadjYw3Yw3adjYb3TYb3(3,3),YeadjYw3Yw3adjYeme2(3,3),YeadjYw3Yw3adjYeYe(3,3),           & 
& YeadjYw3Yw3adjYeTYe(3,3),YeadjYw3Yw3adjYw3Yw3(3,3),YeadjYw3Yw3adjYw3TYw3(3,3)

Complex(dp) :: YeadjYw3TYw3adjYb3Yb3(3,3),YeadjYw3TYw3adjYeYe(3,3),YeadjYw3TYw3adjYw3Yw3(3,3),        & 
& Yumq2adjYdYdadjYu(3,3),Yumq2adjYuYuadjYu(3,3),YuadjYdmd2YdadjYu(3,3),YuadjYdYdmq2adjYu(3,3),& 
& YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYdTYd(3,3),YuadjYdYdadjYumu2(3,3),YuadjYdYdadjYuYu(3,3),& 
& YuadjYdYdadjYuTYu(3,3),YuadjYdTYdadjYdYd(3,3),YuadjYdTYdadjYuYu(3,3),YuadjYdTpYx3CYx3Yd(3,3),& 
& YuadjYdTpYx3CYx3TYd(3,3),YuadjYdTpTYx3CYx3Yd(3,3),YuadjYumu2YuadjYu(3,3),              & 
& YuadjYuYumq2adjYu(3,3),YuadjYuYuadjYumu2(3,3),YuadjYuYuadjYuYu(3,3),YuadjYuYuadjYuTYu(3,3),& 
& YuadjYuTYuadjYuYu(3,3),Yw3ml2adjYb3Yb3adjYw3(3,3),Yw3ml2adjYeYeadjYw3(3,3),            & 
& Yw3ml2adjYw3Yw3adjYw3(3,3),Yw3adjYb3mHb32Yb3adjYw3(3,3),Yw3adjYb3Yb3ml2adjYw3(3,3),    & 
& Yw3adjYb3Yb3adjYb3Yb3(3,3),Yw3adjYb3Yb3adjYb3TYb3(3,3),Yw3adjYb3Yb3adjYw3mHw32(3,3),   & 
& Yw3adjYb3Yb3adjYw3MWM3(3,3),Yw3adjYb3Yb3adjYw3Yw3(3,3),Yw3adjYb3Yb3adjYw3BMWM3(3,3),   & 
& Yw3adjYb3Yb3adjYw3TYw3(3,3),Yw3adjYb3TYb3adjYb3Yb3(3,3),Yw3adjYb3TYb3adjYw3MWM3(3,3),  & 
& Yw3adjYb3TYb3adjYw3Yw3(3,3),Yw3adjYeme2YeadjYw3(3,3),Yw3adjYeYeml2adjYw3(3,3),         & 
& Yw3adjYeYeadjYb3TYb3(3,3),Yw3adjYeYeadjYeYe(3,3),Yw3adjYeYeadjYeTYe(3,3),              & 
& Yw3adjYeYeadjYw3mHw32(3,3),Yw3adjYeYeadjYw3MWM3(3,3),Yw3adjYeYeadjYw3Yw3(3,3),         & 
& Yw3adjYeYeadjYw3BMWM3(3,3),Yw3adjYeYeadjYw3TYw3(3,3),Yw3adjYeTYeadjYb3Yb3(3,3),        & 
& Yw3adjYeTYeadjYeYe(3,3),Yw3adjYeTYeadjYw3MWM3(3,3),Yw3adjYeTYeadjYw3Yw3(3,3),          & 
& Yw3adjYw3mHw32Yw3adjYw3(3,3),Yw3adjYw3Yw3ml2adjYw3(3,3),Yw3adjYw3Yw3adjYb3Yb3(3,3),    & 
& Yw3adjYw3Yw3adjYb3TYb3(3,3),Yw3adjYw3Yw3adjYw3mHw32(3,3),Yw3adjYw3Yw3adjYw3MWM3(3,3),  & 
& Yw3adjYw3Yw3adjYw3Yw3(3,3),Yw3adjYw3Yw3adjYw3BMWM3(3,3),Yw3adjYw3Yw3adjYw3TYw3(3,3),   & 
& Yw3adjYw3TYw3adjYb3Yb3(3,3),Yw3adjYw3TYw3adjYw3MWM3(3,3),Yw3adjYw3TYw3adjYw3Yw3(3,3),  & 
& Yx3md2adjYx3Yx3adjYx3(3,3),Yx3md2CYdTpYdadjYx3(3,3),Yx3adjYx3mHxb32Yx3adjYx3(3,3),     & 
& Yx3adjYx3Yx3md2adjYx3(3,3),Yx3adjYx3Yx3adjYx3mHxb32(3,3),Yx3adjYx3Yx3adjYx3Yx3(3,3),   & 
& Yx3adjYx3Yx3adjYx3TYx3(3,3),Yx3adjYx3TYx3adjYx3Yx3(3,3),Yx3CYdmq2TpYdadjYx3(3,3),      & 
& Yx3CYdTpYdmd2adjYx3(3,3),Yx3CYdTpYdadjYx3mHxb32(3,3),Yx3CYdTpYdadjYx3Yx3(3,3),         & 
& Yx3CYdTpYdadjYx3TYx3(3,3),Yx3CYdTpYdCYdTpYd(3,3),Yx3CYdTpYdCYdTpTYd(3,3),              & 
& Yx3CYdTpYuCYuTpYd(3,3),Yx3CYdTpYuCYuTpTYd(3,3),Yx3CYdTpTYdCYdTpYd(3,3),Yx3CYdTpTYuCYuTpYd(3,3),& 
& Yx3TYdadjYdadjYx3Yx3(3,3),BMBM3CYb3TpYb3CYb3TpYb3(3,3),BMBM3CYb3TpYeCYeTpYb3(3,3),     & 
& BMBM3CYb3TpYw3CYw3TpYb3(3,3),BMWM3CYw3TpYb3CYb3TpYw3(3,3),BMWM3CYw3TpYeCYeTpYw3(3,3),  & 
& BMWM3CYw3TpYw3CYw3TpYw3(3,3),BMXM3CYx3YdadjYdTpYx3(3,3),BMXM3CYx3TpYx3CYx3TpYx3(3,3),  & 
& TYb3adjYb3Yb3adjYb3MBM3(3,3),TYb3adjYb3Yb3adjYb3Yb3(3,3),TYb3adjYb3Yb3adjYw3Yw3(3,3),  & 
& TYb3adjYeYeadjYb3MBM3(3,3),TYb3adjYeYeadjYb3Yb3(3,3),TYb3adjYeYeadjYeYe(3,3),          & 
& TYb3adjYeYeadjYw3Yw3(3,3),TYb3adjYw3Yw3adjYb3MBM3(3,3),TYb3adjYw3Yw3adjYb3Yb3(3,3),    & 
& TYb3adjYw3Yw3adjYw3Yw3(3,3),TYdadjYdYdadjYdYd(3,3),TYdadjYdTpYx3CYx3Yd(3,3),           & 
& TYdadjYuYuadjYdYd(3,3),TYdadjYuYuadjYuYu(3,3),TYeadjYb3Yb3adjYb3Yb3(3,3),              & 
& TYeadjYb3Yb3adjYeYe(3,3),TYeadjYb3Yb3adjYw3Yw3(3,3),TYeadjYeYeadjYeYe(3,3),            & 
& TYeadjYw3Yw3adjYb3Yb3(3,3),TYeadjYw3Yw3adjYeYe(3,3),TYeadjYw3Yw3adjYw3Yw3(3,3),        & 
& TYuadjYdYdadjYdYd(3,3),TYuadjYdYdadjYuYu(3,3),TYuadjYdTpYx3CYx3Yd(3,3),TYuadjYuYuadjYuYu(3,3)

Complex(dp) :: TYw3adjYb3Yb3adjYb3Yb3(3,3),TYw3adjYb3Yb3adjYw3MWM3(3,3),TYw3adjYb3Yb3adjYw3Yw3(3,3),  & 
& TYw3adjYeYeadjYb3Yb3(3,3),TYw3adjYeYeadjYeYe(3,3),TYw3adjYeYeadjYw3MWM3(3,3),          & 
& TYw3adjYeYeadjYw3Yw3(3,3),TYw3adjYw3Yw3adjYb3Yb3(3,3),TYw3adjYw3Yw3adjYw3MWM3(3,3),    & 
& TYw3adjYw3Yw3adjYw3Yw3(3,3),TYx3adjYx3Yx3adjYx3Yx3(3,3),TYx3CYdTpYdadjYx3Yx3(3,3),     & 
& TYx3CYdTpYdCYdTpYd(3,3),TYx3CYdTpYuCYuTpYd(3,3),TpYb3mHb32CYb3TpYb3CYb3(3,3),          & 
& TpYb3mHb32CYb3TpYw3CYw3(3,3),TpYb3CYb3ml2TpYb3CYb3(3,3),TpYb3CYb3ml2TpYw3CYw3(3,3),    & 
& TpYb3CYb3TpYb3mHb32CYb3(3,3),TpYb3CYb3TpYb3CYb3ml2(3,3),TpYb3CYb3TpYw3mHw32CYw3(3,3),  & 
& TpYb3CYb3TpYw3CYw3ml2(3,3),TpYdmd2adjYx3Yx3CYd(3,3),TpYdmd2CYdTpYdCYd(3,3),            & 
& TpYdadjYx3mHxb32Yx3CYd(3,3),TpYdadjYx3Yx3md2CYd(3,3),TpYdadjYx3Yx3CYdmq2(3,3),         & 
& TpYdCYdmq2TpYdCYd(3,3),TpYdCYdTpYdmd2CYd(3,3),TpYdCYdTpYdCYdmq2(3,3),TpYeme2CYeTpYeCYe(3,3),& 
& TpYeCYeml2TpYeCYe(3,3),TpYeCYeTpYeme2CYe(3,3),TpYeCYeTpYeCYeml2(3,3),TpYumu2CYuTpYuCYu(3,3),& 
& TpYuCYumq2TpYuCYu(3,3),TpYuCYuTpYumu2CYu(3,3),TpYuCYuTpYuCYumq2(3,3),TpYw3mHw32CYw3TpYb3CYb3(3,3),& 
& TpYw3mHw32CYw3TpYw3CYw3(3,3),TpYw3CYw3ml2TpYb3CYb3(3,3),TpYw3CYw3ml2TpYw3CYw3(3,3),    & 
& TpYw3CYw3TpYb3mHb32CYb3(3,3),TpYw3CYw3TpYb3CYb3ml2(3,3),TpYw3CYw3TpYw3mHw32CYw3(3,3),  & 
& TpYw3CYw3TpYw3CYw3ml2(3,3),TpYx3mHxb32CYx3TpYx3CYx3(3,3),TpYx3CYx3md2TpYx3CYx3(3,3),   & 
& TpYx3CYx3TpYx3mHxb32CYx3(3,3),TpYx3CYx3TpYx3CYx3md2(3,3),TpYx3CYx3TpYx3CYx3Yd(3,3),    & 
& TpYx3CYx3TpYx3CYx3TYd(3,3),TpYx3CYx3TpTYx3CYx3Yd(3,3),TpTYx3CYx3TpYx3CYx3Yd(3,3)

Complex(dp) :: Trmd2,Trme2,TrmHx32,TrmHxb32,Trml2,Trmq2,Trmu2,TrYb3adjYb3,TrYdadjYd,TrYeadjYe,       & 
& TrYuadjYu,TrYw3adjYw3,TrYx3adjYx3,TradjYb3TYb3,TradjYdTYd,TradjYeTYe,TradjYuTYu,       & 
& TradjYw3TYw3,TradjYx3TYx3,TrCTYb3TpTYb3,TrCTYdTpTYd,TrCTYeTpTYe,TrCTYuTpTYu,           & 
& TrCTYw3TpTYw3,TrCTYx3TpTYx3,Trmd2YdadjYd,Trmd2adjYx3Yx3,Trme2YeadjYe,TrmHb32Yb3adjYb3, & 
& TrmHw32Yw3adjYw3,TrmHxb32Yx3adjYx3,Trml2adjYb3Yb3,Trml2adjYeYe,Trml2adjYw3Yw3,         & 
& Trmq2adjYdYd,Trmq2adjYuYu,Trmu2YuadjYu

Complex(dp) :: TrmHg32,TrmHw32,TrCTYb3TpYb3,TrCTYdTpYd,TrCTYeTpYe,TrCTYuTpYu,TrCTYw3TpYw3,           & 
& TrCTYx3TpYx3,TrYb3adjYb3Yb3adjYb3,TrYb3adjYb3TYb3adjYb3,TrYb3adjYb3TYb3adjTYb3,        & 
& TrYb3adjYeYeadjYb3,TrYb3adjYeTYeadjYb3,TrYb3adjYeTYeadjTYb3,TrYb3adjYw3Yw3adjYb3,      & 
& TrYb3adjYw3TYw3adjYb3,TrYb3adjYw3TYw3adjTYb3,TrYb3adjTYb3TYb3adjYb3,TrYb3adjTYeTYeadjYb3,& 
& TrYb3adjTYw3TYw3adjYb3,TrYb3TpTYb3CTYb3adjYb3,TrYb3TpTYeCTYeadjYb3,TrYb3TpTYw3CTYw3adjYb3,& 
& TrYdadjYdYdadjYd,TrYdadjYdTYdadjYd,TrYdadjYdTYdadjTYd,TrYdadjYdTpYx3CYx3,              & 
& TrYdadjYdTpTYx3CTYx3,TrYdadjYuYuadjYd,TrYdadjYuTYuadjYd,TrYdadjYuTYuadjTYd,            & 
& TrYdadjTYdTYdadjYd,TrYdadjTYuTYuadjYd,TrYdTpTYdCTYdadjYd,TrYdTpTYuCTYuadjYd,           & 
& TrYeadjYb3TYb3adjYe,TrYeadjYb3TYb3adjTYe,TrYeadjYeYeadjYe,TrYeadjYeTYeadjYe,           & 
& TrYeadjYeTYeadjTYe,TrYeadjYw3Yw3adjYe,TrYeadjYw3TYw3adjYe,TrYeadjYw3TYw3adjTYe,        & 
& TrYeadjTYb3TYb3adjYe,TrYeadjTYeTYeadjYe,TrYeadjTYw3TYw3adjYe,TrYeTpTYb3CTYb3adjYe,     & 
& TrYeTpTYeCTYeadjYe,TrYeTpTYw3CTYw3adjYe,TrYuadjYdTYdadjYu,TrYuadjYdTYdadjTYu,          & 
& TrYuadjYuYuadjYu,TrYuadjYuTYuadjYu,TrYuadjYuTYuadjTYu,TrYuadjTYdTYdadjYu,              & 
& TrYuadjTYuTYuadjYu,TrYuTpTYdCTYdadjYu,TrYuTpTYuCTYuadjYu,TrYw3adjYb3TYb3adjYw3,        & 
& TrYw3adjYb3TYb3adjTYw3,TrYw3adjYeTYeadjYw3,TrYw3adjYeTYeadjTYw3,TrYw3adjYw3Yw3adjYw3,  & 
& TrYw3adjYw3TYw3adjYw3,TrYw3adjYw3TYw3adjTYw3,TrYw3adjTYb3TYb3adjYw3,TrYw3adjTYeTYeadjYw3,& 
& TrYw3adjTYw3TYw3adjYw3,TrYw3TpTYb3CTYb3adjYw3,TrYw3TpTYeCTYeadjYw3,TrYw3TpTYw3CTYw3adjYw3,& 
& TrYx3adjYx3Yx3adjYx3,TrYx3adjYx3TYx3adjYx3,TrYx3adjYx3TYx3adjTYx3,TrYx3adjTYx3TYx3adjYx3,& 
& TrYx3CTYdTpTYdadjYx3,TrYx3TpTYx3CTYx3adjYx3,TradjYdTpYx3CYx3TYd,TradjYdTpYx3CTYx3TYd,  & 
& TradjYdTpYx3CTYx3TpTYd,TradjYdTpTYx3CTYx3TpYd,TradjYx3TYx3CYdTpYd,TradjYx3TYx3CTYdTpYd,& 
& Trmd2YdadjYdYdadjYd,Trmd2YdadjYdTpYx3CYx3,Trmd2YdadjYuYuadjYd,Trmd2adjYx3Yx3adjYx3Yx3, & 
& Trmd2TpYx3CYx3YdadjYd,Trme2YeadjYb3Yb3adjYe,Trme2YeadjYeYeadjYe,Trme2YeadjYw3Yw3adjYe, & 
& TrmHb32Yb3adjYb3Yb3adjYb3,TrmHb32Yb3adjYeYeadjYb3,TrmHb32Yb3adjYw3Yw3adjYb3,           & 
& TrmHw32Yw3adjYb3Yb3adjYw3,TrmHw32Yw3adjYeYeadjYw3,TrmHw32Yw3adjYw3Yw3adjYw3,           & 
& TrmHxb32Yx3adjYx3Yx3adjYx3,TrmHxb32Yx3CYdTpYdadjYx3,TrmHxb32CYx3YdadjYdTpYx3,          & 
& Trml2adjYb3Yb3adjYb3Yb3,Trml2adjYb3Yb3adjYeYe,Trml2adjYb3Yb3adjYw3Yw3,Trml2adjYeYeadjYb3Yb3,& 
& Trml2adjYeYeadjYeYe,Trml2adjYeYeadjYw3Yw3,Trml2adjYw3Yw3adjYb3Yb3,Trml2adjYw3Yw3adjYeYe,& 
& Trml2adjYw3Yw3adjYw3Yw3,Trmq2adjYdYdadjYdYd,Trmq2adjYdYdadjYuYu,Trmq2adjYdTpYx3CYx3Yd, & 
& Trmq2adjYuYuadjYdYd,Trmq2adjYuYuadjYuYu,Trmu2YuadjYdYdadjYu,Trmu2YuadjYuYuadjYu

Complex(dp) :: MnuL(3,3), DMnuL(3,3), betaMnuL1(3,3), betaMNuL2(3,3), gammaL1(3,3), &
 & gammaL2(3,3), gammaHd1, gammaHd2

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Real(dp) :: g1p4,g2p4,g3p4

Iname = Iname +1 
NameOfUnit(Iname) = 'rge555' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters555(gy,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3,mue,MXM3,MWM3,MGM3,            & 
& MBM3,TYu,TYd,TYe,TYb3,TYw3,TYx3,Bmue,BMXM3,BMWM3,BMGM3,BMBM3,mq2,ml2,mHd2,             & 
& mHu2,md2,mu2,me2,mHw32,mHg32,mHb32,mHx32,mHxb32,MassB,MassWB,MassG,MnuL)

If (ThresholdCrossed.lt.1) Then 
MXM3(1,:) = 0._dp 
BMXM3(1,:) = 0._dp 
mHx32(1,:) = 0._dp 
mHx32(:,1) = 0._dp 
Yx3(1,:) = 0._dp 
TYx3(1,:) = 0._dp 
MXM3(:,1) = 0._dp 
BMXM3(:,1) = 0._dp 
mHxb32(1,:) = 0._dp 
mHxb32(:,1) = 0._dp 
MGM3(1,:) = 0._dp 
BMGM3(1,:) = 0._dp 
MGM3(:,1) = 0._dp 
BMGM3(:,1) = 0._dp 
mHg32(1,:) = 0._dp 
mHg32(:,1) = 0._dp 
Yb3(1,:) = 0._dp 
TYb3(1,:) = 0._dp 
MBM3(1,:) = 0._dp 
BMBM3(1,:) = 0._dp 
MBM3(:,1) = 0._dp 
BMBM3(:,1) = 0._dp 
mHb32(1,:) = 0._dp 
mHb32(:,1) = 0._dp 
Yw3(1,:) = 0._dp 
TYw3(1,:) = 0._dp 
MWM3(1,:) = 0._dp 
BMWM3(1,:) = 0._dp 
MWM3(:,1) = 0._dp 
BMWM3(:,1) = 0._dp 
mHw32(1,:) = 0._dp 
mHw32(:,1) = 0._dp 
End if 

If (ThresholdCrossed.lt.2) Then 
MXM3(2,:) = 0._dp 
BMXM3(2,:) = 0._dp 
mHx32(2,:) = 0._dp 
mHx32(:,2) = 0._dp 
Yx3(2,:) = 0._dp 
TYx3(2,:) = 0._dp 
MXM3(:,2) = 0._dp 
BMXM3(:,2) = 0._dp 
mHxb32(2,:) = 0._dp 
mHxb32(:,2) = 0._dp 
MGM3(2,:) = 0._dp 
BMGM3(2,:) = 0._dp 
MGM3(:,2) = 0._dp 
BMGM3(:,2) = 0._dp 
mHg32(2,:) = 0._dp 
mHg32(:,2) = 0._dp 
Yb3(2,:) = 0._dp 
TYb3(2,:) = 0._dp 
MBM3(2,:) = 0._dp 
BMBM3(2,:) = 0._dp 
MBM3(:,2) = 0._dp 
BMBM3(:,2) = 0._dp 
mHb32(2,:) = 0._dp 
mHb32(:,2) = 0._dp 
Yw3(2,:) = 0._dp 
TYw3(2,:) = 0._dp 
MWM3(2,:) = 0._dp 
BMWM3(2,:) = 0._dp 
MWM3(:,2) = 0._dp 
BMWM3(:,2) = 0._dp 
mHw32(2,:) = 0._dp 
mHw32(:,2) = 0._dp 
End if 

If (ThresholdCrossed.lt.3) Then 
MXM3(3,:) = 0._dp 
BMXM3(3,:) = 0._dp 
mHx32(3,:) = 0._dp 
mHx32(:,3) = 0._dp 
Yx3(3,:) = 0._dp 
TYx3(3,:) = 0._dp 
MXM3(:,3) = 0._dp 
BMXM3(:,3) = 0._dp 
mHxb32(3,:) = 0._dp 
mHxb32(:,3) = 0._dp 
MGM3(3,:) = 0._dp 
BMGM3(3,:) = 0._dp 
MGM3(:,3) = 0._dp 
BMGM3(:,3) = 0._dp 
mHg32(3,:) = 0._dp 
mHg32(:,3) = 0._dp 
Yb3(3,:) = 0._dp 
TYb3(3,:) = 0._dp 
MBM3(3,:) = 0._dp 
BMBM3(3,:) = 0._dp 
MBM3(:,3) = 0._dp 
BMBM3(:,3) = 0._dp 
mHb32(3,:) = 0._dp 
mHb32(:,3) = 0._dp 
Yw3(3,:) = 0._dp 
TYw3(3,:) = 0._dp 
MWM3(3,:) = 0._dp 
BMWM3(3,:) = 0._dp 
MWM3(:,3) = 0._dp 
BMWM3(:,3) = 0._dp 
mHw32(3,:) = 0._dp 
mHw32(:,3) = 0._dp 
End if 

AbsMassB = Conjg(MassB)*MassB
AbsMassWB = Conjg(MassWB)*MassWB
AbsMassG = Conjg(MassG)*MassG
Call Adjungate(Yu,adjYu)
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(Yb3,adjYb3)
Call Adjungate(Yw3,adjYw3)
Call Adjungate(Yx3,adjYx3)
Call Adjungate(MXM3,adjMXM3)
Call Adjungate(MWM3,adjMWM3)
Call Adjungate(MGM3,adjMGM3)
Call Adjungate(MBM3,adjMBM3)
Call Adjungate(TYu,adjTYu)
Call Adjungate(TYd,adjTYd)
Call Adjungate(TYe,adjTYe)
Call Adjungate(TYb3,adjTYb3)
Call Adjungate(TYw3,adjTYw3)
Call Adjungate(TYx3,adjTYx3)
!Call Adjungate(BMXM3,adjBMXM3)
!Call Adjungate(BMWM3,adjBMWM3)
!Call Adjungate(BMGM3,adjBMGM3)
!Call Adjungate(BMBM3,adjBMBM3)
Call Adjungate(mq2,adjmq2)
Call Adjungate(ml2,adjml2)
Call Adjungate(md2,adjmd2)
Call Adjungate(mu2,adjmu2)
Call Adjungate(me2,adjme2)
Call Adjungate(mHw32,adjmHw32)
Call Adjungate(mHg32,adjmHg32)
Call Adjungate(mHb32,adjmHb32)
Call Adjungate(mHx32,adjmHx32)
Call Adjungate(mHxb32,adjmHxb32)
 md2adjYx3 = Matmul(md2,adjYx3) 
 md2CYd = Matmul(md2,Conjg(Yd)) 
 me2CYe = Matmul(me2,Conjg(Ye)) 
 mHb32CYb3 = Matmul(mHb32,Conjg(Yb3)) 
 mHw32CYw3 = Matmul(mHw32,Conjg(Yw3)) 
 mHxb32CYx3 = Matmul(mHxb32,Conjg(Yx3)) 
 ml2adjYb3 = Matmul(ml2,adjYb3) 
 ml2adjYe = Matmul(ml2,adjYe) 
 ml2adjYw3 = Matmul(ml2,adjYw3) 
 mq2adjYd = Matmul(mq2,adjYd) 
 mq2adjYu = Matmul(mq2,adjYu) 
 mu2CYu = Matmul(mu2,Conjg(Yu)) 
 Yb3adjYb3 = Matmul(Yb3,adjYb3) 
Forall(i2=1:3)  Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3(i2,i2),dp) 
 YdadjYd = Matmul(Yd,adjYd) 
Forall(i2=1:3)  YdadjYd(i2,i2) =  Real(YdadjYd(i2,i2),dp) 
 YeadjYe = Matmul(Ye,adjYe) 
Forall(i2=1:3)  YeadjYe(i2,i2) =  Real(YeadjYe(i2,i2),dp) 
 YuadjYu = Matmul(Yu,adjYu) 
Forall(i2=1:3)  YuadjYu(i2,i2) =  Real(YuadjYu(i2,i2),dp) 
 Yw3adjYw3 = Matmul(Yw3,adjYw3) 
Forall(i2=1:3)  Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3(i2,i2),dp) 
 Yx3adjYx3 = Matmul(Yx3,adjYx3) 
Forall(i2=1:3)  Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3(i2,i2),dp) 
 adjYb3MBM3 = Matmul(adjYb3,MBM3) 
 adjYb3mHb32 = Matmul(adjYb3,mHb32) 
 adjYb3Yb3 = Matmul(adjYb3,Yb3) 
Forall(i2=1:3)  adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3(i2,i2),dp) 
! adjYb3BMBM3 = Matmul(adjYb3,BMBM3) 
 adjYb3TYb3 = Matmul(adjYb3,TYb3) 
 adjYdmd2 = Matmul(adjYd,md2) 
 adjYdYd = Matmul(adjYd,Yd) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYdTYd = Matmul(adjYd,TYd) 
 adjYeme2 = Matmul(adjYe,me2) 
 adjYeYe = Matmul(adjYe,Ye) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYeTYe = Matmul(adjYe,TYe) 
 adjYumu2 = Matmul(adjYu,mu2) 
 adjYuYu = Matmul(adjYu,Yu) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYuTYu = Matmul(adjYu,TYu) 
 adjYw3mHw32 = Matmul(adjYw3,mHw32) 
 adjYw3MWM3 = Matmul(adjYw3,MWM3) 
 adjYw3Yw3 = Matmul(adjYw3,Yw3) 
Forall(i2=1:3)  adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3(i2,i2),dp) 
! adjYw3BMWM3 = Matmul(adjYw3,BMWM3) 
 adjYw3TYw3 = Matmul(adjYw3,TYw3) 
 adjYx3mHxb32 = Matmul(adjYx3,mHxb32) 
 adjYx3Yx3 = Matmul(adjYx3,Yx3) 
Forall(i2=1:3)  adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3(i2,i2),dp) 
 adjYx3TYx3 = Matmul(adjYx3,TYx3) 
 CYb3ml2 = Matmul(Conjg(Yb3),ml2) 
 CYb3TpYb3 = Matmul(Conjg(Yb3),Transpose(Yb3)) 
Forall(i2=1:3)  CYb3TpYb3(i2,i2) =  Real(CYb3TpYb3(i2,i2),dp) 
 CYb3TpTYb3 = Matmul(Conjg(Yb3),Transpose(TYb3)) 
 CYdmq2 = Matmul(Conjg(Yd),mq2) 
 CYdTpYd = Matmul(Conjg(Yd),Transpose(Yd)) 
Forall(i2=1:3)  CYdTpYd(i2,i2) =  Real(CYdTpYd(i2,i2),dp) 
 CYdTpTYd = Matmul(Conjg(Yd),Transpose(TYd)) 
 CYeml2 = Matmul(Conjg(Ye),ml2) 
 CYumq2 = Matmul(Conjg(Yu),mq2) 
 CYw3ml2 = Matmul(Conjg(Yw3),ml2) 
 CYw3TpYw3 = Matmul(Conjg(Yw3),Transpose(Yw3)) 
Forall(i2=1:3)  CYw3TpYw3(i2,i2) =  Real(CYw3TpYw3(i2,i2),dp) 
 CYw3TpTYw3 = Matmul(Conjg(Yw3),Transpose(TYw3)) 
 CYx3md2 = Matmul(Conjg(Yx3),md2) 
 CYx3Yd = Matmul(Conjg(Yx3),Yd) 
 CYx3TYd = Matmul(Conjg(Yx3),TYd) 
 CYx3TpYx3 = Matmul(Conjg(Yx3),Transpose(Yx3)) 
Forall(i2=1:3)  CYx3TpYx3(i2,i2) =  Real(CYx3TpYx3(i2,i2),dp) 
 CYx3TpTYx3 = Matmul(Conjg(Yx3),Transpose(TYx3)) 
 CTYb3TpTYb3 = Matmul(Conjg(TYb3),Transpose(TYb3)) 
 CTYdTpTYd = Matmul(Conjg(TYd),Transpose(TYd)) 
 CTYeTpTYe = Matmul(Conjg(TYe),Transpose(TYe)) 
 CTYuTpTYu = Matmul(Conjg(TYu),Transpose(TYu)) 
 CTYw3TpTYw3 = Matmul(Conjg(TYw3),Transpose(TYw3)) 
 CTYx3TpTYx3 = Matmul(Conjg(TYx3),Transpose(TYx3)) 
 TYb3adjTYb3 = Matmul(TYb3,adjTYb3) 
 TYdadjTYd = Matmul(TYd,adjTYd) 
 TYeadjTYe = Matmul(TYe,adjTYe) 
 TYuadjTYu = Matmul(TYu,adjTYu) 
 TYw3adjTYw3 = Matmul(TYw3,adjTYw3) 
 TYx3adjTYx3 = Matmul(TYx3,adjTYx3) 
 TpYb3CYb3 = Matmul(Transpose(Yb3),Conjg(Yb3)) 
Forall(i2=1:3)  TpYb3CYb3(i2,i2) =  Real(TpYb3CYb3(i2,i2),dp) 
 TpYdCYd = Matmul(Transpose(Yd),Conjg(Yd)) 
Forall(i2=1:3)  TpYdCYd(i2,i2) =  Real(TpYdCYd(i2,i2),dp) 
 TpYeCYe = Matmul(Transpose(Ye),Conjg(Ye)) 
Forall(i2=1:3)  TpYeCYe(i2,i2) =  Real(TpYeCYe(i2,i2),dp) 
 TpYuCYu = Matmul(Transpose(Yu),Conjg(Yu)) 
Forall(i2=1:3)  TpYuCYu(i2,i2) =  Real(TpYuCYu(i2,i2),dp) 
 TpYw3CYw3 = Matmul(Transpose(Yw3),Conjg(Yw3)) 
Forall(i2=1:3)  TpYw3CYw3(i2,i2) =  Real(TpYw3CYw3(i2,i2),dp) 
 TpYx3CYx3 = Matmul(Transpose(Yx3),Conjg(Yx3)) 
Forall(i2=1:3)  TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3(i2,i2),dp) 
 TpTYb3CTYb3 = Matmul(Transpose(TYb3),Conjg(TYb3)) 
 TpTYdCTYd = Matmul(Transpose(TYd),Conjg(TYd)) 
 TpTYeCTYe = Matmul(Transpose(TYe),Conjg(TYe)) 
 TpTYuCTYu = Matmul(Transpose(TYu),Conjg(TYu)) 
 TpTYw3CTYw3 = Matmul(Transpose(TYw3),Conjg(TYw3)) 
 TpTYx3CTYx3 = Matmul(Transpose(TYx3),Conjg(TYx3)) 
 MBM3CYb3TpYb3 = Matmul(MBM3,CYb3TpYb3) 
 MBM3CYb3TpTYb3 = Matmul(MBM3,CYb3TpTYb3) 
 md2YdadjYd = Matmul(md2,YdadjYd) 
 md2adjYx3Yx3 = Matmul(md2,adjYx3Yx3) 
 md2TpYx3CYx3 = Matmul(md2,TpYx3CYx3) 
 me2YeadjYe = Matmul(me2,YeadjYe) 
 mHb32Yb3adjYb3 = Matmul(mHb32,Yb3adjYb3) 
 mHw32Yw3adjYw3 = Matmul(mHw32,Yw3adjYw3) 
 mHxb32Yx3adjYx3 = Matmul(mHxb32,Yx3adjYx3) 
 ml2adjYb3Yb3 = Matmul(ml2,adjYb3Yb3) 
 ml2adjYeYe = Matmul(ml2,adjYeYe) 
 ml2adjYw3Yw3 = Matmul(ml2,adjYw3Yw3) 
 ml2TpYb3CYb3 = Matmul(ml2,TpYb3CYb3) 
 ml2TpYeCYe = Matmul(ml2,TpYeCYe) 
 ml2TpYw3CYw3 = Matmul(ml2,TpYw3CYw3) 
 mq2adjYdYd = Matmul(mq2,adjYdYd) 
 mq2adjYuYu = Matmul(mq2,adjYuYu) 
 mq2TpYdCYd = Matmul(mq2,TpYdCYd) 
 mq2TpYuCYu = Matmul(mq2,TpYuCYu) 
 mu2YuadjYu = Matmul(mu2,YuadjYu) 
 MWM3CYw3TpYw3 = Matmul(MWM3,CYw3TpYw3) 
 MWM3CYw3TpTYw3 = Matmul(MWM3,CYw3TpTYw3) 
 MXM3CYx3TpYx3 = Matmul(MXM3,CYx3TpYx3) 
 MXM3CYx3TpTYx3 = Matmul(MXM3,CYx3TpTYx3) 
 Yb3ml2adjYb3 = Matmul(Yb3,ml2adjYb3) 
 Yb3adjYb3MBM3 = Matmul(Yb3,adjYb3MBM3) 
 Yb3adjYb3mHb32 = Matmul(Yb3,adjYb3mHb32) 
 Yb3adjYb3Yb3 = Matmul(Yb3,adjYb3Yb3) 
! Yb3adjYb3BMBM3 = Matmul(Yb3,adjYb3BMBM3) 
 Yb3adjYb3TYb3 = Matmul(Yb3,adjYb3TYb3) 
 Yb3adjYeYe = Matmul(Yb3,adjYeYe) 
 Yb3adjYeTYe = Matmul(Yb3,adjYeTYe) 
 Yb3adjYw3Yw3 = Matmul(Yb3,adjYw3Yw3) 
 Yb3adjYw3TYw3 = Matmul(Yb3,adjYw3TYw3) 
 Ydmq2adjYd = Matmul(Yd,mq2adjYd) 
 YdadjYdmd2 = Matmul(Yd,adjYdmd2) 
 YdadjYdYd = Matmul(Yd,adjYdYd) 
 YdadjYdTYd = Matmul(Yd,adjYdTYd) 
 YdadjYuYu = Matmul(Yd,adjYuYu) 
 YdadjYuTYu = Matmul(Yd,adjYuTYu) 
 Yeml2adjYe = Matmul(Ye,ml2adjYe) 
 YeadjYb3Yb3 = Matmul(Ye,adjYb3Yb3) 
 YeadjYb3TYb3 = Matmul(Ye,adjYb3TYb3) 
 YeadjYeme2 = Matmul(Ye,adjYeme2) 
 YeadjYeYe = Matmul(Ye,adjYeYe) 
 YeadjYeTYe = Matmul(Ye,adjYeTYe) 
 YeadjYw3Yw3 = Matmul(Ye,adjYw3Yw3) 
 YeadjYw3TYw3 = Matmul(Ye,adjYw3TYw3) 
 Yumq2adjYu = Matmul(Yu,mq2adjYu) 
 YuadjYdYd = Matmul(Yu,adjYdYd) 
 YuadjYdTYd = Matmul(Yu,adjYdTYd) 
 YuadjYumu2 = Matmul(Yu,adjYumu2) 
 YuadjYuYu = Matmul(Yu,adjYuYu) 
 YuadjYuTYu = Matmul(Yu,adjYuTYu) 
 Yw3ml2adjYw3 = Matmul(Yw3,ml2adjYw3) 
 Yw3adjYb3Yb3 = Matmul(Yw3,adjYb3Yb3) 
 Yw3adjYb3TYb3 = Matmul(Yw3,adjYb3TYb3) 
 Yw3adjYeYe = Matmul(Yw3,adjYeYe) 
 Yw3adjYeTYe = Matmul(Yw3,adjYeTYe) 
 Yw3adjYw3mHw32 = Matmul(Yw3,adjYw3mHw32) 
 Yw3adjYw3MWM3 = Matmul(Yw3,adjYw3MWM3) 
 Yw3adjYw3Yw3 = Matmul(Yw3,adjYw3Yw3) 
! Yw3adjYw3BMWM3 = Matmul(Yw3,adjYw3BMWM3) 
 Yw3adjYw3TYw3 = Matmul(Yw3,adjYw3TYw3) 
 Yx3md2adjYx3 = Matmul(Yx3,md2adjYx3) 
 Yx3adjYx3mHxb32 = Matmul(Yx3,adjYx3mHxb32) 
 Yx3adjYx3Yx3 = Matmul(Yx3,adjYx3Yx3) 
 Yx3adjYx3TYx3 = Matmul(Yx3,adjYx3TYx3) 
 Yx3CYdTpYd = Matmul(Yx3,CYdTpYd) 
 Yx3CYdTpTYd = Matmul(Yx3,CYdTpTYd) 
! BMBM3CYb3TpYb3 = Matmul(BMBM3,CYb3TpYb3) 
! BMWM3CYw3TpYw3 = Matmul(BMWM3,CYw3TpYw3) 
! BMXM3CYx3TpYx3 = Matmul(BMXM3,CYx3TpYx3) 
 TYb3adjYb3MBM3 = Matmul(TYb3,adjYb3MBM3) 
 TYb3adjYb3Yb3 = Matmul(TYb3,adjYb3Yb3) 
 TYb3adjYeYe = Matmul(TYb3,adjYeYe) 
 TYb3adjYw3Yw3 = Matmul(TYb3,adjYw3Yw3) 
 TYdadjYdYd = Matmul(TYd,adjYdYd) 
 TYdadjYuYu = Matmul(TYd,adjYuYu) 
 TYeadjYb3Yb3 = Matmul(TYe,adjYb3Yb3) 
 TYeadjYeYe = Matmul(TYe,adjYeYe) 
 TYeadjYw3Yw3 = Matmul(TYe,adjYw3Yw3) 
 TYuadjYdYd = Matmul(TYu,adjYdYd) 
 TYuadjYuYu = Matmul(TYu,adjYuYu) 
 TYw3adjYb3Yb3 = Matmul(TYw3,adjYb3Yb3) 
 TYw3adjYeYe = Matmul(TYw3,adjYeYe) 
 TYw3adjYw3MWM3 = Matmul(TYw3,adjYw3MWM3) 
 TYw3adjYw3Yw3 = Matmul(TYw3,adjYw3Yw3) 
 TYx3adjYx3Yx3 = Matmul(TYx3,adjYx3Yx3) 
 TYx3CYdTpYd = Matmul(TYx3,CYdTpYd) 
 TpYb3mHb32CYb3 = Matmul(Transpose(Yb3),mHb32CYb3) 
 TpYb3CYb3ml2 = Matmul(Transpose(Yb3),CYb3ml2) 
 TpYdmd2CYd = Matmul(Transpose(Yd),md2CYd) 
 TpYdCYdmq2 = Matmul(Transpose(Yd),CYdmq2) 
 TpYeme2CYe = Matmul(Transpose(Ye),me2CYe) 
 TpYeCYeml2 = Matmul(Transpose(Ye),CYeml2) 
 TpYumu2CYu = Matmul(Transpose(Yu),mu2CYu) 
 TpYuCYumq2 = Matmul(Transpose(Yu),CYumq2) 
 TpYw3mHw32CYw3 = Matmul(Transpose(Yw3),mHw32CYw3) 
 TpYw3CYw3ml2 = Matmul(Transpose(Yw3),CYw3ml2) 
 TpYx3mHxb32CYx3 = Matmul(Transpose(Yx3),mHxb32CYx3) 
 TpYx3CYx3md2 = Matmul(Transpose(Yx3),CYx3md2) 
 TpYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3Yd) 
 TpYx3CYx3TYd = Matmul(Transpose(Yx3),CYx3TYd) 
 TpTYx3CYx3Yd = Matmul(Transpose(TYx3),CYx3Yd) 
 Trmd2 = Real(cTrace(md2),dp) 
 Trme2 = Real(cTrace(me2),dp) 
 TrmHx32 = Real(cTrace(mHx32),dp) 
 TrmHxb32 = Real(cTrace(mHxb32),dp) 
 Trml2 = Real(cTrace(ml2),dp) 
 Trmq2 = Real(cTrace(mq2),dp) 
 Trmu2 = Real(cTrace(mu2),dp) 
 TrYb3adjYb3 = Real(cTrace(Yb3adjYb3),dp) 
 TrYdadjYd = Real(cTrace(YdadjYd),dp) 
 TrYeadjYe = Real(cTrace(YeadjYe),dp) 
 TrYuadjYu = Real(cTrace(YuadjYu),dp) 
 TrYw3adjYw3 = Real(cTrace(Yw3adjYw3),dp) 
 TrYx3adjYx3 = Real(cTrace(Yx3adjYx3),dp) 
 TradjYb3TYb3 = Real(cTrace(adjYb3TYb3),dp) 
 TradjYdTYd = Real(cTrace(adjYdTYd),dp) 
 TradjYeTYe = Real(cTrace(adjYeTYe),dp) 
 TradjYuTYu = Real(cTrace(adjYuTYu),dp) 
 TradjYw3TYw3 = Real(cTrace(adjYw3TYw3),dp) 
 TradjYx3TYx3 = Real(cTrace(adjYx3TYx3),dp) 
 TrCTYb3TpTYb3 = Real(cTrace(CTYb3TpTYb3),dp) 
 TrCTYdTpTYd = Real(cTrace(CTYdTpTYd),dp) 
 TrCTYeTpTYe = Real(cTrace(CTYeTpTYe),dp) 
 TrCTYuTpTYu = Real(cTrace(CTYuTpTYu),dp) 
 TrCTYw3TpTYw3 = Real(cTrace(CTYw3TpTYw3),dp) 
 TrCTYx3TpTYx3 = Real(cTrace(CTYx3TpTYx3),dp) 
 Trmd2YdadjYd = Real(cTrace(md2YdadjYd),dp) 
 Trmd2adjYx3Yx3 = Real(cTrace(md2adjYx3Yx3),dp) 
 Trme2YeadjYe = Real(cTrace(me2YeadjYe),dp) 
 TrmHb32Yb3adjYb3 = Real(cTrace(mHb32Yb3adjYb3),dp) 
 TrmHw32Yw3adjYw3 = Real(cTrace(mHw32Yw3adjYw3),dp) 
 TrmHxb32Yx3adjYx3 = Real(cTrace(mHxb32Yx3adjYx3),dp) 
 Trml2adjYb3Yb3 = Real(cTrace(ml2adjYb3Yb3),dp) 
 Trml2adjYeYe = Real(cTrace(ml2adjYeYe),dp) 
 Trml2adjYw3Yw3 = Real(cTrace(ml2adjYw3Yw3),dp) 
 Trmq2adjYdYd = Real(cTrace(mq2adjYdYd),dp) 
 Trmq2adjYuYu = Real(cTrace(mq2adjYuYu),dp) 
 Trmu2YuadjYu = Real(cTrace(mu2YuadjYu),dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 


If (TwoLoopRGE) Then 
 Yb3adjYe = Matmul(Yb3,adjYe) 
 Yb3adjYw3 = Matmul(Yb3,adjYw3) 
 Yb3adjTYb3 = Matmul(Yb3,adjTYb3) 
 Yb3adjTYe = Matmul(Yb3,adjTYe) 
 Yb3adjTYw3 = Matmul(Yb3,adjTYw3) 
 YdadjYu = Matmul(Yd,adjYu) 
 YdadjTYd = Matmul(Yd,adjTYd) 
 YdadjTYu = Matmul(Yd,adjTYu) 
 YeadjYb3 = Matmul(Ye,adjYb3) 
 YeadjYw3 = Matmul(Ye,adjYw3) 
 YeadjTYb3 = Matmul(Ye,adjTYb3) 
 YeadjTYe = Matmul(Ye,adjTYe) 
 YeadjTYw3 = Matmul(Ye,adjTYw3) 
 YuadjYd = Matmul(Yu,adjYd) 
 YuadjTYd = Matmul(Yu,adjTYd) 
 YuadjTYu = Matmul(Yu,adjTYu) 
 Yw3adjYb3 = Matmul(Yw3,adjYb3) 
 Yw3adjYe = Matmul(Yw3,adjYe) 
 Yw3adjTYb3 = Matmul(Yw3,adjTYb3) 
 Yw3adjTYe = Matmul(Yw3,adjTYe) 
 Yw3adjTYw3 = Matmul(Yw3,adjTYw3) 
 Yx3adjTYx3 = Matmul(Yx3,adjTYx3) 
 Yx3CYd = Matmul(Yx3,Conjg(Yd)) 
 Yx3CTYd = Matmul(Yx3,Conjg(TYd)) 
 adjYdadjTYx3 = Matmul(adjYd,adjTYx3) 
 adjYdTpYx3 = Matmul(adjYd,Transpose(Yx3)) 
 adjYdTpTYx3 = Matmul(adjYd,Transpose(TYx3)) 
 adjTYdadjYx3 = Matmul(adjTYd,adjYx3) 
 CYb3TpYw3 = Matmul(Conjg(Yb3),Transpose(Yw3)) 
 CYb3TpTYw3 = Matmul(Conjg(Yb3),Transpose(TYw3)) 
 CYeTpYb3 = Matmul(Conjg(Ye),Transpose(Yb3)) 
 CYeTpYw3 = Matmul(Conjg(Ye),Transpose(Yw3)) 
 CYeTpTYb3 = Matmul(Conjg(Ye),Transpose(TYb3)) 
 CYeTpTYw3 = Matmul(Conjg(Ye),Transpose(TYw3)) 
 CYuTpYd = Matmul(Conjg(Yu),Transpose(Yd)) 
 CYuTpTYd = Matmul(Conjg(Yu),Transpose(TYd)) 
 CYw3TpYb3 = Matmul(Conjg(Yw3),Transpose(Yb3)) 
 CYw3TpTYb3 = Matmul(Conjg(Yw3),Transpose(TYb3)) 
 CTYb3adjYb3 = Matmul(Conjg(TYb3),adjYb3) 
 CTYb3adjYe = Matmul(Conjg(TYb3),adjYe) 
 CTYb3adjYw3 = Matmul(Conjg(TYb3),adjYw3) 
 CTYb3TpYb3 = Matmul(Conjg(TYb3),Transpose(Yb3)) 
 CTYdadjYd = Matmul(Conjg(TYd),adjYd) 
 CTYdadjYu = Matmul(Conjg(TYd),adjYu) 
 CTYdTpYd = Matmul(Conjg(TYd),Transpose(Yd)) 
 CTYeadjYb3 = Matmul(Conjg(TYe),adjYb3) 
 CTYeadjYe = Matmul(Conjg(TYe),adjYe) 
 CTYeadjYw3 = Matmul(Conjg(TYe),adjYw3) 
 CTYeTpYe = Matmul(Conjg(TYe),Transpose(Ye)) 
 CTYuadjYd = Matmul(Conjg(TYu),adjYd) 
 CTYuadjYu = Matmul(Conjg(TYu),adjYu) 
 CTYuTpYu = Matmul(Conjg(TYu),Transpose(Yu)) 
 CTYw3adjYb3 = Matmul(Conjg(TYw3),adjYb3) 
 CTYw3adjYe = Matmul(Conjg(TYw3),adjYe) 
 CTYw3adjYw3 = Matmul(Conjg(TYw3),adjYw3) 
 CTYw3TpYw3 = Matmul(Conjg(TYw3),Transpose(Yw3)) 
 CTYx3adjYx3 = Matmul(Conjg(TYx3),adjYx3) 
 CTYx3TYd = Matmul(Conjg(TYx3),TYd) 
 CTYx3TpYd = Matmul(Conjg(TYx3),Transpose(Yd)) 
 CTYx3TpYx3 = Matmul(Conjg(TYx3),Transpose(Yx3)) 
 CTYx3TpTYd = Matmul(Conjg(TYx3),Transpose(TYd)) 
 TYb3adjYb3 = Matmul(TYb3,adjYb3) 
 TYb3adjYe = Matmul(TYb3,adjYe) 
 TYb3adjYw3 = Matmul(TYb3,adjYw3) 
 TYb3adjTYe = Matmul(TYb3,adjTYe) 
 TYb3adjTYw3 = Matmul(TYb3,adjTYw3) 
 TYdadjYd = Matmul(TYd,adjYd) 
 TYdadjYu = Matmul(TYd,adjYu) 
 TYdadjTYu = Matmul(TYd,adjTYu) 
 TYeadjYb3 = Matmul(TYe,adjYb3) 
 TYeadjYe = Matmul(TYe,adjYe) 
 TYeadjYw3 = Matmul(TYe,adjYw3) 
 TYeadjTYb3 = Matmul(TYe,adjTYb3) 
 TYeadjTYw3 = Matmul(TYe,adjTYw3) 
 TYuadjYd = Matmul(TYu,adjYd) 
 TYuadjYu = Matmul(TYu,adjYu) 
 TYuadjTYd = Matmul(TYu,adjTYd) 
 TYw3adjYb3 = Matmul(TYw3,adjYb3) 
 TYw3adjYe = Matmul(TYw3,adjYe) 
 TYw3adjYw3 = Matmul(TYw3,adjYw3) 
 TYw3adjTYb3 = Matmul(TYw3,adjTYb3) 
 TYw3adjTYe = Matmul(TYw3,adjTYe) 
 TYx3adjYx3 = Matmul(TYx3,adjYx3) 
 TYx3CTYd = Matmul(TYx3,Conjg(TYd)) 
 TpYb3CTYb3 = Matmul(Transpose(Yb3),Conjg(TYb3)) 
 TpYdadjYx3 = Matmul(Transpose(Yd),adjYx3) 
 TpYdadjTYx3 = Matmul(Transpose(Yd),adjTYx3) 
 TpYdCTYd = Matmul(Transpose(Yd),Conjg(TYd)) 
 TpYeCTYe = Matmul(Transpose(Ye),Conjg(TYe)) 
 TpYuCTYu = Matmul(Transpose(Yu),Conjg(TYu)) 
 TpYw3CTYw3 = Matmul(Transpose(Yw3),Conjg(TYw3)) 
 TpYx3CTYx3 = Matmul(Transpose(Yx3),Conjg(TYx3)) 
 TpTYb3CYb3 = Matmul(Transpose(TYb3),Conjg(Yb3)) 
 TpTYdadjYx3 = Matmul(Transpose(TYd),adjYx3) 
 TpTYdCYd = Matmul(Transpose(TYd),Conjg(Yd)) 
 TpTYeCYe = Matmul(Transpose(TYe),Conjg(Ye)) 
 TpTYuCYu = Matmul(Transpose(TYu),Conjg(Yu)) 
 TpTYw3CYw3 = Matmul(Transpose(TYw3),Conjg(Yw3)) 
 TpTYx3CYx3 = Matmul(Transpose(TYx3),Conjg(Yx3)) 
 md2YdadjYu = Matmul(md2,YdadjYu) 
 me2YeadjYb3 = Matmul(me2,YeadjYb3) 
 me2YeadjYw3 = Matmul(me2,YeadjYw3) 
 mHb32Yb3adjYe = Matmul(mHb32,Yb3adjYe) 
 mHb32Yb3adjYw3 = Matmul(mHb32,Yb3adjYw3) 
 mHw32Yw3adjYb3 = Matmul(mHw32,Yw3adjYb3) 
 mHw32Yw3adjYe = Matmul(mHw32,Yw3adjYe) 
 mHxb32Yx3CYd = Matmul(mHxb32,Yx3CYd) 
 mq2TpYdadjYx3 = Matmul(mq2,TpYdadjYx3) 
 mu2YuadjYd = Matmul(mu2,YuadjYd) 
 Yb3ml2adjYe = Matmul(Yb3,ml2adjYe) 
 Yb3ml2adjYw3 = Matmul(Yb3,ml2adjYw3) 
 Yb3adjYeme2 = Matmul(Yb3,adjYeme2) 
 Yb3adjYw3mHw32 = Matmul(Yb3,adjYw3mHw32) 
 Yb3adjYw3MWM3 = Matmul(Yb3,adjYw3MWM3) 
! Yb3adjYw3BMWM3 = Matmul(Yb3,adjYw3BMWM3) 
 Ydmq2adjYu = Matmul(Yd,mq2adjYu) 
 YdadjYdTpYx3 = Matmul(Yd,adjYdTpYx3) 
 YdadjYdTpTYx3 = Matmul(Yd,adjYdTpTYx3) 
 YdadjYumu2 = Matmul(Yd,adjYumu2) 
 YdadjTYdadjYx3 = Matmul(Yd,adjTYdadjYx3) 
 Yeml2adjYb3 = Matmul(Ye,ml2adjYb3) 
 Yeml2adjYw3 = Matmul(Ye,ml2adjYw3) 
 YeadjYb3MBM3 = Matmul(Ye,adjYb3MBM3) 
 YeadjYb3mHb32 = Matmul(Ye,adjYb3mHb32) 
! YeadjYb3BMBM3 = Matmul(Ye,adjYb3BMBM3) 
 YeadjYw3mHw32 = Matmul(Ye,adjYw3mHw32) 
 YeadjYw3MWM3 = Matmul(Ye,adjYw3MWM3) 
! YeadjYw3BMWM3 = Matmul(Ye,adjYw3BMWM3) 
 Yumq2adjYd = Matmul(Yu,mq2adjYd) 
 YuadjYdmd2 = Matmul(Yu,adjYdmd2) 
 Yw3ml2adjYb3 = Matmul(Yw3,ml2adjYb3) 
 Yw3ml2adjYe = Matmul(Yw3,ml2adjYe) 
 Yw3adjYb3MBM3 = Matmul(Yw3,adjYb3MBM3) 
 Yw3adjYb3mHb32 = Matmul(Yw3,adjYb3mHb32) 
! Yw3adjYb3BMBM3 = Matmul(Yw3,adjYb3BMBM3) 
 Yw3adjYeme2 = Matmul(Yw3,adjYeme2) 
 Yx3md2CYd = Matmul(Yx3,md2CYd) 
 Yx3CYdmq2 = Matmul(Yx3,CYdmq2) 
 adjYb3Yb3adjYb3 = Matmul(adjYb3,Yb3adjYb3) 
 adjYb3Yb3adjYe = Matmul(adjYb3,Yb3adjYe) 
 adjYb3Yb3adjYw3 = Matmul(adjYb3,Yb3adjYw3) 
 adjYb3Yb3adjTYb3 = Matmul(adjYb3,Yb3adjTYb3) 
 adjYb3Yb3adjTYe = Matmul(adjYb3,Yb3adjTYe) 
 adjYb3Yb3adjTYw3 = Matmul(adjYb3,Yb3adjTYw3) 
 adjYb3TYb3adjYb3 = Matmul(adjYb3,TYb3adjYb3) 
 adjYb3TYb3adjYe = Matmul(adjYb3,TYb3adjYe) 
 adjYb3TYb3adjYw3 = Matmul(adjYb3,TYb3adjYw3) 
 adjYb3TYb3adjTYb3 = Matmul(adjYb3,TYb3adjTYb3) 
 adjYb3TYb3adjTYe = Matmul(adjYb3,TYb3adjTYe) 
 adjYb3TYb3adjTYw3 = Matmul(adjYb3,TYb3adjTYw3) 
 adjYdYdadjYd = Matmul(adjYd,YdadjYd) 
 adjYdYdadjYu = Matmul(adjYd,YdadjYu) 
 adjYdYdadjTYd = Matmul(adjYd,YdadjTYd) 
 adjYdYdadjTYu = Matmul(adjYd,YdadjTYu) 
 adjYdadjYx3Yx3 = Matmul(adjYd,adjYx3Yx3) 
 adjYdTYdadjYd = Matmul(adjYd,TYdadjYd) 
 adjYdTYdadjYu = Matmul(adjYd,TYdadjYu) 
 adjYdTYdadjTYd = Matmul(adjYd,TYdadjTYd) 
 adjYdTYdadjTYu = Matmul(adjYd,TYdadjTYu) 
 adjYdTpYx3CYx3 = Matmul(adjYd,TpYx3CYx3) 
 adjYdTpYx3CTYx3 = Matmul(adjYd,TpYx3CTYx3) 
 adjYdTpTYx3CTYx3 = Matmul(adjYd,TpTYx3CTYx3) 
 adjYeYeadjYb3 = Matmul(adjYe,YeadjYb3) 
 adjYeYeadjYe = Matmul(adjYe,YeadjYe) 
 adjYeYeadjYw3 = Matmul(adjYe,YeadjYw3) 
 adjYeYeadjTYb3 = Matmul(adjYe,YeadjTYb3) 
 adjYeYeadjTYe = Matmul(adjYe,YeadjTYe) 
 adjYeYeadjTYw3 = Matmul(adjYe,YeadjTYw3) 
 adjYeTYeadjYb3 = Matmul(adjYe,TYeadjYb3) 
 adjYeTYeadjYe = Matmul(adjYe,TYeadjYe) 
 adjYeTYeadjYw3 = Matmul(adjYe,TYeadjYw3) 
 adjYeTYeadjTYb3 = Matmul(adjYe,TYeadjTYb3) 
 adjYeTYeadjTYe = Matmul(adjYe,TYeadjTYe) 
 adjYeTYeadjTYw3 = Matmul(adjYe,TYeadjTYw3) 
 adjYuYuadjYd = Matmul(adjYu,YuadjYd) 
 adjYuYuadjYu = Matmul(adjYu,YuadjYu) 
 adjYuYuadjTYd = Matmul(adjYu,YuadjTYd) 
 adjYuYuadjTYu = Matmul(adjYu,YuadjTYu) 
 adjYuTYuadjYd = Matmul(adjYu,TYuadjYd) 
 adjYuTYuadjYu = Matmul(adjYu,TYuadjYu) 
 adjYuTYuadjTYd = Matmul(adjYu,TYuadjTYd) 
 adjYuTYuadjTYu = Matmul(adjYu,TYuadjTYu) 
 adjYw3Yw3adjYb3 = Matmul(adjYw3,Yw3adjYb3) 
 adjYw3Yw3adjYe = Matmul(adjYw3,Yw3adjYe) 
 adjYw3Yw3adjYw3 = Matmul(adjYw3,Yw3adjYw3) 
 adjYw3Yw3adjTYb3 = Matmul(adjYw3,Yw3adjTYb3) 
 adjYw3Yw3adjTYe = Matmul(adjYw3,Yw3adjTYe) 
 adjYw3Yw3adjTYw3 = Matmul(adjYw3,Yw3adjTYw3) 
 adjYw3TYw3adjYb3 = Matmul(adjYw3,TYw3adjYb3) 
 adjYw3TYw3adjYe = Matmul(adjYw3,TYw3adjYe) 
 adjYw3TYw3adjYw3 = Matmul(adjYw3,TYw3adjYw3) 
 adjYw3TYw3adjTYb3 = Matmul(adjYw3,TYw3adjTYb3) 
 adjYw3TYw3adjTYe = Matmul(adjYw3,TYw3adjTYe) 
 adjYw3TYw3adjTYw3 = Matmul(adjYw3,TYw3adjTYw3) 
 adjYx3Yx3adjYx3 = Matmul(adjYx3,Yx3adjYx3) 
 adjYx3Yx3adjTYx3 = Matmul(adjYx3,Yx3adjTYx3) 
 adjYx3Yx3CYd = Matmul(adjYx3,Yx3CYd) 
 adjYx3Yx3CTYd = Matmul(adjYx3,Yx3CTYd) 
 adjYx3TYx3adjYx3 = Matmul(adjYx3,TYx3adjYx3) 
 adjYx3TYx3adjTYx3 = Matmul(adjYx3,TYx3adjTYx3) 
 adjYx3TYx3CTYd = Matmul(adjYx3,TYx3CTYd) 
 adjTYb3TYb3adjYb3 = Matmul(adjTYb3,TYb3adjYb3) 
 adjTYb3TYb3adjYe = Matmul(adjTYb3,TYb3adjYe) 
 adjTYb3TYb3adjYw3 = Matmul(adjTYb3,TYb3adjYw3) 
 adjTYdTYdadjYd = Matmul(adjTYd,TYdadjYd) 
 adjTYdTYdadjYu = Matmul(adjTYd,TYdadjYu) 
 adjTYeTYeadjYb3 = Matmul(adjTYe,TYeadjYb3) 
 adjTYeTYeadjYe = Matmul(adjTYe,TYeadjYe) 
 adjTYeTYeadjYw3 = Matmul(adjTYe,TYeadjYw3) 
 adjTYuTYuadjYd = Matmul(adjTYu,TYuadjYd) 
 adjTYuTYuadjYu = Matmul(adjTYu,TYuadjYu) 
 adjTYw3TYw3adjYb3 = Matmul(adjTYw3,TYw3adjYb3) 
 adjTYw3TYw3adjYe = Matmul(adjTYw3,TYw3adjYe) 
 adjTYw3TYw3adjYw3 = Matmul(adjTYw3,TYw3adjYw3) 
 adjTYx3TYx3adjYx3 = Matmul(adjTYx3,TYx3adjYx3) 
 CYb3TpYb3CYb3 = Matmul(Conjg(Yb3),TpYb3CYb3) 
 CYb3TpYb3CTYb3 = Matmul(Conjg(Yb3),TpYb3CTYb3) 
 CYb3TpYw3CYw3 = Matmul(Conjg(Yb3),TpYw3CYw3) 
 CYb3TpYw3CTYw3 = Matmul(Conjg(Yb3),TpYw3CTYw3) 
 CYb3TpTYb3CTYb3 = Matmul(Conjg(Yb3),TpTYb3CTYb3) 
 CYb3TpTYw3CTYw3 = Matmul(Conjg(Yb3),TpTYw3CTYw3) 
 CYdTpYdadjYx3 = Matmul(Conjg(Yd),TpYdadjYx3) 
 CYdTpYdadjTYx3 = Matmul(Conjg(Yd),TpYdadjTYx3) 
 CYdTpYdCYd = Matmul(Conjg(Yd),TpYdCYd) 
 CYdTpYdCTYd = Matmul(Conjg(Yd),TpYdCTYd) 
 CYdTpTYdCTYd = Matmul(Conjg(Yd),TpTYdCTYd) 
 CYeTpYeCYe = Matmul(Conjg(Ye),TpYeCYe) 
 CYeTpYeCTYe = Matmul(Conjg(Ye),TpYeCTYe) 
 CYeTpTYeCTYe = Matmul(Conjg(Ye),TpTYeCTYe) 
 CYuTpYuCYu = Matmul(Conjg(Yu),TpYuCYu) 
 CYuTpYuCTYu = Matmul(Conjg(Yu),TpYuCTYu) 
 CYuTpTYuCTYu = Matmul(Conjg(Yu),TpTYuCTYu) 
 CYw3TpYb3CYb3 = Matmul(Conjg(Yw3),TpYb3CYb3) 
 CYw3TpYb3CTYb3 = Matmul(Conjg(Yw3),TpYb3CTYb3) 
 CYw3TpYw3CYw3 = Matmul(Conjg(Yw3),TpYw3CYw3) 
 CYw3TpYw3CTYw3 = Matmul(Conjg(Yw3),TpYw3CTYw3) 
 CYw3TpTYb3CTYb3 = Matmul(Conjg(Yw3),TpTYb3CTYb3) 
 CYw3TpTYw3CTYw3 = Matmul(Conjg(Yw3),TpTYw3CTYw3) 
 CYx3YdadjYd = Matmul(Conjg(Yx3),YdadjYd) 
 CYx3TpYx3CYx3 = Matmul(Conjg(Yx3),TpYx3CYx3) 
 CYx3TpYx3CTYx3 = Matmul(Conjg(Yx3),TpYx3CTYx3) 
 CYx3TpTYx3CTYx3 = Matmul(Conjg(Yx3),TpTYx3CTYx3) 
 CTYb3TpTYb3CYb3 = Matmul(Conjg(TYb3),TpTYb3CYb3) 
 CTYb3TpTYw3CYw3 = Matmul(Conjg(TYb3),TpTYw3CYw3) 
 CTYdTpTYdadjYx3 = Matmul(Conjg(TYd),TpTYdadjYx3) 
 CTYdTpTYdCYd = Matmul(Conjg(TYd),TpTYdCYd) 
 CTYeTpTYeCYe = Matmul(Conjg(TYe),TpTYeCYe) 
 CTYuTpTYuCYu = Matmul(Conjg(TYu),TpTYuCYu) 
 CTYw3TpTYb3CYb3 = Matmul(Conjg(TYw3),TpTYb3CYb3) 
 CTYw3TpTYw3CYw3 = Matmul(Conjg(TYw3),TpTYw3CYw3) 
 CTYx3TpTYx3CYx3 = Matmul(Conjg(TYx3),TpTYx3CYx3) 
 TYb3adjYw3MWM3 = Matmul(TYb3,adjYw3MWM3) 
 TYb3TpYb3CYb3 = Matmul(TYb3,TpYb3CYb3) 
 TYb3TpYw3CYw3 = Matmul(TYb3,TpYw3CYw3) 
 TYdadjYdadjTYx3 = Matmul(TYd,adjYdadjTYx3) 
 TYdadjYdTpYx3 = Matmul(TYd,adjYdTpYx3) 
 TYdTpYdCYd = Matmul(TYd,TpYdCYd) 
 TYeadjYb3MBM3 = Matmul(TYe,adjYb3MBM3) 
 TYeadjYw3MWM3 = Matmul(TYe,adjYw3MWM3) 
 TYeTpYeCYe = Matmul(TYe,TpYeCYe) 
 TYuTpYuCYu = Matmul(TYu,TpYuCYu) 
 TYw3adjYb3MBM3 = Matmul(TYw3,adjYb3MBM3) 
 TYw3TpYb3CYb3 = Matmul(TYw3,TpYb3CYb3) 
 TYw3TpYw3CYw3 = Matmul(TYw3,TpYw3CYw3) 
 TYx3CTYdTpYd = Matmul(TYx3,CTYdTpYd) 
 TYx3TpYx3CYx3 = Matmul(TYx3,TpYx3CYx3) 
 TpYb3CYb3TpYb3 = Matmul(Transpose(Yb3),CYb3TpYb3) 
 TpYb3CYb3TpYw3 = Matmul(Transpose(Yb3),CYb3TpYw3) 
 TpYb3CYb3TpTYb3 = Matmul(Transpose(Yb3),CYb3TpTYb3) 
 TpYb3CYb3TpTYw3 = Matmul(Transpose(Yb3),CYb3TpTYw3) 
 TpYb3CTYb3adjYb3 = Matmul(Transpose(Yb3),CTYb3adjYb3) 
 TpYb3CTYb3adjYe = Matmul(Transpose(Yb3),CTYb3adjYe) 
 TpYb3CTYb3adjYw3 = Matmul(Transpose(Yb3),CTYb3adjYw3) 
 TpYdmd2adjYx3 = Matmul(Transpose(Yd),md2adjYx3) 
 TpYdadjYx3mHxb32 = Matmul(Transpose(Yd),adjYx3mHxb32) 
 TpYdadjYx3Yx3 = Matmul(Transpose(Yd),adjYx3Yx3) 
 TpYdadjYx3TYx3 = Matmul(Transpose(Yd),adjYx3TYx3) 
 TpYdCYdTpYd = Matmul(Transpose(Yd),CYdTpYd) 
 TpYdCYdTpTYd = Matmul(Transpose(Yd),CYdTpTYd) 
 TpYdCTYdadjYd = Matmul(Transpose(Yd),CTYdadjYd) 
 TpYdCTYdadjYu = Matmul(Transpose(Yd),CTYdadjYu) 
 TpYeCYeTpYb3 = Matmul(Transpose(Ye),CYeTpYb3) 
 TpYeCYeTpYw3 = Matmul(Transpose(Ye),CYeTpYw3) 
 TpYeCYeTpTYb3 = Matmul(Transpose(Ye),CYeTpTYb3) 
 TpYeCYeTpTYw3 = Matmul(Transpose(Ye),CYeTpTYw3) 
 TpYeCTYeadjYb3 = Matmul(Transpose(Ye),CTYeadjYb3) 
 TpYeCTYeadjYe = Matmul(Transpose(Ye),CTYeadjYe) 
 TpYeCTYeadjYw3 = Matmul(Transpose(Ye),CTYeadjYw3) 
 TpYuCYuTpYd = Matmul(Transpose(Yu),CYuTpYd) 
 TpYuCYuTpTYd = Matmul(Transpose(Yu),CYuTpTYd) 
 TpYuCTYuadjYd = Matmul(Transpose(Yu),CTYuadjYd) 
 TpYuCTYuadjYu = Matmul(Transpose(Yu),CTYuadjYu) 
 TpYw3CYw3TpYb3 = Matmul(Transpose(Yw3),CYw3TpYb3) 
 TpYw3CYw3TpYw3 = Matmul(Transpose(Yw3),CYw3TpYw3) 
 TpYw3CYw3TpTYb3 = Matmul(Transpose(Yw3),CYw3TpTYb3) 
 TpYw3CYw3TpTYw3 = Matmul(Transpose(Yw3),CYw3TpTYw3) 
 TpYw3CTYw3adjYb3 = Matmul(Transpose(Yw3),CTYw3adjYb3) 
 TpYw3CTYw3adjYe = Matmul(Transpose(Yw3),CTYw3adjYe) 
 TpYw3CTYw3adjYw3 = Matmul(Transpose(Yw3),CTYw3adjYw3) 
 TpYx3CYx3TpYx3 = Matmul(Transpose(Yx3),CYx3TpYx3) 
 TpYx3CYx3TpTYx3 = Matmul(Transpose(Yx3),CYx3TpTYx3) 
 TpYx3CTYx3adjYx3 = Matmul(Transpose(Yx3),CTYx3adjYx3) 
 TpYx3CTYx3TYd = Matmul(Transpose(Yx3),CTYx3TYd) 
 TpYx3CTYx3TpTYd = Matmul(Transpose(Yx3),CTYx3TpTYd) 
 TpTYb3CYb3TpYb3 = Matmul(Transpose(TYb3),CYb3TpYb3) 
 TpTYb3CYb3TpYw3 = Matmul(Transpose(TYb3),CYb3TpYw3) 
 TpTYb3CTYb3adjYb3 = Matmul(Transpose(TYb3),CTYb3adjYb3) 
 TpTYb3CTYb3adjYe = Matmul(Transpose(TYb3),CTYb3adjYe) 
 TpTYb3CTYb3adjYw3 = Matmul(Transpose(TYb3),CTYb3adjYw3) 
 TpTYdCYdTpYd = Matmul(Transpose(TYd),CYdTpYd) 
 TpTYdCTYdadjYd = Matmul(Transpose(TYd),CTYdadjYd) 
 TpTYdCTYdadjYu = Matmul(Transpose(TYd),CTYdadjYu) 
 TpTYeCYeTpYb3 = Matmul(Transpose(TYe),CYeTpYb3) 
 TpTYeCYeTpYw3 = Matmul(Transpose(TYe),CYeTpYw3) 
 TpTYeCTYeadjYb3 = Matmul(Transpose(TYe),CTYeadjYb3) 
 TpTYeCTYeadjYe = Matmul(Transpose(TYe),CTYeadjYe) 
 TpTYeCTYeadjYw3 = Matmul(Transpose(TYe),CTYeadjYw3) 
 TpTYuCYuTpYd = Matmul(Transpose(TYu),CYuTpYd) 
 TpTYuCTYuadjYd = Matmul(Transpose(TYu),CTYuadjYd) 
 TpTYuCTYuadjYu = Matmul(Transpose(TYu),CTYuadjYu) 
 TpTYw3CYw3TpYb3 = Matmul(Transpose(TYw3),CYw3TpYb3) 
 TpTYw3CYw3TpYw3 = Matmul(Transpose(TYw3),CYw3TpYw3) 
 TpTYw3CTYw3adjYb3 = Matmul(Transpose(TYw3),CTYw3adjYb3) 
 TpTYw3CTYw3adjYe = Matmul(Transpose(TYw3),CTYw3adjYe) 
 TpTYw3CTYw3adjYw3 = Matmul(Transpose(TYw3),CTYw3adjYw3) 
 TpTYx3CYx3TpYx3 = Matmul(Transpose(TYx3),CYx3TpYx3) 
 TpTYx3CTYx3adjYx3 = Matmul(Transpose(TYx3),CTYx3adjYx3) 
 TpTYx3CTYx3TpYd = Matmul(Transpose(TYx3),CTYx3TpYd) 
 md2adjYx3Yx3adjYx3 = Matmul(md2,adjYx3Yx3adjYx3) 
 md2adjYx3Yx3CYd = Matmul(md2,adjYx3Yx3CYd) 
 md2CYdTpYdadjYx3 = Matmul(md2,CYdTpYdadjYx3) 
 md2CYdTpYdCYd = Matmul(md2,CYdTpYdCYd) 
 me2CYeTpYeCYe = Matmul(me2,CYeTpYeCYe) 
 mHb32CYb3TpYb3CYb3 = Matmul(mHb32,CYb3TpYb3CYb3) 
 mHb32CYb3TpYw3CYw3 = Matmul(mHb32,CYb3TpYw3CYw3) 
 mHw32CYw3TpYb3CYb3 = Matmul(mHw32,CYw3TpYb3CYb3) 
 mHw32CYw3TpYw3CYw3 = Matmul(mHw32,CYw3TpYw3CYw3) 
 mHxb32CYx3TpYx3CYx3 = Matmul(mHxb32,CYx3TpYx3CYx3) 
 ml2adjYb3Yb3adjYb3 = Matmul(ml2,adjYb3Yb3adjYb3) 
 ml2adjYb3Yb3adjYe = Matmul(ml2,adjYb3Yb3adjYe) 
 ml2adjYb3Yb3adjYw3 = Matmul(ml2,adjYb3Yb3adjYw3) 
 ml2adjYeYeadjYb3 = Matmul(ml2,adjYeYeadjYb3) 
 ml2adjYeYeadjYe = Matmul(ml2,adjYeYeadjYe) 
 ml2adjYeYeadjYw3 = Matmul(ml2,adjYeYeadjYw3) 
 ml2adjYw3Yw3adjYb3 = Matmul(ml2,adjYw3Yw3adjYb3) 
 ml2adjYw3Yw3adjYe = Matmul(ml2,adjYw3Yw3adjYe) 
 ml2adjYw3Yw3adjYw3 = Matmul(ml2,adjYw3Yw3adjYw3) 
 mq2adjYdYdadjYd = Matmul(mq2,adjYdYdadjYd) 
 mq2adjYdYdadjYu = Matmul(mq2,adjYdYdadjYu) 
 mq2adjYuYuadjYd = Matmul(mq2,adjYuYuadjYd) 
 mq2adjYuYuadjYu = Matmul(mq2,adjYuYuadjYu) 
 mq2adjYx3Yx3CYd = Matmul(mq2,adjYx3Yx3CYd) 
 mu2CYuTpYuCYu = Matmul(mu2,CYuTpYuCYu) 
 Yb3adjYb3Yb3adjYb3 = Matmul(Yb3,adjYb3Yb3adjYb3) 
Forall(i2=1:3)  Yb3adjYb3Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3Yb3adjYb3(i2,i2),dp) 
 Yb3adjYb3TYb3adjYb3 = Matmul(Yb3,adjYb3TYb3adjYb3) 
 Yb3adjYb3TYb3adjTYb3 = Matmul(Yb3,adjYb3TYb3adjTYb3) 
 Yb3adjYeYeadjYb3 = Matmul(Yb3,adjYeYeadjYb3) 
Forall(i2=1:3)  Yb3adjYeYeadjYb3(i2,i2) =  Real(Yb3adjYeYeadjYb3(i2,i2),dp) 
 Yb3adjYeTYeadjYb3 = Matmul(Yb3,adjYeTYeadjYb3) 
 Yb3adjYeTYeadjTYb3 = Matmul(Yb3,adjYeTYeadjTYb3) 
 Yb3adjYw3Yw3adjYb3 = Matmul(Yb3,adjYw3Yw3adjYb3) 
Forall(i2=1:3)  Yb3adjYw3Yw3adjYb3(i2,i2) =  Real(Yb3adjYw3Yw3adjYb3(i2,i2),dp) 
 Yb3adjYw3TYw3adjYb3 = Matmul(Yb3,adjYw3TYw3adjYb3) 
 Yb3adjYw3TYw3adjTYb3 = Matmul(Yb3,adjYw3TYw3adjTYb3) 
 Yb3adjTYb3TYb3adjYb3 = Matmul(Yb3,adjTYb3TYb3adjYb3) 
 Yb3adjTYeTYeadjYb3 = Matmul(Yb3,adjTYeTYeadjYb3) 
 Yb3adjTYw3TYw3adjYb3 = Matmul(Yb3,adjTYw3TYw3adjYb3) 
 Yb3TpTYb3CTYb3adjYb3 = Matmul(Yb3,TpTYb3CTYb3adjYb3) 
 Yb3TpTYeCTYeadjYb3 = Matmul(Yb3,TpTYeCTYeadjYb3) 
 Yb3TpTYw3CTYw3adjYb3 = Matmul(Yb3,TpTYw3CTYw3adjYb3) 
 YdadjYdYdadjYd = Matmul(Yd,adjYdYdadjYd) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdTYdadjYd = Matmul(Yd,adjYdTYdadjYd) 
 YdadjYdTYdadjTYd = Matmul(Yd,adjYdTYdadjTYd) 
 YdadjYdTpYx3CYx3 = Matmul(Yd,adjYdTpYx3CYx3) 
 YdadjYdTpTYx3CTYx3 = Matmul(Yd,adjYdTpTYx3CTYx3) 
 YdadjYuYuadjYd = Matmul(Yd,adjYuYuadjYd) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YdadjYuTYuadjYd = Matmul(Yd,adjYuTYuadjYd) 
 YdadjYuTYuadjTYd = Matmul(Yd,adjYuTYuadjTYd) 
 YdadjTYdTYdadjYd = Matmul(Yd,adjTYdTYdadjYd) 
 YdadjTYuTYuadjYd = Matmul(Yd,adjTYuTYuadjYd) 
 YdTpTYdCTYdadjYd = Matmul(Yd,TpTYdCTYdadjYd) 
 YdTpTYuCTYuadjYd = Matmul(Yd,TpTYuCTYuadjYd) 
 YeadjYb3Yb3adjYe = Matmul(Ye,adjYb3Yb3adjYe) 
Forall(i2=1:3)  YeadjYb3Yb3adjYe(i2,i2) =  Real(YeadjYb3Yb3adjYe(i2,i2),dp) 
 YeadjYb3TYb3adjYe = Matmul(Ye,adjYb3TYb3adjYe) 
 YeadjYb3TYb3adjTYe = Matmul(Ye,adjYb3TYb3adjTYe) 
 YeadjYeYeadjYe = Matmul(Ye,adjYeYeadjYe) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYeTYeadjYe = Matmul(Ye,adjYeTYeadjYe) 
 YeadjYeTYeadjTYe = Matmul(Ye,adjYeTYeadjTYe) 
 YeadjYw3Yw3adjYe = Matmul(Ye,adjYw3Yw3adjYe) 
Forall(i2=1:3)  YeadjYw3Yw3adjYe(i2,i2) =  Real(YeadjYw3Yw3adjYe(i2,i2),dp) 
 YeadjYw3TYw3adjYe = Matmul(Ye,adjYw3TYw3adjYe) 
 YeadjYw3TYw3adjTYe = Matmul(Ye,adjYw3TYw3adjTYe) 
 YeadjTYb3TYb3adjYe = Matmul(Ye,adjTYb3TYb3adjYe) 
 YeadjTYeTYeadjYe = Matmul(Ye,adjTYeTYeadjYe) 
 YeadjTYw3TYw3adjYe = Matmul(Ye,adjTYw3TYw3adjYe) 
 YeTpTYb3CTYb3adjYe = Matmul(Ye,TpTYb3CTYb3adjYe) 
 YeTpTYeCTYeadjYe = Matmul(Ye,TpTYeCTYeadjYe) 
 YeTpTYw3CTYw3adjYe = Matmul(Ye,TpTYw3CTYw3adjYe) 
 YuadjYdYdadjYu = Matmul(Yu,adjYdYdadjYu) 
Forall(i2=1:3)  YuadjYdYdadjYu(i2,i2) =  Real(YuadjYdYdadjYu(i2,i2),dp) 
 YuadjYdTYdadjYu = Matmul(Yu,adjYdTYdadjYu) 
 YuadjYdTYdadjTYu = Matmul(Yu,adjYdTYdadjTYu) 
 YuadjYuYuadjYu = Matmul(Yu,adjYuYuadjYu) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 YuadjYuTYuadjYu = Matmul(Yu,adjYuTYuadjYu) 
 YuadjYuTYuadjTYu = Matmul(Yu,adjYuTYuadjTYu) 
 YuadjTYdTYdadjYu = Matmul(Yu,adjTYdTYdadjYu) 
 YuadjTYuTYuadjYu = Matmul(Yu,adjTYuTYuadjYu) 
 YuTpTYdCTYdadjYu = Matmul(Yu,TpTYdCTYdadjYu) 
 YuTpTYuCTYuadjYu = Matmul(Yu,TpTYuCTYuadjYu) 
 Yw3adjYb3Yb3adjYw3 = Matmul(Yw3,adjYb3Yb3adjYw3) 
Forall(i2=1:3)  Yw3adjYb3Yb3adjYw3(i2,i2) =  Real(Yw3adjYb3Yb3adjYw3(i2,i2),dp) 
 Yw3adjYb3TYb3adjYw3 = Matmul(Yw3,adjYb3TYb3adjYw3) 
 Yw3adjYb3TYb3adjTYw3 = Matmul(Yw3,adjYb3TYb3adjTYw3) 
 Yw3adjYeYeadjYw3 = Matmul(Yw3,adjYeYeadjYw3) 
Forall(i2=1:3)  Yw3adjYeYeadjYw3(i2,i2) =  Real(Yw3adjYeYeadjYw3(i2,i2),dp) 
 Yw3adjYeTYeadjYw3 = Matmul(Yw3,adjYeTYeadjYw3) 
 Yw3adjYeTYeadjTYw3 = Matmul(Yw3,adjYeTYeadjTYw3) 
 Yw3adjYw3Yw3adjYw3 = Matmul(Yw3,adjYw3Yw3adjYw3) 
Forall(i2=1:3)  Yw3adjYw3Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3Yw3adjYw3(i2,i2),dp) 
 Yw3adjYw3TYw3adjYw3 = Matmul(Yw3,adjYw3TYw3adjYw3) 
 Yw3adjYw3TYw3adjTYw3 = Matmul(Yw3,adjYw3TYw3adjTYw3) 
 Yw3adjTYb3TYb3adjYw3 = Matmul(Yw3,adjTYb3TYb3adjYw3) 
 Yw3adjTYeTYeadjYw3 = Matmul(Yw3,adjTYeTYeadjYw3) 
 Yw3adjTYw3TYw3adjYw3 = Matmul(Yw3,adjTYw3TYw3adjYw3) 
 Yw3TpTYb3CTYb3adjYw3 = Matmul(Yw3,TpTYb3CTYb3adjYw3) 
 Yw3TpTYeCTYeadjYw3 = Matmul(Yw3,TpTYeCTYeadjYw3) 
 Yw3TpTYw3CTYw3adjYw3 = Matmul(Yw3,TpTYw3CTYw3adjYw3) 
 Yx3adjYx3Yx3adjYx3 = Matmul(Yx3,adjYx3Yx3adjYx3) 
Forall(i2=1:3)  Yx3adjYx3Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3Yx3adjYx3(i2,i2),dp) 
 Yx3adjYx3TYx3adjYx3 = Matmul(Yx3,adjYx3TYx3adjYx3) 
 Yx3adjYx3TYx3adjTYx3 = Matmul(Yx3,adjYx3TYx3adjTYx3) 
 Yx3adjTYx3TYx3adjYx3 = Matmul(Yx3,adjTYx3TYx3adjYx3) 
 Yx3CYdTpYdadjYx3 = Matmul(Yx3,CYdTpYdadjYx3) 
Forall(i2=1:3)  Yx3CYdTpYdadjYx3(i2,i2) =  Real(Yx3CYdTpYdadjYx3(i2,i2),dp) 
 Yx3CTYdTpTYdadjYx3 = Matmul(Yx3,CTYdTpTYdadjYx3) 
 Yx3TYdadjYdadjTYx3 = Matmul(Yx3,TYdadjYdadjTYx3) 
 Yx3TpTYx3CTYx3adjYx3 = Matmul(Yx3,TpTYx3CTYx3adjYx3) 
 adjYb3mHb32Yb3adjYb3 = Matmul(adjYb3,mHb32Yb3adjYb3) 
 adjYb3mHb32Yb3adjYe = Matmul(adjYb3,mHb32Yb3adjYe) 
 adjYb3mHb32Yb3adjYw3 = Matmul(adjYb3,mHb32Yb3adjYw3) 
 adjYb3Yb3ml2adjYb3 = Matmul(adjYb3,Yb3ml2adjYb3) 
 adjYb3Yb3ml2adjYe = Matmul(adjYb3,Yb3ml2adjYe) 
 adjYb3Yb3ml2adjYw3 = Matmul(adjYb3,Yb3ml2adjYw3) 
 adjYb3Yb3adjYb3MBM3 = Matmul(adjYb3,Yb3adjYb3MBM3) 
 adjYb3Yb3adjYb3mHb32 = Matmul(adjYb3,Yb3adjYb3mHb32) 
 adjYb3Yb3adjYb3Yb3 = Matmul(adjYb3,Yb3adjYb3Yb3) 
Forall(i2=1:3)  adjYb3Yb3adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3adjYb3Yb3(i2,i2),dp) 
! adjYb3Yb3adjYb3BMBM3 = Matmul(adjYb3,Yb3adjYb3BMBM3) 
 adjYb3Yb3adjYb3TYb3 = Matmul(adjYb3,Yb3adjYb3TYb3) 
 adjYb3Yb3adjYeme2 = Matmul(adjYb3,Yb3adjYeme2) 
 adjYb3Yb3adjYeYe = Matmul(adjYb3,Yb3adjYeYe) 
 adjYb3Yb3adjYeTYe = Matmul(adjYb3,Yb3adjYeTYe) 
 adjYb3Yb3adjYw3mHw32 = Matmul(adjYb3,Yb3adjYw3mHw32) 
 adjYb3Yb3adjYw3MWM3 = Matmul(adjYb3,Yb3adjYw3MWM3) 
 adjYb3Yb3adjYw3Yw3 = Matmul(adjYb3,Yb3adjYw3Yw3) 
! adjYb3Yb3adjYw3BMWM3 = Matmul(adjYb3,Yb3adjYw3BMWM3) 
 adjYb3Yb3adjYw3TYw3 = Matmul(adjYb3,Yb3adjYw3TYw3) 
 adjYb3TYb3adjYb3MBM3 = Matmul(adjYb3,TYb3adjYb3MBM3) 
 adjYb3TYb3adjYb3Yb3 = Matmul(adjYb3,TYb3adjYb3Yb3) 
 adjYb3TYb3adjYeYe = Matmul(adjYb3,TYb3adjYeYe) 
 adjYb3TYb3adjYw3MWM3 = Matmul(adjYb3,TYb3adjYw3MWM3) 
 adjYb3TYb3adjYw3Yw3 = Matmul(adjYb3,TYb3adjYw3Yw3) 
 adjYdmd2YdadjYd = Matmul(adjYd,md2YdadjYd) 
 adjYdmd2YdadjYu = Matmul(adjYd,md2YdadjYu) 
 adjYdYdmq2adjYd = Matmul(adjYd,Ydmq2adjYd) 
 adjYdYdmq2adjYu = Matmul(adjYd,Ydmq2adjYu) 
 adjYdYdadjYdmd2 = Matmul(adjYd,YdadjYdmd2) 
 adjYdYdadjYdYd = Matmul(adjYd,YdadjYdYd) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYdTYd = Matmul(adjYd,YdadjYdTYd) 
 adjYdYdadjYumu2 = Matmul(adjYd,YdadjYumu2) 
 adjYdYdadjYuYu = Matmul(adjYd,YdadjYuYu) 
 adjYdYdadjYuTYu = Matmul(adjYd,YdadjYuTYu) 
 adjYdTYdadjYdYd = Matmul(adjYd,TYdadjYdYd) 
 adjYdTYdadjYuYu = Matmul(adjYd,TYdadjYuYu) 
 adjYdTpYx3CYx3Yd = Matmul(adjYd,TpYx3CYx3Yd) 
Forall(i2=1:3)  adjYdTpYx3CYx3Yd(i2,i2) =  Real(adjYdTpYx3CYx3Yd(i2,i2),dp) 
 adjYdTpYx3CYx3TYd = Matmul(adjYd,TpYx3CYx3TYd) 
 adjYdTpYx3CTYx3TYd = Matmul(adjYd,TpYx3CTYx3TYd) 
 adjYdTpYx3CTYx3TpTYd = Matmul(adjYd,TpYx3CTYx3TpTYd) 
 adjYdTpTYx3CYx3Yd = Matmul(adjYd,TpTYx3CYx3Yd) 
 adjYdTpTYx3CTYx3TpYd = Matmul(adjYd,TpTYx3CTYx3TpYd) 
 adjYeme2YeadjYb3 = Matmul(adjYe,me2YeadjYb3) 
 adjYeme2YeadjYe = Matmul(adjYe,me2YeadjYe) 
 adjYeme2YeadjYw3 = Matmul(adjYe,me2YeadjYw3) 
 adjYeYeml2adjYb3 = Matmul(adjYe,Yeml2adjYb3) 
 adjYeYeml2adjYe = Matmul(adjYe,Yeml2adjYe) 
 adjYeYeml2adjYw3 = Matmul(adjYe,Yeml2adjYw3) 
 adjYeYeadjYb3MBM3 = Matmul(adjYe,YeadjYb3MBM3) 
 adjYeYeadjYb3mHb32 = Matmul(adjYe,YeadjYb3mHb32) 
 adjYeYeadjYb3Yb3 = Matmul(adjYe,YeadjYb3Yb3) 
! adjYeYeadjYb3BMBM3 = Matmul(adjYe,YeadjYb3BMBM3) 
 adjYeYeadjYb3TYb3 = Matmul(adjYe,YeadjYb3TYb3) 
 adjYeYeadjYeme2 = Matmul(adjYe,YeadjYeme2) 
 adjYeYeadjYeYe = Matmul(adjYe,YeadjYeYe) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYeTYe = Matmul(adjYe,YeadjYeTYe) 
 adjYeYeadjYw3mHw32 = Matmul(adjYe,YeadjYw3mHw32) 
 adjYeYeadjYw3MWM3 = Matmul(adjYe,YeadjYw3MWM3) 
 adjYeYeadjYw3Yw3 = Matmul(adjYe,YeadjYw3Yw3) 
! adjYeYeadjYw3BMWM3 = Matmul(adjYe,YeadjYw3BMWM3) 
 adjYeYeadjYw3TYw3 = Matmul(adjYe,YeadjYw3TYw3) 
 adjYeTYeadjYb3MBM3 = Matmul(adjYe,TYeadjYb3MBM3) 
 adjYeTYeadjYb3Yb3 = Matmul(adjYe,TYeadjYb3Yb3) 
 adjYeTYeadjYeYe = Matmul(adjYe,TYeadjYeYe) 
 adjYeTYeadjYw3MWM3 = Matmul(adjYe,TYeadjYw3MWM3) 
 adjYeTYeadjYw3Yw3 = Matmul(adjYe,TYeadjYw3Yw3) 
 adjYumu2YuadjYd = Matmul(adjYu,mu2YuadjYd) 
 adjYumu2YuadjYu = Matmul(adjYu,mu2YuadjYu) 
 adjYuYumq2adjYd = Matmul(adjYu,Yumq2adjYd) 
 adjYuYumq2adjYu = Matmul(adjYu,Yumq2adjYu) 
 adjYuYuadjYdmd2 = Matmul(adjYu,YuadjYdmd2) 
 adjYuYuadjYdYd = Matmul(adjYu,YuadjYdYd) 
 adjYuYuadjYdTYd = Matmul(adjYu,YuadjYdTYd) 
 adjYuYuadjYumu2 = Matmul(adjYu,YuadjYumu2) 
 adjYuYuadjYuYu = Matmul(adjYu,YuadjYuYu) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYuYuadjYuTYu = Matmul(adjYu,YuadjYuTYu) 
 adjYuTYuadjYdYd = Matmul(adjYu,TYuadjYdYd) 
 adjYuTYuadjYuYu = Matmul(adjYu,TYuadjYuYu) 
 adjYw3mHw32Yw3adjYb3 = Matmul(adjYw3,mHw32Yw3adjYb3) 
 adjYw3mHw32Yw3adjYe = Matmul(adjYw3,mHw32Yw3adjYe) 
 adjYw3mHw32Yw3adjYw3 = Matmul(adjYw3,mHw32Yw3adjYw3) 
 adjYw3Yw3ml2adjYb3 = Matmul(adjYw3,Yw3ml2adjYb3) 
 adjYw3Yw3ml2adjYe = Matmul(adjYw3,Yw3ml2adjYe) 
 adjYw3Yw3ml2adjYw3 = Matmul(adjYw3,Yw3ml2adjYw3) 
 adjYw3Yw3adjYb3MBM3 = Matmul(adjYw3,Yw3adjYb3MBM3) 
 adjYw3Yw3adjYb3mHb32 = Matmul(adjYw3,Yw3adjYb3mHb32) 
 adjYw3Yw3adjYb3Yb3 = Matmul(adjYw3,Yw3adjYb3Yb3) 
! adjYw3Yw3adjYb3BMBM3 = Matmul(adjYw3,Yw3adjYb3BMBM3) 
 adjYw3Yw3adjYb3TYb3 = Matmul(adjYw3,Yw3adjYb3TYb3) 
 adjYw3Yw3adjYeme2 = Matmul(adjYw3,Yw3adjYeme2) 
 adjYw3Yw3adjYeYe = Matmul(adjYw3,Yw3adjYeYe) 
 adjYw3Yw3adjYeTYe = Matmul(adjYw3,Yw3adjYeTYe) 
 adjYw3Yw3adjYw3mHw32 = Matmul(adjYw3,Yw3adjYw3mHw32) 
 adjYw3Yw3adjYw3MWM3 = Matmul(adjYw3,Yw3adjYw3MWM3) 
 adjYw3Yw3adjYw3Yw3 = Matmul(adjYw3,Yw3adjYw3Yw3) 
Forall(i2=1:3)  adjYw3Yw3adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3adjYw3Yw3(i2,i2),dp) 
! adjYw3Yw3adjYw3BMWM3 = Matmul(adjYw3,Yw3adjYw3BMWM3) 
 adjYw3Yw3adjYw3TYw3 = Matmul(adjYw3,Yw3adjYw3TYw3) 
 adjYw3TYw3adjYb3MBM3 = Matmul(adjYw3,TYw3adjYb3MBM3) 
 adjYw3TYw3adjYb3Yb3 = Matmul(adjYw3,TYw3adjYb3Yb3) 
 adjYw3TYw3adjYeYe = Matmul(adjYw3,TYw3adjYeYe) 
 adjYw3TYw3adjYw3MWM3 = Matmul(adjYw3,TYw3adjYw3MWM3) 
 adjYw3TYw3adjYw3Yw3 = Matmul(adjYw3,TYw3adjYw3Yw3) 
 adjYx3mHxb32Yx3adjYx3 = Matmul(adjYx3,mHxb32Yx3adjYx3) 
 adjYx3mHxb32Yx3CYd = Matmul(adjYx3,mHxb32Yx3CYd) 
 adjYx3Yx3md2adjYx3 = Matmul(adjYx3,Yx3md2adjYx3) 
 adjYx3Yx3md2CYd = Matmul(adjYx3,Yx3md2CYd) 
 adjYx3Yx3adjYx3mHxb32 = Matmul(adjYx3,Yx3adjYx3mHxb32) 
 adjYx3Yx3adjYx3Yx3 = Matmul(adjYx3,Yx3adjYx3Yx3) 
Forall(i2=1:3)  adjYx3Yx3adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3adjYx3Yx3(i2,i2),dp) 
 adjYx3Yx3adjYx3TYx3 = Matmul(adjYx3,Yx3adjYx3TYx3) 
 adjYx3Yx3CYdmq2 = Matmul(adjYx3,Yx3CYdmq2) 
 adjYx3TYx3adjYx3Yx3 = Matmul(adjYx3,TYx3adjYx3Yx3) 
 adjYx3TYx3CYdTpYd = Matmul(adjYx3,TYx3CYdTpYd) 
 adjYx3TYx3CTYdTpYd = Matmul(adjYx3,TYx3CTYdTpYd) 
 adjTYb3TYb3TpYb3CYb3 = Matmul(adjTYb3,TYb3TpYb3CYb3) 
 adjTYb3TYb3TpYw3CYw3 = Matmul(adjTYb3,TYb3TpYw3CYw3) 
 adjTYdTYdTpYdCYd = Matmul(adjTYd,TYdTpYdCYd) 
 adjTYeTYeTpYeCYe = Matmul(adjTYe,TYeTpYeCYe) 
 adjTYuTYuTpYuCYu = Matmul(adjTYu,TYuTpYuCYu) 
 adjTYw3TYw3TpYb3CYb3 = Matmul(adjTYw3,TYw3TpYb3CYb3) 
 adjTYw3TYw3TpYw3CYw3 = Matmul(adjTYw3,TYw3TpYw3CYw3) 
 adjTYx3TYx3TpYx3CYx3 = Matmul(adjTYx3,TYx3TpYx3CYx3) 
 CYb3ml2TpYb3CYb3 = Matmul(Conjg(Yb3),ml2TpYb3CYb3) 
 CYb3ml2TpYw3CYw3 = Matmul(Conjg(Yb3),ml2TpYw3CYw3) 
 CYb3TpYb3mHb32CYb3 = Matmul(Conjg(Yb3),TpYb3mHb32CYb3) 
 CYb3TpYb3CYb3ml2 = Matmul(Conjg(Yb3),TpYb3CYb3ml2) 
 CYb3TpYb3CYb3TpYb3 = Matmul(Conjg(Yb3),TpYb3CYb3TpYb3) 
Forall(i2=1:3)  CYb3TpYb3CYb3TpYb3(i2,i2) =  Real(CYb3TpYb3CYb3TpYb3(i2,i2),dp) 
 CYb3TpYb3CYb3TpTYb3 = Matmul(Conjg(Yb3),TpYb3CYb3TpTYb3) 
 CYb3TpYeCYeTpYb3 = Matmul(Conjg(Yb3),TpYeCYeTpYb3) 
Forall(i2=1:3)  CYb3TpYeCYeTpYb3(i2,i2) =  Real(CYb3TpYeCYeTpYb3(i2,i2),dp) 
 CYb3TpYeCYeTpTYb3 = Matmul(Conjg(Yb3),TpYeCYeTpTYb3) 
 CYb3TpYw3mHw32CYw3 = Matmul(Conjg(Yb3),TpYw3mHw32CYw3) 
 CYb3TpYw3CYw3ml2 = Matmul(Conjg(Yb3),TpYw3CYw3ml2) 
 CYb3TpYw3CYw3TpYb3 = Matmul(Conjg(Yb3),TpYw3CYw3TpYb3) 
Forall(i2=1:3)  CYb3TpYw3CYw3TpYb3(i2,i2) =  Real(CYb3TpYw3CYw3TpYb3(i2,i2),dp) 
 CYb3TpYw3CYw3TpTYb3 = Matmul(Conjg(Yb3),TpYw3CYw3TpTYb3) 
 CYb3TpTYb3CYb3TpYb3 = Matmul(Conjg(Yb3),TpTYb3CYb3TpYb3) 
 CYb3TpTYeCYeTpYb3 = Matmul(Conjg(Yb3),TpTYeCYeTpYb3) 
 CYb3TpTYw3CYw3TpYb3 = Matmul(Conjg(Yb3),TpTYw3CYw3TpYb3) 
 CYdmq2TpYdadjYx3 = Matmul(Conjg(Yd),mq2TpYdadjYx3) 
 CYdmq2TpYdCYd = Matmul(Conjg(Yd),mq2TpYdCYd) 
 CYdTpYdmd2adjYx3 = Matmul(Conjg(Yd),TpYdmd2adjYx3) 
 CYdTpYdmd2CYd = Matmul(Conjg(Yd),TpYdmd2CYd) 
 CYdTpYdadjYx3mHxb32 = Matmul(Conjg(Yd),TpYdadjYx3mHxb32) 
 CYdTpYdadjYx3Yx3 = Matmul(Conjg(Yd),TpYdadjYx3Yx3) 
 CYdTpYdadjYx3TYx3 = Matmul(Conjg(Yd),TpYdadjYx3TYx3) 
 CYdTpYdCYdmq2 = Matmul(Conjg(Yd),TpYdCYdmq2) 
 CYdTpYdCYdTpYd = Matmul(Conjg(Yd),TpYdCYdTpYd) 
Forall(i2=1:3)  CYdTpYdCYdTpYd(i2,i2) =  Real(CYdTpYdCYdTpYd(i2,i2),dp) 
 CYdTpYdCYdTpTYd = Matmul(Conjg(Yd),TpYdCYdTpTYd) 
 CYdTpYuCYuTpYd = Matmul(Conjg(Yd),TpYuCYuTpYd) 
Forall(i2=1:3)  CYdTpYuCYuTpYd(i2,i2) =  Real(CYdTpYuCYuTpYd(i2,i2),dp) 
 CYdTpYuCYuTpTYd = Matmul(Conjg(Yd),TpYuCYuTpTYd) 
 CYdTpTYdCYdTpYd = Matmul(Conjg(Yd),TpTYdCYdTpYd) 
 CYdTpTYuCYuTpYd = Matmul(Conjg(Yd),TpTYuCYuTpYd) 
 CYeml2TpYeCYe = Matmul(Conjg(Ye),ml2TpYeCYe) 
 CYeTpYeme2CYe = Matmul(Conjg(Ye),TpYeme2CYe) 
 CYeTpYeCYeml2 = Matmul(Conjg(Ye),TpYeCYeml2) 
 CYumq2TpYuCYu = Matmul(Conjg(Yu),mq2TpYuCYu) 
 CYuTpYumu2CYu = Matmul(Conjg(Yu),TpYumu2CYu) 
 CYuTpYuCYumq2 = Matmul(Conjg(Yu),TpYuCYumq2) 
 CYw3ml2TpYb3CYb3 = Matmul(Conjg(Yw3),ml2TpYb3CYb3) 
 CYw3ml2TpYw3CYw3 = Matmul(Conjg(Yw3),ml2TpYw3CYw3) 
 CYw3TpYb3mHb32CYb3 = Matmul(Conjg(Yw3),TpYb3mHb32CYb3) 
 CYw3TpYb3CYb3ml2 = Matmul(Conjg(Yw3),TpYb3CYb3ml2) 
 CYw3TpYb3CYb3TpYw3 = Matmul(Conjg(Yw3),TpYb3CYb3TpYw3) 
Forall(i2=1:3)  CYw3TpYb3CYb3TpYw3(i2,i2) =  Real(CYw3TpYb3CYb3TpYw3(i2,i2),dp) 
 CYw3TpYb3CYb3TpTYw3 = Matmul(Conjg(Yw3),TpYb3CYb3TpTYw3) 
 CYw3TpYeCYeTpYw3 = Matmul(Conjg(Yw3),TpYeCYeTpYw3) 
Forall(i2=1:3)  CYw3TpYeCYeTpYw3(i2,i2) =  Real(CYw3TpYeCYeTpYw3(i2,i2),dp) 
 CYw3TpYeCYeTpTYw3 = Matmul(Conjg(Yw3),TpYeCYeTpTYw3) 
 CYw3TpYw3mHw32CYw3 = Matmul(Conjg(Yw3),TpYw3mHw32CYw3) 
 CYw3TpYw3CYw3ml2 = Matmul(Conjg(Yw3),TpYw3CYw3ml2) 
 CYw3TpYw3CYw3TpYw3 = Matmul(Conjg(Yw3),TpYw3CYw3TpYw3) 
Forall(i2=1:3)  CYw3TpYw3CYw3TpYw3(i2,i2) =  Real(CYw3TpYw3CYw3TpYw3(i2,i2),dp) 
 CYw3TpYw3CYw3TpTYw3 = Matmul(Conjg(Yw3),TpYw3CYw3TpTYw3) 
 CYw3TpTYb3CYb3TpYw3 = Matmul(Conjg(Yw3),TpTYb3CYb3TpYw3) 
 CYw3TpTYeCYeTpYw3 = Matmul(Conjg(Yw3),TpTYeCYeTpYw3) 
 CYw3TpTYw3CYw3TpYw3 = Matmul(Conjg(Yw3),TpTYw3CYw3TpYw3) 
 CYx3md2TpYx3CYx3 = Matmul(Conjg(Yx3),md2TpYx3CYx3) 
 CYx3YdadjYdTpYx3 = Matmul(Conjg(Yx3),YdadjYdTpYx3) 
Forall(i2=1:3)  CYx3YdadjYdTpYx3(i2,i2) =  Real(CYx3YdadjYdTpYx3(i2,i2),dp) 
 CYx3YdadjYdTpTYx3 = Matmul(Conjg(Yx3),YdadjYdTpTYx3) 
 CYx3TYdadjYdTpYx3 = Matmul(Conjg(Yx3),TYdadjYdTpYx3) 
 CYx3TpYx3mHxb32CYx3 = Matmul(Conjg(Yx3),TpYx3mHxb32CYx3) 
 CYx3TpYx3CYx3md2 = Matmul(Conjg(Yx3),TpYx3CYx3md2) 
 CYx3TpYx3CYx3Yd = Matmul(Conjg(Yx3),TpYx3CYx3Yd) 
 CYx3TpYx3CYx3TYd = Matmul(Conjg(Yx3),TpYx3CYx3TYd) 
 CYx3TpYx3CYx3TpYx3 = Matmul(Conjg(Yx3),TpYx3CYx3TpYx3) 
Forall(i2=1:3)  CYx3TpYx3CYx3TpYx3(i2,i2) =  Real(CYx3TpYx3CYx3TpYx3(i2,i2),dp) 
 CYx3TpYx3CYx3TpTYx3 = Matmul(Conjg(Yx3),TpYx3CYx3TpTYx3) 
 CYx3TpTYx3CYx3Yd = Matmul(Conjg(Yx3),TpTYx3CYx3Yd) 
 CYx3TpTYx3CYx3TpYx3 = Matmul(Conjg(Yx3),TpTYx3CYx3TpYx3) 
 TYb3adjYb3Yb3adjTYb3 = Matmul(TYb3,adjYb3Yb3adjTYb3) 
 TYb3adjYeYeadjTYb3 = Matmul(TYb3,adjYeYeadjTYb3) 
 TYb3adjYw3Yw3adjTYb3 = Matmul(TYb3,adjYw3Yw3adjTYb3) 
 TYb3TpYb3CTYb3adjYb3 = Matmul(TYb3,TpYb3CTYb3adjYb3) 
 TYb3TpYeCTYeadjYb3 = Matmul(TYb3,TpYeCTYeadjYb3) 
 TYb3TpYw3CTYw3adjYb3 = Matmul(TYb3,TpYw3CTYw3adjYb3) 
 TYdadjYdYdadjTYd = Matmul(TYd,adjYdYdadjTYd) 
 TYdadjYdadjYx3Yx3 = Matmul(TYd,adjYdadjYx3Yx3) 
 TYdadjYuYuadjTYd = Matmul(TYd,adjYuYuadjTYd) 
 TYdTpYdCTYdadjYd = Matmul(TYd,TpYdCTYdadjYd) 
 TYdTpYuCTYuadjYd = Matmul(TYd,TpYuCTYuadjYd) 
 TYeadjYb3Yb3adjTYe = Matmul(TYe,adjYb3Yb3adjTYe) 
 TYeadjYeYeadjTYe = Matmul(TYe,adjYeYeadjTYe) 
 TYeadjYw3Yw3adjTYe = Matmul(TYe,adjYw3Yw3adjTYe) 
 TYeTpYb3CTYb3adjYe = Matmul(TYe,TpYb3CTYb3adjYe) 
 TYeTpYeCTYeadjYe = Matmul(TYe,TpYeCTYeadjYe) 
 TYeTpYw3CTYw3adjYe = Matmul(TYe,TpYw3CTYw3adjYe) 
 TYuadjYdYdadjTYu = Matmul(TYu,adjYdYdadjTYu) 
 TYuadjYuYuadjTYu = Matmul(TYu,adjYuYuadjTYu) 
 TYuTpYdCTYdadjYu = Matmul(TYu,TpYdCTYdadjYu) 
 TYuTpYuCTYuadjYu = Matmul(TYu,TpYuCTYuadjYu) 
 TYw3adjYb3Yb3adjTYw3 = Matmul(TYw3,adjYb3Yb3adjTYw3) 
 TYw3adjYeYeadjTYw3 = Matmul(TYw3,adjYeYeadjTYw3) 
 TYw3adjYw3Yw3adjTYw3 = Matmul(TYw3,adjYw3Yw3adjTYw3) 
 TYw3TpYb3CTYb3adjYw3 = Matmul(TYw3,TpYb3CTYb3adjYw3) 
 TYw3TpYeCTYeadjYw3 = Matmul(TYw3,TpYeCTYeadjYw3) 
 TYw3TpYw3CTYw3adjYw3 = Matmul(TYw3,TpYw3CTYw3adjYw3) 
 TYx3YdadjTYdadjYx3 = Matmul(TYx3,YdadjTYdadjYx3) 
 TYx3adjYx3Yx3adjTYx3 = Matmul(TYx3,adjYx3Yx3adjTYx3) 
 TYx3CYdTpYdadjTYx3 = Matmul(TYx3,CYdTpYdadjTYx3) 
 TYx3TpYx3CTYx3adjYx3 = Matmul(TYx3,TpYx3CTYx3adjYx3) 
 TpYb3CYb3TpYb3CYb3 = Matmul(Transpose(Yb3),CYb3TpYb3CYb3) 
Forall(i2=1:3)  TpYb3CYb3TpYb3CYb3(i2,i2) =  Real(TpYb3CYb3TpYb3CYb3(i2,i2),dp) 
 TpYb3CYb3TpYw3CYw3 = Matmul(Transpose(Yb3),CYb3TpYw3CYw3) 
 TpYb3CYb3TpTYb3CTYb3 = Matmul(Transpose(Yb3),CYb3TpTYb3CTYb3) 
 TpYb3CYb3TpTYw3CTYw3 = Matmul(Transpose(Yb3),CYb3TpTYw3CTYw3) 
 TpYb3CTYb3TpTYb3CYb3 = Matmul(Transpose(Yb3),CTYb3TpTYb3CYb3) 
 TpYb3CTYb3TpTYw3CYw3 = Matmul(Transpose(Yb3),CTYb3TpTYw3CYw3) 
 TpYdadjYdTpTYx3CTYx3 = Matmul(Transpose(Yd),adjYdTpTYx3CTYx3) 
 TpYdadjYx3Yx3CYd = Matmul(Transpose(Yd),adjYx3Yx3CYd) 
Forall(i2=1:3)  TpYdadjYx3Yx3CYd(i2,i2) =  Real(TpYdadjYx3Yx3CYd(i2,i2),dp) 
 TpYdadjYx3TYx3CTYd = Matmul(Transpose(Yd),adjYx3TYx3CTYd) 
 TpYdCYdTpYdCYd = Matmul(Transpose(Yd),CYdTpYdCYd) 
Forall(i2=1:3)  TpYdCYdTpYdCYd(i2,i2) =  Real(TpYdCYdTpYdCYd(i2,i2),dp) 
 TpYdCYdTpTYdCTYd = Matmul(Transpose(Yd),CYdTpTYdCTYd) 
 TpYdCTYdTpTYdCYd = Matmul(Transpose(Yd),CTYdTpTYdCYd) 
 TpYeCYeTpYeCYe = Matmul(Transpose(Ye),CYeTpYeCYe) 
Forall(i2=1:3)  TpYeCYeTpYeCYe(i2,i2) =  Real(TpYeCYeTpYeCYe(i2,i2),dp) 
 TpYeCYeTpTYeCTYe = Matmul(Transpose(Ye),CYeTpTYeCTYe) 
 TpYeCTYeTpTYeCYe = Matmul(Transpose(Ye),CTYeTpTYeCYe) 
 TpYuCYuTpYuCYu = Matmul(Transpose(Yu),CYuTpYuCYu) 
Forall(i2=1:3)  TpYuCYuTpYuCYu(i2,i2) =  Real(TpYuCYuTpYuCYu(i2,i2),dp) 
 TpYuCYuTpTYuCTYu = Matmul(Transpose(Yu),CYuTpTYuCTYu) 
 TpYuCTYuTpTYuCYu = Matmul(Transpose(Yu),CTYuTpTYuCYu) 
 TpYw3CYw3TpYb3CYb3 = Matmul(Transpose(Yw3),CYw3TpYb3CYb3) 
 TpYw3CYw3TpYw3CYw3 = Matmul(Transpose(Yw3),CYw3TpYw3CYw3) 
Forall(i2=1:3)  TpYw3CYw3TpYw3CYw3(i2,i2) =  Real(TpYw3CYw3TpYw3CYw3(i2,i2),dp) 
 TpYw3CYw3TpTYb3CTYb3 = Matmul(Transpose(Yw3),CYw3TpTYb3CTYb3) 
 TpYw3CYw3TpTYw3CTYw3 = Matmul(Transpose(Yw3),CYw3TpTYw3CTYw3) 
 TpYw3CTYw3TpTYb3CYb3 = Matmul(Transpose(Yw3),CTYw3TpTYb3CYb3) 
 TpYw3CTYw3TpTYw3CYw3 = Matmul(Transpose(Yw3),CTYw3TpTYw3CYw3) 
 TpYx3CYx3YdadjYd = Matmul(Transpose(Yx3),CYx3YdadjYd) 
 TpYx3CYx3TpYx3CYx3 = Matmul(Transpose(Yx3),CYx3TpYx3CYx3) 
Forall(i2=1:3)  TpYx3CYx3TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3TpYx3CYx3(i2,i2),dp) 
 TpYx3CYx3TpTYx3CTYx3 = Matmul(Transpose(Yx3),CYx3TpTYx3CTYx3) 
 TpYx3CTYx3TpTYx3CYx3 = Matmul(Transpose(Yx3),CTYx3TpTYx3CYx3) 
 TpTYb3CYb3TpYb3CTYb3 = Matmul(Transpose(TYb3),CYb3TpYb3CTYb3) 
 TpTYb3CYb3TpYw3CTYw3 = Matmul(Transpose(TYb3),CYb3TpYw3CTYw3) 
 TpTYdadjYdTpYx3CTYx3 = Matmul(Transpose(TYd),adjYdTpYx3CTYx3) 
 TpTYdadjYx3Yx3CTYd = Matmul(Transpose(TYd),adjYx3Yx3CTYd) 
 TpTYdCYdTpYdCTYd = Matmul(Transpose(TYd),CYdTpYdCTYd) 
 TpTYeCYeTpYeCTYe = Matmul(Transpose(TYe),CYeTpYeCTYe) 
 TpTYuCYuTpYuCTYu = Matmul(Transpose(TYu),CYuTpYuCTYu) 
 TpTYw3CYw3TpYb3CTYb3 = Matmul(Transpose(TYw3),CYw3TpYb3CTYb3) 
 TpTYw3CYw3TpYw3CTYw3 = Matmul(Transpose(TYw3),CYw3TpYw3CTYw3) 
 TpTYx3CYx3TpYx3CTYx3 = Matmul(Transpose(TYx3),CYx3TpYx3CTYx3) 
 MBM3CYb3TpYb3CYb3TpYb3 = Matmul(MBM3,CYb3TpYb3CYb3TpYb3) 
 MBM3CYb3TpYb3CYb3TpTYb3 = Matmul(MBM3,CYb3TpYb3CYb3TpTYb3) 
 MBM3CYb3TpYeCYeTpYb3 = Matmul(MBM3,CYb3TpYeCYeTpYb3) 
 MBM3CYb3TpYeCYeTpTYb3 = Matmul(MBM3,CYb3TpYeCYeTpTYb3) 
 MBM3CYb3TpYw3CYw3TpYb3 = Matmul(MBM3,CYb3TpYw3CYw3TpYb3) 
 MBM3CYb3TpYw3CYw3TpTYb3 = Matmul(MBM3,CYb3TpYw3CYw3TpTYb3) 
 MBM3CYb3TpTYb3CYb3TpYb3 = Matmul(MBM3,CYb3TpTYb3CYb3TpYb3) 
 MBM3CYb3TpTYeCYeTpYb3 = Matmul(MBM3,CYb3TpTYeCYeTpYb3) 
 MBM3CYb3TpTYw3CYw3TpYb3 = Matmul(MBM3,CYb3TpTYw3CYw3TpYb3) 
 md2YdadjYdYdadjYd = Matmul(md2,YdadjYdYdadjYd) 
 md2YdadjYdTpYx3CYx3 = Matmul(md2,YdadjYdTpYx3CYx3) 
 md2YdadjYuYuadjYd = Matmul(md2,YdadjYuYuadjYd) 
 md2adjYx3Yx3adjYx3Yx3 = Matmul(md2,adjYx3Yx3adjYx3Yx3) 
 md2TpYx3CYx3YdadjYd = Matmul(md2,TpYx3CYx3YdadjYd) 
 md2TpYx3CYx3TpYx3CYx3 = Matmul(md2,TpYx3CYx3TpYx3CYx3) 
 me2YeadjYb3Yb3adjYe = Matmul(me2,YeadjYb3Yb3adjYe) 
 me2YeadjYeYeadjYe = Matmul(me2,YeadjYeYeadjYe) 
 me2YeadjYw3Yw3adjYe = Matmul(me2,YeadjYw3Yw3adjYe) 
 mHb32Yb3adjYb3Yb3adjYb3 = Matmul(mHb32,Yb3adjYb3Yb3adjYb3) 
 mHb32Yb3adjYeYeadjYb3 = Matmul(mHb32,Yb3adjYeYeadjYb3) 
 mHb32Yb3adjYw3Yw3adjYb3 = Matmul(mHb32,Yb3adjYw3Yw3adjYb3) 
 mHw32Yw3adjYb3Yb3adjYw3 = Matmul(mHw32,Yw3adjYb3Yb3adjYw3) 
 mHw32Yw3adjYeYeadjYw3 = Matmul(mHw32,Yw3adjYeYeadjYw3) 
 mHw32Yw3adjYw3Yw3adjYw3 = Matmul(mHw32,Yw3adjYw3Yw3adjYw3) 
 mHxb32Yx3adjYx3Yx3adjYx3 = Matmul(mHxb32,Yx3adjYx3Yx3adjYx3) 
 mHxb32Yx3CYdTpYdadjYx3 = Matmul(mHxb32,Yx3CYdTpYdadjYx3) 
 mHxb32CYx3YdadjYdTpYx3 = Matmul(mHxb32,CYx3YdadjYdTpYx3) 
 ml2adjYb3Yb3adjYb3Yb3 = Matmul(ml2,adjYb3Yb3adjYb3Yb3) 
 ml2adjYb3Yb3adjYeYe = Matmul(ml2,adjYb3Yb3adjYeYe) 
 ml2adjYb3Yb3adjYw3Yw3 = Matmul(ml2,adjYb3Yb3adjYw3Yw3) 
 ml2adjYeYeadjYb3Yb3 = Matmul(ml2,adjYeYeadjYb3Yb3) 
 ml2adjYeYeadjYeYe = Matmul(ml2,adjYeYeadjYeYe) 
 ml2adjYeYeadjYw3Yw3 = Matmul(ml2,adjYeYeadjYw3Yw3) 
 ml2adjYw3Yw3adjYb3Yb3 = Matmul(ml2,adjYw3Yw3adjYb3Yb3) 
 ml2adjYw3Yw3adjYeYe = Matmul(ml2,adjYw3Yw3adjYeYe) 
 ml2adjYw3Yw3adjYw3Yw3 = Matmul(ml2,adjYw3Yw3adjYw3Yw3) 
 ml2TpYb3CYb3TpYb3CYb3 = Matmul(ml2,TpYb3CYb3TpYb3CYb3) 
 ml2TpYb3CYb3TpYw3CYw3 = Matmul(ml2,TpYb3CYb3TpYw3CYw3) 
 ml2TpYeCYeTpYeCYe = Matmul(ml2,TpYeCYeTpYeCYe) 
 ml2TpYw3CYw3TpYb3CYb3 = Matmul(ml2,TpYw3CYw3TpYb3CYb3) 
 ml2TpYw3CYw3TpYw3CYw3 = Matmul(ml2,TpYw3CYw3TpYw3CYw3) 
 mq2adjYdYdadjYdYd = Matmul(mq2,adjYdYdadjYdYd) 
 mq2adjYdYdadjYuYu = Matmul(mq2,adjYdYdadjYuYu) 
 mq2adjYdTpYx3CYx3Yd = Matmul(mq2,adjYdTpYx3CYx3Yd) 
 mq2adjYuYuadjYdYd = Matmul(mq2,adjYuYuadjYdYd) 
 mq2adjYuYuadjYuYu = Matmul(mq2,adjYuYuadjYuYu) 
 mq2TpYdCYdTpYdCYd = Matmul(mq2,TpYdCYdTpYdCYd) 
 mq2TpYuCYuTpYuCYu = Matmul(mq2,TpYuCYuTpYuCYu) 
 mu2YuadjYdYdadjYu = Matmul(mu2,YuadjYdYdadjYu) 
 mu2YuadjYuYuadjYu = Matmul(mu2,YuadjYuYuadjYu) 
 MWM3CYw3TpYb3CYb3TpYw3 = Matmul(MWM3,CYw3TpYb3CYb3TpYw3) 
 MWM3CYw3TpYb3CYb3TpTYw3 = Matmul(MWM3,CYw3TpYb3CYb3TpTYw3) 
 MWM3CYw3TpYeCYeTpYw3 = Matmul(MWM3,CYw3TpYeCYeTpYw3) 
 MWM3CYw3TpYeCYeTpTYw3 = Matmul(MWM3,CYw3TpYeCYeTpTYw3) 
 MWM3CYw3TpYw3CYw3TpYw3 = Matmul(MWM3,CYw3TpYw3CYw3TpYw3) 
 MWM3CYw3TpYw3CYw3TpTYw3 = Matmul(MWM3,CYw3TpYw3CYw3TpTYw3) 
 MWM3CYw3TpTYb3CYb3TpYw3 = Matmul(MWM3,CYw3TpTYb3CYb3TpYw3) 
 MWM3CYw3TpTYeCYeTpYw3 = Matmul(MWM3,CYw3TpTYeCYeTpYw3) 
 MWM3CYw3TpTYw3CYw3TpYw3 = Matmul(MWM3,CYw3TpTYw3CYw3TpYw3) 
 MXM3CYx3YdadjYdTpYx3 = Matmul(MXM3,CYx3YdadjYdTpYx3) 
 MXM3CYx3YdadjYdTpTYx3 = Matmul(MXM3,CYx3YdadjYdTpTYx3) 
 MXM3CYx3TYdadjYdTpYx3 = Matmul(MXM3,CYx3TYdadjYdTpYx3) 
 MXM3CYx3TpYx3CYx3TpYx3 = Matmul(MXM3,CYx3TpYx3CYx3TpYx3) 
 MXM3CYx3TpYx3CYx3TpTYx3 = Matmul(MXM3,CYx3TpYx3CYx3TpTYx3) 
 MXM3CYx3TpTYx3CYx3TpYx3 = Matmul(MXM3,CYx3TpTYx3CYx3TpYx3) 
 Yb3ml2adjYb3Yb3adjYb3 = Matmul(Yb3,ml2adjYb3Yb3adjYb3) 
 Yb3ml2adjYeYeadjYb3 = Matmul(Yb3,ml2adjYeYeadjYb3) 
 Yb3ml2adjYw3Yw3adjYb3 = Matmul(Yb3,ml2adjYw3Yw3adjYb3) 
 Yb3adjYb3mHb32Yb3adjYb3 = Matmul(Yb3,adjYb3mHb32Yb3adjYb3) 
 Yb3adjYb3Yb3ml2adjYb3 = Matmul(Yb3,adjYb3Yb3ml2adjYb3) 
 Yb3adjYb3Yb3adjYb3MBM3 = Matmul(Yb3,adjYb3Yb3adjYb3MBM3) 
 Yb3adjYb3Yb3adjYb3mHb32 = Matmul(Yb3,adjYb3Yb3adjYb3mHb32) 
 Yb3adjYb3Yb3adjYb3Yb3 = Matmul(Yb3,adjYb3Yb3adjYb3Yb3) 
! Yb3adjYb3Yb3adjYb3BMBM3 = Matmul(Yb3,adjYb3Yb3adjYb3BMBM3) 
 Yb3adjYb3Yb3adjYb3TYb3 = Matmul(Yb3,adjYb3Yb3adjYb3TYb3) 
 Yb3adjYb3Yb3adjYw3Yw3 = Matmul(Yb3,adjYb3Yb3adjYw3Yw3) 
 Yb3adjYb3Yb3adjYw3TYw3 = Matmul(Yb3,adjYb3Yb3adjYw3TYw3) 
 Yb3adjYb3TYb3adjYb3MBM3 = Matmul(Yb3,adjYb3TYb3adjYb3MBM3) 
 Yb3adjYb3TYb3adjYb3Yb3 = Matmul(Yb3,adjYb3TYb3adjYb3Yb3) 
 Yb3adjYb3TYb3adjYw3Yw3 = Matmul(Yb3,adjYb3TYb3adjYw3Yw3) 
 Yb3adjYeme2YeadjYb3 = Matmul(Yb3,adjYeme2YeadjYb3) 
 Yb3adjYeYeml2adjYb3 = Matmul(Yb3,adjYeYeml2adjYb3) 
 Yb3adjYeYeadjYb3MBM3 = Matmul(Yb3,adjYeYeadjYb3MBM3) 
 Yb3adjYeYeadjYb3mHb32 = Matmul(Yb3,adjYeYeadjYb3mHb32) 
 Yb3adjYeYeadjYb3Yb3 = Matmul(Yb3,adjYeYeadjYb3Yb3) 
! Yb3adjYeYeadjYb3BMBM3 = Matmul(Yb3,adjYeYeadjYb3BMBM3) 
 Yb3adjYeYeadjYb3TYb3 = Matmul(Yb3,adjYeYeadjYb3TYb3) 
 Yb3adjYeYeadjYeYe = Matmul(Yb3,adjYeYeadjYeYe) 
 Yb3adjYeYeadjYeTYe = Matmul(Yb3,adjYeYeadjYeTYe) 
 Yb3adjYeYeadjYw3TYw3 = Matmul(Yb3,adjYeYeadjYw3TYw3) 
 Yb3adjYeTYeadjYb3MBM3 = Matmul(Yb3,adjYeTYeadjYb3MBM3) 
 Yb3adjYeTYeadjYb3Yb3 = Matmul(Yb3,adjYeTYeadjYb3Yb3) 
 Yb3adjYeTYeadjYeYe = Matmul(Yb3,adjYeTYeadjYeYe) 
 Yb3adjYeTYeadjYw3Yw3 = Matmul(Yb3,adjYeTYeadjYw3Yw3) 
 Yb3adjYw3mHw32Yw3adjYb3 = Matmul(Yb3,adjYw3mHw32Yw3adjYb3) 
 Yb3adjYw3Yw3ml2adjYb3 = Matmul(Yb3,adjYw3Yw3ml2adjYb3) 
 Yb3adjYw3Yw3adjYb3MBM3 = Matmul(Yb3,adjYw3Yw3adjYb3MBM3) 
 Yb3adjYw3Yw3adjYb3mHb32 = Matmul(Yb3,adjYw3Yw3adjYb3mHb32) 
 Yb3adjYw3Yw3adjYb3Yb3 = Matmul(Yb3,adjYw3Yw3adjYb3Yb3) 
! Yb3adjYw3Yw3adjYb3BMBM3 = Matmul(Yb3,adjYw3Yw3adjYb3BMBM3) 
 Yb3adjYw3Yw3adjYb3TYb3 = Matmul(Yb3,adjYw3Yw3adjYb3TYb3) 
 Yb3adjYw3Yw3adjYw3Yw3 = Matmul(Yb3,adjYw3Yw3adjYw3Yw3) 
 Yb3adjYw3Yw3adjYw3TYw3 = Matmul(Yb3,adjYw3Yw3adjYw3TYw3) 
 Yb3adjYw3TYw3adjYb3MBM3 = Matmul(Yb3,adjYw3TYw3adjYb3MBM3) 
 Yb3adjYw3TYw3adjYb3Yb3 = Matmul(Yb3,adjYw3TYw3adjYb3Yb3) 
 Yb3adjYw3TYw3adjYw3Yw3 = Matmul(Yb3,adjYw3TYw3adjYw3Yw3) 
 Ydmq2adjYdYdadjYd = Matmul(Yd,mq2adjYdYdadjYd) 
 Ydmq2adjYuYuadjYd = Matmul(Yd,mq2adjYuYuadjYd) 
 Ydmq2adjYx3Yx3CYd = Matmul(Yd,mq2adjYx3Yx3CYd) 
 YdadjYdmd2YdadjYd = Matmul(Yd,adjYdmd2YdadjYd) 
 YdadjYdYdmq2adjYd = Matmul(Yd,adjYdYdmq2adjYd) 
 YdadjYdYdadjYdmd2 = Matmul(Yd,adjYdYdadjYdmd2) 
 YdadjYdYdadjYdYd = Matmul(Yd,adjYdYdadjYdYd) 
 YdadjYdYdadjYdTYd = Matmul(Yd,adjYdYdadjYdTYd) 
 YdadjYdTYdadjYdYd = Matmul(Yd,adjYdTYdadjYdYd) 
 YdadjYdTpYx3CYx3Yd = Matmul(Yd,adjYdTpYx3CYx3Yd) 
 YdadjYdTpYx3CYx3TYd = Matmul(Yd,adjYdTpYx3CYx3TYd) 
 YdadjYdTpTYx3CYx3Yd = Matmul(Yd,adjYdTpTYx3CYx3Yd) 
 YdadjYumu2YuadjYd = Matmul(Yd,adjYumu2YuadjYd) 
 YdadjYuYumq2adjYd = Matmul(Yd,adjYuYumq2adjYd) 
 YdadjYuYuadjYdmd2 = Matmul(Yd,adjYuYuadjYdmd2) 
 YdadjYuYuadjYdYd = Matmul(Yd,adjYuYuadjYdYd) 
 YdadjYuYuadjYdTYd = Matmul(Yd,adjYuYuadjYdTYd) 
 YdadjYuYuadjYuYu = Matmul(Yd,adjYuYuadjYuYu) 
 YdadjYuYuadjYuTYu = Matmul(Yd,adjYuYuadjYuTYu) 
 YdadjYuTYuadjYdYd = Matmul(Yd,adjYuTYuadjYdYd) 
 YdadjYuTYuadjYuYu = Matmul(Yd,adjYuTYuadjYuYu) 
 Yeml2adjYb3Yb3adjYe = Matmul(Ye,ml2adjYb3Yb3adjYe) 
 Yeml2adjYeYeadjYe = Matmul(Ye,ml2adjYeYeadjYe) 
 Yeml2adjYw3Yw3adjYe = Matmul(Ye,ml2adjYw3Yw3adjYe) 
 YeadjYb3mHb32Yb3adjYe = Matmul(Ye,adjYb3mHb32Yb3adjYe) 
 YeadjYb3Yb3ml2adjYe = Matmul(Ye,adjYb3Yb3ml2adjYe) 
 YeadjYb3Yb3adjYb3Yb3 = Matmul(Ye,adjYb3Yb3adjYb3Yb3) 
 YeadjYb3Yb3adjYb3TYb3 = Matmul(Ye,adjYb3Yb3adjYb3TYb3) 
 YeadjYb3Yb3adjYeme2 = Matmul(Ye,adjYb3Yb3adjYeme2) 
 YeadjYb3Yb3adjYeYe = Matmul(Ye,adjYb3Yb3adjYeYe) 
 YeadjYb3Yb3adjYeTYe = Matmul(Ye,adjYb3Yb3adjYeTYe) 
 YeadjYb3Yb3adjYw3Yw3 = Matmul(Ye,adjYb3Yb3adjYw3Yw3) 
 YeadjYb3Yb3adjYw3TYw3 = Matmul(Ye,adjYb3Yb3adjYw3TYw3) 
 YeadjYb3TYb3adjYb3Yb3 = Matmul(Ye,adjYb3TYb3adjYb3Yb3) 
 YeadjYb3TYb3adjYeYe = Matmul(Ye,adjYb3TYb3adjYeYe) 
 YeadjYb3TYb3adjYw3Yw3 = Matmul(Ye,adjYb3TYb3adjYw3Yw3) 
 YeadjYeme2YeadjYe = Matmul(Ye,adjYeme2YeadjYe) 
 YeadjYeYeml2adjYe = Matmul(Ye,adjYeYeml2adjYe) 
 YeadjYeYeadjYeme2 = Matmul(Ye,adjYeYeadjYeme2) 
 YeadjYeYeadjYeYe = Matmul(Ye,adjYeYeadjYeYe) 
 YeadjYeYeadjYeTYe = Matmul(Ye,adjYeYeadjYeTYe) 
 YeadjYeTYeadjYeYe = Matmul(Ye,adjYeTYeadjYeYe) 
 YeadjYw3mHw32Yw3adjYe = Matmul(Ye,adjYw3mHw32Yw3adjYe) 
 YeadjYw3Yw3ml2adjYe = Matmul(Ye,adjYw3Yw3ml2adjYe) 
 YeadjYw3Yw3adjYb3Yb3 = Matmul(Ye,adjYw3Yw3adjYb3Yb3) 
 YeadjYw3Yw3adjYb3TYb3 = Matmul(Ye,adjYw3Yw3adjYb3TYb3) 
 YeadjYw3Yw3adjYeme2 = Matmul(Ye,adjYw3Yw3adjYeme2) 
 YeadjYw3Yw3adjYeYe = Matmul(Ye,adjYw3Yw3adjYeYe) 
 YeadjYw3Yw3adjYeTYe = Matmul(Ye,adjYw3Yw3adjYeTYe) 
 YeadjYw3Yw3adjYw3Yw3 = Matmul(Ye,adjYw3Yw3adjYw3Yw3) 
 YeadjYw3Yw3adjYw3TYw3 = Matmul(Ye,adjYw3Yw3adjYw3TYw3) 
 YeadjYw3TYw3adjYb3Yb3 = Matmul(Ye,adjYw3TYw3adjYb3Yb3) 
 YeadjYw3TYw3adjYeYe = Matmul(Ye,adjYw3TYw3adjYeYe) 
 YeadjYw3TYw3adjYw3Yw3 = Matmul(Ye,adjYw3TYw3adjYw3Yw3) 
 Yumq2adjYdYdadjYu = Matmul(Yu,mq2adjYdYdadjYu) 
 Yumq2adjYuYuadjYu = Matmul(Yu,mq2adjYuYuadjYu) 
 YuadjYdmd2YdadjYu = Matmul(Yu,adjYdmd2YdadjYu) 
 YuadjYdYdmq2adjYu = Matmul(Yu,adjYdYdmq2adjYu) 
 YuadjYdYdadjYdYd = Matmul(Yu,adjYdYdadjYdYd) 
 YuadjYdYdadjYdTYd = Matmul(Yu,adjYdYdadjYdTYd) 
 YuadjYdYdadjYumu2 = Matmul(Yu,adjYdYdadjYumu2) 
 YuadjYdYdadjYuYu = Matmul(Yu,adjYdYdadjYuYu) 
 YuadjYdYdadjYuTYu = Matmul(Yu,adjYdYdadjYuTYu) 
 YuadjYdTYdadjYdYd = Matmul(Yu,adjYdTYdadjYdYd) 
 YuadjYdTYdadjYuYu = Matmul(Yu,adjYdTYdadjYuYu) 
 YuadjYdTpYx3CYx3Yd = Matmul(Yu,adjYdTpYx3CYx3Yd) 
 YuadjYdTpYx3CYx3TYd = Matmul(Yu,adjYdTpYx3CYx3TYd) 
 YuadjYdTpTYx3CYx3Yd = Matmul(Yu,adjYdTpTYx3CYx3Yd) 
 YuadjYumu2YuadjYu = Matmul(Yu,adjYumu2YuadjYu) 
 YuadjYuYumq2adjYu = Matmul(Yu,adjYuYumq2adjYu) 
 YuadjYuYuadjYumu2 = Matmul(Yu,adjYuYuadjYumu2) 
 YuadjYuYuadjYuYu = Matmul(Yu,adjYuYuadjYuYu) 
 YuadjYuYuadjYuTYu = Matmul(Yu,adjYuYuadjYuTYu) 
 YuadjYuTYuadjYuYu = Matmul(Yu,adjYuTYuadjYuYu) 
 Yw3ml2adjYb3Yb3adjYw3 = Matmul(Yw3,ml2adjYb3Yb3adjYw3) 
 Yw3ml2adjYeYeadjYw3 = Matmul(Yw3,ml2adjYeYeadjYw3) 
 Yw3ml2adjYw3Yw3adjYw3 = Matmul(Yw3,ml2adjYw3Yw3adjYw3) 
 Yw3adjYb3mHb32Yb3adjYw3 = Matmul(Yw3,adjYb3mHb32Yb3adjYw3) 
 Yw3adjYb3Yb3ml2adjYw3 = Matmul(Yw3,adjYb3Yb3ml2adjYw3) 
 Yw3adjYb3Yb3adjYb3Yb3 = Matmul(Yw3,adjYb3Yb3adjYb3Yb3) 
 Yw3adjYb3Yb3adjYb3TYb3 = Matmul(Yw3,adjYb3Yb3adjYb3TYb3) 
 Yw3adjYb3Yb3adjYw3mHw32 = Matmul(Yw3,adjYb3Yb3adjYw3mHw32) 
 Yw3adjYb3Yb3adjYw3MWM3 = Matmul(Yw3,adjYb3Yb3adjYw3MWM3) 
 Yw3adjYb3Yb3adjYw3Yw3 = Matmul(Yw3,adjYb3Yb3adjYw3Yw3) 
! Yw3adjYb3Yb3adjYw3BMWM3 = Matmul(Yw3,adjYb3Yb3adjYw3BMWM3) 
 Yw3adjYb3Yb3adjYw3TYw3 = Matmul(Yw3,adjYb3Yb3adjYw3TYw3) 
 Yw3adjYb3TYb3adjYb3Yb3 = Matmul(Yw3,adjYb3TYb3adjYb3Yb3) 
 Yw3adjYb3TYb3adjYw3MWM3 = Matmul(Yw3,adjYb3TYb3adjYw3MWM3) 
 Yw3adjYb3TYb3adjYw3Yw3 = Matmul(Yw3,adjYb3TYb3adjYw3Yw3) 
 Yw3adjYeme2YeadjYw3 = Matmul(Yw3,adjYeme2YeadjYw3) 
 Yw3adjYeYeml2adjYw3 = Matmul(Yw3,adjYeYeml2adjYw3) 
 Yw3adjYeYeadjYb3TYb3 = Matmul(Yw3,adjYeYeadjYb3TYb3) 
 Yw3adjYeYeadjYeYe = Matmul(Yw3,adjYeYeadjYeYe) 
 Yw3adjYeYeadjYeTYe = Matmul(Yw3,adjYeYeadjYeTYe) 
 Yw3adjYeYeadjYw3mHw32 = Matmul(Yw3,adjYeYeadjYw3mHw32) 
 Yw3adjYeYeadjYw3MWM3 = Matmul(Yw3,adjYeYeadjYw3MWM3) 
 Yw3adjYeYeadjYw3Yw3 = Matmul(Yw3,adjYeYeadjYw3Yw3) 
! Yw3adjYeYeadjYw3BMWM3 = Matmul(Yw3,adjYeYeadjYw3BMWM3) 
 Yw3adjYeYeadjYw3TYw3 = Matmul(Yw3,adjYeYeadjYw3TYw3) 
 Yw3adjYeTYeadjYb3Yb3 = Matmul(Yw3,adjYeTYeadjYb3Yb3) 
 Yw3adjYeTYeadjYeYe = Matmul(Yw3,adjYeTYeadjYeYe) 
 Yw3adjYeTYeadjYw3MWM3 = Matmul(Yw3,adjYeTYeadjYw3MWM3) 
 Yw3adjYeTYeadjYw3Yw3 = Matmul(Yw3,adjYeTYeadjYw3Yw3) 
 Yw3adjYw3mHw32Yw3adjYw3 = Matmul(Yw3,adjYw3mHw32Yw3adjYw3) 
 Yw3adjYw3Yw3ml2adjYw3 = Matmul(Yw3,adjYw3Yw3ml2adjYw3) 
 Yw3adjYw3Yw3adjYb3Yb3 = Matmul(Yw3,adjYw3Yw3adjYb3Yb3) 
 Yw3adjYw3Yw3adjYb3TYb3 = Matmul(Yw3,adjYw3Yw3adjYb3TYb3) 
 Yw3adjYw3Yw3adjYw3mHw32 = Matmul(Yw3,adjYw3Yw3adjYw3mHw32) 
 Yw3adjYw3Yw3adjYw3MWM3 = Matmul(Yw3,adjYw3Yw3adjYw3MWM3) 
 Yw3adjYw3Yw3adjYw3Yw3 = Matmul(Yw3,adjYw3Yw3adjYw3Yw3) 
! Yw3adjYw3Yw3adjYw3BMWM3 = Matmul(Yw3,adjYw3Yw3adjYw3BMWM3) 
 Yw3adjYw3Yw3adjYw3TYw3 = Matmul(Yw3,adjYw3Yw3adjYw3TYw3) 
 Yw3adjYw3TYw3adjYb3Yb3 = Matmul(Yw3,adjYw3TYw3adjYb3Yb3) 
 Yw3adjYw3TYw3adjYw3MWM3 = Matmul(Yw3,adjYw3TYw3adjYw3MWM3) 
 Yw3adjYw3TYw3adjYw3Yw3 = Matmul(Yw3,adjYw3TYw3adjYw3Yw3) 
 Yx3md2adjYx3Yx3adjYx3 = Matmul(Yx3,md2adjYx3Yx3adjYx3) 
 Yx3md2CYdTpYdadjYx3 = Matmul(Yx3,md2CYdTpYdadjYx3) 
 Yx3adjYx3mHxb32Yx3adjYx3 = Matmul(Yx3,adjYx3mHxb32Yx3adjYx3) 
 Yx3adjYx3Yx3md2adjYx3 = Matmul(Yx3,adjYx3Yx3md2adjYx3) 
 Yx3adjYx3Yx3adjYx3mHxb32 = Matmul(Yx3,adjYx3Yx3adjYx3mHxb32) 
 Yx3adjYx3Yx3adjYx3Yx3 = Matmul(Yx3,adjYx3Yx3adjYx3Yx3) 
 Yx3adjYx3Yx3adjYx3TYx3 = Matmul(Yx3,adjYx3Yx3adjYx3TYx3) 
 Yx3adjYx3TYx3adjYx3Yx3 = Matmul(Yx3,adjYx3TYx3adjYx3Yx3) 
 Yx3CYdmq2TpYdadjYx3 = Matmul(Yx3,CYdmq2TpYdadjYx3) 
 Yx3CYdTpYdmd2adjYx3 = Matmul(Yx3,CYdTpYdmd2adjYx3) 
 Yx3CYdTpYdadjYx3mHxb32 = Matmul(Yx3,CYdTpYdadjYx3mHxb32) 
 Yx3CYdTpYdadjYx3Yx3 = Matmul(Yx3,CYdTpYdadjYx3Yx3) 
 Yx3CYdTpYdadjYx3TYx3 = Matmul(Yx3,CYdTpYdadjYx3TYx3) 
 Yx3CYdTpYdCYdTpYd = Matmul(Yx3,CYdTpYdCYdTpYd) 
 Yx3CYdTpYdCYdTpTYd = Matmul(Yx3,CYdTpYdCYdTpTYd) 
 Yx3CYdTpYuCYuTpYd = Matmul(Yx3,CYdTpYuCYuTpYd) 
 Yx3CYdTpYuCYuTpTYd = Matmul(Yx3,CYdTpYuCYuTpTYd) 
 Yx3CYdTpTYdCYdTpYd = Matmul(Yx3,CYdTpTYdCYdTpYd) 
 Yx3CYdTpTYuCYuTpYd = Matmul(Yx3,CYdTpTYuCYuTpYd) 
 Yx3TYdadjYdadjYx3Yx3 = Matmul(Yx3,TYdadjYdadjYx3Yx3) 
! BMBM3CYb3TpYb3CYb3TpYb3 = Matmul(BMBM3,CYb3TpYb3CYb3TpYb3) 
! BMBM3CYb3TpYeCYeTpYb3 = Matmul(BMBM3,CYb3TpYeCYeTpYb3) 
! BMBM3CYb3TpYw3CYw3TpYb3 = Matmul(BMBM3,CYb3TpYw3CYw3TpYb3) 
! BMWM3CYw3TpYb3CYb3TpYw3 = Matmul(BMWM3,CYw3TpYb3CYb3TpYw3) 
! BMWM3CYw3TpYeCYeTpYw3 = Matmul(BMWM3,CYw3TpYeCYeTpYw3) 
! BMWM3CYw3TpYw3CYw3TpYw3 = Matmul(BMWM3,CYw3TpYw3CYw3TpYw3) 
! BMXM3CYx3YdadjYdTpYx3 = Matmul(BMXM3,CYx3YdadjYdTpYx3) 
! BMXM3CYx3TpYx3CYx3TpYx3 = Matmul(BMXM3,CYx3TpYx3CYx3TpYx3) 
 TYb3adjYb3Yb3adjYb3MBM3 = Matmul(TYb3,adjYb3Yb3adjYb3MBM3) 
 TYb3adjYb3Yb3adjYb3Yb3 = Matmul(TYb3,adjYb3Yb3adjYb3Yb3) 
 TYb3adjYb3Yb3adjYw3Yw3 = Matmul(TYb3,adjYb3Yb3adjYw3Yw3) 
 TYb3adjYeYeadjYb3MBM3 = Matmul(TYb3,adjYeYeadjYb3MBM3) 
 TYb3adjYeYeadjYb3Yb3 = Matmul(TYb3,adjYeYeadjYb3Yb3) 
 TYb3adjYeYeadjYeYe = Matmul(TYb3,adjYeYeadjYeYe) 
 TYb3adjYeYeadjYw3Yw3 = Matmul(TYb3,adjYeYeadjYw3Yw3) 
 TYb3adjYw3Yw3adjYb3MBM3 = Matmul(TYb3,adjYw3Yw3adjYb3MBM3) 
 TYb3adjYw3Yw3adjYb3Yb3 = Matmul(TYb3,adjYw3Yw3adjYb3Yb3) 
 TYb3adjYw3Yw3adjYw3Yw3 = Matmul(TYb3,adjYw3Yw3adjYw3Yw3) 
 TYdadjYdYdadjYdYd = Matmul(TYd,adjYdYdadjYdYd) 
 TYdadjYdTpYx3CYx3Yd = Matmul(TYd,adjYdTpYx3CYx3Yd) 
 TYdadjYuYuadjYdYd = Matmul(TYd,adjYuYuadjYdYd) 
 TYdadjYuYuadjYuYu = Matmul(TYd,adjYuYuadjYuYu) 
 TYeadjYb3Yb3adjYb3Yb3 = Matmul(TYe,adjYb3Yb3adjYb3Yb3) 
 TYeadjYb3Yb3adjYeYe = Matmul(TYe,adjYb3Yb3adjYeYe) 
 TYeadjYb3Yb3adjYw3Yw3 = Matmul(TYe,adjYb3Yb3adjYw3Yw3) 
 TYeadjYeYeadjYeYe = Matmul(TYe,adjYeYeadjYeYe) 
 TYeadjYw3Yw3adjYb3Yb3 = Matmul(TYe,adjYw3Yw3adjYb3Yb3) 
 TYeadjYw3Yw3adjYeYe = Matmul(TYe,adjYw3Yw3adjYeYe) 
 TYeadjYw3Yw3adjYw3Yw3 = Matmul(TYe,adjYw3Yw3adjYw3Yw3) 
 TYuadjYdYdadjYdYd = Matmul(TYu,adjYdYdadjYdYd) 
 TYuadjYdYdadjYuYu = Matmul(TYu,adjYdYdadjYuYu) 
 TYuadjYdTpYx3CYx3Yd = Matmul(TYu,adjYdTpYx3CYx3Yd) 
 TYuadjYuYuadjYuYu = Matmul(TYu,adjYuYuadjYuYu) 
 TYw3adjYb3Yb3adjYb3Yb3 = Matmul(TYw3,adjYb3Yb3adjYb3Yb3) 
 TYw3adjYb3Yb3adjYw3MWM3 = Matmul(TYw3,adjYb3Yb3adjYw3MWM3) 
 TYw3adjYb3Yb3adjYw3Yw3 = Matmul(TYw3,adjYb3Yb3adjYw3Yw3) 
 TYw3adjYeYeadjYb3Yb3 = Matmul(TYw3,adjYeYeadjYb3Yb3) 
 TYw3adjYeYeadjYeYe = Matmul(TYw3,adjYeYeadjYeYe) 
 TYw3adjYeYeadjYw3MWM3 = Matmul(TYw3,adjYeYeadjYw3MWM3) 
 TYw3adjYeYeadjYw3Yw3 = Matmul(TYw3,adjYeYeadjYw3Yw3) 
 TYw3adjYw3Yw3adjYb3Yb3 = Matmul(TYw3,adjYw3Yw3adjYb3Yb3) 
 TYw3adjYw3Yw3adjYw3MWM3 = Matmul(TYw3,adjYw3Yw3adjYw3MWM3) 
 TYw3adjYw3Yw3adjYw3Yw3 = Matmul(TYw3,adjYw3Yw3adjYw3Yw3) 
 TYx3adjYx3Yx3adjYx3Yx3 = Matmul(TYx3,adjYx3Yx3adjYx3Yx3) 
 TYx3CYdTpYdadjYx3Yx3 = Matmul(TYx3,CYdTpYdadjYx3Yx3) 
 TYx3CYdTpYdCYdTpYd = Matmul(TYx3,CYdTpYdCYdTpYd) 
 TYx3CYdTpYuCYuTpYd = Matmul(TYx3,CYdTpYuCYuTpYd) 
 TpYb3mHb32CYb3TpYb3CYb3 = Matmul(Transpose(Yb3),mHb32CYb3TpYb3CYb3) 
 TpYb3mHb32CYb3TpYw3CYw3 = Matmul(Transpose(Yb3),mHb32CYb3TpYw3CYw3) 
 TpYb3CYb3ml2TpYb3CYb3 = Matmul(Transpose(Yb3),CYb3ml2TpYb3CYb3) 
 TpYb3CYb3ml2TpYw3CYw3 = Matmul(Transpose(Yb3),CYb3ml2TpYw3CYw3) 
 TpYb3CYb3TpYb3mHb32CYb3 = Matmul(Transpose(Yb3),CYb3TpYb3mHb32CYb3) 
 TpYb3CYb3TpYb3CYb3ml2 = Matmul(Transpose(Yb3),CYb3TpYb3CYb3ml2) 
 TpYb3CYb3TpYw3mHw32CYw3 = Matmul(Transpose(Yb3),CYb3TpYw3mHw32CYw3) 
 TpYb3CYb3TpYw3CYw3ml2 = Matmul(Transpose(Yb3),CYb3TpYw3CYw3ml2) 
 TpYdmd2adjYx3Yx3CYd = Matmul(Transpose(Yd),md2adjYx3Yx3CYd) 
 TpYdmd2CYdTpYdCYd = Matmul(Transpose(Yd),md2CYdTpYdCYd) 
 TpYdadjYx3mHxb32Yx3CYd = Matmul(Transpose(Yd),adjYx3mHxb32Yx3CYd) 
 TpYdadjYx3Yx3md2CYd = Matmul(Transpose(Yd),adjYx3Yx3md2CYd) 
 TpYdadjYx3Yx3CYdmq2 = Matmul(Transpose(Yd),adjYx3Yx3CYdmq2) 
 TpYdCYdmq2TpYdCYd = Matmul(Transpose(Yd),CYdmq2TpYdCYd) 
 TpYdCYdTpYdmd2CYd = Matmul(Transpose(Yd),CYdTpYdmd2CYd) 
 TpYdCYdTpYdCYdmq2 = Matmul(Transpose(Yd),CYdTpYdCYdmq2) 
 TpYeme2CYeTpYeCYe = Matmul(Transpose(Ye),me2CYeTpYeCYe) 
 TpYeCYeml2TpYeCYe = Matmul(Transpose(Ye),CYeml2TpYeCYe) 
 TpYeCYeTpYeme2CYe = Matmul(Transpose(Ye),CYeTpYeme2CYe) 
 TpYeCYeTpYeCYeml2 = Matmul(Transpose(Ye),CYeTpYeCYeml2) 
 TpYumu2CYuTpYuCYu = Matmul(Transpose(Yu),mu2CYuTpYuCYu) 
 TpYuCYumq2TpYuCYu = Matmul(Transpose(Yu),CYumq2TpYuCYu) 
 TpYuCYuTpYumu2CYu = Matmul(Transpose(Yu),CYuTpYumu2CYu) 
 TpYuCYuTpYuCYumq2 = Matmul(Transpose(Yu),CYuTpYuCYumq2) 
 TpYw3mHw32CYw3TpYb3CYb3 = Matmul(Transpose(Yw3),mHw32CYw3TpYb3CYb3) 
 TpYw3mHw32CYw3TpYw3CYw3 = Matmul(Transpose(Yw3),mHw32CYw3TpYw3CYw3) 
 TpYw3CYw3ml2TpYb3CYb3 = Matmul(Transpose(Yw3),CYw3ml2TpYb3CYb3) 
 TpYw3CYw3ml2TpYw3CYw3 = Matmul(Transpose(Yw3),CYw3ml2TpYw3CYw3) 
 TpYw3CYw3TpYb3mHb32CYb3 = Matmul(Transpose(Yw3),CYw3TpYb3mHb32CYb3) 
 TpYw3CYw3TpYb3CYb3ml2 = Matmul(Transpose(Yw3),CYw3TpYb3CYb3ml2) 
 TpYw3CYw3TpYw3mHw32CYw3 = Matmul(Transpose(Yw3),CYw3TpYw3mHw32CYw3) 
 TpYw3CYw3TpYw3CYw3ml2 = Matmul(Transpose(Yw3),CYw3TpYw3CYw3ml2) 
 TpYx3mHxb32CYx3TpYx3CYx3 = Matmul(Transpose(Yx3),mHxb32CYx3TpYx3CYx3) 
 TpYx3CYx3md2TpYx3CYx3 = Matmul(Transpose(Yx3),CYx3md2TpYx3CYx3) 
 TpYx3CYx3TpYx3mHxb32CYx3 = Matmul(Transpose(Yx3),CYx3TpYx3mHxb32CYx3) 
 TpYx3CYx3TpYx3CYx3md2 = Matmul(Transpose(Yx3),CYx3TpYx3CYx3md2) 
 TpYx3CYx3TpYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3TpYx3CYx3Yd) 
 TpYx3CYx3TpYx3CYx3TYd = Matmul(Transpose(Yx3),CYx3TpYx3CYx3TYd) 
 TpYx3CYx3TpTYx3CYx3Yd = Matmul(Transpose(Yx3),CYx3TpTYx3CYx3Yd) 
 TpTYx3CYx3TpYx3CYx3Yd = Matmul(Transpose(TYx3),CYx3TpYx3CYx3Yd) 
 TrmHg32 = cTrace(mHg32) 
 TrmHw32 = cTrace(mHw32) 
 TrCTYb3TpYb3 = cTrace(CTYb3TpYb3) 
 TrCTYdTpYd = cTrace(CTYdTpYd) 
 TrCTYeTpYe = cTrace(CTYeTpYe) 
 TrCTYuTpYu = cTrace(CTYuTpYu) 
 TrCTYw3TpYw3 = cTrace(CTYw3TpYw3) 
 TrCTYx3TpYx3 = cTrace(CTYx3TpYx3) 
 TrYb3adjYb3Yb3adjYb3 = cTrace(Yb3adjYb3Yb3adjYb3) 
 TrYb3adjYb3TYb3adjYb3 = cTrace(Yb3adjYb3TYb3adjYb3) 
 TrYb3adjYb3TYb3adjTYb3 = cTrace(Yb3adjYb3TYb3adjTYb3) 
 TrYb3adjYeYeadjYb3 = cTrace(Yb3adjYeYeadjYb3) 
 TrYb3adjYeTYeadjYb3 = cTrace(Yb3adjYeTYeadjYb3) 
 TrYb3adjYeTYeadjTYb3 = cTrace(Yb3adjYeTYeadjTYb3) 
 TrYb3adjYw3Yw3adjYb3 = cTrace(Yb3adjYw3Yw3adjYb3) 
 TrYb3adjYw3TYw3adjYb3 = cTrace(Yb3adjYw3TYw3adjYb3) 
 TrYb3adjYw3TYw3adjTYb3 = cTrace(Yb3adjYw3TYw3adjTYb3) 
 TrYb3adjTYb3TYb3adjYb3 = cTrace(Yb3adjTYb3TYb3adjYb3) 
 TrYb3adjTYeTYeadjYb3 = cTrace(Yb3adjTYeTYeadjYb3) 
 TrYb3adjTYw3TYw3adjYb3 = cTrace(Yb3adjTYw3TYw3adjYb3) 
 TrYb3TpTYb3CTYb3adjYb3 = cTrace(Yb3TpTYb3CTYb3adjYb3) 
 TrYb3TpTYeCTYeadjYb3 = cTrace(Yb3TpTYeCTYeadjYb3) 
 TrYb3TpTYw3CTYw3adjYb3 = cTrace(Yb3TpTYw3CTYw3adjYb3) 
 TrYdadjYdYdadjYd = cTrace(YdadjYdYdadjYd) 
 TrYdadjYdTYdadjYd = cTrace(YdadjYdTYdadjYd) 
 TrYdadjYdTYdadjTYd = cTrace(YdadjYdTYdadjTYd) 
 TrYdadjYdTpYx3CYx3 = cTrace(YdadjYdTpYx3CYx3) 
 TrYdadjYdTpTYx3CTYx3 = cTrace(YdadjYdTpTYx3CTYx3) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYdadjYuTYuadjYd = cTrace(YdadjYuTYuadjYd) 
 TrYdadjYuTYuadjTYd = cTrace(YdadjYuTYuadjTYd) 
 TrYdadjTYdTYdadjYd = cTrace(YdadjTYdTYdadjYd) 
 TrYdadjTYuTYuadjYd = cTrace(YdadjTYuTYuadjYd) 
 TrYdTpTYdCTYdadjYd = cTrace(YdTpTYdCTYdadjYd) 
 TrYdTpTYuCTYuadjYd = cTrace(YdTpTYuCTYuadjYd) 
 TrYeadjYb3TYb3adjYe = cTrace(YeadjYb3TYb3adjYe) 
 TrYeadjYb3TYb3adjTYe = cTrace(YeadjYb3TYb3adjTYe) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYeTYeadjYe = cTrace(YeadjYeTYeadjYe) 
 TrYeadjYeTYeadjTYe = cTrace(YeadjYeTYeadjTYe) 
 TrYeadjYw3Yw3adjYe = cTrace(YeadjYw3Yw3adjYe) 
 TrYeadjYw3TYw3adjYe = cTrace(YeadjYw3TYw3adjYe) 
 TrYeadjYw3TYw3adjTYe = cTrace(YeadjYw3TYw3adjTYe) 
 TrYeadjTYb3TYb3adjYe = cTrace(YeadjTYb3TYb3adjYe) 
 TrYeadjTYeTYeadjYe = cTrace(YeadjTYeTYeadjYe) 
 TrYeadjTYw3TYw3adjYe = cTrace(YeadjTYw3TYw3adjYe) 
 TrYeTpTYb3CTYb3adjYe = cTrace(YeTpTYb3CTYb3adjYe) 
 TrYeTpTYeCTYeadjYe = cTrace(YeTpTYeCTYeadjYe) 
 TrYeTpTYw3CTYw3adjYe = cTrace(YeTpTYw3CTYw3adjYe) 
 TrYuadjYdTYdadjYu = cTrace(YuadjYdTYdadjYu) 
 TrYuadjYdTYdadjTYu = cTrace(YuadjYdTYdadjTYu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYuadjYuTYuadjYu = cTrace(YuadjYuTYuadjYu) 
 TrYuadjYuTYuadjTYu = cTrace(YuadjYuTYuadjTYu) 
 TrYuadjTYdTYdadjYu = cTrace(YuadjTYdTYdadjYu) 
 TrYuadjTYuTYuadjYu = cTrace(YuadjTYuTYuadjYu) 
 TrYuTpTYdCTYdadjYu = cTrace(YuTpTYdCTYdadjYu) 
 TrYuTpTYuCTYuadjYu = cTrace(YuTpTYuCTYuadjYu) 
 TrYw3adjYb3TYb3adjYw3 = cTrace(Yw3adjYb3TYb3adjYw3) 
 TrYw3adjYb3TYb3adjTYw3 = cTrace(Yw3adjYb3TYb3adjTYw3) 
 TrYw3adjYeTYeadjYw3 = cTrace(Yw3adjYeTYeadjYw3) 
 TrYw3adjYeTYeadjTYw3 = cTrace(Yw3adjYeTYeadjTYw3) 
 TrYw3adjYw3Yw3adjYw3 = cTrace(Yw3adjYw3Yw3adjYw3) 
 TrYw3adjYw3TYw3adjYw3 = cTrace(Yw3adjYw3TYw3adjYw3) 
 TrYw3adjYw3TYw3adjTYw3 = cTrace(Yw3adjYw3TYw3adjTYw3) 
 TrYw3adjTYb3TYb3adjYw3 = cTrace(Yw3adjTYb3TYb3adjYw3) 
 TrYw3adjTYeTYeadjYw3 = cTrace(Yw3adjTYeTYeadjYw3) 
 TrYw3adjTYw3TYw3adjYw3 = cTrace(Yw3adjTYw3TYw3adjYw3) 
 TrYw3TpTYb3CTYb3adjYw3 = cTrace(Yw3TpTYb3CTYb3adjYw3) 
 TrYw3TpTYeCTYeadjYw3 = cTrace(Yw3TpTYeCTYeadjYw3) 
 TrYw3TpTYw3CTYw3adjYw3 = cTrace(Yw3TpTYw3CTYw3adjYw3) 
 TrYx3adjYx3Yx3adjYx3 = cTrace(Yx3adjYx3Yx3adjYx3) 
 TrYx3adjYx3TYx3adjYx3 = cTrace(Yx3adjYx3TYx3adjYx3) 
 TrYx3adjYx3TYx3adjTYx3 = cTrace(Yx3adjYx3TYx3adjTYx3) 
 TrYx3adjTYx3TYx3adjYx3 = cTrace(Yx3adjTYx3TYx3adjYx3) 
 TrYx3CTYdTpTYdadjYx3 = cTrace(Yx3CTYdTpTYdadjYx3) 
 TrYx3TpTYx3CTYx3adjYx3 = cTrace(Yx3TpTYx3CTYx3adjYx3) 
 TradjYdTpYx3CYx3TYd = cTrace(adjYdTpYx3CYx3TYd) 
 TradjYdTpYx3CTYx3TYd = cTrace(adjYdTpYx3CTYx3TYd) 
 TradjYdTpYx3CTYx3TpTYd = cTrace(adjYdTpYx3CTYx3TpTYd) 
 TradjYdTpTYx3CTYx3TpYd = cTrace(adjYdTpTYx3CTYx3TpYd) 
 TradjYx3TYx3CYdTpYd = cTrace(adjYx3TYx3CYdTpYd) 
 TradjYx3TYx3CTYdTpYd = cTrace(adjYx3TYx3CTYdTpYd) 
 Trmd2YdadjYdYdadjYd = cTrace(md2YdadjYdYdadjYd) 
 Trmd2YdadjYdTpYx3CYx3 = cTrace(md2YdadjYdTpYx3CYx3) 
 Trmd2YdadjYuYuadjYd = cTrace(md2YdadjYuYuadjYd) 
 Trmd2adjYx3Yx3adjYx3Yx3 = cTrace(md2adjYx3Yx3adjYx3Yx3) 
 Trmd2TpYx3CYx3YdadjYd = cTrace(md2TpYx3CYx3YdadjYd) 
 Trme2YeadjYb3Yb3adjYe = cTrace(me2YeadjYb3Yb3adjYe) 
 Trme2YeadjYeYeadjYe = cTrace(me2YeadjYeYeadjYe) 
 Trme2YeadjYw3Yw3adjYe = cTrace(me2YeadjYw3Yw3adjYe) 
 TrmHb32Yb3adjYb3Yb3adjYb3 = cTrace(mHb32Yb3adjYb3Yb3adjYb3) 
 TrmHb32Yb3adjYeYeadjYb3 = cTrace(mHb32Yb3adjYeYeadjYb3) 
 TrmHb32Yb3adjYw3Yw3adjYb3 = cTrace(mHb32Yb3adjYw3Yw3adjYb3) 
 TrmHw32Yw3adjYb3Yb3adjYw3 = cTrace(mHw32Yw3adjYb3Yb3adjYw3) 
 TrmHw32Yw3adjYeYeadjYw3 = cTrace(mHw32Yw3adjYeYeadjYw3) 
 TrmHw32Yw3adjYw3Yw3adjYw3 = cTrace(mHw32Yw3adjYw3Yw3adjYw3) 
 TrmHxb32Yx3adjYx3Yx3adjYx3 = cTrace(mHxb32Yx3adjYx3Yx3adjYx3) 
 TrmHxb32Yx3CYdTpYdadjYx3 = cTrace(mHxb32Yx3CYdTpYdadjYx3) 
 TrmHxb32CYx3YdadjYdTpYx3 = cTrace(mHxb32CYx3YdadjYdTpYx3) 
 Trml2adjYb3Yb3adjYb3Yb3 = cTrace(ml2adjYb3Yb3adjYb3Yb3) 
 Trml2adjYb3Yb3adjYeYe = cTrace(ml2adjYb3Yb3adjYeYe) 
 Trml2adjYb3Yb3adjYw3Yw3 = cTrace(ml2adjYb3Yb3adjYw3Yw3) 
 Trml2adjYeYeadjYb3Yb3 = cTrace(ml2adjYeYeadjYb3Yb3) 
 Trml2adjYeYeadjYeYe = cTrace(ml2adjYeYeadjYeYe) 
 Trml2adjYeYeadjYw3Yw3 = cTrace(ml2adjYeYeadjYw3Yw3) 
 Trml2adjYw3Yw3adjYb3Yb3 = cTrace(ml2adjYw3Yw3adjYb3Yb3) 
 Trml2adjYw3Yw3adjYeYe = cTrace(ml2adjYw3Yw3adjYeYe) 
 Trml2adjYw3Yw3adjYw3Yw3 = cTrace(ml2adjYw3Yw3adjYw3Yw3) 
 Trmq2adjYdYdadjYdYd = cTrace(mq2adjYdYdadjYdYd) 
 Trmq2adjYdYdadjYuYu = cTrace(mq2adjYdYdadjYuYu) 
 Trmq2adjYdTpYx3CYx3Yd = cTrace(mq2adjYdTpYx3CYx3Yd) 
 Trmq2adjYuYuadjYdYd = cTrace(mq2adjYuYuadjYdYd) 
 Trmq2adjYuYuadjYuYu = cTrace(mq2adjYuYuadjYuYu) 
 Trmu2YuadjYdYdadjYu = cTrace(mu2YuadjYdYdadjYu) 
 Trmu2YuadjYuYuadjYu = cTrace(mu2YuadjYuYuadjYu) 
 g1p4 =g1**4 
 g2p4 =g2**4 
 g3p4 =g3**4 
End If 
 
 
Tr1(1) = (3*(-1._dp*(mHd2) + mHu2 + Trmd2 + Trme2 + 5._dp*(TrmHx32) - 5._dp*(TrmHxb32)& 
&  - Trml2 + Trmq2 - 2._dp*(Trmu2)))/5._dp

If (TwoLoopRGE) Then 
Tr2(1) = (3._dp*(mHd2) + 3._dp*(mHu2) + 2._dp*(Trmd2) + 6._dp*(Trme2) +               & 
&  25._dp*(TrmHx32) + 25._dp*(TrmHxb32) + 3._dp*(Trml2) + Trmq2 + 8._dp*(Trmu2))/10._dp

Tr3(1) = (3*(-20._dp*(Trmd2adjYx3Yx3) - 20._dp*(Trmd2YdadjYd) - 20._dp*(Trme2YeadjYe) & 
&  + 50._dp*(TrmHxb32Yx3adjYx3) + 3._dp*(Trml2adjYb3Yb3) + 10._dp*(Trml2adjYeYe)         & 
&  + 15._dp*(Trml2adjYw3Yw3) - 15*g2p2*(mHd2 - mHu2 - 5._dp*(TrmHx32) + 5._dp*(TrmHxb32) & 
&  + Trml2 - Trmq2) - 10._dp*(Trmq2adjYdYd) - 10._dp*(Trmq2adjYuYu) + (80*g3p2*(Trmd2 +  & 
&  5._dp*(TrmHx32) - 5._dp*(TrmHxb32) + Trmq2 - 2._dp*(Trmu2)))/3._dp - (g1p2*(9._dp*(mHd2)& 
&  - 9._dp*(mHu2) - 4._dp*(Trmd2) - 36._dp*(Trme2) - 125._dp*(TrmHx32) + 125._dp*(TrmHxb32)& 
&  + 9._dp*(Trml2) - Trmq2 + 32._dp*(Trmu2)))/3._dp + 40._dp*(Trmu2YuadjYu)              & 
&  - 3*mHu2*TrYb3adjYb3 + 30*mHd2*TrYdadjYd + 10*mHd2*TrYeadjYe - 30*mHu2*TrYuadjYu -    & 
&  15*mHu2*TrYw3adjYw3 - 30*mHu2*TrYx3adjYx3))/100._dp

Tr2(2) = (mHd2 + mHu2 + 6._dp*(TrmHw32) + 3._dp*(TrmHx32) + 3._dp*(TrmHxb32)          & 
&  + Trml2 + 3._dp*(Trmq2))/2._dp

Tr2(3) = Trmd2/2._dp + 8._dp*(TrmHg32) + TrmHx32 + TrmHxb32 + Trmq2 + Trmu2/2._dp

End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = (g1p3*(66 + 25*NGHx3 + 25*NGHxb3))/10._dp

 
 
If (TwoLoopRGE) Then 
betag12 = (g1p3*(6*(199._dp*(g1p2) + 135._dp*(g2p2) + 440._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -    & 
&  70._dp*(TrYdadjYd) - 90._dp*(TrYeadjYe) - 130._dp*(TrYuadjYu) - 60._dp*(TrYw3adjYw3) -& 
&  190._dp*(TrYx3adjYx3)) + 125*(5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 +& 
&  NGHxb3)))/150._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = (g2p3*(2 + 4*NGHw3 + 3*NGHx3 +               & 
&  3*NGHxb3))/2._dp

 
 
If (TwoLoopRGE) Then 
betag22 = (g2p3*(2*(27._dp*(g1p2) + 375._dp*(g2p2) + 360._dp*(g3p2) - 9._dp*(TrYb3adjYb3) -     & 
&  90._dp*(TrYdadjYd) - 30._dp*(TrYeadjYe) - 90._dp*(TrYuadjYu) - 140._dp*(TrYw3adjYw3) -& 
&  90._dp*(TrYx3adjYx3)) + 720*g2p2*NGHw3 + 15*(5._dp*(g1p2) +            & 
&  21._dp*(g2p2) + 16._dp*(g3p2))*(NGHx3 + NGHxb3)))/30._dp

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = g3p3*(-3 + 3*NGHg3 + NGHx3 + NGHxb3)

 
 
If (TwoLoopRGE) Then 
betag32 = (g3p3*(3*(11._dp*(g1p2) + 45._dp*(g2p2) + 70._dp*(g3p2) - 20._dp*(TrYdadjYd) -        & 
&  20._dp*(TrYuadjYu) - 20._dp*(TrYx3adjYx3)) + 810*g3p2*NGHg3 +          & 
&  5*(5._dp*(g1p2) + 9._dp*(g2p2) + 34._dp*(g3p2))*(NGHx3 +               & 
&  NGHxb3)))/15._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = ((-26._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)    & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yu)               & 
& /30._dp + YuadjYdYd + 3._dp*(YuadjYuYu)

 
 
If (TwoLoopRGE) Then 
betaYu2 = (2._dp*(g1p2)/5._dp - 3._dp*(TrYdadjYd) - TrYeadjYe)*YuadjYdYd + ((4._dp*(g1p2) +     & 
&  60._dp*(g2p2) - 9._dp*(TrYb3adjYb3) - 90._dp*(TrYuadjYu) - 45._dp*(TrYw3adjYw3) -     & 
&  90._dp*(TrYx3adjYx3))*YuadjYuYu)/10._dp - 2*(YuadjYdTpYx3CYx3Yd + YuadjYdYdadjYdYd +  & 
&  YuadjYdYdadjYuYu + 2._dp*(YuadjYuYuadjYuYu)) + (Yu*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - & 
&  540._dp*(TrYb3adjYeYeadjYb3) - 2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(5486._dp*(g1p4) + & 
&  20*g1p2*(45._dp*(g2p2) + 136._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 -       & 
&  32._dp*(g3p4)) - 5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) -        & 
&  1350._dp*(TrYeadjYw3Yw3adjYe) - 8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  8100._dp*(TrYx3adjYx3Yx3adjYx3)) + 60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) +& 
&  65*g1p2*(NGHx3 + NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu +   & 
&  TrYx3adjYx3) + g3p2*(3*NGHg3 + NGHx3 + NGHxb3)) +& 
&  2700*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/1800._dp

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = 2._dp*(TpYx3CYx3Yd) + (-7._dp*(g1p2)/15._dp - 3._dp*(g2p2) -               & 
&  16._dp*(g3p2)/3._dp + 3._dp*(TrYdadjYd) + TrYeadjYe)*Yd + 3._dp*(YdadjYdYd)           & 
&  + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = (4._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYdadjYd) - 3._dp*(TrYeadjYe))*YdadjYdYd +& 
&  (2*TpYx3CYx3Yd*(-3._dp*(TrYb3adjYb3) + 5*(2._dp*(g1p2) + 6._dp*(g2p2) -               & 
&  6._dp*(TrYuadjYu) - 3._dp*(TrYw3adjYw3) - 6._dp*(TrYx3adjYx3))) + (8._dp*(g1p2) -     & 
&  3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*YdadjYuYu)/10._dp -& 
&  2*(TpYx3CYx3TpYx3CYx3Yd + YdadjYdTpYx3CYx3Yd + 2._dp*(YdadjYdYdadjYdYd) +             & 
&  YdadjYuYuadjYdYd + YdadjYuYuadjYuYu) + (Yd*(287._dp*(g1p4) + 90*g1p2*g2p2 +           & 
&  675._dp*(g2p4) + 80*(g1p2 + 9._dp*(g2p2))*g3p2 - 160._dp*(g3p4) - 27._dp*(TrYb3adjYeYeadjYb3) -& 
&  540._dp*(TrYdadjYdTpYx3CYx3) - 810._dp*(TrYdadjYdYdadjYd) - 270._dp*(TrYdadjYuYuadjYd) -& 
&  270._dp*(TrYeadjYeYeadjYe) - 135._dp*(TrYeadjYw3Yw3adjYe) + 135*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3)) + 3*g1p2*(-12*(TrYdadjYd -          & 
&  3._dp*(TrYeadjYe)) + 35*g1p2*(NGHx3 + NGHxb3)) +        & 
&  480*g3p2*(3._dp*(TrYdadjYd) + g3p2*(3*NGHg3 + NGHx3 +   & 
&  NGHxb3))))/90._dp

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = (-9._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)       & 
& *Ye + (3*(YeadjYb3Yb3 + 5*(2._dp*(YeadjYeYe) + YeadjYw3Yw3)))/10._dp

 
 
If (TwoLoopRGE) Then 
betaYe2 = (-36._dp*(YeadjYb3Yb3adjYb3Yb3) - 18*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) +            & 
&  TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(YeadjYb3Yb3 + 5._dp*(YeadjYw3Yw3)) -             & 
&  600*((-2._dp*(g2p2) + 3._dp*(TrYdadjYd) + TrYeadjYe)*YeadjYeYe - 2*g2p2*YeadjYw3Yw3) -& 
&  5*(24._dp*(YeadjYb3Yb3adjYeYe) + 12._dp*(YeadjYb3Yb3adjYw3Yw3) + 160._dp*(YeadjYeYeadjYeYe) +& 
&  9._dp*(YeadjYw3Yw3adjYb3Yb3) + 60*(2._dp*(YeadjYw3Yw3adjYeYe) + YeadjYw3Yw3adjYw3Yw3)) +& 
&  20*Ye*(-4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 12*g1p2*TrYeadjYe + 3*(45._dp*(g1p4) +   & 
&  6*g1p2*g2p2 + 25._dp*(g2p4) - TrYb3adjYeYeadjYb3 - 20._dp*(TrYdadjYdTpYx3CYx3) -      & 
&  30._dp*(TrYdadjYdYdadjYd) - 10._dp*(TrYdadjYuYuadjYd) - 10._dp*(TrYeadjYeYeadjYe) -   & 
&  5._dp*(TrYeadjYw3Yw3adjYe)) + 45*g1p4*NGHx3 + 45*g1p4*NGHxb3 +& 
&  15*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))/200._dp

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yb3 
!-------------------- 
 
betaYb31  = (3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) -              & 
&  2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*Yb3)/10._dp + 9._dp*(Yb3adjYb3Yb3)     & 
& /10._dp + Yb3adjYeYe + 3._dp*(Yb3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYb32 = (18*(4._dp*(g1p2) + 20._dp*(g2p2) - 3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) -        & 
&  15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*Yb3adjYb3Yb3 - 72._dp*(Yb3adjYb3Yb3adjYb3Yb3) +& 
&  40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yb3adjYeYe -               & 
&  30*(3._dp*(TrYb3adjYb3) + 5*(-8._dp*(g2p2) + 6._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3) +& 
&  6._dp*(TrYx3adjYx3)))*Yb3adjYw3Yw3 - 5*(12._dp*(Yb3adjYb3Yb3adjYw3Yw3) +              & 
&  24._dp*(Yb3adjYeYeadjYb3Yb3) + 80._dp*(Yb3adjYeYeadjYeYe) + 45._dp*(Yb3adjYw3Yw3adjYb3Yb3) +& 
&  60._dp*(Yb3adjYw3Yw3adjYw3Yw3)) + Yb3*(3*(276._dp*(g1p4) + 120*g1p2*g2p2 +            & 
&  500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) - 20._dp*(TrYb3adjYeYeadjYb3) -        & 
&  95._dp*(TrYb3adjYw3Yw3adjYb3) - 400._dp*(TrYdadjYdTpYx3CYx3) - 200._dp*(TrYdadjYuYuadjYd) -& 
&  100._dp*(TrYeadjYw3Yw3adjYe) - 600._dp*(TrYuadjYuYuadjYu) - 250._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  600._dp*(TrYx3adjYx3Yx3adjYx3)) + 20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu +            & 
&  20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 + 15*g1p4*NGHx3 + 15*g1p4*NGHxb3) +& 
&  300*g2p2*(4._dp*(TrYw3adjYw3) + g2p2*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/200._dp

 
DYb3 = oo16pi2*( betaYb31 + oo16pi2 * betaYb32 ) 

 
Else 
DYb3 = oo16pi2* betaYb31 
End If 
 
 
!-------------------- 
! Yw3 
!-------------------- 
 
betaYw31  = (-3._dp*(g1p2)/5._dp - 7._dp*(g2p2) + 3._dp*(TrYb3adjYb3)/10._dp +        & 
&  3._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3)/2._dp + 3._dp*(TrYx3adjYx3))*Yw3 +            & 
&  3._dp*(Yw3adjYb3Yb3)/10._dp + Yw3adjYeYe + 5._dp*(Yw3adjYw3Yw3)/2._dp

 
 
If (TwoLoopRGE) Then 
betaYw32 = (-36._dp*(Yw3adjYb3Yb3adjYb3Yb3) + 40*(6._dp*(g1p2) - 15._dp*(TrYdadjYd) -            & 
&  5._dp*(TrYeadjYe))*Yw3adjYeYe + 40*(3._dp*(g1p2) + 25._dp*(g2p2))*Yw3adjYw3Yw3 -      & 
&  6*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(3._dp*(Yw3adjYb3Yb3) +& 
&  25._dp*(Yw3adjYw3Yw3)) - 5*(24._dp*(Yw3adjYb3Yb3adjYw3Yw3) + 80._dp*(Yw3adjYeYeadjYeYe) +& 
&  40._dp*(Yw3adjYeYeadjYw3Yw3) + 9._dp*(Yw3adjYw3Yw3adjYb3Yb3) + 120._dp*(Yw3adjYw3Yw3adjYw3Yw3)) +& 
&  Yw3*(828._dp*(g1p4) + 360*g1p2*g2p2 + 5500._dp*(g2p4) - 54._dp*(TrYb3adjYb3Yb3adjYb3) -& 
&  60._dp*(TrYb3adjYeYeadjYb3) - 285._dp*(TrYb3adjYw3Yw3adjYb3) - 1200._dp*(TrYdadjYdTpYx3CYx3) -& 
&  600._dp*(TrYdadjYuYuadjYd) - 300._dp*(TrYeadjYw3Yw3adjYe) - 1800._dp*(TrYuadjYuYuadjYu) +& 
&  1200*g2p2*TrYw3adjYw3 - 750._dp*(TrYw3adjYw3Yw3adjYw3) - 1800._dp*(TrYx3adjYx3Yx3adjYx3) +& 
&  20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 +       & 
&  15*g1p4*NGHx3 + 15*g1p4*NGHxb3) + 700*g2p4*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3))))/200._dp

 
DYw3 = oo16pi2*( betaYw31 + oo16pi2 * betaYw32 ) 

 
Else 
DYw3 = oo16pi2* betaYw31 
End If 
 
 
!-------------------- 
! Yx3 
!-------------------- 
 
betaYx31  = ((-38._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3)   & 
&  + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))*Yx3)              & 
& /30._dp + 3._dp*(Yx3adjYx3Yx3) + 2._dp*(Yx3CYdTpYd)

 
 
If (TwoLoopRGE) Then 
betaYx32 = (8._dp*(g1p2)/5._dp + 6._dp*(g2p2) - 9._dp*(TrYb3adjYb3)/10._dp - 9._dp*(TrYuadjYu) - & 
&  9._dp*(TrYw3adjYw3)/2._dp - 9._dp*(TrYx3adjYx3))*Yx3adjYx3Yx3 + (2*(g1p2 +            & 
&  15._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yx3CYdTpYd)/5._dp -           & 
&  2*(2._dp*(Yx3adjYx3Yx3adjYx3Yx3) + Yx3CYdTpYdadjYx3Yx3 + Yx3CYdTpYdCYdTpYd +          & 
&  Yx3CYdTpYuCYuTpYd) + (Yx3*(-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2565._dp*(TrYb3adjYw3Yw3adjYb3) + 2*(8246._dp*(g1p4) + 20*g1p2*(153._dp*(g2p2) +      & 
&  232._dp*(g3p2)) + 50*(135._dp*(g2p4) + 144*g2p2*g3p2 - 32._dp*(g3p4)) -               & 
&  5400._dp*(TrYdadjYdTpYx3CYx3) - 2700._dp*(TrYdadjYuYuadjYd) - 1350._dp*(TrYeadjYw3Yw3adjYe) -& 
&  8100._dp*(TrYuadjYuYuadjYu) - 3375._dp*(TrYw3adjYw3Yw3adjYw3) - 8100._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 95*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/1800._dp

 
DYx3 = oo16pi2*( betaYx31 + oo16pi2 * betaYx32 ) 

 
Else 
DYx3 = oo16pi2* betaYx31 
End If 
 
 
!-------------------- 
! mue 
!-------------------- 
 
betamue1  = ((3._dp*(TrYb3adjYb3) + 30._dp*(TrYdadjYd) + 10._dp*(TrYeadjYe)           & 
&  + 30._dp*(TrYuadjYu) + 15._dp*(TrYw3adjYw3) - 6*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))& 
& *mue)/10._dp

 
 
If (TwoLoopRGE) Then 
betamue2 = 6*g2p4*mue*NGHw3 + (mue*(80*(-((g1p2 - 40._dp*(g3p2))*TrYdadjYd) +     & 
&  3*g1p2*TrYeadjYe + 2*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 5*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3) +& 
&  3*(276._dp*(g1p4) + 120*g1p2*g2p2 + 500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) -  & 
&  40._dp*(TrYb3adjYeYeadjYb3) - 95._dp*(TrYb3adjYw3Yw3adjYb3) - 800._dp*(TrYdadjYdTpYx3CYx3) -& 
&  600._dp*(TrYdadjYdYdadjYd) - 400._dp*(TrYdadjYuYuadjYd) - 200._dp*(TrYeadjYeYeadjYe) -& 
&  200._dp*(TrYeadjYw3Yw3adjYe) - 600._dp*(TrYuadjYuYuadjYu) + 400*g2p2*TrYw3adjYw3 -    & 
&  250._dp*(TrYw3adjYw3Yw3adjYw3) - 600._dp*(TrYx3adjYx3Yx3adjYx3)) + 300*(g1p4 +        & 
&  3._dp*(g2p4))*NGHx3 + 300*(g1p4 + 3._dp*(g2p4))*NGHxb3))/200._dp

 
Dmue = oo16pi2*( betamue1 + oo16pi2 * betamue2 ) 

 
Else 
Dmue = oo16pi2* betamue1 
End If 
 
 
!-------------------- 
! MXM3 
!-------------------- 
 
betaMXM31  = -((5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))*MXM3)/3._dp +            & 
&  MXM3CYx3TpYx3

 
 
If (TwoLoopRGE) Then 
betaMXM32 = (-9*(20*(MXM3CYx3TpYx3CYx3TpYx3 + MXM3CYx3YdadjYdTpYx3) + MXM3CYx3TpYx3*(4._dp*(g1p2) +& 
&  3._dp*(TrYb3adjYb3) + 30._dp*(TrYuadjYu) + 15._dp*(TrYw3adjYw3) + 30._dp*(TrYx3adjYx3))) +& 
&  5*MXM3*(223._dp*(g1p4) + 90*g1p2*g2p2 + 135._dp*(g2p4) + 32*(5._dp*(g1p2) +           & 
&  9._dp*(g2p2))*g3p2 - 32._dp*(g3p4) + 108*g2p4*NGHw3 + 3*(25._dp*(g1p4) +& 
&  27._dp*(g2p4))*(NGHx3 + NGHxb3) + 96*g3p4*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)))/90._dp

 
DMXM3 = oo16pi2*( betaMXM31 + oo16pi2 * betaMXM32 ) 

 
Else 
DMXM3 = oo16pi2* betaMXM31 
End If 
 
 
!-------------------- 
! MWM3 
!-------------------- 
 
betaMWM31  = -8*g2p2*MWM3 + MWM3CYw3TpYw3 + Yw3adjYw3MWM3

 
 
If (TwoLoopRGE) Then 
betaMWM32 = (-3._dp*(MWM3CYw3TpYb3CYb3TpYw3) - 10._dp*(MWM3CYw3TpYeCYeTpYw3) - 15._dp*(MWM3CYw3TpYw3CYw3TpYw3) -& 
&  3._dp*(Yw3adjYb3Yb3adjYw3MWM3) - 10._dp*(Yw3adjYeYeadjYw3MWM3) + (6._dp*(g1p2) -      & 
&  10._dp*(g2p2) - 3._dp*(TrYb3adjYb3) - 30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) -     & 
&  30._dp*(TrYx3adjYx3))*(MWM3CYw3TpYw3 + Yw3adjYw3MWM3) - 15._dp*(Yw3adjYw3Yw3adjYw3MWM3) +& 
&  40*g2p4*MWM3*(10 + 4*NGHw3 + 3*NGHx3 + 3*NGHxb3))/10._dp

 
DMWM3 = oo16pi2*( betaMWM31 + oo16pi2 * betaMWM32 ) 

 
Else 
DMWM3 = oo16pi2* betaMWM31 
End If 
 
 
!-------------------- 
! MGM3 
!-------------------- 
 
betaMGM31  = -12*g3p2*MGM3

 
 
If (TwoLoopRGE) Then 
betaMGM32 = 12*g3p4*MGM3*(3 + 3*NGHg3 + NGHx3 + NGHxb3)

 
DMGM3 = oo16pi2*( betaMGM31 + oo16pi2 * betaMGM32 ) 

 
Else 
DMGM3 = oo16pi2* betaMGM31 
End If 
 
 
!-------------------- 
! MBM3 
!-------------------- 
 
betaMBM31  = (3*(MBM3CYb3TpYb3 + Yb3adjYb3MBM3))/5._dp

 
 
If (TwoLoopRGE) Then 
betaMBM32 = (-3*(3._dp*(MBM3CYb3TpYb3CYb3TpYb3) + 10._dp*(MBM3CYb3TpYeCYeTpYb3) + 15._dp*(MBM3CYb3TpYw3CYw3TpYb3) +& 
&  3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) -  & 
&  5._dp*(TrYx3adjYx3)))*(MBM3CYb3TpYb3 + Yb3adjYb3MBM3) + 3._dp*(Yb3adjYb3Yb3adjYb3MBM3) +& 
&  10._dp*(Yb3adjYeYeadjYb3MBM3) + 15._dp*(Yb3adjYw3Yw3adjYb3MBM3)))/50._dp

 
DMBM3 = oo16pi2*( betaMBM31 + oo16pi2 * betaMBM32 ) 

 
Else 
DMBM3 = oo16pi2* betaMBM31 
End If 
 
 
!-------------------- 
! TYu 
!-------------------- 
 
betaTYu1  = TYuadjYdYd + 5._dp*(TYuadjYuYu) + ((26*g1p2*MassB)/15._dp +               & 
&  (32*g3p2*MassG)/3._dp + 6*g2p2*MassWB + 3._dp*(TradjYb3TYb3)/5._dp + 6._dp*(TradjYuTYu)& 
&  + 3._dp*(TradjYw3TYw3) + 6._dp*(TradjYx3TYx3))*Yu + 2._dp*(YuadjYdTYd) +              & 
&  4._dp*(YuadjYuTYu) + ((-26._dp*(g1p2) - 90._dp*(g2p2) - 160._dp*(g3p2) +              & 
&  9._dp*(TrYb3adjYb3) + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3) + 90._dp*(TrYx3adjYx3))& 
& *TYu)/30._dp

 
 
If (TwoLoopRGE) Then 
betaTYu2 = (2700*(8._dp*(g2p2) - TrYb3adjYb3 - 10._dp*(TrYuadjYu) - 5._dp*(TrYw3adjYw3) -        & 
&  10._dp*(TrYx3adjYx3))*TYuadjYuYu + 360*((2._dp*(g1p2) - 15._dp*(TrYdadjYd) -          & 
&  5._dp*(TrYeadjYe))*TYuadjYdYd + 2*(2._dp*(g1p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*YuadjYdTYd -& 
&  2*(2*g1p2*MassB + 15._dp*(TradjYdTYd) + 5._dp*(TradjYeTYe))*YuadjYdYd +               & 
&  6*(g1p2 + 5._dp*(g2p2) - TrYb3adjYb3 - 10._dp*(TrYuadjYu) - 5._dp*(TrYw3adjYw3) -     & 
&  10._dp*(TrYx3adjYx3))*YuadjYuTYu - (4*g1p2*MassB + 60*g2p2*MassWB + 9._dp*(TradjYb3TYb3) +& 
&  90._dp*(TradjYuTYu) + 45._dp*(TradjYw3TYw3) + 90._dp*(TradjYx3TYx3))*YuadjYuYu) -     & 
&  3600*(TYuadjYdTpYx3CYx3Yd + TYuadjYdYdadjYdYd + 2._dp*(TYuadjYdYdadjYuYu) +           & 
&  3._dp*(TYuadjYuYuadjYuYu) + 2._dp*(YuadjYdTpTYx3CYx3Yd) + 2._dp*(YuadjYdTpYx3CYx3TYd) +& 
&  2._dp*(YuadjYdTYdadjYdYd) + 2._dp*(YuadjYdTYdadjYuYu) + 2._dp*(YuadjYdYdadjYdTYd) +   & 
&  YuadjYdYdadjYuTYu + 4._dp*(YuadjYuTYuadjYuYu) + 3._dp*(YuadjYuYuadjYuTYu)) -          & 
&  4*Yu*(486._dp*(TrYb3adjYb3TYb3adjYb3) + 270._dp*(TrYb3adjYeTYeadjYb3) +               & 
&  1215._dp*(TrYb3adjYw3TYw3adjYb3) + 2700._dp*(TrYdadjYuTYuadjYd) + 270._dp*(TrYeadjYb3TYb3adjYe) +& 
&  1350._dp*(TrYeadjYw3TYw3adjYe) + 2700._dp*(TrYuadjYdTYdadjYu) + 16200._dp*(TrYuadjYuTYuadjYu) +& 
&  1215._dp*(TrYw3adjYb3TYb3adjYw3) + 1350._dp*(TrYw3adjYeTYeadjYw3) + 6075._dp*(TrYw3adjYw3TYw3adjYw3) +& 
&  4*(2743*g1p4*MassB + 5*g1p2*(136*g3p2*(MassB + MassG) + 45*g2p2*(MassB +              & 
&  MassWB)) + 25*(-32*g3p4*MassG + 135*g2p4*MassWB + 72*g2p2*g3p2*(MassG +               & 
&  MassWB)) + 1350._dp*(TradjYdTpYx3CYx3TYd) + 1350._dp*(TradjYx3TYx3CYdTpYd) +          & 
&  4050._dp*(TrYx3adjYx3TYx3adjYx3)) + 28800*g3p4*MassG*NGHg3 +           & 
&  60*(6*(-2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +  & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  5*(13*g1p4*MassB + 32*g3p4*MassG)*NGHx3 + 5*(13*g1p4*MassB +           & 
&  32*g3p4*MassG)*NGHxb3) + 2700*g2p2*(-2._dp*(TradjYw3TYw3) +            & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2430._dp*(TrYb3adjYw3Yw3adjYb3) - 10800._dp*(TrYdadjYdTpYx3CYx3) - 5400._dp*(TrYdadjYuYuadjYd) -& 
&  2700._dp*(TrYeadjYw3Yw3adjYe) - 16200._dp*(TrYuadjYuYuadjYu) - 6075._dp*(TrYw3adjYw3Yw3adjYw3) +& 
&  4*(2743._dp*(g1p4) + 450*g1p2*g2p2 + 3375._dp*(g2p4) + 80*(17._dp*(g1p2) +            & 
&  45._dp*(g2p2))*g3p2 - 800._dp*(g3p4) - 4050._dp*(TrYx3adjYx3Yx3adjYx3)) +             & 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 65*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))*TYu)/1800._dp

 
DTYu = oo16pi2*( betaTYu1 + oo16pi2 * betaTYu2 ) 

 
Else 
DTYu = oo16pi2* betaTYu1 
End If 
 
 
!-------------------- 
! TYd 
!-------------------- 
 
betaTYd1  = 4._dp*(TpTYx3CYx3Yd) + 2._dp*(TpYx3CYx3TYd) + 5._dp*(TYdadjYdYd)          & 
&  + TYdadjYuYu + ((14*g1p2*MassB)/15._dp + (32*g3p2*MassG)/3._dp + 6*g2p2*MassWB +      & 
&  6._dp*(TradjYdTYd) + 2._dp*(TradjYeTYe))*Yd + 4._dp*(YdadjYdTYd) + 2._dp*(YdadjYuTYu) & 
&  - ((7._dp*(g1p2) + 45._dp*(g2p2) + 80._dp*(g3p2) - 45._dp*(TrYdadjYd) -               & 
&  15._dp*(TrYeadjYe))*TYd)/15._dp

 
 
If (TwoLoopRGE) Then 
betaTYd2 = 4*g1p2*TpTYx3CYx3Yd + 12*g2p2*TpTYx3CYx3Yd + 2*g1p2*TpYx3CYx3TYd + 6*g2p2*TpYx3CYx3TYd -& 
&  4*g1p2*MassB*TpYx3CYx3Yd - 12*g2p2*MassWB*TpYx3CYx3Yd - (6*TpYx3CYx3Yd*TradjYb3TYb3)/5._dp -& 
&  12*TpYx3CYx3Yd*TradjYuTYu - 6*TpYx3CYx3Yd*TradjYw3TYw3 - 12*TpYx3CYx3Yd*TradjYx3TYx3 -& 
&  (6*TpTYx3CYx3Yd*TrYb3adjYb3)/5._dp - (3*TpYx3CYx3TYd*TrYb3adjYb3)/5._dp -             & 
&  12*TpTYx3CYx3Yd*TrYuadjYu - 6*TpYx3CYx3TYd*TrYuadjYu - 6*TpTYx3CYx3Yd*TrYw3adjYw3 -   & 
&  3*TpYx3CYx3TYd*TrYw3adjYw3 - 12*TpTYx3CYx3Yd*TrYx3adjYx3 - 6*TpYx3CYx3TYd*TrYx3adjYx3 +& 
&  (6*g1p2*TYdadjYdYd)/5._dp + 12*g2p2*TYdadjYdYd - 15*TrYdadjYd*TYdadjYdYd -            & 
&  5*TrYeadjYe*TYdadjYdYd + (4*g1p2*TYdadjYuYu)/5._dp - (3*TrYb3adjYb3*TYdadjYuYu)/10._dp -& 
&  3*TrYuadjYu*TYdadjYuYu - (3*TrYw3adjYw3*TYdadjYuYu)/2._dp - 3*TrYx3adjYx3*TYdadjYuYu +& 
&  (6*g1p2*YdadjYdTYd)/5._dp + 6*g2p2*YdadjYdTYd - 12*TrYdadjYd*YdadjYdTYd -             & 
&  4*TrYeadjYe*YdadjYdTYd - (2*(4*g1p2*MassB + 30*g2p2*MassWB + 45._dp*(TradjYdTYd) +    & 
&  15._dp*(TradjYeTYe))*YdadjYdYd)/5._dp + (8*g1p2*YdadjYuTYu)/5._dp - (3*TrYb3adjYb3*YdadjYuTYu)/5._dp -& 
&  6*TrYuadjYu*YdadjYuTYu - 3*TrYw3adjYw3*YdadjYuTYu - 6*TrYx3adjYx3*YdadjYuTYu -        & 
&  (8*g1p2*MassB*YdadjYuYu)/5._dp - (3*TradjYb3TYb3*YdadjYuYu)/5._dp - 6*TradjYuTYu*YdadjYuYu -& 
&  3*TradjYw3TYw3*YdadjYuYu - 6*TradjYx3TYx3*YdadjYuYu - 2*(2*(TpTYx3CYx3TpYx3CYx3Yd +   & 
&  TpYx3CYx3TpTYx3CYx3Yd) + TpYx3CYx3TpYx3CYx3TYd + TYdadjYdTpYx3CYx3Yd + 3._dp*(TYdadjYdYdadjYdYd) +& 
&  2._dp*(TYdadjYuYuadjYdYd) + TYdadjYuYuadjYuYu + 2._dp*(YdadjYdTpTYx3CYx3Yd) +         & 
&  2._dp*(YdadjYdTpYx3CYx3TYd) + 4._dp*(YdadjYdTYdadjYdYd) + 3._dp*(YdadjYdYdadjYdTYd) + & 
&  2._dp*(YdadjYuTYuadjYdYd) + 2._dp*(YdadjYuTYuadjYuYu) + YdadjYuYuadjYdTYd +           & 
&  2._dp*(YdadjYuYuadjYuTYu)) - (Yd*(2*(287*g1p4*MassB + 5*g1p2*(8*g3p2*(MassB +         & 
&  MassG) + 9*g2p2*(MassB + MassWB)) + 5*(-32*g3p4*MassG + 135*g2p4*MassWB +             & 
&  72*g2p2*g3p2*(MassG + MassWB)) + 270._dp*(TradjYdTpYx3CYx3TYd) + 270._dp*(TradjYx3TYx3CYdTpYd)) +& 
&  27._dp*(TrYb3adjYeTYeadjYb3) + 1620._dp*(TrYdadjYdTYdadjYd) + 270._dp*(TrYdadjYuTYuadjYd) +& 
&  27._dp*(TrYeadjYb3TYb3adjYe) + 540._dp*(TrYeadjYeTYeadjYe) + 135._dp*(TrYeadjYw3TYw3adjYe) +& 
&  270._dp*(TrYuadjYdTYdadjYu) + 135._dp*(TrYw3adjYeTYeadjYw3) + 1080*g2p4*MassWB*NGHw3 +& 
&  30*(-48*g3p2*TradjYdTYd + 48*g3p2*MassG*TrYdadjYd + 96*g3p4*MassG*NGHg3 +& 
&  32*g3p4*MassG*NGHx3 + 27*g2p4*MassWB*NGHx3 +            & 
&  32*g3p4*MassG*NGHxb3 + 27*g2p4*MassWB*NGHxb3) +         & 
&  6*g1p2*(6*(TradjYdTYd - 3._dp*(TradjYeTYe)) + MassB*(-6*(TrYdadjYd - 3._dp*(TrYeadjYe)) +& 
&  35*g1p2*(NGHx3 + NGHxb3)))))/45._dp + ((287._dp*(g1p4) +& 
&  90*g1p2*g2p2 + 675._dp*(g2p4) + 80*(g1p2 + 9._dp*(g2p2))*g3p2 - 160._dp*(g3p4) -      & 
&  27._dp*(TrYb3adjYeYeadjYb3) - 540._dp*(TrYdadjYdTpYx3CYx3) - 810._dp*(TrYdadjYdYdadjYd) -& 
&  270._dp*(TrYdadjYuYuadjYd) - 270._dp*(TrYeadjYeYeadjYe) - 135._dp*(TrYeadjYw3Yw3adjYe) +& 
&  135*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3)) +& 
&  3*g1p2*(-12*(TrYdadjYd - 3._dp*(TrYeadjYe)) + 35*g1p2*(NGHx3 +         & 
&  NGHxb3)) + 480*g3p2*(3._dp*(TrYdadjYd) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)))*TYd)/90._dp

 
DTYd = oo16pi2*( betaTYd1 + oo16pi2 * betaTYd2 ) 

 
Else 
DTYd = oo16pi2* betaTYd1 
End If 
 
 
!-------------------- 
! TYe 
!-------------------- 
 
betaTYe1  = 3._dp*(TYeadjYb3Yb3)/10._dp + 5._dp*(TYeadjYeYe) + 3._dp*(TYeadjYw3Yw3)   & 
& /2._dp + ((18*g1p2*MassB)/5._dp + 6*g2p2*MassWB + 6._dp*(TradjYdTYd) + 2._dp*(TradjYeTYe))& 
& *Ye + 3._dp*(YeadjYb3TYb3)/5._dp + 4._dp*(YeadjYeTYe) + 3._dp*(YeadjYw3TYw3)           & 
&  + ((-9._dp*(g1p2) - 15._dp*(g2p2) + 15._dp*(TrYdadjYd) + 5._dp*(TrYeadjYe))           & 
& *TYe)/5._dp

 
 
If (TwoLoopRGE) Then 
betaTYe2 = (-36._dp*(TYeadjYb3Yb3adjYb3Yb3) - 240._dp*(TYeadjYb3Yb3adjYeYe) - 45._dp*(TYeadjYb3Yb3adjYw3Yw3) -& 
&  1200._dp*(TYeadjYeYeadjYeYe) - 45._dp*(TYeadjYw3Yw3adjYb3Yb3) - 1200._dp*(TYeadjYw3Yw3adjYeYe) -& 
&  225._dp*(TYeadjYw3Yw3adjYw3Yw3) - 72._dp*(YeadjYb3TYb3adjYb3Yb3) - 240._dp*(YeadjYb3TYb3adjYeYe) -& 
&  90._dp*(YeadjYb3TYb3adjYw3Yw3) - 72._dp*(YeadjYb3Yb3adjYb3TYb3) - 120._dp*(YeadjYb3Yb3adjYeTYe) -& 
&  90._dp*(YeadjYb3Yb3adjYw3TYw3) - 1600._dp*(YeadjYeTYeadjYeYe) - 1200*(2*g2p2*MassWB + & 
&  3._dp*(TradjYdTYd) + TradjYeTYe)*YeadjYeYe - 1200._dp*(YeadjYeYeadjYeTYe) -           & 
&  18*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))*(TYeadjYb3Yb3 +& 
&  5._dp*(TYeadjYw3Yw3) + 2._dp*(YeadjYb3TYb3) + 10._dp*(YeadjYw3TYw3)) - 90._dp*(YeadjYw3TYw3adjYb3Yb3) -& 
&  1200._dp*(YeadjYw3TYw3adjYeYe) - 450._dp*(YeadjYw3TYw3adjYw3Yw3) - 36*(TradjYb3TYb3 + & 
&  5*(2._dp*(TradjYuTYu) + TradjYw3TYw3 + 2._dp*(TradjYx3TYx3)))*(YeadjYb3Yb3 +          & 
&  5._dp*(YeadjYw3Yw3)) + 40*(-6*g1p2*TYeadjYeYe + 60*g2p2*TYeadjYeYe - 75*TrYdadjYd*TYeadjYeYe -& 
&  25*TrYeadjYe*TYeadjYeYe + 30*g2p2*TYeadjYw3Yw3 + (6._dp*(g1p2) + 30._dp*(g2p2) -      & 
&  60._dp*(TrYdadjYd) - 20._dp*(TrYeadjYe))*YeadjYeTYe + 60*g2p2*YeadjYw3TYw3 -          & 
&  60*g2p2*MassWB*YeadjYw3Yw3) - 90._dp*(YeadjYw3Yw3adjYb3TYb3) - 600._dp*(YeadjYw3Yw3adjYeTYe) -& 
&  450._dp*(YeadjYw3Yw3adjYw3TYw3) - 40*Ye*(-4*(-((g1p2 - 40._dp*(g3p2))*TradjYdTYd) +   & 
&  3*g1p2*TradjYeTYe + (g1p2*MassB - 40*g3p2*MassG)*TrYdadjYd - 3*g1p2*MassB*TrYeadjYe) +& 
&  3*(90*g1p4*MassB + 50*g2p4*MassWB + 6*g1p2*g2p2*(MassB + MassWB) + 20._dp*(TradjYdTpYx3CYx3TYd) +& 
&  20._dp*(TradjYx3TYx3CYdTpYd) + TrYb3adjYeTYeadjYb3 + 60._dp*(TrYdadjYdTYdadjYd) +     & 
&  10._dp*(TrYdadjYuTYuadjYd) + TrYeadjYb3TYb3adjYe + 20._dp*(TrYeadjYeTYeadjYe) +       & 
&  5._dp*(TrYeadjYw3TYw3adjYe) + 10._dp*(TrYuadjYdTYdadjYu) + 5._dp*(TrYw3adjYeTYeadjYw3)) +& 
&  90*g1p4*MassB*NGHx3 + 90*g1p4*MassB*NGHxb3 +            & 
&  30*g2p4*MassWB*(4*NGHw3 + 3*(NGHx3 + NGHxb3))) +& 
&  20*(-4*(g1p2 - 40._dp*(g3p2))*TrYdadjYd + 12*g1p2*TrYeadjYe + 3*(45._dp*(g1p4) +      & 
&  6*g1p2*g2p2 + 25._dp*(g2p4) - TrYb3adjYeYeadjYb3 - 20._dp*(TrYdadjYdTpYx3CYx3) -      & 
&  30._dp*(TrYdadjYdYdadjYd) - 10._dp*(TrYdadjYuYuadjYd) - 10._dp*(TrYeadjYeYeadjYe) -   & 
&  5._dp*(TrYeadjYw3Yw3adjYe)) + 45*g1p4*NGHx3 + 45*g1p4*NGHxb3 +& 
&  15*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))*TYe)/200._dp

 
DTYe = oo16pi2*( betaTYe1 + oo16pi2 * betaTYe2 ) 

 
Else 
DTYe = oo16pi2* betaTYe1 
End If 
 
 
!-------------------- 
! TYb3 
!-------------------- 
 
betaTYb31  = (12*(2*g1p2*MassB + 10*g2p2*MassWB + TradjYb3TYb3 + 10._dp*(TradjYuTYu)  & 
&  + 5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3))*Yb3 + 24._dp*(Yb3adjYb3TYb3)          & 
&  + 5*(6._dp*(TYb3adjYb3Yb3) + 4._dp*(TYb3adjYeYe) + 12._dp*(TYb3adjYw3Yw3)             & 
&  + 8._dp*(Yb3adjYeTYe) + 15._dp*(Yb3adjYw3TYw3)) + 6*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) & 
&  + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*TYb3)/20._dp

 
 
If (TwoLoopRGE) Then 
betaTYb32 = (144*g1p2*TYb3adjYb3Yb3 + 720*g2p2*TYb3adjYb3Yb3 - 90*TrYb3adjYb3*TYb3adjYb3Yb3 -     & 
&  900*TrYuadjYu*TYb3adjYb3Yb3 - 450*TrYw3adjYw3*TYb3adjYb3Yb3 - 900*TrYx3adjYx3*TYb3adjYb3Yb3 -& 
&  108._dp*(TYb3adjYb3Yb3adjYb3Yb3) - 135._dp*(TYb3adjYb3Yb3adjYw3Yw3) + 240*g1p2*TYb3adjYeYe -& 
&  600*TrYdadjYd*TYb3adjYeYe - 200*TrYeadjYe*TYb3adjYeYe - 240._dp*(TYb3adjYeYeadjYb3Yb3) -& 
&  400._dp*(TYb3adjYeYeadjYeYe) - 300._dp*(TYb3adjYeYeadjYw3Yw3) + 180*g1p2*TYb3adjYw3Yw3 +& 
&  2100*g2p2*TYb3adjYw3Yw3 - 180*TrYb3adjYb3*TYb3adjYw3Yw3 - 1800*TrYuadjYu*TYb3adjYw3Yw3 -& 
&  900*TrYw3adjYw3*TYb3adjYw3Yw3 - 1800*TrYx3adjYx3*TYb3adjYw3Yw3 - 405._dp*(TYb3adjYw3Yw3adjYb3Yb3) -& 
&  675._dp*(TYb3adjYw3Yw3adjYw3Yw3) + 72*g1p2*Yb3adjYb3TYb3 + 360*g2p2*Yb3adjYb3TYb3 -   & 
&  72*TrYb3adjYb3*Yb3adjYb3TYb3 - 720*TrYuadjYu*Yb3adjYb3TYb3 - 360*TrYw3adjYw3*Yb3adjYb3TYb3 -& 
&  720*TrYx3adjYx3*Yb3adjYb3TYb3 - 144._dp*(Yb3adjYb3TYb3adjYb3Yb3) - 180._dp*(Yb3adjYb3TYb3adjYw3Yw3) -& 
&  36*(4*g1p2*MassB + 20*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYuTYu) +      & 
&  15._dp*(TradjYw3TYw3) + 30._dp*(TradjYx3TYx3))*Yb3adjYb3Yb3 - 108._dp*(Yb3adjYb3Yb3adjYb3TYb3) -& 
&  135._dp*(Yb3adjYb3Yb3adjYw3TYw3) + 480*g1p2*Yb3adjYeTYe - 1200*TrYdadjYd*Yb3adjYeTYe -& 
&  400*TrYeadjYe*Yb3adjYeTYe - 240._dp*(Yb3adjYeTYeadjYb3Yb3) - 800._dp*(Yb3adjYeTYeadjYeYe) -& 
&  300._dp*(Yb3adjYeTYeadjYw3Yw3) - 480*g1p2*MassB*Yb3adjYeYe - 1200*TradjYdTYd*Yb3adjYeYe -& 
&  400*TradjYeTYe*Yb3adjYeYe - 120._dp*(Yb3adjYeYeadjYb3TYb3) - 800._dp*(Yb3adjYeYeadjYeTYe) -& 
&  150._dp*(Yb3adjYeYeadjYw3TYw3) + 90*g1p2*Yb3adjYw3TYw3 + 2850*g2p2*Yb3adjYw3TYw3 -    & 
&  225*TrYb3adjYb3*Yb3adjYw3TYw3 - 2250*TrYuadjYu*Yb3adjYw3TYw3 - 1125*TrYw3adjYw3*Yb3adjYw3TYw3 -& 
&  2250*TrYx3adjYx3*Yb3adjYw3TYw3 - 450._dp*(Yb3adjYw3TYw3adjYb3Yb3) - 900._dp*(Yb3adjYw3TYw3adjYw3Yw3) -& 
&  180*g1p2*MassB*Yb3adjYw3Yw3 - 3300*g2p2*MassWB*Yb3adjYw3Yw3 - 270*TradjYb3TYb3*Yb3adjYw3Yw3 -& 
&  2700*TradjYuTYu*Yb3adjYw3Yw3 - 1350*TradjYw3TYw3*Yb3adjYw3Yw3 - 2700*TradjYx3TYx3*Yb3adjYw3Yw3 -& 
&  270._dp*(Yb3adjYw3Yw3adjYb3TYb3) - 675._dp*(Yb3adjYw3Yw3adjYw3TYw3) - 4*Yb3*(40*(-    & 
& 2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +           & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  3*(18._dp*(TrYb3adjYb3TYb3adjYb3) + 10._dp*(TrYb3adjYeTYeadjYb3) + 45._dp*(TrYb3adjYw3TYw3adjYb3) +& 
&  100._dp*(TrYdadjYuTYuadjYd) + 10._dp*(TrYeadjYb3TYb3adjYe) + 50._dp*(TrYeadjYw3TYw3adjYe) +& 
&  100._dp*(TrYuadjYdTYdadjYu) + 600._dp*(TrYuadjYuTYuadjYu) + 45._dp*(TrYw3adjYb3TYb3adjYw3) +& 
&  50._dp*(TrYw3adjYeTYeadjYw3) + 225._dp*(TrYw3adjYw3TYw3adjYw3) + 4*(69*g1p4*MassB +   & 
&  125*g2p4*MassWB + 15*g1p2*g2p2*(MassB + MassWB) + 50._dp*(TradjYdTpYx3CYx3TYd) +      & 
&  50._dp*(TradjYx3TYx3CYdTpYd) + 150._dp*(TrYx3adjYx3TYx3adjYx3))) + 300*g1p4*MassB*NGHx3 +& 
&  300*g1p4*MassB*NGHxb3 + 300*g2p2*(-2._dp*(TradjYw3TYw3) +              & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (3*(276._dp*(g1p4) + 120*g1p2*g2p2 + 500._dp*(g2p4) -     & 
&  18._dp*(TrYb3adjYb3Yb3adjYb3) - 20._dp*(TrYb3adjYeYeadjYb3) - 90._dp*(TrYb3adjYw3Yw3adjYb3) -& 
&  400._dp*(TrYdadjYdTpYx3CYx3) - 200._dp*(TrYdadjYuYuadjYd) - 100._dp*(TrYeadjYw3Yw3adjYe) -& 
&  600._dp*(TrYuadjYuYuadjYu) - 225._dp*(TrYw3adjYw3Yw3adjYw3) - 600._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu + 20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 +       & 
&  15*g1p4*NGHx3 + 15*g1p4*NGHxb3) + 300*g2p2*(4._dp*(TrYw3adjYw3) +& 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))*TYb3)/200._dp

 
DTYb3 = oo16pi2*( betaTYb31 + oo16pi2 * betaTYb32 ) 

 
Else 
DTYb3 = oo16pi2* betaTYb31 
End If 
 
 
!-------------------- 
! TYw3 
!-------------------- 
 
betaTYw31  = 9._dp*(TYw3adjYb3Yb3)/10._dp + TYw3adjYeYe + 3._dp*(TYw3adjYw3Yw3)       & 
&  + ((6*g1p2*MassB)/5._dp + 14*g2p2*MassWB + 3._dp*(TradjYb3TYb3)/5._dp +               & 
&  6._dp*(TradjYuTYu) + 3._dp*(TradjYw3TYw3) + 6._dp*(TradjYx3TYx3))*Yw3 +               & 
&  9._dp*(Yw3adjYb3TYb3)/10._dp + 2._dp*(Yw3adjYeTYe) + 15._dp*(Yw3adjYw3TYw3)           & 
& /4._dp + ((-6._dp*(g1p2) - 70._dp*(g2p2) + 3._dp*(TrYb3adjYb3) + 30._dp*(TrYuadjYu)    & 
&  + 15._dp*(TrYw3adjYw3) + 30._dp*(TrYx3adjYx3))*TYw3)/10._dp

 
 
If (TwoLoopRGE) Then 
betaTYw32 = (72*g1p2*TYw3adjYb3Yb3 - 120*g2p2*TYw3adjYb3Yb3 - 54*TrYb3adjYb3*TYw3adjYb3Yb3 -      & 
&  540*TrYuadjYu*TYw3adjYb3Yb3 - 270*TrYw3adjYw3*TYw3adjYb3Yb3 - 540*TrYx3adjYx3*TYw3adjYb3Yb3 -& 
&  72._dp*(TYw3adjYb3Yb3adjYb3Yb3) - 135._dp*(TYw3adjYb3Yb3adjYw3Yw3) + 240*g1p2*TYw3adjYeYe -& 
&  600*TrYdadjYd*TYw3adjYeYe - 200*TrYeadjYe*TYw3adjYeYe - 120._dp*(TYw3adjYeYeadjYb3Yb3) -& 
&  400._dp*(TYw3adjYeYeadjYeYe) - 300._dp*(TYw3adjYeYeadjYw3Yw3) + 180*g1p2*TYw3adjYw3Yw3 +& 
&  900*g2p2*TYw3adjYw3Yw3 - 180*TrYb3adjYb3*TYw3adjYw3Yw3 - 1800*TrYuadjYu*TYw3adjYw3Yw3 -& 
&  900*TrYw3adjYw3*TYw3adjYw3Yw3 - 1800*TrYx3adjYx3*TYw3adjYw3Yw3 - 225._dp*(TYw3adjYw3Yw3adjYb3Yb3) -& 
&  675._dp*(TYw3adjYw3Yw3adjYw3Yw3) + 36*g1p2*Yw3adjYb3TYb3 - 60*g2p2*Yw3adjYb3TYb3 -    & 
&  54*TrYb3adjYb3*Yw3adjYb3TYb3 - 540*TrYuadjYu*Yw3adjYb3TYb3 - 270*TrYw3adjYw3*Yw3adjYb3TYb3 -& 
&  540*TrYx3adjYx3*Yw3adjYb3TYb3 - 108._dp*(Yw3adjYb3TYb3adjYb3Yb3) - 180._dp*(Yw3adjYb3TYb3adjYw3Yw3) -& 
&  24*(3*g1p2*MassB - 5*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYuTYu) +       & 
&  15._dp*(TradjYw3TYw3) + 30._dp*(TradjYx3TYx3))*Yw3adjYb3Yb3 - 90._dp*(Yw3adjYb3Yb3adjYb3TYb3) -& 
&  135._dp*(Yw3adjYb3Yb3adjYw3TYw3) + 480*g1p2*Yw3adjYeTYe - 1200*TrYdadjYd*Yw3adjYeTYe -& 
&  400*TrYeadjYe*Yw3adjYeTYe - 120._dp*(Yw3adjYeTYeadjYb3Yb3) - 800._dp*(Yw3adjYeTYeadjYeYe) -& 
&  300._dp*(Yw3adjYeTYeadjYw3Yw3) - 480*g1p2*MassB*Yw3adjYeYe - 1200*TradjYdTYd*Yw3adjYeYe -& 
&  400*TradjYeTYe*Yw3adjYeYe - 60._dp*(Yw3adjYeYeadjYb3TYb3) - 800._dp*(Yw3adjYeYeadjYeTYe) -& 
&  150._dp*(Yw3adjYeYeadjYw3TYw3) + 90*g1p2*Yw3adjYw3TYw3 + 2250*g2p2*Yw3adjYw3TYw3 -    & 
&  225*TrYb3adjYb3*Yw3adjYw3TYw3 - 2250*TrYuadjYu*Yw3adjYw3TYw3 - 1125*TrYw3adjYw3*Yw3adjYw3TYw3 -& 
&  2250*TrYx3adjYx3*Yw3adjYw3TYw3 - 270._dp*(Yw3adjYw3TYw3adjYb3Yb3) - 900._dp*(Yw3adjYw3TYw3adjYw3Yw3) -& 
&  180*g1p2*MassB*Yw3adjYw3Yw3 - 2100*g2p2*MassWB*Yw3adjYw3Yw3 - 270*TradjYb3TYb3*Yw3adjYw3Yw3 -& 
&  2700*TradjYuTYu*Yw3adjYw3Yw3 - 1350*TradjYw3TYw3*Yw3adjYw3Yw3 - 2700*TradjYx3TYx3*Yw3adjYw3Yw3 -& 
&  180._dp*(Yw3adjYw3Yw3adjYb3TYb3) - 675._dp*(Yw3adjYw3Yw3adjYw3TYw3) - 4*Yw3*(54._dp*(TrYb3adjYb3TYb3adjYb3) +& 
&  30._dp*(TrYb3adjYeTYeadjYb3) + 135._dp*(TrYb3adjYw3TYw3adjYb3) + 300._dp*(TrYdadjYuTYuadjYd) +& 
&  30._dp*(TrYeadjYb3TYb3adjYe) + 150._dp*(TrYeadjYw3TYw3adjYe) + 300._dp*(TrYuadjYdTYdadjYu) +& 
&  1800._dp*(TrYuadjYuTYuadjYu) + 135._dp*(TrYw3adjYb3TYb3adjYw3) + 150._dp*(TrYw3adjYeTYeadjYw3) +& 
&  675._dp*(TrYw3adjYw3TYw3adjYw3) + 4*(207*g1p4*MassB + 1375*g2p4*MassWB +              & 
&  45*g1p2*g2p2*(MassB + MassWB) + 150._dp*(TradjYdTpYx3CYx3TYd) + 150._dp*(TradjYx3TYx3CYdTpYd) +& 
&  450._dp*(TrYx3adjYx3TYx3adjYx3)) + 300*g1p4*MassB*NGHx3 +              & 
&  20*(2*(-2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +  & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  15*g1p4*MassB*NGHxb3) + 100*g2p2*(-6._dp*(TradjYw3TYw3) +              & 
&  6*MassWB*TrYw3adjYw3 + 7*g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (828._dp*(g1p4) + 360*g1p2*g2p2 + 5500._dp*(g2p4) -       & 
&  54._dp*(TrYb3adjYb3Yb3adjYb3) - 60._dp*(TrYb3adjYeYeadjYb3) - 270._dp*(TrYb3adjYw3Yw3adjYb3) -& 
&  1200._dp*(TrYdadjYdTpYx3CYx3) - 600._dp*(TrYdadjYuYuadjYd) - 300._dp*(TrYeadjYw3Yw3adjYe) -& 
&  1800._dp*(TrYuadjYuYuadjYu) + 1200*g2p2*TrYw3adjYw3 - 675._dp*(TrYw3adjYw3Yw3adjYw3) -& 
&  1800._dp*(TrYx3adjYx3Yx3adjYx3) + 20*(8*(g1p2 + 20._dp*(g3p2))*TrYuadjYu +            & 
&  20*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3 + 15*g1p4*NGHx3 + 15*g1p4*NGHxb3) +& 
&  700*g2p4*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))*TYw3)/200._dp

 
DTYw3 = oo16pi2*( betaTYw31 + oo16pi2 * betaTYw32 ) 

 
Else 
DTYw3 = oo16pi2* betaTYw31 
End If 
 
 
!-------------------- 
! TYx3 
!-------------------- 
 
betaTYx31  = 4._dp*(TYx3adjYx3Yx3) + 2._dp*(TYx3CYdTpYd) + ((38*g1p2*MassB)           & 
& /15._dp + (32*g3p2*MassG)/3._dp + 6*g2p2*MassWB + 3._dp*(TradjYb3TYb3)/5._dp +         & 
&  6._dp*(TradjYuTYu) + 3._dp*(TradjYw3TYw3) + 6._dp*(TradjYx3TYx3))*Yx3 +               & 
&  5._dp*(Yx3adjYx3TYx3) + 4._dp*(Yx3CYdTpTYd) + ((-38._dp*(g1p2) - 90._dp*(g2p2)        & 
&  - 160._dp*(g3p2) + 9._dp*(TrYb3adjYb3) + 90._dp*(TrYuadjYu) + 45._dp*(TrYw3adjYw3)    & 
&  + 90._dp*(TrYx3adjYx3))*TYx3)/30._dp

 
 
If (TwoLoopRGE) Then 
betaTYx32 = (180*(12*(g1p2 + 5._dp*(g2p2) - TrYb3adjYb3 - 10._dp*(TrYuadjYu) - 5._dp*(TrYw3adjYw3) -& 
&  10._dp*(TrYx3adjYx3))*TYx3adjYx3Yx3 + 4*(g1p2 + 15._dp*(g2p2) - 15._dp*(TrYdadjYd) -  & 
&  5._dp*(TrYeadjYe))*TYx3CYdTpYd + 3*(12._dp*(g1p2) + 40._dp*(g2p2) - 5._dp*(TrYb3adjYb3) -& 
&  50._dp*(TrYuadjYu) - 25._dp*(TrYw3adjYw3) - 50._dp*(TrYx3adjYx3))*Yx3adjYx3TYx3 -     & 
&  2*(16*g1p2*MassB + 60*g2p2*MassWB + 9._dp*(TradjYb3TYb3) + 90._dp*(TradjYuTYu) +      & 
&  45._dp*(TradjYw3TYw3) + 90._dp*(TradjYx3TYx3))*Yx3adjYx3Yx3 + 8*(g1p2 +               & 
&  15._dp*(g2p2) - 15._dp*(TrYdadjYd) - 5._dp*(TrYeadjYe))*Yx3CYdTpTYd - 8*(g1p2*MassB + & 
&  15*g2p2*MassWB + 15._dp*(TradjYdTYd) + 5._dp*(TradjYeTYe))*Yx3CYdTpYd) -              & 
&  3600*(3._dp*(TYx3adjYx3Yx3adjYx3Yx3) + 2._dp*(TYx3CYdTpYdadjYx3Yx3) + TYx3CYdTpYdCYdTpYd +& 
&  TYx3CYdTpYuCYuTpYd + 4._dp*(Yx3adjYx3TYx3adjYx3Yx3) + 3._dp*(Yx3adjYx3Yx3adjYx3TYx3) +& 
&  2._dp*(Yx3CYdTpTYdCYdTpYd) + 2._dp*(Yx3CYdTpTYuCYuTpYd) + Yx3CYdTpYdadjYx3TYx3 +      & 
&  2._dp*(Yx3CYdTpYdCYdTpTYd) + 2._dp*(Yx3CYdTpYuCYuTpTYd) + 2._dp*(Yx3TYdadjYdadjYx3Yx3)) -& 
&  4*Yx3*(486._dp*(TrYb3adjYb3TYb3adjYb3) + 270._dp*(TrYb3adjYeTYeadjYb3) +              & 
&  1215._dp*(TrYb3adjYw3TYw3adjYb3) + 2700._dp*(TrYdadjYuTYuadjYd) + 270._dp*(TrYeadjYb3TYb3adjYe) +& 
&  1350._dp*(TrYeadjYw3TYw3adjYe) + 2700._dp*(TrYuadjYdTYdadjYu) + 16200._dp*(TrYuadjYuTYuadjYu) +& 
&  1215._dp*(TrYw3adjYb3TYb3adjYw3) + 1350._dp*(TrYw3adjYeTYeadjYw3) + 6075._dp*(TrYw3adjYw3TYw3adjYw3) +& 
&  4*(4123*g1p4*MassB + 5*g1p2*(232*g3p2*(MassB + MassG) + 153*g2p2*(MassB +             & 
&  MassWB)) + 25*(-32*g3p4*MassG + 135*g2p4*MassWB + 72*g2p2*g3p2*(MassG +               & 
&  MassWB)) + 1350._dp*(TradjYdTpYx3CYx3TYd) + 1350._dp*(TradjYx3TYx3CYdTpYd) +          & 
&  4050._dp*(TrYx3adjYx3TYx3adjYx3)) + 28800*g3p4*MassG*NGHg3 +           & 
&  60*(6*(-2*(g1p2 + 20._dp*(g3p2))*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 +  & 
&  2*(g1p2*MassB + 20*g3p2*MassG)*TrYuadjYu + 5*(g1p2*MassB + 8*g3p2*MassG)*TrYx3adjYx3) +& 
&  5*(19*g1p4*MassB + 32*g3p4*MassG)*NGHx3 + 5*(19*g1p4*MassB +           & 
&  32*g3p4*MassG)*NGHxb3) + 2700*g2p2*(-2._dp*(TradjYw3TYw3) +            & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))) + (-486._dp*(TrYb3adjYb3Yb3adjYb3) - 540._dp*(TrYb3adjYeYeadjYb3) -& 
&  2430._dp*(TrYb3adjYw3Yw3adjYb3) - 10800._dp*(TrYdadjYdTpYx3CYx3) - 5400._dp*(TrYdadjYuYuadjYd) -& 
&  2700._dp*(TrYeadjYw3Yw3adjYe) - 16200._dp*(TrYuadjYuYuadjYu) - 6075._dp*(TrYw3adjYw3Yw3adjYw3) +& 
&  4*(4123._dp*(g1p4) + 1530*g1p2*g2p2 + 3375._dp*(g2p4) + 80*(29._dp*(g1p2) +           & 
&  45._dp*(g2p2))*g3p2 - 800._dp*(g3p4) - 4050._dp*(TrYx3adjYx3Yx3adjYx3)) +             & 
&  60*g1p2*(24._dp*(TrYuadjYu) + 60._dp*(TrYx3adjYx3) + 95*g1p2*(NGHx3 +  & 
&  NGHxb3)) + 9600*g3p2*(3*(TrYuadjYu + TrYx3adjYx3) + g3p2*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 2700*g2p2*(4._dp*(TrYw3adjYw3) +     & 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3))))*TYx3)/1800._dp

 
DTYx3 = oo16pi2*( betaTYx31 + oo16pi2 * betaTYx32 ) 

 
Else 
DTYx3 = oo16pi2* betaTYx31 
End If 
 
 
!-------------------- 
! Bmue 
!-------------------- 
 
betaBmue1  = ((6*g1p2*MassB + 30*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYdTYd)& 
&  + 10._dp*(TradjYeTYe) + 30._dp*(TradjYuTYu) + 15._dp*(TradjYw3TYw3) + 30._dp*(TradjYx3TYx3))& 
& *mue)/5._dp + (-3._dp*(g1p2)/5._dp - 3._dp*(g2p2) + 3._dp*(TrYb3adjYb3)/10._dp +       & 
&  3._dp*(TrYdadjYd) + TrYeadjYe + 3._dp*(TrYuadjYu) + 3._dp*(TrYw3adjYw3)               & 
& /2._dp + 3._dp*(TrYx3adjYx3))*Bmue

 
 
If (TwoLoopRGE) Then 
betaBmue2 = (-48*(69*g1p4*MassB + 125*g2p4*MassWB + 15*g1p2*g2p2*(MassB + MassWB))*mue +          & 
&  Bmue*(80*(-((g1p2 - 40._dp*(g3p2))*TrYdadjYd) + 3*g1p2*TrYeadjYe + 2*(g1p2 +          & 
&  20._dp*(g3p2))*TrYuadjYu + 5*(g1p2 + 8._dp*(g3p2))*TrYx3adjYx3) + 3*(276._dp*(g1p4) + & 
&  120*g1p2*g2p2 + 500._dp*(g2p4) - 18._dp*(TrYb3adjYb3Yb3adjYb3) - 40._dp*(TrYb3adjYeYeadjYb3) -& 
&  90._dp*(TrYb3adjYw3Yw3adjYb3) - 800._dp*(TrYdadjYdTpYx3CYx3) - 600._dp*(TrYdadjYdYdadjYd) -& 
&  400._dp*(TrYdadjYuYuadjYd) - 200._dp*(TrYeadjYeYeadjYe) - 200._dp*(TrYeadjYw3Yw3adjYe) -& 
&  600._dp*(TrYuadjYuYuadjYu) - 225._dp*(TrYw3adjYw3Yw3adjYw3) - 600._dp*(TrYx3adjYx3Yx3adjYx3)) +& 
&  300*g1p4*NGHx3 + 300*g1p4*NGHxb3 + 300*g2p2*(4._dp*(TrYw3adjYw3) +& 
&  g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))) -  & 
&  4*mue*(40*(g1p2*TradjYdTYd - 40*g3p2*TradjYdTYd - 3*g1p2*TradjYeTYe - 2*g1p2*TradjYuTYu -& 
&  40*g3p2*TradjYuTYu - 5*(g1p2 + 8._dp*(g3p2))*TradjYx3TYx3 + (-(g1p2*MassB) +          & 
&  40*g3p2*MassG)*TrYdadjYd + 3*g1p2*MassB*TrYeadjYe + 2*g1p2*MassB*TrYuadjYu +          & 
&  40*g3p2*MassG*TrYuadjYu + 5*g1p2*MassB*TrYx3adjYx3 + 40*g3p2*MassG*TrYx3adjYx3) +     & 
&  3*(18._dp*(TrYb3adjYb3TYb3adjYb3) + 5*(4._dp*(TrYb3adjYeTYeadjYb3) + 9._dp*(TrYb3adjYw3TYw3adjYb3) +& 
&  120._dp*(TrYdadjYdTYdadjYd) + 40._dp*(TrYdadjYuTYuadjYd) + 4._dp*(TrYeadjYb3TYb3adjYe) +& 
&  40._dp*(TrYeadjYeTYeadjYe) + 20._dp*(TrYeadjYw3TYw3adjYe) + 40._dp*(TrYuadjYdTYdadjYu) +& 
&  120._dp*(TrYuadjYuTYuadjYu) + 9._dp*(TrYw3adjYb3TYb3adjYw3) + 5*(4._dp*(TrYw3adjYeTYeadjYw3) +& 
&  9._dp*(TrYw3adjYw3TYw3adjYw3) + 8*(2*(TradjYdTpYx3CYx3TYd + TradjYx3TYx3CYdTpYd) +    & 
&  3._dp*(TrYx3adjYx3TYx3adjYx3))))) + 300*g1p4*MassB*NGHx3 +             & 
&  300*g1p4*MassB*NGHxb3 + 300*g2p2*(-2._dp*(TradjYw3TYw3) +              & 
&  2*MassWB*TrYw3adjYw3 + g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 +& 
&  NGHxb3)))))/200._dp

 
DBmue = oo16pi2*( betaBmue1 + oo16pi2 * betaBmue2 ) 

 
Else 
DBmue = oo16pi2* betaBmue1 
End If 
 
 
!-------------------- 
! BMXM3 
!-------------------- 
 
!betaBMXM31  = (2*(5*g1p2*MassB + 16*g3p2*MassG + 9*g2p2*MassWB)*MXM3 + 3*(BMXM3CYx3TpYx3 +& 
!&  2._dp*(MXM3CYx3TpTYx3)) - (5._dp*(g1p2) + 9._dp*(g2p2) + 16._dp*(g3p2))               & 
!& *BMXM3)/3._dp

 
 
! If (TwoLoopRGE) Then 
! betaBMXM32 = -2*(BMXM3CYx3TpYx3CYx3TpYx3 + BMXM3CYx3YdadjYdTpYx3 + 2._dp*(MXM3CYx3TpTYx3CYx3TpYx3) +& 
! &  2._dp*(MXM3CYx3TpYx3CYx3TpTYx3) + 2._dp*(MXM3CYx3TYdadjYdTpYx3) + 2._dp*(MXM3CYx3YdadjYdTpTYx3)) -& 
! &  (2*MXM3*(223*g1p4*MassB - 32*g3p4*MassG + 135*g2p4*MassWB + 144*g2p2*g3p2*(MassG +    & 
! &  MassWB) + 5*g1p2*(16*g3p2*(MassB + MassG) + 9*g2p2*(MassB + MassWB)) + 288*g3p4*MassG*NGHg3 +& 
! &  108*g2p4*MassWB*NGHw3 + 75*g1p4*MassB*(NGHx3 +          & 
! &  NGHxb3) + 3*(32*g3p4*MassG + 27*g2p4*MassWB)*(NGHx3 +   & 
! &  NGHxb3)))/9._dp + (18*MXM3CYx3TpYx3*(4*g1p2*MassB - 3._dp*(TradjYb3TYb3) -& 
! &  30._dp*(TradjYuTYu) - 15._dp*(TradjYw3TYw3) - 30._dp*(TradjYx3TYx3)) - 9*(BMXM3CYx3TpYx3 +& 
! &  2._dp*(MXM3CYx3TpTYx3))*(4._dp*(g1p2) + 3._dp*(TrYb3adjYb3) + 30._dp*(TrYuadjYu) +    & 
! &  15._dp*(TrYw3adjYw3) + 30._dp*(TrYx3adjYx3)) + 5*BMXM3*(223._dp*(g1p4) +              & 
! &  90*g1p2*g2p2 + 135._dp*(g2p4) + 32*(5._dp*(g1p2) + 9._dp*(g2p2))*g3p2 -               & 
! &  32._dp*(g3p4) + 108*g2p4*NGHw3 + 3*(25._dp*(g1p4) + 27._dp*(g2p4))*(NGHx3 +& 
! &  NGHxb3) + 96*g3p4*(3*NGHg3 + NGHx3 +     & 
! &  NGHxb3)))/90._dp
! 
!  
! DBMXM3 = oo16pi2*( betaBMXM31 + oo16pi2 * betaBMXM32 ) 
! 
!  
! Else 
! DBMXM3 = oo16pi2* betaBMXM31 
! End If 
  DBMXM3 = 0._dp
 
!-------------------- 
! BMWM3 
!-------------------- 
 
! betaBMWM31  = BMWM3CYw3TpYw3 + 2._dp*(MWM3CYw3TpTYw3) + 2._dp*(TYw3adjYw3MWM3)        & 
! &  + Yw3adjYw3BMWM3 + 8*g2p2*(2*MassWB*MWM3 - BMWM3)
! 
!  
!  
! If (TwoLoopRGE) Then 
! betaBMWM32 = (-3._dp*(BMWM3CYw3TpYb3CYb3TpYw3) - 10._dp*(BMWM3CYw3TpYeCYeTpYw3) - 15._dp*(BMWM3CYw3TpYw3CYw3TpYw3) +& 
! &  6*BMWM3CYw3TpYw3*g1p2 - 10*BMWM3CYw3TpYw3*g2p2 - 6._dp*(MWM3CYw3TpTYb3CYb3TpYw3) -    & 
! &  20._dp*(MWM3CYw3TpTYeCYeTpYw3) - 30._dp*(MWM3CYw3TpTYw3CYw3TpYw3) - 6._dp*(MWM3CYw3TpYb3CYb3TpTYw3) -& 
! &  20._dp*(MWM3CYw3TpYeCYeTpTYw3) - 30._dp*(MWM3CYw3TpYw3CYw3TpTYw3) - 2*MWM3CYw3TpYw3*(6*g1p2*MassB -& 
! &  10*g2p2*MassWB + 3._dp*(TradjYb3TYb3) + 30._dp*(TradjYuTYu) + 15._dp*(TradjYw3TYw3) + & 
! &  30._dp*(TradjYx3TYx3)) - 3*BMWM3CYw3TpYw3*TrYb3adjYb3 - 30*BMWM3CYw3TpYw3*TrYuadjYu - & 
! &  15*BMWM3CYw3TpYw3*TrYw3adjYw3 - 30*BMWM3CYw3TpYw3*TrYx3adjYx3 - MWM3CYw3TpTYw3*(-     & 
! & 12._dp*(g1p2) + 20._dp*(g2p2) + 6._dp*(TrYb3adjYb3) + 60._dp*(TrYuadjYu) +             & 
! &  30._dp*(TrYw3adjYw3) + 60._dp*(TrYx3adjYx3)) - 6._dp*(TYw3adjYb3Yb3adjYw3MWM3) -      & 
! &  20._dp*(TYw3adjYeYeadjYw3MWM3) + 2*(6._dp*(g1p2) - 10._dp*(g2p2) - 3._dp*(TrYb3adjYb3) -& 
! &  30._dp*(TrYuadjYu) - 15._dp*(TrYw3adjYw3) - 30._dp*(TrYx3adjYx3))*TYw3adjYw3MWM3 -    & 
! &  30._dp*(TYw3adjYw3Yw3adjYw3MWM3) - 6._dp*(Yw3adjYb3TYb3adjYw3MWM3) - 3._dp*(Yw3adjYb3Yb3adjYw3BMWM3) -& 
! &  20._dp*(Yw3adjYeTYeadjYw3MWM3) - 10._dp*(Yw3adjYeYeadjYw3BMWM3) + 6*g1p2*Yw3adjYw3BMWM3 -& 
! &  10*g2p2*Yw3adjYw3BMWM3 - 3*TrYb3adjYb3*Yw3adjYw3BMWM3 - 30*TrYuadjYu*Yw3adjYw3BMWM3 - & 
! &  15*TrYw3adjYw3*Yw3adjYw3BMWM3 - 30*TrYx3adjYx3*Yw3adjYw3BMWM3 - 12*g1p2*MassB*Yw3adjYw3MWM3 +& 
! &  20*g2p2*MassWB*Yw3adjYw3MWM3 - 6*TradjYb3TYb3*Yw3adjYw3MWM3 - 60*TradjYuTYu*Yw3adjYw3MWM3 -& 
! &  30*TradjYw3TYw3*Yw3adjYw3MWM3 - 60*TradjYx3TYx3*Yw3adjYw3MWM3 - 30._dp*(Yw3adjYw3TYw3adjYw3MWM3) -& 
! &  15._dp*(Yw3adjYw3Yw3adjYw3BMWM3))/10._dp - 16*g2p4*MassWB*MWM3*(10 + 4*NGHw3 +& 
! &  3*NGHx3 + 3*NGHxb3) + 4*g2p4*BMWM3*(10 + 4*NGHw3 +& 
! &  3*NGHx3 + 3*NGHxb3)
! 
!  
! DBMWM3 = oo16pi2*( betaBMWM31 + oo16pi2 * betaBMWM32 ) 
! 
!  
! Else 
! DBMWM3 = oo16pi2* betaBMWM31 
! End If 
 DBMWM3 = 0._dp
 
!-------------------- 
! BMGM3 
!-------------------- 
 
! betaBMGM31  = 12*g3p2*(2*MassG*MGM3 - BMGM3)
! 
!  
!  
! If (TwoLoopRGE) Then 
! betaBMGM32 = -12*g3p4*(4*MassG*MGM3 - BMGM3)*(3 + 3*NGHg3 + NGHx3 +  & 
! &  NGHxb3)
! 
!  
! DBMGM3 = oo16pi2*( betaBMGM31 + oo16pi2 * betaBMGM32 ) 
! 
!  
! Else 
! DBMGM3 = oo16pi2* betaBMGM31 
! End If 
 DBMGM3 = 0._dp
 
!-------------------- 
! BMBM3 
!-------------------- 
 
! betaBMBM31  = (3*(BMBM3CYb3TpYb3 + 2._dp*(MBM3CYb3TpTYb3) + 2._dp*(TYb3adjYb3MBM3)    & 
! &  + Yb3adjYb3BMBM3))/5._dp
! 
!  
!  
! If (TwoLoopRGE) Then 
! betaBMBM32 = (-3*(3._dp*(BMBM3CYb3TpYb3CYb3TpYb3) + 10._dp*(BMBM3CYb3TpYeCYeTpYb3) +               & 
! &  15._dp*(BMBM3CYb3TpYw3CYw3TpYb3) + 6._dp*(MBM3CYb3TpTYb3CYb3TpYb3) + 20._dp*(MBM3CYb3TpTYeCYeTpYb3) +& 
! &  30._dp*(MBM3CYb3TpTYw3CYw3TpYb3) + 6._dp*(MBM3CYb3TpYb3CYb3TpTYb3) + 20._dp*(MBM3CYb3TpYeCYeTpTYb3) +& 
! &  30._dp*(MBM3CYb3TpYw3CYw3TpTYb3) + 6*MBM3CYb3TpYb3*(2*g1p2*MassB + 10*g2p2*MassWB +   & 
! &  TradjYb3TYb3 + 10._dp*(TradjYuTYu) + 5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3)) +  & 
! &  6*MBM3CYb3TpTYb3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) -            & 
! &  2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3))) + 6._dp*(TYb3adjYb3Yb3adjYb3MBM3) +    & 
! &  20._dp*(TYb3adjYeYeadjYb3MBM3) + 30._dp*(TYb3adjYw3Yw3adjYb3MBM3) + 3*((TrYb3adjYb3 + & 
! &  10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*(BMBM3CYb3TpYb3 +& 
! &  2._dp*(TYb3adjYb3MBM3) + Yb3adjYb3BMBM3) + 2*(2*g1p2*MassB + 10*g2p2*MassWB +         & 
! &  TradjYb3TYb3 + 10._dp*(TradjYuTYu) + 5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3))*Yb3adjYb3MBM3) +& 
! &  6._dp*(Yb3adjYb3TYb3adjYb3MBM3) + 3._dp*(Yb3adjYb3Yb3adjYb3BMBM3) + 20._dp*(Yb3adjYeTYeadjYb3MBM3) +& 
! &  10._dp*(Yb3adjYeYeadjYb3BMBM3) + 30._dp*(Yb3adjYw3TYw3adjYb3MBM3) + 15._dp*(Yb3adjYw3Yw3adjYb3BMBM3)))/50._dp
! 
!  
! DBMBM3 = oo16pi2*( betaBMBM31 + oo16pi2 * betaBMBM32 ) 
! 
!  
! Else 
! DBMBM3 = oo16pi2* betaBMBM31 
! End If 
 DBMBM3 = 0._dp
 
!-------------------- 
! mq2 
!-------------------- 
 
betamq21  = mq2TpYdCYd + mq2TpYuCYu + 2._dp*(TpTYdCTYd) + 2._dp*(TpTYuCTYu)           & 
&  + 2*mHd2*TpYdCYd + TpYdCYdmq2 + 2._dp*(TpYdmd2CYd) + 2*mHu2*TpYuCYu + TpYuCYumq2 +    & 
&  2._dp*(TpYumu2CYu) - (id3R*(90*AbsMassWB*g2p2 + 160*AbsMassG*g3p2 + g1p2*(2._dp*(AbsMassB)& 
&  - 5*Tr1(1))))/15._dp

 
 
If (TwoLoopRGE) Then 
betamq22 = (2*AbsMassWB*g1p2*g2p2*id3R)/5._dp + 33*AbsMassWB*g2p4*id3R + 32*AbsMassWB*g2p2*g3p2*id3R +& 
&  (2*g1p2*mq2TpYdCYd)/5._dp + (4*g1p2*mq2TpYuCYu)/5._dp + (4*g1p2*TpTYdCTYd)/5._dp +    & 
&  (8*g1p2*TpTYuCTYu)/5._dp - 4*mHd2*TpYdadjYx3Yx3CYd - 4*mHu2*TpYdadjYx3Yx3CYd -        & 
&  (4*g1p2*MassB*TpYdCTYd)/5._dp + (4*g1p2*mHd2*TpYdCYd)/5._dp + (2*g1p2*TpYdCYdmq2)/5._dp -& 
&  8*mHd2*TpYdCYdTpYdCYd + (4*g1p2*TpYdmd2CYd)/5._dp - (8*g1p2*MassB*TpYuCTYu)/5._dp +   & 
&  (8*g1p2*mHu2*TpYuCYu)/5._dp + (4*g1p2*TpYuCYumq2)/5._dp - 8*mHu2*TpYuCYuTpYuCYu +     & 
&  (8*g1p2*TpYumu2CYu)/5._dp - (3*TpYuCTYu*TradjYb3TYb3)/5._dp - 6*TpYdCTYd*TradjYdTYd - & 
&  2*TpYdCTYd*TradjYeTYe - 6*TpYuCTYu*TradjYuTYu - 3*TpYuCTYu*TradjYw3TYw3 -             & 
&  6*TpYuCTYu*TradjYx3TYx3 - (3*TpYuCYu*TrCTYb3TpTYb3)/5._dp - 6*TpYdCYd*TrCTYdTpTYd -   & 
&  2*TpYdCYd*TrCTYeTpTYe - 2*TpTYdCYd*(3._dp*(TrCTYdTpYd) + TrCTYeTpYe) - 6*TpYuCYu*TrCTYuTpTYu -& 
&  3*TpYuCYu*TrCTYw3TpTYw3 - 6*TpYuCYu*TrCTYx3TpTYx3 - (3*TpTYuCYu*(TrCTYb3TpYb3 +       & 
&  5*(2._dp*(TrCTYuTpYu) + TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3))))/5._dp - 6*TpYuCYu*Trmd2adjYx3Yx3 -& 
&  6*TpYdCYd*Trmd2YdadjYd - 2*TpYdCYd*Trme2YeadjYe - (3*TpYuCYu*TrmHb32Yb3adjYb3)/5._dp -& 
&  3*TpYuCYu*TrmHw32Yw3adjYw3 - 6*TpYuCYu*TrmHxb32Yx3adjYx3 - (3*TpYuCYu*Trml2adjYb3Yb3)/5._dp -& 
&  2*TpYdCYd*Trml2adjYeYe - 3*TpYuCYu*Trml2adjYw3Yw3 - 6*TpYdCYd*Trmq2adjYdYd -          & 
&  6*TpYuCYu*Trmq2adjYuYu - 6*TpYuCYu*Trmu2YuadjYu - (3*mq2TpYuCYu*TrYb3adjYb3)/10._dp - & 
&  (3*TpTYuCTYu*TrYb3adjYb3)/5._dp - (6*mHu2*TpYuCYu*TrYb3adjYb3)/5._dp - (3*TpYuCYumq2*TrYb3adjYb3)/10._dp -& 
&  (3*TpYumu2CYu*TrYb3adjYb3)/5._dp - 3*mq2TpYdCYd*TrYdadjYd - 6*TpTYdCTYd*TrYdadjYd -   & 
&  12*mHd2*TpYdCYd*TrYdadjYd - 3*TpYdCYdmq2*TrYdadjYd - 6*TpYdmd2CYd*TrYdadjYd -         & 
&  mq2TpYdCYd*TrYeadjYe - 2*TpTYdCTYd*TrYeadjYe - 4*mHd2*TpYdCYd*TrYeadjYe -             & 
&  TpYdCYdmq2*TrYeadjYe - 2*TpYdmd2CYd*TrYeadjYe - 3*mq2TpYuCYu*TrYuadjYu  
betamq22 =  betamq22- 6*TpTYuCTYu*TrYuadjYu - 12*mHu2*TpYuCYu*TrYuadjYu - 3*TpYuCYumq2*TrYuadjYu -          & 
&  6*TpYumu2CYu*TrYuadjYu - (3*mq2TpYuCYu*TrYw3adjYw3)/2._dp - 3*TpTYuCTYu*TrYw3adjYw3 - & 
&  6*mHu2*TpYuCYu*TrYw3adjYw3 - (3*TpYuCYumq2*TrYw3adjYw3)/2._dp - 3*TpYumu2CYu*TrYw3adjYw3 -& 
&  3*mq2TpYuCYu*TrYx3adjYx3 - 6*TpTYuCTYu*TrYx3adjYx3 - 12*mHu2*TpYuCYu*TrYx3adjYx3 -    & 
&  3*TpYuCYumq2*TrYx3adjYx3 - 6*TpYumu2CYu*TrYx3adjYx3 - 2*(2._dp*(adjTYdTYdTpYdCYd) +   & 
&  2._dp*(adjTYuTYuTpYuCYu) + mq2TpYdCYdTpYdCYd + mq2TpYuCYuTpYuCYu + 2._dp*(TpTYdadjYdTpYx3CTYx3) +& 
&  2._dp*(TpTYdadjYx3Yx3CTYd) + 2._dp*(TpTYdCYdTpYdCTYd) + 2._dp*(TpTYuCYuTpYuCTYu) +    & 
&  2._dp*(TpYdadjYdTpTYx3CTYx3) + 2._dp*(TpYdadjYx3mHxb32Yx3CYd) + 2._dp*(TpYdadjYx3TYx3CTYd) +& 
&  TpYdadjYx3Yx3CYdmq2 + 2._dp*(TpYdadjYx3Yx3md2CYd) + 2._dp*(TpYdCTYdTpTYdCYd) +        & 
&  2._dp*(TpYdCYdmq2TpYdCYd) + 2._dp*(TpYdCYdTpTYdCTYd) + TpYdCYdTpYdCYdmq2 +            & 
&  2._dp*(TpYdCYdTpYdmd2CYd) + 2._dp*(TpYdmd2adjYx3Yx3CYd) + 2._dp*(TpYdmd2CYdTpYdCYd) + & 
&  2._dp*(TpYuCTYuTpTYuCYu) + 2._dp*(TpYuCYuTpTYuCTYu) + TpYuCYuTpYuCYumq2 +             & 
&  2*(TpYuCYumq2TpYuCYu + TpYuCYuTpYumu2CYu + TpYumu2CYuTpYuCYu) + Ydmq2adjYx3Yx3CYd) +  & 
&  (g1p2*g2p2*id3R*MassB*Conjg(MassWB))/5._dp + 16*g2p2*g3p2*id3R*MassG*Conjg(MassWB) +  & 
&  36*AbsMassWB*g2p4*id3R*NGHw3 + 27*AbsMassWB*g2p4*id3R*NGHx3 +& 
&  27*AbsMassWB*g2p4*id3R*NGHxb3 + (16*g3p2*id3R*Conjg(MassG)*(g1p2*(MassB +& 
&  2._dp*(MassG)) + 15*(-8*g3p2*MassG + 3*g2p2*(2._dp*(MassG) + MassWB)) +               & 
&  90*g3p2*MassG*(3*NGHg3 + NGHx3 + NGHxb3)))/45._dp +& 
&  (g1p2*Conjg(MassB)*(-180*(TpTYdCYd + 2._dp*(TpTYuCYu)) + 360*MassB*(TpYdCYd +         & 
&  2._dp*(TpYuCYu)) + id3R*(597*g1p2*MassB + 5*(16*g3p2*(2._dp*(MassB) + MassG) +        & 
&  9*g2p2*(2._dp*(MassB) + MassWB)) + 225*g1p2*MassB*(NGHx3 +             & 
&  NGHxb3))))/225._dp + (2*g1p4*id3R*Tr2(1))/15._dp + 6*g2p4*id3R*Tr2(2)  
betamq22 =  betamq22+ (32*g3p4*id3R*Tr2(3))/3._dp + (4*g1p2*id3R*Tr3(1))/3._dp

 
Dmq2 = oo16pi2*( betamq21 + oo16pi2 * betamq22 ) 

 
Else 
Dmq2 = oo16pi2* betamq21 
End If 
 
 
Forall(i1=1:3) Dmq2(i1,i1) =  Real(Dmq2(i1,i1),dp) 
!-------------------- 
! ml2 
!-------------------- 
 
betaml21  = 3._dp*(ml2TpYb3CYb3)/10._dp + ml2TpYeCYe + 3._dp*(ml2TpYw3CYw3)           & 
& /2._dp + 3._dp*(TpTYb3CTYb3)/5._dp + 2._dp*(TpTYeCTYe) + 3._dp*(TpTYw3CTYw3)           & 
&  + (3*mHu2*TpYb3CYb3)/5._dp + 3._dp*(TpYb3CYb3ml2)/10._dp + 3._dp*(TpYb3mHb32CYb3)     & 
& /5._dp + 2*mHd2*TpYeCYe + TpYeCYeml2 + 2._dp*(TpYeme2CYe) + 3*mHu2*TpYw3CYw3 +         & 
&  3._dp*(TpYw3CYw3ml2)/2._dp + 3._dp*(TpYw3mHw32CYw3) - (id3R*(30*AbsMassWB*g2p2 +      & 
&  g1p2*(6._dp*(AbsMassB) + 5*Tr1(1))))/5._dp

 
 
If (TwoLoopRGE) Then 
betaml22 = (6*g1p2*ml2TpYeCYe)/5._dp + 6*g2p2*ml2TpYw3CYw3 + (12*g1p2*TpTYeCTYe)/5._dp +         & 
&  12*g2p2*TpTYw3CTYw3 - (18*mHu2*TpYb3CYb3TpYb3CYb3)/25._dp - (9*mHu2*TpYb3CYb3TpYw3CYw3)/10._dp -& 
&  (12*g1p2*MassB*TpYeCTYe)/5._dp + (12*g1p2*mHd2*TpYeCYe)/5._dp + (6*g1p2*TpYeCYeml2)/5._dp -& 
&  8*mHd2*TpYeCYeTpYeCYe + (12*g1p2*TpYeme2CYe)/5._dp - 12*g2p2*MassWB*TpYw3CTYw3 +      & 
&  12*g2p2*mHu2*TpYw3CYw3 + 6*g2p2*TpYw3CYw3ml2 - (9*mHu2*TpYw3CYw3TpYb3CYb3)/10._dp -   & 
&  (9*mHu2*TpYw3CYw3TpYw3CYw3)/2._dp + 12*g2p2*TpYw3mHw32CYw3 - (9*TpYw3CTYw3*TradjYb3TYb3)/10._dp -& 
&  6*TpYeCTYe*TradjYdTYd - 2*TpYeCTYe*TradjYeTYe - 9*TpYw3CTYw3*TradjYuTYu -             & 
&  (9*TpYw3CTYw3*TradjYw3TYw3)/2._dp - 9*TpYw3CTYw3*TradjYx3TYx3 - (9*TpYw3CYw3*TrCTYb3TpTYb3)/10._dp -& 
&  6*TpYeCYe*TrCTYdTpTYd - 2*TpYeCYe*TrCTYeTpTYe - 9*TpYw3CYw3*TrCTYuTpTYu -             & 
&  (9*TpYw3CYw3*TrCTYw3TpTYw3)/2._dp - 9*TpYw3CYw3*TrCTYx3TpTYx3 - 9*TpYw3CYw3*Trmd2adjYx3Yx3 -& 
&  6*TpYeCYe*Trmd2YdadjYd - 2*TpYeCYe*Trme2YeadjYe - (9*TpYw3CYw3*TrmHb32Yb3adjYb3)/10._dp -& 
&  (9*TpYw3CYw3*TrmHw32Yw3adjYw3)/2._dp - 9*TpYw3CYw3*TrmHxb32Yx3adjYx3 - (9*TpYw3CYw3*Trml2adjYb3Yb3)/10._dp -& 
&  2*TpYeCYe*Trml2adjYeYe - (9*TpYw3CYw3*Trml2adjYw3Yw3)/2._dp - 6*TpYeCYe*Trmq2adjYdYd -& 
&  9*TpYw3CYw3*Trmq2adjYuYu - 9*TpYw3CYw3*Trmu2YuadjYu - (9*ml2TpYb3CYb3*TrYb3adjYb3)/100._dp -& 
&  (9*ml2TpYw3CYw3*TrYb3adjYb3)/20._dp - (9*TpTYb3CTYb3*TrYb3adjYb3)/50._dp -            & 
&  (9*TpTYw3CTYw3*TrYb3adjYb3)/10._dp - (9*TpYb3CYb3ml2*TrYb3adjYb3)/100._dp -           & 
&  (9*TpYb3mHb32CYb3*TrYb3adjYb3)/50._dp - (9*mHu2*TpYw3CYw3*TrYb3adjYb3)/5._dp -        & 
&  (9*TpYw3CYw3ml2*TrYb3adjYb3)/20._dp - (9*TpYw3mHw32CYw3*TrYb3adjYb3)/10._dp -         & 
&  3*ml2TpYeCYe*TrYdadjYd - 6*TpTYeCTYe*TrYdadjYd - 12*mHd2*TpYeCYe*TrYdadjYd -          & 
&  3*TpYeCYeml2*TrYdadjYd - 6*TpYeme2CYe*TrYdadjYd - ml2TpYeCYe*TrYeadjYe -              & 
&  2*TpTYeCTYe*TrYeadjYe - 4*mHd2*TpYeCYe*TrYeadjYe - TpYeCYeml2*TrYeadjYe  
betaml22 =  betaml22- 2*TpYeme2CYe*TrYeadjYe - (9*ml2TpYb3CYb3*TrYuadjYu)/10._dp - (9*ml2TpYw3CYw3*TrYuadjYu)/2._dp -& 
&  (9*TpTYb3CTYb3*TrYuadjYu)/5._dp - 9*TpTYw3CTYw3*TrYuadjYu - (9*TpYb3CYb3ml2*TrYuadjYu)/10._dp -& 
&  (9*TpYb3mHb32CYb3*TrYuadjYu)/5._dp - 18*mHu2*TpYw3CYw3*TrYuadjYu - (9*TpYw3CYw3ml2*TrYuadjYu)/2._dp -& 
&  9*TpYw3mHw32CYw3*TrYuadjYu - (9*ml2TpYb3CYb3*TrYw3adjYw3)/20._dp - (9*ml2TpYw3CYw3*TrYw3adjYw3)/4._dp -& 
&  (9*TpTYb3CTYb3*TrYw3adjYw3)/10._dp - (9*TpTYw3CTYw3*TrYw3adjYw3)/2._dp -              & 
&  (9*TpYb3CYb3ml2*TrYw3adjYw3)/20._dp - (9*TpYb3mHb32CYb3*TrYw3adjYw3)/10._dp -         & 
&  9*mHu2*TpYw3CYw3*TrYw3adjYw3 - (9*TpYw3CYw3ml2*TrYw3adjYw3)/4._dp - (9*TpYw3mHw32CYw3*TrYw3adjYw3)/2._dp -& 
&  (9*ml2TpYb3CYb3*TrYx3adjYx3)/10._dp - (9*ml2TpYw3CYw3*TrYx3adjYx3)/2._dp -            & 
&  (9*TpTYb3CTYb3*TrYx3adjYx3)/5._dp - 9*TpTYw3CTYw3*TrYx3adjYx3 - (9*TpYb3CYb3ml2*TrYx3adjYx3)/10._dp -& 
&  (9*TpYb3mHb32CYb3*TrYx3adjYx3)/5._dp - 18*mHu2*TpYw3CYw3*TrYx3adjYx3 - (9*TpYw3CYw3ml2*TrYx3adjYx3)/2._dp -& 
&  9*TpYw3mHw32CYw3*TrYx3adjYx3 + (-72._dp*(adjTYb3TYb3TpYb3CYb3) - 90._dp*(adjTYb3TYb3TpYw3CYw3) -& 
&  800._dp*(adjTYeTYeTpYeCYe) - 90._dp*(adjTYw3TYw3TpYb3CYb3) - 450._dp*(adjTYw3TYw3TpYw3CYw3) -& 
&  36._dp*(ml2TpYb3CYb3TpYb3CYb3) - 45._dp*(ml2TpYb3CYb3TpYw3CYw3) - 400._dp*(ml2TpYeCYeTpYeCYe) -& 
&  45._dp*(ml2TpYw3CYw3TpYb3CYb3) - 225._dp*(ml2TpYw3CYw3TpYw3CYw3) - 72._dp*(TpTYb3CYb3TpYb3CTYb3) -& 
&  90._dp*(TpTYb3CYb3TpYw3CTYw3) - 800._dp*(TpTYeCYeTpYeCTYe) - 90._dp*(TpTYw3CYw3TpYb3CTYb3) -& 
&  450._dp*(TpTYw3CYw3TpYw3CTYw3) - 72._dp*(TpYb3CTYb3TpTYb3CYb3) - 90._dp*(TpYb3CTYb3TpTYw3CYw3) -& 
&  72._dp*(TpYb3CYb3ml2TpYb3CYb3) - 90._dp*(TpYb3CYb3ml2TpYw3CYw3) - 72._dp*(TpYb3CYb3TpTYb3CTYb3) -& 
&  90._dp*(TpYb3CYb3TpTYw3CTYw3) - 36._dp*(TpYb3CYb3TpYb3CYb3ml2) - 72._dp*(TpYb3CYb3TpYb3mHb32CYb3) -& 
&  72._dp*(TpYb3mHb32CYb3TpYb3CYb3) - 90._dp*(TpYb3mHb32CYb3TpYw3CYw3) - 800._dp*(TpYeCTYeTpTYeCYe) -& 
&  800._dp*(TpYeCYeTpTYeCTYe) - 90._dp*(TpYw3CTYw3TpTYb3CYb3) - 450._dp*(TpYw3CTYw3TpTYw3CYw3) -& 
&  90._dp*(TpYw3CYw3TpTYb3CTYb3) - 450._dp*(TpYw3CYw3TpTYw3CTYw3) - 5*(9._dp*(TpYb3CYb3TpYw3CYw3ml2) +& 
&  18._dp*(TpYb3CYb3TpYw3mHw32CYw3) + 160._dp*(TpYeCYeml2TpYeCYe) + 80._dp*(TpYeCYeTpYeCYeml2) +& 
&  160._dp*(TpYeCYeTpYeme2CYe) + 160._dp*(TpYeme2CYeTpYeCYe) + 9*(2._dp*(TpYw3CYw3ml2TpYb3CYb3) +& 
&  10._dp*(TpYw3CYw3ml2TpYw3CYw3) + TpYw3CYw3TpYb3CYb3ml2 + 2._dp*(TpYw3CYw3TpYb3mHb32CYb3) +& 
&  5._dp*(TpYw3CYw3TpYw3CYw3ml2) + 10._dp*(TpYw3CYw3TpYw3mHw32CYw3) + 2._dp*(TpYw3mHw32CYw3TpYb3CYb3) +& 
&  10._dp*(TpYw3mHw32CYw3TpYw3CYw3))) - 36*TpYb3CTYb3*(TradjYb3TYb3 + 5*(2._dp*(TradjYuTYu) +& 
&  TradjYw3TYw3 + 2._dp*(TradjYx3TYx3))) - 400*TpTYeCYe*(3._dp*(TrCTYdTpYd) +            & 
&  TrCTYeTpYe) - 36*(TpTYb3CYb3 + 5._dp*(TpTYw3CYw3))*(TrCTYb3TpYb3 + 5*(2._dp*(TrCTYuTpYu) +& 
&  TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3))) - 36*TpYb3CYb3*(TrCTYb3TpTYb3 + 10._dp*(TrCTYuTpTYu) +& 
&  5._dp*(TrCTYw3TpTYw3) + 10._dp*(TrCTYx3TpTYx3) + 10._dp*(Trmd2adjYx3Yx3) +            & 
&  TrmHb32Yb3adjYb3 + 5._dp*(TrmHw32Yw3adjYw3) + 10._dp*(TrmHxb32Yx3adjYx3) +            & 
&  Trml2adjYb3Yb3 + 5*(Trml2adjYw3Yw3 + 2*(Trmq2adjYuYu + Trmu2YuadjYu)) +               & 
&  2*mHu2*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3)))))/200._dp  
betaml22 =  betaml22+ (3*g1p2*Conjg(MassB)*(-20._dp*(TpTYeCYe) + 40*MassB*TpYeCYe + 3*id3R*(69*g1p2*MassB + & 
&  5*g2p2*(2._dp*(MassB) + MassWB) + 25*g1p2*MassB*(NGHx3 +               & 
&  NGHxb3))))/25._dp + (3*g2p2*Conjg(MassWB)*(-20._dp*(TpTYw3CYw3) +      & 
&  40*MassWB*TpYw3CYw3 + id3R*(55*g2p2*MassWB + 3*g1p2*(MassB + 2._dp*(MassWB)) +        & 
&  15*g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/5._dp +& 
&  (6*g1p4*id3R*Tr2(1))/5._dp + 6*g2p4*id3R*Tr2(2) - 4*g1p2*id3R*Tr3(1)

 
Dml2 = oo16pi2*( betaml21 + oo16pi2 * betaml22 ) 

 
Else 
Dml2 = oo16pi2* betaml21 
End If 
 
 
Forall(i1=1:3) Dml2(i1,i1) =  Real(Dml2(i1,i1),dp) 
!-------------------- 
! mHd2 
!-------------------- 
 
betamHd21  = 2*(-3*AbsMassWB*g2p2 + 3._dp*(TrCTYdTpTYd) + TrCTYeTpTYe +               & 
&  3._dp*(Trmd2YdadjYd) + Trme2YeadjYe + Trml2adjYeYe + 3._dp*(Trmq2adjYdYd)             & 
&  + 3*mHd2*TrYdadjYd + mHd2*TrYeadjYe) - (g1p2*(6._dp*(AbsMassB) + 5*Tr1(1)))/5._dp

 
 
If (TwoLoopRGE) Then 
betamHd22 = (g1p2*Conjg(MassB)*(9*(69*g1p2*MassB + 5*g2p2*(2._dp*(MassB) + MassWB)) +             & 
&  20._dp*(TradjYdTYd) - 60._dp*(TradjYeTYe) + 5*MassB*(-8*(TrYdadjYd - 3._dp*(TrYeadjYe)) +& 
&  45*g1p2*(NGHx3 + NGHxb3))) + 5*(-4*g1p2*TrCTYdTpTYd +   & 
&  160*g3p2*TrCTYdTpTYd + 4*g1p2*MassB*TrCTYdTpYd - 160*g3p2*MassG*TrCTYdTpYd +          & 
&  12*g1p2*TrCTYeTpTYe - 12*g1p2*MassB*TrCTYeTpYe - 4*g1p2*Trmd2YdadjYd + 160*g3p2*Trmd2YdadjYd +& 
&  12*g1p2*Trme2YeadjYe + 12*g1p2*Trml2adjYeYe - 4*g1p2*Trmq2adjYdYd + 160*g3p2*Trmq2adjYdYd -& 
&  3*mHd2*TrYb3adjYeYeadjYb3 - 3*mHu2*TrYb3adjYeYeadjYb3 - 4*(-80*AbsMassG*g3p2 +        & 
&  (g1p2 - 40._dp*(g3p2))*mHd2)*TrYdadjYd - 60*mHd2*TrYdadjYdTpYx3CYx3 - 60*mHu2*TrYdadjYdTpYx3CYx3 -& 
&  180*mHd2*TrYdadjYdYdadjYd - 30*mHd2*TrYdadjYuYuadjYd - 30*mHu2*TrYdadjYuYuadjYd +     & 
&  12*g1p2*mHd2*TrYeadjYe - 60*mHd2*TrYeadjYeYeadjYe - 15*(mHd2 + mHu2)*TrYeadjYw3Yw3adjYe -& 
&  3*(20._dp*(TradjYdTpYx3CTYx3TYd) + 20._dp*(TradjYx3TYx3CTYdTpYd) + 20._dp*(Trmd2TpYx3CYx3YdadjYd) +& 
&  20._dp*(Trmd2YdadjYdTpYx3CYx3) + 60._dp*(Trmd2YdadjYdYdadjYd) + 10._dp*(Trmd2YdadjYuYuadjYd) +& 
&  Trme2YeadjYb3Yb3adjYe + 20._dp*(Trme2YeadjYeYeadjYe) + 5._dp*(Trme2YeadjYw3Yw3adjYe) +& 
&  TrmHb32Yb3adjYeYeadjYb3 + 5._dp*(TrmHw32Yw3adjYeYeadjYw3) + 20._dp*(TrmHxb32Yx3CYdTpYdadjYx3) +& 
&  Trml2adjYb3Yb3adjYeYe + Trml2adjYeYeadjYb3Yb3 + 5*(4._dp*(Trml2adjYeYeadjYeYe) +      & 
&  Trml2adjYeYeadjYw3Yw3 + Trml2adjYw3Yw3adjYeYe + 2*(2._dp*(Trmq2adjYdTpYx3CYx3Yd) +    & 
&  6._dp*(Trmq2adjYdYdadjYdYd) + Trmq2adjYdYdadjYuYu + Trmq2adjYuYuadjYdYd +             & 
&  Trmu2YuadjYdYdadjYu)) + TrYb3adjYeTYeadjTYb3 + TrYb3TpTYeCTYeadjYb3 + 30._dp*(TrYdadjTYdTYdadjYd) +& 
&  10._dp*(TrYdadjTYuTYuadjYd) + 20._dp*(TrYdadjYdTpTYx3CTYx3) + 60._dp*(TrYdadjYdTYdadjTYd) +& 
&  10._dp*(TrYdadjYuTYuadjTYd) + 30._dp*(TrYdTpTYdCTYdadjYd) + TrYeadjTYb3TYb3adjYe +    & 
&  10._dp*(TrYeadjTYeTYeadjYe) + 5._dp*(TrYeadjTYw3TYw3adjYe) + TrYeadjYb3TYb3adjTYe +   & 
&  20._dp*(TrYeadjYeTYeadjTYe) + 5._dp*(TrYeadjYw3TYw3adjTYe) + 10._dp*(TrYeTpTYeCTYeadjYe) +& 
&  10._dp*(TrYuadjYdTYdadjTYu) + 10._dp*(TrYuTpTYdCTYdadjYu) + 5._dp*(TrYw3adjYeTYeadjTYw3) +& 
&  5._dp*(TrYw3TpTYeCTYeadjYw3) + 20._dp*(TrYx3CTYdTpTYdadjYx3)) - 160*g3p2*TradjYdTYd*Conjg(MassG) +& 
&  3*g2p2*Conjg(MassWB)*(55*g2p2*MassWB + 3*g1p2*(MassB + 2._dp*(MassWB)) +              & 
&  15*g2p2*MassWB*(4*NGHw3 + 3*(NGHx3 + NGHxb3))) +& 
&  6*g1p4*Tr2(1) + 30*g2p4*Tr2(2) - 20*g1p2*Tr3(1)))/25._dp

 
DmHd2 = oo16pi2*( betamHd21 + oo16pi2 * betamHd22 ) 

 
Else 
DmHd2 = oo16pi2* betamHd21 
End If 
 
 
!-------------------- 
! mHu2 
!-------------------- 
 
betamHu21  = (-6*AbsMassB*g1p2)/5._dp + (3*(-10*AbsMassWB*g2p2 + TrCTYb3TpTYb3 +      & 
&  10._dp*(TrCTYuTpTYu) + 5._dp*(TrCTYw3TpTYw3) + 10._dp*(TrCTYx3TpTYx3) +               & 
&  10._dp*(Trmd2adjYx3Yx3) + TrmHb32Yb3adjYb3 + 5._dp*(TrmHw32Yw3adjYw3) +               & 
&  10._dp*(TrmHxb32Yx3adjYx3) + Trml2adjYb3Yb3 + 5*(Trml2adjYw3Yw3 + 2*(Trmq2adjYuYu +   & 
&  Trmu2YuadjYu)) + mHu2*TrYb3adjYb3 + 5*mHu2*(2._dp*(TrYuadjYu) + TrYw3adjYw3 +         & 
&  2._dp*(TrYx3adjYx3))))/5._dp + g1p2*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHu22 = (8*g1p2*TrCTYuTpTYu)/5._dp + 32*g3p2*TrCTYuTpTYu - (8*g1p2*MassB*TrCTYuTpYu)/5._dp -  & 
&  32*g3p2*MassG*TrCTYuTpYu + 12*g2p2*TrCTYw3TpTYw3 - 12*g2p2*MassWB*TrCTYw3TpYw3 +      & 
&  4*g1p2*TrCTYx3TpTYx3 + 32*g3p2*TrCTYx3TpTYx3 - 4*g1p2*MassB*TrCTYx3TpYx3 -            & 
&  32*g3p2*MassG*TrCTYx3TpYx3 + 4*g1p2*Trmd2adjYx3Yx3 + 32*g3p2*Trmd2adjYx3Yx3 +         & 
&  12*g2p2*TrmHw32Yw3adjYw3 + 4*g1p2*TrmHxb32Yx3adjYx3 + 32*g3p2*TrmHxb32Yx3adjYx3 +     & 
&  12*g2p2*Trml2adjYw3Yw3 + (8*g1p2*Trmq2adjYuYu)/5._dp + 32*g3p2*Trmq2adjYuYu +         & 
&  (8*g1p2*Trmu2YuadjYu)/5._dp + 32*g3p2*Trmu2YuadjYu - (27*mHu2*TrYb3adjYb3Yb3adjYb3)/25._dp -& 
&  (3*mHu2*TrYb3adjYeYeadjYb3)/5._dp - (27*mHu2*TrYb3adjYw3Yw3adjYb3)/5._dp -            & 
&  12*mHu2*TrYdadjYdTpYx3CYx3 - 6*mHu2*TrYdadjYuYuadjYd - 3*mHu2*TrYeadjYw3Yw3adjYe +    & 
&  64*AbsMassG*g3p2*TrYuadjYu + (8*g1p2*mHu2*TrYuadjYu)/5._dp + 32*g3p2*mHu2*TrYuadjYu - & 
&  36*mHu2*TrYuadjYuYuadjYu + 12*g2p2*mHu2*TrYw3adjYw3 - (27*mHu2*TrYw3adjYw3Yw3adjYw3)/2._dp +& 
&  64*AbsMassG*g3p2*TrYx3adjYx3 + 4*g1p2*mHu2*TrYx3adjYx3 + 32*g3p2*mHu2*TrYx3adjYx3 -   & 
&  36*mHu2*TrYx3adjYx3Yx3adjYx3 - (3*(18._dp*(TrYb3adjTYb3TYb3adjYb3) + 20._dp*(TrYb3adjTYeTYeadjYb3) +& 
&  45._dp*(TrYb3adjTYw3TYw3adjYb3) + 36._dp*(TrYb3adjYb3TYb3adjTYb3) + 20._dp*(TrYb3adjYeTYeadjTYb3) +& 
&  90._dp*(TrYb3adjYw3TYw3adjTYb3) + 18._dp*(TrYb3TpTYb3CTYb3adjYb3) + 45._dp*(TrYb3TpTYw3CTYw3adjYb3) +& 
&  200._dp*(TrYdadjYuTYuadjTYd) + 200._dp*(TrYdTpTYuCTYuadjYd) + 20._dp*(TrYeadjYb3TYb3adjTYe) +& 
&  100._dp*(TrYeadjYw3TYw3adjTYe) + 20*mHd2*(TrYb3adjYeYeadjYb3 + 5*(4._dp*(TrYdadjYdTpYx3CYx3) +& 
&  2._dp*(TrYdadjYuYuadjYd) + TrYeadjYw3Yw3adjYe)) + 20._dp*(TrYeTpTYb3CTYb3adjYe) +     & 
&  100._dp*(TrYeTpTYw3CTYw3adjYe) + 200._dp*(TrYuadjTYdTYdadjYu) + 600._dp*(TrYuadjTYuTYuadjYu) +& 
&  200._dp*(TrYuadjYdTYdadjTYu) + 1200._dp*(TrYuadjYuTYuadjTYu) + 600._dp*(TrYuTpTYuCTYuadjYu) +& 
&  45._dp*(TrYw3adjTYb3TYb3adjYw3) + 100._dp*(TrYw3adjTYeTYeadjYw3) + 225._dp*(TrYw3adjTYw3TYw3adjYw3) +& 
&  90._dp*(TrYw3adjYb3TYb3adjTYw3) + 100._dp*(TrYw3adjYeTYeadjTYw3) + 450._dp*(TrYw3adjYw3TYw3adjTYw3) +& 
&  45._dp*(TrYw3TpTYb3CTYb3adjYw3) + 225._dp*(TrYw3TpTYw3CTYw3adjYw3) + 2*(200._dp*(TradjYdTpTYx3CTYx3TpYd) +& 
&  200._dp*(TradjYdTpYx3CTYx3TpTYd) + 200._dp*(TradjYx3TYx3CTYdTpYd) + 600._dp*(Trmd2adjYx3Yx3adjYx3Yx3) +& 
&  200._dp*(Trmd2TpYx3CYx3YdadjYd) + 200._dp*(Trmd2YdadjYdTpYx3CYx3) + 100._dp*(Trmd2YdadjYuYuadjYd) +& 
&  10._dp*(Trme2YeadjYb3Yb3adjYe) + 50._dp*(Trme2YeadjYw3Yw3adjYe) + 18._dp*(TrmHb32Yb3adjYb3Yb3adjYb3) +& 
&  10._dp*(TrmHb32Yb3adjYeYeadjYb3) + 45._dp*(TrmHb32Yb3adjYw3Yw3adjYb3) +               & 
&  45._dp*(TrmHw32Yw3adjYb3Yb3adjYw3) + 50._dp*(TrmHw32Yw3adjYeYeadjYw3) +               & 
&  225._dp*(TrmHw32Yw3adjYw3Yw3adjYw3) + 200._dp*(TrmHxb32CYx3YdadjYdTpYx3) +            & 
&  600._dp*(TrmHxb32Yx3adjYx3Yx3adjYx3) + 18._dp*(Trml2adjYb3Yb3adjYb3Yb3) +             & 
&  5*(2._dp*(Trml2adjYb3Yb3adjYeYe) + 9._dp*(Trml2adjYb3Yb3adjYw3Yw3) + 2._dp*(Trml2adjYeYeadjYb3Yb3) +& 
&  10._dp*(Trml2adjYeYeadjYw3Yw3) + 9._dp*(Trml2adjYw3Yw3adjYb3Yb3) + 5*(2._dp*(Trml2adjYw3Yw3adjYeYe) +& 
&  9._dp*(Trml2adjYw3Yw3adjYw3Yw3) + 4*(2._dp*(Trmq2adjYdTpYx3CYx3Yd) + Trmq2adjYdYdadjYuYu +& 
&  Trmq2adjYuYuadjYdYd + 6._dp*(Trmq2adjYuYuadjYuYu) + Trmu2YuadjYdYdadjYu +             & 
&  6._dp*(Trmu2YuadjYuYuadjYu)))) + 300._dp*(TrYx3adjTYx3TYx3adjYx3) + 600._dp*(TrYx3adjYx3TYx3adjTYx3) +& 
&  200._dp*(TrYx3CTYdTpTYdadjYx3) + 300._dp*(TrYx3TpTYx3CTYx3adjYx3))))/100._dp  
betamHu22 =  betamHu22- 32*g3p2*TradjYuTYu*Conjg(MassG) - 32*g3p2*TradjYx3TYx3*Conjg(MassG) &
& + (g1p2*Conjg(MassB)*(9*(69*g1p2*MassB +& 
&  5*g2p2*(2._dp*(MassB) + MassWB)) - 40._dp*(TradjYuTYu) - 100._dp*(TradjYx3TYx3) +     & 
&  5*MassB*(8*(2._dp*(TrYuadjYu) + 5._dp*(TrYx3adjYx3)) + 45*g1p2*(NGHx3 +& 
&  NGHxb3))))/25._dp + (3*g2p2*Conjg(MassWB)*(55*g2p2*MassWB +            & 
&  3*g1p2*(MassB + 2._dp*(MassWB)) - 20._dp*(TradjYw3TYw3) + 5*MassWB*(8._dp*(TrYw3adjYw3) +& 
&  3*g2p2*(4*NGHw3 + 3*(NGHx3 + NGHxb3)))))/5._dp +& 
&  (6*g1p4*Tr2(1))/5._dp + 6*g2p4*Tr2(2) + 4*g1p2*Tr3(1)

 
DmHu2 = oo16pi2*( betamHu21 + oo16pi2 * betamHu22 ) 

 
Else 
DmHu2 = oo16pi2* betamHu21 
End If 
 
 
!-------------------- 
! md2 
!-------------------- 
 
betamd21  = 2*(md2TpYx3CYx3 + md2YdadjYd + 2._dp*(TpTYx3CTYx3) + 2*mHu2*TpYx3CYx3 +   & 
&  TpYx3CYx3md2 + 2._dp*(TpYx3mHxb32CYx3) + 2._dp*(TYdadjTYd) + 2*mHd2*YdadjYd +         & 
&  YdadjYdmd2 + 2._dp*(Ydmq2adjYd)) - (2*id3R*(80*AbsMassG*g3p2 + g1p2*(4._dp*(AbsMassB) & 
&  - 5*Tr1(1))))/15._dp

 
 
If (TwoLoopRGE) Then 
betamd22 = 2*g1p2*md2TpYx3CYx3 + 6*g2p2*md2TpYx3CYx3 + (2*g1p2*md2YdadjYd)/5._dp +               & 
&  6*g2p2*md2YdadjYd + 4*g1p2*TpTYx3CTYx3 + 12*g2p2*TpTYx3CTYx3 - 4*g1p2*MassB*TpYx3CTYx3 -& 
&  12*g2p2*MassWB*TpYx3CTYx3 + 24*AbsMassWB*g2p2*TpYx3CYx3 + 4*g1p2*mHu2*TpYx3CYx3 +     & 
&  12*g2p2*mHu2*TpYx3CYx3 + 2*g1p2*TpYx3CYx3md2 + 6*g2p2*TpYx3CYx3md2 - 8*mHu2*TpYx3CYx3TpYx3CYx3 +& 
&  4*g1p2*TpYx3mHxb32CYx3 + 12*g2p2*TpYx3mHxb32CYx3 - (6*TpYx3CTYx3*TradjYb3TYb3)/5._dp -& 
&  12*TpYx3CTYx3*TradjYuTYu - 6*TpYx3CTYx3*TradjYw3TYw3 - 12*TpYx3CTYx3*TradjYx3TYx3 -   & 
&  (6*TpYx3CYx3*TrCTYb3TpTYb3)/5._dp - (6*TpTYx3CYx3*TrCTYb3TpYb3)/5._dp -               & 
&  12*TpYx3CYx3*TrCTYuTpTYu - 12*TpTYx3CYx3*TrCTYuTpYu - 6*TpYx3CYx3*TrCTYw3TpTYw3 -     & 
&  6*TpTYx3CYx3*TrCTYw3TpYw3 - 12*TpYx3CYx3*TrCTYx3TpTYx3 - 12*TpTYx3CYx3*TrCTYx3TpYx3 - & 
&  12*TpYx3CYx3*Trmd2adjYx3Yx3 - (6*TpYx3CYx3*TrmHb32Yb3adjYb3)/5._dp - 6*TpYx3CYx3*TrmHw32Yw3adjYw3 -& 
&  12*TpYx3CYx3*TrmHxb32Yx3adjYx3 - (6*TpYx3CYx3*Trml2adjYb3Yb3)/5._dp - 6*TpYx3CYx3*Trml2adjYw3Yw3 -& 
&  12*TpYx3CYx3*Trmq2adjYuYu - 12*TpYx3CYx3*Trmu2YuadjYu - (3*md2TpYx3CYx3*TrYb3adjYb3)/5._dp -& 
&  (6*TpTYx3CTYx3*TrYb3adjYb3)/5._dp - (12*mHu2*TpYx3CYx3*TrYb3adjYb3)/5._dp -           & 
&  (3*TpYx3CYx3md2*TrYb3adjYb3)/5._dp - (6*TpYx3mHxb32CYx3*TrYb3adjYb3)/5._dp -          & 
&  6*md2YdadjYd*TrYdadjYd - 2*md2YdadjYd*TrYeadjYe - 6*md2TpYx3CYx3*TrYuadjYu -          & 
&  12*TpTYx3CTYx3*TrYuadjYu - 24*mHu2*TpYx3CYx3*TrYuadjYu - 6*TpYx3CYx3md2*TrYuadjYu -   & 
&  12*TpYx3mHxb32CYx3*TrYuadjYu - 3*md2TpYx3CYx3*TrYw3adjYw3 - 6*TpTYx3CTYx3*TrYw3adjYw3 -& 
&  12*mHu2*TpYx3CYx3*TrYw3adjYw3 - 3*TpYx3CYx3md2*TrYw3adjYw3 - 6*TpYx3mHxb32CYx3*TrYw3adjYw3 -& 
&  6*md2TpYx3CYx3*TrYx3adjYx3 - 12*TpTYx3CTYx3*TrYx3adjYx3 - 24*mHu2*TpYx3CYx3*TrYx3adjYx3 -& 
&  6*TpYx3CYx3md2*TrYx3adjYx3 - 12*TpYx3mHxb32CYx3*TrYx3adjYx3 + (4*g1p2*TYdadjTYd)/5._dp +& 
&  12*g2p2*TYdadjTYd - 12*TrYdadjYd*TYdadjTYd - 4*TrYeadjYe*TYdadjTYd - 12*TrCTYdTpYd*TYdadjYd  
betamd22 =  betamd22- 4*TrCTYeTpYe*TYdadjYd - (4*g1p2*MassB*YdadjTYd)/5._dp - 12*g2p2*MassWB*YdadjTYd -     & 
&  12*TradjYdTYd*YdadjTYd - 4*TradjYeTYe*YdadjTYd + 24*AbsMassWB*g2p2*YdadjYd +          & 
&  (4*g1p2*mHd2*YdadjYd)/5._dp + 12*g2p2*mHd2*YdadjYd - 12*TrCTYdTpTYd*YdadjYd -         & 
&  4*TrCTYeTpTYe*YdadjYd - 12*Trmd2YdadjYd*YdadjYd - 4*Trme2YeadjYe*YdadjYd -            & 
&  4*Trml2adjYeYe*YdadjYd - 12*Trmq2adjYdYd*YdadjYd - 24*mHd2*TrYdadjYd*YdadjYd -        & 
&  8*mHd2*TrYeadjYe*YdadjYd + (2*g1p2*YdadjYdmd2)/5._dp + 6*g2p2*YdadjYdmd2 -            & 
&  6*TrYdadjYd*YdadjYdmd2 - 2*TrYeadjYe*YdadjYdmd2 - 8*mHd2*YdadjYdYdadjYd -             & 
&  4*mHd2*YdadjYuYuadjYd - 4*mHu2*YdadjYuYuadjYd + (4*g1p2*Ydmq2adjYd)/5._dp +           & 
&  12*g2p2*Ydmq2adjYd - 12*TrYdadjYd*Ydmq2adjYd - 4*TrYeadjYe*Ydmq2adjYd -               & 
&  2*(2._dp*(adjTYx3TYx3TpYx3CYx3) + md2TpYx3CYx3TpYx3CYx3 + md2YdadjYdYdadjYd +         & 
&  md2YdadjYuYuadjYd + 2._dp*(TpTYx3CYx3TpYx3CTYx3) + 2._dp*(TpYx3CTYx3TpTYx3CYx3) +     & 
&  2._dp*(TpYx3CYx3TpTYx3CTYx3) + TpYx3CYx3TpYx3CYx3md2 + 2*(TpYx3CYx3md2TpYx3CYx3 +     & 
&  TpYx3CYx3TpYx3mHxb32CYx3 + TpYx3mHxb32CYx3TpYx3CYx3) + 2._dp*(TYdadjYdYdadjTYd) +     & 
&  2._dp*(TYdadjYuYuadjTYd) + 2._dp*(TYdTpYdCTYdadjYd) + 2._dp*(TYdTpYuCTYuadjYd) +      & 
&  2._dp*(YdadjTYdTYdadjYd) + 2._dp*(YdadjTYuTYuadjYd) + 2._dp*(YdadjYdmd2YdadjYd) +     & 
&  2._dp*(YdadjYdTYdadjTYd) + YdadjYdYdadjYdmd2 + 2._dp*(YdadjYdYdmq2adjYd) +            & 
&  2._dp*(YdadjYumu2YuadjYd) + 2._dp*(YdadjYuTYuadjTYd) + YdadjYuYuadjYdmd2 +            & 
&  2._dp*(YdadjYuYumq2adjYd) + 2._dp*(Ydmq2adjYdYdadjYd) + 2._dp*(Ydmq2adjYuYuadjYd)) -  & 
&  12*g2p2*TpTYx3CYx3*Conjg(MassWB) - 12*g2p2*TYdadjYd*Conjg(MassWB) + (32*g3p2*id3R*Conjg(MassG)*(-& 
& 60*g3p2*MassG + 2*g1p2*(MassB + 2._dp*(MassG)) + 45*g3p2*MassG*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)))/45._dp + (4*g1p2*Conjg(MassB)*(-      & 
& 45*(5._dp*(TpTYx3CYx3) + TYdadjYd - 2*MassB*(5._dp*(TpYx3CYx3) + YdadjYd)) +           & 
&  id3R*(606*g1p2*MassB + 80*g3p2*(2._dp*(MassB) + MassG) + 225*g1p2*MassB*(NGHx3 +& 
&  NGHxb3))))/225._dp + (8*g1p4*id3R*Tr2(1))/15._dp + (32*g3p4*id3R*Tr2(3))/3._dp  
betamd22 =  betamd22+ (8*g1p2*id3R*Tr3(1))/3._dp

 
Dmd2 = oo16pi2*( betamd21 + oo16pi2 * betamd22 ) 

 
Else 
Dmd2 = oo16pi2* betamd21 
End If 
 
 
Forall(i1=1:3) Dmd2(i1,i1) =  Real(Dmd2(i1,i1),dp) 
!-------------------- 
! mu2 
!-------------------- 
 
betamu21  = 2*(mu2YuadjYu + 2._dp*(TYuadjTYu) + 2*mHu2*YuadjYu + YuadjYumu2 +         & 
&  2._dp*(Yumq2adjYu)) - (4*id3R*(40*AbsMassG*g3p2 + g1p2*(8._dp*(AbsMassB)              & 
&  + 5*Tr1(1))))/15._dp

 
 
If (TwoLoopRGE) Then 
betamu22 = (-2*g1p2*mu2YuadjYu)/5._dp + 6*g2p2*mu2YuadjYu - (3*mu2YuadjYu*TrYb3adjYb3)/5._dp -   & 
&  6*mu2YuadjYu*TrYuadjYu - 3*mu2YuadjYu*TrYw3adjYw3 - 6*mu2YuadjYu*TrYx3adjYx3 -        & 
&  (4*g1p2*TYuadjTYu)/5._dp + 12*g2p2*TYuadjTYu - (6*TrYb3adjYb3*TYuadjTYu)/5._dp -      & 
&  12*TrYuadjYu*TYuadjTYu - 6*TrYw3adjYw3*TYuadjTYu - 12*TrYx3adjYx3*TYuadjTYu -         & 
&  (6*TrCTYb3TpYb3*TYuadjYu)/5._dp - 12*TrCTYuTpYu*TYuadjYu - 6*TrCTYw3TpYw3*TYuadjYu -  & 
&  12*TrCTYx3TpYx3*TYuadjYu + (4*g1p2*MassB*YuadjTYu)/5._dp - 12*g2p2*MassWB*YuadjTYu -  & 
&  (6*TradjYb3TYb3*YuadjTYu)/5._dp - 12*TradjYuTYu*YuadjTYu - 6*TradjYw3TYw3*YuadjTYu -  & 
&  12*TradjYx3TYx3*YuadjTYu - 4*mHu2*YuadjYdYdadjYu + 24*AbsMassWB*g2p2*YuadjYu -        & 
&  (4*g1p2*mHu2*YuadjYu)/5._dp + 12*g2p2*mHu2*YuadjYu - (6*TrCTYb3TpTYb3*YuadjYu)/5._dp -& 
&  12*TrCTYuTpTYu*YuadjYu - 6*TrCTYw3TpTYw3*YuadjYu - 12*TrCTYx3TpTYx3*YuadjYu -         & 
&  12*Trmd2adjYx3Yx3*YuadjYu - (6*TrmHb32Yb3adjYb3*YuadjYu)/5._dp - 6*TrmHw32Yw3adjYw3*YuadjYu -& 
&  12*TrmHxb32Yx3adjYx3*YuadjYu - (6*Trml2adjYb3Yb3*YuadjYu)/5._dp - 6*Trml2adjYw3Yw3*YuadjYu -& 
&  12*Trmq2adjYuYu*YuadjYu - 12*Trmu2YuadjYu*YuadjYu - (12*mHu2*TrYb3adjYb3*YuadjYu)/5._dp -& 
&  24*mHu2*TrYuadjYu*YuadjYu - 12*mHu2*TrYw3adjYw3*YuadjYu - 24*mHu2*TrYx3adjYx3*YuadjYu -& 
&  (2*g1p2*YuadjYumu2)/5._dp + 6*g2p2*YuadjYumu2 - (3*TrYb3adjYb3*YuadjYumu2)/5._dp -    & 
&  6*TrYuadjYu*YuadjYumu2 - 3*TrYw3adjYw3*YuadjYumu2 - 6*TrYx3adjYx3*YuadjYumu2 -        & 
&  8*mHu2*YuadjYuYuadjYu - (4*g1p2*Yumq2adjYu)/5._dp + 12*g2p2*Yumq2adjYu -              & 
&  (6*TrYb3adjYb3*Yumq2adjYu)/5._dp - 12*TrYuadjYu*Yumq2adjYu - 6*TrYw3adjYw3*Yumq2adjYu -& 
&  12*TrYx3adjYx3*Yumq2adjYu - 2*(mu2YuadjYdYdadjYu + mu2YuadjYuYuadjYu + 2._dp*(TYuadjYdYdadjTYu) +& 
&  2._dp*(TYuadjYuYuadjTYu) + 2._dp*(TYuTpYdCTYdadjYu) + 2._dp*(TYuTpYuCTYuadjYu) +      & 
&  2._dp*(YuadjTYdTYdadjYu) + 2._dp*(YuadjTYuTYuadjYu) + 2._dp*(YuadjYdmd2YdadjYu) +     & 
&  2._dp*(YuadjYdTYdadjTYu) + 2*mHd2*YuadjYdYdadjYu + YuadjYdYdadjYumu2 + 2._dp*(YuadjYdYdmq2adjYu) +& 
&  2._dp*(YuadjYuTYuadjTYu) + YuadjYuYuadjYumu2 + 2*(YuadjYumu2YuadjYu + YuadjYuYumq2adjYu) +& 
&  2._dp*(Yumq2adjYdYdadjYu) + 2._dp*(Yumq2adjYuYuadjYu)) - 12*g2p2*TYuadjYu*Conjg(MassWB)  
betamu22 =  betamu22+ (32*g3p2*id3R*Conjg(MassG)*(-60*g3p2*MassG + 8*g1p2*(MassB + 2._dp*(MassG)) +         & 
&  45*g3p2*MassG*(3*NGHg3 + NGHx3 + NGHxb3)))/45._dp +& 
&  (4*g1p2*Conjg(MassB)*(45*(TYuadjYu - 2*MassB*YuadjYu) + 4*id3R*(642*g1p2*MassB +      & 
&  80*g3p2*(2._dp*(MassB) + MassG) + 225*g1p2*MassB*(NGHx3 +              & 
&  NGHxb3))))/225._dp + (32*g1p4*id3R*Tr2(1))/15._dp + (32*g3p4*id3R*Tr2(3))/3._dp -& 
&  (16*g1p2*id3R*Tr3(1))/3._dp

 
Dmu2 = oo16pi2*( betamu21 + oo16pi2 * betamu22 ) 

 
Else 
Dmu2 = oo16pi2* betamu21 
End If 
 
 
Forall(i1=1:3) Dmu2(i1,i1) =  Real(Dmu2(i1,i1),dp) 
!-------------------- 
! me2 
!-------------------- 
 
betame21  = 2*(me2YeadjYe + 2._dp*(TYeadjTYe) + 2*mHd2*YeadjYe + YeadjYeme2 +         & 
&  2._dp*(Yeml2adjYe)) + (2*g1p2*id3R*(-12._dp*(AbsMassB) + 5*Tr1(1)))/5._dp

 
 
If (TwoLoopRGE) Then 
betame22 = (12*g1p2*Conjg(MassB)*(5*(TYeadjYe - 2*MassB*YeadjYe) + 3*g1p2*id3R*MassB*(78 +       & 
&  25*NGHx3 + 25*NGHxb3)))/25._dp + (-3._dp*(me2YeadjYb3Yb3adjYe) -& 
&  10._dp*(me2YeadjYeYeadjYe) - 15._dp*(me2YeadjYw3Yw3adjYe) - 6._dp*(TYeadjYb3Yb3adjTYe) -& 
&  20._dp*(TYeadjYeYeadjTYe) - 30._dp*(TYeadjYw3Yw3adjTYe) - 6._dp*(TYeTpYb3CTYb3adjYe) -& 
&  20._dp*(TYeTpYeCTYeadjYe) - 30._dp*(TYeTpYw3CTYw3adjYe) - 6._dp*(YeadjTYb3TYb3adjYe) -& 
&  20._dp*(YeadjTYeTYeadjYe) - 30._dp*(YeadjTYw3TYw3adjYe) - 6._dp*(YeadjYb3mHb32Yb3adjYe) -& 
&  6._dp*(YeadjYb3TYb3adjTYe) - 3._dp*(YeadjYb3Yb3adjYeme2) - 6._dp*(YeadjYb3Yb3ml2adjYe) -& 
&  4*(-30*AbsMassWB*g2p2 + 3*g1p2*mHd2 - 15*g2p2*mHd2 + 15._dp*(TrCTYdTpTYd) +           & 
&  5._dp*(TrCTYeTpTYe) + 15._dp*(Trmd2YdadjYd) + 5._dp*(Trme2YeadjYe) + 5._dp*(Trml2adjYeYe) +& 
&  15._dp*(Trmq2adjYdYd) + 10*mHd2*(3._dp*(TrYdadjYd) + TrYeadjYe))*YeadjYe -            & 
&  20._dp*(YeadjYeTYeadjTYe) - 30._dp*(YeadjYw3TYw3adjTYe) - 2*(3*(mHd2 + mHu2)*YeadjYb3Yb3adjYe +& 
&  20*mHd2*YeadjYeYeadjYe + 15*(mHd2 + mHu2)*YeadjYw3Yw3adjYe) - 5*(4._dp*(YeadjYeme2YeadjYe) +& 
&  2._dp*(YeadjYeYeadjYeme2) + 4._dp*(YeadjYeYeml2adjYe) + 3._dp*(YeadjYw3Yw3adjYeme2) + & 
&  6*(YeadjYw3mHw32Yw3adjYe + YeadjYw3Yw3ml2adjYe)) - 6._dp*(Yeml2adjYb3Yb3adjYe) -      & 
&  20._dp*(Yeml2adjYeYeadjYe) - 30._dp*(Yeml2adjYw3Yw3adjYe) - 2*(10*(3._dp*(TrCTYdTpYd) +& 
&  TrCTYeTpYe)*TYeadjYe + (-6*g1p2*MassB + 30*g2p2*MassWB + 30._dp*(TradjYdTYd) +        & 
&  10._dp*(TradjYeTYe))*YeadjTYe + (3*(g1p2 - 5._dp*(g2p2) + 5._dp*(TrYdadjYd)) +        & 
&  5._dp*(TrYeadjYe))*(me2YeadjYe + 2._dp*(TYeadjTYe) + YeadjYeme2 + 2._dp*(Yeml2adjYe)) +& 
&  30*g2p2*TYeadjYe*Conjg(MassWB)) + 8*id3R*(3*g1p4*Tr2(1) + 5*g1p2*Tr3(1)))/5._dp

 
Dme2 = oo16pi2*( betame21 + oo16pi2 * betame22 ) 

 
Else 
Dme2 = oo16pi2* betame21 
End If 
 
 
Forall(i1=1:3) Dme2(i1,i1) =  Real(Dme2(i1,i1),dp) 
!-------------------- 
! mHw32 
!-------------------- 
 
betamHw321  = -16*AbsMassWB*g2p2*id3R + mHw32Yw3adjYw3 + 2._dp*(TYw3adjTYw3)          & 
&  + 2*mHu2*Yw3adjYw3 + Yw3adjYw3mHw32 + 2._dp*(Yw3ml2adjYw3)

 
 
If (TwoLoopRGE) Then 
betamHw322 = -3._dp*(mHw32Yw3adjYb3Yb3adjYw3)/10._dp - mHw32Yw3adjYeYeadjYw3 + (3*g1p2*mHw32Yw3adjYw3)/5._dp -& 
&  g2p2*mHw32Yw3adjYw3 - 3._dp*(mHw32Yw3adjYw3Yw3adjYw3)/2._dp - (3*mHw32Yw3adjYw3*TrYb3adjYb3)/10._dp -& 
&  3*mHw32Yw3adjYw3*TrYuadjYu - (3*mHw32Yw3adjYw3*TrYw3adjYw3)/2._dp - 3*mHw32Yw3adjYw3*TrYx3adjYx3 +& 
&  (6*g1p2*TYw3adjTYw3)/5._dp - 2*g2p2*TYw3adjTYw3 - (3*TrYb3adjYb3*TYw3adjTYw3)/5._dp - & 
&  6*TrYuadjYu*TYw3adjTYw3 - 3*TrYw3adjYw3*TYw3adjTYw3 - 6*TrYx3adjYx3*TYw3adjTYw3 -     & 
&  3._dp*(TYw3adjYb3Yb3adjTYw3)/5._dp - 2._dp*(TYw3adjYeYeadjTYw3) - (3*TrCTYb3TpYb3*TYw3adjYw3)/5._dp -& 
&  6*TrCTYuTpYu*TYw3adjYw3 - 3*TrCTYw3TpYw3*TYw3adjYw3 - 6*TrCTYx3TpYx3*TYw3adjYw3 -     & 
&  3._dp*(TYw3adjYw3Yw3adjTYw3) - 3._dp*(TYw3TpYb3CTYb3adjYw3)/5._dp - 2._dp*(TYw3TpYeCTYeadjYw3) -& 
&  3._dp*(TYw3TpYw3CTYw3adjYw3) - 3._dp*(Yw3adjTYb3TYb3adjYw3)/5._dp - 2._dp*(Yw3adjTYeTYeadjYw3) -& 
&  (6*g1p2*MassB*Yw3adjTYw3)/5._dp + 2*g2p2*MassWB*Yw3adjTYw3 - (3*TradjYb3TYb3*Yw3adjTYw3)/5._dp -& 
&  6*TradjYuTYu*Yw3adjTYw3 - 3*TradjYw3TYw3*Yw3adjTYw3 - 6*TradjYx3TYx3*Yw3adjTYw3 -     & 
&  3._dp*(Yw3adjTYw3TYw3adjYw3) - 3._dp*(Yw3adjYb3mHb32Yb3adjYw3)/5._dp - 3._dp*(Yw3adjYb3TYb3adjTYw3)/5._dp -& 
&  (6*mHu2*Yw3adjYb3Yb3adjYw3)/5._dp - 3._dp*(Yw3adjYb3Yb3adjYw3mHw32)/10._dp -          & 
&  3._dp*(Yw3adjYb3Yb3ml2adjYw3)/5._dp - 2._dp*(Yw3adjYeme2YeadjYw3) - 2._dp*(Yw3adjYeTYeadjTYw3) -& 
&  2*mHd2*Yw3adjYeYeadjYw3 - 2*mHu2*Yw3adjYeYeadjYw3 - Yw3adjYeYeadjYw3mHw32 -           & 
&  2._dp*(Yw3adjYeYeml2adjYw3) + (12*AbsMassB*g1p2*Yw3adjYw3)/5._dp + (6*g1p2*mHu2*Yw3adjYw3)/5._dp -& 
&  2*g2p2*mHu2*Yw3adjYw3 - (3*TrCTYb3TpTYb3*Yw3adjYw3)/5._dp - 6*TrCTYuTpTYu*Yw3adjYw3 - & 
&  3*TrCTYw3TpTYw3*Yw3adjYw3 - 6*TrCTYx3TpTYx3*Yw3adjYw3 - 6*Trmd2adjYx3Yx3*Yw3adjYw3 -  & 
&  (3*TrmHb32Yb3adjYb3*Yw3adjYw3)/5._dp - 3*TrmHw32Yw3adjYw3*Yw3adjYw3 - 6*TrmHxb32Yx3adjYx3*Yw3adjYw3 -& 
&  (3*Trml2adjYb3Yb3*Yw3adjYw3)/5._dp - 3*Trml2adjYw3Yw3*Yw3adjYw3 - 6*Trmq2adjYuYu*Yw3adjYw3 -& 
&  6*Trmu2YuadjYu*Yw3adjYw3 - (6*mHu2*TrYb3adjYb3*Yw3adjYw3)/5._dp - 12*mHu2*TrYuadjYu*Yw3adjYw3  
betamHw322 =  betamHw322- 6*mHu2*TrYw3adjYw3*Yw3adjYw3 - 12*mHu2*TrYx3adjYx3*Yw3adjYw3 + (3*g1p2*Yw3adjYw3mHw32)/5._dp -& 
&  g2p2*Yw3adjYw3mHw32 - (3*TrYb3adjYb3*Yw3adjYw3mHw32)/10._dp - 3*TrYuadjYu*Yw3adjYw3mHw32 -& 
&  (3*TrYw3adjYw3*Yw3adjYw3mHw32)/2._dp - 3*TrYx3adjYx3*Yw3adjYw3mHw32 - 3._dp*(Yw3adjYw3mHw32Yw3adjYw3) -& 
&  3._dp*(Yw3adjYw3TYw3adjTYw3) - 6*mHu2*Yw3adjYw3Yw3adjYw3 - 3._dp*(Yw3adjYw3Yw3adjYw3mHw32)/2._dp -& 
&  3._dp*(Yw3adjYw3Yw3ml2adjYw3) - 3._dp*(Yw3ml2adjYb3Yb3adjYw3)/5._dp - 2._dp*(Yw3ml2adjYeYeadjYw3) +& 
&  (6*g1p2*Yw3ml2adjYw3)/5._dp - 2*g2p2*Yw3ml2adjYw3 - (3*TrYb3adjYb3*Yw3ml2adjYw3)/5._dp -& 
&  6*TrYuadjYu*Yw3ml2adjYw3 - 3*TrYw3adjYw3*Yw3ml2adjYw3 - 6*TrYx3adjYx3*Yw3ml2adjYw3 -  & 
&  3._dp*(Yw3ml2adjYw3Yw3adjYw3) - (6*g1p2*TYw3adjYw3*Conjg(MassB))/5._dp +              & 
&  2*g2p2*Conjg(MassWB)*(TYw3adjYw3 - 2*MassWB*Yw3adjYw3 + 4*g2p2*id3R*MassWB*(26 +      & 
&  12*NGHw3 + 9*NGHx3 + 9*NGHxb3)) +        & 
&  16*g2p4*id3R*Tr2(2)

 
DmHw32 = oo16pi2*( betamHw321 + oo16pi2 * betamHw322 ) 

 
Else 
DmHw32 = oo16pi2* betamHw321 
End If 
 
 
Forall(i1=1:3) DmHw32(i1,i1) =  Real(DmHw32(i1,i1),dp) 
!-------------------- 
! mHg32 
!-------------------- 
 
betamHg321  = -24*AbsMassG*g3p2*id3R

 
 
If (TwoLoopRGE) Then 
betamHg322 = 24*g3p4*id3R*(3*AbsMassG*(2 + 3*NGHg3 + NGHx3 +         & 
&  NGHxb3) + Tr2(3))

 
DmHg32 = oo16pi2*( betamHg321 + oo16pi2 * betamHg322 ) 

 
Else 
DmHg32 = oo16pi2* betamHg321 
End If 
 
 
Forall(i1=1:3) DmHg32(i1,i1) =  Real(DmHg32(i1,i1),dp) 
!-------------------- 
! mHb32 
!-------------------- 
 
betamHb321  = (3*(mHb32Yb3adjYb3 + 2._dp*(TYb3adjTYb3) + 2*mHu2*Yb3adjYb3 +           & 
&  Yb3adjYb3mHb32 + 2._dp*(Yb3ml2adjYb3)))/5._dp

 
 
If (TwoLoopRGE) Then 
betamHb322 = (3*(-3._dp*(mHb32Yb3adjYb3Yb3adjYb3) - 10._dp*(mHb32Yb3adjYeYeadjYb3) -               & 
&  15._dp*(mHb32Yb3adjYw3Yw3adjYb3) - 3*mHb32Yb3adjYb3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) +& 
&  5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3))) - 6*(TrYb3adjYb3 +& 
&  10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) - 5._dp*(TrYx3adjYx3)))*TYb3adjTYb3 -& 
&  6._dp*(TYb3adjYb3Yb3adjTYb3) - 20._dp*(TYb3adjYeYeadjTYb3) - 30._dp*(TYb3adjYw3Yw3adjTYb3) -& 
&  6._dp*(TYb3TpYb3CTYb3adjYb3) - 20._dp*(TYb3TpYeCTYeadjYb3) - 30._dp*(TYb3TpYw3CTYw3adjYb3) -& 
&  6*(2*g1p2*MassB + 10*g2p2*MassWB + TradjYb3TYb3 + 10._dp*(TradjYuTYu) +               & 
&  5._dp*(TradjYw3TYw3) + 10._dp*(TradjYx3TYx3))*Yb3adjTYb3 - 6._dp*(Yb3adjTYb3TYb3adjYb3) -& 
&  20._dp*(Yb3adjTYeTYeadjYb3) - 30._dp*(Yb3adjTYw3TYw3adjYb3) + 6*(4*AbsMassB*g1p2 +    & 
&  20*AbsMassWB*g2p2 + 2*g1p2*mHu2 + 10*g2p2*mHu2 - TrCTYb3TpTYb3 - 10._dp*(TrCTYuTpTYu) -& 
&  5._dp*(TrCTYw3TpTYw3) - 10._dp*(TrCTYx3TpTYx3) - 10._dp*(Trmd2adjYx3Yx3) -            & 
&  TrmHb32Yb3adjYb3 - 5._dp*(TrmHw32Yw3adjYw3) - 10._dp*(TrmHxb32Yx3adjYx3) -            & 
&  Trml2adjYb3Yb3 - 5._dp*(Trml2adjYw3Yw3) - 10._dp*(Trmq2adjYuYu) - 10._dp*(Trmu2YuadjYu) -& 
&  2*mHu2*(TrYb3adjYb3 + 5*(2._dp*(TrYuadjYu) + TrYw3adjYw3 + 2._dp*(TrYx3adjYx3))))*Yb3adjYb3 -& 
&  3*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) -  & 
&  5._dp*(TrYx3adjYx3)))*Yb3adjYb3mHb32 - 6._dp*(Yb3adjYb3mHb32Yb3adjYb3) -              & 
&  6._dp*(Yb3adjYb3TYb3adjTYb3) - 3._dp*(Yb3adjYb3Yb3adjYb3mHb32) - 6._dp*(Yb3adjYb3Yb3ml2adjYb3) -& 
&  20._dp*(Yb3adjYeme2YeadjYb3) - 20._dp*(Yb3adjYeTYeadjTYb3) - 20*mHd2*Yb3adjYeYeadjYb3 -& 
&  10._dp*(Yb3adjYeYeadjYb3mHb32) - 20._dp*(Yb3adjYeYeml2adjYb3) - 30._dp*(Yb3adjYw3mHw32Yw3adjYb3) -& 
&  30._dp*(Yb3adjYw3TYw3adjTYb3) - 4*mHu2*(3._dp*(Yb3adjYb3Yb3adjYb3) + 5*(Yb3adjYeYeadjYb3 +& 
&  3._dp*(Yb3adjYw3Yw3adjYb3))) - 15._dp*(Yb3adjYw3Yw3adjYb3mHb32) - 30._dp*(Yb3adjYw3Yw3ml2adjYb3) -& 
&  6*(TrYb3adjYb3 + 10._dp*(TrYuadjYu) + 5._dp*(TrYw3adjYw3) - 2*(g1p2 + 5._dp*(g2p2) -  & 
&  5._dp*(TrYx3adjYx3)))*Yb3ml2adjYb3 - 6._dp*(Yb3ml2adjYb3Yb3adjYb3) - 20._dp*(Yb3ml2adjYeYeadjYb3) -& 
&  30._dp*(Yb3ml2adjYw3Yw3adjYb3) - 6*TYb3adjYb3*(TrCTYb3TpYb3 + 5*(2._dp*(TrCTYuTpYu) + & 
&  TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3)) + 2*g1p2*Conjg(MassB) + 10*g2p2*Conjg(MassWB))))/50._dp

 
DmHb32 = oo16pi2*( betamHb321 + oo16pi2 * betamHb322 ) 

 
Else 
DmHb32 = oo16pi2* betamHb321 
End If 
 
 
Forall(i1=1:3) DmHb32(i1,i1) =  Real(DmHb32(i1,i1),dp) 
!-------------------- 
! mHx32 
!-------------------- 
 
betamHx321  = -(id3R*(18*AbsMassWB*g2p2 + 32*AbsMassG*g3p2 + 5*g1p2*(2._dp*(AbsMassB) & 
&  - Tr1(1))))/3._dp

 
 
If (TwoLoopRGE) Then 
betamHx322 = (id3R*(g1p2*Conjg(MassB)*(669*g1p2*MassB + 5*(16*g3p2*(2._dp*(MassB) + MassG) +       & 
&  9*g2p2*(2._dp*(MassB) + MassWB)) + 225*g1p2*MassB*(NGHx3 +             & 
&  NGHxb3)) + 16*g3p2*Conjg(MassG)*(5*g1p2*(MassB + 2._dp*(MassG)) +      & 
&  3*(-8*g3p2*MassG + 3*g2p2*(2._dp*(MassG) + MassWB)) + 18*g3p2*MassG*(3*NGHg3 +& 
&  NGHx3 + NGHxb3)) + 3*(3*g2p2*Conjg(MassWB)*(33*g2p2*MassWB +& 
&  5*g1p2*(MassB + 2._dp*(MassWB)) + 16*g3p2*(MassG + 2._dp*(MassWB)) + 9*g2p2*MassWB*(4*NGHw3 +& 
&  3*(NGHx3 + NGHxb3))) + 2*(5*g1p4*Tr2(1) +               & 
&  9*g2p4*Tr2(2) + 16*g3p4*Tr2(3) + 10*g1p2*Tr3(1)))))/9._dp

 
DmHx32 = oo16pi2*( betamHx321 + oo16pi2 * betamHx322 ) 

 
Else 
DmHx32 = oo16pi2* betamHx321 
End If 
 
 
Forall(i1=1:3) DmHx32(i1,i1) =  Real(DmHx32(i1,i1),dp) 
!-------------------- 
! mHxb32 
!-------------------- 
 
betamHxb321  = mHxb32Yx3adjYx3 + 2._dp*(TYx3adjTYx3) + 2*mHu2*Yx3adjYx3 +             & 
&  Yx3adjYx3mHxb32 + 2._dp*(Yx3md2adjYx3) - (id3R*(18*AbsMassWB*g2p2 + 32*AbsMassG*g3p2 +& 
&  5*g1p2*(2._dp*(AbsMassB) + Tr1(1))))/3._dp

 
 
If (TwoLoopRGE) Then 
betamHxb322 = 10*AbsMassWB*g1p2*g2p2*id3R + 33*AbsMassWB*g2p4*id3R + 32*AbsMassWB*g2p2*g3p2*id3R -  & 
&  (2*g1p2*mHxb32Yx3adjYx3)/5._dp - (3*mHxb32Yx3adjYx3*TrYb3adjYb3)/10._dp -             & 
&  3*mHxb32Yx3adjYx3*TrYuadjYu - (3*mHxb32Yx3adjYx3*TrYw3adjYw3)/2._dp - 3*mHxb32Yx3adjYx3*TrYx3adjYx3 -& 
&  (4*g1p2*TYx3adjTYx3)/5._dp - (3*TrYb3adjYb3*TYx3adjTYx3)/5._dp - 6*TrYuadjYu*TYx3adjTYx3 -& 
&  3*TrYw3adjYw3*TYx3adjTYx3 - 6*TrYx3adjYx3*TYx3adjTYx3 - (3*(TrCTYb3TpYb3 +            & 
&  5*(2._dp*(TrCTYuTpYu) + TrCTYw3TpYw3 + 2._dp*(TrCTYx3TpYx3)))*TYx3adjYx3)/5._dp +     & 
&  (4*g1p2*MassB*Yx3adjTYx3)/5._dp - (3*TradjYb3TYb3*Yx3adjTYx3)/5._dp - 6*TradjYuTYu*Yx3adjTYx3 -& 
&  3*TradjYw3TYw3*Yx3adjTYx3 - 6*TradjYx3TYx3*Yx3adjTYx3 - (4*g1p2*mHu2*Yx3adjYx3)/5._dp -& 
&  (3*TrCTYb3TpTYb3*Yx3adjYx3)/5._dp - 6*TrCTYuTpTYu*Yx3adjYx3 - 3*TrCTYw3TpTYw3*Yx3adjYx3 -& 
&  6*TrCTYx3TpTYx3*Yx3adjYx3 - 6*Trmd2adjYx3Yx3*Yx3adjYx3 - (3*TrmHb32Yb3adjYb3*Yx3adjYx3)/5._dp -& 
&  3*TrmHw32Yw3adjYw3*Yx3adjYx3 - 6*TrmHxb32Yx3adjYx3*Yx3adjYx3 - (3*Trml2adjYb3Yb3*Yx3adjYx3)/5._dp -& 
&  3*Trml2adjYw3Yw3*Yx3adjYx3 - 6*Trmq2adjYuYu*Yx3adjYx3 - 6*Trmu2YuadjYu*Yx3adjYx3 -    & 
&  (6*mHu2*TrYb3adjYb3*Yx3adjYx3)/5._dp - 12*mHu2*TrYuadjYu*Yx3adjYx3 - 6*mHu2*TrYw3adjYw3*Yx3adjYx3 -& 
&  12*mHu2*TrYx3adjYx3*Yx3adjYx3 - (2*g1p2*Yx3adjYx3mHxb32)/5._dp - (3*TrYb3adjYb3*Yx3adjYx3mHxb32)/10._dp -& 
&  3*TrYuadjYu*Yx3adjYx3mHxb32 - (3*TrYw3adjYw3*Yx3adjYx3mHxb32)/2._dp - 3*TrYx3adjYx3*Yx3adjYx3mHxb32 -& 
&  8*mHu2*Yx3adjYx3Yx3adjYx3 - 4*mHu2*Yx3CYdTpYdadjYx3 - (4*g1p2*Yx3md2adjYx3)/5._dp -   & 
&  (3*TrYb3adjYb3*Yx3md2adjYx3)/5._dp - 6*TrYuadjYu*Yx3md2adjYx3 - 3*TrYw3adjYw3*Yx3md2adjYx3 -& 
&  6*TrYx3adjYx3*Yx3md2adjYx3 - 2*(mHxb32Yx3adjYx3Yx3adjYx3 + mHxb32Yx3CYdTpYdadjYx3 +   & 
&  2._dp*(TYx3adjYx3Yx3adjTYx3) + 2._dp*(TYx3CYdTpYdadjTYx3) + 2._dp*(TYx3TpYx3CTYx3adjYx3) +& 
&  2._dp*(TYx3YdadjTYdadjYx3) + 2._dp*(Yx3adjTYx3TYx3adjYx3) + 2._dp*(Yx3adjYx3mHxb32Yx3adjYx3) +& 
&  2._dp*(Yx3adjYx3TYx3adjTYx3) + Yx3adjYx3Yx3adjYx3mHxb32 + 2._dp*(Yx3adjYx3Yx3md2adjYx3) +& 
&  2._dp*(Yx3CTYdTpTYdadjYx3) + 2*mHd2*Yx3CYdTpYdadjYx3 + Yx3CYdTpYdadjYx3mHxb32 +       & 
&  2*(Yx3CYdmq2TpYdadjYx3 + Yx3CYdTpYdmd2adjYx3) + 2._dp*(Yx3md2adjYx3Yx3adjYx3) +       & 
&  2._dp*(Yx3md2CYdTpYdadjYx3) + 2._dp*(Yx3TYdadjYdadjTYx3)) + 5*g1p2*g2p2*id3R*MassB*Conjg(MassWB)  
betamHxb322 =  betamHxb322+ 16*g2p2*g3p2*id3R*MassG*Conjg(MassWB) + 36*AbsMassWB*g2p4*id3R*NGHw3 + & 
&  27*AbsMassWB*g2p4*id3R*NGHx3 + 27*AbsMassWB*g2p4*id3R*NGHxb3 +& 
&  (16*g3p2*id3R*Conjg(MassG)*(5*g1p2*(MassB + 2._dp*(MassG)) + 3*(-8*g3p2*MassG +       & 
&  3*g2p2*(2._dp*(MassG) + MassWB)) + 18*g3p2*MassG*(3*NGHg3 +            & 
&  NGHx3 + NGHxb3)))/9._dp + (g1p2*Conjg(MassB)*(36*(TYx3adjYx3 -& 
&  2*MassB*Yx3adjYx3) + 5*id3R*(669*g1p2*MassB + 5*(16*g3p2*(2._dp*(MassB) +             & 
&  MassG) + 9*g2p2*(2._dp*(MassB) + MassWB)) + 225*g1p2*MassB*(NGHx3 +    & 
&  NGHxb3))))/45._dp + (10*g1p4*id3R*Tr2(1))/3._dp + 6*g2p4*id3R*Tr2(2) + & 
&  (32*g3p4*id3R*Tr2(3))/3._dp - (20*g1p2*id3R*Tr3(1))/3._dp

 
DmHxb32 = oo16pi2*( betamHxb321 + oo16pi2 * betamHxb322 ) 

 
Else 
DmHxb32 = oo16pi2* betamHxb321 
End If 
 
 
Forall(i1=1:3) DmHxb32(i1,i1) =  Real(DmHxb32(i1,i1),dp) 
!-------------------- 
! MassB 
!-------------------- 
 
betaMassB1  = (g1p2*MassB*(66 + 25*NGHx3 + 25*NGHxb3))/5._dp

 
 
If (TwoLoopRGE) Then 
betaMassB2 = (g1p2*(6*(398*g1p2*MassB + 5*(88*g3p2*(MassB + MassG) + 27*g2p2*(MassB +              & 
&  MassWB)) + 9._dp*(TradjYb3TYb3) + 70._dp*(TradjYdTYd) + 90._dp*(TradjYeTYe) +         & 
&  130._dp*(TradjYuTYu) + 60._dp*(TradjYw3TYw3) + 190._dp*(TradjYx3TYx3) -               & 
&  9*MassB*TrYb3adjYb3 - 10*MassB*(7._dp*(TrYdadjYd) + 9._dp*(TrYeadjYe) +               & 
&  13._dp*(TrYuadjYu) + 6._dp*(TrYw3adjYw3) + 19._dp*(TrYx3adjYx3))) + 125*(10*g1p2*MassB +& 
&  16*g3p2*(MassB + MassG) + 9*g2p2*(MassB + MassWB))*NGHx3 +             & 
&  125*(10*g1p2*MassB + 16*g3p2*(MassB + MassG) + 9*g2p2*(MassB + MassWB))*NGHxb3))/75._dp

 
DMassB = oo16pi2*( betaMassB1 + oo16pi2 * betaMassB2 ) 

 
Else 
DMassB = oo16pi2* betaMassB1 
End If 
 
 
!-------------------- 
! MassWB 
!-------------------- 
 
betaMassWB1  = g2p2*MassWB*(2 + 4*NGHw3 + 3*NGHx3 +     & 
&  3*NGHxb3)

 
 
If (TwoLoopRGE) Then 
betaMassWB2 = (g2p2*(2*(9._dp*(TradjYb3TYb3) + 90._dp*(TradjYdTYd) + 30._dp*(TradjYeTYe) +          & 
&  90._dp*(TradjYuTYu) + 140._dp*(TradjYw3TYw3) + 3*(9*g1p2*(MassB + MassWB) +           & 
&  10*(25*g2p2*MassWB + 12*g3p2*(MassG + MassWB)) + 30._dp*(TradjYx3TYx3))) -            & 
&  18*MassWB*TrYb3adjYb3 - 180*MassWB*TrYdadjYd - 60*MassWB*TrYeadjYe - 180*MassWB*TrYuadjYu -& 
&  280*MassWB*TrYw3adjYw3 - 180*MassWB*TrYx3adjYx3 + 1440*g2p2*MassWB*NGHw3 +& 
&  15*(5*g1p2*(MassB + MassWB) + 2*(21*g2p2*MassWB + 8*g3p2*(MassG + MassWB)))*NGHx3 +& 
&  75*g1p2*MassB*NGHxb3 + 240*g3p2*MassG*NGHxb3 +          & 
&  75*g1p2*MassWB*NGHxb3 + 630*g2p2*MassWB*NGHxb3 +        & 
&  240*g3p2*MassWB*NGHxb3))/15._dp

 
DMassWB = oo16pi2*( betaMassWB1 + oo16pi2 * betaMassWB2 ) 

 
Else 
DMassWB = oo16pi2* betaMassWB1 
End If 
 
 
!-------------------- 
! MassG 
!-------------------- 
 
betaMassG1  = 2*g3p2*MassG*(-3 + 3*NGHg3 + NGHx3 +      & 
&  NGHxb3)

 
 
If (TwoLoopRGE) Then 
betaMassG2 = (2*g3p2*(3*(11*g1p2*(MassB + MassG) + 5*(28*g3p2*MassG + 9*g2p2*(MassG +              & 
&  MassWB)) + 20._dp*(TradjYdTYd) + 20._dp*(TradjYuTYu) + 20._dp*(TradjYx3TYx3)) +       & 
&  1620*g3p2*MassG*NGHg3 + 5*(68*g3p2*MassG + 5*g1p2*(MassB +             & 
&  MassG) + 9*g2p2*(MassG + MassWB))*NGHx3 + 5*(-12*MassG*(TrYdadjYd +    & 
&  TrYuadjYu + TrYx3adjYx3) + (68*g3p2*MassG + 5*g1p2*(MassB + MassG) + 9*g2p2*(MassG +  & 
&  MassWB))*NGHxb3)))/15._dp

 
DMassG = oo16pi2*( betaMassG1 + oo16pi2 * betaMassG2 ) 

 
Else 
DMassG = oo16pi2* betaMassG1 
End If 

!-------------------- 
! MnuL 
!-------------------- 
 

Do i1 = 1,3
Do i2 = 1,3
gammaL1(i1,i2) =  3._dp*adjYb3Yb3(i1,i2)/10._dp + adjYeYe(i1,i2) + 3._dp*adjYw3Yw3(i1,i2)/2._dp

If(i1.eq.i2) Then
 gammaL1(i1,i2) = + gammaL1(i1,i2) -3._dp*g1p2/10._dp - 3._dp*g2p2/2._dp
End If

End Do
End Do

gammaHd1 = -3._dp*g1p2/10._dp - 3._dp*g2p2/2._dp + 3._dp*TrYb3adjYb3/10._dp + &
& 3._dp*TrYuadjYu + 3._dp*TrYw3adjYw3/2._dp + 3._dp*TrYx3adjYx3

betaMnuL1 = MatMul(Transpose(gammaL1),MnuL) + MatMul(MnuL,gammaL1) + 2._dp*MnuL*gammaHd1

 
If (TwoLoopRGE) Then 

Do i1 = 1,3
Do i2 = 1,3
gammaL2(i1,i2) = -9._dp*TrYb3adjYb3*adjYb3Yb3(i1,i2)/100._dp - &
& 9._dp*TrYuadjYu*adjYb3Yb3(i1,i2)/10._dp - &
& 9._dp*TrYw3adjYw3*adjYb3Yb3(i1,i2)/20._dp - &
& 9._dp*TrYx3adjYx3*adjYb3Yb3(i1,i2)/10._dp - &
& 9._dp*adjYb3Yb3adjYb3Yb3(i1,i2)/50._dp - &
& 3._dp*adjYb3Yb3adjYw3Yw3(i1,i2)/10._dp + 6._dp*g1p2*adjYeYe(i1,i2)/5._dp - &
& 3._dp*TrYdadjYd*adjYeYe(i1,i2) - TrYeadjYe*adjYeYe(i1,i2) - &
& 2._dp*adjYeYeadjYeYe(i1,i2) + 6._dp*g2p2*adjYw3Yw3(i1,i2) - &
& 9._dp*TrYb3adjYb3*adjYw3Yw3(i1,i2)/20._dp - &
& 9._dp*TrYuadjYu*adjYw3Yw3(i1,i2)/2._dp - &
& 9._dp*TrYw3adjYw3*adjYw3Yw3(i1,i2)/4._dp - &
& 9._dp*TrYx3adjYx3*adjYw3Yw3(i1,i2)/2._dp - &
& 9._dp*adjYw3Yw3adjYb3Yb3(i1,i2)/40._dp - 3._dp*adjYw3Yw3adjYw3Yw3(i1,i2)/2._dp

If(i1.eq.i2) Then
 gammaL2(i1,i2) =  gammaL2(i1,i2) + 657._dp*g1**4/100._dp + 9._dp*g1p2*g2p2/10._dp &
	& + 105._dp*g2**4/4._dp
End If

End Do
End Do

gammaHd2 = 657._dp*g1**4/100._dp + 9._dp*g1p2*g2p2/10._dp + 105._dp*g2**4/4._dp - &
& 27._dp*TrYb3adjYb3Yb3adjYb3/100._dp - 3._dp*TrYb3adjYeYeadjYb3/10._dp - &
& 3._dp*TrYb3adjYw3Yw3adjYb3 /4._dp - 6._dp*TrYdadjYdTpYx3CYx3 - &
& 3._dp*TrYdadjYuYuadjYd + 4._dp*g1p2*TrYuadjYu/5._dp + 16._dp*g3p2*TrYuadjYu &
& - 9._dp*TrYuadjYuYuadjYu - 27._dp*TrYb3adjYw3Yw3adjYb3 /40._dp - &
& 3._dp*TrYeadjYw3Yw3adjYe/2._dp + 6._dp*g2p2*TrYw3adjYw3 - &
& 15._dp*TrYw3adjYw3Yw3adjYw3/4._dp + 2._dp*g1p2*TrYx3adjYx3 + &
& 16._dp*g3p2*TrYx3adjYx3 - 9._dp*TrYx3adjYx3Yx3adjYx3

betaMnuL2 = MatMul(Transpose(gammaL2),MnuL) + MatMul(MnuL,gammaL2) + 2._dp*MnuL*gammaHd2


 
DMnuL = oo16pi2*( betaMnuL1 + oo16pi2 * betaMnuL2 ) 

 
Else 
DMnuL = oo16pi2* betaMnuL1
End If 
 
 
If (ThresholdCrossed.lt.1) Then 
DMXM3(1,:) = 0._dp 
DBMXM3(1,:) = 0._dp 
DmHx32(1,:) = 0._dp 
DmHx32(:,1) = 0._dp 
DYx3(1,:) = 0._dp 
DTYx3(1,:) = 0._dp 
DMXM3(:,1) = 0._dp 
DBMXM3(:,1) = 0._dp 
DmHxb32(1,:) = 0._dp 
DmHxb32(:,1) = 0._dp 
DMGM3(1,:) = 0._dp 
DBMGM3(1,:) = 0._dp 
DMGM3(:,1) = 0._dp 
DBMGM3(:,1) = 0._dp 
DmHg32(1,:) = 0._dp 
DmHg32(:,1) = 0._dp 
DYb3(1,:) = 0._dp 
DTYb3(1,:) = 0._dp 
DMBM3(1,:) = 0._dp 
DBMBM3(1,:) = 0._dp 
DMBM3(:,1) = 0._dp 
DBMBM3(:,1) = 0._dp 
DmHb32(1,:) = 0._dp 
DmHb32(:,1) = 0._dp 
DYw3(1,:) = 0._dp 
DTYw3(1,:) = 0._dp 
DMWM3(1,:) = 0._dp 
DBMWM3(1,:) = 0._dp 
DMWM3(:,1) = 0._dp 
DBMWM3(:,1) = 0._dp 
DmHw32(1,:) = 0._dp 
DmHw32(:,1) = 0._dp 
End if 

If (ThresholdCrossed.lt.2) Then 
DMXM3(2,:) = 0._dp 
DBMXM3(2,:) = 0._dp 
DmHx32(2,:) = 0._dp 
DmHx32(:,2) = 0._dp 
DYx3(2,:) = 0._dp 
DTYx3(2,:) = 0._dp 
DMXM3(:,2) = 0._dp 
DBMXM3(:,2) = 0._dp 
DmHxb32(2,:) = 0._dp 
DmHxb32(:,2) = 0._dp 
DMGM3(2,:) = 0._dp 
DBMGM3(2,:) = 0._dp 
DMGM3(:,2) = 0._dp 
DBMGM3(:,2) = 0._dp 
DmHg32(2,:) = 0._dp 
DmHg32(:,2) = 0._dp 
DYb3(2,:) = 0._dp 
DTYb3(2,:) = 0._dp 
DMBM3(2,:) = 0._dp 
DBMBM3(2,:) = 0._dp 
DMBM3(:,2) = 0._dp 
DBMBM3(:,2) = 0._dp 
DmHb32(2,:) = 0._dp 
DmHb32(:,2) = 0._dp 
DYw3(2,:) = 0._dp 
DTYw3(2,:) = 0._dp 
DMWM3(2,:) = 0._dp 
DBMWM3(2,:) = 0._dp 
DMWM3(:,2) = 0._dp 
DBMWM3(:,2) = 0._dp 
DmHw32(2,:) = 0._dp 
DmHw32(:,2) = 0._dp 
End if 

If (ThresholdCrossed.lt.3) Then 
DMXM3(3,:) = 0._dp 
DBMXM3(3,:) = 0._dp 
DmHx32(3,:) = 0._dp 
DmHx32(:,3) = 0._dp 
DYx3(3,:) = 0._dp 
DTYx3(3,:) = 0._dp 
DMXM3(:,3) = 0._dp 
DBMXM3(:,3) = 0._dp 
DmHxb32(3,:) = 0._dp 
DmHxb32(:,3) = 0._dp 
DMGM3(3,:) = 0._dp 
DBMGM3(3,:) = 0._dp 
DMGM3(:,3) = 0._dp 
DBMGM3(:,3) = 0._dp 
DmHg32(3,:) = 0._dp 
DmHg32(:,3) = 0._dp 
DYb3(3,:) = 0._dp 
DTYb3(3,:) = 0._dp 
DMBM3(3,:) = 0._dp 
DBMBM3(3,:) = 0._dp 
DMBM3(:,3) = 0._dp 
DBMBM3(:,3) = 0._dp 
DmHb32(3,:) = 0._dp 
DmHb32(:,3) = 0._dp 
DYw3(3,:) = 0._dp 
DTYw3(3,:) = 0._dp 
DMWM3(3,:) = 0._dp 
DBMWM3(3,:) = 0._dp 
DMWM3(:,3) = 0._dp 
DBMWM3(:,3) = 0._dp 
DmHw32(3,:) = 0._dp 
DmHw32(:,3) = 0._dp 
End if 

Call ParametersToG555(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYb3,DYw3,DYx3,Dmue,DMXM3,              & 
& DMWM3,DMGM3,DMBM3,DTYu,DTYd,DTYe,DTYb3,DTYw3,DTYx3,DBmue,DBMXM3,DBMWM3,DBMGM3,         & 
& DBMBM3,Dmq2,Dml2,DmHd2,DmHu2,Dmd2,Dmu2,Dme2,DmHw32,DmHg32,DmHb32,DmHx32,               & 
& DmHxb32,DMassB,DMassWB,DMassG,DMnuL,f)

Iname = Iname - 1 
 
End Subroutine rge555

#endif SEESAWIII


 Subroutine Set_Decoupling_Heavy_States(set)
  Implicit none
  Logical, Intent(in) :: set

   decoupling_heavy_states = set

 End Subroutine Set_Decoupling_Heavy_States


End Module RGEs

